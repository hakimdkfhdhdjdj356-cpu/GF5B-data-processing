"""大气/辐射归一化后端注册表与调度。

四档后端，产物等级互不混淆：
- none       只交付辐亮度，不做任何归一化。
- iarr       纯 numpy/rasterio 的场景内相对反射率（默认）。零外部依赖。
- sixs_lut   本机开源 6SV1.1 + Py6S 生成 LUT，反演模型估计地表反射率。
- envi_quac  调用本机 ENVI QUAC 任务引擎（可选商业插件，需许可）。

说明：sixs_lut 的 6S 是 conda-forge 开源本机二进制（无商业许可），
因此 requires_external_software=False 指的是"无商业/许可软件依赖"；
它仍依赖外部二进制，见 requires_external_binary 与 external_dependencies。

设计约束：
- 只有真跑成功并产出文件，才允许流程声称有对应等级的产品。
- 相对反射率绝不叫"地表反射率"；QUAC 估计反射率也不声称经独立精度验证。
- 后端不可用时按配置降级到辐亮度层，不伪造任何反射率。
- 不读取、不修改任何许可文件。
"""
from __future__ import annotations

import json
import logging
import os
import subprocess
import tempfile
from pathlib import Path

from . import relative

log = logging.getLogger(__name__)

METHOD_INFO = {
    "none": {"product_grade": "radiance_only", "produces_analysis_product": False,
             "requires_external_software": False, "external_dependencies": []},
    "iarr": {"product_grade": "relative_reflectance_NOT_surface_reflectance",
             "produces_analysis_product": True, "requires_external_software": False,
             "external_dependencies": []},
    "sixs_lut": {"product_grade": "sixs_model_estimated_surface_reflectance_NOT_independently_validated",
                 "produces_analysis_product": True, "requires_external_software": False,
                 "requires_external_binary": True,
                 "external_dependencies": [
                     "6SV1.1 (conda-forge sixs, open-source local binary, no license)",
                     "Py6S (pip, LGPL)"]},
    "envi_quac": {"product_grade": "quac_estimated_reflectance_NOT_independently_validated",
                  "produces_analysis_product": True, "requires_external_software": True,
                  "external_dependencies": ["ENVI 5.6+ taskengine.exe", "ENVI/IDL license"]},
}


class BackendUnavailable(RuntimeError):
    pass


def _py6s_importable() -> bool:
    try:
        import Py6S  # noqa: F401
        return True
    except ImportError:
        return False


def detect(cfg: dict) -> dict:
    """探测当前配置下后端是否可用；iarr 恒可用，两种物理/估计后端依赖外部程序。"""
    a = cfg["atmospheric"]
    method = a["method"]
    info = METHOD_INFO[method]
    engine = str(a.get("envi_taskengine") or "")
    sixs_exe = str((a.get("sixs") or {}).get("executable") or "")
    if method == "envi_quac":
        available = bool(engine) and Path(engine).is_file()
    elif method == "sixs_lut":
        available = bool(sixs_exe) and Path(sixs_exe).is_file() and _py6s_importable()
    else:
        available = True
    notes = {
        "iarr": "无需外部软件，可无人值守重跑",
        "none": "不做大气校正，只交付辐亮度",
        "sixs_lut": ("本机开源 6SV1.1 + Py6S；无商业许可，但水汽/气溶胶/臭氧为"
                     "假设情景，产物为模型估计反射率，未经独立精度验证"),
        "envi_quac": "需本机安装 ENVI 并取得许可",
    }
    return {"method": method, "engine": engine or None,
            "sixs_executable": sixs_exe or None, "available": available,
            "sensor": a.get("sensor") if method == "envi_quac" else None,
            "note": notes[method], **info}


def _ascii_junction(target: Path, work: Path) -> Path:
    """IDL 对非 ASCII 路径不可靠：建一个临时目录联接指向任务目录，只传 ASCII 路径。"""
    link = work / "prepared"
    result = subprocess.run(["cmd.exe", "/d", "/c", "mklink", "/J", str(link), str(target)],
                            check=True, capture_output=True)
    if not link.is_dir():
        raise BackendUnavailable("目录联接创建失败: " + (result.stdout or b"").decode("utf-8", "replace"))
    return link


def _validate_output(output_dat: Path, input_dat: Path) -> dict:
    """用产物本身校验 QUAC 是否真的完成，而不是只信进程退出码。"""
    import rasterio
    info = {"exists": output_dat.is_file(), "readable": False, "bands": None,
            "size": None, "scale_factor": None, "bbl_ok": False, "nonzero_data": False}
    if not info["exists"]:
        return info
    try:
        with rasterio.open(output_dat) as q, rasterio.open(input_dat) as src:
            envi = q.tags(ns="ENVI")
            info["readable"] = True
            info["bands"] = q.count
            info["size"] = [q.width, q.height]
            info["expected_bands"] = src.count
            info["expected_size"] = [src.width, src.height]
            info["scale_factor"] = float(envi.get("reflectance_scale_factor", 0)) or None
            bbl = [int(x) for x in envi.get("bbl", "").strip("{} ").replace(",", " ").split()]
            info["bbl_ok"] = len(bbl) == q.count
            import numpy as np
            a = q.read(1, out_shape=(32, 32))
            info["nonzero_data"] = bool(np.any(a != 0))
    except Exception as exc:
        info["error"] = f"{type(exc).__name__}: {exc}"
    info["valid"] = bool(info["readable"] and info["bands"] == info.get("expected_bands")
                         and info["size"] == info.get("expected_size")
                         and info["scale_factor"] and info["bbl_ok"] and info["nonzero_data"])
    return info


def run_quac(cfg: dict, input_dat: Path, output_dat: Path, log_dir: Path) -> dict:
    """对单景执行 QUAC。返回结构化结果，异常不吞。"""
    status = detect(cfg)
    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    if not status["available"]:
        return {"successful": False, "backend_method": "envi_quac",
                "reason": "envi_taskengine_not_found", "engine": status["engine"]}

    engine = Path(status["engine"])
    input_dat = Path(input_dat).resolve()
    output_dat = Path(output_dat).resolve()
    output_dat.parent.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="ahsi_quac_"))
    link = _ascii_junction(output_dat.parents[1], work)
    # link 指向 <scene_dir>；ENVI_input / QUAC_v2 为其子目录
    rel_in = input_dat.relative_to(output_dat.parents[1])
    rel_out = output_dat.relative_to(output_dat.parents[1])
    ascii_in, ascii_out = link / rel_in, link / rel_out

    payload = {"taskName": "QUAC", "inputParameters": {
        "INPUT_RASTER": {"factory": "URLRaster", "url": str(ascii_in)},
        "SENSOR": cfg["atmospheric"]["sensor"],
        "OUTPUT_RASTER_URI": str(ascii_out)}}
    (log_dir / "task_parameters.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

    env = os.environ.copy()
    env["IDL_STARTUP"] = ""
    env["IDL_TMPDIR"] = str(work)
    record = {"method": f"envi_quac ({cfg['atmospheric']['sensor']})",
              "backend_method": "envi_quac", "engine": str(engine),
              "input": str(input_dat), "output": str(output_dat)}
    try:
        result = subprocess.run(
            [str(engine), f"--license-timeout={int(cfg['atmospheric']['license_timeout_s'])}", "ENVI"],
            input=json.dumps(payload), text=True, capture_output=True, cwd=work, env=env,
            timeout=float(cfg["atmospheric"]["timeout_s"]))
        (log_dir / "engine_stdout.txt").write_text(result.stdout or "", encoding="utf-8")
        (log_dir / "engine_stderr.txt").write_text(result.stderr or "", encoding="utf-8")
        integrity = _validate_output(output_dat, input_dat)
        accept = bool(cfg["atmospheric"].get("accept_nonzero_exit_with_valid_output", False))
        exit_ok = result.returncode == 0
        record.update(exit_code=result.returncode, output_file_created=output_dat.is_file(),
                      output_integrity=integrity,
                      engine_exit_nonzero_warning=bool(not exit_ok and integrity["valid"]),
                      successful=bool(integrity["valid"] and (exit_ok or accept)))
        if not record["successful"]:
            record["reason"] = ("engine_exit_nonzero_output_valid_but_rejected" if integrity["valid"]
                                else "engine_failed_invalid_output" if not exit_ok else "output_missing")
            record["stderr_head"] = (result.stderr or "")[:600]
        (log_dir / "output_integrity.json").write_text(
            json.dumps(integrity, ensure_ascii=False, indent=2), encoding="utf-8")
    except subprocess.TimeoutExpired:
        record.update(successful=False, reason="timeout",
                      timeout_s=cfg["atmospheric"]["timeout_s"])
    return record


def run_backend(cfg: dict, scene, cube_result: dict, scene_dir: Path,
                backend: dict, envi_input: Path | None = None) -> dict:
    """大气/归一化阶段统一入口。返回的 record 会被写入 stage_quac.json。"""
    scene_dir = Path(scene_dir)
    method = cfg["atmospheric"]["method"]
    if method == "none":
        return {"successful": False, "backend_method": "none",
                "reason": "method_none_deliver_radiance_only",
                "note": "配置选择只交付辐亮度立方体；这是有意行为，不是失败"}
    if method == "iarr":
        rec = relative.compute_factors(scene, cfg, cube_result, scene_dir / "IARR")
        rec["backend_method"] = "iarr"
        return rec
    if method == "sixs_lut":
        from . import sixs_lut
        if not backend["available"]:
            return {"successful": False, "backend_method": "sixs_lut",
                    "reason": "sixs_executable_or_py6s_unavailable",
                    "sixs_executable": backend.get("sixs_executable"),
                    "external_dependencies": backend["external_dependencies"]}
        rec = sixs_lut.build_lut(cfg, scene, cube_result, scene_dir / "SIXS_LUT")
        rec["backend_method"] = "sixs_lut"
        return rec
    if method == "envi_quac":
        if not backend["available"]:
            return {"successful": False, "backend_method": "envi_quac",
                    "reason": "envi_taskengine_not_found", "engine": backend["engine"],
                    "external_dependencies": backend["external_dependencies"]}
        if not envi_input or not Path(envi_input).is_file():
            return {"successful": False, "backend_method": "envi_quac",
                    "reason": "envi_input_missing", "envi_input": str(envi_input)}
        quac_dir = scene_dir / "QUAC_v2"
        output_dat = quac_dir / f"{scene.scene_id}_QUAC_reflectance.dat"
        rec = run_quac(cfg, Path(envi_input), output_dat, quac_dir)
        rec["backend_method"] = "envi_quac"
        return rec
    raise ValueError(f"未实现的大气后端: {method}")


def finalize(cfg: dict, scene, scene_dir: Path, cube_result: dict,
             stage_record: dict) -> dict:
    """成果收尾阶段统一入口：按后端方法产出不同等级、不同命名的产品。"""
    method = stage_record.get("backend_method") or cfg["atmospheric"]["method"]
    if method == "iarr":
        return relative.finalize_relative_reflectance(scene, cfg, scene_dir, cube_result, stage_record)
    if method == "sixs_lut":
        from . import sixs_lut
        return sixs_lut.finalize_sixs_reflectance(scene, cfg, scene_dir,
                                                  cube_result, stage_record)
    if method == "envi_quac":
        from . import core
        return core.finalize_reflectance(scene, cfg, scene_dir, Path(stage_record["output"]), cube_result)
    raise ValueError(f"该后端无成果收尾步骤: {method}")
