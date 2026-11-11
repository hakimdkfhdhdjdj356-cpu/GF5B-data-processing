"""阶段编排：run_id 隔离、断点续跑、输入清单、统一质量报告。

一个 run_id 对应一个自包含输出目录；换一景数据 / 换研究区只需改配置后新开 run_id，
旧结果永不被覆盖，也不依赖当前工作目录。
"""
from __future__ import annotations

import json
import logging
import platform
import shutil
import sys
from datetime import datetime
from pathlib import Path

from . import atmos, core
from .config import public_config
from .qc import evaluate_scene, overall, write_qc_report
from .scene import Scene, discover_scenes, file_inventory, scene_brief

log = logging.getLogger(__name__)

STAGES = ["radiance", "envi", "quac", "reflectance", "qc"]


def _stage_path(scene_dir: Path, stage: str) -> Path:
    return scene_dir / f"stage_{stage}.json"


def _done(scene_dir: Path, stage: str) -> bool:
    return _stage_path(scene_dir, stage).is_file()


def _mark(scene_dir: Path, stage: str, payload: dict) -> None:
    core.save_json(_stage_path(scene_dir, stage),
                   {"stage": stage, "finished_utc": datetime.utcnow().isoformat() + "Z", **payload})


def _load(stage_file: Path) -> dict:
    return json.loads(Path(stage_file).read_text(encoding="utf-8"))


def _free_gb(path: Path) -> float:
    path.mkdir(parents=True, exist_ok=True)
    return shutil.disk_usage(str(path)).free / (1024 ** 3)


def _environment(cfg: dict) -> dict:
    import numpy, rasterio
    versions = {"python": sys.version.split()[0], "platform": platform.platform(),
                "numpy": numpy.__version__, "rasterio": rasterio.__version__}
    try:
        from osgeo import gdal
        versions["gdal"] = gdal.__version__
    except Exception:
        pass
    versions["config_sha256"] = cfg.get("_config_sha256")
    versions["sensors_sha256"] = cfg.get("_sensors_sha256")
    return versions


def _validate_resume(run_dir: Path, cfg: dict, max_bands) -> None:
    """Do not overwrite the manifest and reuse stages from a different input/config or SMOKE."""
    p = run_dir / 'manifest.json'
    if p.is_file():
        old = _load(p)
        if old.get('config') != public_config(cfg) or old.get('smoke_max_bands') != max_bands:
            raise ValueError('Existing run has different input/config or smoke_max_bands; choose a new run_id')


def run(cfg: dict, run_id: str | None = None, scenes_filter: list[str] | None = None,
        stages: list[str] | None = None, max_bands: int | None = None,
        force: bool = False, progress=log.info) -> dict:
    stages = stages or list(STAGES)
    unknown = [s for s in stages if s not in STAGES]
    if unknown:
        raise ValueError(f"未知阶段 {unknown}；可选 {STAGES}")
    run_id = run_id or datetime.now().strftime("%Y%m%d_%H%M%S")
    if Path(run_id).name != run_id or run_id in ('.', '..'):
        raise ValueError('run_id must be a directory name, not a path')
    run_dir = Path(cfg["output_root"]) / run_id
    _validate_resume(run_dir, cfg, max_bands)
    run_dir.mkdir(parents=True, exist_ok=True)

    scenes, errors = discover_scenes(cfg["input_root"], cfg)
    if scenes_filter:
        keys = [k.upper() for k in scenes_filter]
        scenes = [s for s in scenes if any(k in s.scene_id.upper() or k == s.satellite.upper() for k in keys)]
    if not scenes:
        raise RuntimeError(f"未发现可处理场景。错误: {errors}")

    backend = atmos.detect(cfg)
    notes = [f"发现 {len(scenes)} 景；解析失败 {len(errors)} 项",
             f"大气校正后端 {backend['method']} 可用={backend['available']}"]
    if max_bands:
        notes.append(f"SMOKE 模式：每景仅输出前 {max_bands} 个波段，产物不可用于交付")
    if errors:
        notes.append(f"被拒绝的场景：{errors}")

    manifest = {"run_id": run_id, "created_utc": datetime.utcnow().isoformat() + "Z",
                "pipeline": "ahsi_pipeline v1", "stages_requested": stages, "smoke_max_bands": max_bands,
                "config": public_config(cfg), "environment": _environment(cfg),
                "backend": backend, "scenes": [], "discovery_errors": errors}
    core.save_json(run_dir / "manifest.json", manifest)

    results: dict[str, dict] = {}
    all_checks = []
    for scene in scenes:
        scene_dir = run_dir / "scenes" / scene.scene_id
        scene_dir.mkdir(parents=True, exist_ok=True)
        core.save_json(scene_dir / "scene_manifest.json",
                       {"scene": scene_brief(scene),
                        "input_files": file_inventory(scene, bool(cfg["runtime"]["hash_inputs"])),
                        "input_hashes_computed": bool(cfg["runtime"]["hash_inputs"])})
        r: dict = {}
        progress(f"[{scene.scene_id}] 开始")

        # --- 辐亮度 ---
        if "radiance" in stages:
            free = _free_gb(run_dir)
            if free < float(cfg["runtime"]["min_free_gb"]):
                raise RuntimeError(f"磁盘可用 {free:.1f} GB < 配置下限 {cfg['runtime']['min_free_gb']} GB")
            if _done(scene_dir, "radiance") and not force:
                r["radiance"] = _load(_stage_path(scene_dir, "radiance"))["result"]
                progress(f"[{scene.scene_id}] radiance 已完成，跳过")
            else:
                r["radiance"] = core.orthorectify(scene, cfg, scene_dir, max_bands, progress)
                _mark(scene_dir, "radiance", {"result": r["radiance"]})
        elif _done(scene_dir, "radiance"):
            r["radiance"] = _load(_stage_path(scene_dir, "radiance"))["result"]

        rad = r.get("radiance")
        # --- ENVI 输入（仅 envi_quac 后端需要；iarr / none 直接读 GeoTIFF 立方体） ---
        if "envi" in stages and rad and backend["method"] == "envi_quac":
            if _done(scene_dir, "envi") and not force:
                r["envi"] = _load(_stage_path(scene_dir, "envi"))["result"]
            else:
                p = core.export_envi(scene, cfg, scene_dir, rad)
                r["envi"] = {"path": str(p) if p else None, "status": "ok" if p else "skipped_no_candidate_vrt"}
                _mark(scene_dir, "envi", {"result": r["envi"]})
        elif _done(scene_dir, "envi"):
            r["envi"] = _load(_stage_path(scene_dir, "envi"))["result"]

        # --- 大气 / 归一化后端（阶段名 quac 保留以兼容既有 run 的断点续跑） ---
        method = backend["method"]
        stage_record = None
        if "quac" in stages and cfg["atmospheric"]["enabled"] and rad and method != "none":
            if _done(scene_dir, "quac") and not force:
                cached = _load(_stage_path(scene_dir, "quac"))["result"]
                old = cached.get("backend_method")
                if old and old != method:
                    progress(f"[{scene.scene_id}] 后端由 {old} 改为 {method}，旧阶段标记作废并重算")
                elif cached.get("successful"):
                    stage_record = cached
            if stage_record is None:
                envi_path = r.get("envi", {}).get("path") if method == "envi_quac" else None
                stage_record = atmos.run_backend(cfg, scene, rad, scene_dir, backend, envi_path)
                _mark(scene_dir, "quac", {"result": stage_record})
            r["quac"] = stage_record
            if not stage_record.get("successful"):
                reason = stage_record.get("reason")
                if cfg["atmospheric"]["on_unavailable"] == "fail":
                    raise RuntimeError(f"大气校正后端失败（{reason}），且 on_unavailable=fail")
                notes.append(f"{scene.scene_id}: 后端 {method} 未产出（{reason}），停在辐亮度层")
        elif method != "none" and _done(scene_dir, "quac"):
            stage_record = _load(_stage_path(scene_dir, "quac"))["result"]
            r["quac"] = stage_record

        # --- 成果收尾 ---
        if "reflectance" in stages and rad and stage_record and stage_record.get("successful"):
            grade = atmos.METHOD_INFO[method]["product_grade"]
            cached = None
            if _done(scene_dir, "reflectance") and not force:
                c = _load(_stage_path(scene_dir, "reflectance"))["result"]
                if c.get("product_grade", "quac_estimated_reflectance_NOT_independently_validated") == grade:
                    cached = c
            if cached is not None:
                r["reflectance"] = cached
            else:
                r["reflectance"] = atmos.finalize(cfg, scene, scene_dir, rad, stage_record)
                _mark(scene_dir, "reflectance", {"result": r["reflectance"]})
        elif (method != "none" and (stage_record is None or stage_record.get("successful"))
              and _done(scene_dir, "reflectance")):
            r["reflectance"] = _load(_stage_path(scene_dir, "reflectance"))["result"]

        # --- 质量检查 ---
        checks = evaluate_scene(scene, cfg, r.get("radiance"), r.get("reflectance"),
                                r.get("quac"), backend)
        all_checks.extend(checks)
        core.save_json(scene_dir / "scene_qc.json",
                       {"overall": overall(checks), "checks": [c.__dict__ for c in checks]})
        results[scene.scene_id] = r
        progress(f"[{scene.scene_id}] 完成，QC={overall(checks)}")

    status = "NOT_RUN"
    if "qc" in stages:
        status = write_qc_report(run_dir / "qc_report.md", cfg, run_id, scenes, results,
                                 all_checks, backend, notes)
    manifest["scenes"] = [scene_brief(s) for s in scenes]
    manifest["overall_status"] = status
    manifest["notes"] = notes
    core.save_json(run_dir / "manifest.json", manifest)
    core.save_json(run_dir / "run_results.json", results)
    return {"run_id": run_id, "run_dir": str(run_dir), "status": status,
            "scenes": [s.scene_id for s in scenes], "discovery_errors": errors,
            "notes": notes, "results": results}
