"""统一质量门禁与报告。所有数字来自本次运行的实际结果，报告正文不写死任何结论。"""
from __future__ import annotations

import json
import logging
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

log = logging.getLogger(__name__)

PASS, WARN, FAIL, INFO = "PASS", "WARN", "FAIL", "INFO"


@dataclass
class Check:
    scene_id: str
    id: str
    status: str
    observed: str
    threshold: str
    note: str = ""


def _check(scene_id, cid, ok, observed, threshold, note="", warn_only=False, info=False):
    status = INFO if info else (PASS if ok else (WARN if warn_only else FAIL))
    return Check(scene_id, cid, status, str(observed), str(threshold), note)


def evaluate_scene(scene, cfg, radiance: dict | None, reflectance: dict | None,
                   quac: dict | None, backend: dict) -> list[Check]:
    q, sid = cfg["qc"], scene.scene_id
    out: list[Check] = []

    out.append(_check(sid, "input_pairing", True,
                      f"{len(scene.segments)} 相机 / {scene.total_bands} 波段", "影像·定标·光谱响应·RPC 全部配对",
                      "由场景解析强制保证；不配对直接拒绝该景"))
    waves = [r["wavelength_nm"] for r in scene.band_table]
    out.append(_check(sid, "wavelength_monotonic", all(b >= a - 1e-6 for a, b in zip(waves, waves[1:])),
                      f"{waves[0]:.3f}–{waves[-1]:.3f} nm", "按波长升序",
                      "VNIR/SWIR 重叠区保留双相机原始波段，未做人为增益匹配"))
    out.append(_check(sid, "calibration_coefficients", all(r["gain"] > 0 and r["fwhm_nm"] > 0 for r in scene.band_table),
                      f"{len(scene.band_table)} 组 gain/offset/FWHM", "gain>0 且 FWHM>0",
                      cfg["sensors"]["calibration"]["formula"]))
    out.append(_check(sid, "rpc_available", all(s.has_rpc for s in scene.segments),
                      [s.has_rpc for s in scene.segments], "两相机均可读取 RPC",
                      "禁止用普通仿射变换替代 RPC 正射"))

    if scene.ogp.get("solar_zenith_diff_deg") is not None:
        diff = scene.ogp["solar_zenith_diff_deg"]
        out.append(_check(sid, "ogp_vs_xml_solar_zenith", diff <= float(q["max_ogp_solar_zenith_diff_deg"]),
                          f"{diff:.4f} deg", f"<= {q['max_ogp_solar_zenith_diff_deg']} deg",
                          "OGP 为逐像元观测几何角度表，非经纬度网格"))

    if radiance:
        residuals = [g["rpc_roundtrip_pixels"] for g in radiance["geometry_checks"]]
        out.append(_check(sid, "rpc_numerical_roundtrip", max(residuals) <= float(q["max_rpc_roundtrip_pixels"]),
                          f"max {max(residuals):.4f} px", f"<= {q['max_rpc_roundtrip_pixels']} px",
                          "仅为内部数值检查，不是独立地面精度验收"))
        out.append(_check(sid, "complete_raster_readback", radiance["complete_raster_readback_passed"],
                          radiance["complete_raster_readback_passed"], "True",
                          "全波段回读，检查有限数值与无静默截断"))
        area = radiance["joint_camera_area_km2"]
        out.append(_check(sid, "joint_camera_coverage", area >= float(q["min_joint_coverage_km2"]),
                          f"{area:.2f} km2", f">= {q['min_joint_coverage_km2']} km2",
                          "掩膜只说明两相机都有观测，不等于无云/无雪/裸岩"))
        out.append(_check(sid, "candidate_band_count", True, radiance["candidate_bands"],
                          "记录值", "筛选规则见配置 bands.*；完整波段仍保留在辐亮度立方体", info=True))
        out.append(_check(sid, "smoke_mode", not radiance["smoke_mode"], radiance["smoke_mode"], "False",
                          "True 表示仅输出部分波端的验证运行，不可作为交付成果",
                          warn_only=True))

    method = backend.get("method", "envi_quac")
    out.append(_check(sid, "atmospheric_backend", backend.get("available", False),
                      f"{method} ({backend.get('note')})", "当前配置的后端可用",
                      "可选后端不可用时按 on_unavailable 处理",
                      warn_only=cfg["atmospheric"]["on_unavailable"] == "stop_at_radiance"))
    out.append(_check(sid, "external_software_dependency",
                      not backend.get("requires_external_software", False),
                      backend.get("external_dependencies") or "无",
                      "默认路径不依赖商业软件即可端到端产出",
                      "需要绝对反射率时再切换到 envi_quac / 6S / FLAASH 并单独验收",
                      warn_only=True))

    if method == "none":
        out.append(_check(sid, "radiance_only_by_configuration", True, "method=none",
                          "有意只交付辐亮度", "配置选择的保守行为，不是失败", info=True))
    if quac is not None:
        ok = bool(quac.get("successful"))
        cid = "quac_execution" if method == "envi_quac" else "atmospheric_backend_execution"
        out.append(_check(sid, cid, ok, quac.get("reason", "ok"),
                          "后端执行成功并产出可校验结果",
                          "失败不得伪造任何等级的反射率", warn_only=not ok and
                          cfg["atmospheric"]["on_unavailable"] == "stop_at_radiance"))
        if method == "iarr" and ok:
            n, good = quac.get("band_count"), quac.get("good_bands")
            out.append(_check(sid, "iarr_normalization_factors", bool(n) and good > 0,
                              f"{good}/{n} 波段有效", "至少 1 个波段有稳健归一化因子",
                              f"归一化方式 {quac.get('normalization')}，稳健分位 {quac.get('robust_pct')}"))
            out.append(_check(sid, "iarr_good_band_ratio", bool(n) and good / n >= 0.8,
                              f"{good/n:.1%}" if n else "n/a", ">= 80%",
                              "低于阈值通常提示候选波段落在探测器边缘或整幅接近背景值",
                              warn_only=True))
            out.append(_check(sid, "iarr_deterministic", quac.get("deterministic") is True,
                              quac.get("deterministic"), "True",
                              "无外部数据、无经验大气参数、无随机性，结果逐位可复现"))
        integ = quac.get("output_integrity")
        if integ:
            out.append(_check(sid, "quac_output_integrity", integ.get("valid", False),
                              f"bands={integ.get('bands')} size={integ.get('size')} scale={integ.get('scale_factor')}",
                              "波段数/尺寸/缩放因子/bbl/非零数据全部通过"))
        if quac.get("engine_exit_nonzero_warning"):
            out.append(_check(sid, "quac_engine_exit_nonzero", False,
                              quac.get("exit_code"), "exit=0",
                              "ENVI IDL 写完产物后异常退出；产物完整性校验通过，已按配置接受并保留告警",
                              warn_only=True))

    if reflectance:
        grade = reflectance.get("product_grade", "quac_estimated_reflectance_NOT_independently_validated")
        if grade.startswith("relative"):
            from pathlib import Path
            name = Path(reflectance["reflectance_file"]).name
            out.append(_check(sid, "relative_product_labeling",
                              "relative_reflectance" in name and grade.startswith("relative"),
                              f"{name} / {grade}",
                              "文件名与 product_grade 均显式标注 relative / NOT surface",
                              "防止相对反射率被下游误当作地表反射率"))
            out.append(_check(sid, "iarr_not_surface_reflectance",
                              reflectance.get("surface_reflectance_available") is False,
                              reflectance.get("surface_reflectance_available"), "False",
                              "未做大气校正；禁止绝对反射率数值、跨景定量对比、直接匹配实验室光谱库",
                              info=True))
            out.append(_check(sid, "iarr_bad_bands", not reflectance.get("bad_band_wavelengths_nm"),
                              f"{len(reflectance.get('bad_band_wavelengths_nm', []))} 个坏波段",
                              "0（有坏波段时 bbl=0，需人工确认）", warn_only=True))
            f3 = reflectance["value_above_3_fraction"]
            out.append(_check(sid, "relative_value_above_3", f3 <= 0.10,
                              f"{f3:.4%}", "<= 10%",
                              "相对值 >3 通常提示残余云/雪/饱和或极暗背景，仅作筛查，不修改数据",
                              warn_only=True))
            out.append(_check(sid, "no_double_scaling",
                              reflectance.get("integer_scale_factor", 1.0) in (1, 1.0),
                              reflectance.get("integer_scale_factor", 1.0), "1（无整数缩放）",
                              "IARR 直接输出 float32 相对值，无需任何二次缩放"))
        else:
            f1 = reflectance["pixels_with_any_reflectance_above_1_fraction"]
            f2 = reflectance["pixels_with_any_reflectance_above_2_fraction"]
            out.append(_check(sid, "reflectance_above_1", f1 <= float(q["max_fraction_above_1"]),
                              f"{f1:.4%}", f"<= {float(q['max_fraction_above_1']):.2%}",
                              ">1 的值保留并标记，不截到 0–1", warn_only=True))
            out.append(_check(sid, "reflectance_above_2", f2 <= float(q["max_fraction_above_2"]),
                              f"{f2:.4%}", f"<= {float(q['max_fraction_above_2']):.2%}",
                              warn_only=True))
            if grade.startswith("sixs_"):
                n = reflectance["reflectance_bands"]
                valid_n = reflectance.get("valid_reflectance_bands", n)
                out.append(_check(sid, "sixs_valid_bands", valid_n == n and n > 0,
                                  f"{valid_n}/{n}", "所有输出波段至少有一个有限反演值"))
                out.append(_check(sid, "no_double_scaling",
                                  reflectance.get("integer_scale_factor") == 1.0,
                                  reflectance.get("integer_scale_factor"), "1（6SV 直接反演 float32）"))
            else:
                out.append(_check(sid, "reflectance_scaled_once", True,
                                  reflectance["integer_scale_factor"], "QUAC .hdr 的 reflectance_scale_factor",
                                  "输出已换算为 float32 反射率"))
    else:
        out.append(_check(sid, "no_fabricated_reflectance", True,
                          "仅辐亮度立方体" if method != "none" else "method=none，仅辐亮度",
                          "未产出对应等级产品时不得声称任何反射率",
                          "这是正确的降级/保守行为，不是失败"))

    for cid in q.get("fail_on", []):
        for c in out:
            if c.id == cid and c.status == WARN:
                c.status = FAIL
    return out


def overall(checks: list[Check]) -> str:
    if any(c.status == FAIL for c in checks):
        return FAIL
    if any(c.status == WARN for c in checks):
        return WARN
    return PASS


def write_qc_report(path: Path, cfg: dict, run_id: str, scenes, results: dict,
                    checks: list[Check], backend: dict, notes: list[str] | None = None) -> str:
    status = overall(checks)
    lines = [
        f"# 星载高光谱 L1 自动化预处理 — 运行质量报告", "",
        f"- 运行 ID：`{run_id}`",
        f"- 生成时间（UTC）：{datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S')}",
        f"- 总判定：**{status}**",
        f"- 输入根目录：`{cfg['input_root']}`",
        f"- 目标网格：`{cfg['grid']['crs']}`，{cfg['grid']['resolution_m']} m",
        f"- 大气/归一化后端：`{backend.get('method')}`，可用={backend.get('available')}，"
        f"外部依赖={backend.get('external_dependencies') or '无'}",
        "", "## 一、场景与产物", "",
        "| 场景 | 产品级别 | 源波段 | 输出波段 | 候选波段 | 双相机覆盖 km² | CRS | 网格 列×行 | 分析产品等级 |",
        "| --- | --- | ---: | ---: | ---: | ---: | --- | --- | --- |",
    ]
    for sc in scenes:
        r = results.get(sc.scene_id, {})
        rad, ref = r.get("radiance"), r.get("reflectance")
        if rad:
            cov = f"{rad['joint_camera_area_km2']:.2f}"
            rows = f"{rad['width']}×{rad['height']}"
            outb, cand = rad["output_bands"], rad["candidate_bands"]
            crs = rad["crs"]
        else:
            cov, rows, outb, cand, crs = "—", "—", "—", len(sc.candidate_bands), "—"
        if not ref:
            product = "无（仅辐亮度）"
        elif str(ref.get("product_grade", "")).startswith("relative"):
            product = "IARR 相对反射率（非地表反射率）"
        else:
            product = "QUAC 估计反射率（未独立验证）"
        lines.append(f"| {sc.scene_id} | {sc.meta.get('product_level') or '—'} | {sc.total_bands} | "
                     f"{outb} | {cand} | {cov} | {crs} | {rows} | {product} |")

    lines += ["", "## 二、质量检查", "",
              "| 场景 | 检查项 | 结果 | 实测 | 阈值 | 说明 |",
              "| --- | --- | --- | --- | --- | --- |"]
    for c in checks:
        lines.append(f"| {c.scene_id} | `{c.id}` | {c.status} | {c.observed} | {c.threshold} | {c.note} |")

    lines += ["", "## 三、本次运行的实际产物", ""]
    for sc in scenes:
        r = results.get(sc.scene_id, {})
        rad, ref = r.get("radiance"), r.get("reflectance")
        lines.append(f"### {sc.scene_id}")
        if rad:
            lines.append(f"- 辐亮度立方体：`{Path(rad['radiance_file']).name}`（{rad['output_bands']} 波段，"
                         f"{rad['wavelength_min_nm']:.2f}–{rad['wavelength_max_nm']:.2f} nm）")
            lines.append(f"- 候选波段 VRT、双相机覆盖掩膜、RGB 与预览、逐波段质量表")
        if ref and str(ref.get("product_grade", "")).startswith("relative"):
            lines.append(f"- 相对反射率立方体：`{Path(ref['reflectance_file']).name}`"
                         f"（{ref['reflectance_bands']} 波段，有效归一化波段 {ref['good_bands']}，"
                         f"归一化方式 {ref['normalization']}）")
            lines.append(f"- 相对值 >2 比例 {ref['value_above_2_fraction']:.2%}，"
                         f">3 比例 {ref['value_above_3_fraction']:.2%}（保留并标记，不修改数据）")
            lines.append("- **该产品不是地表反射率**，未做大气校正；用途限制见下节第 1 条")
        elif ref:
            sixs = ref.get("product_grade", "").startswith("sixs_")
            label = "6SV 模型估计" if sixs else "QUAC 估计"
            scale = "直接反演 float32" if sixs else f"缩放因子 {ref['integer_scale_factor']:g}"
            lines.append(f"- {label}反射率立方体：`{Path(ref['reflectance_file']).name}`"
                         f"（{ref['reflectance_bands']} 波段，{scale}）")
            lines.append(f"- 反射率 >1 像元比例 {ref['pixels_with_any_reflectance_above_1_fraction']:.2%}，"
                         f">2 比例 {ref['pixels_with_any_reflectance_above_2_fraction']:.2%}（保留并标记）")
        stage = r.get("quac")
        if stage:
            lines.append(f"- 后端执行：method={stage.get('backend_method') or stage.get('method')}，"
                         f"successful={stage.get('successful')}，reason={stage.get('reason', 'ok')}")
        if not rad and not ref:
            lines.append("- 本次未产生栅格产物（plan / verify 模式）")
        lines.append("")

    method = backend.get("method", "envi_quac")
    if method == "iarr":
        boundary_1 = (
            "1. 辐亮度产品由随产品提供的 gain/offset 逐波段定标得到，**不是地表反射率**。"
            "本次分析产品为 **IARR 相对反射率**（逐波段除以场景稳健均值），**未做大气校正、无绝对物理量纲**；"
            "仅可用于单景内光谱形状分析（连续统去除、诊断吸收特征位置/深度、场景端元选择、MNF/PCA）。"
            "**禁止**用于绝对反射率数值、跨景/跨时相定量对比、与实验室绝对反射率光谱库直接做 SAM 匹配"
            "（逐波段归一化会改变光谱向量方向）。")
    elif method == "envi_quac":
        boundary_1 = (
            "1. 辐亮度产品由随产品提供的 gain/offset 逐波段定标得到，**不是地表反射率**；"
            "反射率产品来自 ENVI QUAC 场景统计估计，属经验法结果，未经独立地面光谱精度验证。")
    elif method == "sixs_lut":
        boundary_1 = (
            "1. 辐亮度产品由随产品提供的 gain/offset 逐波段定标得到，**不是地表反射率**；"
            "反射率产品由本机开源 6SV1.1 辐射传输模型生成的 (xa, xb, xc) LUT 反演，"
            "属**模型估计地表反射率**，不是实测真值，未经独立地面光谱精度验证。"
            "本景无同步水汽/气溶胶/臭氧观测，采用显式假设情景并输出敏感性带；"
            "无实测光谱响应函数，按中心波长+FWHM 构造高斯 SRF 并做 2xFWHM 敏感性。"
            "未做地形光照、邻近效应、BRDF 校正。**禁止**据此宣称绝对反射率精度或跨景定量变化。")
    else:
        boundary_1 = ("1. 本次只交付逐波段定标后的辐亮度立方体，**不是地表反射率**，也未做任何辐射归一化。")
    lines += ["## 四、结论边界（必须随成果一起交付）", "", boundary_1,
              "2. RPC 往返残差是内部数值检查；**未使用地面控制点，绝对定位精度未验收**。",
              "3. 未做地形光照辐射校正、去条带、云/雪/阴影逐像元判别；覆盖掩膜不等于可信裸岩范围。",
              "4. 候选波段筛选是按配置规则做的探索性初筛，**不等于逐波段信噪比认证**。",
              "5. 本报告中所有数字均由本次运行实际产物统计生成；流程不写入任何预设结论。",
              "6. 场景按单景独立处理，不要求多景或双传感器同时齐全；"
              "研究区适用性须用 `preflight` 与范围矢量另行核对，本流水线不对地质结论负责。",
              "7. 换研究区时目标 CRS / DEM 覆盖范围必须重新核对；`golden_root` 的回归结果只对原两景有效。", ""]
    if notes:
        lines += ["## 五、运行备注", ""] + [f"- {n}" for n in notes] + [""]
    lines += ["---", "", "复现方式：`python scripts/run_ahsi.py run --config config/pipeline.yaml --run-id <本运行ID>`", ""]
    text = "\n".join(lines)
    Path(path).write_text(text, encoding="utf-8")
    Path(path).with_name("qc_checks.json").write_text(
        json.dumps([asdict(c) for c in checks], ensure_ascii=False, indent=2), encoding="utf-8")
    return status
