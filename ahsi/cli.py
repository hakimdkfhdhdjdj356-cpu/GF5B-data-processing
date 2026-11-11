"""单入口 CLI：plan / run / verify。"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import sys
from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from ahsi import atmos, core, pipeline, preflight  # noqa: E402
from ahsi.config import load_config  # noqa: E402
from ahsi.scene import discover_scenes, scene_brief  # noqa: E402

log = logging.getLogger("ahsi")


def _default_config() -> str:
    return str(PACKAGE_ROOT / "config" / "pipeline.yaml")


def _setup(verbose: bool) -> None:
    logging.basicConfig(level=logging.DEBUG if verbose else logging.INFO,
                        format="%(asctime)s %(levelname)-7s %(message)s", datefmt="%H:%M:%S")


def cmd_plan(args) -> int:
    cfg = load_config(args.config)
    scenes, errors = discover_scenes(cfg["input_root"], cfg)
    backend = atmos.detect(cfg)
    brief = [scene_brief(s) for s in scenes]
    print(f"输入根目录 : {cfg['input_root']}")
    print(f"发现场景   : {len(scenes)}    解析失败: {len(errors)}")
    crs_note = (f"auto_utm 按景推导（基准 {cfg['grid']['crs']}）"
                if cfg["grid"]["auto_utm"] else cfg["grid"]["crs"])
    print(f"目标网格   : {crs_note} @ {cfg['grid']['resolution_m']} m")
    print(f"大气后端   : {backend['method']} 可用={backend['available']} "
          f"外部依赖={backend['external_dependencies'] or '无'}")
    print("-" * 108)
    print(f"{'场景':<28} {'级别':<9} {'源波段':>6} {'候选':>5} {'VNIR/SWIR':>12} {'波长 nm':>18} {'中心经纬度':>20}")
    for s in brief:
        seg = s["segments"]
        pair = f"{seg['VNIR']['bands']}/{seg['SWIR']['bands']}"
        wl = f"{s['wavelength_range_nm'][0]:.1f}-{s['wavelength_range_nm'][1]:.1f}"
        ctr = f"{s['center_lon_lat'][0]}, {s['center_lon_lat'][1]}"
        print(f"{s['scene_id']:<28} {str(s['product_level']):<9} {s['source_bands']:>6} "
              f"{s['candidate_bands']:>5} {pair:>12} {wl:>18} {ctr:>20}")
    for s in brief:
        if s["warnings"]:
            print(f"  [WARN] {s['scene_id']}: {'; '.join(s['warnings'])}")
    for e in errors:
        print(f"  [ERROR] {e['xml']}: {e['error']}")
    print("-" * 108)
    print("提示：中心经纬度用于核对是否覆盖研究区；本流水线不替你判断地质适用性。")
    out = Path(args.out) if args.out else Path(cfg["output_root"]) / "plan.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    core.save_json(out, {"config": cfg["project_root"], "input_root": cfg["input_root"],
                         "grid": cfg["grid"], "backend": backend, "scenes": brief, "errors": errors})
    print(f"计划已写入: {out}")
    return 0 if scenes and not errors else 1


def cmd_run(args) -> int:
    cfg = load_config(args.config)
    stages = [s.strip() for s in args.stages.split(",")] if args.stages else None
    scenes = [s.strip() for s in args.scenes.split(",")] if args.scenes else None
    result = pipeline.run(cfg, run_id=args.run_id, scenes_filter=scenes, stages=stages,
                          max_bands=args.max_bands, force=args.force, progress=log.info)
    print(json.dumps({k: v for k, v in result.items() if k != "results"},
                     ensure_ascii=False, indent=2))
    return 0 if result["status"] in ("PASS", "WARN", "NOT_RUN") else 2


def _read_golden_csv(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def cmd_verify(args) -> int:
    """与已验证的人工成果做非破坏性回归：只比元数据与抽样像元，不重算 23 GB。"""
    import numpy as np
    import rasterio

    cfg = load_config(args.config)
    golden = Path(cfg["golden_root"])
    scenes, errors = discover_scenes(cfg["input_root"], cfg)
    items = []

    def add(scene, name, ok, expected, observed, note=""):
        items.append({"scene": scene, "check": name, "status": "PASS" if ok else "FAIL",
                      "expected": str(expected), "observed": str(observed), "note": note})

    if not golden.is_dir():
        print(f"提示：未找到 golden_root: {golden}；换数据后属正常，相关检查记 SKIP")

    for sc in scenes:
        sid, sat = sc.scene_id, sc.satellite
        gdir = golden / sat
        summary_path = gdir / "processing_summary.json"
        band_csv = gdir / "band_parameters_and_quality.csv"
        if not gdir.is_dir():
            # 换数据 / 换卫星名后没有 golden 参考是正常情况：记 SKIP，不伪造通过也不算失败
            items.append({"scene": sid, "check": "golden_applicability", "status": "SKIP",
                          "expected": "n/a", "observed": "该景无 golden 参考",
                          "note": "golden 只对 2026-11-06 原两景有效；不拿旧结论给新数据背书"})
            continue
        if not summary_path.is_file() or not band_csv.is_file():
            add(sid, "golden_artifacts_present", False, "summary + band csv", "missing",
                "golden 目录存在但文件不全，检查 golden_root 配置")
            continue
        summ = json.loads(summary_path.read_text(encoding="utf-8"))
        grows = _read_golden_csv(band_csv)

        add(sid, "source_band_count", sc.total_bands == int(summ["source_bands"]),
            summ["source_bands"], sc.total_bands)
        add(sid, "candidate_band_count", len(sc.candidate_bands) == int(summ["candidate_bands"]),
            summ["candidate_bands"], len(sc.candidate_bands))
        gw = [float(r["wavelength_nm"]) for r in grows]
        ow = [r["wavelength_nm"] for r in sc.band_table]
        add(sid, "wavelength_sequence", len(gw) == len(ow) and np.allclose(gw, ow, atol=1e-4),
            f"{len(gw)} bands {gw[0]:.4f}..{gw[-1]:.4f}",
            f"{len(ow)} bands {ow[0]:.4f}..{ow[-1]:.4f}")
        gg = [float(r["gain"]) for r in grows]
        og = [r["gain"] for r in sc.band_table]
        add(sid, "gain_sequence", len(gg) == len(og) and np.allclose(gg, og, atol=1e-6),
            "golden gain array", "discovered gain array")
        gc = [r["analysis_candidate"].strip().lower() == "true" for r in grows]
        oc = [r["analysis_candidate"] for r in sc.band_table]
        add(sid, "candidate_flags", gc == oc, f"{sum(gc)} True", f"{sum(oc)} True")
        add(sid, "grid_crs", summ["crs"] == cfg["grid"]["crs"], summ["crs"], cfg["grid"]["crs"])
        add(sid, "configured_resolution", float(summ["resolution_m"]) == float(cfg["grid"]["resolution_m"]),
            summ["resolution_m"], cfg["grid"]["resolution_m"])

        rad = gdir / f"{sat}_radiance_RPC_DEM_30m.tif"
        if rad.is_file():
            with rasterio.open(rad) as d:
                add(sid, "golden_grid_size",
                    [d.width, d.height] == [int(summ["width"]), int(summ["height"])],
                    f"{summ['width']}x{summ['height']}", f"{d.width}x{d.height}")
                env = d.tags(ns="ENVI")
                w = [float(x) for x in env.get("wavelength", "").strip("{}").split(",") if x]
                add(sid, "golden_radiance_band_count", d.count == len(ow),
                    len(ow), d.count)
                add(sid, "golden_radiance_wavelength_metadata",
                    len(w) == len(ow) and np.allclose(w, ow, atol=1e-4),
                    "ENVI wavelength tag == discovered", f"{len(w)} values")
                add(sid, "golden_radiance_crs", str(d.crs.to_epsg()) == str(cfg["grid"]["crs"]).split(":")[-1],
                    cfg["grid"]["crs"], d.crs)
                win = ((d.height // 2, d.height // 2 + 32), (d.width // 2, d.width // 2 + 32))
                for b in (1, d.count // 2, d.count):
                    a = d.read(b, window=win)
                    finite = np.isfinite(a).all()
                    add(sid, f"golden_radiance_sample_band_{b}", finite,
                        "finite values", f"min={float(a.min()):.5g} max={float(a.max()):.5g}")

        refs = sorted((gdir / "reflectance_final").glob("*_reflectance_float32_30m.tif"))
        if refs:
            with rasterio.open(refs[0]) as d:
                add(sid, "golden_reflectance_band_count", d.count == int(summ["candidate_bands"]),
                    summ["candidate_bands"], d.count)
                a = d.read(1, window=((d.height // 2, d.height // 2 + 32),
                                      (d.width // 2, d.width // 2 + 32)))
                v = a[a != d.nodata]
                ok = v.size > 0 and np.isfinite(v).all() and float(v.min()) >= 0 and float(v.max()) < 50
                add(sid, "golden_reflectance_scale_sanity", ok, "0 <= refl < 50 (already /10000)",
                    f"min={float(v.min()):.5g} max={float(v.max()):.5g}" if v.size else "no valid pixels")

        # 若本次新 run 存在，额外做数值级抽样比对
        if args.run_id:
            newcube = Path(cfg["output_root"]) / args.run_id / "scenes" / sid / f"{sid}_radiance.tif"
            if newcube.is_file() and rad.is_file():
                with rasterio.open(newcube) as n, rasterio.open(rad) as g:
                    win = ((g.height // 2, g.height // 2 + 16), (g.width // 2, g.width // 2 + 16))
                    for b in (1, min(n.count, g.count)):
                        a, c = n.read(b, window=win).astype("float64"), g.read(b, window=win).astype("float64")
                        both = (a != cfg["grid"]["nodata"]) & (c != -9999)
                        diff = float(np.nanmax(np.abs(a[both] - c[both]))) if both.any() else None
                        add(sid, f"new_vs_golden_radiance_band_{b}",
                            diff is not None and diff < 1e-4, "max|diff| < 1e-4", diff,
                            "同数值内核的回归证据")
            else:
                # 该 run 本来就没处理这一景（--scenes 过滤）：跳过，不算失败
                scene_run = Path(cfg["output_root"]) / args.run_id / "scenes" / sid
                if (scene_run / "stage_radiance.json").is_file():
                    add(sid, "new_run_available", False, "radiance.tif", "missing", "该景标记已运行但找不到产物")
                else:
                    items.append({"scene": sid, "check": "new_vs_golden_sampling", "status": "SKIP",
                                  "expected": "n/a", "observed": "scene not in this run",
                                  "note": "--scenes 过滤，非失败"})

    failed = [i for i in items if i["status"] == "FAIL"]
    out = Path(args.out) if args.out else Path(cfg["output_root"]) / "verify_report.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    lines = ["# 与已验证成果的回归比对", "",
             f"- golden_root: `{golden}`", f"- 检查项 {len(items)}，失败 {len(failed)}",
             f"- 总判定：**{'PASS' if not failed else 'FAIL'}**", "",
             "| 场景 | 检查项 | 结果 | 期望 | 实测 | 说明 |", "| --- | --- | --- | --- | --- | --- |"]
    for i in items:
        lines.append(f"| {i['scene']} | `{i['check']}` | {i['status']} | {i['expected']} | "
                     f"{i['observed']} | {i['note']} |")
    lines += ["", "说明：本比对只读取已保存的元数据与小窗口抽样像元，不重算大立方体；"
              "通过表示新流水线的场景发现、波段表与筛选规则同已验证成果一致，"
              "不代表绝对定位精度或反射率物理精度通过验收。",
              "`SKIP` 表示该景没有对应的 golden 参考（换数据/换卫星名后属正常），"
              "既不伪造通过，也不阻断生产流程；新数据的正确性请用 `preflight` 与 `run` 的 QC 报告判定。", ""]
    out.write_text("\n".join(lines), encoding="utf-8")
    core.save_json(out.with_suffix(".json"), {"items": items, "discovery_errors": errors,
                                              "failed": len(failed)})
    print(f"回归比对: {len(items)} 项，失败 {len(failed)}；报告 {out}")
    for f in failed:
        print(f"  FAIL {f['scene']} {f['check']}: expected {f['expected']}, got {f['observed']}")
    return 0 if not failed else 2


def cmd_preflight(args) -> int:
    """换数据/换研究区前，把"哪些能自动、哪些必须人工准备"一次列清楚。"""
    cfg = load_config(args.config)
    preview = bool(args.input_root)
    if preview:
        cfg["input_root"] = str(Path(args.input_root).resolve())
    try:
        scenes, errors = discover_scenes(cfg["input_root"], cfg)
    except FileNotFoundError as exc:
        print(f"[FAIL] {exc}")
        return 2
    backend = atmos.detect(cfg)
    result = preflight.run_preflight(cfg, scenes, errors, backend)
    result["preview_input_root_only"] = preview

    print("=== 换数据 / 换研究区可迁移性检查 ===")
    if preview:
        print(f"[预演] 临时 input_root={cfg['input_root']}（配置文件未修改）")
    print(f"输入根目录 : {cfg['input_root']}")
    print(f"发现场景   : {result['scenes']}（单景合法，不要求多景 / 双传感器齐全）")
    print(f"大气后端   : {result['method']}   目标投影: "
          f"{'auto_utm 按景自动推导' if result['auto_utm'] else result['configured_crs']}")
    print("-" * 100)
    for i in result["items"]:
        print(f"[{i['status']:<4}] {i['id']:<46} {i['observed']}")
        if i["status"] == "FAIL" and i["action"]:
            print(f"         -> {i['action']}")
    print("-" * 100)
    print(f"检查 {len(result['items'])} 项：" +
          "  ".join(f"{k}={v}" for k, v in result["counts"].items()))
    print(f"判定：{'READY，可直接 run' if result['ready'] else 'NOT READY，请先处理阻塞项'}")
    print("\n【换数据清单】")
    for c in result["checklist"]:
        print(f"  [{c['level']}] {c['item']}\n      {c['detail']}")

    out = Path(args.out) if args.out else Path(cfg["output_root"]) / "preflight.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    core.save_json(out, result)
    lines = ["# 换数据 / 换研究区可迁移性检查", "",
             f"- 判定：**{'READY' if result['ready'] else 'NOT READY'}**",
             f"- 输入根目录：`{cfg['input_root']}`",
             f"- 发现场景：{result['scenes']}（单景合法，不要求多景 / 双传感器齐全）",
             f"- 大气后端：`{result['method']}`，目标投影："
             f"{'auto_utm 按景自动推导' if result['auto_utm'] else result['configured_crs']}",
             f"- 检查 {len(result['items'])} 项：" +
             "，".join(f"{k}={v}" for k, v in result["counts"].items()), "",
             "## 一、检查结果", "",
             "| 检查项 | 结果 | 实测 | 要求 | 处理建议 |",
             "| --- | --- | --- | --- | --- |"]
    for i in result["items"]:
        lines.append(f"| `{i['id']}` | {i['status']} | {i['observed']} | "
                     f"{i['requirement']} | {i['action']} |")
    lines += ["", "## 二、换数据清单", "",
              "| 类别 | 项目 | 说明 |", "| --- | --- | --- |"]
    for c in result["checklist"]:
        lines.append(f"| {c['level']} | {c['item']} | {c['detail']} |")
    lines += ["", "## 三、边界", "",
              "1. 本检查只读输入与辅助数据，不修改配置、不写原始数据、不产出影像。",
              "2. DEM 只核对场景中心点是否落在范围内；**整景覆盖仍需人工确认**。",
              "3. `golden_root` 的回归结论只对 2026-11-06 原两景有效，换数据后自动记 SKIP。",
              "4. 判定为 READY 只表示流水线可以开跑，不代表绝对定位精度或反射率物理精度通过验收。", ""]
    md = out.with_suffix(".md")
    md.write_text("\n".join(lines), encoding="utf-8")
    print(f"\n报告：{out} / {md}")
    return 0 if result["ready"] else 2


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="run_ahsi", description="星载高光谱 L1 自动化预处理流水线")
    p.add_argument("--config", default=_default_config(), help="主配置 yaml 路径")
    p.add_argument("-v", "--verbose", action="store_true")
    sub = p.add_subparsers(dest="command", required=True)

    pl = sub.add_parser("plan", help="扫描输入、校验配对、输出场景计划（秒级，不写大文件）")
    pl.add_argument("--out")
    pl.set_defaults(func=cmd_plan)

    rn = sub.add_parser("run", help="执行完整流水线")
    rn.add_argument("--run-id", help="指定运行 ID；复用即断点续跑")
    rn.add_argument("--scenes", help="逗号分隔，按卫星名或场景 ID 过滤")
    rn.add_argument("--stages", help=f"逗号分隔，可选 {','.join(pipeline.STAGES)}")
    rn.add_argument("--max-bands", type=int, help="仅输出前 N 个波段（端到端 SMOKE，非交付）")
    rn.add_argument("--force", action="store_true", help="忽略已完成阶段标记，强制重算")
    rn.set_defaults(func=cmd_run)

    pf = sub.add_parser("preflight", help="换数据/换研究区前的可迁移性检查与人工待办清单")
    pf.add_argument("--input-root", help="预演：临时覆盖 input_root，不修改配置文件")
    pf.add_argument("--out", help="输出 json 路径（同时写同名 .md）")
    pf.set_defaults(func=cmd_preflight)

    vf = sub.add_parser("verify", help="与已验证成果做非破坏性回归比对（无 golden 时记 SKIP）")
    vf.add_argument("--run-id", help="可选：同时比对新 run 的辐亮度抽样像元")
    vf.add_argument("--out")
    vf.set_defaults(func=cmd_verify)

    args = p.parse_args(argv)
    _setup(args.verbose)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
