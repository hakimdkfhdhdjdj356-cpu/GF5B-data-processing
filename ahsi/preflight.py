"""换数据 / 换研究区前的可迁移性检查。

回答一个具体问题："把 input_root 换成新数据，能不能一条命令直接跑？"
不假设多景齐全、不假设传感器型号不变、不假设 DEM 覆盖新研究区。

本模块只读，不写原始输入、不改配置；仅输出判定与人工待办清单。
"""
from __future__ import annotations

import shutil
from pathlib import Path

import rasterio
from rasterio.warp import transform_bounds

from .config import aux_path

PASS, WARN, FAIL, SKIP, INFO = "PASS", "WARN", "FAIL", "SKIP", "INFO"


def _item(cid, status, observed, requirement, action="") -> dict:
    return {"id": cid, "status": status, "observed": str(observed),
            "requirement": str(requirement), "action": action}


def _utm_crs(lon: float, lat: float) -> str:
    zone = int((lon + 180) // 6) + 1
    return f"EPSG:{32600 + zone if lat >= 0 else 32700 + zone}"


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def migration_checklist(cfg: dict) -> list[dict]:
    """明确区分：必须人工准备 / 必须改配置 / 按需改 / 完全不用改。"""
    return [
        {"level": "必须人工准备", "item": "覆盖新研究区的 DEM",
         "detail": f"当前 auxiliary_root={cfg['auxiliary_root']} 下的 "
                   f"{cfg['grid']['rpc_dem']} 用于 RPC 正射；"
                   "新研究区若不在其范围内会直接失败。"},
        {"level": "必须人工准备", "item": "新产品目录结构",
         "detail": "每个产品目录内需齐全：1 个 XML、*_VN.tif / *_SW.tif、"
                   "同名 *.rpb、VNIR/SWIR 各一份定标 .raw 与光谱响应 .raw、OGP 角度表。"
                   "场景按 XML 逐个独立发现，单景即可运行。"},
        {"level": "必须改配置", "item": "input_root",
         "detail": "指向新数据根目录（可用 preflight --input-root 预演）"},
        {"level": "必须改配置", "item": "grid.rpc_dem / auxiliary_root",
         "detail": "换成覆盖新区 DEM 的文件名与所在目录；EGM2008 geoid 网格是全局的，不用换。"},
        {"level": "按需改配置", "item": "grid.crs",
         "detail": "auto_utm=true 时按每景中心经度自动推导，无需改；"
                   "但请核对 XML 的 CenterLongitude/CenterLatitude 是否可解析。"},
        {"level": "按需改配置", "item": "bands.* 与 sensors.yaml.file_patterns / known_sensors",
         "detail": "换传感器级别时波长范围、水汽吸收带、相机边缘、文件命名可能不同；"
                   "卫星名不在 known_sensors 也能跑（波段表从数据自读），只是少一层完整性校验。"},
        {"level": "按需改配置", "item": "golden_root",
         "detail": "只对 2026-11-06 那两景的已验证成果有效。换数据后 verify 会自动记 SKIP，"
                   "不会伪造通过，也不会阻断生产流程。"},
        {"level": "完全不用改", "item": "所有 .py 代码",
         "detail": "场景发现、定标、RPC 正射、IARR 归一化、QC 门禁、阶段编排、报告生成"
                   "全部由配置与输入数据驱动，无硬编码卫星名/路径/结论。"},
    ]


def run_preflight(cfg: dict, scenes, errors: list[dict], backend: dict) -> dict:
    items: list[dict] = []

    # 1. 输入与场景发现：单景合法，不以"两景齐全"为前提
    root = Path(cfg["input_root"])
    items.append(_item("input_root_exists", PASS if root.is_dir() else FAIL,
                       root, "目录存在且可读"))
    items.append(_item("scene_discovery", PASS if scenes and not errors else FAIL,
                       f"发现 {len(scenes)} 景，解析失败 {len(errors)}",
                       ">=1 景且全部可解析（单景合法，不要求多景/双传感器齐全）",
                       "解析失败时检查 XML 字段别名与文件配对规则" if errors else ""))
    for e in errors:
        items.append(_item("scene_parse_error", FAIL, f"{Path(e['xml']).name}: {e['error']}",
                           "可解析", "补全侧车文件或在 sensors.yaml 增加别名"))

    # 2. 每景：中心经纬度、UTM 带、传感器规格
    known = cfg["sensors"].get("known_sensors", {})
    auto_utm = bool(cfg["grid"]["auto_utm"])
    configured_crs = str(cfg["grid"]["crs"])
    for sc in scenes:
        sid = sc.scene_id
        lon, lat = _num(sc.meta.get("center_lon")), _num(sc.meta.get("center_lat"))
        items.append(_item(f"{sid}:center_coordinates", PASS if lon is not None and lat is not None else FAIL,
                           f"{sc.meta.get('center_lon')}, {sc.meta.get('center_lat')}",
                           "XML 可解析出 CenterLongitude / CenterLatitude",
                           "auto_utm 与 DEM 覆盖检查都依赖它；缺字段时在 sensors.yaml 补别名"))
        if lon is not None and lat is not None:
            expected = _utm_crs(lon, lat)
            if auto_utm:
                items.append(_item(f"{sid}:utm_zone", INFO,
                                   f"中心 ({lon:.3f}, {lat:.3f}) -> 自动使用 {expected}",
                                   "auto_utm=true，按景自动选带",
                                   f"配置里 grid.crs={configured_crs} 仅作回退/参考，不生效"))
            else:
                items.append(_item(f"{sid}:utm_zone", PASS if expected == configured_crs else FAIL,
                                   f"应有 {expected}，配置 {configured_crs}",
                                   "目标 CRS 与场景所在 UTM 带一致",
                                   "改 grid.crs，或设 auto_utm=true" if expected != configured_crs else ""))
        spec = known.get(sc.satellite)
        items.append(_item(f"{sid}:sensor_spec",
                           PASS if spec else WARN,
                           f"卫星名 {sc.satellite}，{sc.total_bands} 波段 "
                           f"(VNIR {sc.segments[0].count} / SWIR {sc.segments[1].count})",
                           "known_sensors 中有对应规格（缺失仍可运行）",
                           "建议在 sensors.yaml 补一条规格用于完整性校验" if not spec else ""))
        if sc.warnings:
            items.append(_item(f"{sid}:scene_warnings", WARN, "; ".join(sc.warnings),
                               "无告警", "逐条核对，不影响流程继续"))

    # 3. DEM / geoid：换研究区最常踩的坑
    dem = aux_path(cfg, "rpc_dem")
    items.append(_item("dem_file_exists", PASS if dem.is_file() else FAIL,
                       dem, "RPC 正射所需 DEM 存在",
                       "准备覆盖新研究区的 DEM 并放到 auxiliary_root"))
    if dem.is_file():
        try:
            with rasterio.open(dem) as d:
                b = transform_bounds(d.crs, "EPSG:4326", *d.bounds) if d.crs else d.bounds
                shape = [d.width, d.height]
                dem_crs = str(d.crs)
                dem_bounds = d.bounds
            items.append(_item("dem_readable", PASS,
                               f"{dem_crs} {shape[0]}x{shape[1]} bounds={tuple(round(x, 3) for x in b)}",
                               "可被 GDAL 读取且带 CRS"))
            for sc in scenes:
                lon, lat = _num(sc.meta.get("center_lon")), _num(sc.meta.get("center_lat"))
                if lon is None or lat is None:
                    continue
                requested = cfg["grid"].get("bounds_m")
                if requested is not None:
                    target = _utm_crs(lon, lat) if auto_utm else configured_crs
                    bounds = transform_bounds(dem_crs, target, *dem_bounds)
                    tolerance = 0.01
                    inside = (bounds[0] - tolerance <= requested[0]
                              and bounds[1] - tolerance <= requested[1]
                              and bounds[2] + tolerance >= requested[2]
                              and bounds[3] + tolerance >= requested[3])
                    items.append(_item(f"{sc.scene_id}:dem_covers_requested_bounds",
                                       PASS if inside else FAIL, f"{target} DEM {bounds} vs {requested}",
                                       "DEM 外包框覆盖配置范围",
                                       "像元级覆盖由正射输出 QA 核实" if inside else "扩展 DEM 或核对目标范围"))
                    continue
                inside = b[0] <= lon <= b[2] and b[1] <= lat <= b[3]
                items.append(_item(f"{sc.scene_id}:dem_covers_center", PASS if inside else FAIL,
                                   f"中心 ({lon:.3f}, {lat:.3f}) vs DEM {tuple(round(x, 3) for x in b)}",
                                   "场景中心落在 DEM 范围内",
                                   "中心点只是最低要求；换研究区请确认整景范围都被 DEM 覆盖"
                                   if inside else "必须换成覆盖新研究区的 DEM，否则 RPC 正射会失败或高程失真"))
        except Exception as exc:
            items.append(_item("dem_readable", FAIL, f"{type(exc).__name__}: {exc}",
                               "可被 GDAL 读取", "检查 DEM 格式与 CRS"))
    if cfg["grid"]["apply_vdatum_shift"]:
        geoid = aux_path(cfg, "geoid_grid")
        items.append(_item("geoid_grid_exists", PASS if geoid.is_file() else FAIL,
                           geoid, "EGM2008 高程网格存在（全局文件，换区不用换）"))

    # 4. 后端依赖
    method = backend["method"]
    if backend.get("requires_external_software") or backend.get("requires_external_binary"):
        executable = backend.get("sixs_executable") if method == "sixs_lut" else backend.get("engine")
        items.append(_item("external_backend", PASS if backend["available"] else FAIL,
                           f"{method}: {executable} 可用={backend['available']}",
                           ", ".join(backend["external_dependencies"]),
                           "检查已选后端的可执行文件及运行依赖" if not backend["available"] else ""))
    else:
        items.append(_item("external_backend", PASS,
                           f"{method}: 无外部软件依赖", "零外部依赖即可端到端产出"))

    # 5. golden 回归的适用性
    golden = Path(cfg["golden_root"]) if cfg.get("golden_root") else None
    if golden and golden.is_dir():
        matched = [sc.satellite for sc in scenes if (golden / sc.satellite / "processing_summary.json").is_file()]
        status = PASS if matched else SKIP
        items.append(_item("golden_applicability", status,
                           f"匹配 {matched or '无'}（golden 只对 2026-11-06 原两景有效）",
                           "新数据集无 golden 时记 SKIP，不算失败、不伪造通过",
                           "不要拿旧两景的回归结论给新数据背书"))
    else:
        items.append(_item("golden_applicability", SKIP, "未配置或不存在",
                           "换数据后 verify 自动跳过该景的数值比对"))

    # 6. 输出目录与磁盘
    out = Path(cfg["output_root"])
    try:
        out.mkdir(parents=True, exist_ok=True)
        free = shutil.disk_usage(str(out)).free / (1024 ** 3)
        need = float(cfg["runtime"]["min_free_gb"])
        items.append(_item("disk_space", PASS if free >= need else FAIL,
                           f"可用 {free:.1f} GB", f">= {need:g} GB",
                           "全量运行会产出数十 GB，SMOKE 可用 --max-bands 限波端"))
    except Exception as exc:
        items.append(_item("output_writable", FAIL, f"{type(exc).__name__}: {exc}", "输出根目录可写"))

    blocking = [i for i in items if i["status"] == FAIL]
    return {"ready": not blocking, "blocking": blocking,
            "scenes": len(scenes), "input_root": str(root), "method": method,
            "auto_utm": auto_utm, "configured_crs": configured_crs,
            "items": items, "checklist": migration_checklist(cfg),
            "counts": {s: sum(i["status"] == s for i in items) for s in (PASS, WARN, FAIL, SKIP, INFO)}}
