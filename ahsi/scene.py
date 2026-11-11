"""场景发现、产品元数据解析、波段表构建。

支持 GF5B AHSI 与 ZY1E AHSI 两类 L1 产品；差异全部由 config/sensors.yaml 描述。
只读原始目录，从不修改、解压或移动输入文件。
"""
from __future__ import annotations

import io
import logging
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import rasterio

log = logging.getLogger(__name__)


@dataclass
class Segment:
    name: str                  # VNIR / SWIR
    image: Path
    rpc_file: Path
    radcal_file: Path
    spectral_file: Path
    cal: np.ndarray            # (n,2) gain, offset
    spec: np.ndarray           # (n,2) center_wavelength_nm, fwhm_nm
    count: int
    width: int
    height: int
    dtype: str
    has_rpc: bool


@dataclass
class Scene:
    scene_id: str
    satellite: str
    directory: Path
    meta: dict
    files: dict
    segments: list[Segment]
    band_table: list[dict] = field(default_factory=list)
    ogp: dict = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)

    @property
    def total_bands(self) -> int:
        return sum(s.count for s in self.segments)

    @property
    def candidate_bands(self) -> list[dict]:
        return [r for r in self.band_table if r["analysis_candidate"]]


def numeric_table(path: Path) -> np.ndarray:
    """定标/光谱响应表：GF5B 逗号分隔，ZY1E 空格分隔，统一兼容。"""
    text = Path(path).read_text("utf-8-sig").replace(",", " ")
    data = np.loadtxt(io.StringIO(text), ndmin=2)
    if data.shape[1] != 2 or not np.isfinite(data).all():
        raise ValueError(f"辅助参数表格式错误（应为两列数值）: {path}")
    return data


def _flatten_xml(path: Path) -> dict:
    root = ET.parse(path).getroot()
    return {n.tag: (n.text or "").strip() for n in root.iter() if not list(n)}


def _alias(meta: dict, aliases: dict, key: str, default=None):
    for name in aliases.get(key, []):
        if name in meta and meta[name] != "":
            return meta[name]
    return default


def _one(directory: Path, patterns: list[str], label: str) -> Path:
    hits = []
    for pattern in patterns:
        hits.extend(p for p in directory.glob(pattern) if p.is_file())
    hits = sorted(set(hits))
    if len(hits) != 1:
        raise ValueError(f"{label}: 期望 1 个文件，实际 {len(hits)} 个 -> {[p.name for p in hits][:5]}")
    return hits[0]


def _by_keywords(directory: Path, keywords: list[str], label: str) -> Path:
    keys = [k.lower() for k in keywords]
    hits = [p for p in directory.glob("*.raw")
            if p.is_file() and all(k in p.name.lower() for k in keys)]
    hits = sorted(set(hits))
    if len(hits) != 1:
        raise ValueError(f"{label}: 期望 1 个 .raw，实际 {len(hits)} 个 -> {[p.name for p in hits][:5]}")
    return hits[0]


def _read_ogp(path: Path, scale: float, order: list[str]) -> dict:
    with rasterio.open(path) as ds:
        if ds.count != 4:
            raise ValueError(f"OGP 应为 4 波段角度表，实际 {ds.count}: {path}")
        center = ds.read([1, 2, 3, 4], window=((ds.height // 2, ds.height // 2 + 1),
                                               (ds.width // 2, ds.width // 2 + 1)))
        vals = [float(v[0, 0]) / scale for v in center]
        sample = ds.read(1, out_shape=(64, 64)).astype("float64") / scale
        return {
            "file": path.name,
            "scale": scale,
            "center": dict(zip(order, vals)),
            "solar_zenith_range_deg": [float(sample.min()), float(sample.max())],
            "shape": [ds.height, ds.width],
            "interpretation": "逐像元观测几何角度查找表；非地理编码网格、非经纬度",
        }


def _screening_reasons(wave: float, segment: str, cfg: dict) -> list[str]:
    b = cfg["bands"]
    reasons = []
    if wave < b["wavelength_min_nm"] or wave > b["wavelength_max_nm"]:
        reasons.append("sensor_edge_candidate_exclusion")
    for low, high in b["water_absorption_nm"]:
        if low <= wave <= high:
            reasons.append("atmospheric_absorption_screening")
            break
    if segment == "VNIR" and wave > b["vnir_max_nm"]:
        reasons.append("VNIR_SWIR_edge_overlap_screening")
    if segment == "SWIR" and wave < b["swir_min_nm"]:
        reasons.append("VNIR_SWIR_edge_overlap_screening")
    return reasons


def build_band_table(scene: "Scene", cfg: dict) -> list[dict]:
    records = []
    for seg in scene.segments:
        for i, (wave, fwhm) in enumerate(seg.spec):
            gain, offset = seg.cal[i]
            reasons = _screening_reasons(float(wave), seg.name, cfg)
            records.append({
                "segment": seg.name,
                "source_band": i + 1,
                "source_file": seg.image.name,
                "wavelength_nm": float(wave),
                "fwhm_nm": float(fwhm),
                "gain": float(gain),
                "offset": float(offset),
                "analysis_candidate": not reasons,
                "screening_reason": ";".join(reasons),
            })
    records.sort(key=lambda r: (r["wavelength_nm"], r["segment"]))
    for i, r in enumerate(records, 1):
        r["output_band"] = i

    waves = np.array([r["wavelength_nm"] for r in records])
    if cfg["qc"]["require_wavelength_monotonic"] and len(waves) > 1:
        # 两相机重叠区允许相等/极近波长，但不能整体倒序
        if np.any(np.diff(waves) < -1e-6):
            raise ValueError("波长排序后出现倒序，检查 SpectralResponse 文件")
    if not all(r["gain"] > 0 for r in records):
        raise ValueError("存在非正定标增益")
    if not all(r["fwhm_nm"] > 0 for r in records):
        raise ValueError("存在非正 FWHM")
    return records


def parse_scene(xml_path: Path, cfg: dict) -> Scene:
    patterns = cfg["sensors"]["file_patterns"]
    aliases = cfg["sensors"]["meta_aliases"]
    directory = xml_path.parent
    raw_meta = _flatten_xml(xml_path)

    meta = {key: _alias(raw_meta, aliases, key) for key in aliases}
    satellite = str(meta.get("satellite") or "UNKNOWN")
    product_id = str(meta.get("product_id") or directory.name)
    date = str(meta.get("start_time") or "")[:10].replace("-", "") or "unknown"
    scene_id = f"{satellite}_{product_id}_{date}"

    files = {"metadata": xml_path}
    files["vn_image"] = _one(directory, patterns["vn_image"], "VNIR 影像")
    files["sw_image"] = _one(directory, patterns["sw_image"], "SWIR 影像")
    files["vn_rpc"] = _one(directory, patterns["vn_rpc"], "VNIR RPC")
    files["sw_rpc"] = _one(directory, patterns["sw_rpc"], "SWIR RPC")
    files["radcal_vnir"] = _by_keywords(directory, patterns["radcal_vnir_keywords"], "VNIR 定标")
    files["radcal_swir"] = _by_keywords(directory, patterns["radcal_swir_keywords"], "SWIR 定标")
    files["spectral_vnir"] = _by_keywords(directory, patterns["spectral_vnir_keywords"], "VNIR 光谱响应")
    files["spectral_swir"] = _by_keywords(directory, patterns["spectral_swir_keywords"], "SWIR 光谱响应")
    ogp_hits = sorted(set(p for pat in patterns["ogp"] for p in directory.glob(pat) if p.is_file()))
    files["ogp"] = ogp_hits[0] if ogp_hits else None
    shp_hits = sorted(set(p for pat in patterns["footprint"] for p in directory.glob(pat) if p.is_file()))
    files["footprint"] = shp_hits[0] if shp_hits else None

    segments = []
    for seg_name, img_key, cal_key, spec_key, rpc_key in [
        ("VNIR", "vn_image", "radcal_vnir", "spectral_vnir", "vn_rpc"),
        ("SWIR", "sw_image", "radcal_swir", "spectral_swir", "sw_rpc"),
    ]:
        cal = numeric_table(files[cal_key])
        spec = numeric_table(files[spec_key])
        with rasterio.open(files[img_key]) as ds:
            count, width, height, dtype = ds.count, ds.width, ds.height, ds.dtypes[0]
            has_rpc = bool(ds.rpcs)
        if not (len(cal) == len(spec) == count):
            raise ValueError(f"{scene_id} {seg_name}: 影像 {count} 波段 / 定标 {len(cal)} 行 / 光谱响应 {len(spec)} 行 不一致")
        if np.any(np.diff(spec[:, 0]) <= 0):
            raise ValueError(f"{scene_id} {seg_name}: 中心波长非严格递增")
        segments.append(Segment(seg_name, files[img_key], files[rpc_key], files[cal_key],
                                files[spec_key], cal, spec, count, width, height, dtype, has_rpc))

    scene = Scene(scene_id, satellite, directory, meta, files, segments)
    scene.band_table = build_band_table(scene, cfg)

    expected = cfg["sensors"].get("known_sensors", {}).get(satellite)
    if expected:
        if scene.total_bands != expected.get("total_bands"):
            scene.warnings.append(f"总波段数 {scene.total_bands} 与已知规格 {expected['total_bands']} 不符")
        if segments[0].count != expected.get("vnir_bands") or segments[1].count != expected.get("swir_bands"):
            scene.warnings.append("VNIR/SWIR 分相机波段数与已知规格不符")
    if not all(s.has_rpc for s in segments):
        raise ValueError(f"{scene_id}: 影像缺少可读取的 RPC（检查 .rpb 侧车文件是否与 .tif 同名）")
    if (segments[0].width, segments[0].height) != (segments[1].width, segments[1].height):
        scene.warnings.append("两相机影像尺寸不同，正射后统一到同一目标网格")

    if files["ogp"]:
        scale = float(meta.get("ogp_scale_factor") or 10000)
        scene.ogp = _read_ogp(files["ogp"], scale, cfg["sensors"]["ogp_band_order"])
        xml_sza = meta.get("solar_zenith")
        if xml_sza:
            scene.ogp["xml_center_solar_zenith_deg"] = float(xml_sza)
            scene.ogp["solar_zenith_diff_deg"] = abs(scene.ogp["center"]["solar_zenith_deg"] - float(xml_sza))
    else:
        scene.warnings.append("缺少 OGP 角度表，逐像元观测几何不可用")
    return scene


def discover_scenes(input_root: str | Path, cfg: dict) -> tuple[list[Scene], list[dict]]:
    """扫描数据根目录下所有产品 XML；不完整场景记录为错误而不是静默跳过。"""
    root = Path(input_root)
    if not root.is_dir():
        raise FileNotFoundError(f"输入根目录不存在: {root}")
    xmls = sorted(p for p in root.rglob("*.xml") if p.is_file())
    scenes, errors = [], []
    for xml in xmls:
        try:
            scene = parse_scene(xml, cfg)
            scenes.append(scene)
        except Exception as exc:
            errors.append({"xml": str(xml), "error": f"{type(exc).__name__}: {exc}"})
            log.warning("场景解析失败 %s: %s", xml, exc)
    scenes.sort(key=lambda s: s.scene_id)
    return scenes, errors


def file_inventory(scene: Scene, with_hash: bool = False) -> list[dict]:
    items = []
    for key, p in scene.files.items():
        if p is None or not Path(p).is_file():
            continue
        p = Path(p)
        rec = {"role": key, "path": str(p), "bytes": p.stat().st_size, "mtime_utc": p.stat().st_mtime}
        if with_hash:
            from .config import _sha256
            rec["sha256"] = _sha256(p)
        items.append(rec)
    return items


def scene_brief(scene: Scene) -> dict:
    return {
        "scene_id": scene.scene_id,
        "satellite": scene.satellite,
        "product_level": scene.meta.get("product_level"),
        "start_time": scene.meta.get("start_time"),
        "center_lon_lat": [scene.meta.get("center_lon"), scene.meta.get("center_lat")],
        "gsd_m": scene.meta.get("gsd_m"),
        "source_bands": scene.total_bands,
        "segments": {s.name: {"bands": s.count, "size": [s.width, s.height], "dtype": s.dtype,
                              "rpc": s.has_rpc} for s in scene.segments},
        "candidate_bands": len(scene.candidate_bands),
        "wavelength_range_nm": [scene.band_table[0]["wavelength_nm"], scene.band_table[-1]["wavelength_nm"]],
        "ogp_center_angles": scene.ogp.get("center"),
        "warnings": scene.warnings,
        "directory": str(scene.directory),
    }
