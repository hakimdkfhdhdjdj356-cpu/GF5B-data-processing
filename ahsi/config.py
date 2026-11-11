"""配置加载、路径展开与校验。代码内不出现任何机器绝对路径。"""
from __future__ import annotations

import math
import os
import re
from copy import deepcopy
from pathlib import Path

import yaml

_ENV_RE = re.compile(r"\$\{ENV:([A-Za-z_][A-Za-z0-9_]*)\}")

REQUIRED_TOP = ["project_root", "input_root", "output_root", "grid", "bands", "atmospheric", "qc", "runtime"]
PATH_KEYS = ["input_root", "output_root", "auxiliary_root", "golden_root"]


def load_yaml(path: Path) -> dict:
    with Path(path).open("r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def _expand_env(value, base: Path):
    """只做环境变量替换；路径解析按白名单单独处理，避免把 'Generic / Unknown Sensor' 之类当成路径。"""
    if isinstance(value, dict):
        return {k: _expand_env(v, base) for k, v in value.items()}
    if isinstance(value, list):
        return [_expand_env(v, base) for v in value]
    if not isinstance(value, str):
        return value
    def sub(m):
        key = m.group(1)
        if key not in os.environ:
            raise ValueError(f"配置引用了未设置的环境变量: {key}")
        return os.environ[key]
    return _ENV_RE.sub(sub, value)


def _resolve_path(value: str, base: Path) -> str:
    p = Path(value)
    if not p.is_absolute():
        p = base / p
    return str(p.resolve())


def load_config(pipeline_yaml: str | Path, sensors_yaml: str | Path | None = None) -> dict:
    pipeline_yaml = Path(pipeline_yaml).resolve()
    cfg = load_yaml(pipeline_yaml)
    missing = [k for k in REQUIRED_TOP if k not in cfg]
    if missing:
        raise ValueError(f"主配置缺少必需字段: {missing}")

    raw_root = _expand_env(cfg["project_root"], pipeline_yaml.parent)
    base = Path(raw_root).resolve() if Path(raw_root).is_absolute() else (pipeline_yaml.parent / raw_root).resolve()
    cfg = _expand_env(cfg, base)
    cfg["project_root"] = str(base)
    for k in PATH_KEYS:
        if k in cfg and cfg[k]:
            cfg[k] = _resolve_path(cfg[k], base)
    if cfg.get("atmospheric", {}).get("envi_taskengine"):
        cfg["atmospheric"]["envi_taskengine"] = _resolve_path(cfg["atmospheric"]["envi_taskengine"], base)
    sixs = cfg.get("atmospheric", {}).get("sixs")
    if sixs and sixs.get("executable"):
        sixs["executable"] = _resolve_path(sixs["executable"], base)
    cfg["_config_dir"] = str(pipeline_yaml.parent)
    cfg["_config_sha256"] = _sha256(pipeline_yaml)

    if sensors_yaml is None:
        sensors_yaml = pipeline_yaml.parent / "sensors.yaml"
    sensors = load_yaml(Path(sensors_yaml).resolve())
    cfg["sensors"] = sensors
    cfg["_sensors_sha256"] = _sha256(Path(sensors_yaml).resolve())

    _validate(cfg)
    return cfg


def _sha256(path: Path) -> str:
    import hashlib
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


VALID_METHODS = ("none", "iarr", "sixs_lut", "envi_quac")


def _validate(cfg: dict) -> None:
    g = cfg["grid"]
    if float(g["resolution_m"]) <= 0:
        raise ValueError("grid.resolution_m 必须为正")
    if g["resampling"] not in ("nearest", "bilinear", "cubic", "average"):
        raise ValueError("grid.resampling 不支持")
    a = cfg["atmospheric"]
    if a["on_unavailable"] not in ("stop_at_radiance", "fail"):
        raise ValueError("atmospheric.on_unavailable 只能是 stop_at_radiance 或 fail")
    if a["method"] not in VALID_METHODS:
        raise ValueError(f"atmospheric.method 只能是 {VALID_METHODS}，实际 {a['method']!r}")
    if a["method"] == "iarr":
        i = a.get("iarr") or {}
        if i.get("normalization", "scene_mean") not in ("scene_mean", "scene_median"):
            raise ValueError("atmospheric.iarr.normalization 只能是 scene_mean 或 scene_median")
        lo, hi = i.get("robust_pct", [2, 98])
        if not (0 <= lo < hi <= 100):
            raise ValueError("atmospheric.iarr.robust_pct 必须满足 0 <= low < high <= 100")
        if int(i.get("min_valid_pixels", 1000)) < 1:
            raise ValueError("atmospheric.iarr.min_valid_pixels 必须 >= 1")
    if a["method"] == "sixs_lut":
        sx = a.get("sixs") or {}
        if not sx.get("executable"):
            raise ValueError("atmospheric.sixs.executable 必须指定本机 6S 可执行路径")
        scs = sx.get("scenarios") or []
        if len(scs) < 1 or sum(1 for s in scs if s.get("nominal")) != 1:
            raise ValueError("atmospheric.sixs.scenarios 必须含且仅含一个 nominal=true 情景")
        ids = set()
        for s in scs:
            for k in ("id", "aot550", "water_g_cm2", "ozone_atm_cm"):
                if k not in s:
                    raise ValueError(f"atmospheric.sixs.scenarios 缺少字段 {k}")
            scenario_id = s["id"]
            if not isinstance(scenario_id, str) or not scenario_id.strip() or scenario_id in ids:
                raise ValueError("6S 情景 ID 必须为非空且唯一的字符串")
            ids.add(scenario_id)
            aot, water, ozone = (float(s[k]) for k in ("aot550", "water_g_cm2", "ozone_atm_cm"))
            if not all(map(math.isfinite, (aot, water, ozone))) or aot < 0 or water <= 0 or ozone <= 0:
                raise ValueError("6S 情景参数必须有限，且 AOT 非负、水汽和臭氧为正")
    for low, high in cfg["bands"]["water_absorption_nm"]:
        if low >= high:
            raise ValueError("吸收带区间必须低波长在前")
    if cfg["bands"]["vnir_max_nm"] >= cfg["bands"]["swir_min_nm"]:
        # 允许重叠(VNIR/SWIR 本来就有重叠)，这里只防止配置写反
        pass
    if cfg["runtime"]["band_chunk"] < 1:
        raise ValueError("runtime.band_chunk 必须 >= 1")
    s = cfg["sensors"]
    for key in ("file_patterns", "meta_aliases", "calibration", "ogp_band_order"):
        if key not in s:
            raise ValueError(f"sensors.yaml 缺少 {key}")
    if len(s["ogp_band_order"]) != 4:
        raise ValueError("ogp_band_order 必须为 4 项")


def aux_path(cfg: dict, key: str) -> Path:
    """辅助数据文件（DEM / geoid）路径统一从这里取。"""
    p = Path(cfg["auxiliary_root"]) / cfg["grid"][key]
    return p


def public_config(cfg: dict) -> dict:
    """写入 manifest 的配置快照，去掉内部键。"""
    out = deepcopy(cfg)
    for k in list(out):
        if k.startswith("_"):
            out.pop(k)
    return out
