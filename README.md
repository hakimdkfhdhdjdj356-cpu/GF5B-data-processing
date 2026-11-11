# ahsi-pipeline

[![CI](https://github.com/hakimdkfhdhdjdj356-cpu/ahsi-pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/hakimdkfhdhdjdj356-cpu/ahsi-pipeline/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](pyproject.toml)

配置驱动的国产星载高光谱 **L1 预处理工程底座**：把 L1 DN 处理成定标、正射、质量受控的辐亮度／反射率立方体，并为后续的岩性、蚀变矿物识别提供物理上可追溯的输入。

> **English.** A configuration-driven Python pipeline for spaceborne hyperspectral L1 products (validated on GF-5B AHSI and ZY1E AHSI). It performs radiometric calibration, RPC + DEM orthorectification, four interchangeable atmospheric / normalisation backends, band screening and a unified quality gate, plus optional SCS+C terrain illumination correction. It is an engineering base, **not** a new classification algorithm, and it deliberately keeps product grades (radiance / relative reflectance / model-estimated surface reflectance) strictly separated.

**定位**：工程底座，不是新算法。相对早期人工链路，改进的是路径与阈值配置化、run_id 隔离、断点续跑、QC 判定、输入清单、golden 回归，以及「产品等级绝不混淆」。

**这个仓库不含任何影像数据**：L1 产品、DEM、参考影像与运行产物都需自备；`config/` 里只有模板和占位路径。

---

## 1. 能力

| 环节 | 实现 |
| --- | --- |
| L1 解析 | 按 `config/sensors.yaml` 配对 VN/SW 影像、RadCal、SpectralResponse、RPC、XML、OGP；校验波段数、波长单调、gain > 0 |
| 辐射定标 | 逐波段 `L = DN × gain + offset`，全波段回读校验 |
| 正射 | 显式 RPC + DEM（需自备椭球高 DSM）+ EGM2008 大地水准面；按景自动推导 UTM 带；双相机联合覆盖掩膜 |
| 大气 / 归一化 | `iarr`（默认，零外部依赖）、`sixs_lut`（6SV 查表）、`envi_quac`（可选商业软件）、`none` |
| 波段筛选 | 420–2400 nm、水汽吸收带、VNIR/SWIR 相机边缘，保留完整辐亮度立方体 |
| 质量门禁 | 输入配对、OGP 太阳角核对、RPC 往返残差、全波段回读、覆盖面积、后端可用性、按产品等级分派的检查；统一 PASS/WARN/FAIL |
| 地形校正（独立阶段） | DEM → 坡度/坡向 → `cos(i)` → 逐波段回归求 C → SCS+C / C-correction / cosine，含低光照 QA |
| 同点光谱链对比（独立阶段） | RPC 反算输出格网点到原始相机像素，`output_band → (segment, source_band)` 重排 DN，输出 DN / 辐亮度 / 大气校正 / 地形校正四阶段配对曲线 |
| 资产对接 | `ahsi/assets.py` 通过 catalog 向下游工具提供按产品等级的读取层 |

**默认整条链路不依赖 ENVI / IDL 或任何商业许可。** ENVI QUAC 只是可选插件。

---

## 2. 安装

```bash
git clone https://github.com/hakimdkfhdhdjdj356-cpu/ahsi-pipeline.git
cd ahsi-pipeline
python -m venv .venv && source .venv/bin/activate     # Windows: .venv\Scripts\activate
pip install -e .
```

`rasterio` / `fiona` 的官方 wheel 已内置 GDAL，通常不需要单独编译。使用 `sixs_lut` 后端时另需本机安装 6S V1.1 可执行文件与 `Py6S`（`pip install -e ".[sixs]"`）。

---

## 3. 快速开始

仓库里的 `config/pipeline.yaml` 是**不含本机路径的模板**，数据目录 `data/l1` 默认为空。请先复制示例并填写自己的路径：

```bash
cp config/pipeline.local.yaml.example config/pipeline.local.yaml
# 编辑该文件，或用 AHSI_* 环境变量提供路径
python scripts/run_ahsi.py --config config/pipeline.local.yaml preflight
python scripts/run_ahsi.py --config config/pipeline.local.yaml plan
python scripts/run_ahsi.py --config config/pipeline.local.yaml run --run-id demo_001
python scripts/run_ahsi.py --config config/pipeline.local.yaml verify
```

安装后也可直接用 `ahsi preflight` / `ahsi plan` / `ahsi run` / `ahsi verify`（等价入口）。

| 命令 | 作用 |
| --- | --- |
| `preflight` | 秒级。换数据/换研究区前列出哪些能自动、哪些必须人工准备 |
| `plan` | 秒级。扫描输入、校验文件配对、输出场景计划 |
| `run` | 正式跑全波段。约 30–60 分钟、20 GB 量级输出，建议后台执行 |
| `verify` | 与自备 golden 成果做非破坏性回归；无 golden 时记 `SKIP`，**不伪造通过** |

常用参数：`--run-id`（相同输入/配置可断点续跑）、`--scenes`、`--stages`、`--max-bands`（仅用于 SMOKE）、`--force`。代码、配置或 SMOKE 波段数改变后必须使用**新的** `--run-id`，流水线会在写 manifest 前拒绝复用。

独立阶段：

```bash
python scripts/run_terrain_correction.py --help
python scripts/spectral_pair_qc.py --help
python scripts/spectral_chain_comparison.py --help
```

---

## 4. 处理链（每个 stage 落盘、可独立复现）

| 阶段 | 关键处理 | 主要产物 |
| --- | --- | --- |
| 发现 | 产品目录配对与元数据校验 | `plan.json`、`scene_manifest.json` |
| `radiance` | 定标 + RPC/DEM 正射 + 覆盖掩膜 + 全波段回读 | `*_radiance.tif`、`band_parameters_and_quality.csv`、预览 |
| `envi` | 仅 `envi_quac` 需要：导出标准 ENVI 辐亮度 | `ENVI_input/*` |
| `quac` | 执行所选大气/归一化后端 | `IARR/iarr_factors.json`、`SIXS_LUT/sixs_lut.csv` 或 `QUAC_v2/*` |
| `reflectance` | 归一化 / 6SV 反演 / QUAC 缩放 | float32 光谱立方体、质量 CSV、flags |
| `qc` | 统一判定与按产品等级生成结论边界 | `qc_report.md`、`qc_checks.json`、`manifest.json` |

> 阶段名 `quac` 是历史遗留，实际含义是「大气/归一化后端执行」，保留是为了让既有 run_id 的断点续跑标记继续有效。

一个 run_id 自包含：`manifest.json`（配置快照 + 环境版本 + 输入清单 + 后端状态与产品等级）、`qc_report.md` / `qc_checks.json`、`run_results.json`、`scenes/<scene_id>/`。

---

## 5. 四种后端与产品等级（**不要混用**）

`config/pipeline.yaml` 的 `atmospheric.method`：

| method | 外部依赖 | 产品等级 | 典型用途 |
| --- | --- | --- | --- |
| `iarr`（默认） | 无 | **相对反射率，不是地表反射率** | 单景内连续统去除、吸收特征位置/深度、场景端元选择、MNF/PCA、相对异常检测 |
| `sixs_lut` | 本机 6SV1.1 + Py6S | 假设大气条件下的**模型估计**地表反射率 | 需要估计地表反射率的探索性分析 |
| `envi_quac` | ENVI 5.6+ taskengine + 许可 | QUAC 场景统计估计反射率 | 与既有 ENVI 成果对比 |
| `none` | 无 | 辐亮度 | 最保守交付、后续自行处理 |

各后端的**文件名、GeoTIFF tag、ENVI bbl、QC 检查项、报告结论边界互相独立**，任何情况下都不会把相对反射率标成地表反射率。

- **IARR**：逐波段除以场景有效像元的稳健均值（默认裁掉 2–98 分位外极值）。无外部数据、无经验大气参数、无随机性，同一输入逐位可复现。逐波段缩放会改变光谱向量方向，**不能**直接拿实验室绝对反射率光谱库做 SAM，也**不能**用于跨景/跨时相定量对比。
- **6SV**：从本景 XML 读取观测几何，用配置 DEM 中位高程和显式大气情景生成逐波段 LUT（`y = xa·L − xb; ρ = y/(1 + xc·y)`），高斯近似响应采样于 6S 从 0.25 µm 起的 2.5 nm 网格。AOT / 水汽 / 臭氧是**假设情景**，不是影像反演结果；生成 LUT 需要本地 6S，续跑可复用已生成的 LUT。

---

## 6. 换数据 / 换研究区

先跑 `preflight`，它会给出判定与人工待办清单。要点：

| 类别 | 项目 |
| --- | --- |
| 必须自备 | 覆盖研究区的椭球高 DEM（RPC 正射必需）、EGM2008 网格 |
| 必须自备 | 符合 `sensors.yaml` 匹配规则的 L1 产品目录（XML、VN/SW 影像、`.rpb`、定标与光谱响应文件、OGP 角度表） |
| 必须改配置 | `input_root`、`auxiliary_root`、`grid.rpc_dem` / `geoid_grid` |
| 按需改配置 | `grid.crs`（`auto_utm: true` 时自动推导）、`bands.*`、`sensors.yaml`、`golden_root` |
| **不用改** | 所有 `.py` 代码：无硬编码卫星名、研究区、路径或结论 |

场景按 XML 逐个独立发现，**不要求多景齐全**；新增传感器型号只需改 `sensors.yaml` 的文件匹配模式、XML 字段别名与已知波段数。`golden_root` 只对验证过的那两景有效，换数据后 `verify` 自动记 `SKIP`。

---

## 7. 已知边界（必须随成果一起说明）

1. 辐亮度**不是**地表反射率；默认 IARR 产品**不是**地表反射率；6SV 产品是模型估计值，**未经独立地面光谱精度验证**。
2. RPC 往返残差只用于内部数值检查。参考影像精化正射 / GCP / CE95 仍未实现，与参考影像的相对偏移未校正，因此**不能报告绝对定位精度**。
3. 地形校正阶段只做一阶太阳直射归一化，**未做**天空漫射、邻近地形反射、投射阴影、BRDF、逐像元 OGP 观测角；也未做去条带、云/雪/阴影逐像元判别。覆盖掩膜 ≠ 可信裸岩范围。逐波段 C 参数会改变光谱向量方向，用于 SAM / 光谱库匹配前须重新评估。
4. 候选波段是按配置规则的探索性初筛，不等于逐波段信噪比认证。
5. 6SV 大气情景为假设值；水汽/AOT 影像反演、波长重校准、光谱平滑均未实现。参考光谱库若按其它传感器 SRF 采样，需按本传感器 SRF 重卷积后再匹配。
6. GDAL 对 `RPC_DEM_SRS` 的 PROJ 字符串可能打印非致命 `proj_create: unrecognized format` 日志；严格的垂直基准验收需单独做。
7. 本仓库只提供预处理底座，**不含**任何岩性/蚀变分类算法，也不对下游分类精度作承诺。

---

## 8. 与下游「创新层」的接口

按产品等级分别读取：

- 默认：`relative_reflectance_final/<scene_id>_relative_reflectance.tif`（float32，NoData = −9999，带 `wavelength` / `fwhm` / `bbl` 标签）
- 可选：`reflectance_final/<scene_id>_surface_reflectance.tif`、`sixs_reflectance_final/<scene_id>_6SV_surface_reflectance.tif`、地形校正产品

配套 `relative_reflectance_band_quality.csv` 与 `quality_flags.tif`。用相对反射率时端元必须来自影像本身（或先做连续统去除），不能直接套实验室绝对反射率光谱库。

---

## 9. 测试与 CI

```bash
python -m unittest discover -s tests -v
python -m compileall -q ahsi scripts
```

测试分两类：

- **纯单元测试**（IARR 数值与确定性、后端契约、6SV 反演公式与 SRF、地形几何与 SCS+C）不需要任何外部数据，在 CI 中始终执行。
- **golden 回归**（`tests/test_golden.py`）需要自备 L1 场景数据；未配置 `input_root` 时整类 `skip`，**不伪造通过**。作者本地在有数据的环境下全量 49 项通过。

CI 见 [`.github/workflows/ci.yml`](.github/workflows/ci.yml)，在 Python 3.10 / 3.11 / 3.12 上运行。

---

## 10. 仓库结构

```
ahsi/          核心包（配置、场景发现、定标正射、后端、QC、地形、资产读取）
scripts/       CLI 入口与独立阶段脚本
config/        主配置模板、传感器适配规则、本机配置示例
tests/         unittest 套件
docs/          预处理链与结论边界的补充说明
.github/       Issue / PR 模板与 CI
```

`scripts/build_scene_asset.py` 与 `scripts/plot_scene_qa.py` 是**示例适配器**：它们演示如何把产物与自备的 EMIT / ASTER / Sentinel-2 参考数据做交叉比较，需要 `scene_assets` 配置（见 `config/scene_assets.example.yaml`），不是主流水线的必需环节。其中的对比状态描述来自作者的一次具体实跑，换数据后应按自己的几何检查结果修改，不要当作通用结论。

---

## 11. 数据与隐私

仓库中**没有**任何卫星影像、DEM、研究成果、钻孔或地质图数据，也没有本机绝对路径、人名或研究区标识。所有路径都通过配置文件或环境变量提供；`config/*.local.yaml`、`.env`、`data/`、`outputs/` 等均已被 `.gitignore` 排除。请勿把原始数据或凭据提交到仓库。

---

## 12. 许可与引用

MIT，见 [LICENSE](LICENSE)。若用于学术论文，请引用本仓库（见 [CITATION.cff](CITATION.cff)）并说明所用后端的产品等级与上述结论边界。
