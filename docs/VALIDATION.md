# 实现与验证记录

## A. 自动化测试（CI 持续执行）

| 项目 | 内容 | 证据 |
| --- | --- | --- |
| 纯单元测试 | IARR 数值与确定性、坏波段与旧因子拒绝、后端注册与契约、6SV 反演公式与 SRF 采样、坡度/坡向/光照/SCS+C 几何、配置校验 | `tests/`、`.github/workflows/ci.yml` |
| golden 回归 | 自备 L1 场景数据时执行；未配置时整类跳过 | `tests/test_golden.py` |
| 字节编译 | `ahsi/` 与 `scripts/` 全部模块 | CI |
| 元数据检查 | 配置解析、默认配置无本机路径、仓库未跟踪栅格数据 | CI |

CI 在 Python 3.10 / 3.11 / 3.12、Ubuntu 上运行；`rasterio` / `fiona` 依赖官方 wheel 内置的 GDAL。作者本地在有 L1 数据的环境下全量 49 项通过；无数据的 CI 中 golden 类自动跳过，其余 41 项执行。

## B. 处理链实现记录

| 环节 | 实现位置 | 本机实跑输出 |
| --- | --- | --- |
| L1 解析与辐射定标 | `ahsi/scene.py`、`ahsi/core.py` | 330 波段辐亮度立方体、逐波段质量 CSV |
| RPC + DEM 正射 | `ahsi/core.py`、`scripts/prepare_rpc_dem.py` | 30 m UTM 格网、联合覆盖掩膜、RGB 预览 |
| IARR 归一化 | `ahsi/relative.py` | 相对反射率立方体、因子 JSON、确定性复算 |
| 6SV 大气校正 | `ahsi/sixs_lut.py` | 762 行 LUT、254 波段反射率、三档情景敏感性表 |
| ENVI QUAC（可选） | `ahsi/atmos.py` | QUAC 产品与阶段记录 |
| 地形光照校正 | `ahsi/terrain.py`、`scripts/run_terrain_correction.py` | SCS+C 产品、QA 位、逐波段 C 参数 |
| 同点光谱链对比 | `scripts/spectral_chain_comparison.py`、`scripts/spectral_pair_qc.py` | DN / 辐亮度 / 大气校正 / 地形校正四阶段配对曲线与 CSV |
| 质量门禁 | `ahsi/qc.py` | `qc_report.md`、`qc_checks.json`、`manifest.json` |
| 资产对接 | `ahsi/assets.py`、`scripts/build_scene_asset.py` | 按产品等级注册的 catalog |

## C. 可扩展环节

以下处理可按研究区与数据条件继续接入，配置与 stage 接口已留出位置：

- 基于 GCP / 参考影像的 RPC 精化正射与定位精度验收
- 水汽与气溶胶的影像反演，丰富当前由配置给定的大气情景
- 天空漫射、邻近地形反射、投射阴影与 BRDF 的二阶地形辐射校正
- 去条带、云/雪/阴影逐像元判别、光谱平滑与波长重校准
- 按 GF5B SRF 重卷积的参考光谱库与影像端元分类

## D. 运行记录

每个 `run_id` 的 `manifest.json` 保存配置快照、环境版本、输入清单、后端与产品等级；`qc_checks.json` 保存逐项检查结果；`stage_*.json` 保存各阶段完成标记。换数据或改配置时使用新的 `run_id`，即可完整复现该次处理。
