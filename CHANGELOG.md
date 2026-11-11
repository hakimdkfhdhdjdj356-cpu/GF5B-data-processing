# 更新日志

本文件记录本仓库（公开版本）的显著变更。格式参考 [Keep a Changelog](https://keepachangelog.com/zh-CN/1.1.0/)，版本号遵循语义化版本。

## [2.0.0] - 2026-11-11

首个公开发布版本。由内部工程底座整理而来，**未改动任何数值内核**；公开化只涉及脱敏、可移植性与仓库规范化。

### 新增

- 独立仓库与 GitHub 社区标准文件：README、LICENSE (MIT)、`.gitignore`、CONTRIBUTING、SECURITY、CODE_OF_CONDUCT、CITATION.cff、Issue / PR 模板、GitHub Actions CI、`pyproject.toml` 与 `requirements.txt`。
- `config/pipeline.yaml` 改为不含本机绝对路径的通用模板；新增 `config/pipeline.local.yaml.example` 与 `config/scene_assets.example.yaml`。
- `tests/test_golden.py` 在缺少本地 L1 数据时整类跳过，而不是报错或伪造通过；CI 因此可在无数据环境运行。

### 变更

- 移除配置、文档与测试中的本机绝对路径、人名、研究区标识与内部报告引用。
- 三个场景适配脚本重命名为通用名称：`build_scene_asset.py`、`plot_scene_qa.py`、`prepare_rpc_dem.py`；场景资产配置键统一为 `scene_assets`。
- `tests/test_backends.py` 的"无硬编码研究区"检查改为通用的"无本机绝对路径"检查。
- `ahsi/__init__.py` 文档字符串同步为四种后端（`iarr` / `sixs_lut` / `envi_quac` / `none`）。

### 已知限制

见 [README §7](README.md#7-已知边界必须随成果一起说明)：无 GCP 精配准、6SV 大气参数为假设情景、地形校正仅一阶直射归一化、模型估计反射率未经独立验证。

---

## 内部沿革（早于本仓库）

- 默认大气后端由 `envi_quac` 改为 `iarr`，整条链路脱离 ENVI/IDL 许可。
- 新增 `sixs_lut` 后端（6SV 查表）、SCS+C 地形光照校正、同点四阶段光谱链对比与资产读取层。
- 增加 `preflight` 可迁移性检查、run_id 隔离与断点续跑、统一质量门禁、golden 回归。
