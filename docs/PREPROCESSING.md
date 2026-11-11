# 预处理链与运行语义

本文补充 README 中未展开的配置与 run_id 语义。

## 配置分层

| 文件 | 作用 | 是否入库 |
| --- | --- | --- |
| `config/pipeline.yaml` | 仓库默认模板，全部为相对路径与占位值 | ✅ |
| `config/pipeline.local.yaml.example` | 本机配置示例，展示 `${ENV:...}` 与 6SV 情景写法 | ✅ |
| `config/pipeline.local.yaml` | 你自己填写的本机配置 | ❌（`.gitignore`） |
| `config/sensors.yaml` | 传感器文件匹配模式、XML 字段别名、已知波段数 | ✅ |
| `config/scene_assets.example.yaml` | 可选资产适配层示例 | ✅ |

路径规则：

- 相对路径相对于 `project_root` 解析；`project_root` 本身相对于**配置文件所在目录**解析，模板中写 `..` 即仓库根目录。
- `${ENV:VAR}` 只在 `load_config()` 中展开；**未设置该环境变量时直接抛错**，不会回退到错误目录。
- 代码里不得出现绝对路径；`tests/test_backends.py` 会检查代码与配置中是否混入本机路径标记。

## stage 与断点续跑

`run` 的阶段顺序为 `radiance → envi → quac → reflectance → qc`，每个阶段完成后写入 `stage_<name>.json` 标记。相同 `run_id` 且相同输入/配置/SMOKE 波段数时可续跑，已完成的阶段跳过。

以下情况**必须使用新的 run_id**，流水线会在写 manifest 之前拒绝复用：

- 输入产品或输入目录变化
- 配置文件任一字段变化（配置会写入 sha256）
- 后端变化
- `--max-bands` 等 SMOKE 参数变化

`--force` 用于显式强制重算相同配置；`--stages` 可只跑指定阶段。

## 产品等级的执行约束

QC 报告正文不含硬编码结论：它读取本次运行的 `product_grade`，据此生成结论边界。修改后端时如果只改数值、不改 tag 与 QC 分派，会被视为缺陷而不是"小改动"。

| 后端 | 输出目录 | 关键 tag |
| --- | --- | --- |
| `iarr` | `relative_reflectance_final/` | 明确标注相对反射率 |
| `sixs_lut` | `sixs_reflectance_final/` | 模型估计地表反射率，标注未独立验证 |
| `envi_quac` | `reflectance_final/` | QUAC 统计估计 |
| `none` | 仅 `*_radiance.tif` | 辐亮度 |

## 独立阶段

主 CLI（`scripts/run_ahsi.py`）**默认不做**地形光照校正。以下脚本相互独立、可单独重跑，各自写入指定输出目录：

- `scripts/run_terrain_correction.py`：DEM → 坡度/坡向 → `cos(i)` → 逐波段 C 回归 → SCS+C；低光照、DEM 无效或 C 参数退化的像元保留上游反射率并写 QA 位。
- `scripts/spectral_pair_qc.py`：同点光谱质检，使用上凸包做连续统去除。
- `scripts/spectral_chain_comparison.py`：RPC 反算输出格网点到原始相机像素，按 `output_band → (segment, source_band)` 重排 DN，输出四阶段配对曲线。VNIR/SWIR 在输出立方体中是**交错**的，不能按 `VN 1..N + SW 1..M` 顺序拼接原始 DN。

## 换研究区清单

1. `preflight --input-root <新目录>` 预演，读 `preflight.json/.md` 的人工待办。
2. 准备覆盖新区的椭球高 DSM 与 EGM2008 网格，放到 `auxiliary_root`。
3. 核对 `grid.crs`（`auto_utm: true` 时自动推导）、`grid.rpc_dem`、`grid.geoid_grid`。
4. 核对 `sensors.yaml` 的文件匹配模式与 XML 别名；必要时补 `known_sensors`。
5. 用新 `run_id` 跑；`golden_root` 不匹配时 `verify` 自动记 `SKIP`。

详细已验证/未验证事项见 [VALIDATION.md](VALIDATION.md)。
