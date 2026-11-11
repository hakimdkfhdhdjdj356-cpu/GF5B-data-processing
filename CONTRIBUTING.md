# 贡献指南

欢迎提 Issue 和 Pull Request。这是一个科研工程底座，最看重的是**可复现**和**产品等级不被混淆**，其次才是新功能。

## 环境

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e .
python -m unittest discover -s tests -v
```

纯单元测试不依赖任何外部数据；`tests/test_golden.py` 需要自备 L1 场景数据，没有时会自动跳过。

## 提交前检查

- [ ] `python -m unittest discover -s tests -v` 通过（允许 golden 类 skip）
- [ ] `python -m compileall -q ahsi scripts` 无错误
- [ ] 没有把影像、DEM、运行产物或本机绝对路径加入提交
- [ ] 没有在代码里写死卫星名、研究区、路径或结论
- [ ] 新增/修改后端时，同步更新 `METHOD_INFO`、QC 检查项与结论边界措辞
- [ ] 任何数值改动都给出可复现的输入、配置与对比结果

## 代码约定

- 机器相关路径**只能**出现在配置文件；支持 `${ENV:变量名}`，未设置时报错而不是静默回退。
- 所有落盘产物都要写产品等级标签；相对反射率不得标注为地表反射率。
- QC 报告正文不得硬编码结论；结论边界应随产品等级自动改写。
- 新增传感器型号只改 `config/sensors.yaml`，不改 Python 代码。
- 保持 `run_id` 的断点续跑语义：输入、配置或 SMOKE 波段数变化必须换 run_id。
- 测试用 `unittest`，不引入额外测试框架依赖。

## Pull Request

请说明：改动动机、影响哪些 stage / 后端、如何验证（命令与结果）、是否改变已有产物的数值或标签。若涉及数值内核，请附与既有成果的抽样对比。

## 分支与提交

`main` 保持可用。提交信息用简短的祈使句，例如 `fix: reject non-finite 6S coefficients before LUT write`。

## 行为准则

见 [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md)。安全问题请按 [SECURITY.md](SECURITY.md) 私下报告，不要公开提 Issue。
