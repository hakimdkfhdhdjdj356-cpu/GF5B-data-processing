## 改动内容

<!-- 一句话说明改了什么 -->

## 动机

<!-- 为什么需要这个改动；关联 Issue：# -->

## 影响范围

- 影响阶段：<!-- preflight / plan / radiance / envi / quac / reflectance / qc / terrain / spectral chain / 无 -->
- 影响后端：<!-- iarr / sixs_lut / envi_quac / none / 无 -->
- 是否改变已有产物的数值或 GeoTIFF tag：**是 / 否**
- 是否新增配置项或改变配置语义：**是 / 否**（如是，请同步更新 `config/*.example` 与 README）

## 验证

```bash
python -m unittest discover -s tests -v
python -m compileall -q ahsi scripts
```

<!-- 贴出结果摘要；涉及数值内核时请附与既有成果的抽样对比 -->

## 检查清单

- [ ] 测试通过（允许 golden 类因缺少本地数据而 skip）
- [ ] 没有把影像、DEM、运行产物、凭据或本机绝对路径加入提交
- [ ] 代码里没有写死卫星名、研究区、路径或结论
- [ ] 新增/修改后端时，`METHOD_INFO`、QC 检查项与产品等级标签已同步
- [ ] 相对反射率没有被标注为地表反射率
- [ ] `CHANGELOG.md` 已更新（用户可见变更）
