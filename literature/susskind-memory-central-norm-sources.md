# 中央相干范数的来源

2026-11-06实际在线打开并读取Kota–Vyas，https://arxiv.org/pdf/1504.00759 ，§7.1前后的粒子移除定义、式(101)–(106)及相邻particle-space矩接口。这一文献使用与H独立的随机transition operator；我们的A由同一H构造，不能照搬其独立性假设。复用的是有限sector/矩形算符和U(N)协方差的框架。

项目已有一般charge的charged covariance与norm transfer是本轮直接输入：

- `manuscript/susskind-memory-coherent-cross/main.tex` 第2节的T_a(t)、两charge相干关系。
- `manuscript/susskind-memory-charge-lift/main.tex` 的精确ell_r²传递。
- `manuscript/susskind-memory-residual-norm/main.tex` 的一般charge Casimir交换插入及半填充norm相消；本轮显式重算n−m项，未机械延用半填充式。
- `results/susskind-memory-central-window/central-moments.json` 的中央紧窗六阶output矩。
- `results/susskind-memory-two-charge/symbolic-filter.json` 的gH精确前因子，以及上一轮central-scale的gA紧窗展开。

新增计算：一般charge的coherent A norm整理、中央x的p²/p相消、实际A输入尾、三方向Gram和强gA替换，再组装W2core。含output-deficit的完整Gram与动态谱仍待下一步；不作原创性优先权结论。
