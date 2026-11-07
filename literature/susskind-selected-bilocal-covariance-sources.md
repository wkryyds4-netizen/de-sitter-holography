# 制备协方差：已有数学与本模型接口

2026-11-08。本阶段沿用此前核读的原文；收尾再核 arXiv 书目信息。原文数学、项目先前结果与本次两插入推导分别列明。

- Bożejko–Kümmerer–Speicher，[q-Gaussian processes: non-commutative and classical aspects](https://arxiv.org/abs/funct-an/9604010)，Commun. Math. Phys. 185 (1997),129–154，§§2、4：q-Fock、second quantization、Markov/Mehler 结构。条件积分的正性、收缩与半群来自这些已有结构；不把两步 Mehler 视为新数学。
- Effros–Popa，[Feynman Diagrams and Wick Products Associated with q-Fock Spaces](https://arxiv.org/abs/math/0303045)：whole-Wick 禁配对及 crossing 计数的前驱。实际 D/L 源如何耦合仍由本项目局域 CAR 计算决定，不能仅凭抽象 Wick 法则指定。
- 上一阶段[来源记录](susskind-selected-bilocal-mean-sources.md)已经核对 Casper–Kolb–Yakimov 的双变量 q-Hermite（r=q）及 de Boor 的 Genocchi–Hermite 插值余项。本轮完整 beta 和 factorial 尾直接继承，不重复列为新公式。

模型证明输入：

1. [多源局域性](../notes/susskind-bilocal-multisource-proof.md) §§2–6：nilpotent 源、有限 N 支配及本 D/L 类高源累积量消失。
2. [高阶 Wick 与 selected 密度](../notes/susskind-slater-higher-wick.md)：真实单选态 rho、W-W 同弦及 D-W 零连接。
3. [单插入投影](../notes/susskind-selected-bilocal-mean.md)：one-U 投影、完整 Xi_j、实际双局域时间滤波与均值。
4. [中性制备响应](../notes/susskind-neutral-selected-response.md)：相邻四 U 的两步 Mehler 恒等式。
5. [共同回归与谱](../notes/susskind-memory-joint-regression-proof.md) §§5–10、[有序后代](../notes/susskind-memory-wick-descendants.md)：同一原始耦合的 charged/neutral/mixed 正矩阵谱与绝对收敛时间级数。
6. [固定电荷实际双局域](../notes/susskind-sector-bilocal.md)：有限 N 精确 pinching、半填充、actual/local 强桥及制备下的时间极限。

本次模型内推导把 rho−f 的两个边界源单独展开，得到完整两场协方差的谱基线加显式制备修正；继而得共同中心时间退相位、t⁶/t⁷ 及符号反例。Riemann–Lebesgue、L² 稠密性和 Cauchy–Schwarz 均是标准工具。未声称优先权查新已经完备，也未把模型内结果识别成有限 G 引力动力学。
