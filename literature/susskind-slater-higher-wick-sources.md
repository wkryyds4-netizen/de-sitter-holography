# 高阶 Wick 选态接口：前驱与使用范围

2026-11-08。本轮先回读已封存source、selected密度、charged谱与tail证明，再做新计算。

- Berkooz–Narovlansky–Raj，Complex Sachdev-Ye-Kitaev model in the double scaling limit，[arXiv:2006.13983v2](https://arxiv.org/html/2006.13983v2)，§3及§5。原文提供complex弦图与fixed-charge投影背景；不是本文same-J selected高阶loss源或pinching的现成证明。本次重新打开原文，核§5 (5.1)–(5.2)按电荷投影的定义。
- Ben Efraim–Lust-Piquard，Poincaré type inequalities on the discrete cube and in the CAR algebra，[math/0702233v1](https://arxiv.org/html/math/0702233v1)，§6.1与§6.4 Theorem6.8(1),(3)。本次重新核原文：CAR Schatten/normalized trace上的reverse gradient inequality；未给本项目fixed-charge条件迹版本。不能只将sector trace写到同一不等式下。旧完整fullFock尾使用该定理仍有效。
- 项目slater-composite及selected-state：local D/L矩阵、normalized-defect收敛与single-selected rho1；本轮推广mark长度和subblock，不将旧结论改写为新发现。
- 项目memory-joint-regression-proof §§5–7：J_d回归、无桥分组相消、外部Hermite矩的有限带宽及q衰减、L²两能量谱和exact pinching。此次先用allG源相等独立接central矩，再移植谱论论证；没有使用其fullcharge→central自动推论。
- 项目slater-polynomial-density：固定Gaussian degree的orbitdensity矩界；本轮应用于J_d及固定selected多项式制备，常数不可不经估计用于degree→∞。
- 项目bilocal-fourth-tail-diagram-draft与operator-draft：pure/repeated tagged loss、normalized sector CP path、coupling OU/Chebyshev估计已逐central-sector成立。ordinary charged-proof的Clifford Riesz乘积环节仍需新的central处理。

本轮新内容：四端点局部压缩与完整有限代数证书、任意固定whole-Wick D/L子块的boundary恒等式、charged塔限制态f(h)，以及据此转接已有charged谱到central/selected逐阶actual pinching。没有宣称新q-Mehler公式、所有高阶场有界、普通内部pair已解决或完整deSitter字典完成。
