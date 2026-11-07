# 中央中性场：原文与继承范围

2026-11-08。本次重新打开 [BNR 2006.13983v2](https://arxiv.org/html/2006.13983v2) 与 [Effros–Popa math/0303045](https://arxiv.org/abs/math/0303045)。它们分别提供固定电荷 complex SYK 与 q-Wick 图组合背景，不替代本稿的同样本强误差和精确能量投影。

## 本轮实际核过的项目接口

- [wick-rank 稿](../manuscript/susskind-memory-wick-rank/main.tex) §§5–6：固定 sector 归一 HS、正 rank transition、不同输入 rank 的确定性交叉正交、creation sym contraction、annihilation adjoint。不是 fulltrace 平均界。
- [central-window 稿](../manuscript/susskind-memory-central-window/main.tex) §§2–3：QQ[A,B,z] 六矩证书、抽出 epsilon 幂后非零分母、紧集一致性。fixed k 的 r≤kp 只扩大紧集，不将 k 作为增长量。
- [neutral-marked-series](../manuscript/susskind-memory-neutral-marked-series/main.tex) §§2–6：实际 eta recurrence、fixed-degree Dfrak²Z 与 Dfrak²H^n、graph Cauchy–Schwarz。原稿还有完整 rank-proxy 重求和；本轮不需要重走这一段。
- [covariance tangent](../manuscript/susskind-sector-covariance-tangent/main.tex)：固定 v 的 A_v same-copy O(p^-1)，保留 finite coefficients。
- [operator draft](../notes/susskind-bilocal-fourth-tail-operator-draft.md) §§5–7 与 [gradient tail](../manuscript/susskind-sector-gradient-tail/main.tex)：逐 sector 实际 eta/T_j 和 local proxy 的统一次数包络。有限次 Holder 先展开递推，阶数不随次数翻倍。
- [wick-descendants](../manuscript/susskind-memory-wick-descendants/main.tex)：通用局域 primitive、完整回归、无桥成组抵消、两能量密度及 Portmanteau。不是仅用若干低阶矩猜无对角质量。
- [higher-wick](../manuscript/susskind-slater-higher-wick/main.tex)：允许整个纯 Wick 源块、其中的 L 子标记、普通外部 H 与整体 chaos；all-G 局部源逐图相同给 central trace，极端 F 边界给 sigma1。标记区内部 Gaussian pairs 不在该定理范围，本轮 primitive 不含这种 pair。
- [polynomial-density](../manuscript/susskind-slater-polynomial-density/main.tex)：固定 Gaussian degree 的有限N轨道密度界；仅用于固定制备，不以其常数支配增长制备次数。

本轮组合的新关键是 (Psi−Theta)pE 的图范数插值。没有把标准 spectral Cauchy–Schwarz、CP 收缩、q-Hermite 回归等已有公式称为原创。
