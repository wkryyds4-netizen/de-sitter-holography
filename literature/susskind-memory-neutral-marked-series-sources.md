# neutral marked series：原文与可复用接口

核读日期：2026-11-07。本轮先核原文与本地已有证书，再确定新计算；不把标准 Wick/Jacobi 公式列为原创发现。

1. Edward G. Effros、Mihai Popa，*Feynman Diagrams and Wick Products Associated with q-Fock Space*，[原文 PDF](https://arxiv.org/pdf/math/0303045)。Theorems 3.2–3.3，PDF第8–11页：多block q-Wick 期望的配对不得连接同block内部点，按crossing计权；乘积规则另计singleton在弧下的贡献。可用于之后的marked contractions。不能据此直接把物理rank缺额D等同q-Fock number。此前内存中的36个一变量线性化检查未落盘，本阶段不将其计为新验证证据。

2. Micha Berkooz、Mikhail Isachenkov、Vladimir Narovlansky、Genis Torrents，*Towards a full solution of the large N double-scaled SYK model*，[原文](https://arxiv.org/html/1811.02584)。§4.1 Eq.(4.8)/(4.9)为uncrossed双probe，§4.2/§5给crossed vertices与R matrix，Appendix C.4 Eq.(C.8)给独立probe权重的多点表达。对权重求导是可探索的计算方式，但必须先证明导数插入与本模型Uj或charged defect一致。

3. Jiuci Xu，*On Chord Dynamics and Complexity Growth in Double-Scaled SYK*，[v3原文](https://arxiv.org/html/2411.04251)，2025-03-22。§3.4 Eq.(101)、(106)–(110)给chord-number生成函数与crossed四点接口。这里的chord number不是D=pNcal−Rcal。直接微分Eq.(101)可见其k阶矩为(−1)^k∂μ^kI|0；Eq.(102)所印额外1/k!与定义不一致。本轮未依赖该式证明我们的结论。

本地前驱：uniform-bound的有限参数Jacobi系数、rank阶段的一次commutator及全k二矩、central-window的Q[A,B,z]六阶相消证书、two-component的A+B及低阶cross。新稿显式区分调用的旧结论和本轮推广。尚未在以上原文找到直接给出本项目共同J rank-marked/charged联合涨落核的定理。

下一步可借这些标准公式组织固定-degree图，但先核真实rank mark及Gaussian Wick block如何映射；零频质量与精确pinching要单独分析。
