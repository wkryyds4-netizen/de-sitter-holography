# 统一二阶界：已有工具与本轮桥接

2026-11-07。先查 Gaussian 矩、Jacobi 与 q-Wick 文献，再写本轮定量推导。以下原文已在本轮阅读；主稿不以检索未发现同式作为原创优先权证据。

- J. A. Tropp, *Second-Order Matrix Concentration Inequalities*，作者站[原文](https://tropp.caltech.edu/papers/Tro16-Second-Order-Matrix-preprint.pdf)，Prop.2.3/(2.4)、§7。使用标准 Gaussian Schatten 矩工具。本稿针对 `Σ B_a²=CI` 用 Gaussian integration by parts 与谱基 weighted Young 自行写出 `(2q−1)!! C^q` 矩界，从而得到同副本 OU force 的 `ε³|s|⁴√11!! C³/24` 余项；该 force 公式并非原文现成结论。
- NIST DLMF [§18.2](https://dlmf.nist.gov/18.2)，三项递推、正交多项式及 Jacobi 结构。`Q_{j+1}=xQ_j−Ca_jQ_{j−1}` 属标准框架。本模型 `a_j=Σ μ(vp)`、与 coupling chaos 的定量差、全复平面 ODE 延拓是这里具体展开的推导。辅助 Jacobi 测度不等于有限 N 的实际 DOS。
- M. Bożejko, B. Kümmerer and R. Speicher, *q-Gaussian processes: non-commutative and classical aspects*，[funct-an/9604010 原文](https://arxiv.org/pdf/funct-an/9604010)，Prop.2.9、Thm.2.11。标准 Wick tensor/条件期望可求真正 q-Gaussian 极限中的领先能量背景。本轮有限 N 统一界不依赖将该极限公式代入有限 N。后续若使用带外部 u 的条件期望，必须保留两外侧与内部的配对。
- J. Tall and S. Tomsovic, *Double-scaled bosonic and fermionic embedded ensembles, complex SYK, and the dual Hilbert space*，[2604.14522v2 §3.1](https://arxiv.org/html/2604.14522v2#S3.SS1)，(3.1)–(3.2)及后续 Wick 结构。原文 spectral Wick polynomials 是 annealed 谱测度的正交多项式，不是有限 N Gaussian coupling 最高 chaos 投影。不能把两个定义直接视为恒等；本轮 Jacobi 比较是为补这个接口。

本地依赖依次为 `susskind-memory-wick-rank` 的全次数缺额与 number 估计、`susskind-memory-charged-transfer` 的正反缺陷及 polar 传递、`susskind-memory-return-background` 的首返守恒强近似。前三份证据在本轮 recorder 中重新核 immutable hashes。新稿证明指定 microscopic field 的统一二阶界及两点替换；实际有限时间极限谱与体引力字典仍开放。
