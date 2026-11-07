# 固定电荷梯度尾：来源与继承接口

2026-11-08。先核原文适用范围与已有项目公式，再写新的逐扇区估计。

## 外部原文

- L. Ben Efraim–F. Lust-Piquard, *Poincaré type inequalities on the discrete cube and in the CAR algebra*, [arXiv:math/0702233v1](https://arxiv.org/html/math/0702233v1)，本轮重开并定位 §6.4 Theorem 6.8，亦核其 CAR 代数和归一矩阵迹背景。该定理给 full-CAR Riesz/symmetrized-gradient 控制；并非固定 charge 的定理。本轮不宣称扩展了任意场的 sector Riesz 定理，而是对实际需要的字段建立直接梯度包络。
- M. Berkooz–V. Narovlansky–H. Raj, *Complex Sachdev–Ye–Kitaev model in the double scaling limit*, [arXiv:2006.13983v2](https://arxiv.org/html/2006.13983v2)，本轮重开 §5，确认 (5.1)–(5.2) 固定电荷投影的背景。它提供 complex SYK 固定电荷文献入口，不替代当前有限 N 的矩形梯度估计。

## 已有项目证明，按原范围继承

- [纯 tag/受限 covariance](../notes/susskind-bilocal-fourth-tail-diagram-draft.md) §§1–4.1：任意交叉 Gaussian contraction 的 ancilla 范数、有限 membership 条件的 covariance、纯 Wick no-internal-pair 双标签条件、repeated deficits。新的内容是在梯度矩中明确保留 r 个成对 gradient index，并以 graded half-unit shift 保住 membership。
- [operator draft](../notes/susskind-bilocal-fourth-tail-operator-draft.md) §§1–5：归一 charge channels、精确 covariance path、OU/Chebyshev top-chaos 界、A_v telescope。没有将这些已有公式改称新发现。
- [旧 charged 证明](../notes/susskind-bilocal-fourth-tail-charged-proof.md) §§2–3、5–7：D/L 与通道交换子、A_v导数、实际eta recurrence、inverse Jacobi coefficient差分及ordinary五项重组。旧§4使用fulltrace Riesz的步骤由本轮直接梯度+矩形CP估计替换。
- [central ordinary](../notes/susskind-central-ordinary.md)：半填充 fixed-degree 源与同J强桥，接 exact pinching 与固定多项式制备；本轮提供它所缺的actual infinite-series tail。
- [固定多项式轨道密度](../notes/susskind-slater-polynomial-density.md)：固定 Gaussian degree 制备的一致有限 N 密度矩，不允许把制备次数随着展开次数增长后仍沿用同一常数。

## 本轮证明与边界

CAR分次导数移位、归一通道梯度主项收缩、受限tagged gradient计数和normalized-trace rectangular CP Holder都是正文逐步给出的有限矩阵推导。不以文献搜索不见同样写法作为原创性证据。最终中央窗口保留charge差分中的1+|X_N|限制，完整selected bilocal的neutral/rank/memory-return尚需独立检查。
