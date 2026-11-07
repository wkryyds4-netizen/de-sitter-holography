# Joint memory：来源、继承结果与本轮接口

核读与整理：2026-11-07。先使用已核前驱与原文，然后计算其未覆盖的共同电荷回归。没有进行无关的宽泛文献扫查。

1. [Berkooz–Narovlansky–Raj, arXiv:2006.13983v2](https://arxiv.org/html/2006.13983v2)，本轮核对 §4.2 的 (4.4)–(4.9)：化学势给局部乘积态、single-site trace 与 double-scaling 源计数。原文讨论的 1/sqrt(N) 化学势缩放对应本项目 1/p 的电荷源尺度；本项目还保留 whole-field loss 的 nilpotent source 与有限 singleton 修正。原文 CAR=2，本项目 CAR=1，不照搬其能量归一。标量 charge Gaussian 极限、局部迹/Poisson 方法属于既有工具；g_n(X,H)、内部弦charge回归及其对实际记忆场的应用在本项目内推导。

2. [Effros–Popa, math/0303045](https://arxiv.org/abs/math/0303045)，前驱已核 Theorems 3.2–3.3 的 multiblock q-Wick 规则。本轮用它消除外部 q-Hermite 块的内部配对；仅能删除不属于 marked incident set 的弦。不能把它当作任意矩阵系数的交换律。Mahonian 平均、q-Hermite 正交与 compact moment determinacy 也是既有数学工具。

3. [oriented transfer](../notes/susskind-memory-oriented-transfer.md) 给 exact D 接口、full-time weighted E-G 替换、ordinary K_n 尾、DD/DL 局部源与charge背景。新证明先消singleton再取Poisson极限；不对内部匹配权 C^i q^c_gamma 再加 q 导数。混合源使用较小 neutral mark 子段的 incident 集 I_L。

4. [neutral descendants](../notes/susskind-memory-wick-descendants.md) 给 V_j 递推、确定性回归、固定场双能量 L2 密度与有限时间 BB。此次复用其谱界组织；补齐 charged AA 和 mixed AB/BA。正矩阵测度的 Cauchy–Schwarz、闭集 Portmanteau 及 total variation 闭包是标准工具，不声称发明新的谱论。

5. [two-component](../notes/susskind-memory-two-component.md)、[uniform-bound](../notes/susskind-memory-uniform-bound.md) 和 [bilocal-channel](../notes/susskind-bilocal-channel.md) 给实际有限矩阵的 M(s)、F=A+B+o(L2)、共同耦合与精确 pinching。这里先使用有限 N 普通级数尾，随后才做极限 Hermite 重求和；没有仅凭低阶矩重建整个函数。A0B2 的非零 Mehler 校准是前驱结果，本轮独立恢复它。

后续物理入口另见 [memory-to-bilocal 来源笔记](susskind-memory-to-bilocal-sources.md)。之前已核的 Susskind/FLZ/Itakura bulk 结果仍是匹配目标，完整 memory covariance 不是已经完成的 de Sitter 质量字典。本轮不作原创性优先权宣称，也不新发对外稿件。
