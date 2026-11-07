# 输出 rank、费米子 size 与共同耦合弦图：原文接口

2026-11-07。本轮实际核读下列原文，并对照[最终阶段稿](../manuscript/susskind-memory-size-bridge/main.tex)。这里记录可复用的公式及其适用条件，不作新结果优先权声明。本项目使用标准 CAR、完整重叠 complex Gaussian p-body、同一组静态耦合，先固定次数，再取 p²/N→λ；W_j 始终是 coupling-Gaussian chaos 投影。

## 1. Tarakemeh–Xu：sector rank 的 SU(2) 表示

[M. Tarakemeh、S. Xu，arXiv:2511.16634v1](https://arxiv.org/html/2511.16634)，正文 Eq. (3)、(5)，Appendix C Eq. (23)–(26)。Eq. (5) 的辅助自旋超算符为

\[
J_+X=\sum_i c_i^\dagger Xc_i,\qquad
J_-X=\sum_i c_iXc_i^\dagger,\qquad
J_z=(Q_L+Q_R-N)/2.
\]

原文定义非阿贝尔 size 为 s=N/2−J。**匹配本项目 Casimir 归一后**，neutral 空间满足

\[
\mathcal C_x=N(N+2)/4-\mathbf J^2,
\qquad \mathcal R=N/2-J,
\qquad \mathcal L=(D+\overline D)/2=N/2-J_x.
\]

这里 R 是 charge-sector 的 U(N) 表示 rank；L 是 full-Fock Clifford 长度的一半，两者不是同一个有限 N 算符。它们的差为 J−J_x。前两项与本项目 Casimir 的匹配、后续 ladder/overlap 估计属于本阶段的接口推导，不把它们误引成原文已证明的 DSSYK 结论。

Appendix C 对含 m_c 个 c/c†、m_n 个 2n−I 的局部基元，以 Clebsch–Gordan 系数证明固定 support 的大 N 关系 s→m_c/2+m_n。该值是 Majorana 长度的一半。应用时以 Eq. (24)/(25) 为准；Eq. (26) 的最小 J 下标漏印了 1/2，末段关于局部计数的文字也须按上述整体因子二解读。

原文的固定 support 极限不能直接覆盖本项目 O(p)=O(√N) support，逐基元结论也不能自动用于任意相干叠加：Q−N/2 是 rank 0、Clifford 长度 2 的反例。本阶段用共同耦合的 overlap 与 ladder 界处理 Z_j。原文 Brownian 模型的 master equation 不用于本项目静态耦合动力学。

## 2. Carlen–Lieb：物理 Clifford 长度的生成映射

[E. A. Carlen、E. H. Lieb，arXiv:hep-th/9209053](https://arxiv.org/pdf/hep-th/9209053)，Eq. (1.3)/(1.4)，PDF 第 3 页。其 fermionic oscillator semigroup 在 Clifford 单项上满足

\[
P_t^{\rm Cl}\Gamma_A=e^{-t|A|}\Gamma_A.
\]

这是已建立的长度乘子，不是新的弦图核。本阶段据此定义

\[
\mathcal M_j(t)=e^{tjp}P_{t/2}^{\rm Cl}Z_j,
\qquad \mathcal M_j'(0)=(jp-\mathcal L)Z_j.
\]

将前因子按原始费米子 incidences 分配后，singleton 权重精确抵消；双出现 occupation 的局部映射为 n→n+(e^t−1)I/2，1−n→1−n+(e^t−1)I/2。这些是本阶段对既有 semigroup 的代数应用。这里 t 是 source 参数；Clifford OU 与作用于随机耦合的 Gaussian OU 不同。长度映射须作用于 full-Fock lift，再限制 charge；不能仅从一个 sector 矩阵重建它。

## 3. Berkooz–Narovlansky–Raj：有向 CAR word 与 DS 交叠规则

[M. Berkooz、V. Narovlansky、H. Raj，arXiv:2006.13983v2](https://arxiv.org/html/2006.13983v2)，§3.2 Eq. (3.11)/(3.12)，§3.3 Eq. (3.13)–(3.15)。前两式在有限 N、一般 index-set lengths 下给出有向 word 的交叉符号、局部交替条件及多重交叠因子；非交替的局部 word 迹为零。§3.3 随后才使用固定数量 p-set 的 pair-Poisson 极限及高重叠消失。

本阶段先保留这套有限 CAR 局部结构，再插入 Clifford-size source；符号在局部取迹之前计算，不能因零 source 的局部迹为零就删掉 marked 项。原文费米子规范与本项目标准 CAR 的尺度差已由 C_N 归一匹配；这里采用 λ=p²/N、q=e^(−4λ)。

可直接复用的是原始 word 规则及固定图的交叠极限。两 mark 的 K(I)、同一 site 的 λK_12 项、共同耦合 Z_j 缺额分布及其 rank 强桥接是本阶段组合和估计的对象，不归给原文现成公式。中央固定 charge 窗的定量转换继续使用本项目旧证书；独立随机 probe 的相关函数不能代替这些 common-J 复合场。

## 4. Effros–Popa：多块谱 q-Wick 的既有组合框架

[E. G. Effros、M. Popa，arXiv:math/0303045](https://arxiv.org/pdf/math/0303045)，Theorem 3.2（PDF 第 8–9 页）、Theorem 3.3 与 Eq. (11)（第 10–11 页）。前者将多块 q-Wick 期望写成禁止块内连线的完美配对和，带 q^crossing 与协方差积；后者给出允许 singleton 的算符乘积展开，其指数还包括弧下 singleton 的退化交叉数。

该文提供成熟的多块配对框架与谱 q-Wick 乘法。它的 Wick ordering 是 q-Fock/spectral ordering，**不是**本项目有限 N 的 coupling-Gaussian W_j；不能据此直接断言二者相等。本阶段 coupling-Wick 块内部禁止收缩来自原始耦合 Gaussian Wick 定理，物理 size mark 来自上述 CAR/Clifford 接口。既有 q-Hermite、Mehler 与 Mahonian 公式均按既有工具使用。

## 本阶段与后续边界

[neutral-marked-series](../notes/susskind-memory-neutral-marked-series.md)提供固定次数 rank 高矩及全时间级数尾；[central-spectrum](../notes/susskind-memory-central-spectrum.md)提供固定 marked graph 的中央迹转换；[oriented-kernel](../notes/susskind-memory-oriented-kernel.md)提供 D/Db 局部 mark 的前驱。

本阶段稿把这些输入用于 rank–half-length 强桥接、固定次数的一/二 mark 弦图矩、固定 j 的离散缺额律及两能量矩问题。有限系数与局部矩阵枚举验证恒等式，不能替代其 uniform overlap、ladder 和矩确定性论证。一般 j 的零频质量、实际 samplewise pinching 的极限、ordered Wick 后代及 charged reverse-source 的共同核仍待处理；本阶段没有得到 finite-G 引力四点动力学或 de Sitter 体质量/位置字典。
