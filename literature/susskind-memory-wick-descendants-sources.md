# Ordered Wick 后代、生成元 masa 与零频质量：文献接口

核读日期：2026-11-07。适用背景沿用 size-bridge：同一组静态 Gaussian 耦合、标准 CAR、β=μ=0 归一迹、固定次数先于双标度极限，λ>0、q=e^(−4λ)∈(0,1)。本文件只记录原文定理、由它们导出的可用判据及尚待证明的接口，不宣布 common-J 后代的无对角质量已经解决。

## 1. 可直接调用的 q-Gaussian 定理

[É. Ricard, *Factoriality of q-Gaussian von Neumann algebras*, math/0311413](https://arxiv.org/pdf/math/0311413)，§3 Theorem 1；此前第二量子化讨论给生成元子代数的条件期望。对 tracial q-Gaussian 代数 M=Γ_q(K_R)，单位向量 e 对应 h=s_q(e)，A=W*(h) 是 masa。条件期望在 Fock 向量上投影到 F_q(Re)。实际核读 PDF 第 4–6 页的定理与证明；网页对 ≥ 字符的提取不可靠，维数条件以后一文明确的 dim K_R≥2 为准。

[P. Bikram、K. Mukherjee, *Generator masas in q-deformed Araki-Woods von Neumann algebras and factoriality*, 1606.04752](https://arxiv.org/html/1606.04752)，Eq. (17)、Theorems 5.3–5.4。取 U_t=I 即 tracial 情形。固定 q∈(−1,1)，生成元 masa 的补空间具有对二维 Lebesgue 测度绝对连续的 left-right measure；不要求小 q。有限 Wick words x,y 若 E_A(x)=E_A(y)=0，则

\[
T_{x,y}:L^2(A)\longrightarrow L^2(A),\qquad
T_{x,y}(a)=E_A(xay)
\]

是 Hilbert–Schmidt。证明的归一 q-Hermite 尾含 |q|^(2n) 和固定阶 q-factorial 比，相邻比趋 |q|²。有限 Wick 稠密性随后扩展至整个补空间。这里引用的是原文既有结果；任意补空间向量的密度只保证 L¹，不能将有限 Wick words 的 L² 密度强加给任意 L² 向量。

**对本项目的条件推论。** 若某 marked 残量的左右能量作用，已等距识别为 L²(M)⊖L²(A) 中的向量，则它的两能量谱无 E=E′ 对角质量。仅需无对角结论时，masa 已足够：标准 L²(M) 中 h-central 向量恰为 L²(A)，故 ker(L_h−R_h)=L²(A)。后一等式还可由有限迹下 affiliated L² 算符的谱截断得到，不是对任意外部 bimodule 的结论。

当前 ordinary H 的 q-Gaussian 极限，尚不等于所有 common-J rank-marked fluctuation fields 都在这个标准 L²(M) 中。可接受的实现证明是：构造候选向量，核对所有左右能量多项式生成向量的共同 Gram 矩阵，并把能量回归与 E_A 对齐。有限阶矩吻合不够。

## 2. 不依赖 Fock 实现的平方可和判据

以下是适用于本项目正两能量测度的直接数学推论，并非上述论文已经证明我们的 marked 系数满足它。

设 ρ=ρ_(q,C) 是方差 C>0 的 q-Gaussian 能量分布，支撑为紧区间 I。采用 monic 多项式

\[
H_{n+1}^{(q,C)}(E)=E H_n^{(q,C)}(E)-C[n]_qH_{n-1}^{(q,C)}(E),\qquad
p_n(E)=\frac{H_n^{(q,C)}(E)}{\sqrt{C^n[n]_q!}}.
\]

设 ν_R 是一个残量的有限正两能量测度，支撑 I×I。若

\[
M_{mn}=\int p_m(E)p_n(E')\,d\nu_R(E,E'),\qquad
\sum_{m,n\ge0}|M_{mn}|^2<\infty,
\]

则 ν_R=f_R(ρ⊗ρ)，f_R∈L²(ρ⊗ρ)，且 ν_R(Δ)=0，Δ={(E,E′):E=E′}。

证明：平方可和系数在 L²(ρ⊗ρ) 中定义 f_R；与 ν_R 对所有双变量多项式的积分相同。紧支撑上多项式稠密于连续函数，故两有限测度相等。正性给 f_R≥0，ρ 无原子给 (ρ⊗ρ)(Δ)=0。对实际代数向量 R，M 是 T_(R*,R) 的矩阵，至多交换两指标，因此这与 Hilbert–Schmidt 判据相同。

可实施的目标是从真实 common-J marked diagrams 推出归一系数的全阶可和界，例如有限带宽与 |M_mn|≤K(1+m+n)^d r^(min(m,n))，其中 0<r<1、|m−n|≤D，K,d,D 对固定后代有限。**该示例界在本文件未获证明**。固定数量的矩、有限截断平方和及未归一的 q-Hermite 系数都不能代替此条件。

所有常数可依赖固定 λ；λ→0 导致 q→1 的统一控制是另一个问题，不能从固定 q 定理直接取得。

## 3. 有限后代与全时间级数的闭包

在共同正谱空间中记 P_Δ=1_Δ(L_h,R_h)。若已证明 P_Δξ_a=0，则有限和

\[
\sum_{a,b} f_{ab}(L_h)g_{ab}(R_h)\xi_a
\]

也无对角质量，只需 f_ab,g_ab 在能量支撑上有界。这是 P_Δ 与左右能量函数对易的直接结果；多项式在紧支撑上自动有界。共同矩阵谱测度的正性同时控制不同 ξ_a 的 cross 对角项。若后代级数在同一 Hilbert 空间中范数收敛，P_Δ 的有界性把结论传至级数极限。

这里有两个必须实际完成的步骤：先证基本残量无对角质量，再证 ordered coupling-Gaussian Wick 操作确实落入上述共同谱空间与闭包。不能把 diamond 当成普通左右能量乘法；也不能从“某个同律 q-Gaussian 变量存在”跳到 marked 后代的共同实现。

“weakly contained in coarse”本身不是本文件使用的判据。所需的是实际 normal coarse 实现、标准 q-Gaussian 补空间实现，或直接的绝对连续性/平方可和证明；仅有弱包含不能排除奇异谱质量。

## 4. 零对角质量如何传递 exact pinching

这是正测度弱收敛的直接应用，可避免额外交换无限时间平均与双标度极限。

对有限 N 的随机矩阵 R_N、H_N，定义 annealed 正测度

\[
\nu_N^R(B\times D)=\mathbb E\,\tau_N\!\left(R_N^*P_N(B)R_NP_N(D)\right).
\]

若 \(\mathcal P_{H_N}R_N=\sum_E P_N(\{E\})R_NP_N(\{E\})\) 是包含整个简并能量块的 exact pinching，则

\[
\mathbb E\|\mathcal P_{H_N}R_N\|_{2,\tau_N}^2=\nu_N^R(\Delta).
\]

若有限正测度 ν_N^R 弱收敛到 ν_R，且 ν_R(Δ)=0，闭集 Δ 的 Portmanteau 不等式立即给 pinching 的 annealed L² 范数趋零。若有已确定的守恒能量回归 g_N(H_N)，应先令 R_N=X_N−g_N(H_N) 再用此式。这里“精确”指每个样本上使用真实谱投影；收敛结论是 annealed L²，并推出概率收敛，不自动推出几乎处处收敛。

## 5. EP 与 BINT 能复用什么

[E. G. Effros、M. Popa, math/0303045](https://arxiv.org/pdf/math/0303045)，Theorems 3.2–3.3：谱 q-Wick 多块期望使用块内禁配对；算符乘积还带 singleton 位于弧下的权重。此前已核原文 PDF 第 8–11 页，本轮复核本地 [neutral-marked-series 来源](susskind-memory-neutral-marked-series-sources.md)与 [size-bridge 来源](susskind-memory-size-bridge-sources.md)。这是成熟的组合框架；其谱 Wick 次序不是本项目 coupling-Gaussian W_j。

[BINT, 1811.02584v2](https://arxiv.org/html/1811.02584v2)，本轮重新核 §2.2 Eq. (2.6)–(2.8)：probe 随机耦合与 Hamiltonian 的耦合独立，同 flavor 才相互配对。§4.1、§4.2 的核及 Appendix C.4 的多点接口沿用本地 [OTOC 来源](susskind-singlet-otoc-sources.md)、[首阶来源](susskind-first-order-sources.md)的已核定位。

独立 probe 权重的导数只改变其允许图上的权重，不会补回 common-J 复合场与背景 H 之间本来允许的耦合收缩。因此必须先证明图集与权重匹配，才可使用该导数作为我们的 mark；独立 probe 的 AC 核不是 U_j 或 ordered descendants 的自动替代。

当前可推进的接口是：保留真实 coupling-Wick blocks 的禁配对与同 J 收缩，用已经建立的 CAR/size 局部 mark 生成两能量矩；随后证明第 2 节的全阶界，或构造第 1 节要求的共同实现。一般残量的无对角质量、以及由此导出的 exact pinching 收敛，在该接口完成前仍是待证任务。

## 6. 同轮后续：已完成的实际接口

以上段落保留本轮先查文献时的条件边界。随后在[正式稿](../manuscript/susskind-memory-wick-descendants/main.tex)中完成了直接图证明，无需假设标准 q-Fock 实现：对H^r(U_v diamond Z_w)固定field内部partial matching，再对开放端点置换求和，得到确定性能量回归。两个field之间无桥时，两个置换和独立，K12=0，中心化后整组抵消；q^N只给每个内部图组乘q^d，不破坏抵消。外部q-Hermite块禁止内部pair，余图至少一条field桥与全部energy桥交叉，给有限带宽与poly(n)q^n界，满足第2节的HS判据。

因此本阶段已经证明每个固定后代的回归残量两能量谱有L²(ρ⊗ρ)密度。由第4节传递exact pinching，再用已有全时间范数尾得到完整neutral核的收敛构造；全时间密度仅需L¹，未声明统一L²。谱N极限只对固定外部时间，级数截断尾才对全部实外部时间一致。charged defect及其cross不包含在这项结论内。原文定理仍作为成熟方法引用，不被改写为已有文献已经计算了我们的common-J核。
