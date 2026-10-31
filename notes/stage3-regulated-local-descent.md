# 原局部模式核的受调节弱下降：生成 FP 空间上的创建与伴随

日期：2026-10-29。模型为 dS₂、L=1、G→0、ν>0 的自由实物质与主系列量子 observer；固定 proper-time Gaussian 宽度 σ>0、实际左右 observer 下界 s*>0，逐次保留质量投影。先处理 proper SO⁺(2,1)，采用 KL 正文模式及相位约定。

**本稿给出一个明确的受调节弱场实现：先由原 (5.4) 已核的运动学模式核定义共同紧群 cutoff 下的矩阵元，再证明极限存在、在真实产品生成的 FP 闭包中给出唯一输出，并与 normalized multiplicity 映射一致。创建映射不必保持强轨道域 D_orb。**

论证采用[全有限粒子数生成空间](stage3-all-n-generated-space.md) N11 已建立的各 J_N，作为独立输入；本文不替代该识别的证明。结合下文的弱下降，最终得到这个生成空间上、完整有限粒子数域中的本质自伴场。这里不宣称未调节位置积分逐点存在，也不把生成闭包自动等同任意形式不变分布构成的最大空间。

## 1. 输入和精确前提

主要原式是 [KL 2411.08090v1](https://arxiv.org/html/2411.08090v1) (3.14)、(4.8)—(4.10)、(4.42)、(5.4)—(5.6) 与附录 E。模型的局部 observer 矩阵元为

\[
\langle a,n_a|\phi_R|b,n_b\rangle
=\int d^2x\sqrt{-g}\,\overline{\Psi^a_{n_a}(x)}\,
\varphi(x)\Psi^b_{n_b}(x),
\qquad
\varphi^+(x)=\sum_n\overline{\psi_n^\nu(x)}a_n^\dagger.
\tag{L1}
\]

每次 Gaussian 涂抹与原负号时间演化给创建模式因子

\[
h_{\sigma,\tau_0}(a,b)
=e^{-\sigma^2[m(a)-m(b)]^2/2}
 e^{-i[m(a)-m(b)]\tau_0},\qquad m(s)=\sqrt{s^2+1/4}.
\]

下文的运动学输入，是[创建审计](stage3-creation-descent-audit.md) C3—C15 已给出的固定模式识别和有界延拓：

\[
W_N=H_{R,\ge s_*}\otimes\operatorname{Sym}^N H_\nu,
\qquad
C_N:W_N\longrightarrow W_{N+1},
\]
\[
C_N=\sqrt{N+1}(1_R\otimes P_{\mathrm{Sym}^{N+1}})
(T_{\sigma,\tau_0}\otimes1),
\quad
\|C_N\|\le b_N:=\sqrt{(N+1)\mathfrak M}.
\tag{L2}
\]

C_N 是 SO⁺(2,1) 有界等变映射，且在有限原模式、紧质量弱矩阵元上与 (L1) 同一创建核相等。该输入采用经原模式积分重建的顶角，不把 literal A.7 的不相容常数默默当作已刊正确公式。C_N* 与相同实时间窗的原湮灭模式核相匹配。

令 K=SO(2)，W_N⁰ 为右 observer 加全部物质的总 K=0 子空间。使用[全 N 生成空间稿](stage3-all-n-generated-space.md) N2—N3 给出的一致正文相位下的完整酉分解

\[
W_N^0\simeq\int_{d>0}^{\oplus}dd\,\mathcal M_{N,d},
\qquad
\mathscr H_N=\int_{d\ge s_*}^{\oplus}\rho_{\rm Pl}(d)dd\,
\mathcal M_{N,d},
\tag{L3}
\]

并且 C_N 在这些分解中为可测 fiber 映射 c_N(d)，满足

\[
\|c_N(d)\|\le b_N\quad\text{几乎处处}.
\tag{L4}
\]

所有内部 principal/离散通道及 Bose 投影均包含在 M_N,d 中。总 discrete 没有 K=0，不影响 (L3)。N=0 时总标签与实际右 observer 标签相同，不能人为再引入独立内部质量积分；本稿将该扇区的既有独立识别列为输入。

全 N 生成空间稿的 N11 给出等距满射

\[
J_N:\mathscr H_N\longrightarrow\mathcal H_N^{gen},
\tag{L5}
\]

其中 H_N^gen **定义为所选真实产品种子群平均像在位置 FP 范数中的生成闭包**。具体要求如下：

- 原种子为左磁指标 0 的 C∞紧质量包，乘右总 K=0 的有限原模式、C∞紧实际右质量包；Bose 化只含有限项。这些种子已验证属于位置强轨道域 D_orb。
- 它们的 GA 幅度在谱空间 𝓗_N 中稠密，且 J_N 把这些幅度映到原位置群平均像，FP 内积等于 GA 内积。
- 在固定左包下，任意 W_N⁰ 向量都可由上述右原模式包作运动学逼近。此项是普通右 Hilbert 空间的稠密性；它仅与下一节的明确连续界联合使用。

N=2 的幅度与密度证明见[产品种子稿](stage3-product-seed-ga-density.md)，位置充分等距见[FP 稿](stage3-fp-two-particle-inner-product.md)。全 N 生成空间稿对 N≥1 以 PP/PD 完整酉分解归纳，并重复同一 K=0 密度证明；N=0 以 d=a 独立处理。真实 observer 产品的轨道条件由允许任意酉物质 Hilbert 空间的[正频共同域定理](stage3-positive-frequency-orbit-core.md)提供。因此这里不是仅从 N=2 的结果自动外推。

## 2. 扩大的可控 GA 测试类

记完整运动学空间为 K_N=H_L⊗W_N，其对角群表示为 U_N(g)。定义

\[
\mathcal E_N=\operatorname{span}\{
|f,0\rangle_L\otimes v:
 f\in C_c^\infty((s_*,\infty)),\ v\in W_N^0\}.
\tag{L6}
\]

这里右向量 v 可以有无限原磁指标及非紧内部谱，不要求属于 D_orb。对 S=f⊗v，已有波包版 Schur 公式给 GA 幅度

\[
\mathcal A_NS(d)=K_P(d)f(d)v(d),
\quad K_P(d)=\rho_{\rm Pl}(d)^{-1},
\quad
\|\mathcal A_N(f\otimes v)\|_{\mathscr H_N}
\le B_f\|v\|,
\]
\[
B_f^2=\operatorname*{ess\,sup}_{d>0}K_P(d)|f(d)|^2<\infty.
\tag{L7}
\]

对有限线性组合按线性延拓。若写成不同的张量和，已证 GA 交叉配对公式保证 A_N 表示同一元素，而不是逐分解任意赋值。

本节所需的绝对收敛是：[产品种子稿](stage3-product-seed-ga-density.md) P3—P10 的统一球函数界 |D₀₀^d|≤Ξ₀≤C(1+λ) exp(−λ/2)，结合左质量先积分后的任意阶分部积分，给

\[
\int_G |L_{f_1f_2}(g)|\Xi_0(g)dg<\infty.
\]

右总谱密度 R₁₂(d)=〈v₁(d),v₂(d)〉属于 L¹(dd)，故所有 E_N 中向量的 GA 交叉配对绝对收敛，并且

\[
Q_{GA}(S_1,S_2)=
\int_G\langle S_1,U_N(g)S_2\rangle dg
=\langle\mathcal A_NS_1,\mathcal A_NS_2\rangle_{\mathscr H_N}.
\tag{L8}
\]

这个论证只使用总 K=0 主系列和完整 multiplicity，形式上与 N 无关，但仍要求 (L3) 成立。

因为 C_N 与群作用对易，它与 K 对易，所以 C_NW_N⁰⊂W_N+1⁰；其伴随也保持相反方向的 K=0 子空间。因此

\[
\mathcal C_N:=1_L\otimes C_N:
\mathcal E_N\to\mathcal E_{N+1},
\qquad
\mathcal C_N^*\mathcal E_{N+1}\subset\mathcal E_N.
\tag{L9}
\]

这不是 D_orb 不变性。E_N 的作用是提供可控的弱测试配对；不要求它等于全运动学位置测试函数空间。

## 3. 先按局部核和共同 cutoff 定义弱矩阵元

取紧支撑可测函数 0≤χ_R≤1，χ_R(g⁻¹)=χ_R(g)，且 χ_R→1。设

\[
\eta_{N,R}S:=\int_G\chi_R(g)U_N(g)S\,dg.
\tag{L10}
\]

每个有限 R 都是普通运动学 Hilbert 空间中的 Bochner 积分。η_N,R 有界且自伴，通常不把种子变成严格 singlet；本稿不对这些有限 R 非不变态套用位置 FP 不变态公式。

**定义。** 对 Φ∈E_N+1、Ψ∈E_N，局部创建弱矩阵元取

\[
B_N^{loc}(\Phi,\Psi)
:=\lim_{R\to\infty}
\langle\Phi,\mathcal C_N\eta_{N,R}\Psi\rangle_{kin}.
\tag{L11}
\]

这里 C_N 来自 (L1) 的已核运动学弱核。式 (L11) 是定义的起点，不以乘法空间中的 c_N RHS 作为定义。没有假定未调节时空积分可以直接作用于群平均后的分布。

还可以把这一起点直接落实为有限模式矩阵。由原有限磁指标、C∞紧实际质量包向量经Gram–Schmidt，分别在W_N与W_N+1选递增有限秩正交投影P_j、Q_j，强收敛到恒等算符。每个Q_j C_N P_j的有限矩阵元均由同一Gaussian、同一逐次质量截止下的原(L1)弱积分给出，且

\[
Q_jC_NP_j\to C_N\ \text{强收敛},\qquad
\|Q_jC_NP_j\|\le b_N.
\tag{L11a}
\]

在完整运动学空间使用1_L⊗(Q_j C_N P_j)。固定R先取j→∞，由普通Hilbert连续性恢复L11中的有限R配对，再取R→∞。有限j的矩阵一般不等变，故只有恢复C_N之后才使用L12；这里没有交换两个极限，也没有把有限模式截断后的非不变态代入FP切片式。


有限 R 时，bounded intertwining 与 Bochner 积分直接给

\[
\mathcal C_N\eta_{N,R}\Psi
=\eta_{N+1,R}\mathcal C_N\Psi.
\tag{L12}
\]

由于 Φ、C_NΨ 都属于 E_N+1，(L8) 的绝对收敛适用。于是 (L11) 存在，并且由支配收敛对任意上述共同 cutoff 独立：

\[
B_N^{loc}(\Phi,\Psi)
=\int_G\langle\Phi,U_{N+1}(g)\mathcal C_N\Psi\rangle dg
=Q_{GA}(\Phi,\mathcal C_N\Psi).
\tag{L13}
\]

这里交换的是已经有界的运动学映射与有限 Haar 积分，再取绝对收敛的标量极限；并未交换原无穷模式和、时空积分与群规约的全部次序。

## 4. 计算该弱型，并得到物理范数界

对单项 f⊗v，等变分解和 (L7) 给

\[
\mathcal A_{N+1}\mathcal C_N(f\otimes v)(d)
=K_P(d)f(d)c_N(d)v(d)
=c_N(d)\mathcal A_N(f\otimes v)(d).
\]

它也对有限和成立。因此从 (L13) 独立计算出

\[
\boxed{B_N^{loc}(\Phi,\Psi)
=\int_{d\ge s_*}\rho_{\rm Pl}(d)
\langle\mathcal A_{N+1}\Phi(d),
 c_N(d)\mathcal A_N\Psi(d)\rangle\,dd.}
\tag{L14}
\]

这才是与 multiplicity RHS 的识别步骤。由 (L4)，

\[
|B_N^{loc}(\Phi,\Psi)|
\le b_N\|\mathcal A_{N+1}\Phi\|\,
\|\mathcal A_N\Psi\|.
\tag{L15}
\]

(L15) 是物理 GA/FP 范数的界。它依赖 fiber c_N 的统一界，**不是**把 (L7) 中依赖左包的 B_f 误当作一般物理算符范数。

若 AinΨ=0，(L14) 给该输入的全部弱矩阵元为零；若 AoutΦ=0，测试端也为零。所以局部弱型下降到 GA 零范数商，不依赖种子的不同运动学代表。

由 (L5) 和物理测试像的稠密性，Riesz 定理给唯一输出 Y_Ψ∈H_N+1^gen，满足

\[
\langle J_{N+1}\mathcal A_{N+1}\Phi,Y_\Psi\rangle_{FP}
=B_N^{loc}(\Phi,\Psi),
\qquad
\|Y_\Psi\|_{FP}\le b_N\|J_N\mathcal A_N\Psi\|_{FP}.
\tag{L16}
\]

因此可定义 C_N^weak(J_NA_NΨ)=Y_Ψ，再唯一有界延拓到整个 H_N^gen。

## 5. 输出也可由真实强域产品在 FP 范数中构造

这个构造解释为何不需要 C_N 保持 D_orb。先取 Ψ=f⊗v。令 w=C_Nv∈W_N+1⁰，选右原模式产品包有限和 w_j→w，于普通右 Hilbert 范数。根据 (L5) 的种子前提，Γ_j=f⊗w_j 都是已验证的强域种子；它们的实际位置群平均 ηΓ_j 已存在。

由 (L7)及位置 FP–GA 等距，

\[
\|\eta\Gamma_j-\eta\Gamma_k\|_{FP}
=\|\mathcal A_{N+1}(\Gamma_j-\Gamma_k)\|
\le B_f\|w_j-w_k\|\longrightarrow0.
\tag{L17}
\]

所以它们在 H_N+1^gen 中有极限。另一逼近列给相同极限，因二者与 w 的差都由 B_f 控制。与任意真实产品测试的 FP 配对，再利用 (L8) 在固定左包下的连续性，极限恰等于 (L13)。由物理测试稠密性，该极限就是 (L16) 的 Y_Ψ。

有限多个左包的和可逐项逼近；误差上界是相应 B_f‖w_j−w‖ 的有限和。不同张量分解、不同运动学代表的独立性最终由 (L14)—(L16)给出，不靠逐项 B_f 界猜测。

这一证明只在强域中使用 ηΓ_j；没有声称输出种子 C_NΨ 本身满足逐轨道 Bochner-L¹ 条件，也没有把未经证明的强域不变性藏进记号。

## 6. 伴随、物理创建界与谱表示

为湮灭定义同一共同 cutoff 的弱型

\[
D_N^{loc}(\Psi,\Phi)
:=\lim_R\langle\Psi,\mathcal C_N^*\eta_{N+1,R}\Phi\rangle_{kin}.
\tag{L18}
\]

C_N* 有界等变且保持 E，所以该极限同样存在。有限 R 时，χ_R 逆元对称使 η_N,R 自伴；由 (L12)的伴随关系，

\[
\langle\Psi,\mathcal C_N^*\eta_{N+1,R}\Phi\rangle
=\overline{\langle\Phi,\mathcal C_N\eta_{N,R}\Psi\rangle}.
\]

取极限并使用同一 FP 稠密测试空间，得到完整物理块上的

\[
D_N^{weak}=(C_N^{weak})^*,
\qquad
\|C_N^{weak}\|,\ \|D_N^{weak}\|
\le\sqrt{(N+1)\mathfrak M}.
\tag{L19}
\]

完成弱型和输出构造之后，才把算符写成等价的谱表达

\[
C_N^{weak}
=J_{N+1}\left(\int^{\oplus}c_N(d)\,\rho_{\rm Pl}(d)dd\right)J_N^{-1}.
\tag{L20}
\]

这不是用 RHS 先指定一个算符再宣称它等于原局部核；RHS 是从 (L11) 计算得到的表示。

(L3)—(L5) 的全 N 识别由前述独立生成空间稿提供；下面据此实施[创建审计](stage3-creation-descent-audit.md) C19 原先保留为条件命题的解析向量步骤。结论始终指这个生成空间和 (L11) 的弱定义。

## 7. 完整有限粒子数域上的本质自伴场

令

\[
\mathcal H^{gen}=\bigoplus_{N\ge0}\mathcal H_N^{gen},
\qquad
\mathcal D_{fin}=\bigoplus_{N\ge0}^{alg}\mathcal H_N^{gen}.
\tag{L21}
\]

D_fin 在每个已出现扇区中包含所有 Hilbert 向量，只限制非零扇区数有限。它不是原小种子域 S_N，也不是有限通道紧谱域；本稿不要求那些小域在创建作用下不变。

在 D_fin 上定义 C 把第 N 层送到 C_N^weak 的第 N+1 层，D 把第 N+1 层送到 (C_N^weak)* 的第 N 层，并令

\[
\Phi_{\sigma,\tau_0}^{weak}=C+D\quad\text{于 }\mathcal D_{fin}.
\tag{L22}
\]

因为每个扇区块已在全部输入扇区上有界，D_fin 是不变稠密域；(L19)使 Φ 在该域对称。若 ψ 的物质粒子数支撑不超过 r，则

\[
\|C\psi\|\le\sqrt{(r+1)\mathfrak M}\|\psi\|,
\qquad
\|D\psi\|\le\sqrt{r\mathfrak M}\|\psi\|,
\]
\[
\|\Phi\psi\|\le2\sqrt{\mathfrak M}\sqrt{r+1}\|\psi\|.
\tag{L23}
\]

每次作用至多增加一个粒子，故归纳给

\[
\|\Phi^k\psi\|\le
(2\sqrt{\mathfrak M})^k
\sqrt{\frac{(r+k)!}{r!}}\|\psi\|.
\tag{L24}
\]

对任意 t≥0，Σ_k t^k‖Φ^kψ‖/k! 的控制项相邻比为 2t√𝔐√(r+k+1)/(k+1)→0。因此 D_fin 中每个向量都是整解析向量。由 [Nelson 解析向量定理](https://doi.org/10.2307/1970331)（准确表述亦见[已缓存 Neeb 2011](../sources/neeb-2011-aif2660.pdf) Theorem 4.8），

\[
\boxed{\Phi_{\sigma,\tau_0}^{weak}|_{\mathcal D_{fin}}
\text{ 在 }\mathcal H^{gen}\text{ 上本质自伴。}}
\tag{L25}
\]

所以该明确受调节弱实现具有唯一自伴闭包。这里的唯一性针对这个空间上 (L22) 的稠密对称算符，不是对尚未识别的更大 FP 分布空间作出的唯一性结论。√N 增长也没有被改写为完整场有界。


## 8. 得到的是什么，仍未得到什么

得到的是共同群 cutoff 下的局部模式弱型、其 cutoff 独立极限、FP 生成闭包中的输出、代表独立性、伴随、√N 块界，以及这个生成空间上完整有限粒子数域的本质自伴场。整个链条是

\[
\text{原模式弱核}\ \longrightarrow\
\text{有限群cutoff的运动学配对}\ \longrightarrow\
\text{绝对GA弱极限}\ \longrightarrow\
\text{FP生成闭包输出}\ \longrightarrow\
\text{multiplicity表示}.
\]

尚未额外声称以下两件事：

1. 对每一个已规约的位置不变波函数，直接执行原 (5.4) 的**未调节**时空积分，逐点或以某个更强分布拓扑收敛，并且与本稿弱极限相同。若需要这个更强表述，须另给相应位置测试空间、积分控制及极限交换证明。
2. H_N^gen 已穷尽所有可想象的、满足形式约束且具有有限 FP 积分的不变分布。若物理空间采用更大的定义，仍须证明那些额外候选态在生成闭包中；本稿不在其正交补任意补零，也不在那里宣称唯一场实现。

有界等变性本身并不赋予任意分布上的连续作用。本稿之所以能下降，是额外使用了已建立的 J_N、固定左包 GA 连续界、完整 fiber 物理界和指定的共同弱调节。D_orb 不变性可以不用，但这些条件不能一起略去。
