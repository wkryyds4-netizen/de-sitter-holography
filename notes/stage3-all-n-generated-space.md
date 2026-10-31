# 全有限粒子数的正频种子生成空间

2026-10-29。模型为dS₂、L=1、G→0、固定ν>0自由实标量、两主系列observer与同一正谱下限s_*>0；群为proper SO⁺(2,1)。本文将[两粒子产品GA密度](stage3-product-seed-ga-density.md)的构造推广到每个有限物质粒子数，并明确物理空间采用**真实正频产品种子经群规约后，在原位置FP范数中的闭包**。

所证明的是这个具体生成空间的全通道识别。若另行把所有形式不变分布都定义成更大的物理空间，仍需比较该最大定义与这里的闭包；不能仅凭下述结果宣布两者等同。原局部场的受调节弱作用另见[弱下降稿](stage3-regulated-local-descent.md)。

## 1. 各粒子数的运动学表示

记

\[
H_R=\int_{a\ge s_*}^{\oplus}da\,P_a,
\quad K_N=\operatorname{Sym}^N P_\nu,
\quad W_N=H_R\otimes K_N.
\tag{N1}
\]

K₀=ℂ，表示物质真空；右总旋转磁指标0的子空间另记W_N^{(0)}，避免将两个0混淆。H_L同样是ℓ≥s_*的正频主系列直积分。

需要的既有表示论输入为KL附录E中完整酉分解：P×P的两重连续P与D⁺/D⁻，以及P×D^±的一重连续P与同号离散系列。normalized权和相位见[两粒子CG树](../literature/stage3-two-particle-cg-map.md)。特别不把某条局部创建顶角只含S误当作全树只含S。

## 2. 对任意N≥1的分解存在性

先固定a，不作物质Bose投影。N=1时，P_a⊗P_ν的KL完整分解直接给Lebesgue总principal谱dd及正、负离散支。

若N时的分解已由principal及discrete组成，继续张量一条P_ν：每条principal支使用PP定理；每条discrete支交换张量次序后使用PD定理。先在各纤维作完整normalized酉变换，再作Hilbert直积分，得到N+1时的完整分解。principal输出均可选裸dd；所有此前连续标签、通道和离散路径转入multiplicity空间。有限N下空间可分，涉及可数路径与标准连续谱；Fubini对Hilbert范数的非负积分合法。该归纳不产生complementary系列。

随后施加物质对称化P_Sym。它是与群作用对易的正交投影，因此在不重复参数化的不可约直积分上形如1⊗p_N(a,d)，把multiplicity纤维换成其闭投影像。它可以混合树路径，但不会产生新的表示谱。最终可写

\[
W_{N,a}^{P}\simeq_G\int_{d>0}^{\oplus}dd\,
 [P_d\otimes\mathcal M_{N,d,a}],\qquad N\ge1,
\tag{N2}
\]

并另保留全部最终discrete支。这一步只断言存在一个完整normalized酉分解，不声称已列出任意N的闭式Racah核。

每个P_d含唯一单位K=0向量，而整数离散系列没有K=0。因此限制总磁指标0后，N2给整个该子空间的酉映射

\[
U_{N,a}:W_{N,a}^{(0)}\longrightarrow
\int_{d>0}^{\oplus}dd\,\mathcal M_{N,d,a}.
\tag{N3}
\]

内部discrete路径仍在M中。原磁指标基使a上的Hilbert场平凡；KL显式CG及其可去退化极限给联合可测分析变换。有限步复合保持可测；Bose投影变换后为可测可分解投影，其纤维像构成可测Hilbert场。由此可选联合可测的U_N,a并在共同(a,d)零集意义下作a直积分，记M_N,d=∫_{a≥s_*}⊕da M_N,d,a。

## 3. 小产品域及其位置FP意义

对每个N≥1，选W_N,a^{(0)}的原粒子Bose正交基e_{N,i}(a)，i可数，且每个向量只含有限原磁指标。定义小种子域S_N为以下向量的有限线性包：

\[
S_{i;f,k}=|f,0\rangle_L\otimes k(a)e_{N,i}(a),
\qquad f,k\in C_c^\infty((s_*,\infty)).
\tag{N4}
\]

它在完整运动学空间中并不稠密，也不需要这个性质。原粒子基的右observer磁指标固定；余下有限Bose态是一个物质Hilbert向量。因此N4正是[真实正频共同域定理](stage3-positive-frequency-orbit-core.md)覆盖的产品种子，或其有限和。该定理允许任意酉物质K，故对每个K_N原样适用。

其轨道Bochner积分ηS满足位置FP＝GA。定义

\[
\mathcal H_N^{gen}:=
\overline{\{\eta S:S\in\mathcal S_N\}}^{L^2(Z,d\mu;K_N)}.
\tag{N5}
\]

这里六切片测度与原KL(4.42)一致；每个元素是原不变向量值函数序列的FP极限。先取零范数商再完备化给相同空间。

## 4. 同一估计给全部N的GA幅度与密度

N3中定义u_N,i(a,d)=U_N,a e_N,i(a)，于是∫dd‖u_N,i(a,d)‖²=1。对于任意右总K=0向量v∈W_N^{(0)}，其群矩阵系数仍为

\[
\langle v_1,U_N(g)v_2\rangle
=\int dd\,D^d_{00}(g)\langle v_1(d),v_2(d)\rangle_{M_{N,d}}.
\tag{N6}
\]

[产品GA稿](stage3-product-seed-ga-density.md)的统一球函数界和左质量先积分估计只使用这一形式及∫|〈v₁(d),v₂(d)〉|dd≤‖v₁‖‖v₂‖，不使用物质数为2。故同一绝对Fubini与波包Schur证明给

\[
Q_{GA}(f_1\otimes v_1,f_2\otimes v_2)
=\int dd\,K_P(d)\overline{f_1(d)}f_2(d)
 \langle v_1(d),v_2(d)\rangle,
\]

\[
A_N(f\otimes v)(d)=K_P(d)f(d)v(d),
\quad\|A_N(f\otimes v)\|^2_{spec}
\le B_f^2\|v\|^2,
\quad B_f^2=\operatorname*{ess\,sup}K_P|f|^2.
\tag{N7}
\]

N≥1的目标为

\[
\mathcal H_N^{spec}=
\int_{d\ge s_*}^{\oplus}\rho_{Pl}(d)dd\,
\mathcal M_{N,d},\qquad \rho_{Pl}=K_P^{-1}.
\tag{N8}
\]

要证S_N的幅度像稠密，设h与全部K_P f(d)k(a)u_N,i(a,d)正交。ρK_P=1给所有光滑乘积测试的配对为零；在紧a、d区间，局部L¹界为

\[
\int_{A\times D}|\langle h,u_{N,i}\rangle|da\,dd
\le\|h\|_{spec}\sqrt{|A|\sup_D\rho_{Pl}^{-1}}.
\]

所以〈h,u_N,i〉=0几乎处处，对可数i同时成立。将h在紧A×D截断、在d<s_*处零延拓后，得到裸da dd的L²向量。对几乎处处a作U_N,a逆变换，全部原基e_N,i系数为零，故该向量为零。可数紧区间穷竭给h=0。这与两粒子稿P14—P17相同，既不依赖无界对角迹的连续性，也不忽略内部通道。

## 5. N=0必须单独处理

N=0时固定a的右表示就是P_a，总标签d=a；不能错误地写成固定a下的Lebesgue dd谱。正确的是

\[
W_0^{(0)}\simeq L^2([s_*,\infty),dd),
\qquad \mathcal M_{0,d}=\mathbb C.
\tag{N9}
\]

取S₀=|f,0〉_L⊗|k,0〉_R，仍属于真实轨道共同域。N7变为

\[
g_{f,k}(d)=K_P(d)f(d)k(d),\qquad
\|\eta S_0\|_{FP}^2=\int K_P(d)|f(d)k(d)|^2dd.
\tag{N10}
\]

给定任意g∈Cc∞((s_*,∞))，选f∈Cc∞在supp g上恒1，并取k=ρ_Pl g，则g_f,k=g。故像在L²(ρ_Pl dd)稠密。这里不需要N≥1的固定a CG基总性论证。

## 6. 各N的明确等距识别

对所有有限N，N4或N10的真实种子同时有位置FP等距与谱GA范数。由已证的像稠密性，唯一延拓得到满射等距

\[
\boxed{J_N:\mathcal H_N^{spec}\longrightarrow\mathcal H_N^{gen},
\qquad J_N A_N(S)=\eta S.}
\tag{N11}
\]

N5中的闭包与N8的完整multiplicity匹配；这是一项带具体映射的等距识别，不是抽象同维论证。更大产品种子域的比较使用[有限K延拓稿](stage3-finite-k-ga-extension.md)的同一论证。该稿显式写N=2；推广到这里每个N时，N2已保证仍只有principal与discrete，固定总磁指标仍只容许有限个最终离散标签。故统一Ξ界、左质量先积分、跨系列为零和正文对偶相位完全不变。总旋转投影仅作用于K_N物质向量，仍保持真实产品强域。任意有限左磁指标产品种子的幅度为g=K_PΣ_j(−1)^{m_j}f_j v^P_{j,−m_j}；其每态范数及与小域的交叉配对都是相同谱内积。由S_N幅度像稠密，FP距离等于谱距离而趋于零，故全部这些真实产品的FP闭包也等于H_N^gen。N=0时右侧没有discrete，使用N9的谱分解同样成立。这没有声称任意非产品运动学向量都可未经域检验群平均。

同样，N11不把最大形式分布空间纳入结论。后续可用N11和原模式有界等变创建块构造[规约弱场](stage3-regulated-local-descent.md)，并检验物理伴随；该操作的定义域与调节顺序须明确写出。有限G、微观全息字典、反射与反线性时间反演的进一步约束也不由此得到。
