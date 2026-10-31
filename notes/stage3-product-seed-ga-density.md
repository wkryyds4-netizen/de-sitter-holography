# 原粒子产品种子的 GA 幅度与完整 multiplicity 稠密性

日期：2026-10-29。模型为 dS₂、L=1、G→0、ν>0 的两个相同自由玻色子与两主系列 observer，实际 observer 谱下界 s*>0；只处理 proper SO⁺(2,1)。本稿不修改已有文件，不重跑数值计算。

**结果：仅取左 observer 磁指标 0、右 observer 加两物质的总磁指标 0，已经足以产生完整 N=2 主系列 multiplicity 空间的稠密 GA 像。关键不是运动学稠密性自动传到 GA，而是先证明固定左光滑紧质量包时，GA 对任意右总 K=0 Hilbert 向量有明确连续界，再用原粒子产品测试和 CG 酉性证明正交补为零。**

这一 GA 结果不需要先假定产品种子属于位置轨道强域 D_orb。若另一条证明确认这些真实产品种子属于 D_orb，则可与位置 FP–GA 充分等距定理结合，识别其 FP 生成闭包。它仍不独自证明任意另行允许的不变分布都落在该闭包中，也不证明局部创建算符已下降。

输入为 [KL 2411.08090v1 附录 E](https://arxiv.org/html/2411.08090v1#A5)、[MY 2610.24073v1](https://arxiv.org/html/2610.24073v1) (3.2)、(3.8)、(3.24)，以及已经列清全部通道及测度的[normalized CG 树](../literature/stage3-two-particle-cg-map.md)。本文应用这些既有表示论结果，不主张新的群正交定理。

## 1. 为什么总 K=0 已经足够

令

\[
W_2=\left(\int_{a\ge s_*}^{\oplus}da\,P_a\right)
\otimes\operatorname{Sym}^2P_\nu,
\qquad W_2^0=\ker L_{0,\mathrm{total}}.
\]

这里上标 0 指总旋转磁指标，不是物质真空。所有出现的 principal P_d 都恰有一个 K=0 向量；D_Δ⁺ 的磁指标从 Δ≥1 开始，D_Δ⁻ 从 −Δ 开始，所以没有最终总 discrete 的 K=0 向量。完整 CG 分解因此给

\[
W_2^0\simeq\int_{a\ge s_*}^{\oplus}da
\int_{d>0}^{\oplus}dd\,\mathcal M_{d,a}.
\tag{P1}
\]

M_d,a 包含物质优先树中的内部 principal q>0、外层 S/A 两通道，以及内部偶数 Δ 的 D⁺/D⁻。总 K=0 没有最终 discrete，并不删除这些内部 discrete 路径。原分解不出现 complementary 系列，是这里使用的明确模型输入。

对固定 a，取原三因子 K=0 子空间的可数正规正交基 e_i(a)：observer 磁指标 n_R 与两物质磁指标满足 n_R+n₁+n₂=0，物质部分作标准归一化 Bose 对称。每个 e_i 只含有限个原粒子磁指标；例如 n₁≠n₂ 时物质向量为 (|n₁,n₂〉+|n₂,n₁〉)/√2。

固定 a 的 CG 分析酉映射记 U_a，定义

\[
u_i(a,d,\xi)=(U_ae_i(a))(d,\xi),
\qquad
\int_0^\infty dd\,\|u_i(a,d)\|_{\mathcal M_{d,a}}^2=1.
\tag{P2}
\]

U_a 使用合成系数的复共轭；u_i 不是未共轭的合成列。ξ 含 q、通道及离散指标，相应积分/求和压入纤维内积。整族 U_a 可测，由已给显式 CG 系数及谱合成构造取得。

## 2. 一个对全部 principal 参数一致的球函数界

在 KL 正文或附录基中，磁指标 0 的矩阵元不受半相位重定义影响。由 MY (3.8) 或其超几何变换，

\[
D^s_{00}(k_1a_\lambda k_2)
=D^s_{00}(a_\lambda)
=P_{-1/2+is}(\cosh\lambda),\qquad \lambda\ge0.
\]

Legendre 的圆周积分表示（[DLMF 14.12.7](https://dlmf.nist.gov/14.12.E7)，取 m=0、x=coshλ，再将半圆积分写成整圆平均）给

\[
D^s_{00}(a_\lambda)
=\frac1{2\pi}\int_0^{2\pi}
(\cosh\lambda+\sinh\lambda\cos\theta)^{-1/2+is}d\theta.
\]

底数为正，故 s 只给模一相位。因此对每个 s>0，

\[
|D^s_{00}(a_\lambda)|\le
\Xi_0(\lambda):=P_{-1/2}(\cosh\lambda)
\le C(1+\lambda)e^{-\lambda/2}.
\tag{P3}
\]

第二个不等式可直接估计，不需要对 s→0 的渐近做不合法交换。将圆周平移并折叠后，

\[
\Xi_0(\lambda)=\frac2\pi\int_0^{\pi/2}
\frac{du}{\sqrt{e^{-\lambda}+2\sinh\lambda\sin^2u}}.
\]

λ≥1 时令 c₀=1−e⁻²，则 2sinhλ≥c₀e^λ、sin u≥2u/π，所以

\[
\Xi_0(\lambda)
\le \frac{e^{-\lambda/2}}{\sqrt{c_0}}
\operatorname{arsinh}(\sqrt{c_0}e^\lambda)
\le C(1+\lambda)e^{-\lambda/2}.
\]

0≤λ≤1 由连续性补上统一常数。这里只需要球函数 n=0；无须为所有磁指标另建一致界。

## 3. 左质量先积分，才交换右总谱

取 f₁,f₂∈C∞c((s*,∞))，左 observer 波包为 |f_j,0〉。它们的左矩阵系数是

\[
L_{12}(g)=\langle f_1,0|U_L(g)|f_2,0\rangle
=\int d\ell\,\overline{f_1(\ell)}f_2(\ell)D^\ell_{00}(g).
\tag{P4}
\]

质量积分在这里已经完成。紧谱远离 0，MY (3.8) 的连接公式给

\[
D^\ell_{00}(a_\lambda)
=e^{-\lambda/2}\sum_{\epsilon=\pm}
 A_\epsilon(\ell)e^{i\epsilon\ell\lambda}
+O(e^{-3\lambda/2}),
\]

其中 A_ε 在该紧集上光滑，余项可允许固定次 λ 多项式。对 leading 项按 ℓ 分部积分 r 次，所有边界项为零，得到

\[
|L_{12}(a_\lambda)|
\le C_r e^{-\lambda/2}(1+\lambda)^{-r}
+C'(1+\lambda)^p e^{-3\lambda/2}.
\tag{P5}
\]

MY Haar 是 dg=dθ₁dθ₂ sinhλ dλ。结合 (P3)，r≥3 时有

\[
J_{12}:=\int_G|L_{12}(g)|\Xi_0(\lambda(g))dg<\infty.
\tag{P6}
\]

例如 leading 径向上界为常数乘 (1+λ)^(1−r)。这不是宣称固定左质量的矩阵元与 Ξ₀ 之积群可积；左光滑质量积分不可撤回。

现在令 v₁,v₂ 是任意 W₂⁰ Hilbert 向量，不要求其总 d 或内部 ξ 支撑紧，也不要求有限原磁指标。按 (P1) 把其系数记为 v_j(d)∈M_d，其中 M_d=∫⊕da M_d,a。由于两向量只有总 K=0 分量，

\[
\langle v_1,U_W(g)v_2\rangle
=\int_0^\infty dd\,D^d_{00}(g)R_{12}(d),
\quad
R_{12}(d)=\langle v_1(d),v_2(d)\rangle_{\mathcal M_d},
\]
\[
\int_0^\infty |R_{12}(d)|dd\le\|v_1\|\|v_2\|.
\tag{P7}
\]

(P3)、(P6)、(P7)给真正的绝对 Fubini 控制：

\[
\int_Gdg\int_0^\infty dd\,
|L_{12}(g)R_{12}(d)D^d_{00}(g)|
\le J_{12}\|v_1\|\|v_2\|<\infty.
\tag{P8}
\]

所以左波包乘任意右 K=0 Hilbert 向量的群平均交叉配对已绝对收敛，而且现在可以交换右 d 积分与群积分。这一估计避开了对非紧右谱密度直接形式套用 Plancherel 的识别缺口。

## 4. 波包版 Schur 恒等式及真实产品幅度

定义

\[
I_{12}(d)=\int_G L_{12}(g)D^d_{00}(g)dg.
\]

(P6)保证普通积分存在，且其统一支配与 D_d 的连续性使 I₁₂(d)在 d>0 连续。对任意光滑紧谱测试函数 q(d)，可由 (P6)交换 q 的积分，再使用 MY (3.24) 的分布 Schur 关系。磁指标为 0 时 D⁽ᵈ⁾₀₀(g⁻¹)=D⁽ᵈ⁾₀₀(g)，没有对偶符号或半相位残留。得到作为分布、从而作为连续函数的等式

\[
I_{12}(d)=K_P(d)\overline{f_1(d)}f_2(d),
\qquad
K_P(d)=\frac{4\pi^2}{d\tanh\pi d}
=\rho_{\rm Pl}(d)^{-1}.
\tag{P9}
\]

这里使用已有 Schur 分布恒等式；没有声称固定两个质量的群积分绝对收敛。

对 S_j=|f_j,0〉⊗v_j，(P8)—(P9)给

\[
Q_{GA}(S_1,S_2)
=\int_0^\infty dd\,K_P(d)\overline{f_1(d)}f_2(d)
 \langle v_1(d),v_2(d)\rangle.
\tag{P10}
\]

所以在已经建立的 GA 不变幅度坐标中，

\[
g_{f,v}(d)=K_P(d)f(d)v(d).
\tag{P11}
\]

这与[光滑谱 GA 桥](stage3-ga-spectral-bridge.md)的 g=K_Pb 完全相同，并非重新定义群平均。具体兼容性可这样取得：先选 χ∈C∞c((s*,∞))，使 χ 在 supp f 上恒为1；若同时比较多个左包，则在它们支撑的并集上恒为1。按已独立证明的 (P10)，v 与 χ(d)v 有相同 GA 像，因为二者之差的 GA 范数为零。然后仅把 χv 在归一化谱表示中用光滑紧 d、有限内部通道测试函数逼近，所有 d 支撑均取在 (s*,∞) 内；这些向量属于既有谱测试域，幅度由原 G8 给出。由 (P10)和下一节界，逼近在 GA 范数中收敛到 (P11)，且普通群积分也由 (P8)对固定左包连续。这里无需、也没有在该受限测试域中作全部 v 的运动学逼近；原 v 可能含 d<s* 的总谱。非紧右总谱因而没有引入另一个延拓。

取右原粒子产品种子

\[
v_i^k(a)=k(a)e_i(a),\qquad k\in C_c^\infty((s_*,\infty)),
\]

其精确 GA 幅度为

\[
\boxed{g_{i;f,k}(d,a,\xi)
=K_P(d)f(d)k(a)u_i(a,d,\xi).}
\tag{P12}
\]

(C∞紧左包)×(C∞紧右包)×(有限原粒子磁指标)正是这里的真实产品种子；Bose 基向量只有有限项。u_i 的全部内部谱尾由 CG 酉范数控制，没有截掉任何 A 或内部离散路径。

## 5. 固定左包的连续界与无界迹反例不冲突

令

\[
B_f^2=\operatorname*{ess\,sup}_{d>0}K_P(d)|f(d)|^2<\infty.
\]

由 f 的紧支撑远离谱端点，常数有限；(P10)给

\[
\|g_{f,v}\|_{L^2(\rho_{\rm Pl}dd;\mathcal M_d)}^2
=\int dd\,K_P(d)|f(d)|^2\|v(d)\|^2
\le B_f^2\|v\|_{W_2}^2.
\tag{P13}
\]

更一般地，交叉配对不超过 ‖K_P overline(f₁)f₂‖∞‖v₁‖‖v₂‖。这只是**固定左波包**后的右 Hilbert 连续性。常数依赖 f 的加权上确界，不能仅用 ‖f‖₂统一控制，故不与先前对角迹的全运动学无界反例冲突。

该界是将 CG 谱逼近传到 GA 的实际估计；本稿没有使用“运动学稠密，所以 GA 稠密”的错误推理。

## 6. 产品 GA 像的正交补严格为零

目标空间记为

\[
\mathcal H_{spec}
=\int_{d\ge s_*}^{\oplus}\rho_{\rm Pl}(d)dd
\int_{a\ge s_*}^{\oplus}da\,\mathcal M_{d,a}.
\]

设 h∈H_spec 与所有 (P12) 正交。因为 ρ_Pl K_P=1，对全部 i 和 f,k∈C∞c((s*,∞))，

\[
0=\langle h,g_{i;f,k}\rangle_{spec}
=\int dd\,da\,f(d)k(a)
 \langle h(d,a),u_i(a,d)\rangle_{\mathcal M_{d,a}}.
\tag{P14}
\]

这里是合法的局部 L¹ 配对。确切地，对紧区间 A,D⊂(s*,∞)，Cauchy–Schwarz 与 (P2)给

\[
\begin{aligned}
\int_{A\times D}|\langle h,u_i\rangle|\,da\,dd
&\le\left(\int_{A\times D}\rho_{\rm Pl}(d)\|h\|^2\,da\,dd\right)^{1/2}
\left(\int_{A\times D}\rho_{\rm Pl}(d)^{-1}\|u_i\|^2\,da\,dd\right)^{1/2}\\
&\le\|h\|_{spec}\sqrt{|A|\,\sup_{d\in D}\rho_{\rm Pl}(d)^{-1}}<\infty.
\end{aligned}
\tag{P15}
\]

光滑乘积测试的有限线性组合在局部测试函数空间中稠密，所以 (P14)推出

\[
\langle h(d,a),u_i(a,d)\rangle=0
\quad\text{对几乎处处 }(a,d),\text{对全部 }i.
\tag{P16}
\]

i 可数，能取共同零测集。最后一步不需要在单个 d 上贸然把广义 CG 列称为正交基。固定紧 A,D，令

\[
w_{A,D}(a,d)=\mathbf1_A(a)\mathbf1_D(d)h(d,a),
\]

并在 d<s* 处零延拓。由于 ρ_Pl 在 D 上有正下界，w_A,D 属于裸 da dd 的右 K=0 谱 Hilbert 空间。对几乎处处 a∈A，逆 CG 向量的原基系数为

\[
\langle e_i(a),U_a^{-1}w_{A,D}(a)\rangle
=\int_D dd\,\langle u_i(a,d),h(d,a)\rangle=0.
\tag{P17}
\]

等式可由固定 a 的 L² Cauchy–Schwarz解释；Fubini与可数性给共同的 a 零测集。e_i(a)在原 K=0 子空间是完整正规正交基，故 U_a⁻¹w_A,D(a)=0，进而 w_A,D=0。取 A,D 的可数紧区间穷竭，得到 h=0。

因此所有真实产品幅度 (P12) 的有限线性包在 H_spec 中稠密。证明使用的是：明确幅度公式、局部 L¹ 控制、产品测试的总性和 CG 完整酉变换；不经过无界对角取值映射的连续性假设。

## 7. 与位置 FP 强域结果怎样连接

[位置 FP 充分定理](stage3-fp-two-particle-inner-product.md)的假设是实际位置种子 ψ 满足

\[
a_\psi(z)=\int_G\|V(g^{-1})\psi(g\cdot s(z))\|dg\in L^2(Z,d\mu),
\]

而不只是标量矩阵系数的 L¹(G) 性。本文的 (P8)并未证明这个更强的轨道条件。

若独立证明确认所有本稿使用的产品种子 S_i=|f,0〉⊗k(a)e_i(a)都属于该 D_orb，并属于原正频 observer 子空间，则位置稿给它们及其有限线性组合的

\[
\langle\eta S_1,\eta S_2\rangle_{FP}
=Q_{GA}(S_1,S_2)=\langle g_{S_1},g_{S_2}\rangle_{spec}.
\tag{P18}
\]

结合已证的 GA 幅度稠密性，得到准确的完成化结论

\[
\boxed{\mathcal H_{spec}\ \simeq\
\overline{\operatorname{span}\{\eta S_i\}}^{\,FP}.}
\tag{P19}
\]

即真实产品种子生成的 FP 闭包与完整 normalized multiplicity 空间等距同构；不存在某个遗漏的 S/A 或内部 discrete multiplicity 正交补。

(P19)并不单独排除一种更大定义：若把所有满足某些形式约束的任意不变分布都纳入“完整位置 FP 空间”，还要证明它们可由上述种子逼近，或先明确物理 Hilbert 空间按该种子规约闭包定义。不能从两个空间抽象同维、或群正则表示的测度零 trivial fiber 取值推出这一覆盖性。

原 (5.4) 的局部模式创建与乘法空间 fiber 算符的弱矩阵元兼容，仍须按[创建下降审计](stage3-creation-descent-audit.md) C17 独立证明。全 N 归纳和本质自伴性也不由本稿 N=2 密度结果自动完成。

本轮新增的是对非紧右总谱仍有效的 (P8)—(P13)，以及无需无界迹连续性的 (P14)—(P17) 稠密性证明。位置强域、任意不变分布覆盖和局部场识别各自保持明确边界。
