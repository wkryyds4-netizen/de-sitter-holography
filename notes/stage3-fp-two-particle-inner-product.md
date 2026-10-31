# 两粒子 FP 内积：可证明的等距步骤、群平均字典与受控充分域

后续更新：本稿保留最初充分域证明；后续[共同域定理](stage3-positive-frequency-orbit-core.md)与[完整产品像密度](stage3-product-seed-ga-density.md)已补齐真实种子生成闭包的覆盖，[全N稿](stage3-all-n-generated-space.md)给相应推广。它们不声称覆盖任意另行定义的最大不变分布空间。

日期：2026-10-29。模型为 Kolchmeyer–Liu 的 dS₂、L=1、G_N→0、正实主系列观察者和自由物质。先在 proper SO⁺(2,1) 空间工作，Bose 对称与空间反射分别处理。本稿承接[FP multiplicity 路线](stage3-fp-multiplicity-bridge.md)，不修改旧绑定文件。

**本轮能严格得到三件事：物质自身的 CG 酉变换可直接穿过位置 FP 内积；满足明确混合 L¹/L² 条件的种子，其群平均与位置 FP 内积精确相等；已知 Schur 关系给出的群平均幅度转换恰产生 ρ_Pl 权。尚未证明这类种子在原正频观察者物理空间中的像稠密，故完整 N=2 酉桥没有在这里闭合。**

来源：KL [2411.08090v1](https://arxiv.org/html/2411.08090v1) (4.10)、(4.22)、(4.25)—(4.45)、(4.49)—(4.50) 与附录 E；Marolf–Yu [2610.24073v1](https://arxiv.org/html/2610.24073v1) §3.1.3、§3.2、§3.3，特别 (3.19)、(3.24)、(3.35)—(3.38) 及该节末的积分次序讨论。[本地新文](../sources/marolf-yu-2610.24073v1.html)。群平均的 tempered 正性与 CG reduction 是这些已有工作的结果，不是本项目新定理。以下轨道等距证明是标准 Haar 测度分解及 Bochner–Fubini 定理在本模型的应用。

## 1. 不需要先跨越 FP 难点的物质 CG 分解

令 K₂=Sym²Hν。采用附录 E 的完整归一化 CG 后，物质自身可酉分解为主系列连续谱及离散序列；对于相同 ν 的两个玻色子，内部 principal 只保留 S，内部 discrete 只保留偶数 Δ 的 DΔ⁺、DΔ⁻。这里以归一化后的 dq 和离散计数写成 K₂≅∫⊕K_j dμ_m(j)，并将全部所需多重度包含在 j 中。

原 (4.40) 的物理波函数在每个观察者位置对 y=(x_L,x_R) 取值于 K₂。若 V_m:K₂→∫⊕K_jdμ_m 是这一固定的、只作用于物质的酉映射，定义 Ψ_j(y)=(V_mΨ(y))_j。它不混合观察者位置，故对几乎处处的 y，Parseval 给

\[
\|\Psi(y)\|_{K_2}^{2}=\int\|\Psi_j(y)\|_{K_j}^{2}d\mu_m(j).
\]

在 KL 的六条规约切片上用 Tonelli 积分，即得到允许两边为∞的等式

\[
\boxed{\|\Psi\|_{\rm FP}^{2}=\int d\mu_m(j)\,\|\Psi_j\|_{{\rm FP},j}^{2}.}
\tag{F1}
\]

有限范数情形再由极化得到交叉内积。这一步已成立，不要求对 observer+matter 的非局域 CG 变换与 GF 测度作交换。因此 N=2 的剩余问题可以先降为固定有效物质 irrep j 的三腿 FP Gram 矩阵。

不同 j 的交叉项在 (F1) 中严格正交。相同 principal j 内的 outer S/A 是另一层通道，不能由物质 Parseval 自动消掉。若两个相应位置波函数已属于 FP 的 L² 域，显式 S 偶、A 奇及六切片测度的空间反射不变性给〈S,A〉FP=0；Cauchy–Schwarz 保证这个反射换元合法。该论证只确定交叉项，不证明两条对角范数相等，更不先行证明 A 波函数的 FP 可积性。原脚注33的文字张力沿用[顶角审查](stage3-vertex-analytic-audit.md)的显式约定处理。

## 2. 六切片的实际轨道测度

设 X=dS₂，Y=X×X，G=SO⁺(2,1)。排除两点嵌入内积 c=±1 的退化集合。该集合对普通位置测度为零；本稿不把集中在该集合上的分布态自动纳入 FP Hilbert 空间。

非退化 Y 分成 KL (4.26)—(4.28) 的六个自由群轨道族。每族可选截面

\[
x_L=(t,\theta)=(0,0),\quad
x_R=(\pm\tau,0),(\pm\tau,\pi)\ (\tau>0),
\quad\hbox{或}\quad x_R=(0,\pm\alpha)\ (0<\alpha<\pi).
\]

令 Z 表示这六条切片，y=g·s(z)。有

\[
dy=dg\,d\mu(z),\qquad
d\mu=\sinh\tau\,d\tau\quad\hbox{或}\quad\sin\alpha\,d\alpha.
\tag{F2}
\]

这正是原文的 d(cosh τ) 和绝对值意义的 d(cos α)。常数可直接核，而不是靠δ(0)：取嵌入坐标 X=(sinh t,cosh t cosθ,cosh t sinθ)，群单位元附近参数为旋转 θ_g 与两boost β₁,β₂。对 x_R=(τ,0)，四坐标的Jacobian绝对值是 tanh τ，乘位置体积 cosh τ 得 sinh τ；对 x_R=(0,α)，Jacobian绝对值是 |sin α|。Haar 在单位元的密度取 dθ_g dβ₁dβ₂，与 MY (3.2) 的 dθ_L dθ_R sinhλ dλ 相容。群不变性把该局部计算推广到整个非退化轨道。

## 3. 群平均等于位置 FP 的一个严格充分命题

本节允许任意可分 Hilbert 空间 K 及其酉 G 表示 V，包括 K₂及任意固定物质 irrep。环境运动学空间为 L²(Y;K)，群作用取

\[
[U(k)\psi](y)=V(k)\psi(k^{-1}y).
\]

通过 (F2) 作 untwist：

\[
F_\psi(g,z)=V(g^{-1})\psi(g\cdot s(z)).
\tag{F3}
\]

它是 L²(Y;K)→L²(G×Z;K) 的酉映射；群作用变成 F_{U(k)ψ}(g,z)=Fψ(k⁻¹g,z)。原两观察者正频空间由 (4.10) 等距嵌入这个环境空间，但下面的充分条件仍须在其种子上单独检查。

定义

\[
a_\psi(z)=\int_G\|F_\psi(g,z)\|_Kdg,\qquad
\mathcal D_{\rm orb}=\{\psi\in L^2(Y;K):a_\psi\in L^2(Z,d\mu)\}.
\tag{F4}
\]

**命题。** 若 φ,ψ∈D_orb，则普通群平均 ∫Gdk〈φ,U(k)ψ〉绝对收敛；Bψ(z)=∫G Fψ(g,z)dg 是 Bochner 积分并属于 L²(Z;K)。定义不变量

\[
(\eta\psi)(g\cdot s(z))=V(g)B_\psi(z).
\tag{F5}
\]

则该函数在分布意义下是原种子的群平均，而且

\[
\boxed{\int_Gdk\,\langle\phi,U(k)\psi\rangle
=\int_Zd\mu(z)\langle B_\phi(z),B_\psi(z)\rangle_K
=\langle\eta\phi,\eta\psi\rangle_{\rm FP}.}
\tag{F6}
\]

**证明。** 单模性和 Haar 换元 h=k⁻¹g 给

\[
\begin{aligned}
\int_Gdk\,|\langle\phi,U(k)\psi\rangle|
&\le\int_Zd\mu\int_Gdg\int_Gdk\,
\|F_\phi(g,z)\|\|F_\psi(k^{-1}g,z)\|\\
&=\int_Za_\phi a_\psi\,d\mu
\le\|a_\phi\|_2\|a_\psi\|_2<\infty.
\end{aligned}
\]

故 Fubini 合法，三重积分分别成为两个 Bochner 轨道积分的内积；这给第一等号。第二等号由原 (4.42) 与 (F2)、(F5) 直接得到。对紧支撑位置测试函数重复同一Fubini，也验证 (F5) 与 ∫U(k)ψdk 作为分布相同。证明没有使用无限群体积的相消。

**受控调节。** 取逆元对称的紧支撑可测 cutoff 0≤χ_R≤1，χ_R(g⁻¹)=χ_R(g)，χ_R↑1，例如逐渐扩大的 Cartan 球。令 Bψ,R=∫χ_R(g)Fψ(g,z)dg。逆元对称使它直接等于普通截断群平均 η_Rψ=∫χ_R(k)U(k)ψ dk 在截面 g=e 的值；若不取逆元对称，必须在 Bψ,R 中改用 χ_R(g⁻¹)。由

\[
\|B_{\psi,R}-B_\psi\|_{L^2(Z;K)}
\le\left\|\int_G(1-\chi_R(g))\|F_\psi(g,z)\|dg\right\|_{L^2(Z)}\longrightarrow0
\tag{F7}
\]

取得 FP 范数中的极限；最后一步由 aψ² 的可积支配。该 cutoff 是证明中的辅助轨道调节。有限R的普通截断群平均 η_Rψ 未必群不变；若 ψ 原已在观察者正频闭子空间，则 η_Rψ 仍在该子空间，因为紧支撑L¹群积分保持群不变闭子空间。另一个对象——从切片 Bψ,R 人为作出的不变延拓 V(g)Bψ,R(z)——不保证仍满足原正频谱约束。不能混淆这两个有限R对象；本命题中的物理识别指 (F5) 的极限。

## 4. 群平均侧的 universal ρ 权来自何处

MY (3.24) 是已知的分布式 Schur 关系；其常数由 (3.19) 为

\[
K_P(d)=\frac{4\pi^2}{d\tanh\pi d}=\rho_{\rm Pl}(d)^{-1},
\qquad w(d)=\frac{\pi^3}{K_P(d)}=\pi^3\rho_{\rm Pl}(d).
\tag{F8}
\]

这不是简单把 GA 范数的 K_P 改名为ρ。先把 W_N 的总 principal 部分用 normalized CG 写为 ∫⊕dd(H_d⊗M_{N,d})。再保留左 observer 的独立连续标签 ℓ。对具有明确光滑双谱测试核的种子，令 bψ(d,ξ) 为 ℓ=d 上与单位模对偶张量收缩的系数；在 KL 位置基的相容约定下，对偶张量可写为Σ_n(−1)^n|d,n〉|d,−n〉。Schur 公式给

\[
g_\psi(d,\xi)=K_P(d)b_\psi(d,\xi),\qquad
\langle\phi,\psi\rangle_{\rm GA}
=\int dd\,d\xi\,K_P(d)\overline{b_\phi}b_\psi
=\int dd\,d\xi\,\rho_{\rm Pl}(d)\overline{g_\phi}g_\psi.
\tag{F9}
\]

这里 g 是**群平均之后的不变分布幅度**，b 是种子对角收缩；二者相差 K_P。与原 (4.22)、(4.45) 的校准一致。MY §3.3 已证明相应 tempered CG reduction 的正性，并说明分布测试及积分次序；不能把这些结果归为我们的新发现。

光滑双谱种子的普通 GA 合法性和 (F9) 的完整测试核约定见[谱 GA 桥](stage3-ga-spectral-bridge.md)：有限总K磁指标、远离零端点的联合光滑紧谱核，以及 multiplicity 积分后的谱导数控制，给出绝对收敛。此处绝对收敛指**先完成质量积分**再积 |〈φ,U(g)ψ〉|；不宣称逐质量绝对值与群积分的交换。它证明 GA 谱种子合法，不直接证明 (F4) 的位置轨道强L¹条件；该谱稿的满射构造与无界对角迹例也说明，运动学稠密性本身不能替代物理范数下的稠密性。

对同时满足 (F4) 且其模式分布确由上述光滑谱种子表达的态，(F6)+(F9) 就得到目标 FP 等距式。若这样的不变幅度 g 在全部 multiplicity L²(ρdd) 中稠密，再通过完备化才得到完整酉桥。

## 5. 已覆盖什么，具体缺口还在哪里

环境空间中，untwist后形如 F(g,z)=χ(g)v(z) 的种子（χ∈Cc(G)，v∈L²(Z;K)，可取z紧支撑避开退化端点）明确满足 (F4)，且其有限和在环境L²中稠密。因此充分域非空，调节证明不是形式操作。

但这些环境种子通常不在两个 observer 的正频谱投影像内。**尚未证明**平滑紧质量包、有限原粒子磁指标的物理种子自动属于 D_orb，也未证明 D_orb与物理种子空间交集的群平均像覆盖所需完整通道。质量积分可能把大时间振荡变成任意多项式衰减；仍须控制群K角的退化方向及对z的L²积分，不能只在一条generic boost线上估计后认定成立。

同样不能采用以下捷径：untwist之后正频观察者投影虽是有界等变投影，但它在正则表示的 Plancherel 分解中只定义几乎处处的fibers。非amenable G的trivial表示不在正则Plancherel支撑上；没有连续核/测试域定理，不能把该投影直接“取trivial fiber”，从而声称它保持FP或与群平均交换。

原 (4.50) 的抽象 L² 直和描述也不替代 (F6) 在所需物理稠密域的适用性。当前缺口已缩小为可检验任务：证明足够多的真实观察者谱种子满足 (F4)，或者用比 (F4) 更弱但仍可控制的调节证明位置FP与GA一致，并证明这些态的multiplicity幅度稠密。

## 6. N=2 三腿核的明确验收值

经 (F1) 的物质自身归一化变换后，固定物质 principal j=q 时，outer raw PP→P的正交权为 w(d)，故若完整桥成立，raw系数 f 的位置FP范数应为

\[
\int dd\,da\,\frac{\rho_{\rm Pl}(d)}{w(d)}|f|^2
=\frac1{\pi^3}\int dd\,da\,|f|^2
\]

对S、A各自成立且二者交叉为零。固定物质 discrete j=DΔ± 时，outer raw PD→P的权为

\[
v_\Delta(d)=\frac{d\tanh\pi d}{\Delta-1/2},
\qquad\frac{\rho_{\rm Pl}(d)}{v_\Delta(d)}=\frac{\Delta-1/2}{4\pi^2}.
\tag{F10}
\]

所以其预测 raw FP 对角权为 (Δ−1/2)/(4π²)，不能沿用1/π³。本节是待位置FP验证的明确归一化目标；并未仅据GA计算就把它们标作原 (4.42) 的全通道证明。

Bose投影与群作用对易；采用物质优先树可先在 K₂完成它，再用 (F1)。D⁺与D⁻受空间反射交换，反射投影应在proper群空间上另做；时间反演的反线性实结构再另列。局部创建算符的矩阵元识别、全物理定义域和自伴性仍不由本稿自动推出。