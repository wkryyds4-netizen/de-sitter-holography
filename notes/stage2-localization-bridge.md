# 第二阶段字典 A：由静态管测试函数逼近世界线场

日期：2026-10-29。模型固定为 dS₄ 上实、自由、无质量、共形耦合标量，Bunch–Davies（BD）准自由态，半径 L>0，中心测地观察者，β=2πL。本文证明裸世界线时间涂抹场的局域代数归属，并给出供共同钟导数接口使用的范数控制。结论不包含有限 G、相互作用、反冲、复合场，亦不声称一个世界线场生成整个静态补丁代数。

**结论。** 对实 C∞c 时间窗口，以及本文所用实 cos² 紧支撑窗口和 Gaussian 窗口，存在静态补丁内的实 C∞c 时空测试函数序列 j_n，使其 BD 单粒子向量收敛到该世界线场的向量 v_f。相应 Weyl 算符在 BD 表示中强收敛，因此 W(λv_f) 属于静态补丁 von Neumann 代数，所有实 λ 均成立；Segal 场 A(v_f) 的自伴闭包附属于该代数。对紧支撑窗口，可把所有充分后期的测试函数放进包含其世界线支撑段的任意给定开管。Gaussian 的归属结论使用包含整条世界线的开管或整个静态补丁，不把它说成有限时间局域算符。

## 1 原始依据与本证明的位置

- [Ali–Bhattacharya–Lochan 2003.11046v2](https://arxiv.org/html/2003.11046v2)，式 (1)、(2) 给平坦共形坐标，式 (8) 及其后的 y 定义、§3.1 式 (13) 给共形 BD 完整两点函数。这里保留空间变量，不仅使用已经限制到观察者世界线的式子。
- [Pinamonti 0806.0803v2](https://arxiv.org/pdf/0806.0803v2)，2009-02-03，§2 Lemma 2.1 的共形 Klein–Gordon 变换和 §3.1 Proposition 3.2：四维该基本场的共形权为 1，协变测试函数对应权为 3。下文仍显式计算体积因子，避免把两种权混用。
- [Witten 2303.02837v2](https://arxiv.org/html/2303.02837v2)，§2—3 讨论时间涂抹场及其有界函数产生的代数。本证明给出当前自由模型的具体近似，不以一般性讨论代替极限计算。

原论文提供 Wightman 函数和共形协变输入；下面的局域化构造、误差公式与序列选择是据此写出的直接证明，未作首创性声明。这里不需要把任意正谱函数猜测成一个局域场论。

## 2 坐标、完整两点函数和测试函数权

平坦共形图中的度规为

\[
ds^2=a(\eta)^2(-d\eta^2+d\mathbf x^2),\qquad
\eta<0,\quad a(\eta)=-\frac L\eta.
\tag{1}
\]

中心曲线 γ(t)=(η(t),0)，其中 t 为物理固有时，

\[
\eta(t)=-Le^{-t/L},\qquad
 t(\eta)=-L\log(-\eta/L),\qquad dt=a(\eta)d\eta.
\tag{2}
\]

该中心观察者的静态补丁在此图中是

\[
D=\{\eta<0,\ |\mathbf x|<-\eta\}.
\tag{3}
\]

例如物理面积半径 r=a|x|，条件 (3) 等价于 r<L。中心静态时间平移 δ 对应 (η,x)↦e^(−δ/L)(η,x)，其 Killing 场为 −(η∂η+x·∂x)/L，在 D 中类时，且在 x=0 上即 ∂t。

完整共形 BD 两点函数可写为边界值分布

\[
W_{\rm dS}(X,X')
=\frac{a(\eta)^{-1}a(\eta')^{-1}}
 {4\pi^2[|\mathbf x-\mathbf x'|^2-(\eta-\eta'-i0)^2]}
=a(\eta)^{-1}a(\eta')^{-1}W_M(X,X').
\tag{4}
\]

它等价于原文的 H²/(4π²y)，H=1/L。不能把 i0 丢掉后作普通奇异函数积分。因场是自由准自由态，(4) 同时确定所有 Wick 关联函数；在该坐标域的 Weyl 代数上给出 BD 态与共形 Minkowski 真空实现之间的 GNS 识别。

令 φdS(j)=∫dvol_dS jφdS，采用实标量测试函数 j。由 dvol_dS=a⁴dηd³x 和 φdS=a⁻¹φM，

\[
\phi_{\rm dS}(j)=\phi_M(a^3j).
\tag{5}
\]

因此如果需要 Minkowski 坐标测试**密度** h(η,x)，de Sitter 协变测试**函数**是 j=a⁻³h。这里的 a⁻³ 不能改成 a⁻¹ 或 a⁻⁴。

## 3 世界线时间窗在共形坐标中的准确系数

先用物理时间窗口 F(t)，目标为

\[
A_F=\int_{\mathbb R}dt\,F(t)\phi_{\rm dS}(\eta(t),0).
\tag{6}
\]

由于 dt=a dη 与场的 a⁻¹ 恰好抵消，

\[
A_F=\int_{-\infty}^0d\eta\,g(\eta)\phi_M(\eta,0),
\qquad \boxed{g(\eta)=F[-L\log(-\eta/L)]}.
\tag{7}
\]

本项目模窗口的定义为 A_f(s₀)=∫du f(u)φdS(β(s₀+u),0)，所以

\[
F(t)=\frac1\beta f(t/\beta-s_0),\qquad
\boxed{g(\eta)=\frac1\beta
 f\!\left[-\frac{L}{\beta}\log(-\eta/L)-s_0\right]}.
\tag{8}
\]

此处 1/β 必须保留。若 f 无量纲且归一化，F 和 g 的量纲为长度⁻¹；j 的量纲为长度⁻⁴，(5)—(8) 得到的 A 仍具有场的长度⁻¹ 量纲。

若 F∈C∞c(R)，其支撑映到某个 I=[η₁,η₂]⊂(−∞,0)，g 零延拓后属于 C∞c(R)。这是下面空间逼近最先证明的情形。

## 4 显式静态管近似与完整交叉范数

选 χ∈C∞c(R³)，实、非负、球对称，∫χd³x=1，suppχ⊂{|x|≤Rχ}。定义 χε(x)=ε⁻³χ(x/ε)，以及

\[
\boxed{j_\varepsilon(\eta,\mathbf x)
 =a(\eta)^{-3}g(\eta)\chi_\varepsilon(\mathbf x)}.
\tag{9}
\]

若 εRχ<−η₂，则 jε∈C∞c(D)。若另指定开管 O⊂D，且 O 包含紧集 {(η,0):η∈suppg}，由紧性可进一步缩小 ε，使 suppjε⊂O。这不是只在每个时间点分别选宽度后忽略整体支撑，而是利用整个紧支撑段与 O 的补集之间的正距离。

用下列 Fourier 约定：

\[
\widehat g(k)=\int d\eta\,g(\eta)e^{ik\eta},\qquad
\widehat\chi(\mathbf q)=\int d^3x\,\chi(\mathbf x)e^{-i\mathbf q\cdot\mathbf x}.
\tag{10}
\]

Minkowski 模展开取 [a(k),a†(k')]=δ³(k−k') 及系数 [(2π)³2|k|]⁻¹ᐟ²。作用于真空的创建向量为

\[
v_\varepsilon(\mathbf k)
=\frac{\widehat g(|\mathbf k|)\widehat\chi(\varepsilon\mathbf k)}
 {(2\pi)^{3/2}\sqrt{2|\mathbf k|}},\qquad
v_0(\mathbf k)
=\frac{\widehat g(|\mathbf k|)}
 {(2\pi)^{3/2}\sqrt{2|\mathbf k|}}.
\tag{11}
\]

向量空间采用 L²(R³,d³k)。由于 χ 球对称，下式的 χ̂(εk) 表示其径向值。完整 Wightman 函数 (4) 对两份时空测试函数配对，逐一消去各自的 a⁴、a⁻³、a⁻¹，给出

\[
\begin{aligned}
\langle v_\varepsilon,v_\delta\rangle
 &=\iint dvol_Xdvol_{X'}\,
   j_\varepsilon(X)W_{\rm dS}(X,X')j_\delta(X')\\
 &=\int\frac{d^3k}{(2\pi)^3\,2|\mathbf k|}
   |\widehat g(|\mathbf k|)|^2
   \overline{\widehat\chi(\varepsilon\mathbf k)}
   \widehat\chi(\delta\mathbf k).
\end{aligned}
\tag{12}
\]

因此本证明确实控制两份逼近的交叉项，不只是比较方差 Z。直接得到

\[
\|v_0\|^2=\frac1{4\pi^2}\int_0^\infty
 k|\widehat g(k)|^2dk,
\tag{13}
\]

\[
\boxed{\|v_\varepsilon-v_0\|^2
=\frac1{4\pi^2}\int_0^\infty k|\widehat g(k)|^2
 |\widehat\chi(\varepsilon k)-1|^2dk\ \longrightarrow\ 0.}
\tag{14}
\]

证明：|χ̂|≤1、χ̂(0)=1，所以误差因子≤4 并逐点趋零；g∈C∞c 使 ĝ 快速下降，k=0 端的 k|ĝ(k)|² 也可积。支配收敛定理适用。这个空间极限与 i0→0 的边界值已通过正频模表示正确处理，不要求在光锥奇点处交换普通绝对积分。

还有可选的显式速率。令 m₂=∫|x|²χ(x)d³x。球对称性使一阶矩为零，|χ̂(q)−1|≤m₂|q|²/2，故

\[
\|v_\varepsilon-v_0\|^2
\le \frac{m_2^2\varepsilon^4}{16\pi^2}
\int_0^\infty k^5|\widehat g(k)|^2dk.
\tag{15}
\]

紧支撑光滑 g 保证右端有限。该常数依赖具体时间窗口；Gaussian 的截断长度改变时，不把它当作统一常数。

## 5 与已用静态模谱的核对

把 (4) 限制到中心线并用 (2)，得到

\[
W_t(t-t')=-\frac1{16\pi^2L^2
 \sinh^2[(t-t'-i0)/(2L)]}.
\tag{16}
\]

它的物理频率谱是 r(ω)=ω/[4π²(1−exp(−βω))]，W_t(t)=∫r(ω)e^(−iωt)dω。于是对实光滑紧支撑 F，

\[
\begin{aligned}
\|v_F\|^2
 &=\frac1{4\pi^2}\int_0^\infty k|\widehat g(k)|^2dk\\
 &=\int_{\mathbb R}r(\omega)|\widehat F(\omega)|^2d\omega\\
 &=\frac1{4\pi^2}\int_0^\infty
 \omega\coth(\beta\omega/2)|\widehat F(\omega)|^2d\omega.
\end{aligned}
\tag{17}
\]

这三个式子是同一完整二点分布在不同变量下的配对。第一行中的 k≥0 是 Minkowski 共形时间 η 的正频率；第二行中的 ω 是静态固有时 t 的双侧频率，两者没有被等同。

对本项目模时间窗口 f，可写成

\[
\boxed{\|v_f\|^2
=\frac1{4\pi^2\beta^2}\int_0^\infty
 E\coth(E/2)|\widehat f(E)|^2dE.}
\tag{18}
\]

平移 s₀ 不改变该范数。因 f 实，|f̂|² 为偶函数；(18) 相当于全线谱 ρ(E)=E/[4π²β²(1−e^(−E))] 的配对。对 E≥0，

max(E,2)≤E coth(E/2)≤E+2。

故在固定 β 下，实窗口的该单粒子范数与 H^(1/2)(R) 范数等价（Fourier 权重为 sqrt(1+E²)）。这为下面两个非 C∞c 窗口提供具体的逼近拓扑。

## 6 cos² 窗口：先在谱范数中光滑化

设

\[
f_T(u)=T^{-1}\cos^2(\pi u/(2T))\,\mathbf1_{|u|\le T}.
\tag{19}
\]

其零延拓为 C¹、H²，故属于 H^(1/2)；既有 Fourier 式 f̂T(E)=π²sin(ET)/[ET(π²−E²T²)] 也直接显示 (18) 有限。

选实、非负、偶的一维 m∈C∞c，∫m=1，令 fδ=mδ*fT。则 fδ∈C∞c，实、归一化，支撑仅比 [−T,T] 扩大 O(δ)，且

\[
\|v_{f_\delta}-v_{f_T}\|^2
=\int_{\mathbb R}\rho(E)|\widehat f_T(E)|^2
 |\widehat m(\delta E)-1|^2dE\ \longrightarrow\ 0.
\tag{20}
\]

这里可先在光滑窗口差的范数中构造完备极限；右端就是该极限的范数，避免预先假设非光滑点场已存在。|m̂|≤1 给出支配收敛与 ||v_fδ||≤||v_fT||。其极限的二点和高阶 Wick 矩与原来用实际谱定义的窗口完全一致。

对每个 δ 先用 (8)—(14) 构造 jδ,ε。指定 δn→0，再选 εn 使

\[
\|v_{j_{\delta_n,\varepsilon_n}}-v_{f_{\delta_n}}\|<1/n,
\quad \operatorname{supp}j_{\delta_n,\varepsilon_n}\subset O.
\tag{21}
\]

只要 O 包含原窗口的紧世界线支撑段，它也包含充分小 δ 所扩展的支撑段。三角不等式结合 (20)、(21) 给 v_jn→v_fT。这些 jn 都是实 C∞c 时空测试函数；并没有把 C¹ 窗口直接冒充测试函数空间中的光滑元素。

## 7 Gaussian：先截断时间，再对角选择空间宽度

对 fσ(u)=exp[−u²/(2σ²)]/(sqrt(2π)σ)，选 c∈C∞c(R)，0≤c≤1，c=1 于 [−1,1]，支撑于 [−2,2]，并取实偶 c。令

\[
\widetilde f_R(u)=c(u/R)f_\sigma(u),\qquad
f_R=\widetilde f_R/\int\widetilde f_R.
\tag{22}
\]

归一化分母趋 1。未归一化差 hR=(1−c(u/R))fσ 满足

hR'=(1−c(u/R))fσ'−R⁻¹c'(u/R)fσ。

Gaussian 及其导数的 L² 尾趋零，故 hR→0 于 H¹；归一化后仍 fR→fσ 于 H¹，从而也于 H^(1/2) 及 (18) 的单粒子范数。这里给出的控制强于“Gaussian 方差有限”，保证它由已经证明局域的紧时间窗逼近。

每个 R 的物理时间支撑包含于 t∈β[s₀−2R,s₀+2R]，所以 η 支撑的右端是

η₂,R=−L exp[−β(s₀+2R)/L]。

为了使空间支撑处处留在静态补丁，必须选择

\[
\varepsilon_R R_\chi<-\eta_{2,R}.
\tag{23}
\]

固定非零 ε 然后令 R→∞ 一般会违反 (23)。正确做法是选 Rn→∞ 后，再依每个 Rn 选 εn，使 (23)、所要求的管支撑和

||v_j(Rn,εn)−v_fRn||<1/n

同时成立。因此 v_jn→v_fσ。若要固定物理管半径 r₀<L，可额外要求 εnRχ<r₀ exp[−β(s₀+2Rn)/L]，依然可以完成该对角选择。

Gaussian 覆盖整条中心线。所得 Weyl 算符属于包含整条世界线的开管代数，尤其属于 A₀=A(D)；本证明不赋予它任意固定有限时间段的局域性。

## 8 从单粒子收敛到实际局域代数归属

在所选 BD 自由场表示中，定义

\[
\mathcal A(O)=\{e^{i\phi_{\rm dS}(j)}:
 j\in C_c^\infty(O;\mathbb R)\}^{\prime\prime}.
\tag{24}
\]

对于上面构造的 jn，φdS(jn)=A(vn)，A(v)=a(v)+a†(v)，且 vn→vf。令 W(v)=exp(iA(v))。Fock Weyl 关系与真空期望

⟨Ω,W(w)Ω⟩=exp(−||w||²/2)

给出 W(vn)Ω→W(vf)Ω；在稠密的 Weyl 相干向量 W(w)Ω 上同样成立，因为额外的 Weyl 相位随 vn 连续。W(vn) 都是范数为 1 的酉算符，故收敛扩展到所有向量：

\[
W(\lambda v_n)\xrightarrow{\rm strong}W(\lambda v_f),
\qquad\lambda\in\mathbb R.
\tag{25}
\]

自由 BD 态的共形实现与其局域 Weyl GNS 表示由相同的准自由关联函数识别；若采用该局域代数的循环子空间，(25) 也限制到该子空间。因 A(O) 强闭，得

\[
\boxed{W(\lambda v_f)\in\mathcal A(O)\subseteq\mathcal A_0
\quad\text{对所有实 }\lambda.}
\tag{26}
\]

Segal 场在有限粒子域上本质自伴，(26) 表明其自伴闭包的全部谱投影属于 A(O)，即 A(vf) **附属于** A(O)。这与“无界 A(vf) 本身是有界 von Neumann 代数元素”是不同命题。

紧窗口取任意包含其支撑段的开管 O；Gaussian 按第 7 节取整条世界线的开管。若只需要 A₀ 归属，所有近似始终在 D 内已经足够。该证明给出包含关系，没有证明 A(γ)=A(D) 的反向包含，也没有证明 type 分类。

## 9 传给共同钟接口的统一范数条件

令 N 是此自由 Fock 表示中的粒子数。创建/湮灭算符的基本界给出

\[
\|[A(v_n)-A(v_f)](N+1)^{-1/2}\|
\le 2\|v_n-v_f\|\longrightarrow0,
\qquad \sup_n\|v_n\|<\infty.
\tag{27}
\]

这是实际裸场的图范数控制，强于只比较一个二点矩。

静态等距平移 δ 在该共形动量表示中可实现为

\[
(U_1(\delta)v)(\mathbf k)
=e^{-3\delta/(2L)}v(e^{-\delta/L}\mathbf k).
\tag{28}
\]

它是 L²(d³k) 上的酉群。对 (7) 的中心时间平移 Fδ(t)=F(t−δ)，对应 gδ(η)=g(e^(δ/L)η)，直接 Fourier 变换可核对 (28)。因此

||U₁(δ)(vn−vf)||=||vn−vf||

对所有实 δ 成立。其生成元是静态缩放生成元，具有双侧谱，不能把它替换成 (11) 中正的共形频率 |k|。第二量子化 U(δ)=Γ(U₁(δ)) 保持 N；在与等价 BD 热 GNS 实现作识别时应一并保持对应的实际静态演化，而不凭相同二点函数另选任意 Hamiltonian。

(27) 加上这项同时平移的一致性，是对全线钟关系式酉变换逐纤维取极限所需的输入。它本身不等于证明压缩 Weyl 产品的混合导数交换；后者属于字典 B。尤其 Ππ₀(exp(iλA))Π 一般仍不等于 exp(iλΠπ₀(A)Π)。

## 10 已补齐与仍未主张的内容

本阶段字典 A 在以下明确条件下成立：真实 dS₄ 自由共形标量、BD 准自由态、中心测地线、实窗口，以及由此场的时空 Weyl 算符定义的局域 von Neumann 代数。光滑紧窗、本文 cos² 窗和 Gaussian 已分别给出完整近似。体积因子、交叉范数、空间支撑、时间平滑/截断和极限次序均已指定。

因此，先前审查中“裸世界线窗口是否能从静态管测试函数得到”这一缺口，在此模型中已有直接肯定答案。这两项输入现已接入[字典 B](stage2-weyl-bridge.md)的自伴性、附属关系与有界关系式 Weyl 产品混合导数证明。有限 G 字典、观察者反冲、相互作用场及重整化复合算符都未由这里推出。
