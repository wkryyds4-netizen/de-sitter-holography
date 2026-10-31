# 真实正频观察者种子的轨道共同域

2026-10-29。模型取dS₂、L=1、proper SO⁺(2,1)，两个observer采用KL正频主系列、固定谱下界s_*>0。物质空间K可为任意可分Hilbert空间上的酉群表示；特别包括本项目每个有限粒子数扇区。本文证明真实紧质量包产品种子满足[位置FP稿](stage3-fp-two-particle-inner-product.md) F4的强轨道条件，从而其位置FP内积精确等于群平均。

证明由原模式两端估计、[径向散射定理](stage3-radial-scattering-core.md)和一个非负Tonelli恒等式组成。它不以运动学稠密代替物理稠密；后者另见[产品种子的GA幅度与密度](stage3-product-seed-ga-density.md)。完整局部场作用还须另核。

## 1. 种子域与定理

记单observer包为

\[
u_{n,f}(t,\theta)=\int f(s)\Psi_n^s(t,\theta)ds,
\quad n\in\mathbb Z,\quad f\in C_c^\infty((s_*,\infty)).
\tag{P1}
\]

定义C_orb为有限线性组合

\[
\psi(x_L,x_R)=u_{n_L,f_L}(x_L)u_{n_R,f_R}(x_R)\otimes v,
\qquad v\in K.
\tag{P2}
\]

这些是普通运动学Hilbert向量；有限和在所选observer谱投影后的运动学正频空间中稠密。若研究完整s>0观察者空间，应对s_*>0取并集。这个运动学稠密结论本身不用于推断物理像稠密。

**定理。** C_orb包含于D_orb。对φ、ψ∈C_orb，群平均绝对收敛，轨道Bochner积分定义不变的物质向量值位置波函数ηψ，并满足

\[
Q_{GA}(\phi,\psi)=\langle\eta\phi,\eta\psi\rangle_{FP}.
\tag{P3}
\]

对逆元对称的紧Cartan cutoff，截面值在FP范数中收敛到ηψ的截面值。这里使用的是原六切片的正测度，没有无限群体积相消。下面逐项验证F4的假设。

## 2. 原模式的紧质量包具有快速下降尾

KL [2411.08090v1](https://arxiv.org/html/2411.08090v1) 的(3.11)、(4.9)、A.5—A.6固定本段模式和半相位。固定n，令ℓ=|n|、h=ℓ+1/2。先用(3.11)将径向模式的指标取为−ℓ，避免正整数n时下参数的表观退化。对t的径向部分有

\[
v_s(t):=\sqrt{\cosh t}\,\Psi_n^s(t,0)
=N_\ell(s)\cosh^h t\,
\mathbf F(h-is,h+is;\ell+1;z(t)),
\quad z(t)=\frac{1-i\sinh t}{2}.
\tag{P4}
\]

其中\(\mathbf F={}_2F_1/\Gamma(c)\)，N_ℓ为原A.5前因子乘(4.9)归一化，并选原模式连续半相位。它在任意紧正谱区间光滑、无零点；以下估计只需其各阶导数在该紧集有界，不另改原模式的归一化。由KG方程直接核对

\[
[-\partial_t^2+(1/4-n^2)\operatorname{sech}^2t]v_s=s^2v_s.
\tag{P5}
\]

取s∈I=[s_0,s_1]⊂(s_*,∞)，t=ετ、ε=±1、τ>0、x=e^(−τ)。有精确代数恒等式

\[
-z(\epsilon\tau)=\frac{i\epsilon}{4x}(1+i\epsilon x)^2,
\qquad z(\epsilon\tau)^{-1}=\frac{4i\epsilon x}{(1+i\epsilon x)^2}.
\tag{P6}
\]

由[DLMF 15.8.2](https://dlmf.nist.gov/15.8.E2)的连接公式，两个上参数差2is远离0，得到两端各有两支：

\[
v_s(\epsilon\tau)=\sum_{\sigma=\pm1}
 A_{\epsilon\sigma}(s)e^{i\sigma s\tau}M_{\epsilon\sigma}(x,s),
\quad q_\sigma=h-i\sigma s,
\tag{P7}
\]

\[
A_{\epsilon\sigma}(s)=N_\ell(s)2^{2q_\sigma-h}
 e^{-i\epsilon\pi q_\sigma/2}
\frac{\Gamma(2i\sigma s)}{\Gamma(h+i\sigma s)\Gamma(1/2+i\sigma s)},
\]

\[
M_{\epsilon\sigma}(x,s)=(1+x^2)^h(1+i\epsilon x)^{-2q_\sigma}
{}_2F_1\left(q_\sigma,\tfrac12-i\sigma s;1-2i\sigma s;
\frac{4i\epsilon x}{(1+i\epsilon x)^2}\right).
\tag{P8}
\]

两端−z的相位趋近±π/2，使用与P4连续的复幂分支。没有删除任一频率支。

在I的足够小复邻域及|x|<1/16上，超几何变量模小于1，参数不触及极点；[DLMF 15.2](https://dlmf.nist.gov/15.2)的Gauss级数和参数解析性给联合解析。P8中两个一次项分别是−2iεq_σx与+2iεq_σx，因此严格相消。故

\[
M_{\epsilon\sigma}(x,s)=1+x^2Q_{\epsilon\sigma}(x,s)
\]

且对任意固定j,r及τ≥log32，

\[
\sup_{s\in I}|\partial_s^j\partial_\tau^r(M-1)|
\le C_{n,I,j,r}e^{-2\tau}.
\tag{P9}
\]

若把相位也放进余项，s导数会产生τ，因此完整余项的界应为C(1+τ)^j e^(−2τ)，不能遗漏该因子。

对f∈Cc∞(I°)，P7的主项质量积分是紧支撑光滑函数的Fourier变换；反复分部积分没有端点项，余项及其固定阶t导数指数下降。有限t区间的模式及其导数连续有界。于是

\[
\boxed{\sqrt{\cosh t}\,u_{n,f}(t,0)\in\mathcal S(\mathbb R),}
\]

\[
\sup_t(1+|t|)^N\left|\partial_t^r
 [\sqrt{\cosh t}\,u_{n,f}(t,0)]\right|
\le C_{n,I,N,r}\sum_{j=0}^N\|f^{(j)}\|_{L^1(I)}.
\tag{P10}
\]

常数依赖固定磁指标和紧谱集，不主张对无限磁指标或支撑触及s=0的一致性。取模不必保留光滑性；P10已保证

\[
Q_u(t):=\sqrt{2\pi\cosh t}\,|u_{n,f}(t,0)|\in L^2(dt),
\qquad M_u:=\int(1+|t|)Q_u(t)dt<\infty.
\tag{P11}
\]

√(2π)来自θ∈[0,2π)的完整位置测度，与径向散射稿的酉约定相同。

## 3. 所需自相关在群上平方可积

设f_u(x)=|u_{n,f}(x)|。因单磁指标的角向相位模为1，f_u对K=SO(2)不变。令

\[
A_{f_u}(k)=\langle f_u,U_X(k)f_u\rangle_{L^2(X)}
=\int_X f_u(x)f_u(k^{-1}x)dx\ge0.
\tag{P12}
\]

[径向散射稿](stage3-radial-scattering-core.md)R1—R16证明：H₀=−∂t²+sech²t/4没有束缚态或零能共振；两通道酉ds变换给q(s)，模式阈值抑制抵消群Schur权的s⁻²。其结论应用于P11给

\[
\|A_{f_u}\|_{L^2(G)}\le C M_u\|Q_u\|_2<\infty.
\tag{P13}
\]

该论证允许f_u的径向谱抵达0。observer原来的质量下限不能转嫁给取模之后的谱。

## 4. 六轨道整体的非负恒等式

记F(y)=f_L(x_L)f_R(x_R)、y=(x_L,x_R)，其中f_i=|u_i|。由[FP稿](stage3-fp-two-particle-inner-product.md)F2，非退化位置对可写为y=g s(z)，dy=dg dμ(z)，Z为六条切片。物质表示V酉，故untwist之后

\[
a_\psi(z)=\|v\|\int_G F(g s(z))dg.
\]

全部 integrand 非负，可先用Tonelli而不预设有限性。令第二个群变量h=kg，使用群的单模性和轨道测度，逐步得到

\[
\begin{aligned}
\|a_\psi\|_2^2
&=\|v\|^2\int_Zd\mu(z)\int_Gdg\int_Gdh\,
 F(g s(z))F(h s(z))\\
&=\|v\|^2\int_Gdk\int_Ydy\,F(y)F(ky)\\
&=\|v\|^2\int_Gdk\,A_{f_L}(k^{-1})A_{f_R}(k^{-1})\\
&=\|v\|^2\int_Gdk\,A_{f_L}(k)A_{f_R}(k).
\end{aligned}
\tag{P14}
\]

最后一步仅用Haar逆元不变。因f_i非负实函数，各A亦非负；也可用A(k⁻¹)=A(k)直接核对。退化轨道在原位置测度中为零，本定理不增加集中于其上的分布态。

P13与Cauchy–Schwarz给

\[
\boxed{\|a_\psi\|_2^2
\le C\|v\|^2M_{u_L}M_{u_R}\|Q_{u_L}\|_2\|Q_{u_R}\|_2<\infty.}
\tag{P15}
\]

这控制的是六轨道整体，不只某一条boost渐近。有限和满足a_(Σψ_j)≤Σa_ψ_j，故整个C_orb包含于D_orb。P3及cutoff收敛随即由FP稿F6—F7推出。

## 5. 得到的物理空间与尚需区别的完成化

C_orb上的FP半范数与GA半范数相同，零范数商相同。其像在六切片L²(Z;K)中的闭包因此与C_orb的GA完成化等距。该结论对任意上述酉物质K成立；它没有使用某个N特有的物质衰减假设。

实际产品种子的幅度公式与物理像稠密由[产品种子谱稿](stage3-product-seed-ga-density.md)在N=2、左K=0小域中证明；[有限K延拓](stage3-finite-k-ga-extension.md)另核全部C_orb态的范数与交叉配对，证明其像也在同一小域闭包中。因此N=2的整个C_orb生成闭包已与完整normalized-CG multiplicity空间等距。每个有限N的对应推广见[全N生成空间稿](stage3-all-n-generated-space.md)。还应区分由这些原模式种子定义的完成化，与“所有满足某种分布谱约束且FP范数有限的函数”构成的最大空间：后者如果作为独立定义采用，须另证没有额外扇区，不能用抽象Hilbert同构略过。

本稿没有证明KL原(5.4)局部模式核在物理空间上的全通道弱作用，也没有宣布全N物理场本质自伴。下一步从已经可控制的共同域推进这些识别。已有位置FP充分定理与散射理论是明确输入，本文为其在真实KL制备态上的应用，不作未经查新的原创性声明。


## 6. 环境紧支撑测试与最大空间比较的起点

以上径向结果还给一个不依赖正频限制的引理。对任意φ∈Cc∞(Y;K)，固定其紧支撑C，选非负径向h_L,h_R∈Cc∞(ℝ)，在C的两个时间投影上各≥1。角变量本来紧，故

\[
\|\phi(y)\|\le M h_L(t_L)h_R(t_R),\qquad
M=\sup_C\|\phi\|.
\]

h_i的Q_i光滑紧支撑，P13—P14原样适用，得到

\[
\|a_\phi\|_{L^2(Z)}\le C_C\sup_C\|\phi\|.
\tag{P16}
\]

因此环境空间的紧支撑光滑测试函数都在D_orb，且强域半范数受固定紧集上的一致范数控制。

对任意k∈L²(Z;K)，在非退化位置对上定义等变延拓Rk(gs(z))=V(g)k(z)。单模性及非负Tonelli给

\[
\int_Y\|Rk(y)\|\,\|\phi(y)\|dy
=\int_Z\|k(z)\|a_\phi(z)d\mu(z)
\le\|k\|_2\|a_\phi\|_2.
\tag{P17}
\]

取φ=χv₀，其中χ≥1_C紧支撑、v₀为固定单位K向量，得∫_C‖Rk‖≤C_C‖k‖₂。于是

\[
R:L^2(Z;K)\longrightarrow L^1_{loc}(Y;K)
\tag{P18}
\]

连续，且给出零阶位置分布。若Rk为零，用dy=dg dμ和任意有限正Haar测度的群集合可得k=0几乎处处；注入性不依靠对L¹函数作切片取值。

这为比较更大的位置态类提供了明确起点：FP有限截面确实有合法位置分布。尚未证明正频谱投影在这一拓扑中连续，也未证明满足正频约束的所有此类分布都可由C_orb逼近。两项不能由P18直接推出。

## 7. 有限点公式诊断

[诊断代码](../code/validate_orbit_core_identities.py)及[70位原始结果](../results/stage3/orbit-core/orbit-core-identity-70dps.json)保存51项检查：36个双端连接式点、9个Jost微分方程点、3个Wronskian符号和3个透射模平方。全部通过预设1e−55容差，最大缩放差约6.04×10⁻⁷¹。原式与连接式使用同一mpmath库，可能共享内部算法；它们是有限点实现诊断，不是独立库复验、区间认证或共同域/稠密性证明。解析论证在上文及所引径向稿中。
