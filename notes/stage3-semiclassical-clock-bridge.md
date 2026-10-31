# dS₂ 大质量极限与共同钟非交叉尾核的桥接

日期：2026-10-29。原始输入是 Kolchmeyer–Liu [2411.08090v1，式 (6.18)—(6.21)](https://arxiv.org/html/2411.08090v1) 的大质量渐近，以及 [四点通道证明](stage3-fourpoint-domains.md) 的精确非交叉表达。本稿核定系数、变量、投影阈值和极限次序，并证明这两个通道对全部实时间的绝对一致收敛。它只建立同一 dS₂ 模型中两类算符矩的接口；不建立 dS₂ 与 dS₄ 的跨维度等价，也不处理交叉 ABAB 的 6j 反冲极限。

## 1 保持不变的物理任务

固定：物质主系列参数 ν>0；proper-time Gaussian 宽度 σ>0；归一化平坦钟包 ζ(p)，p≥0；每个插入间的观察者质量下界 m≥m₀。令 m₀→∞，保持 p=2π(m−m₀) 而不是保持整个质量 m 不变。

相应物理外态随 m₀ 变化，按此前等距映射为

\[
f_{\zeta,m_0}(s)=\mathbf1_{m(s)\ge m_0}
\sqrt{\frac{2\pi s}{m(s)\rho_{Pl}(s)}}\,
\zeta(2\pi[m(s)-m_0]),
\quad m(s)=\sqrt{s^2+\tfrac14}.
\tag{1}
\]

两侧真空质量相等，外部平坦概率 P(p)=|ζ(p)|²始终保持不变。所有右插入都保留 Qm₀=1[HR≥m₀]1[HL≥m₀]，不是只对初态取一次下界。时间采用 L=1 的 proper time τ；模时间为 τ/(2π)。

本轮紧支撑双峰满足以下证明全部条件。证明实际上适用于任意**固定**正规化 P∈L¹(R₊)；如果 P 也随 m₀ 展宽、平移或改变尾部，需要另加控制，不能直接套用。

## 2 单个顶角的极限与精确 Jacobian

令 a=s(m′)、b=s(m)，其中 s(m)=√(m²−1/4)，质量差 δ=m′−m 固定。八 Gamma 核可以分解成两组四 Gamma，分别依赖 a−b 和 a+b。定义

\[
\Gamma_4(\nu,x)=\prod_{\epsilon,\epsilon'=\pm1}
\Gamma\!\left(\frac14+\frac i2[\epsilon\nu+\epsilon' x]\right)
=\left|\Gamma\!\left(\frac14+\frac i2(\nu+x)\right)\right|^2
 \left|\Gamma\!\left(\frac14+\frac i2(\nu-x)\right)\right|^2.
\tag{2}
\]

于是 Γ₈(ν,a,b)=Γ₄(ν,a−b)Γ₄(ν,a+b)。Stirling 公式在 S=a+b→∞、ν固定时给

\[
\Gamma_4(\nu,S)\sim\frac{8\pi^2}{S}e^{-\pi S}.
\tag{3}
\]

把它代入精确 Vν(a|b)，使用 cosh²(πb)∼e²πb/4、a/(a+b)→1/2、a−b→δ，得到

\[
\boxed{V_\nu(s(m+\delta)|s(m))\frac{ds(m+\delta)}{dm'}
\longrightarrow
v_\nu(\delta):=\frac{e^{-\pi\delta}}{16\pi^3}\Gamma_4(\nu,\delta).}
\tag{4}
\]

这里 ds/dm′=m′/s(m′)→1。系数 1/(16π³) 与原文 (6.21) 一致；没有将未归一化 HH 外态的 ZHH 因子错误地再除一次。v 是质量差的谱权，不是归一化跃迁概率。

真正适合证明全部实时间收敛的变量是每段平坦钟能量，而非用 a−b 近似精确质量差。令

\[
m_p=m_0+\frac p{2\pi},\quad s_p=s(m_p),\quad
W_{m_0}(q|p)
=\frac{m_q}{2\pi s_q}V_\nu(s_q|s_p)
\exp\!\left[-\left(\frac{\sigma}{2\pi}\right)^2(q-p)^2\right].
\tag{5}
\]

q 是插入之后的钟能量，积分范围 q≥0。因 da=(m_q/(2πs_q))dq，(5) 包含完整 Jacobian。其逐点极限为平移核

\[
\boxed{W_\infty(q|p)=w_\sigma(q-p),\qquad
w_\sigma(z)=\frac1{2\pi}v_\nu\!\left(\frac z{2\pi}\right)
\exp\!\left[-\left(\frac{\sigma}{2\pi}\right)^2z^2\right].}
\tag{6}
\]

同时每段相位**精确**等于 exp[i(q−p)u/(2π)]，对任何 m₀ 都成立。保留 m(s) 而不以 s 替换它，正是后面时间一致性的要点。

## 3 与 m₀ 无关的可积包络

[独立 Gamma 核核验](../literature/stage3-observer-kernel-audit.md) 给出固定 ν 的全质量上界

\[
V_\nu(a|b)\le C_\nu\frac{e^{-2\pi(a-b)_+}}{1+|a-b|}.
\tag{7}
\]

为了避免 s=0 的 Jacobian 端点，只需在极限尾段取 m₀≥1。因为 ds/dm=m/s≥1，故

\[
|s_q-s_p|\ge\frac{|q-p|}{2\pi},\qquad
(s_q-s_p)_+\ge\frac{(q-p)_+}{2\pi},\qquad
\frac{m_q}{s_q}\le\frac2{\sqrt3}.
\tag{8}
\]

由 (5)—(8)，存在与 m₀、p、q 无关的有限常数，使

\[
0\le W_{m_0}(q|p)\le H(q-p),
\]

\[
\boxed{H(z)=\frac{C_\nu}{\pi\sqrt3}
\frac{e^{-z_+}}{1+|z|/(2\pi)}
\exp\!\left[-\left(\frac{\sigma}{2\pi}\right)^2z^2\right]
\in L^1(\mathbb R).}
\tag{9}
\]

这个常数只用于存在性证明，当前没有提供可用于区间认证的最优显式值。W∞同样受 H 支配。因此所有 p≥0 的核质量均有一致上界 ||H||₁。

## 4 两条通道的 L¹ 收敛及全部实时间一致性

时间标记沿用从右到左 t₁→t₄，即显示乘积为 φ(t₄)φ(t₃)φ(t₂)φ(t₁)。在平坦钟变量下，两个非交叉谱密度分别是

\[
F^{adj}_{m_0}(p,q,r)
=P(p)W_{m_0}(q|p)W_{m_0}(r|p),
\tag{10}
\]

\[
F^{nest}_{m_0}(p,q,r)
=P(p)W_{m_0}(q|p)W_{m_0}(r|q),
\tag{11}
\]

三个变量均≥0。它们的可积控制为 P(p)H(q−p)H(r−p) 和 P(p)H(q−p)H(r−q)，总积分各≤||H||₁²。因此支配收敛给出

\[
\boxed{\|F^{adj}_{m_0}-F^{adj}_\infty\|_{L^1}\to0,
\qquad
\|F^{nest}_{m_0}-F^{nest}_\infty\|_{L^1}\to0.}
\tag{12}
\]

相邻通道相位为

\[
\exp\!\left\{\frac i{2\pi}[(q-p)(t_2-t_1)+(r-p)(t_4-t_3)]\right\},
\]

嵌套通道相位为

\[
\exp\!\left\{\frac i{2\pi}[(q-p)(t_4-t_1)+(r-q)(t_3-t_2)]\right\}.
\]

这些相位模长为 1，且与 m₀ 无关。故

\[
\boxed{\sup_{\mathbf t\in\mathbb R^4}
|C^{adj}_{m_0}(\mathbf t)-C^{adj}_\infty(\mathbf t)|
\le\|F^{adj}_{m_0}-F^{adj}_\infty\|_1\to0,}
\tag{13}
\]

嵌套通道有完全相同的式子。二点也由 P(p)W(q|p) 的 L¹ 收敛得到全部实时间绝对一致极限。

这是固定质量包形状、固定物质与 Gaussian 任务的**绝对**一致性，不是任意时刻分母可能很小的相对误差控制。证明没有给出一般 m₀ 收敛速率。若在有限质量阶段先把质量相位 m(a)−m(b) 换成 a−b，则由 t 放大的相位误差不再自动满足这项证明。

可供实现的误差分解是定义

\[
\delta_{m_0}(p)=\int_0^\infty dq\,|W_{m_0}(q|p)-W_\infty(q|p)|.
\tag{14}
\]

若 M=||H||₁，则相邻通道 L¹ 差≤2M∫P(p)δm₀(p)dp；嵌套通道可界为

\[
M\int P(p)\delta_{m_0}(p)dp
+\int P(p)dp\int_0^\infty dq\,W_\infty(q|p)\delta_{m_0}(q).
\tag{15}
\]

这展示了为什么嵌套通道必须控制中间质量处的下一顶角，不能只验证初态支持区间内一个二点数值。

## 5 回到场模能量：共同钟阈值从哪里来

设场模能量 E=−2πδ。因此一次跃迁使钟从 p 变成 p−E。定义实际 dS₂ 场谱

\[
\rho_\nu(E)=\frac1{2\pi}v_\nu\!\left(-\frac E{2\pi}\right)
=\frac{e^{E/2}}{32\pi^4}\Gamma_4\!\left(\nu,\frac E{2\pi}\right),
\quad \rho_\nu(-E)=e^{-E}\rho_\nu(E),
\tag{16}
\]

及涂抹后有限正测度

\[
d\mu_{\nu,\sigma}(E)=\rho_\nu(E)
 e^{-(\sigma/(2\pi))^2E^2}\,dE,
\qquad Z_{\nu,\sigma}=\int d\mu_{\nu,\sigma}<\infty.
\tag{17}
\]

这里 μ 尚未除以 Z；与先前使用归一化单场 a=A/√Z 的记号比较时，四点应相应除 Z²。KMS 符号与 proper/mod 时间变换都已固定：exp[iδu]=exp[−iE u/(2π)]。

令 K_P(x)=∫x∞P(p)dp，x≥0。相邻图的质量路径是 p→p−E→p→p−F→p。每段钟能量非负要求 p≥max(0,E,F)，故

\[
\boxed{C^{adj}_\infty(\mathbf t)
=\iint d\mu(E)d\mu(F)\,
K_P(\max(0,E,F))
 e^{-\frac i{2\pi}[E(t_2-t_1)+F(t_4-t_3)]}.}
\tag{18}
\]

嵌套图的路径是 p→p−E→p−E−F→p−E→p，阈值为 max(0,E,E+F)=[E+F₊]₊，从而

\[
\boxed{C^{nest}_\infty(\mathbf t)
=\iint d\mu(E)d\mu(F)\,
K_P(\max(0,E,E+F))
 e^{-\frac i{2\pi}[E(t_4-t_1)+F(t_3-t_2)]}.}
\tag{19}
\]

二点同理给 ∫dμ(E)K_P(E₊)exp[−iEτ/(2π)]。两个四点阈值恰好是先前共同钟非交叉通道的 q₁₂ 和 q₁₄；此处它们从量子化有限质量顶角、实际 Jacobian 和逐次投影的大质量极限推导而来，没有先把旧尾核乘到有限质量公式上。

式 (18)、(19) 也解释了外质量包与中间投影怎样成为同一积分中的资源条件。先只截外态而不在各插入间放 Qm₀，会保留负 q 或 r 的中间区间，不能得到这些 K_P 阈值。

## 6 成果边界及数值比较要求

这条桥建立了同一 dS₂、G→0 模型中两种完整混合 species 四点（AABB、ABBA）的大质量算符矩接口。没有证明一个世界线生成整个静态补丁代数，没有识别完整观测者代数，也没有把 dS₂ 主系列谱当成阶段二 dS₄ 共形谱。共同钟阈值结构可复用，具体场谱与物理量不能跨维度直接替换。

ABAB 的交叉核含 6j 再耦合，其大质量极限还出现奇异 delta 结构；这里的 L¹ 控制没有覆盖它。因此不从 (13) 推断完整 OTOC 对全部实时间一致，更不否定 τ∼log m₀ 的反冲增强。有限 G 也未进入此推导。

数值比较应固定同一个 P、ν、proper σ、时间顺序及逐次质量投影，并保留 m(s) 的精确相位。分别比较有限质量的平坦 q 核与 (6)、三变量通道与 (18)—(19) 的旧式尾核表示。至少在零时和非零时检查；有限阶求积与理论 m₀ 极限误差分开记录。不能仅凭 C₂ 接近就断言两条四点都同样接近，也不能用经验拟合的幂律替代尚未给出的收敛速率证明。