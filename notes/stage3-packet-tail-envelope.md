# 有限外包正二跳积分的解析高谱尾包络

2026-10-29。本文核查并补充本项目 dS₂、G→0 模型的正谱尾估计。对于 σ=1、外谱标签 d∈[.7,.9]，内部谱 cap R=20 之外的二步态范数平方遗漏满足严格数字界 **ε₂₀<10⁻³⁸**。这个结论不依赖有限 cap 内的数值求积，也不把该求积变成区间认证。

对象是 [既定正核 Vν(a∣d)](../code/quantized_observer.py) 与其二跳积分，定义见[四点域记录](stage3-fourpoint-domains.md)、[历史截止正交性](stage3-crossing-integral-audit.md)。本文不重新证明其与原文全部多粒子 CG 域、闭场的识别，也不覆盖有限 G。以下是对该已定义正核的解析补充，未作首创性声明。

## 1. 对象与可用的统一 Gamma 界

固定实 principal-series ν>0，σ>0，下谱投影 s*>0，精确质量

\[
 m(s)=\sqrt{s^2+\tfrac14},\qquad M=m(d).
\]

定义

\[
 N_R(d)=\int_{s_*}^{R}da\int_{s_*}^{R}db\,
 V_\nu(a|d)V_\nu(b|a)
 e^{-\sigma^2[(m(a)-M)^2+(m(b)-m(a))^2]},
 \quad\varepsilon_R(d)=N_\infty(d)-N_R(d).
\tag{T1}
\]

这里 a,b,d,R 都是谱标签，R 对应的物理质量为 K=m(R)，不是 K=R。两个 Gaussian 来自二步态振幅取模平方，故指数系数为 −σ²；它们不是交叉四点振幅中 −σ²/2 的四边窗。

正核按既定归一化为

\[
 V_\nu(a|d)=\frac{a\tanh(\pi a)\cosh^2(\pi d)}{16\pi^5}
 \prod_{\epsilon_1,\epsilon_2,\epsilon_3=\pm1}
 \Gamma\!\left(\tfrac14+\tfrac i2(\epsilon_1\nu+
 \epsilon_2a+\epsilon_3d)\right).
\tag{T2}
\]

八 Gamma 成四个共轭对，乘积实非负。由 Euler 积分，对于实 y，

\[
 |\Gamma(\tfrac14+iy)|
 \le\int_0^\infty t^{-3/4}e^{-t}dt=\Gamma(\tfrac14).
\]

再用 tanh(πa)≤1，得到对所有实 ν 一致的粗界

\[
 0\le V_\nu(a|d)\le C_0\,a\cosh^2(\pi d),
 \qquad C_0=\frac{\Gamma(1/4)^8}{16\pi^5}.
\tag{T3}
\]

它很宽松，但足够控制大 cap 的余尾。无需使用振荡抵消或未经证明的 Gamma 渐近。

## 2. 质量变换与补平方

令 x=m(a)、y=m(b)。因为 a da=x dx、b db=y dy，两核乘积及测度满足

\[
 V(a|d)V(b|a)\,da\,db
 \le C_0^2\cosh^2(\pi d)\,xy\cosh^2(\pi a)\,dx\,dy.
\]

在 x,y>0 的原域上，cosh²(πa)≤e²πᵃ≤e²πˣ 且 xy≤(x²+y²)/2。设 s=σ²、

\[
 \mu=M+\frac\pi s,\qquad Q=2\pi M+\frac{\pi^2}s.
\]

直接展开可核对

\[
 -s[(x-M)^2+(y-x)^2]+2\pi x
 =Q-s\,(x-\mu,y-\mu)
 \begin{pmatrix}2&-1\\-1&1\end{pmatrix}
 \binom{x-\mu}{y-\mu}.
\tag{T4}
\]

该 Gaussian 的全平面积分是 π/s。归一化后视为随机向量 Z=(X,Y)，有

\[
 \mathbb E Z=(\mu,\mu),\qquad
 \Sigma=\frac1{2s}\begin{pmatrix}1&1\\1&2\end{pmatrix}.
\tag{T5}
\]

Σ 是协方差，指数矩阵 s[[2,−1],[−1,1]] 等于 Σ⁻¹/2；不得遗漏这个 2。

原谱尾域是 x,y≥m(s*) 且 max(x,y)>K。先用非负二次式 H=(x²+y²)/2 替换 xy，再作 union bound、延至整个平面。因此

\[
 \varepsilon_R(d)\le A_d
 \sum_{j=1}^2\mathbb E[H(Z)\mathbf1_{Z_j\ge K}],
 \qquad A_d=\frac{\pi C_0^2}{s}\cosh^2(\pi d)e^Q.
\tag{T6}
\]

union 会重复算 x,y 都超过 K 的部分，方向仍是上界。必须先替换成 H：若直接把 xy 延至全平面，负象限贡献会破坏该比较。

## 3. Chernoff 二次矩：两项系数的来源

要求 K>μ，记 Δ=K−μ>0。对 j=1,2 取 λ_j=Δ/Σ_jj。则

\[
 \mathbf1_{Z_j\ge K}\le e^{\lambda_j(Z_j-K)},
\]

而 Gaussian 倾斜公式给

\[
 \mathbb E[H(Z)e^{\lambda_j(Z_j-K)}]
 =e^{-\Delta^2/(2\Sigma_{jj})}
 \mathbb E_j H(Z),
\tag{T7}
\]

其中倾斜后协方差仍为 Σ，均值是 (μ,μ)+λ_jΣe_j。

| 被界定的尾 | λ_j | 倾斜均值 | E_j H |
|---|---|---|---|
| X≥K | 2sΔ | (K,K) | K²+3/(4s) |
| Y≥K | sΔ | ((K+μ)/2,K) | K²/2+(K+μ)²/8+3/(4s) |

这就得到

\[
\boxed{\begin{aligned}
 B_R(d)=\frac{\pi C_0^2}{\sigma^2}\cosh^2(\pi d)e^Q
 \bigg[&\left(K^2+\frac{3}{4\sigma^2}\right)
 e^{-\sigma^2(K-\mu)^2}\\
 &+\left(\frac{K^2}{2}+\frac{(K+\mu)^2}{8}
 +\frac{3}{4\sigma^2}\right)
 e^{-\frac12\sigma^2(K-\mu)^2}\bigg],
 \qquad\varepsilon_R(d)\le B_R(d).
\end{aligned}}
\tag{T8}
\]

Chernoff 的 λ 优化这里只最小化指数部分，未声称对带 H 的整体矩达到最优。合法上界不需要这一额外最优性。

## 4. 正规化外包与 d 的单调性

固定 σ、R，并令外标签 0≤d≤d_max。若 K>μ(d_max)，则 T8 的各正因子均随 d 增加：M、μ、Q 和 cosh²πd 增加；K−μ保持正且减小，故两指数增加；第二多项式也增加。因此

\[
 \sup_{0\le d\le d_{\max}}\varepsilon_R(d)
 \le B_R(d_{\max}).
\tag{T9}
\]

任意正规化 P(d)dd 支撑于此区间时，

\[
 \varepsilon_{R,P}=\int P(d)\varepsilon_R(d)dd
 \le B_R(d_{\max}).
\tag{T10}
\]

这包括当前 [.7,.9] 上的 sin⁴ 包，且不需要外包求积来证明 T10。这里乘的是已经正规化的 P(d)dd，不再乘 Plancherel 密度。

## 5. σ=1 的高精度评价与严格粗数字界

在 d_max=.9，50 位 mpmath 评价给

\[
 C_0\simeq6.09788534362667575987,\quad
 M\simeq1.02956301409870003158,\quad
 \mu\simeq4.17115566768849327004.
\]

| 谱 cap R | 物理质量 K | T8 的 50 位数值评价 |
|---:|---:|---:|
| 16 | 16.0078105935821217162 | 7.09752404005008966462×10⁻¹⁸ |
| 20 | 20.0062490237425566241 | 1.01870686617260239720×10⁻⁴¹ |

可将数值评价保守展示为约 7.10×10⁻¹⁸、1.02×10⁻⁴¹，但普通高精度浮点不提供有向舍入证书。严格解析结论首先是 T8–T10。下面另给不用浮点舍入认证的较粗严格数字界。

由 logΓ 在 [1,2] 上严格凸且 Γ(1)=Γ(2)=1，Γ(5/4)<1，故 Γ(1/4)=4Γ(5/4)<4。再用 3<π<22/7，

\[
 C_0<\frac{4^8}{16\,3^5}=\frac{4096}{243}<17.
\]

对 d≤9/10，有 M²≤53/50，因而 M<26/25。于是

\[
 \cosh^2(\pi d)<e^6,\quad
 Q<2\frac{22}{7}\frac{26}{25}+\left(\frac{22}{7}\right)^2
 =\frac{20108}{1225}<17,
\]

\[
 \mu<\frac{26}{25}+\frac{22}{7}=\frac{732}{175}<5.
\]

R=20 时 K>20，故

\[
 \frac{(K-\mu)^2}{2}>
 \frac12\left(\frac{2768}{175}\right)^2
 =\frac{3830912}{30625}>125.
\]

又 K²=400.25<401、K<21，T8 的两个多项式分别小于 402、286。将第一个更快的 e⁻Δ² 也放宽至 e⁻¹²⁵，得到

\[
 B_{20}(d)<4\cdot17^2\cdot(402+286)e^{6+17-125}
 =795328\,e^{-102}.
\]

Euler 常数 e 的前五项 Taylor 和给 e>65/24>27/10，且精确整数比较

\[
 795328\,10^{140}<27^{102}
\]

成立。因此

\[
 \boxed{\sup_{d\in[.7,.9]}\varepsilon_{20}(d)
 \le B_{20}(.9)<10^{-38},\qquad\varepsilon_{20,P}<10^{-38}.}
\tag{T11}
\]

最后的整数不等式已用整数运算核验，完全独立于上述 mp50 小数评价。严格粗界比锐闭式评价宽松约三阶，仍足以排除 cap=20 外尾在当前精度上的影响。

## 6. 与 Λ=3 的实际截断误差怎样合并

R 是正核尾积分的远端 cap；Λ=3 是当前交叉积分的内部截止。二者不可混同。对准确数学积分，

\[
 \varepsilon_\Lambda=(N_R-N_\Lambda)+(N_\infty-N_R),
 \quad0\le N_\infty-N_R\le B_R.
\tag{T12}
\]

所以 T11 认证的是有限 cap 之后的遗漏；Λ 到 R 之间的正尾仍需求积。已有历史通道配对界可写为

\[
 \sup_t|X_{P,\infty}(t)-X_{P,\Lambda}(t)|
 \le\sqrt{\varepsilon_{\Lambda,P}}
 \left(\sqrt{N_{\infty,P}}+\sqrt{N_{\Lambda,P}}\right).
\tag{T13}
\]

将 T12 的解析界代入 T13，在使用准确 N_R、N_Λ 时是严格的；若只使用普通自适应求积的估计值，有限区间求积误差仍未被严格上界覆盖。当前能准确描述为：**cap 外尾已由解析式控制，有限域交叉求积和 cap 内正核求积仍是经验数值检查。** 不能据 T11 宣称整个 Λ=3 截断误差或最终四点总误差已获区间认证。

## 7. 可选的收紧方式

T6 的 Gaussian 半空间二次矩也可直接计算，省去 Chernoff 的放宽。令 z₁=√(2s)(K−μ)、z₂=√s(K−μ)，φ(z)=e⁻ᶻ²ᐟ²/√(2π)，Q_G(z)=erfc(z/√2)/2，C=μ²+3/(4s)。则

\[
 \mathbb E[H\mathbf1_{X\ge K}]
 =C Q_G(z_1)+\frac{K+\mu}{\sqrt{2s}}\phi(z_1),
\]

\[
 \mathbb E[H\mathbf1_{Y\ge K}]
 =C Q_G(z_2)+\frac{5K+7\mu}{8\sqrt s}\phi(z_2).
\tag{T14}
\]

这由条件 Gaussian 二次矩及 ∫_z^∞u²φ(u)du=Q_G(z)+zφ(z)直接得到，仍保留 union 的重复部分。固定 d=.9、σ=1 时，其右侧乘 A_d 的 mp50 评价在 R=16、20 分别约为 2.39046×10⁻¹⁹、2.56581×10⁻⁴³。当前统一外包严格结论只使用已逐项证明单调性的 T8 与粗数字界 T11，不需要依赖这个更紧表达。
## 8. 不借助数值范数的严格配对尾数字界

T6 的同一全平面包络若不加尾事件，就给

\[
 N_\infty(d)\le
 \frac{\pi C_0^2}{\sigma^2}\cosh^2(\pi d)e^Q
 \left(\mu^2+\frac{3}{4\sigma^2}\right).
\tag{T15}
\]

σ=1、d≤.9 时，前面已经严格得到 prefactor<1156e²³、μ<5。又 e<3（由 n!≥2ⁿ⁻¹ 的级数比较，且 n≥3 时严格），于是

\[
 N_\infty(d)<1156\,3^{23}\left(25+\frac34\right)
 =2802360004143309<3\times10^{15}.
\tag{T16}
\]

虽然这个全范数上界极为宽松，它足以与 T11 合用，完全避开未认证的数值 N₂₀。以 R=20 作为交叉三内标签共同上界，使用[交叉积分审查](stage3-crossing-integral-audit.md)中已证的历史截止配对不等式，有

\[
 \begin{aligned}
 \sup_{t\in\mathbb R}|X_{d,\infty}(t)-X_{d,20}(t)|
 &\le\sqrt{\varepsilon_{20}(d)}
 \left(\sqrt{N_\infty(d)}+\sqrt{N_{20}(d)}\right)\\
 &\le2\sqrt{\varepsilon_{20}(d)N_\infty(d)}
 <2\times10^{-11}.
 \end{aligned}
\tag{T17}
\]

最后一步仅用 ε₂₀N∞<3×10⁻²³<10⁻²²。它对当前时间配置 (t,0,t,0) 的全部实 t 成立；时间相位不改变相应二步态范数或其投影尾。

对于当前正规化有限包，∫P(d)dd=1，ε₂₀,P<10⁻³⁸ 且 N∞,P<3×10¹⁵，所以同一证明给

\[
 \boxed{\sup_{t\in\mathbb R}
 |X_{P,\infty}(t)-X_{P,20}(t)|<2\times10^{-11}.}
\tag{T18}
\]

这里 X_P,20 是准确数学意义上整个内部 [s*,20]³ 的积分，不是声称已经在该大盒完成数值求积。实际生产交叉网格仍为 Λ=3。T18 严格控制 cap=20 外的配对遗漏；它不认证 [3,20] 的正核尾积分，不认证 [s*,3]³ 的交叉求积，也不提供整个生产结果的严格总误差条。