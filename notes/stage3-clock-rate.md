# 非交叉大质量极限：全域速率与首阶修正

日期：2026-10-29。模型固定为 Kolchmeyer–Liu [2411.08090v1 §6.2—6.3](https://arxiv.org/html/2411.08090v1) 的 dS₂、L=1、G→0、自由主系列物质；两 species 同质量。本文是从已有精确顶角出发的误差定理与首项推导，不作首创性声明。此前[无速率证明](stage3-semiclassical-clock-bridge.md)保持原样；[四点定义](stage3-fourpoint-domains.md)仍控制物理任务及定义域边界。数学输入是 [DLMF 5.11.1—2](https://dlmf.nist.gov/5.11) 的 Stirling/digamma 渐近及 [5.7.6](https://dlmf.nist.gov/5.7#E6) 的 digamma 级数。

**结论。** 固定 ν>0、proper-time Gaussian 宽度 σ>0，任意固定 P≥0、∫₀∞P(p)dp=1，保留每次插入间的质量下界投影。两点及完整混合 species AABB、ABBA 满足全实时间绝对误差 O(m₀⁻¹)，并有全实时间一致首阶展开 C_m₀=C∞+C₁/m₀+o(m₀⁻¹)。不要求 ∫pP 有限。首阶四点核是在相应极限尾核上乘 −(E+F)/(4π)。若另有 ∫pP<∞，一阶展开余项可加强到 O(m₀⁻²)。这些结论均不包含 ABAB 交叉 6j、相对误差、有限 G 或代数类型的收敛。

## 1 精确核及符号

以下 M=m₀≥1，p,q≥0，m=M+p/(2π)，m′=M+q/(2π)，a=√(m′²−1/4)，b=√(m²−1/4)，δ=m′−m=(q−p)/(2π)，d=a−b，S=a+b。δ 是精确总质量差，d 是主系列谱标签差，不混用。

\[
\Gamma_4(\nu,x)=\left|\Gamma\left(\tfrac14+\tfrac i2(\nu+x)\right)\right|^2
\left|\Gamma\left(\tfrac14+\tfrac i2(\nu-x)\right)\right|^2,
\qquad v_\nu(x)=\frac{e^{-\pi x}\Gamma_4(\nu,x)}{16\pi^3}.
\tag{1}
\]

\[
W_M(q|p)=\frac{m'}{2\pi a}V_\nu(a|b)e^{-\sigma^2\delta^2},
\qquad w(z)=\frac1{2\pi}v_\nu\left(\frac z{2\pi}\right)e^{-\sigma^2(z/(2\pi))^2}.
\tag{2}
\]

v、W、w 均为未归一化正谱权。定义比值 R_M(p,q)=W_M(q|p)/w(q−p)，以及

\[
T_\nu(S)=\frac{S e^{\pi S}}{8\pi^2}\Gamma_4(\nu,S).
\tag{3}
\]

用 Γ₈=Γ₄(ν,d)Γ₄(ν,S)，精确地得到

\[
\boxed{R_M(p,q)=\frac{2m'}S\,
\tanh(\pi a)(1+e^{-2\pi b})^2\,
T_\nu(S)\,\frac{v_\nu(d)}{v_\nu(\delta)}.}
\tag{4}
\]

特别是 Jacobian 消去了原顶角中的 a；不能把前因子误写成 2a/S。σ 在该比值中消去，但在积分控制中不可省略。

## 2 全域估计：不把固定增量展开直接用于整个积分域

### 2.1 Stirling 因子

固定实 α 的模方 Gamma 渐近为

\[
|\Gamma(\alpha+iy)|^2=2\pi |y|^{2\alpha-1}e^{-\pi|y|}
[1+O_\alpha(|y|^{-2})].
\tag{5}
\]

由 log Γ 的 Stirling 展开取实部可见没有 |y|⁻¹ 项；不能只由最低阶等价式推出这一点。在 α=1/4、S→∞、ν 固定时，两模方的频率为 (S±ν)/2，因此

\[
T_\nu(S)=\frac{S}{\sqrt{S^2-\nu^2}}[1+O_\nu(S^{-2})]
=1+O_\nu(S^{-2}).
\tag{6}
\]

第一式只用于 S>ν 的渐近区域。全实轴无 Gamma 极点，所以 T 连续正；通过紧区间与渐近尾拼接，可定义有限常数

\[
B_\nu:=\sup_{S\ge\sqrt3}S^2|T_\nu(S)-1|<\infty.
\tag{7}
\]

故对所有 M≥1、p,q≥0，S≥√3M，|T−1|≤Bν/(3M²)，且 T≤1+Bν/3。本稿没有计算 Bν 的区间认证数值；以下 O 常数是理论有限常数，不是现成的机器误差条。

### 2.2 谱标签与真实质量差

设 h(m)=m−√(m²−1/4)=1/[4(m+√(m²−1/4))]。m≥1 时

\[
0<h(m)\le\frac1{4m},\qquad |h'(m)|=\frac{h(m)}{\sqrt{m^2-1/4}}
\le\frac1{2\sqrt3\,m^2}.
\]

h 递减，因此 d−δ=h(m)−h(m′)，并有同时有效的两个界

\[
\boxed{|d-\delta|\le\min\left(\frac1{4M},\frac{|\delta|}{2\sqrt3 M^2}\right).}
\tag{8}
\]

又有

\[
0<\frac{2m'}S\le\frac4{\sqrt3},\qquad
\left|\frac{2m'}S-1\right|
\le\frac{|\delta|}{\sqrt3 M}+\frac1{2\sqrt3 M^2}.
\tag{9}
\]

这些都是全域估计，允许 p、q 随 M 任意增大。

### 2.3 v 的对数导数有界

由 digamma 级数，y≥0 时

\[
\operatorname{Im}\psi(\tfrac14+iy)
=y\sum_{n=0}^\infty\frac1{(n+1/4)^2+y^2}
\le2+\frac\pi2.
\]

这里 n=0 项至多 2，n≥1 的单调和由从 0 开始的积分控制；y<0 用奇性。于是可取与 ν 无关的常数 L=4+2π，使

\[
|\partial_x\log v_\nu(x)|\le L\quad(x\in\mathbb R).
\tag{10}
\]

结合 (8)，G=v(d)/v(δ) 满足 G≤e^(L/4)，且

\[
|G-1|\le e^{L/4}L\min\left(\frac1{4M},\frac{|\delta|}{2\sqrt3 M^2}\right).
\tag{11}
\]

热因子 J=tanh(πa)(1+e⁻²πb)² 则有 0<J≤4，|J−1|≤11e⁻π√3M≤11/M²。将这些界逐项代入 (4)，得到存在只依赖 ν 的有限常数 Cν,Rν，使

\[
\boxed{0<R_M\le R_\nu,\qquad
|R_M-1|\le C_\nu\left(\frac{|\delta|}M+\frac1{M^2}\right).}
\tag{12}
\]

例如 Rν=16e^(L/4)(1+Bν/3)/√3 可用；Cν 可由 (7)、(9)、(11) 和 J 的上述界取和。常数无需依赖 P、p、q、σ 或实时间。

## 3 一跳积分速率与首项

令

\[
J_k=\int_{\mathbb R}|z/(2\pi)|^k w(z)dz<\infty\quad(k=0,1,2,\ldots).
\tag{13}
\]

Gamma 渐近与 σ>0 保证全部增量矩有限。由 (12)，对每个 p≥0，

\[
\int_0^\infty|W_M(q|p)-w(q-p)|dq
\le\frac{C_\nu}M(J_1+J_0),
\tag{14}
\]

且 ∫W_M(q|p)dq≤RνJ₀。这一速率对 p 一致。

对于固定 p,q，m,m′=M+O(1)。T−1、G−1 以及 J−1 均无 M⁻¹ 项，(9) 的前因子给

\[
\boxed{M(R_M-1)\longrightarrow\frac\delta2=\frac{q-p}{4\pi}.}
\tag{15}
\]

(15) 是点态结论；将它升级为 L¹ 需使用 (12)：

\[
M|R_M-1|\le C_\nu(1+|\delta|).
\tag{16}
\]

对任意固定正规化 P，P(p)w(q−p)(1+|δ|) 在 R₊² 上可积，其积分≤J₀+J₁。所以

\[
\left\|M P(p)[W_M(q|p)-w(q-p)]
-P(p)w(q-p)\frac\delta2\right\|_1\longrightarrow0.
\tag{17}
\]

P 只通过 ∫P=1 进入支配条件，不需要任何 p 矩。P(p)=(1+p)⁻² 就是允许但平均能量发散的例子。

## 4 两条四点谱密度及全部实时间

采用显示乘积 φ(t₄)φ(t₃)φ(t₂)φ(t₁)。相邻图取 δ₁=(q−p)/(2π)、δ₂=(r−p)/(2π)，(u,v)=(t₂−t₁,t₄−t₃)；嵌套图取 δ₁=(q−p)/(2π)、δ₂=(r−q)/(2π)，(u,v)=(t₄−t₁,t₃−t₂)。分别定义

\[
F_M^{adj}=P(p)W_M(q|p)W_M(r|p),\qquad
F_M^{nest}=P(p)W_M(q|p)W_M(r|q),
\tag{18}
\]

F∞ 将两 W 换成 w，p,q,r≥0。两通道都有

\[
F_M=F_\infty R_1R_2,\qquad
M|R_1R_2-1|\le C_\nu'(1+|\delta_1|+|\delta_2|).
\tag{19}
\]

后一估计用 |R₁R₂−1|≤|R₁−1|R₂+|R₂−1| 及 R₂≤Rν，不能把两个未受控的点态展开直接相乘。把两个增量的积分域扩到 R²，(19) 的总支配积分至多 Cν′(J₀²+2J₀J₁)，与 P 的尾部无关。因此两通道均满足

\[
\|F_M-F_\infty\|_1\le\frac{C_{\nu,\sigma}}M,
\quad
\boxed{\left\|M(F_M-F_\infty)-F_1\right\|_1\to0,
\qquad F_1=\tfrac12(\delta_1+\delta_2)F_\infty.}
\tag{20}
\]

相位 e^[iδ₁u+iδ₂v] 对全部 M 完全相同，模长为 1。令 C_M(u,v) 是相应四点，则

\[
\boxed{\sup_{u,v\in\mathbb R}|C_M-C_\infty|\le\frac{C_{\nu,\sigma}}M,
\qquad
\sup_{u,v\in\mathbb R}\left|M(C_M-C_\infty)-C_1\right|\to0,}
\tag{21}
\]

\[
C_1(u,v)=\iiint F_1(p,q,r)e^{i(\delta_1u+\delta_2v)}dp\,dq\,dr.
\tag{22}
\]

等价于 C_M=C∞+C₁/M+o(M⁻¹)，余项对全部实时间一致。二点由 (17) 得同样结论，系数为单增量 δ/2 的积分。该一致性不是非零下界，因此不保证相对误差一致。

速率常数在所有正规化 P 中统一，所以即使比较随 M 变化的 P_M，仍可对使用**同一个 P_M**的有限质量和极限表达应用 O(M⁻¹) 界。但 (17)、(20) 的固定首系数极限要求 P 固定；它们不能无条件用于向 p≈cM 漂移的制备族。

## 5 明确的首阶四点修正及可测任务

令 E=−2πδ₁、F=−2πδ₂，dμ(E)=vν(−E/(2π))e⁻(σE/(2π))²dE/(2π)，K_P(x)=∫x∞P(p)dp。两通道阈值分别为

\[
Q_{adj}(E,F)=\max(0,E,F),\qquad
Q_{nest}(E,F)=\max(0,E,E+F).
\tag{23}
\]

于是首阶系数明确为

\[
\boxed{C_1^{j}(u,v)=-\frac1{4\pi}\iint(E+F)\,
K_P(Q_j(E,F))e^{-i(Eu+Fv)/(2\pi)}d\mu(E)d\mu(F),
\quad j=adj,nest.}
\tag{24}
\]

它不是拟合参数，也没有额外 HH 归一化常数。两点对应单 E、单阈值 E₊。全部增量矩有限允许对 u,v 微分，故还有独立可检查的表示

\[
\boxed{C_1^j(u,v)=\frac1{2i}(\partial_u+\partial_v)C_\infty^j(u,v).}
\tag{25}
\]

这里 u,v 是两条收缩对应的独立时差；不能将 (25) 改写成共同时间平移导数，后者为零。

物理上，AABB 与 ABBA 是两独立同质量自由场共享量子观察者的完整混合四点。对不同 M 的同一 P、ν、σ、Q 制备，可比较 M[C_M−C∞] 与 (24)，或比较两质量值 M、λM：

\[
\frac{M[C_M-C_{\lambda M}]}{1-\lambda^{-1}}\longrightarrow C_1
\quad(\lambda>1\text{ 固定}),
\tag{26}
\]

收敛仍对全部实时间绝对一致。零时系数实，非零时可为复数；其厄米和反厄米组合分别给关联任务的实、虚响应。这里“可测”指本模型的指定关联任务，并非已有实验方案。若四点整体除固定 Zνσ²，C₁也除同一个常数；若改用随 M 的 C₂(M)²归一化，首项须同时展开分母，不能直接沿用 (24)。

## 6 什么时候可以进一步写 O(M⁻²) 余项

从精确前因子直接算

\[
\frac{2m'}S-1-\frac\delta{2m}
=\frac{\delta[-\delta+h(m)+h(m')]}{2mS}
+\frac{h(m)+h(m')}S.
\]

结合第 2 节各因子，得到全域局部质量展开

\[
\left|R_M-1-\frac\delta{2m_p}\right|
\le\frac{C_\nu''(1+\delta^2)}{M^2}.
\tag{27}
\]

换成固定基准 M 会多出

\[
\left|\frac\delta{2m_p}-\frac\delta{2M}\right|
\le\frac{|\delta|p}{4\pi M^2}.
\tag{28}
\]

所以若 ∫pP(p)dp<∞，(27)—(28) 可在两点积分中控制一阶展开余项。四点中用有界 R 展开两个因子；嵌套第二跳的初态钟能量 q≤p+2π|δ₁|，同一 P 一阶矩和有限增量矩仍足够。于是

\[
\boxed{\sup_{\mathbf t}\left|C_M-C_\infty-\frac{C_1}M\right|
\le\frac{C_{\nu,\sigma}(1+\int pP(p)dp)}{M^2}.}
\tag{29}
\]

这是充分条件，不宣称有限 P 一阶矩是必要条件。任意固定正规化 P 已有 (21) 的 o(M⁻¹)，不能仅据此宣称 O(M⁻²)；本稿没有为所有无限矩 P 给出反向下界或最优余项速率。

## 7 证明、核验与边界

核心全域比值、digamma 界、两增量支配、任意固定 P 无矩条件、首项负号、时差导数及一阶 P 矩充分条件已经由另一研究代理独立复核。源公式、当前定理与数值趋势分开：定理依赖上述全域估计，数值只能检查实现、首项符号和选定参数，不能认证全部实时间或所有重尾波包。旧绑定代码和无速率证明均未修改。

本结果只处理两点与两种非交叉混合场四点。ABAB 的 6j 极限仍有独立的奇异结构；闭场乘积域仍遵循此前限制。不得由 (21) 推断完整 OTOC、τ≈log M 的交叉反冲、有限 G 或静态补丁完整代数的收敛。ν、σ 固定；σ→0、ν随M增长、P随M漂移的首项均是不同问题。
### 7.1 独立实现诊断

新增 [validate_clock_rate.py](../code/validate_clock_rate.py)，结果分别保存于 [study.json](../results/stage3/clock-rate/study.json) 与 [validation.json](../results/stage3/clock-rate/validation.json)。58 项检查全部通过：60→85 位精度与精确比值分解；固定增量首系数；大标签全域分解；两点和两种四点在零时及非零时的首修正；质量 M=10、20、40、80、160；求积加阶/截止扩大及独立时差有限差分。未把点态样本或曲线收敛当成定理证明。

标准双峰包、ν=1、σ=.5 的零时首系数为：

| 任务 | C∞ | C₁ |
|---|---:|---:|
| C₂ | 0.06000295441548 | −0.02579858544693 |
| AABB | 0.005964890741940 | −0.005514181603426 |
| ABBA | 0.001154650370282 | −0.0005760666263351 |

M=160 时去掉 C₁/M 后的绝对余项依次约 1.76×10⁻⁶、4.16×10⁻⁷、5.55×10⁻⁸；这些是指定紧支撑包的数值，不是对所有 P 的常数。验证输出绑定代码、导入模块与本笔记 SHA-256，并记录依赖版本。精度细化与截止扩展仅为经验检查，没有给 Bν、谱尾或全部实时间做区间认证。