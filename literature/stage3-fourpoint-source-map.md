# 阶段三：四点源公式、质量链与可编码字典

日期：2026-10-29。唯一主来源为 Kolchmeyer–Liu，[arXiv:2411.08090v1](https://arxiv.org/html/2411.08090v1)，固定版本 2024-11-12。已读取本项目 NEXT、[量子观察者字典](stage3-quantized-observer-dictionary.md)、原文 §5、§6.3、附录 A.5 与 B.21—B.27。适用模型仍是 dS₂、L=1、严格 G_N→0、有限质量量子双观察者和自由 principal-series 物质场。

**结论：** 两个非交叉图与交叉 6j 图均有原文精确输入。非交叉图可化为已有 Q 核的两次积分；交叉图可用四个复参数 ₄F₃(1) 加 Gamma 因子计算。下面区分原文输入与本项目的换元、外包推广和时间涂抹推导。非交叉图之和不是同种物质的完整四点；采用两独立 species 后，AABB 和 ABBA 分别是只有一个非交叉图的完整四点任务。

## 1 原文定位及最容易混淆的测度

| 原文位置 | 输入 | 实现含义 |
|---|---|---|
| (5.6) | φ_R(τ)=exp(−iH_Rτ)φ_R exp(+iH_Rτ) | 每个跃迁使用精确 m(s)=√(s²+1/4)，保留原负号约定 |
| (5.13) | ⟨f₂\|P_d φ_R P_c φ_R P_b φ_R P_a φ_R\|f₁⟩ | 外质量是 d；从右往左的质量链是 d→a→b→c→d |
| (6.10) | ρ(s)=ρ_Pl(s)sech²(πs) | 这是未归一化 HH 密度，不是 ρ_Pl |
| (6.23)、(6.24) | 两个不交叉 chord 图 | 分别约束 c=a、b=d，delta 相对于裸 ds |
| (6.25) | Γ^ν_ab | 就是此前的 D_ν(a,b)，不另乘谱密度 |
| (6.26)、(A.14) | crossed kernel 和十六 Gamma 乘积 | 完整 cross 不具循环性；其中的 ℳS 具循环/反射对称性 |
| (B.21)—(B.24) | 6j 的 G、B、F 表达 | 两个 G，每个含两个 ₄F₃；B 不是 Euler beta |
| (B.25)—(B.27) | 允许不同物质质量的版本 | 本稿为可复用现有 D 的任务固定两 species 同 ν |
| (6.29)、脚注 44 | 双 species 的纯交叉 OTOC、外态微正则投影 | 原文外态投影不等于每一段均加半线 Q |
| (B.46)、§B.1.1 | 大质量极限产生 delta | 须在积分/分布意义检验，不能逐点把有限质量 cross 换成 delta |

记

\[
\rho_H(s):=\rho_{\rm Pl}(s)\operatorname{sech}^2(\pi s),\quad
\rho_{\rm Pl}(s)=\frac{s\tanh\pi s}{4\pi^2},\quad
Z_H=\int\rho_H(s)ds=\frac1{8\pi^4}.
\tag{F1}
\]

原文 ρ 即这里 ρ_H。所有下述四质量密度相对于 da db dc dd；ρ_Pl 不能再重复乘进去。P_{R,s} 是谱密度记号，不是正规化单质量本征态的投影。

外部真空物理包 f 的归一化密度为 P_f(d)=ρ_Pl(d)|f(d)|²。左观察者的质量保持不变，且真空外态两侧质量相等，因此左右端均为 d。一般 f₂、f₁ 的四点密度由未归一化 HH 图乘

\[
f_2(d)^*f_1(d)\cosh^2(\pi d).
\tag{F2}
\]

同态时该因子也等于 P_f(d)/ρ_H(d)。这来自外质量谱保持及 f_HH=sech(πd)，不额外假设经典混合。若使用归一化 HH，P_f=ρ_H/Z_H，整个原图只除一次 Z_H。d=0 处按相应连续/加权积分极限理解，不在代码中先形成 0/0。

## 2 两个非交叉图：精确质量链

沿用已有正谱核

\[
D_\nu(a,b)=\frac{a\tanh\pi a\,b\tanh\pi b}{64\pi^7}\Gamma_8(\nu,a,b),\qquad
Q_b^\nu(a)=\frac{D_\nu(a,b)}{\rho_H(b)}.
\tag{F3}
\]

Q 的积分不是 1。定义乘积从左至右为 φ₁φ₂φ₃φ₄；它的质量路径从右至左仍是 d→a→b→c→d。

原 (6.23) 是嵌套配对 (14)(23)：

\[
N_{\rm nest}^{HH}(a,b,c,d)
=\frac{D_\nu(a,b)D_\nu(a,d)}{\rho_H(a)}\delta(c-a).
\tag{F4}
\]

原 (6.24) 是相邻配对 (12)(34)：

\[
N_{\rm adj}^{HH}(a,b,c,d)
=\frac{D_\nu(b,c)D_\nu(b,a)}{\rho_H(b)}\delta(d-b).
\tag{F5}
\]

乘 (F2) 后消去 delta，一般归一化外包的两种密度化为

\[
\boxed{N_{\rm adj}:\quad P_f(d)Q_d^\nu(c)Q_d^\nu(a),\qquad b=d,}
\tag{F6}
\]

\[
\boxed{N_{\rm nest}:\quad P_f(d)Q_d^\nu(a)Q_a^\nu(b),\qquad c=a.}
\tag{F7}
\]

(F7) 第二次跃迁以中间质量 a 为条件，不能换成 Q_d(b)。即使 (F6)，同一个 d 也须在最后共同平均；不能把它拆成两个已经独立平均的二点之积。物质 Wick 定理在本模型中产生这些 chord 图，但观察者顶角与中间质量仍保留在图内。

若每次插入都压缩 Φ_R=Qφ_RQ，Q=1[H_R≥m₀]1[H_L≥m₀]，令 s₀=√(m₀²−1/4)，则外包支持 d≥s₀，各未消去中间积分也均从 s₀ 起。若只制备外包而不压缩算符，中间积分仍从 0 起。纯右链中左侧投影可移至真空外端，理由见[核审查第 6、9 节](stage3-observer-kernel-audit.md)。

## 3 左至右时间标签与四个 Gaussian 因子

本稿严格规定

\[
F(\tau_1,\tau_2,\tau_3,\tau_4)
=\langle\Psi_f|\Phi_1(\tau_1)\Phi_2(\tau_2)
\Phi_3(\tau_3)\Phi_4(\tau_4)|\Psi_f\rangle,
\tag{F8}
\]

下标按书写顺序从左至右，不按先作用次序。令 m_x=m(x)。四个跃迁的“输出减输入”分别为

\[
(\Delta_1,\Delta_2,\Delta_3,\Delta_4)
=(m_d-m_c,m_c-m_b,m_b-m_a,m_a-m_d).
\tag{F9}
\]

由原 (5.6)，总相位及各插入独立 Gaussian 涂抹为

\[
\boxed{T_{\boldsymbol\sigma,\boldsymbol\tau}
=\exp\!\left[-i\sum_{j=1}^4\Delta_j\tau_j
-\frac12\sum_{j=1}^4\sigma_j^2\Delta_j^2\right].}
\tag{F10}
\]

σ_j 是 proper-time 宽度；如果来自此前模时间宽度，必须乘 2π。四个 σ 因子不可先成对随意合并；只有实施相应图的质量 delta 后才能化简。总时间平移不改变 (F10)，因为 ∑Δ_j=0。

相邻图 b=d 的相位为

\[
e^{i(m_c-m_d)(\tau_1-\tau_2)+i(m_a-m_d)(\tau_3-\tau_4)},
\tag{F11}
\]

其涂抹指数为 −[(σ₁²+σ₂²)(m_c−m_d)²+(σ₃²+σ₄²)(m_a−m_d)²]/2。

嵌套图 c=a 的相位为

\[
e^{i(m_a-m_d)(\tau_1-\tau_4)+i(m_b-m_a)(\tau_2-\tau_3)},
\tag{F12}
\]

其涂抹指数为 −[(σ₁²+σ₄²)(m_a−m_d)²+(σ₂²+σ₃²)(m_b−m_a)²]/2。所有 σ_j=σ 时，两式分别化为 −σ²[(m_c−m_d)²+(m_a−m_d)²] 与 −σ²[(m_a−m_d)²+(m_b−m_a)²]。

## 4 两独立 species 给出可直接完成的四点任务

原 (6.29) 已使用两独立自由 species A、B，令 ⟨ϕ_Aϕ_B⟩=0，以只保留交叉收缩。采用同一原理而保留量子观察者顶角：

- 左至右 AABB 只保留相邻配对 (12)(34)。
- 左至右 ABBA 只保留嵌套配对 (14)(23)。
- 左至右 ABAB 只保留交叉配对 (13)(24)。

这是物质种类选择，不是任意删除同种四点里的图。两 species 的物质质量可相同，仍是相互独立的场；此时全部复用同一个 D_ν。裸场的独立性也不意味着穿衣后的右观察者算符可任意交换。

为便于编码，所有 σ_j=σ、两 species 同 ν 时，定义

\[
K_{\sigma,t}(y|x)=Q_x^\nu(y)
 e^{-\sigma^2[m_y-m_x]^2}e^{i(m_y-m_x)t},\qquad
C_x(t)=\int_{s_0}^{\infty}K_{\sigma,t}(y|x)dy.
\tag{F13}
\]

中间不压缩时以 0 替换 s₀。完整两种任务为

\[
\boxed{F_{AABB}=\int dd\,P_f(d)
C_d(\tau_1-\tau_2)C_d(\tau_3-\tau_4),}
\tag{F14}
\]

\[
\boxed{F_{ABBA}=\int dd\,P_f(d)\int_{s_0}^{\infty}da\,
K_{\sigma,\tau_1-\tau_4}(a|d)
C_a(\tau_2-\tau_3).}
\tag{F15}
\]

更一般质量 ν_A、ν_B 时，(F14) 左对使用 ν_A、右对使用 ν_B；(F15) 外跳 Q_d(a) 使用 ν_A，内跳 Q_a(b) 使用 ν_B。本文交叉实现暂限同 ν，不能把 (A.14) 不加更改用于不同质量。

已证明对固定 ν、σ>0，有 M_σ=sup_x∫₀∞Q_x(y)e^{−σ²(m_y−m_x)²}dy<∞。因此 (F14)、(F15) 对任意归一化 P_f 均绝对收敛，且模不超过 M_σ²；半线投影只减小这个正包络。该结论确立本谱积分任务的有限性，不自动提供完整 Hilbert 空间上所有原始无界乘积的闭扩张。

等时、同窗下，F_AABB=∫P_f(d)V_σ(d)²dd；它与 [∫P_f V_σ]² 的差恰为 V_σ 在 P_f 下的方差。这是共同质量平均的诊断，不能称作外包相位敏感性。F_ABBA 则保持两跳条件质量，不能用上述平方替代。

## 5 交叉图的可编码精确表达

### 5.1 主文公式与一般外包

原 (6.26) 给出未归一化 HH 密度

\[
X_H(a,b,c,d)=a\sinh\pi a\,b\tanh\pi b\,
c\sinh\pi c\,d\tanh\pi d\,\mathcal S^\nu(a,b,c,d).
\tag{F16}
\]

一般归一化外包的纯交叉任务为

\[
F_{ABAB}=\int dd\,P_f(d)
\int da\,db\,dc\,
\frac{X_H(a,b,c,d)}{\rho_H(d)}
T_{\boldsymbol\sigma,\boldsymbol\tau}(a,b,c,d),
\tag{F17}
\]

积分或正则化意义须另外核验。这里不把 X_H 当作概率或非负核。若保留每段半线投影，a、b、c 都必须≥s₀。未进行这项高维积分时不能把点值核验证称为完整 ABAB 算例。

### 5.2 十六 Gamma 因子的紧凑实现

令 r=(r₁,r₂,r₃,r₄)、r₅=r₁，η=(−1,−1,+1,+1)。原 (A.14) 等价于

\[
\mathcal P_\nu(r)=
\prod_{j=1}^{4}\prod_{\epsilon=\pm1}
\Gamma\!\left(\frac14+\frac i2[\eta_j\nu+
\epsilon(r_j+r_{j+1})]\right)
\Gamma\!\left(\frac14+\frac i2[\epsilon\nu+r_j-r_{j+1}]\right),
\tag{F18}
\]

\[
\boxed{\mathcal S^\nu(r)
=\frac{S(r_1,r_2,r_3,r_4,\nu,\nu)\,\mathcal P_\nu(r)}{4096\pi^{11}}.}
\tag{F19}
\]

(F18) 是十六个 Gamma，不是四份 Γ₈；不能分别取模平方或实部。只有完整 (F19) 应为实数。复 Gamma 可以用 loggamma 相加后一次指数化，但所有复相位须保留。

### 5.3 B.21 的六个参数与四个 ₄F₃

为避免和谱核 D 混淆，以下把原 (B.23) 的 B 写为 ℬ：

\[
\mathcal B(u,v)=
\frac{\sqrt\pi\,
\Gamma((1-u)/2)\Gamma((1-v)/2)\Gamma((u+v-1)/2)}
{\Gamma(u/2)\Gamma(v/2)\Gamma((2-u-v)/2)}.
\tag{F20}
\]

它对 u、v 对称，但**不是** Euler beta(u,v)。原 (B.24) 为

\[
\mathcal F(h_1,h_2,h_3,h_4,h,k)=
{}_4F_3\!\left(
\begin{matrix}
1+k-h_1-h_4,\ k+h_1-h_4,\ 2-h_3-h_4-h,\ 1+h_3-h_4-h\\
2-2h_4,\ 2+k-h_2-h_4-h,\ 1+k+h_2-h_4-h
\end{matrix};1\right).
\tag{F21}
\]

原 (B.22) 等价于

\[
\begin{aligned}
\mathcal G(h_1,h_2,h_3,h_4,h,k)={}&
\mathcal B(2-h-h_3-h_4,2h_4)
\mathcal B(1+h-h_1-h_2,1+h_1-h_4-k)\\
&\times\mathcal B(h_2+h_3-k,1+h-h_3-h_4)
\mathcal F(h_1,h_2,h_3,1-h_4,h,k)\\
+{}&\mathcal B(h_1+h_4+k-1,1-h_1+h_4-k)
\mathcal B(1+h-h_1-h_2,h_1+h_4-k)\\
&\times\mathcal B(h_2+h_3-k,h-h_3+h_4)
\mathcal F(h_1,h_2,h_3,h_4,h,k).
\end{aligned}
\tag{F22}
\]

对 r₁,r₂,r₃,r₄,ν_A,ν_B 定义

\[
\begin{aligned}
h_1&=\tfrac12+\tfrac i2(-r_2+r_4-\nu_A-\nu_B),&
h_2&=\tfrac12+\tfrac i2(r_1+r_3+\nu_A-\nu_B),\\
h_3&=\tfrac12+\tfrac i2(-r_2+r_4+\nu_A+\nu_B),&
h_4&=\tfrac12+\tfrac i2(-r_1-r_3+\nu_A-\nu_B),\\
h&=\tfrac12+\tfrac i2(r_1+r_2-r_3+r_4),&
k&=\tfrac12+\tfrac i2(-r_1+r_2+r_3+r_4).
\end{aligned}
\tag{F23}
\]

则原 (B.21) 的短表达是

\[
\boxed{S(r_1,r_2,r_3,r_4,\nu_A,\nu_B)
=\pi^{-3}\left[\mathcal G(h_1,h_2,h_3,h_4,h,k)
+\mathcal G(h_3,h_4,h_1,h_2,1-h,1-k)\right].}
\tag{F24}
\]

交叉实现用 ν_A=ν_B=ν，并将 r=(a,b,c,d) 代入 (F19)、(F16)。这已经包含 6j，不另乘任何“标准归一化 6j”或 Plancherel 密度。

### 5.4 主文与附录的标签核对

附录 B.25 的外质量是 s₁，prefactor (B.26) 对 s₁、s₃ 用 tanh、对 s₂、s₄ 用 sinh。与主文外质量 d 的同一链对照为

\[
(s_1,s_2,s_3,s_4)=(d,a,b,c).
\tag{F25}
\]

相同 species 质量时，完整 ℳS 具有循环对称性，所以主文可写 ℳS(a,b,c,d)；不能因此把未装配的 S 或 Gamma 乘积单独循环重排。原 B.25 前的文字明确要求乘 sech²(πs₁) 后才得到 HH 图，这也独立核实 (F2) 应依赖 d。

## 6 交叉实现的验证与数值边界

从 (F21) 可直接算出下参数之和减上参数之和恰等于 1。因此泛型非奇异参数的 ₄F₃(1) 级数绝对收敛，但简单截项的尾可能收敛很慢。各 G 项之间也可能严重抵消，尤其大质量时。必须记录精度细化、原始虚部与相对残差；不能先取 real 再声称实性通过。

(F20) 的某些 Gamma 在 r₁=r₃ 等超曲面有单项极点，而完整 6j 可存在有限的可去极限。首轮只取四个质量互异、两个对角差都非零的点。不能把单项极点设为零，也不能把任意小扰动的结果当作极限已证明。

最低限度的点核检查：

1. 独立比较 (F18) 的分组乘积与逐项 A.14 展开。
2. 在更高 mp.dps 下重复完整核，保留复值、收敛差和精度记录。
3. 检验 ℳS 的四次循环和反射 r→(r₄,r₃,r₂,r₁)，以及完整虚部残差。
4. 检验原 (6.27)：X_H(a,b,c,d)/X_H(d,c,b,a)=cosh(πa)cosh(πc)/[cosh(πd)cosh(πb)]。接近零核时使用乘法残差而非比值。
5. 一个循环也交换奇偶位置，故 X_H(a,b,c,d)=X_H(b,c,d,a)cosh(πa)cosh(πc)/[cosh(πb)cosh(πd)]；循环两步的 full cross 则相同。这些检查不能独立固定整体常数，常数仍来自原式。

源公式 (F16)—(F24) 已完整，没有缺失的符号或测度；尚待完成的是稳定实现、对退化参数的极限处理、以及完整 cross 波包积分的误差控制。点核通过不是高维积分绝对收敛证明。

## 7 HH、大质量与域：哪些已经能检验

归一化 HH 的直接替换是 P=ρ_H/Z_H。对非交叉图，(F6)、(F7) 还原原 (6.23)、(6.24) 后统一除 Z_H；这可作为独立常数核验。等时 AABB 还可与逐质量 V_σ² 的独立积分比较，ABBA 则与真正两跳条件核比较。

原文微正则极限以外态质量窗 Λ 附近制备，再取 Λ→∞；§B.1.1 的 delta 是分布极限，只有完成相应频谱积分与整体归一化后才与固定世界线 QFT 比较。有限 σ、有限包和每段半线 Q 是本项目任务的明确选择，不把这些选择默默替换成原文无窗的 OTOC。

对于 ABAB，还可从物理两步态给出不同于“四质量核绝对 L¹”的有限性路线：两独立 species 的 Φ_BσΦ_AσΨ_f 位于一个 A 粒子、一个 B 粒子的扇区，其范数由唯一嵌套图 (F7) 给出；对同 ν、同 σ，平方范数≤M_σ²。谱截断尾可由该图的 delta 正交结构与正核控制，构造相应 L² 极限。ABAB 是两种有序两步态的内积，Cauchy–Schwarz 因而控制这个物理 pairing。此路线仍须明确选择的截断/闭扩张，并不自动证明 (F17) 可任意换积分次序或为绝对收敛的普通四重积分。

当前最稳妥的交付层级是：完整 AABB、ABBA 波包四点任务；精确 cross 的少量非退化点值验证；把 cross 高维积分、可去极点和完整算符域问题保留为明确后续项。全部有限纯右任务仍对外包逐质量相位不敏感，研究对象是质量分布、质量链与投影，不是靠提高阶数读取同一真空包的相位。
