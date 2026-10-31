# 两点谱校准、制备协议与归一化自由度

2026-10-31。自查推导与指定模型计算；[LaTeX](../manuscript/recoil-spectral/main.tex)、[PDF](../output/pdf/de-sitter-recoil-spectral.pdf)、[来源](../literature/stage4-recoil-spectral-sources.md)。承接[直接C₀因子](stage4-recoil-poles.md)。

## 1. 原文已有的校准，不重复当作发现

半径ℓ=1，dS₃自由主系列标量，ν>0，场质量平方1+ν²。CLP(3.38)在无亏角α=2π时，与CSTY(1.16)–(1.17)的谱核一致：
\[
k_\nu(\omega)=\frac1{8\pi[\cosh(\pi\omega)+\cosh(\pi\nu)]},\qquad
g_\nu(t)=\int_\mathbb R d\omega\,e^{-i\omega t}k_\nu(\omega)
=\frac{\sin(\nu t)}{4\pi\sinh(\pi\nu)\sinh t}.
\]
g(0)=ν/[4πsinh(πν)]，无额外1/(2π)。这是θ=π+it的反极点两侧关联函数，**不是同侧Wightman**。CSTY固定平均质量的混合表示和带阈值的迹已有这一差别；下面追踪我们的固定质量C因子及有限质量记账。

仍先固定质量κ→∞，再取大质量或时间极限。真实有限G态、内积和完整四点插入尚未确定；本轮不假定CLP的𝒩在固定质量极限已知。

## 2. 加上观察者权重后，能确定哪一部分

CLP(3.10)与α=2π(1−m/κ)给，在去掉与m无关的常数后
\[
z(m)=m^2e^{-2\pi m},\qquad m>0.
\]
上一轮直接C配对D(m₁,m₂,ν)已知。令m₁,₂=M±ω/2，A=cosh(2πM)、B=cosh(πω)、C=cosh(πν)。使用sinh乘积与cosh乘积身份，严格得到
\[
\frac{z(m_1)z(m_2)D(m_1,m_2,\nu)}{z(M)}
=16\pi^2\frac{e^{\pi\nu}\sinh\pi\nu}{\nu}\,k_\nu(\omega)R_M(\omega),
\]
\[
R_M(\omega)=\left(1-\frac{\omega^2}{4M^2}\right)
\frac{A-B}{A+C}.
\]

在上述已去常数的规范中，选择n₀(ν)=νe⁻πν/[16π²sinhπν]能使有界ω、大M时对准k。它不是原文完整𝒩的绝对值。把尚未知的剩余插入函数记为n(M,ω)，则**候选**能量核为
\[
\mathcal Q_n(m_1,m_2)=z(M)k_\nu(\omega)R_M(\omega)n(M,\omega).
\]
物理两点正性和两标签对称性建议n在实域正且对ω偶；能否由实际算符和内积实现，仍须另证。本轮n=1只表示一个明确的诊断模型，不宣布它就是完整物理两点。

对实|ω|<2M，0<R<1。若质量下界m*>0，固定M积分区域是|ω|≤D_M=2(M−m*)。n=1时对所有实t有一致绝对界
\[
\left|\int_{-D_M}^{D_M}e^{-i\omega t}k_\nu R_M\,d\omega-g_\nu(t)\right|
\le\frac{g_\nu(0)(1+\nu^2)}{12M^2}
+\frac{D_M}{4\pi(A+C)}+\frac{e^{-\pi D_M}}{2\pi^2}.
\]
证明：1−R≤ω²/(4M²)+(B+C)/(A+C)，∫ω²k=g(0)(1+ν²)/3，且k≤e⁻π|ω|/(4π)。界不在g的零点除首项，也不是全时间相对界。固定有限M的硬截断另有2k(D_M)R(D_M)sin(D_Mt)/t端点项；因此固定M晚时与先M→∞不同。

## 3. 大质量校准不能唯一固定有限质量修正

对任意λ≥0，取
\[
n_\lambda(M,\omega)=1+\lambda\frac{\omega^2}{4M^2}.
\]
它在实质量域正、偶，nλ(M,0)=1，且固定ω时nλ→1；所有λ都通过同一个大M谱校准。然而
\[
R_M n_\lambda
=1+(\lambda-1)\frac{\omega^2}{4M^2}+O(M^{-4})+O(e^{-2\pi M})
\]
在固定ν、ω下改变首修正，λ=1甚至消去该阶。这是具体非唯一性例子：即使要求零能差值及领先QFT归一一致，仍不能从现有两点校准推出唯一有限M修正。不可据n=1的−ω²/(4M²)修正宣布物理预言。有限G四点的歧义更不可能由一个谱形状自动消除。

## 4. 固定初态质量与同侧Wightman

采用CSTY的能量符号E=−m，但明确设m₁=m是初态、m₂=m−ω是末态，所以ω=E_final−E_initial，W(t)W(0)的相位是e⁻ⁱωt。定义候选条件谱
\[
w_{m,n}(\omega)=\mathcal Q_n(m,m-\omega)/z(m),\quad \omega\le m-m_*.
\]
同一对质量交换得到精确反向关系
\[
z(m)w_{m,n}(\omega)=z(m-\omega)w_{m-\omega,n}(-\omega).
\]
此关系连接**不同初态**，不是同一有限m态的KMS条件。任意对称n均满足。n=1时M=m−ω/2，两个多项式精确约成1−ω/m：
\[
w_{m,1}(\omega)=e^{\pi\omega}k_\nu(\omega)
\left(1-\frac\omega m\right)
\frac{\cosh\pi(2m-\omega)-\cosh\pi\omega}
{\cosh\pi(2m-\omega)+\cosh\pi\nu}.
\]
若n在固定ω、大m时趋1，则
\[
w_\infty(\omega)=e^{\pi\omega}k_\nu(\omega),\qquad
\frac{w_\infty(-\omega)}{w_\infty(\omega)}=e^{-2\pi\omega}.
\]
这是标准温度1/(2π)的详细平衡，非新温度。密度反向比为z(m−ω)/z(m)=(1−ω/m)²e²πω，局部熵斜率2π−2/m来自z的已知多项式；不能把这个跨初态比当成同一态的精确热性。

严格的时域关系在解析条带中是
\[
g_\nu(t+i(\pi-\epsilon))
=\int d\omega\,e^{-i\omega t}e^{(\pi-\epsilon)\omega}k_\nu(\omega),
\quad0<\epsilon<2\pi.
\]
ε↓0给同侧Wightman分布边界值，t→0有−i/[4π(t−i0)]短距奇性。Wightman正频谱趋常数，不可直接使用两侧谱的L¹尾界或无调节绝对求积；本轮只数值核条带内部，并用谱式证明详细平衡。

## 5. 全质量迹：一个可完全计算的指定模型

现在明确选择n=1，取m₁,m₂≥m₀>0。以
\[
Z_0=\int_{m_0}^\infty z(m)dm
=e^{-\beta m_0}(m_0^2/\beta+2m_0/\beta^2+2/\beta^3),\quad\beta=2\pi
\]
归一谱积分形式的迹模型，得到偶谱（尚未从引力Hilbert空间证明该积分确为算符迹）
\[
s_{m_0}(\omega)=\frac{k_\nu(\omega)}{Z_0}
\int_{m_0+|\omega|/2}^\infty dM\,z(M)R_M(\omega).
\]
它为正且偶，属于β=0的迹协议，并不等于上述Bunch–Davies条件谱。与CSTY忽略缓变质量前因子的结果相比，这里保留了m²及R的多项式。

令c₀=(m₀+1/β)/(m₀²+2m₀/β+2/β²)。去掉R的指数小部分后，质量积分恰给
\[
s^{(0)}_{m_0}(\omega)
=k_\nu(\omega)e^{-\pi|\omega|}(1+c_0|\omega|),
\quad
0\le1-\frac{s_{m_0}}{s^{(0)}_{m_0}}
\le2(1+\cosh\pi\nu)e^{-2\pi m_0}.
\]
这个相对界对全部实ω一致；它来自M≥m₀+|ω|/2，并非把ω固定后随意积分。

要算指数修正，若m₀>ν/2，令bₙ=(n+2)β、a=|ω|/2、L=m₀+a、qₙ=(−1)ⁿsinh((n+1)πν)/sinh(πν)，则绝对收敛级数为
\[
s_{m_0}(\omega)=s^{(0)}_{m_0}(\omega)
-\frac1{4\pi Z_0}\sum_{n\ge0}q_n e^{-b_nL}
\left[\frac{m_0^2+m_0|\omega|}{b_n}
+\frac{2m_0+|\omega|}{b_n^2}+\frac2{b_n^3}\right].
\]
不满足m₀>ν/2时仍可使用原质量积分，不能使用此级数。

## 6. 迹的长时间尾与可复查数值

移动下界使s在零频出现角点：s′(0+)=−kν(0)z(m₀)R_m₀(0)/(2Z₀)。在半轴上四阶导数可积，分部积分给
\[
T_{m_0}(t):=\int d\omega e^{-i\omega t}s_{m_0}(\omega)
=\frac{k_\nu(0)z(m_0)R_{m_0}(0)}{Z_0t^2}+O(t^{-4}).
\]
这是硬质量阈值和迹协议的结果，不是有限G混沌指数。CSTY已经给出相应的Cauchy卷积/长尾机制；本轮保留当前n=1模型的质量前因子与精确R。对于一般光滑、偶的n，角点系数还乘n(m₀,0)，并须额外保证可积性；故也不能不加条件地声称它是唯一物理尾。

数值把s⁽⁰⁾作实际Fourier积分，将级数每一项的exp(−λ|ω|)乘一次多项式精确变换为有理函数。不是用已知t⁻²系数生成曲线。ν=.6，m₀分别为1、3，t取8、16、32、64、128，末点与领先尾之比分别1.00028507、.99994557。

质量积分与级数直接比较，20/24项细化；Fourier采用28/40阶、45/65位、步长.1、频率截断12，最大t另扩到14。基线尾有显式指数界；截断N项级数的全时域误差界为
\[
\frac{e^{-2\beta m_0}\rho^N}{2\pi Z_0(1-e^{-2\pi\nu})(1-\rho)}
\left[\frac{2m_0^2}{b_N^2}+\frac{8m_0}{b_N^3}+\frac8{b_N^4}\right],
\quad\rho=e^{\pi\nu-\beta m_0}<1.
\]
本次解析截断界/信号最大2.05e−28；求积加阶最大约2.61e−43。两者不等于完整区间误差认证：求积舍入和离散误差仅有经验核查。保留了完整复数条带积分及真实小实信号，没有以大平台值掩盖误差。

最终[71项诊断](../results/stage4/recoil-spectral/diagnostics.json)通过。[初版59项](../results/stage4/recoil-spectral/initial-diagnostics.json)及[对应源码](../results/stage4/recoil-spectral/initial-recoil-spectral.py)保留，新增12项为解析截断界和截止扩展，未降低容差。最初一次执行因mpmath不接受字符串'.0'而在保存前中断，改为数值0后继续；不是科学检验失败。

## 7. 四点的必要验收与下一任务

在G=0、大观察者质量、两种独立自由标量W、V的明确QFT基准下，反极点OTOC的断开归一应为
\[
F_{\rm disc}=g_\nu(0)g_{\nu'}(0)
=\frac{\nu\nu'}{16\pi^2\sinh\pi\nu\sinh\pi\nu'}.
\]
这来自Wick收缩，不是某个四面体图在ω=0的值。相同标量还须保留其他Wick配对。有限M的共同质量约束也不自动因式分解。

两点校准仅给必要条件：未来四点插入规则必须先解释n及谱测度，再恢复这个基准；仍需交叉核的独立信息，不能由两点唯一确定。下一步优先核已有Virasoro/复群缝合和identity-leg极限的真正归一规则，或给出可检验的条件性候选。identity插入不是ν→0（ν=0仍对应质量平方1）；若需要Δ→0，涉及ν=±i的非主系列继续和分布极限，必须从原式论证。有限G、完整全息对偶仍开放。
