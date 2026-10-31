# 纯反冲极限下的合并顶角

2026-10-30；作者自查证明草稿，没有独立认证或优先权声明。
[四页PDF](../output/pdf/de-sitter-recoil-vertex.pdf)、[LaTeX](../manuscript/recoil-vertex/main.tex)、[原文核查](../literature/stage4-recoil-vertex-sources.md)。承接[Gram与横向涨落](stage4-gram-prefactor.md)。

## 对象与精确约分

取三维、de Sitter半径1、CLP固定拓扑四面体图。并未定义新的物理态或代数；该图和CSTY固定平均质量的物理混合表示仍有归一化及轮廓接口待补。定义Q=b+b⁻¹、pₐ=iPₐ、xσ=Q/2+Σσₐpₐ，E为三个符号乘积为+1的四组符号，O=−E。CLP (3.17)、(3.46)—(3.47)的C、V给精确恒等式

\[
U_b=\frac{V^2}{C}=K_b^{-1}
\frac{\prod_a\Gamma_b(Q+2p_a)\Gamma_b(Q-2p_a)}
{\prod_{\sigma\in O}\Gamma_b(x_\sigma)^2},
\qquad K_b=\frac{\Gamma_b(2Q)}{\sqrt2\Gamma_b(Q)^3}.
\]

这是先做S_b(x)=Γ_b(x)/Γ_b(Q−x)的精确约分，没有丢掉order-one项。实Euclidean分支的双手征合并顶角为T_b=|U_b|。完整图的内部积分表示是Z_tet=(∏四顶角T_b)|J_b|²，其中J_b为原(3.46)的内部double-sine积分。

## 固定质量下的纯反冲极限

b²=i/κ，pₐ=−1/(2b)−ibmₐ/2，L=b⁻²=−iκ；先固定正实三质量，κ→∞。记S=Σmₐ、z₀=(1+iS)/2、zₐ=[1−i(S−2mₐ)]/2。对2/b+bz使用两次精确逆b移位：

\[
\Gamma_b(2/b+bz)=\frac{2\pi b^{1-L-2z}}{\Gamma(z)\Gamma(L+z)}\Gamma_b(bz).
\]

剩余双Gamma平衡商E_b含六个1±imₐ、系数+1，四个zⱼ、系数−2。这组带符号参数的零、一、二次矩为−2、2、4，均与质量无关。

令fκ(z)=logΓ_b(bz)−logΓ(z)。由Barnes乘积/双zeta定义的三阶导数，

\[
f_\kappa'''(z)=-2\sum_{q\ge1,n\ge0}(z+n-iq\kappa)^{-3}.
\]

在Re z≥1/2、|Im z|≤R、κ≥2R上，绝对值≤Kκ=16ζ(3)/κ³+8ζ(2)/κ²。以z=1作二阶Taylor展开，任意质量三元组与固定参考三元组之间的多项式项全部抵消；log误差≤Kκ/6乘两组三次绝对矩Σ|c||z−1|³。故不用假定实b渐近可直接一致延拓到本射线，就能证明此特定归一商的极限。

普通大Gamma比B_L=Γ(L+z₀)²/∏Γ(L+1+imₐ)~Γ(L)⁻¹L⁻²，其领先部分也与质量无关。因此

\[
\frac{U_b(\mathbf m)}{U_b(\mathbf r)}\longrightarrow
\frac{F(\mathbf m)}{F(\mathbf r)},\qquad
F(\mathbf m)=\frac{\prod_a\Gamma(1-im_a)}
{\prod_a\Gamma(\tfrac12-\tfrac i2(S-2m_a))^2}.
\]

参考质量r固定；只消去与质量无关的归一化，未固定完整物理振幅的绝对值。取实Euclidean模长，得

\[
T(\mathbf m)=\pi^{-3/2}
\sqrt{\frac{m_1m_2m_3}{\prod_a\sinh(\pi m_a)}}
\prod_a\cosh\frac{\pi(S-2m_a)}2.
\]

双Gamma余项已给显式界，但B_L联合增长质量误差未在本轮统一定量；所以仍采用先κ→∞、再大质量，不能宣布有限G精度。

## 横向区域内的一致变化

严格质量三角条件ℓₐ=S−2mₐ>0时，精确分解

\[
T=\frac{\sqrt{m_1m_2m_3}}{2\sqrt2\pi^{3/2}}\,
\frac{\prod_a(1+e^{-\pi\ell_a})}
{\prod_a(1-e^{-2\pi m_a})^{1/2}}
\equiv\frac{\sqrt{m_1m_2m_3}}{2\sqrt2\pi^{3/2}}\mathcal E.
\]

延拓实分支后的极限函数，在Reℓₐ≥d>0、Remₐ≥a>0时满足
|log E|≤3e⁻πd/(1−e⁻πd)+(3/2)e⁻²πa/(1−e⁻²πa)。这不是在复质量上重新取模。

等探针ν，四顶角为(mᵢ,mᵢ₊₁,ν)，其中m₁,₃=M+ω/2±u、m₂,₄=M−ω/2±v。四顶角乘积相对ω=u=v=0的比是

\[
\mathcal R_V=P\exp(\Delta\log\mathcal E_{\rm four}),\qquad
P=\frac{[(M+\omega/2)^2-u^2][(M-\omega/2)^2-v^2]}{M^4}.
\]

八个顶角的上述界控制log修正。取ω=iA₀y、A₀=ν²/(2M)、有界实y；u、v为实且|u|、|v|≤L_box√M。若ν→∞、ν/M→0及L_box√M=o(ν)，三角间隙至少ν−2L_box√M，而

\[
\mathcal R_V=1+O_Y(\nu^4/M^4+L_{\rm box}^2/M)
+O(e^{-\pi(\nu-2L_{\rm box}\sqrt M)}).
\]

先前ν^(4/3)≪M≪ν²层级容许这样的缓慢增长盒。这给极限顶角在盒内的一致控制；没有控制盒外尾，没有证明完整作用量的高斯余项。

CLP的已知Z_obs在固定Σmᵢ=4M下，四个观察者因子相对中心的乘积趋P²。这里只核该已知构件，未凭此规定物理四点必须含这个乘积。两点N₁₂;W的领先α依赖趋常数，也不确定未知修正或四点插入规则。

## 与上一轮的关系及下一步

必须区分两种分解：上一轮Gram是假设完整6j的一圈式；本轮V已经是那个6j的一部分。不能再把R_V直接乘到上一轮Gram答案上。当前证明的是内部积分表示中V/√C的局部平坦性；下一步须在同一表示中求J_b的退化、内部涨落、相位及轮廓，才能比较两套prefactor。这一记账要求避免潜在重复计数，不将局部结果宣传为完整匹配。

[42项检查](../results/stage4/recoil-vertex/diagnostics.json)通过，参数及依赖/哈希由[脚本](../code/recoil_vertex.py)保存。含3个符号矩、4组Gamma/双曲恒等式及界、6组大Gamma比、4组格点三阶导数界、9组复杂能量横向盒及3组观察者比。65位工作；6组Gamma比和3组观察者比95位复核。格点q截40、n用polygamma求和，余尾有显式绝对界。没有直接计算整个双Gamma商或6j，也没有数值认证一致定理。

例：M=ν^1.6、u=√M/2、v=−√M/3、y=0时，ν=100、1000、10000的顶角乘积相对1的差约2.27835×10⁻⁴、5.72322×10⁻⁶、1.43761×10⁻⁷。极小指数修正用log1p计算，其显式界与多项式误差分开保存；未用这些点代替证明。
