# 顶角标量的解析重建：Legendre 轮廓、圆周三角积分与反射通道

日期：2026-10-29。observer_literature 定位 Bros 等的原始积分文献并提出轮廓路线；clock_math_audit 独立检查跳跃符号、归一化、相位及两通道选择。本稿不修改已冻结的[A.7 源疑点记录](stage3-cg-domain-investigation.md)。

**已能推进到：对全部 s₁,s₂,s₃>0，原 A.5 的 n=0 三模积分解析确定了 D.8 已指定 S 通道的标量。结合原附录 E 的两通道分解及空间反射，可以排除另一 A 通道，并在固定模式的三线性积分层面重建顶角。所得表达式恰为 literal A.7 除以 Γ(1/2+is₃)²。该结果不等于作者确认的勘误，更不解决 FP 物理 Hilbert 空间的全场定义域。**

原始来源：

1. [Kolchmeyer–Liu 2411.08090v1](https://arxiv.org/html/2411.08090v1)，(3.10)—(3.13)、(4.9)、(A.5)—(A.7)、(D.8)、(E.14)—(E.24)。
2. [Bros–Epstein–Gaudin–Moschella–Pasquier 0901.4223v1](https://arxiv.org/html/0901.4223v1)，(1)—(2)、(64)—(66)、(124)。[本地 HTML](../sources/bros-et-al-0901.4223v1.html)、[本地 PDF](../sources/bros-et-al-0901.4223v1.pdf)。
3. [NIST DLMF 14.24.1](https://dlmf.nist.gov/14.24.E1)，Legendre 主支延拓。[本地页面](../sources/dlmf-14.24.html)。

源 URL、下载时间、字节数、SHA256 见[analytic-sources.json](../results/stage3/cg-domain-source-audit/analytic-sources.json)。

## 1. 符号与零磁指标模式

设 cⱼ=cosh(πsⱼ)、rⱼ=sⱼtanh(πsⱼ)、νⱼ=−1/2+isⱼ，并以

\[
P_j(z):=P_{\nu_j}(z)={}_2F_1(1/2-is_j,1/2+is_j;1;(1-z)/2)
\]

定义支：Pⱼ在复平面的割线 (−∞,−1] 外解析，z=1 正则。A.5/A.6 在 n=0、半相位取 e^{iπ/4} 时给

\[
\Psi_0^{s_j}(t)=\frac{\sqrt{r_j}}{2\sqrt{\pi c_j}}P_j(i\sinh t).
\tag{A1}
\]

P_{−1/2+is}=P_{−1/2−is}，故正实 s 的实解析结构给 overline(Pⱼ(i sinh t))=Pⱼ(−i sinh t)。为避免与三模积分混淆，定义无归一化积分

\[
J=\int_{\mathbb R}\cosh t\,P_1(-i\sinh t)P_2(i\sinh t)P_3(i\sinh t)dt
=\frac1i\int_{-i\infty}^{i\infty}P_1(-z)P_2(z)P_3(z)dz.
\tag{A2}
\]

## 2. 轮廓及跳跃符号

对 u>1，DLMF 14.24.1 取 μ=0、其整数延拓参数为±1。两式中 Q 项相同，相减得

\[
P_{\nu}(-u-i0)-P_{\nu}(-u+i0)
=-2i\sin(\pi\nu)P_\nu(u)=2i\cosh(\pi s)P_\nu(u).
\tag{A3}
\]

竖直轮廓向上，在右半平面闭合；仅 P₁(−z) 在 [1,∞) 有割线。上唇的 P₁ 自变量是 −u−i0，下唇是 −u+i0，所以竖直积分等于上唇减下唇的积分。于是

\[
J=2c_1\,h(s_1,s_2,s_3),\qquad
h=\int_1^\infty P_1(u)P_2(u)P_3(u)du.
\tag{A4}
\]

端点和无穷远项不能略去：z=1 只有 P₁(−z) 的对数奇性，小圆为 O(ε|log ε|)；固定 sⱼ>0 时右半平面大圆的各 P 分支为 O(R⁻¹ᐟ²)，乘积 O(R⁻³ᐟ²)，弧积分 O(R⁻¹ᐟ²) 消失。正标签排除了 s=0 的退化渐近；本稿不对标签趋零给一致估计。

Bros (1)—(2) 取 d=2 给

\[
h=\frac{\Gamma_8(s_1,s_2,s_3)}{4\pi^2\prod_{j=1}^3\Gamma(1/2+is_j)\Gamma(1/2-is_j)}
=\frac{\Gamma_8c_1c_2c_3}{4\pi^5},
\tag{A5}
\]

其中 Γ₈是所有八个±号的 Γ[1/4+i(±s₁±s₂±s₃)/2] 乘积。包括 (A1) 的三个归一化与 θ 积分的2π，真正 D.8 左边为

\[
T_{000}=\frac{c_1\sqrt{r_1r_2r_3}}{2\sqrt\pi\sqrt{c_1c_2c_3}}h
=\frac{c_1\sqrt{r_1r_2r_3}\sqrt{c_1c_2c_3}}{8\pi^{11/2}}\Gamma_8.
\tag{A6}
\]

它严格为正实数，且对 s₂,s₃ 交换对称。

## 3. Bros 圆周三角平均映射到 E.19/E.24

先以一般 Δ₁,Δ₂,Δ₃ 书写 E.15，令 α=Δ₁+Δ₃−Δ₂、β=Δ₂+Δ₃−Δ₁、γ=Δ₁+Δ₂−Δ₃。C_{S,000} 是该正弦三角乘积的归一化三角平均，归一化为 (2π)⁻³；不是未归一化的角积分。

Bros 的 J 在 (64)—(66) 是弦长平方的平均。单位圆上 |nᵢ−nⱼ|²=4sin²[(θᵢ−θⱼ)/2]，故代入 a₁=−β/2、a₂=−α/2、a₃=−γ/2，a=a₁+a₂+a₃，得到 C_{S,000}=4^{-a}J_Bros。其 (124) 在 d=2 的前因子是 4^a/π^{3/2}，因此

\[
C_{S,000}=
\frac{\Gamma(1-(\Delta_1+\Delta_2+\Delta_3)/2)
\Gamma((1-\alpha)/2)\Gamma((1-\beta)/2)\Gamma((1-\gamma)/2)}
{\pi^{3/2}\Gamma(1-\Delta_1)\Gamma(1-\Delta_2)\Gamma(1-\Delta_3)}.
\tag{A7}
\]

本问题 Re aⱼ=−1/4，单边碰撞幂可积；Re a=−3/4>−1，三点同时碰撞也可积，因此不靠越过发散域的任意延拓来定义角平均。

按原 E.19 取 Δ=(1/2−is₃,1/2−is₂,1/2+is₁)。记

\[
G_4=\Gamma\!\left(\frac14+\frac i2(-s_1+s_2+s_3)\right)
\Gamma\!\left(\frac14+\frac i2(-s_1-s_2+s_3)\right)
\Gamma\!\left(\frac14+\frac i2(-s_1+s_2-s_3)\right)
\Gamma\!\left(\frac14+\frac i2(s_1+s_2+s_3)\right),
\]
\[
D=\Gamma(1/2-is_1)\Gamma(1/2+is_2)\Gamma(1/2+is_3),\quad
C_0=\frac{G_4}{\pi^{3/2}D},\quad C^{phys}_0=e^{-i\pi/4}C_0.
\tag{A8}
\]

最后一个相位来自 E.24 的三个 n=0 半相位，不是拟合。Gamma 无零点且这里所有实部为1/4或1/2，所以 C₀在全部正实标签上非零。

## 4. 标量重建与 literal A.7 的关系

由 T₀₀₀/overline(C₀^{phys}) 得到

\[
\boxed{\mathcal V^{s_1}_{s_2,s_3,\,rec}
=e^{-i\pi/4}\frac{c_1\sqrt{r_1r_2r_3}}{8\pi\sqrt{c_1c_2c_3}}\frac{G_4}{D}
=\frac{\mathcal V^{s_1}_{s_2,s_3,\,literal}}{\Gamma(1/2+is_3)^2}.}
\tag{A9}
\]

独立检验代数：|C₀|²=Γ₈c₁c₂c₃/π⁶，且 (A9)=e^{-iπ/4}√π c₁√(r₁r₂r₃)C₀/[8√(c₁c₂c₃)]。乘 overline(C₀^{phys})，精确得到 (A6)。所以这一步是全部正实标签的 n=0 解析等式，超出了单点数值拟合。

## 5. n=0 能确定什么，反射又补上什么

n=0 本身对 A 通道盲：C_{A,000}=0。不能只凭 (A9) 就声称独立确定了所有磁指标。

这里使用原附录 E.1 开头的表示分解陈述、E.1.1 的 S/A 两份 principal 通道，以及 E.1.3 的完备关系给定的 principal×principal→principal 两重度分解：三线性不变量在 S、A 两通道展开。对固定各正实标签及有限磁指标，A.5 的模式在两端按 O(e^{−|t|/2}) 衰减，三模积分绝对收敛。时空积分因体积元与模式协变性定义 SO⁺(2,1) 不变量。

进一步，原 (3.10)—(3.11) 给 Ψₙ(t,−θ)=Ψ₋ₙ(t,θ)。对原三模积分作 θ→−θ，得磁指标同时取反时不变。另一方面，E.15 的 I_S 在三角同时反射时不变；E.16 的 I_A 含三个正弦的符号积，在反射时变号。因此其 Fourier 系数 C_S 同时反号为偶、C_A 为奇；A.6 明确给的是 e^{iφₙ}=e^{iφ₋ₙ}；半相位沿与 (3.11) 相容的约定，对 n 与 −n 选相同分支，于是 E.24/E.25 的相位变换保持此奇偶性。

**原文脚注张力须单列。** §4.3 在 (4.42) 后、(4.43) 前的脚注33，PDF印刷页26，写道：“This may be derived using the fact that (E.15) and (E.16) are invariant under θ→−θ.” PDF 已视觉复核，见[原页截图](../results/stage3/cg-domain-source-audit/v1-page26.png)。该句显示的是纯 θ 反射，没有共轭或 PT；接下来另说时间反演要求系数实。按 E.16 明写的三个正弦符号积，同步 θᵢ→−θᵢ 会变号，故这句文字与显式公式有张力。本稿没有静默改写脚注；全 n 的陈述特指采用 (3.11)、E.15/E.16、E.24/E.25 及上述半相位分支的一致 S/A 约定重建。只凭 n=0 的标量结果不需要解决脚注中的 parity 文字。SO⁺(2,1) 的两通道分解与额外反射约束在此分开使用，也不把时间反演当作空间反射。
故在该完整的两通道 CG 分解中，物理三模积分的 A 系数必须为零。S 系数又由非零 C₀唯一确定为 (A9)。**这是“先用反射排除 A，再用 n=0 定 S 标量”；不是把 n=0 的匹配当成两通道完备性证明。** 两通道分解来自引用的表示论结果；本稿没有重新证明该分解。

若把直积分的两重度陈述先理解为标签几乎处处成立，也不由此直接断言每一点。固定磁指标时，A.5 的大时间展开在正标签紧集上一致，质量导数至多带多项式时间因子，仍有可积包络；E.15/E.16 的圆周积分在这些参数上也连续。故几乎处处的顶角恒等式可由连续性延拓到所有正实标签。这里不包括 s=0 的端点。

结论在固定模式的矩阵元/三线性分布层面成立，且相位约定严格沿用 A.5、E.24。它不是完整场的自伴性定理。

## 6. Bubble 障碍消除及仍需补的物理定义域

由 (A9)

\[
|\mathcal V^{s_1}_{s_2,s_3,rec}|^2
=\frac{c_1^2}{64\pi^5}r_1r_2r_3\Gamma_8.
\]

沿[A.7 审查](stage3-cg-domain-investigation.md)的正确创建方向，以原 E.22 收缩，得到

\[
\frac{4|\mathcal V^b_{\nu,a,rec}|^2}{\nu\tanh\pi\nu\,b\tanh\pi b}
=\frac{a\tanh\pi a\cosh^2\pi b}{16\pi^5}\Gamma_8(\nu,a,b)
=V_\nu(a|b).
\tag{A10}
\]

因此主文谱权与从模式积分重建的运动学 CG 顶角相容；Gaussian 窗的相应创建映射范数≤√𝔐。这修复的是新证明尝试中的源归一化障碍，不需要改已直接实现正文核/A.14的旧代码。

FP 物理态仍是规约后的分布式 singlet，不是运动学 Hilbert 空间的闭子空间。还需证明全物理稠密测试域、各内部表示/多重度通道上的规约积分交换，以及创建/湮灭在 FP 内积中的伴随关系。只有这些补齐，才能回到[闭平方域条件](stage3-closed-domain-next.md)；本稿不宣称它们已解决。

## 7. 数值证据与共同偏差的保留

[数值记录](stage3-vertex-source-check.md)的原始时域积分独立于 Bros 解析式。T=80 时与 (A6) 的差9.78×10⁻²⁰，小于解析双尾界8.49×10⁻¹⁹；求积阶数/精度变化约2.40×10⁻³¹。

两种4F3的CG表达式虽彼此一致到约40位，却共同偏离本稿纯Gamma参考约1.28×10⁻²⁷。故没有把二者互相符合当作真实40位准确。早期数值、反事实比较、纯Gamma补充结果分别保留；[纯Gamma参考](../results/stage3/cg-domain-source-audit/vertex-analytic-reference.json)记录该偏差。60位计算的 (A6) 与 (A9)×overline(C₀^{phys})差2.36×10⁻⁶⁴，仅作解析代数的数值核查。

本研究将 (A9) 称为依原约定的“重建顶角”；检索没有找到作者确认的勘误，未联系作者。旧 A.7 必须继续按原貌保存，不能把我们的重建悄悄写成原文已刊公式。