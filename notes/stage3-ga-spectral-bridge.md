# 两粒子群平均的谱系数与物理权重

2026-10-29。本稿在完整normalized-CG分解之后，构造一个明确的群平均测试域，证明其不变幅度的内积为ρ_Pl加权直积分。它是已有群Schur正交关系的应用，并未证明该完备化等于KL位置式(4.42)的全部物理空间。后者由[位置FP稿](stage3-fp-two-particle-inner-product.md)单独处理。

模型取dS₂、L=1、ν>0、proper SO⁺(2,1)，两observer标签均≥s*>0，两物质玻色子。内部分解包含所有允许principal及discrete通道，不能对每个内部标签施加observer下界。来源为[Marolf–Yu 2610.24073v1 §3](https://arxiv.org/html/2610.24073v1)及[KL附录E](https://arxiv.org/html/2411.08090v1#A5)。完整通道和normalized系数见[两粒子映射](../literature/stage3-two-particle-cg-map.md)。

## 运动学空间与明确测试域

记右observer和两个物质的空间为

\[
W_2=\left(\int_{a\ge s_*}^{\oplus}da\,P_a\right)
       \otimes\operatorname{Sym}^2 P_\nu.
\]

其总principal部分由normalized-CG合成作酉分解

\[
W_2^P\simeq\int_{d>0}^{\oplus}dd\,
(P_d\otimes\mathcal M_{2,d}).\tag{G1}
\]

物质优先的乘法空间基测度可取da dq（外顶角S/A两通道）加da及偶Δ、±的离散计数；a≥s*而q>0。记其标签统称ξ，基测度dν(ξ)。总discrete分量仍属于W₂，但不能与左principal形成不变二腿配对；这不同于从CG树删除内部discrete路径。

加入左observer后，在总principal部分用系数

\[
F_{mn}(\ell,d,\xi),\qquad
\|F\|^2_{kin}=\int d\ell\,dd\,d\nu(\xi)
\sum_{m,n}|F_{mn}(\ell,d,\xi)|^2.\tag{G2}
\]

测试域取：m,n只有有限项；ℓ、d及ξ中的连续标签为光滑紧支撑，离散ξ只有有限项，支持避开正谱端点。ℓ、d的支撑严格位于(s*,∞)内部，内部q的紧支撑避开0；特别是下界s*不是可留非零边界值的积分端点。各分支采用连续相位约定。至少要求完成ξ积分后的任意双线性谱核及其ℓ、d各阶导数可积，且有共同紧支撑。

这里m,n是**左observer及W₂总表示**的磁指标，不是原来三个W₂粒子磁指标的逐个有限截断。测试向量由酉CG合成回运动学空间，因此是普通Hilbert向量；没有将形式singlet作为运动学正规化态。

## 先质量积分再群积分的绝对收敛

对固定有限K型及紧正谱，MY(3.8)的超几何函数在λ→∞给

\[
D^{s}_{mm'}(a_\lambda)
=e^{-\lambda/2}\sum_{\epsilon=\pm1}
 A^{\epsilon}_{mm'}(s)e^{i\epsilon s\lambda}
 +O(e^{-3\lambda/2}),\tag{G3}
\]

允许余项带固定次多项式λ；紧谱上的s导数满足同样控制。由超几何连接公式可直接取得此展开，s远离0排除两频率合并。KAK两个旋转角只给有限的模一因子。

因此先按(G2)完成质量与ξ积分后，每个leading项形如

\[
e^{-\lambda}\int d\ell\,dd\,
 B_{\epsilon\epsilon'}(\ell,d)
 e^{i\lambda(\epsilon\ell+\epsilon'd)}.\tag{G4}
\]

对ℓ分部积分r次，边界项由紧支撑消失，得到O(e⁻λ(1+λ)⁻ʳ)。相位对ℓ的导数恒为±1，即使ℓ=d也没有驻相。余项为O(e⁻²λ)乘固定次多项式。Haar径向因子sinhλ只补回一个e^λ，故r≥2便可积；有限λ处连续且有紧谱控制。

这证明该谱测试域上的⟨F₁,U(g)F₂⟩属于L¹(G)。它不声称把每个固定质量矩阵元先取绝对值再积分也收敛，更不声称位置表示的逐轨道Bochner-L¹条件已经满足。

## Schur系数与不变幅度

采用MY(3.24)的Haar约定，主系列的分布Schur系数为

\[
K_P(d)=\frac{4\pi^2}{d\tanh\pi d}
      =\rho_{\rm Pl}(d)^{-1}.\tag{G5}
\]

在KL正文半相位基底，可用形式不变配对

\[
I_d=\sum_n(-1)^n|d,-n\rangle_L\otimes|d,n\rangle_W.
\]

若使用MY附录基底，须相应插入其(3.25)—(3.26)相位；以下用I_d自身定义收缩，因而不混用基底。令

\[
b_F(d,\xi)=\sum_n(-1)^nF_{-n,n}(d,d,\xi).\tag{G6}
\]

对上述光滑谱测试域，分布Schur关系给出

\[
Q_{GA}(F_1,F_2)
=\int_{d\ge s_*}dd\,d\nu(\xi)\,
 K_P(d)\overline{b_{F_1}(d,\xi)}b_{F_2}(d,\xi).\tag{G7}
\]

同一关系说明群平均后的形式不变态系数是

\[
g_F(d,\xi)=K_P(d)b_F(d,\xi),\tag{G8}
\]

而不是b_F。因此

\[
\boxed{Q_{GA}(F_1,F_2)
=\int\rho_{\rm Pl}(d)dd\,d\nu(\xi)\,
\overline{g_{F_1}(d,\xi)}g_{F_2}(d,\xi).}\tag{G9}
\]

这里没有δ(0)或群体积比值。群平均前的种子系数与不变态系数是不同坐标，K_P和ρ_Pl互为倒数正源于(G8)。

对任意光滑紧支撑的g(d,ξ)，令χ∈C∞c(ℝ)、χ(0)=1，并取窄到不越过ℓ下界的支撑，设唯一非零系数

\[
F_{00}(\ell,d,\xi)=\rho_{\rm Pl}(d)g(d,\xi)\chi(\ell-d).
\tag{G10}
\]

它属于上述测试域且g_F=g。因此取零范数商并完备化，得到整个(G9)加权multiplicity直积分，而非其中一个单通道像。这是**群平均构造**中的满射；对KL位置FP构造的满射仍未证明。

## 零磁模式的独立归一化交叉检查

[上一轮Bros解析核查](stage3-vertex-analytic-audit.md)已给全部正标签的

\[
h=\int_1^\infty P_{-1/2+is_1}(u)P_{-1/2+is_2}(u)P_{-1/2+is_3}(u)du
=\frac{\Gamma_8\prod_j\cosh\pi s_j}{4\pi^5},
\qquad |C_{000}|^2=\frac{\Gamma_8\prod_j\cosh\pi s_j}{\pi^6}.
\]

KAK两角各积分2π、u=coshλ，故零磁三主系列的直接群积分为4π²h=π³|C₀₀₀|²，正好复核MY(3.36)的系数。A通道零磁值由奇性为零。这个检查只校准三腿GA的归一化，不替代全两粒子位置FP等距证明；它复用已核的解析公式，不新增数值扫描。

## 为什么普通Hilbert范数不足以传递规约极限

取g紧支撑于d的正内部，χ(0)≠0。令

\[
F^{(\varepsilon)}_{00}(\ell,d,\xi)
=\varepsilon^{-1/2}\chi((\ell-d)/\varepsilon)g(d,\xi).
\tag{G11}
\]

对充分小ε，运动学范数不随ε变化：它等于‖χ‖²₂∫|g|²dd dν。但

\[
Q_{GA}(F^{(\varepsilon)},F^{(\varepsilon)})
=\varepsilon^{-1}|\chi(0)|^2
 \int K_P(d)|g(d,\xi)|^2dd\,d\nu.\tag{G12}
\]

因此从运动学空间到群平均空间的映射无界。这个显式对角收缩例说明：运动学范数收敛、甚至运动学有界的算符，都不能单凭这些性质保证物理规约极限存在。等变分解及所需测试域必须另核。

不过若一个已构造的运动学等变映射在(G1)中为1⊗c(d)，且sup‖c(d)‖≤B，则(G9)直接给不变幅度上的算符范数≤B；这使用了同一c(d)对(G8)的作用，不能推广到未识别的局部核。关于实际Gaussian创建及其剩余兼容条件，见[创建下降审查](stage3-creation-descent-audit.md)。

## 与位置FP任务的边界

(G1)—(G12)补齐群平均侧的测试域、权重和满射。下一步仍须在KL(4.40)—(4.42)中证明这些不变幅度产生可积的物质向量值位置波函数，并与(G9)等距；或证明[轨道Fubini充分类](stage3-fp-two-particle-inner-product.md)的群平均像在这个空间中稠密且给出同一闭包。原局部模式创建还须在该共同域满足相同弱矩阵元。

目前不能用抽象等维Hilbert空间、a.e.正则表示分解中的trivial点值，或单个真空范数匹配替代上述识别。有限质量交叉数值保持原有二次型解释；本稿没有把全物理闭场问题宣布完成。
