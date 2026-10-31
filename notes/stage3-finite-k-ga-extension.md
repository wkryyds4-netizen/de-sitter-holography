# 有限左磁指标的 GA 延拓与大小种子域闭包

2026-10-29。dS₂、L=1、G_N→0、proper SO⁺(2,1)，两个主系列observer有固定s_*>0截止，物质为Sym²Pν。本稿补齐[真实正频共同域](stage3-positive-frequency-orbit-core.md)中全部C_orb产品种子，与[产品GA密度稿](stage3-product-seed-ga-density.md)中左磁指标0的小生成域S，二者FP闭包相同这一环节。不改变最大不变分布空间、局部场识别和全N问题的既有边界。

输入为KL正文[4.20](https://arxiv.org/html/2411.08090v1)的singlet相位，MY [2610.24073v1](https://arxiv.org/html/2610.24073v1)的固定K型渐近及Schur正交，和已建立的C_orb⊂D_orb、FP=GA、S的完整multiplicity稠密性。圆周模型可由Sun [2111.04591v2 §3.2](https://arxiv.org/html/2111.04591v2)的诱导实线模型作立体投影得到；下面直接给出所需作用和界，不移植该文与KL不同的相位。

## 1. 对全部主系列磁指标的一致界

在L²(S¹,dθ/(2π))上，令τ_g(θ)=g⁻¹·θ，j_g(θ)=dτ_g/dθ>0。主系列的圆周实现取

\[
[U_s(g)F](\theta)=j_g(\theta)^{1/2+is}F(\tau_g(\theta)).
\tag{K1}
\]

复幂使用正j的实log，群合成由Jacobian链式法则验证；ReΔ=1/2使该作用酉。单位Fourier基ε_n(θ)=e^{-inθ}逐点模为1，故

\[
|D^s_{mn}(g)|\le\frac1{2\pi}\int_0^{2\pi}j_g(\theta)^{1/2}d\theta
=\Xi_0(\lambda(g)),\qquad
\Xi_0(\lambda)\le C(1+\lambda)e^{-\lambda/2}.
\tag{K2}
\]

对boost可取j_a(θ)=[coshλ+sinhλ cosθ]⁻¹；外侧两旋转不影响右端。此界对所有整数m,n、所有s>0成立，不含随磁指标或质量增长的常数。每个K型一维，转为任意单位归一化KL主系列K基只乘单位相位，因此界不变。

## 2. 任意固定左右总K的群谱交换

记W=H_R⊗Sym²Pν，按完整normalized-CG分为总principal和总discrete。取固定m,m′及左紧谱包f₁,f₂，定义

\[
L_{12}^{mm'}(g)=\int d\ell\,\overline{f_1(\ell)}f_2(\ell)D^\ell_{mm'}(g).
\]

固定K型、紧正谱上的两支渐近及按ℓ分部积分，同[谱GA稿](stage3-ga-spectral-bridge.md)G3，给

\[
|L_{12}^{mm'}(a_\lambda)|
\le C_r e^{-\lambda/2}(1+\lambda)^{-r}
+C(1+\lambda)^p e^{-3\lambda/2}.
\]

KAK角度仅贡献模一因子。取r≥3即有

\[
J_{12}^{mm'}=\int_G|L_{12}^{mm'}(g)|\Xi_0(\lambda(g))dg<\infty.
\tag{K3}
\]

令v_N∈W、w_N′∈W分别只有总磁指标N、N′，其principal谱系数为v_N(d)、w_N′(d)∈M_d。矩阵系数的principal部分是

\[
\int_0^\infty dd\,D^d_{NN'}(g)
\langle v_N(d),w_{N'}(d)\rangle.
\]

纤维内积R(d)满足‖R‖₁≤‖v_N‖‖w_N′‖。由K2—K3，

\[
\int_Gdg\int_0^\infty dd\,
|L_{12}^{mm'}(g)D^d_{NN'}(g)R(d)|
\le J_{12}^{mm'}\|v_N\|\|w_{N'}\|<\infty.
\tag{K4}
\]

所以任意非紧右principal谱的交换合法，不要求原matter向量有有限磁指标。这是小域稿P8的有限左K推广；左质量积分仍必须先做。

固定N、N′的总discrete交叉项，只可能来自同号且

\[
1\le\Delta\le\min(|N|,|N'|)
\]

的有限组DΔ±；若N、N′异号或其中一个为0，就没有这些项。每个固定矩阵元满足|D^Δ_NN′(aλ)|≤C e⁻Δλ。结合|L|≤‖overline(f₁)f₂‖₁Ξ₀和Haar增长，径向尾至多为多项式乘e^{-(Δ-1/2)λ}，故普通群积分绝对收敛。multiplicity内积只贡献有限Hilbert系数。已知principal–discrete跨系列Schur正交使每项为零；这里固定左质量的交换也绝对合法，因为ℓ支撑紧且K2一致。**为零的是最终总discrete，不能删去M_d内的内部discrete路径。**

## 3. KL正文相位下的波包Schur与幅度

全篇把W输出的主系列K基选为与左observer同一个KL正文半相位规范。原(4.20)给

\[
I_d=\sum_n(-1)^n|d,-n\rangle_L\otimes|d,n\rangle_W.
\]

因此与MY(3.24)相同Haar归一化的Schur关系为

\[
\int_Gdg\,D^\ell_{mm'}(g)D^d_{NN'}(g)
=K_P(d)\delta(\ell-d)(-1)^{N+N'}
\delta_{m+N,0}\delta_{m'+N',0}.
\tag{K5}
\]

K5按谱分布理解。先完成左包，再以右紧光滑谱函数测试，并用K4与d上的支配连续性，即得普通波包积分版本；最后用K4覆盖任意右Hilbert谱系数。不存在从固定质量非绝对积分直接交换的步骤。若改回KL附录E基，须恢复MY(3.25)—(3.26)相位；不能在附录基里照抄(-1)^n。

对总旋转匹配的种子S=|f,m〉_L⊗v_-m，和S′=|f′,m′〉_L⊗v′_-m′，得到

\[
Q_{GA}(S,S')=\int_{d\ge s_*}dd\,K_P(d)
(-1)^{m+m'}\overline{f(d)}f'(d)
\langle v^P_{-m}(d),v'^P_{-m'}(d)\rangle.
\tag{K6}
\]

所有最终discrete项为零。因此，在已经固定的H_spec=L²(ρ_Pl dd;M_d)坐标中，

\[
g_S(d)=K_P(d)(-1)^m f(d)v^P_{-m}(d),
\qquad \|g_S\|_{spec}\le B_f\|v_{-m}\|,
\quad B_f^2=\|K_P|f|^2\|_\infty<\infty.
\tag{K7}
\]

有限和的幅度必须相干相加，不能把不同左m先当作物理正交通道：

\[
g_\psi(d)=K_P(d)\sum_j(-1)^{m_j}f_j(d)v^P_{j,-m_j}(d).
\tag{K8}
\]

K6的交叉配对证明其范数及所有交叉内积确为该g的谱内积；不只是给出与小域的配对。因此不会遗漏一个与小域正交、但自身范数非零的额外分量。

## 4. 任意C_orb种子可先投总旋转零模

C_orb的单项为|f_L,m〉⊗|f_R,n_R〉⊗v，其中v是任意matter Hilbert向量。以P_rot表示总旋转零模的正规化紧群平均，则

\[
P_{rot}\psi
=|f_L,m\rangle\otimes|f_R,n_R\rangle
\otimes P^{mat}_{-m-n_R}v.
\tag{K9}
\]

右侧仍是C_orb单项：只把matter向量替换为有界磁投影后的向量。它的右总磁指标为−m，因此适用K6—K8。

位置D_orb上的untwist将U(k)作用化为群变量左平移，故a_U(k)ψ(z)=a_ψ(z)。Bochner轨道积分在紧K平均下可交换，Haar不变性给

\[
\eta\psi=\eta P_{rot}\psi.
\tag{K10}
\]

更具体地，积分绝对范数被∫_Kdk,a_ψ(z)=a_ψ(z)控制；紧群平均用总质量1。故这里既保持强轨道域，也没有把单纯运动学投影误认为未经证明的物理连续算符。有限线性组合逐项成立。

## 5. 大小生成闭包相同

每个ψ∈C_orb先由K9—K10匹配旋转，再由K8取得g_ψ∈H_spec；对任意ψ、φ∈C_orb，强域FP定理与K6给

\[
\langle\eta\phi,\eta\psi\rangle_{FP}
=Q_{GA}(\phi,\psi)
=\langle g_\phi,g_\psi\rangle_{spec}.
\tag{K11}
\]

小域S包含在C_orb中，其幅度像已由[产品GA密度稿](stage3-product-seed-ga-density.md)P14—P17证明在H_spec稠密。对任意ψ∈C_orb，选择S_j∈S使g_Sj→g_ψ，则

\[
\|\eta\psi-\eta S_j\|_{FP}^2
=\|g_\psi-g_{S_j}\|_{spec}^2\longrightarrow0.
\tag{K12}
\]

因此

\[
\boxed{
\overline{\eta\mathcal S}^{\,FP}
=\overline{\eta\mathcal C_{orb}}^{\,FP}
\cong\mathcal H_{spec}.}
\tag{K13}
\]

这里S、C_orb均指其有限线性包。K13完成大小产品种子域的闭包识别，不依赖全运动学对角迹有界。它仍不是“所有允许不变分布”的最大位置空间覆盖定理，也没有自动识别KL原(5.4)的局部场作用。