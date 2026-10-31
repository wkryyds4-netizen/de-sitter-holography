# K₀ 径向散射与轨道绝对范数：O2/O3 的独立证明

2026-10-29。本稿仅处理 [共同域候选稿](stage3-fp-orbit-core-candidate.md) 的径向谱步骤，采用 X=dS₂、G=SO⁺(2,1)、L=1。它不使用 KL 正负频模式之间未经核对的 cross Gram，也不宣称物理像稠密或完整 N=2 酉桥完成。

**证明范围。** 对任意 K=SO(2) 不变的 f∈L²(X)，若 Q(t)=√(2π cosh t)f(t) 满足 (1+|t|)Q∈L¹(ℝ)，则 A_f(g)=〈f,U_X(g)f〉 属于 L²(G)，且有精确谱范数与定量界

\[
\|A_f\|_2^2=\int_0^\infty\frac{4\pi^2}{s\tanh\pi s}q_f(s)^2ds,
\qquad
\|A_f\|_2\le C\|(1+|t|)Q\|_1\|Q\|_2.
\tag{R1}
\]

其中 q_f 是下述酉 ds 双通道谱变换的平方模之和。C 只依赖固定径向势和 Haar 约定。真实 observer 的紧质量包是否满足该位置矩条件，应由对应包尾定理单独提供；原模型固定的 observer 下限 s_*>0 不给 |u_n| 的径向谱提供同一下限。

## 使用的已知数学输入

- Germain–Pusateri–Rousset，*The nonlinear Schrödinger equation with a potential*，Ann. IHP C **35** (2018), 1477–1530：[出版社全文](https://ems.press/content/serial-article-files/15932)，§2.1—2.3、(2.13)—(2.14)、Theorem 2.3（印刷页1487）。该定理明确给无束缚态短程实势的 distorted Fourier 酉变换、逆变换和 dk 测度。此处仅使用线性散射部分。
- Goldberg–Schlag，[math/0306108v2](https://arxiv.org/pdf/math/0306108v2)，§2：Jost 解、Wronskian 与零能共振的关系。下文对本势另给直接能量及 Volterra 论证。
- Marolf–Yu，[2610.24073v1](https://arxiv.org/html/2610.24073v1)，§3.1.2 的 SO⁺(2,1) 酉不可约表示约定，(3.2)、(3.12)、(3.19)、(3.24) 的 Haar、单位 K₀ 球函数与 Schur 常数。使用标准酉直积分分解，不声称这是本项目新结果。
- [DLMF 14.3.1](https://dlmf.nist.gov/14.3.E1)、[15.8.4](https://dlmf.nist.gov/15.8.E4) 用于显式 Jost 解与透射幅的独立校准；[14.20.11—12](https://dlmf.nist.gov/14.20.vi) 的 μ=0 Mehler–Fock 变换再次核对球函数权。

## 1. 位置测度与自伴径向算子

沿用位置测度

\[
dx=\cosh t\,dt\,d\theta,\qquad 0\le\theta<2\pi.
\]

因此，K 不变函数 f(t) 的酉径向变换是

\[
\mathcal Rf(t)=Q(t)=\sqrt{2\pi\cosh t}\,f(t),
\qquad \|Q\|_{L^2(dt)}=\|f\|_{L^2(X)}.
\tag{R2}
\]

若先写 f(t,θ)=r(t)/√(2π)，则等价地 Q=√cosh(t)r；两个约定不能混用。候选稿只写 √cosh 的形式不影响微分算子，但在 O3 定义 q 时必须补上这一角向常数。

位置 Casimir 采用正质量平方约定

\[
C_X=\square_X=-\partial_t^2-\tanh t\,\partial_t
+\operatorname{sech}^2t\,\partial_\theta^2.
\]

在第 n 个角 Fourier 模上，直接微分给

\[
\mathcal R_n C_X\mathcal R_n^{-1}
=-\partial_t^2+\frac14+\left(\frac14-n^2\right)\operatorname{sech}^2t.
\tag{R3}
\]

尤其

\[
H_0:=\mathcal R C_X\mathcal R^{-1}-\frac14
=-\partial_t^2+V_0(t),\qquad V_0(t)=\frac14\operatorname{sech}^2t.
\tag{R4}
\]

V₀为有界光滑实势，故 H₀ 在 H²(ℝ) 自伴，Cc∞(ℝ)为本质自伴核心。每个固定 n 亦如此；有限 Fourier 模的核心直和给全 Casimir 的自伴闭包。Cc∞(X) 在位置群作用下保持，且其上 Casimir 与 U_X(g) 对易，因而自伴闭包与群作用强对易。这允许将径向 Casimir 谱与群表示谱对应，而不必先假设 KL 模式的完备性。

## 2. 无束缚态及无零能共振

二次型 ∫(|Q′|²+V₀|Q|²)≥0 排除负本征值；零能 L² 本征函数由同一二次型只能为零。V₀指数衰减，标准一维短程谱定理排除正能嵌入本征值及奇异连续谱，正谱纯绝对连续。

零能共振的排除还可直接做：设 h 是有界解 −h″+V₀h=0。因为 V₀h∈L¹，h′在两端有极限；h有界迫使两个极限都为0。于是

\[
h'(t)=-\int_t^\infty V_0(u)h(u)du\quad(t\to+\infty)
\]

及负端的对应式给 h′=O(e^{-2|t|})。乘以 h 的共轭、积分并取实部，边界项消失，故

\[
\int_{\mathbb R}(|h'|^2+V_0|h|^2)dt=0.
\]

V₀严格正，故 h=0，不存在非零有界零能解。这同时排除零能共振；不是仅凭“势非负”略去边界项。

## 3. Jost 约定、精确透射幅与分支

用 F_±(t,s) 表示 Jost 解，s>0：

\[
(-\partial_t^2+V_0)F_\pm=s^2F_\pm,
\quad F_+(t,s)\sim e^{ist}\ (t\to+\infty),
\quad F_-(t,s)\sim e^{-ist}\ (t\to-\infty).
\]

固定 Wronskian 的顺序为

\[
W(s)=F_-F_+'-F_-'F_+,\qquad T(s)=\frac{2is}{W(s)}.
\tag{R5}
\]

若使用相反 W 排序，T 公式应补负号。本势偶对称，F_-(t,s)=F_+(-t,s)。直接代入微分方程或用 DLMF 的 Ferrers 超几何表示，可得

\[
F_+(t,s)=e^{ist}\,{}_2F_1\left(\frac12,\frac12;1-is;z(t)\right),
\qquad z(t)=\frac{1-\tanh t}{2}\in(0,1).
\tag{R6}
\]

所有超几何变量均在实区间(0,1)，连接公式中的 log(1−z)取实值；本步骤不存在虚轴 Legendre 边界选择问题。z→1 时的连接公式给

\[
F_+(t,s)\sim a(s)e^{ist}+b(s)e^{-ist}\quad(t\to-\infty),
\]

\[
a(s)=\frac{\Gamma(1-is)\Gamma(-is)}{\Gamma(1/2-is)^2},
\qquad b(s)=\frac{\Gamma(1-is)\Gamma(is)}{\Gamma(1/2)^2}.
\]

这里 1/T=a、R_-/T=b，因此

\[
\boxed{T(s)=\frac{\Gamma(1/2-is)^2}{\Gamma(1-is)\Gamma(-is)},
\qquad |T(s)|^2=\tanh^2(\pi s).}
\tag{R7}
\]

特别 T(s)=−iπs+O(s²)，W(0)=−2/π。这给上一节的非共振结论一个独立显式校准。一般实势的通量关系 |T|²+|R_±|²=1 在本约定下成立。

## 4. 全位置、全能量的模式界

令

\[
e_\pm(t,s)=(2\pi)^{-1/2}T(s)F_\pm(t,s),\qquad s>0.
\tag{R8}
\]

它们是下一节采用的 ds 归一化左右入射模式。证明存在固定 C，使

\[
\boxed{|e_\pm(t,s)|\le C\min\{1,s(1+|t|)\}.}
\tag{R9}
\]

先证低能量的线性位置界。F₊满足 Volterra 方程

\[
F_+(t,s)=e^{ist}+\int_t^\infty\frac{\sin(s(u-t))}{s}
V_0(u)F_+(u,s)du.
\]

在 t≥0，用 |sin(sv)|/s≤v 及 ∫₀∞uV₀(u)du<∞，Volterra/Gronwall 给 |F₊(t,s)|≤C，常数对实 s≥0 一致。微分一次后，在0≤s≤1有 |F₊′(0,s)|≤C。向 t≤0 延伸 Cauchy 解，并令 h_s(t)=|F₊(t,s)|/(1+|t|)，得到

\[
h_s(t)\le C+\int_t^0(1+|u|)V_0(u)h_s(u)du.
\]

反向 Gronwall 给 |F₊(t,s)|≤C(1+|t|)，0≤s≤1；F₋同理。由 |T(s)|≤πs 即得 R9 的低能量部分。

再证不随 t、s 增长的一致界。t≥0 时使用 |T|≤1 和 F₊的好半轴界。t≤0 时改用精确散射关系

\[
T(s)F_+(t,s)=\overline{F_-(t,s)}+R_-(s)F_-(t,s).
\]

F₋在该半轴一致有界，且 |R_-|≤1，故 |TF₊|≤C。F₋入射模式对称处理。两种界合并便得 R9。该推导同时控制低能量、大|t|的联动区域；没有先固定t后把阈值极限当作一致估计。

## 5. ds 双通道酉谱变换

Germain–Pusateri–Rousset Theorem 2.3 适用于本势：V₀∈L¹((1+|t|)dt)，并且没有束缚态。将其 k∈ℝ 变换拆为 k=s>0 与 k=−s<0，定义

\[
a_\pm(s)=\int_{\mathbb R}\overline{e_\pm(t,s)}Q(t)dt.
\tag{R10}
\]

对 Q∈L¹∩L²，积分按通常意义定义，并与 L² 延拓一致。定理给

\[
\|Q\|_2^2=\int_0^\infty(|a_+(s)|^2+|a_-(s)|^2)ds,
\]

\[
Q(t)=\int_0^\infty[e_+(t,s)a_+(s)+e_-(t,s)a_-(s)]ds
\quad\hbox{在 }L^2\hbox{ 中},
\tag{R11}
\]

以及 H₀对角化为乘 s²。**没有额外 2s、s tanh(πs) 或 K_P 因子放进这一 ds 谱测度。** 若改用能量 E=s²，才须相应更换密度。

设 M₁=∫(1+|t|)|Q(t)|dt<∞，q_f(s)=|a₊(s)|²+|a₋(s)|²。由 R9 和 R10，

\[
q_f(s)\le C M_1^2\min(s^2,1),\qquad
\int_0^\infty q_f(s)ds=\|Q\|_2^2.
\tag{R12}
\]

这一步没有假定 |u_n| 光滑，也没有假定其径向谱远离零。

## 6. 为什么相同 q_f 正好进入群球函数

令 H_sph 是位置表示中由所有 K 不变向量生成的闭 G 子空间。只讨论这个子空间；**不声称整个 L²(dS₂) 没有离散系列。**

标准 SO⁺(2,1) 酉直积分分解中，含 K₀的非平凡不可约表示只有球主系列和互补系列；每个主系列的 K₀空间为一维，取单位向量 |0〉_s。离散系列的角动量严格同号，不含 K₀；平凡表示与互补系列的 Casimir 小于1/4。R4、R11说明 K₀空间上的 Casimir 谱恰为 1/4+s²、s>0，纯绝对连续且重数二，因而 H_sph 的表示分解可以取为

\[
H_{\rm sph}\cong_G\int_0^\infty{}^\oplus(P_s\otimes\mathbb C^2)\,ds.
\tag{R13}
\]

该结论用到了标准酉表示分解和已知不可约分类；不是仅由一个微分方程就猜出完整群表示。径向散射基与群分解的 multiplicity 基可以相差可测 U(2) 矩阵，但它们对角化的是同一自伴 Casimir，使用同一酉 ds 归一化。因此对 K 不变 f，在群分解中其向量为 |0〉_s⊗v_f(s)，且

\[
\|v_f(s)\|_{\mathbb C^2}^2=q_f(s)\quad\hbox{几乎处处}.
\]

不用求出这一个可测 U(2) 矩阵；范数平方不受它影响。这正是绕过 KL 两组模式 cross Gram 的环节。

单位 K₀球函数为 MY (3.12)

\[
\varphi_s(g)={}_s\langle0|U(g)|0\rangle_s,
\qquad \varphi_s(a_\lambda)=P_{-1/2+is}(\cosh\lambda),
\quad\varphi_s(e)=1.
\]

从而

\[
\boxed{A_f(g)=\int_0^\infty q_f(s)\varphi_s(g)ds.}
\tag{R14}
\]

q_f∈L¹且|φ_s(g)|≤1，所以这个积分逐点绝对收敛，并且给连续矩阵系数。

## 7. 群 Schur 常数与合法的 L² 极限

沿用 MY (3.2) Haar：dg=dθ_L dθ_R sinhλ dλ，两个角均长2π。其单位 K₀ Schur 关系给

\[
\int_G\overline{\varphi_s(g)}\varphi_{s'}(g)dg
=K_P(s)\delta(s-s'),\qquad
K_P(s)=\frac{4\pi^2}{s\tanh\pi s}.
\tag{R15}
\]

也可用 Mehler–Fock 变换核对：μ=0 时其谱权是 s tanh(πs)，x=coshλ 给 dx=sinhλ dλ，两个角积分提供4π²。故这里既没有缺少角体积，也没有把散射变换的 ds 归一化和群 Plancherel 权混在一起。

R15是分布式关系，不能直接把任意 q 的非绝对多重积分交换。先注意

\[
\sup_{s>0}K_P(s)\min(s^2,1)<\infty.
\]

结合 R12，已经独立证明

\[
\int_0^\infty K_P(s)q_f(s)^2ds
\le C M_1^2\int_0^\infty q_f(s)ds
=C M_1^2\|Q\|_2^2<\infty.
\tag{R16}
\]

现在取 q_j∈Cc∞((0,∞))，使 q_j→q_f 同时在 L¹(ds) 和 L²(K_P ds) 中成立；先截去两个谱端点，再在紧正区间光滑逼近即可。对这些测试核，R15给

\[
\left\|\int q_j(s)\varphi_s(\cdot)ds\right\|_{L^2(G)}^2
=\int K_P(s)|q_j(s)|^2ds.
\]

右端使球函数积分列在L²(G)中Cauchy，左侧积分又因L¹逼近与|φ_s|≤1而在G上一致收敛到R14。于是两个极限相同，R1的精确等式及定量界成立。这一延拓不使用 δ(0)，也没有把分布式Schur公式直接当作任意态的Fubini定理。

## 8. 对轨道充分域的实际贡献与边界

若 ψ(x_L,x_R)=u_L(x_L)u_R(x_R)⊗v，u_L、u_R各是单磁指标包，则 f_i=|u_i| 为K不变非负函数。假设其 Q_i=√(2π cosh t)f_i 满足上面的矩条件；利用候选稿 O1 的非负 Tonelli 恒等式，有

\[
\|a_\psi\|_{L^2(Z)}^2
=\|v\|^2\int_G A_{f_L}(g)A_{f_R}(g)dg
\le C\|v\|^2 M_{1,L}M_{1,R}\|Q_L\|_2\|Q_R\|_2.
\tag{R17}
\]

因此这些产品种子进入 [FP稿](stage3-fp-two-particle-inner-product.md) 的 D_orb；有限和由轨道范数的三角不等式处理。这里物质向量只要求在其Hilbert空间内，酉群作用保持‖v‖，无需为它另设位置衰减。

本稿独立完成的是 O2、径向谱界及 O3 的规范化证明，并把所需原包假设缩成明确的加权L¹条件。原observer谱包的两端均匀估计另行验收；即使已满足，也仍然**不能从这些种子在运动学空间中稠密推出其GA像在物理范数中稠密**。全通道FP完成化、局部创建算符识别和全N下降都没有由此自动完成。