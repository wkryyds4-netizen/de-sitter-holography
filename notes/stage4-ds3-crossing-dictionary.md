# dS₃ 交叉重耦合：物理通道与独立 Euclidean 检验

2026-10-31。Alan Day / Chengheng。内部自查推导。接续[两物质范数](stage4-ds3-two-matter.md)，[来源](../literature/stage4-ds3-crossing-dictionary-sources.md)。

## 1. 模型与要回答的问题

ℓ=1、Newton G=0、proper SO⁺(3,1)；两标量 BD 观察者标签≥s>0；不同自由 KG 物质 species p,q>0；局部耦合1。h(x)=√(1+x²)，ρ(x)=2x²/π，K_xyz 为此前固定 Haar 概率的三体群积分。所有公式是该模型的规约物理内积，不是有限 G 动力学。

固定外部 e=(a,q,c,p)，s 通道为 (a,q)|(c,p)，t 通道按 LPRS 是 (c,q)|(a,p)。实际相反创建次序使用 (a,p)|(c,q)，其与上述 t 通道相差交换两对及内部 shadow。必须先除去 0000 分量，再利用反射等价，才可把两者视为相同物理通道。

全四体 CG 谱包含全部整数自旋；两种局部创建像分别在自己配对的 J=0 通道。所需 U_00 是完整变换的一个子块，不是独立酉变换。

## 2. 球周期、反射及物理规范

取 Δ_x=1+ix，B(x,y,z) 为 all-ket 边界三点核的球面 000 分量，即代码 boundary_ll0([x,y,z],0,(True,True,True))。负标签只表示 shadow 等价代表。若 α_i=−1/2−i(S−2ν_i)/2，A=Σα_i，

\[
B=(4\pi)^{3/2}4^A\frac{\Gamma(A+2)\prod_i\Gamma(\alpha_i+1)}{\prod_i\Gamma(A-\alpha_i+2)},\qquad
|B|^2=\frac{\pi^4}{D}\prod_i\frac{\sinh\pi\nu_i}{\nu_i}.
\]

在球面正交基中，将两个外部 l=0 投影到标量内部表示，只留下内部 l=0。因此标准 shadow 积分 Ψ 的 0000 分量恰为

\[
N_s(b)=B(a,q,b)B(c,p,-b),\quad N_t(d)=B(c,q,d)B(a,p,-d).
\]

无额外 4π：两个三点投影各留下一个归一 Y00，内部积分是其配对1。平面到单位球面的酉共形映射在内部两互补维数之和2时恰吸收面积 Jacobian。由固定通道 CG 唯一性与 0000=1，J_s=Ψ_s/N_s，J_t=Ψ_t/N_t。

反射并非 Ψ_b=Ψ_−b。写 h=(1+ib)/2、h_i=(1+ie_i)/2，标准 shadow 系数满足

\[
\frac{N_s(b)}{N_s(-b)}=
\frac{K_{1-h}^{34}}{K_{1-h}^{12}},\qquad
\frac{\Psi_b}{N_s(b)}=\frac{\Psi_{-b}}{N_s(-b)}.
\]

可由上式 Gamma 周期和两 shadow block 的系数直接证明。交换两对将 shadow 积分的两三点因子交换，并翻转内部标签，故归一后的 J 不变。

## 3. quotient 常数与物理核

KKS 的标量三点 quotient 为1/4，本项目同一未归一边界核物理范数 |B|²/K=8π³，故两种 Haar 规范的比例为32π³。这是固定 Haar 体积比例，对不同点数共同适用；没有拟合自由常数。

四体逆范数沿用 W_s=ρ(b)K_abq K_bcp；W_t=ρ(d)K_cdq K_dap。由球周期直接消去外部标签，

\[
|N_s|^2/W_s=32\pi^7/b^2.
\]

独立地，KKS shadow-square=π²/b²，μ_KKS=b²/(8π³)，其 (39) 给 (1/4)²·2π/μ_KKS=π⁴/b²。这与32π³比例吻合。正 b 已选一个反射代表，没有额外2。

令 e_b^s=√W_s J_s 为 db 的 delta 归一基，LPRS 核为 𝒥₂=〈Ψ_s,Ψ_t〉_LPRS，则

\[
\langle J_s(b),J_t(d)\rangle_{phys}
=\frac{32\pi^3\mathcal J_2(b,d;e)}{\overline{N_s(b)}N_t(d)},\qquad
U_{00}(b,d;e)=\frac{bd}{\pi^4}\frac{N_s/|N_s|}{N_t/|N_t|}\mathcal J_2(b,d;e).
\]

不能把 N_s/N_t 换成模方，不能交换 t 的两对却保留原有相位。这个字典不修复早前 H/Z 候选的四点矩形反例；使用的是实际局部 CG 链。

## 4. 独立 Euclidean 面积积分

固定点 (0,z,1,∞)。记 k_h(z)=z^h·₂F₁(h−h12,h+h34;2h;z)。物理 scalar 意味着 hbar=h，不是 conj(h)。

\[
F_s(z)=|z|^{-2(h_1+h_2)}[K_{1-h}^{34}k_h(z)k_h(\bar z)+K_h^{12}k_{1-h}(z)k_{1-h}(\bar z)].
\]

t 波用 (h3,h2,h1,h4)、内部 d 及 1−z。于是 𝒥₂=(1/4)∫_C conj(F_s)F_t d²z。这个积分没有调用旧 omega_hyper 或 omega_quad，因而独立于原 Lorentzian inversion 的 ₄F₃/Euler 乘积实现；仍共享标准三点 shadow 系数。

取基域 |z|<1、|1−z|>1，六个 Möbius 图 z,1−z,1/z,1/(1−z),z/(z−1),(z−1)/z。它是按三种距离 |z|、|1−z|、1 排序所得六个域的代表，除边界外不重叠覆盖复平面。上半域 z=r exp(iθ)，acos(r/2)<θ<π，下半因 scalar 两手征交换而给相同值，连同1/4后系数为1/2。用 r=t^P、P=4或8平滑端点，Gauss–Legendre分别积分t和角度。逐图保存复值与绝对贡献和。

## 5. 不依赖质量标签的 puncture 尾界

由绝对值作用于原始 shadow 积分，所有表示维数实部为1，定义

\[
I(z)=\int_{\mathbb C}\frac{d^2x}{|x||x-z||x-1|}.
\]

有 |F_s|≤I/|z|、|F_t|≤I/|1−z|。r=|z|≤1/4时，将x平面分为 |x|≤2r、2r<|x|≤1/2、1/2<|x|≤2、|x|>2：各贡献分别≤40π、8πlog(1/(4r))、48π、16π/7。因此 I(z)≤8π[12+log(1/r)]。第一块缩放x=ry后，以 |y|+|y−e|≥1 把1/(|y||y−e|)控制为2(1/|y|+1/|y−e|)，在半径2及3圆盘积分即可。

z→1−z给1附近同界；x=zy给 I(z)=|z|⁻¹I(1/z)，控制无穷远。删去六个图全部 r<ε 区域，其像包含在0、1半径δ=ε/(1−ε)圆盘与 |z|>1/δ 的并集中。δ≤1/4，V=12+log(1/δ)，C=128π³/3，则原始6j遗漏绝对值≤

\[
T(\epsilon)=C\{2\delta(V^2+2V+2)+\tfrac12\delta^2(V^2+V+\tfrac12)\}.
\]

这是保守解析界，证明交叉内积绝对收敛且可撤去 puncture cutoff。它不是 Gauss 网格误差界，更不是同通道 delta 范数的普通积分公式。一般同通道内积不能按此绝对收敛处理。

## 6. 实际创建态的交叉泛函

真空系数f(c)的范数∫ρ(c)|f(c)|²dc，局部000振幅 E_xyk 沿用上一轮。p先q后的通道振幅

\[
v_s(a,b,c)=\frac{f(c)E_{abq}E_{bcp}}{\sqrt{W_s}}
e^{-it_p[h(b)-h(c)]-it_q[h(a)-h(b)]}.
\]

q先p后的 t 振幅为

\[
v_t(a,d,c)=\frac{f(c)E_{adp}E_{dcq}}{\sqrt{W_t}}
e^{-iu_q[h(d)-h(c)]-iu_p[h(a)-h(d)]}.
\]

所有观察者标签≥s。两步态重叠 X=∫da dc db dd conj(v_s) U_00 v_t，按紧谱截断的 Hilbert 空间极限理解，不假定绝对四重积分。未调节时∫ρ|f|²c²有限是充分域；两步均取半周期调节后，振幅多共同因子 exp{π[h(a)−h(c)]/2}，对任意初态L²有界。

完整 Racah 变换酉，而 P_s^0 U P_t^0 为 contraction。因此 |X|≤||v_s||||v_t||。对任意所选紧谱投影，

\[
|X-X_{L}|\le\|v_s-v_{s,L}\|\,\|v_t\|
+\|v_{s,L}\|\,\|v_t-v_{t,L}\|.
\]

所需尾范数来自已知正密度 |f(c)|²ρ(c)R_q(a|b)R_p(b|c) 及其交换，调节时改用tilde R。这样无需以 U 的逐点绝对值粗暴控制多重谱尾，也没有声称 scalar 子块酉。当前尚未实际完成该四重波包积分；同species Bose干涉、全Fock域、有限κ/G仍留待后续。

## 7. 数值证据入口

[45位规范诊断](../results/stage4/ds3-crossing-dictionary/diagnostics.json)保存三组标签的33项检查：三点/四点范数、相位、反射和配对交换；全部通过。第一个测试点 e=(.4,.7,.93,.6)、b=.65、d=.57，原始𝒥₂≈190.9971396918+3.7321669491i，物理U≈.7266051675149，舍入级虚部保留在数据内。

[P=4完整积分](../results/stage4/ds3-crossing-dictionary/euclidean.json)在8、16、24、40阶的复相对差依次约1.13e−2、8.93e−5、5.13e−6、1.01e−7。[P=8加密及puncture检查](../results/stage4/ds3-crossing-dictionary/euclidean-refined.json)保留第二组标签、各图复值、相消因子与截断变化。具体最终数值汇总另存；不能用一个实部小误差替代完整复值检验，也不能以这些点代表全谱精度。


最终加密：第一组P=8的24/40阶复相对差为2.80867e−7、8.73209e−11；第二组e=(1.1,.43,.82,1.3)、b=.37、d=1.17为3.37010e−5、4.58403e−8。最终绝对和/复和模约2.171与3.662。第一组ε=.001/.0001的截断32阶相对差为5.06092e−4、9.55716e−5，它们包含网格误差，不能称为纯尾误差。解析T(1e−12)=4.37066e−6（原始核绝对单位）。

[几何及尾原函数检查](../results/stage4/ds3-crossing-dictionary/geometry.json)10项通过；以(1+|z|²)⁻²完整面积π检验图的重数/Jacobian，32阶相对误差约1.24e−40。发现初始code注释/文稿把基域写成Re z<1/2；实际代码acos(r/2)一直对应|1−z|>1。已修正文稿和注释，保留[运行时源码](../results/stage4/ds3-crossing-dictionary/executed-kernel.py.txt)及[修正记录](../results/stage4/ds3-crossing-dictionary/source-correction.json)，剔除docstring后的AST完全相同，旧数值未重写。
