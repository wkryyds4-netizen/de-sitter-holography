# dS₃ 两物质扇区、物理范数与非交叉极限

2026-10-31。Alan Day / Chengheng。内部自查推导；[六页LaTeX](../manuscript/ds3-two-matter/main.tex)、[PDF](../output/pdf/de-sitter-ds3-two-matter.pdf)、[来源核查](../literature/stage4-ds3-two-matter-sources.md)。接续[群配对](stage4-ds3-group-pairing.md)，所有旧科学文件保持原字节。

## 1. 模型及通道

ℓ=1、Newton G=0、proper SO⁺(3,1)，两个可区分BD标量观察者质量标签≥s>0；两个不同的BD自由物质species，标签p,q>0；局部耦合1。物质KG归一，观察者时空δ归一。沿用ρ(b)=2b²/π，μ=b tanhπb/π，σ(b)=ρ(b)sech²πb，真实静质量h(b)=√(1+b²)，三体rigging系数K_xyz和正两点核Q_p(a,b)=μ(a)μ(b)/(8D_abp)。D及K见上轮。

全标量张量积有内部主系列(J,ν)，J∈Z。局部创建b→a+q的输入观察者b是标量，所以在(a,q)配对中其像落在J=0。**这是实际产生子空间，不是全四体空间的spin0截断。** 改成(a,p)配对要保留完整自旋；同种粒子Bose干涉也涉及此改变。

## 2. CG测度与四体规约范数

取C_aq;b的000分量为1。其标量谱测度
\[
w_{aq}(b)=\rho(b)K_{aqb},\qquad b>0.
\]
一种推导是在a⊗q分解后和b做Schur群平均：K₃=wK₂，K₂=1/ρ。另一种完全独立的常数核查用DS原始标量ρ_DS=b²/(4π⁴)，乘酉球面基中边界核
\[
|B_{000}|^2=\pi^4D_{aqb}^{-1}\prod_{n=a,q,b}\frac{\sinh\pi n}{n}.
\]
得到同一w。反射代表已选正ν，不能加第二个2。数学CG谱b>0，不施加物理观察者下界；实际创建的像才有b≥s。

令J_b=(C_aq;b†⊗1_c⊗1_p)I_bcp。J_b的四个球分量全0时等于1。四体有限K群平均绝对收敛，酉CG谱分解给标量部分
\[
\eta_4^{(0)}=\int_0^\infty db\,W(a,b,c;p,q)|J_b\rangle\langle J_b|,
\quad W=\rho(b)K_{abq}K_{bcp}.
\]
其他J仍在全rigging form中。由此，分布∫F(b)J_bdb的系数范数是∫|F|²/W db；先在紧谱测试上论证，再由直接积分完成。此处不使用δ²或形式群体积。

一次物质后的系数f₁(b,c)，第二次创建给F(a,b,c)=E_abq f₁(b,c)。于是
\[
\|A_q f_1\|_2^2
=\int_{b,c\ge s}\frac{|f_1(b,c)|^2}{K_{bcp}}\mathcal S_q(b)\,dbdc,
\quad R_p(a|b)=Q_p(a,b)/\sigma(b),\quad \mathcal S_q(b)=\int_s^\infty R_q(a|b)da.
\]
该纤维乘法映射闭，最大域由精确S加权可积性给定；紧外部质量支撑稠密。反向块为伴随，对像闭包的正交补为零。本轮没有把它当作全Fock场。

任意真空f(c)两步系数为f(c)E_bcp E_abq。对线性配对态f=sechπc，物理平方范数密度约分成
\[
\mathcal W_{pq}(a,b,c)=Q_q(a,b)Q_p(b,c)/\sigma(b).
\]
这类共享区域逆密度已有KL dS₂(6.23)–(6.25)前驱。当前常数和三标签核由本模型推得。

## 3. 未调节作用域

已有S_q(b)≤b/(2π)+c_s，c_s=1/(2π²)+1/(4π³s)。从
R_p(b|c)≤b exp[−2π(b−c)_+]/(πc)分开b≤c、b>c，令β=2π：
\[
\int_s^\infty bR_p(b|c)db\le
c^2/(3\pi)+c/(\pi\beta)+2/(\pi\beta^2)+2/(\pi\beta^3s).
\]
因此ρ(c)|f(c)|²的有限二阶质量矩足以保证两步范数有限。线性态符合条件。其非交叉两步态重叠可用W乘
exp{−it_p[h(b)−h(c)]−it_q[h(a)−h(b)]}三重积分表示。

未调节条件率在大质量下有宽的向下能量区域，不能把其同点范数直接要求等于反极点的有限QFT基准。

## 4. 真正静质量的半周期延拓

在紧质量测试上定义
\[
\widetilde A_p=e^{\pi H_{\rm out}/2}A_pe^{-\pi H_{\rm in}/2}.
\]
时间约定A_p(t)=e^(−iHout t)A_p e^(iHin t)，故这是t=iπ/2。平方范数条件核
\[
\widetilde R_p(a|b)=e^{\pi[h(a)-h(b)]}R_p(a|b)
=\frac{g(a)}{g(b)}\widehat R_p(a,b),
\quad g(x)=\frac{x e^{\pi[h(x)-x]}}{1+e^{-2\pi x}},
\]
\[
\widehat R_p(a,b)=\frac{A-B}{8\pi(A+C)(B+C)},\quad
A=\cosh\pi(a+b),\ B=\cosh\pi(a-b),\ C=\cosh\pi p.
\]
这是精确恒等式；h不能替换成标签。定义k_p(u)=[8π(coshπu+coshπp)]⁻¹，B_s=2exp[π(√(1+s²)−s)]。有
\[
0<\widehat R_p(a,b)\le k_p(a-b),\quad
g(a)/g(b)\le B_s(1+|a-b|/s),\quad
k_p(u)\le e^{-\pi|u|}/(4\pi).
\]
所以sup_b∫tildeR da≤B_s[1/(2π²)+1/(2sπ³)]，上述0→1和1→2调节块均有有界延拓。形式指数算符不必在任意输入上分别存在。不因此声称有限质量已经具有精确几何反极解释。

## 5. 条件化非交叉四点与自由场极限

真空质量概率ρ|f|²取趋向δ(c−M)的归一光滑包；连续纤维极限定义
\[
F_M(t_p,t_q)=\int_{a,b\ge s}da\,db\,
\widetilde R_q(a|b)\widetilde R_p(b|M)
e^{-it_p[h(b)-h(M)]-it_q[h(a)-h(b)]}.
\]
这是两次创建态的重叠，限定operator ordering。并非未条件化的HH期望，也不是交叉OTOC。两步g比望远镜相消为g(a)/g(M)。取u=b−M、v=a−b，固定u,v和s时M→∞：
g比→1，两个hatR→k_p(u)k_q(v)，真正能差→u,v。

全局支配B_s[1+(|u|+|v|)/s]k_pk_q可积，阈值排除域逐点消失，故紧实时间范围一致
\[
F_M\longrightarrow L_p(t_p)L_q(t_q),\quad
L_p(t)=\frac{\sin(pt)}{4\pi\sinh(\pi p)\sinh t},\quad
L_p(0)=\frac p{4\pi\sinh\pi p}.
\]
零时恢复pp′/[16π²sinhπp sinhπp′]（这里p′=q）。本轮实际检验了所需独立自由探针反极断开基准，适用条件包含半周期延拓、固定s、外部质量条件化、大M和不同species。不能省略这些条件，更不包括移动阈值或未延拓同点。

## 6. 明确误差界及相位首项

M≥2s，T=max(|t_p|,|t_q|)，记
d_M=2/M+2π/M²+2πe^(−πM)，C_M=exp(Md_M/2)，
ε_M=2(2+coshπp+coshπq)e^(−πM)。定义
\[
\mathcal T_s(R)=\frac{e^{-\pi R}}{4\pi^2}
\left[(1+B_s)(R/\pi+1/\pi^2)
+\frac{B_s}{s}(R^2/\pi+2R/\pi^2+2/\pi^3)\right].
\]
则
\[
|F_M-L_pL_q|\le
\frac{C_M(d_M+2T/M^2)}{2\pi^5}
+\epsilon_M L_p(0)L_q(0)+\mathcal T_s(M/2).
\]
证明：在r=|u|+|v|≤M/2菱形，全部质量≥M/2，log g导数≤d_M，g比偏差≤C_Md_Mr；两hatR相对k的亏损和≤ε_M。由1−h′≤1/(2x²)，相位误差≤2Tr/M²。再用∫r k_pk_q≤1/(2π⁵)。外部实际核加极限核以[1+B_s(1+r/s)]e^(−πr)/(16π²)控制，四象限积分给T_s。

这是保守的全积分绝对O(1/M)界；有限M时可能远大于信号。相对界仅可在L_pL_q不为零且远离零点时除以其模，不声称已认证M=32的百分比精度。

进一步在内菱形Taylor展开g比=1+(u+v)/M+O((r+r²)/M²)，hatR修正指数小，真实能差相位误差O_T(r/M²)，外部含线性减项仍指数小。所有多项式矩有限，故
\[
F_M=L_pL_q+\frac iM(L_p'L_q+L_pL_q')+O_T(M^{-2}).
\]
首项是虚相位。保存数据中的second_formal只作探索，未作二阶定理或验证阈值。

## 7. 独立核查和误差等级

[主诊断](../results/stage4/ds3-two-matter/diagnostics.json)56项，[非球面](../results/stage4/ds3-two-matter/angular.json)10项，均通过。源码在[实现](../code/ds3_two_matter.py)及两个独立驱动中。群/顶角50位，时间NumPy双精度64/96节点；保存实际复值、Python/mpmath/NumPy版本及源哈希。

三组DS球周期规范比对；两组(a,q)、r=0,.7,2的球函数乘积谱重建，b截到32。三组四球函数积分分别做直接径向、谱积分和闭式：
\[
\int_0^\infty \sinh^2r\prod_{j=1}^4\phi_{n_j}(r)dr
=\frac1{16\prod n_j}\sum_{\epsilon_j=\pm1}(\prod\epsilon_j)
f(\sum\epsilon_j n_j),\quad f(k)=1-\frac{\pi k}{2}\coth\frac{\pi k}{2},\ f(0)=0.
\]
由四正弦展开和∫(cos kr−1)/sinh²r dr=f(k)推得；包含全标签1的零频退化样本。径向T=40尾≤2/[∏n(expm1(2T))]；谱L=18尾≤∏(sinhπn/n)e^(−2πL)。

非球面取(a,q)的l=1旋转singlet，I(S₁)=√3(1+a²+q²−b²)/[2√((1+a²)(1+q²))]。六个实际局部顶角检查此多项式；两组直接boost四体群积分对谱积分W I(S₁)，含各自解析尾。证明依靠完整CG酉理论，数值不是全自旋穷举。

时间样本s=.2、p=.6、q=1.2，|u|,|v|≤14，精确切开阈值角点，外盒误差≤T_.2(14)=1.74578e−17。盒内求积未区间认证。

| M | F_M(0,0) | 相对极限误差 | (.5,1)复数相对误差 |
|---|---:|---:|---:|
|4|.00006571136262|.0050433|.16709|
|8|.00006540451475|.00035012|.10797|
|16|.00006538310306|.000022634|.060391|
|32|.00006538171794|.0000014491|.031824|

零时极限.0000653816231972；(.5,1)减去首相位后残差4.44e−6、1.11e−6、2.79e−7、6.98e−8。零时很小的误差不代表完整时间关联同样精确。

## 8. 下一步

交叉物理四点需要(a,q)与(a,p)完整酉重耦合，包括自旋、反射、相位、真实作用域及谱尾。先构建正确物理Racah矩阵元和独立内积检查，不能再将旧正核H当完整酉核。两物质同species的Bose交叉、全N域、复合gluing、有限κ、有限G及完整全息目标未完成。
