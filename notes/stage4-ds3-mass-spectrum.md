# dS₃时空谱、BD投影与观察者归一

2026-10-31。Alan Day / Chengheng。[LaTeX](../manuscript/ds3-mass-spectrum/main.tex)、[PDF](../output/pdf/de-sitter-ds3-mass-spectrum.pdf)。承接[局部顶角](stage4-ds3-local-vertex.md)，先核已有[LMW/MM/KL来源](../literature/stage4-ds3-mass-spectrum-sources.md)。以下是指定运动学模型中的自查推导；不是有限G物理Hilbert空间或完整对偶。

## 1. 全时空算符与完整谱

半径ℓ=1，G=0，dx=cosh²t dt dΩ，球谐单位归一。对F=Σf_lm(t)Y_lm，酉映射RF=(cosh t f_lm)把L²(dS₃,dx)送到⊕_lm L²(R,dt)。正质量约定的Casimir为□，
\[
R\Box R^{-1}=\bigoplus_{lm}(1+H_l),\qquad
H_l=-\partial_t^2-l(l+1)\operatorname{sech}^2t.\tag{1}
\]
每个H_l在H²(R)自伴，Cc∞为核心，因为势有界实值。全算符定义域为各径向H²向量且Σ||(1+H_l)χ_lm||²有限；有限球谐、紧径向核心给其自伴闭包。Cc∞(dS₃)是包含此核心的对称定义域，故同一闭包与保持它的proper orthochronous等距作用强对易。

完整性可以由Darboux极分解直接证明，不必把无束缚态版本的散射定理误套在这里。令D_j=∂t−j tanh t，
\[
D_j^\dagger D_j=H_{j-1}+j^2,\quad
D_jD_j^\dagger=H_j+j^2,\quad
W_j=D_j(H_{j-1}+j^2)^{-1/2}.\tag{2}
\]
归纳中H_{j−1}最低点≥−(j−1)²，逆平方根有界。W_j是等距算符，范围为ker(D_j†)的正交补；核由sech^j t张成，能量−j²。W_j intertwine前后H，所以每步只增加一个归一束缚态，并把其余谱酉搬运。由H₀的普通Fourier变换，得到：

- 每个l连续谱ν²>0有两个通道，纯绝对连续；无奇异连续谱。
- 束缚谱为−k²，k=1,…,l，各径向本征值简单。
- 零阈值有有界解P_l(tanh t)，不是L²态，无零能原子。l=0也有常数共振。

特别整个□有点谱1−k²、k≥1，每个含所有l≥k及m，故全空间中每个点的重数无限。这里不将它称作有限重数的算符离散谱，也不在未做表示分类前叫作群的离散系列。原自由物质的ν>0主系列与这些L²束缚扇区不是同一概念。

定义d_l(ν)=∏_{j=1}^l√(ν²+j²)，
\[
e_{\nu l}^{\mp}(t)=\frac{i^l}{\sqrt{2\pi}d_l(\nu)}
D_l\cdots D_1 e^{\mp i\nu t}.\tag{3}
\]
两通道ν>0各以δ(ν−ν′)归一，彼此正交。束缚模可取实的
\[
b_{lk}(t)=\frac{D_l\cdots D_{k+1}\operatorname{sech}^{k}t}
 {\sqrt{B(1/2,k)}\prod_{j=k+1}^{l}\sqrt{j^2-k^2}}.\tag{4}
\]
由(2)其范数为1。若a_±(ν)=〈e_νl^±,χ〉、c_k=〈b_lk,χ〉，则
\[
\|\chi\|^2=\int_0^\infty (|a_-|^2+|a_+|^2)d\nu+
\sum_{k=1}^l|c_k|^2,\quad
\chi=\int_0^\infty(e^-a_-+e^+a_+)d\nu+\sum b_{lk}c_k.\tag{5}
\]
逆变换在L²成立，δ核是分布。不能删掉第二通道和束缚部分后仍把结果写成全位置δ。

## 2. KG与时空Gram不是同一个内积

前稿KG归一BD模式为u_νlm=f_νlY_lm。写χ_B=cosh(t)f_νl，q=e^(−πν)，N=[2ν(1−q²)]^(−1/2)，则
\[
\chi_B=\sqrt{2\pi}N(e^- -qe^+).\tag{6}
\]
空间反极映射A:(t,n)↦(−t,−n)是L²中的自伴酉算符，保持质量，且与proper等距群对易。对固定lm，A后的径向模χ_A=(−1)^lχ_B(−t)=(−1)^lχ_B(t)*，故
\[
\chi_A=\sqrt{2\pi}N(e^+ -qe^-).\tag{7}
\]
对于复球谐，这里保持Y_lm；完整u*另含m反转，不可混淆。

由(5)—(7)，在(B,A)双通道基中
\[
\langle u_{\nu lm}^{r},u_{\nu' l'm'}^{s}\rangle_{L^2}
=G_{rs}(\nu)\delta(\nu-\nu')\delta_{ll'}\delta_{mm'},
\quad G=\frac\pi\nu
\begin{pmatrix}\coth\pi\nu&-\operatorname{csch}\pi\nu\\
-\operatorname{csch}\pi\nu&\coth\pi\nu\end{pmatrix}.\tag{8}
\]
两个本征值为(π/ν)tanh(πν/2)、(π/ν)coth(πν/2)，detG=(π/ν)²。
这不同于固定质量的KG签名(+1,−1)；时空Gram正定，交叉项不为零。若以原χ_B*为第二径向基，交叉项额外乘(−1)^l。

KL dS₂(4.9)—(4.11)给其特定反极模式的正交性。这里从dS₃无反射模式重推，得到非零交叉项，故不能跨维照搬(4.11)。虽然对角归一函数恰与KL形式相同，这不意味着两个谱分解相同，也不是指控KL的dS₂公式错误。

两条非正交线若都保留，完整连续投影须用G⁻¹：∫u^r(G⁻¹)_rs〈u^s,·〉dν。分别除以各自范数后直接相加不是恒等算符。

## 3. 明确的BD观察者子空间

取质量下界ν_*>0，定义
\[
\mu(\nu)=\frac{\nu\tanh\pi\nu}{\pi},\quad
v_{\nu lm}=\sqrt{\mu(\nu)}u_{\nu lm},\quad
\mathcal H_{\rm obs}=\int_{\nu_*}^\infty{}^\oplus\ell^2(l,m)\,d\nu.\tag{9}
\]
映射Ja=Σ∫v_νlm a_lm(ν)dν为等距嵌入时空L²，范数∫Σ|a|²dν。在e通道中，BD单位向量和其真实正交补为
\[
v_B=\frac{(1,-q)^T}{\sqrt{1+q^2}},\qquad
v_\perp=\frac{(q,1)^T}{\sqrt{1+q^2}},\qquad
P_B=\frac1{1+q^2}\begin{pmatrix}1&-q\\-q&q^2\end{pmatrix}.\tag{10}
\]
P_obs在连续谱为1_{ν≥ν_*}P_B，在所有束缚模为0。因此J J†=P_obs，J†J=1；它是范数1的正交投影，不是全位置恒等。

两个e通道的boost矩阵相同：沿z轴的iLie_B在相邻l上系数
\[
\beta_{j,m}(\nu)=\sqrt{\frac{(j^2-m^2)(j^2+\nu^2)}{(2j-1)(2j+1)}}.\tag{11}
\]
从(3)及球谐递推可直接得到；与MM(20)—(21)在d=2、σ=−1+iν一致。系数与通道符号无关，旋转也相同，故(10)中只依赖ν的矩阵与proper群作用对易。也可由BD子空间不变、A与群对易及二者在ν>0张成两个通道证明。质量下界是Casimir谱截断，同样保持群。

反极单位向量Av_B=(−q,1)/√(1+q²)满足
\[
\langle v_B,Av_B\rangle=-\operatorname{sech}\pi\nu,\quad
v_\perp=\frac{Av_B+\operatorname{sech}(\pi\nu)v_B}{\tanh\pi\nu}.
\tag{12}
\]
所以P_B与AP_BA不互补。ν→0时这两条线近乎平行，双基Gram条件数coth²(πν/2)发散；P_B本身仍有有界极限，μ~ν²。ν_*用于本轮模型，不借低能极限偷偷删去它。ν→∞时反极重叠指数小，μ~ν/π。

恢复半径ℓ：KG模u_ℓ=ℓ^(−1/2)u₁，时空Gram=ℓ²G，v_ℓ=ℓ^(−3/2)v₁。若用物理有效质量m=√(1+ν²)/ℓ，则dν/dm=ℓ²m/ν，δ(m−m′)模为√(ℓ²m/ν)v_ν。若用λ=ℓ²m²=1+ν²，则模为v_ν/√(2ν)。这些是谱变量转换，不是额外的群Plancherel密度。

## 4. 局部顶角与受控两点的接口

采用(9)的观察者、前稿的局部耦合1与KG归一物质腿，三标签系数变为
\[
C_\eta^{\rm obs}(\nu_o,\nu_i,\nu_p)
=\sqrt{\mu(\nu_o)\mu(\nu_i)}C_\eta^{\rm KG}.\tag{13}
\]
只重标两个观察者腿，物质场的创造湮灭CCR仍采用KG模式。ℓ恢复后整体量纲ℓ^(−1/2)，另按物理耦合恢复单位。式(13)固定了这一运动学实现的腿归一，尚未固定旧两点n：那还包含规约态、Haar体积/Schur常数、反射配对及有限κ内积。

可先给一个严格的有界算符接口。对g∈L∞(dS₃)，定义
\[
A_g=J^\dagger M_gJ,\quad\|A_g\|\le\|g\|_\infty,
\quad A_g^\dagger=A_{\bar g}.\tag{14}
\]
在位置子空间中即P M_g P。它对proper等距作用协变。中间质量求和是ν≥ν_*的dν及所有lm，不可换成无投影位置δ。对任意g,h，
\[
A_gA_h=J^\dagger M_gP M_hJ,\quad
A_{gh}-A_gA_h=J^\dagger M_g(1-P)M_hJ.\tag{15}
\]
实g时后者为[(1−P)M_gJ]†[(1−P)M_gJ]≥0。漏掉中间投影会改变有限质量二点，即便原乘法算符彼此对易。

每个固定物质球谐模式u_p有界。若先保留有限个物质模式p∈S，
O_S=Σ_p(A_{u_p}⊗a_p+A_{u_p}†⊗a_p†)在有限粒子核心本质自伴（有限个有界系数的Segal估计给解析向量）。真空二点严格为
\[
\langle\psi,0|O_S^2|\psi,0\rangle
=\sum_{p\in S}\|A_{u_p}^\dagger\psi\|^2
=\sum_{p\in S}\sum_{lm}\int_{\nu_*}^\infty
 |\langle\nu lm|A_{u_p}^\dagger|\psi\rangle|^2d\nu.\tag{16}
\]
有限物质截断不保持完整de Sitter不变性；撤去截断的域、紫外收敛和group averaging仍须证明。式(16)是运动学受控基准，不是所求物理HH两点；不把它替代最终目标。反极点大质量断开四点基准及有限G要求保留。

## 5. 诊断与继续位置

[实现](../code/ds3_mass_spectrum.py)、[验证脚本](../code/validate_ds3_mass_spectrum.py)、[原始诊断](../results/stage4/ds3-mass-spectrum/diagnostics.json)记录50位矩阵/模式/束缚谱计算，sech^p测试函数的真实连续谱积分加束缚权重、独立时积分Fourier变换，以及紧支撑质量包的双精度时空积分。质量包同时检BD范数、反极非零重叠和真实正交补，不能用单个ν的矩阵代数充当δ归一的直接检验。有限求积均非区间认证；谱完备性的依据是(2)的极分解证明。

下一步优先构造两观察者的SO⁺(3,1)规约内积及态，核Haar和DS反射/测度；再把(13)、(16)接到有明确物理态的两点谱。撤去物质截断需独立的正范数/域估计。现在不以ν²型群密度覆盖dν，也不以运动学μ唯一便声称旧n唯一。保留有限G及完整全息开放目标。

最终150项诊断通过，另完整保留[初始146项及两项失败](../results/stage4/ds3-mass-spectrum/initial-diagnostics.json)。失败来自紧质量包h∝(ν−.6)^4(1.8−ν)^4在T=20的有限时盒：范数亏损3.67314e−6、反极重叠误差2.76710e−7，均超过原2e−7目标。T=40、80范数亏损降至5.46933e−9和1.37811e−11；反极重叠趋向−.050551528061。真实正交补重叠趋零。没有放宽原目标后伪称T=20达标；[补充脚本](../code/complete_ds3_mass_diagnostics.py)另给尾估计并保留失败历史，不重跑其余计算。

尾部处理：令b±=h(1,−q)±/√(1+q²)，零延拓的前3阶边界导数消失，Plancherel给∫_{|t|>T}|χ|²≤T⁻⁸E₄，E₄=∫(|b_-''''|²+|b_+''''|²)dν。反极与正交补波包有相同尾范数界，由Cauchy–Schwarz也控制交叉尾。E₄的64/96点求积约3.26727403e6；T20估计1.27628e−4。该不等式解析成立，数值E₄未作区间包络，最终尾检查仍为经验验证，不称认证。

连续谱及束缚权重的四组Parseval测试分别使用(p,l)=(.5,0)、(.5,1)、(1.3,2)、(2.2,3)。后者连续权重仅约.00227920，束缚权重约1.26198698；漏束缚态的范数错误巨大，测试针对完整性主张而非仅重复定义。代码环境Python3.14.5、mpmath1.3.0、NumPy2.4.6、SymPy1.14.0；波包用双精度，其余50位。
