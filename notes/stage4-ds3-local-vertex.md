# dS₃ 的局部场三点矩阵元与全局 CG 系数

2026-10-31。Alan Day / Chengheng。接续[全局Racah检验](stage4-global-racah.md)。[原文定位](../literature/stage4-ds3-local-vertex-sources.md)、[实现](../code/ds3_local_vertex.py)、[诊断](../results/stage4/ds3-local-vertex/diagnostics.json)、[与C₀比较](../results/stage4/ds3-local-vertex/clp-comparison.json)。

本轮在**明确选择的G=0局部乘法插入模型**中，恢复了包含完整三质量依赖和相位的标量CG系数。它比只匹配三点协变性多了一步：系数由KG归一的BD模式和实际时空积分确定。该构造尚未确定物理观察者质量直积分、规约Hilbert空间或CLP有限κ插入，不能据此宣布两点n已唯一。

## 1. 态、模式与作用域

半径ℓ=1，取global dS₃度规
\[
ds^2=-dt^2+\cosh^2t\,d\Omega_2^2,\qquad m_{\rm eff}^2=1+\nu^2,\quad\nu>0.
\]
标准球谐归一为∫|Y_lm|²dΩ=1，Y₀₀=1/√(4π)。这里l是SO(3)的K-type角动量，**不是**全局SL(2,C)主系列自旋J；单个标量J=0表示包含所有l。采用BD真空、[a_lm,a†_l′m′]=δ，φ=Σ(a u+a†u*)。

设uνlm=fνl(t)Ylm，χνl=cosh(t)fνl，则
\[
\chi''+\big[\nu^2+l(l+1)\operatorname{sech}^2t\big]\chi=0,\qquad
i(\chi^*\chi'-\chi'^*\chi)=1. \tag{1}
\]
这是LMW的Pöschl–Teller重写；其无反射in/out基不等于BD基。由下半欧氏球正则性及KG归一，
\[
\chi_{\nu0}=\frac{\sinh(\pi\nu/2-i\nu t)}
 {\sqrt{\nu\sinh\pi\nu}}
=N_\nu(e^{-i\nu t}-q_\nu e^{i\nu t}),\quad
q_\nu=e^{-\pi\nu},\quad
N_\nu=[2\nu(1-q_\nu^2)]^{-1/2}. \tag{2}
\]
在t=−iπ/2分子消零，除以cosh后为正则解。全局in/out纯指数不具有同一正则性。

令D_j=∂t−j tanh t，d_l=∏_{j=1}^l√(ν²+j²)。取固定相位
\[
f_{\nu l}(t)=\frac{i^l}{d_l\cosh t}D_l\cdots D_1\chi_{\nu0}(t). \tag{3}
\]
D_j将l−1方程映到l方程，Wronskian乘以ν²+j²，故(3)保KG归一；欧氏正则支也保持。i^l使其与MM(16)–(17)约定相同：
\[
f_{\nu l}=
\sqrt{\Gamma(l+1-i\nu)\Gamma(l+1+i\nu)/2}\,
(\cosh t)^{-1/2}
P_{-1/2+i\nu}^{-l-1/2}(i\sinh t). \tag{4}
\]
分支采用type-2 Legendre P。恢复ℓ时每个KG归一模式有ℓ^{-1/2}，无量纲ν对应m_eff²ℓ²=1+ν²；以下所有数值都在ℓ=1。

## 2. 实际插入，而非仅由协变性猜系数

仿照KL(5.4)的结构，**定义本轮dS₃模型**
\[
\langle\nu_o,l_om_o|\mathcal O_\phi|\nu_i,l_im_i\rangle
:=\int_{\mathrm{dS}_3}dV\,u_{\nu_o,l_om_o}^*(X)
 \phi(X)u_{\nu_i,l_im_i}(X). \tag{5}
\]
它当前是固定质量模式上的形式矩阵元；不先声称它在某质量直积分上闭或自伴。局部耦合取1是模型输入，不是由对称性导出。

场的湮灭和创建系数分别为
\[
T_{\rm ann}=\int dV\,u_o^*u_i u_p,\qquad
T_{\rm em}=\int dV\,u_o^*u_i u_p^*. \tag{6}
\]
下文用η_j=+1表示u_j、η_j=−1表示u_j*，以及Tη=∫dV∏u_j^{(η_j)}统一记号。有限球谐组合的每一模式在两端O(e^{−|t|})，故体积乘三模式为O(e^{−|t|})，绝对收敛。不是无限Fock算符域的证明。三个KG模式和体积给ℓ^{3/2}量纲；定义物理O还需所选局部耦合及观察者谱约定。

Killing导数作用于三模式积后是散度，角向无边界，时间通量为O(e^{−|t|})，所以此三线性泛函不变。有限群变换的光滑共形因子在球面上有上下界，仍有相同端点衰减；可将此论证用于有限模式的群轨道。对一般主系列参数，DS(17)–(19)的三点协变核唯一到一个整体三标签系数。该标准表示论事实不提供物理系数，我们下面直接算它。

## 3. 所有有限模式的Beta求和与真实尾界

对σ=±1定义多项式
\[
P_{0,\sigma}=1,\quad
P_{l,\sigma}=(1-y^2)P'_{l-1,\sigma}
 +(i\sigma\nu-ly)P_{l-1,\sigma},\qquad y=\tanh t. \tag{7}
\]
(3)为两项P e^{iσνt}乘sech t，含系数i^lN/d_l及−q。三项乘积只有八个频率分支，角积分为通常Gaunt系数（一般m时保留共轭与−m相位）。剩余每个单项用
\[
I_n(\omega)=\int_{\mathbb R}dt\,\operatorname{sech}t\,\tanh^nt\,e^{i\omega t}
=\sum_{k=0}^n {n\choose k}2^k(-1)^{n-k}
B\!\left(k+\frac{1+i\omega}{2},\frac{1-i\omega}{2}\right). \tag{8}
\]
证明仅用y=tanh t、x=(1+y)/2。n=0给π/cosh(πω/2)。

若第j个模式展开为sech t Σσ e^{iωσt}p_{jσ}(tanh t)，令b_j为全部p系数绝对值之和。则截断|t|≤T的绝对余项
\[
|E_T|\le4e^{-T}|\mathcal G_{l_1m_1,l_2m_2,l_3m_3}|
\prod_j b_j. \tag{9}
\]
这只认证解析尾界，不认证盒内求积。程序保存Beta项的相消比，不能把小相对残差不加解释地当所有相位抵消的总误差。

## 4. 球对称矩阵元闭式

记S=ν₁+ν₂+ν₃，
\[
\mathcal D=\cosh\frac{\pi S}{2}\prod_{j=1}^3
\cosh\frac{\pi(S-2\nu_j)}2,\qquad
B_*=\frac{\pi}{4\sqrt{4\pi}\mathcal D}
\sqrt{\prod_{j=1}^3\frac{\sinh\pi\nu_j}{\nu_j}}. \tag{10}
\]
用(2)八项展开和I₀可直接代数约分：

| 共轭模式数 | Tη,000 |
|---|---|
| 0或3 | −B_* |
| 1，唯一共轭腿为r | B_* cosh(πν_r) |
| 2，唯一未共轭腿为r | B_* cosh(πν_r) |

于是ann的r为out，em的r为in。后者与交换in/out后的ann共轭一致，符合Hermitian场的矩阵元要求。cosh因子来自BD两频率混合；不是可凭直觉删除的热因子。

## 5. 把KG模式接到边界表示，并固定shadow相位

定义
\[
c_l(\nu)=i^l\prod_{j=1}^l\frac{-i\nu-j}{\sqrt{\nu^2+j^2}}
=(-i)^l\frac{(1+i\nu)_l}{\sqrt{\prod_{j=1}^l(\nu^2+j^2)}} ,
\qquad |c_l|=1. \tag{11}
\]
对有限模式组合a_lm，归一未来数据为Aν=Σa_lm c_l(ν)Y_lm，因而此映射在KG模系数空间与L²(S²)间保范数。晚时展开为
\[
u(t,n)=2N_\nu e^{-t}
\left[e^{-i\nu t}A_\nu(n)
-q_\nu e^{i\nu t}\mathcal S_\nu A_\nu(n)\right]+O(e^{-3t}), \tag{12}
\]
\[
\mathcal S_\nu Y_{lm}=
\frac{(1-i\nu)_l}{(1+i\nu)_l}Y_{lm}. \tag{13}
\]
这是单位模shadow相位，s₀=1，S_{−ν}Sν=1。其球面核为
iν/(π4^{iν}) |n−n′|^{-2+2iν}，须按解析延拓的分布理解；在Re指数=−1不能误作绝对收敛的普通积分。两个晚时分支在BD空间中已相关，不是两个独立Hilbert空间。

给每条有共轭选择的边取w_j=1+iη_jν_j。三线性积分所需核的权重为2−w_j，故
\[
a_i=-\frac12+\frac i2(\eta_j\nu_j+\eta_k\nu_k-\eta_i\nu_i),
\quad
\mathscr K_\eta(n_1,n_2,n_3)=
\prod_{\rm cyc}|n_j-n_k|^{2a_i}. \tag{14}
\]
在球面测度下定义Bη(A₁,A₂,A₃)=∫∏dΩ_j A_j^{(η_j)}Kη。不要把核权重2−w与波函数权重w混用。

立体投影下Ω(z)=2/(1+|z|²)，dΩ₂=Ω²d²z，
|n_i−n_j|²=Ω_iΩ_j|z_i−z_j|²，A_plane=Ω^w A_sphere。每个积分点的Ω幂为2−w_i+a_j+a_k=0，因此平面核恰为∏|z_j−z_k|^{2a_i}，无遗漏常数。

ann重排为输入(in,probe)、输出out，DS标签a_DS=−iν：此平面核正是DS(17)的W；em是输入(out,probe)、输出in之W的共轭。对一般参数，Tη=Cη Bη，三标签系数由下一节的球面周期固定。这里的匹配是局部标量CG，不是上一轮被否定的H→完整6j识别。

## 6. 球面周期及整体系数

令A=a₁+a₂+a₃。标准球面三角积分的**归一平均**
\[
J(a_1,a_2,a_3)=\frac1{(4\pi)^3}\int(d\Omega)^3
\prod_{\rm cyc}|n_j-n_k|^{2a_i}
=4^A\frac{\Gamma(A+2)\prod_i\Gamma(a_i+1)}
 {\prod_i\Gamma(A-a_i+2)}. \tag{15}
\]
绝对收敛域含Re a_i>−1、Re A>−2；本轮Rea_i=−1/2、ReA=−3/2位于其中。

为避免原文整体常数的不一致，独立核归一：(2−2x)^a球面卷积的第l本征值为
\[
\lambda_l(a)=\frac{4\pi\,4^a}{a+1}\frac{(-a)_l}{(a+2)_l}.
\]
三卷积迹除以(4π)³得到Σ_l(2l+1)∏(-a_i)_l/(a_i+2)_l乘4^A/∏(a_i+1)。DLMF16.4.9的Rogers–Dougall ₅F₄公式立刻给(15)；先在平滑核域论证迹，再解析延拓到上述绝对收敛域。级数项为O(l^{−2ReA−5})，在本轮是O(l^{-2})。J(0,0,0)=1，单边J(a,0,0)=4^a/(a+1)也直接积分可得。此公式是已有三角积分的规范化复推，不是新数学恒等式。

因为所有归一l=0数据为Y₀₀，Bη,000=(4π)^{3/2}J(a)。故完整复系数为
\[
\boxed{C_\eta=\frac{T_{\eta,000}}{(4\pi)^{3/2}J(a)}}. \tag{16}
\]
Gamma反射给出简洁的模方：
\[
|C_\eta|^2=\frac{F_\eta}{64\pi^3\mathcal D},\qquad
F_\eta=
\begin{cases}1,&0或3条共轭腿,\\
\cosh^2(\pi\nu_r),&1或2条共轭腿，r如表.\end{cases}\tag{17}
\]
更高角动量不是以l=0值直接冒充：对(l,l,0)三腿，旋转平均P_l(n₁·n₂)将边界矩阵元化为(15)中a₃的有限整数移位，并乘c_l(ν₁)^{(η₁)}c_l(ν₂)^{(η₂)}。本轮l=1,2,3的独立有限模式积分通过了这一匹配；另直接检查非平凡boost Ward关系。

## 7. 与此前直接C₀配对的关系

冻结的[直接C₀](stage4-recoil-poles.md)质量依赖是
\[
D_{\rm CLP}=\frac{\pi e^{\pi S}}{\mathcal D}
\prod_j\frac{\sinh\pi\nu_j}{\nu_j}.
\]
由(17)得
\[
\frac{|C_\eta|^2}{D_{\rm CLP}}
=\frac{F_\eta}{64\pi^4}
\prod_j\frac{\nu_j e^{-\pi\nu_j}}{\sinh\pi\nu_j}. \tag{18}
\]
这一次差别确实是**明确可分离的腿因子**，并有由ann/em方向决定的cosh²。因此G=0直接三点配对存在一个局部规范接口；它不修复上一轮H/Z的四点比值差异。不能把(18)直接设为旧两点n：旧谱还含观察者z、质量积分协议、三点基的逆范数、态与物理缝合。

## 8. 验证、来源异常与继续位置

[249项诊断](../results/stage4/ds3-local-vertex/diagnostics.json)与[24项C₀比较](../results/stage4/ds3-local-vertex/clp-comparison.json)通过。前者：KG/ODE/独立Legendre各36，l=0闭式及系数模方各24，高l边界匹配27，boost Ward27，直接时积分及加密8，Dougall级数2，零/单边周期4，shadow16，加精度9。ν三元组(.7,1.1,.6)、(1.4,.9,1.2)、(3.2,2.7,.5)，40/60位。最大记录scaled error约6.34e−38；四个直接积分使用T=96、步长3及2，最大相对差约1.55e−41；相应解析尾/信号≤2.46e−40。高l相消比最大632.41。盒内积分未区间认证；249不是独立物理预测数。

原文核对发现两个需要显式记录的地方：MM HTML(14)的角项印为cosh²，与度规(4)及模式(16)不符，本轮从度规重建为sech²；BEGMP的HTML(93),(127),(128)在d=3合用会给J(0)=16而非定义要求的1，对应旧PDF(74),(101),(102)也有该常数。这不是作者确认的勘误；本轮不用该常数，改用(15)的独立归一推导。旧原文未改。

下一步应区分KG固定质量归一与观察者的时空L²/质量δ归一：核dS₃全局径向连续谱的双重数、BD两分支、可能的离散束缚谱和物理质量下界，再选择并说明观察者投影。随后将(16)与DS正交及反射测度一起缝合，先完成自由场两点及独立探针断开四点基准。有限G与完整全息对偶仍未完成。
