# 缝合测度、真空标签与约化两点的分布极限

2026-10-31。Alan Day / Chengheng。配套[原文定位](../literature/stage4-recoil-identity-sources.md)、[代码](../code/recoil_identity.py)、[数值记录](../results/stage4/recoil-identity/diagnostics.json)。本轮为来源核查与自查数学推导，不是物理四点构造。

## 1. 实参数理论究竟缝合什么

以下 CEZ 指 arXiv:2304.13650v3；其 HTML 采用连续编号，不能套用早期 PDF 的章节编号。实 b>0、正实内部 P 的标准共形块规范中，(23)–(25)给

\[
\langle s|s'\rangle=\frac{\delta(s-s')}{\rho_s(s)},\quad
\rho_s(s)=\rho_0(s)C_0(P_1,P_2,s)C_0(P_3,P_4,s),\quad
\rho_0(P)=4\sqrt2\sinh(2\pi bP)\sinh(2\pi P/b).
\tag{1}
\]

因此四孔球的分解恒等算符是 \(\int_0^\infty ds\,\rho_s(s)|s\rangle\langle s|\)。若先以两个 C₀ 的平方根重标基矢，测度才变成单独的 ρ₀；这正是原文(27)–(28)，不是可随意省略顶角的理由。

取融合变换 \(|s\rangle=\int_0^\infty dt\,F_{st}|t\rangle\)，内积保持要求
\[
\int_0^\infty dt\,\frac{F_{st}\overline{F_{s't}}}{\rho_t(t)}
=\frac{\delta(s-s')}{\rho_s(s)}.\tag{2}
\]
这由(1)直接推出；与 CEZ(55)的交换关系及附录(122)融合逆关系相容，但我们未数值计算完整 F。逆变换的双线性恒等式与带复共轭的酉性恒等式不可混用。

若基矢变成 \(a_s|s\rangle\)，同时令 \(\rho'_s=\rho_s/|a_s|^2\)、\(F'_{st}=a_sF_{st}/a_t\)，则(2)保持。这只是基规范自由；并不证明任意改变物理插入 n 后都能满足固定 F 的 crossing。

CEZ(29)的恒等腿是分布：
\[
\lim_{P_3\to iQ/2}C_0(P_1,P_2,P_3)
=\rho_0(P_1)^{-1}\delta(P_1-P_2),\qquad Q=b+b^{-1}.\tag{3}
\]
P>0 约定吸收了反射冗余；在全实线上另有反射 delta，见附录(134)、(140)。原文(46)区分外部/入射恒等表示与输出连续谱；不能把输出 P 积分中某个普通点当作恒等态。恒等块本身非正规化。尤其不能将(3)直接取模平方，得到未定义的 δ²。

## 2. de Sitter 延拓尚缺的接口

CLP §6明确区分固定拓扑振幅的解析延拓与新的 dS 缝合理论。CEM arXiv:2501.01486v1 的(20)、(24)、§5则从晚时切片的模空间/Teichmüller积分定义另一内积；普通共形块不将其对角化。故(1)中的ρ₀以及观察者权重 z(m)都不能未经推导认作该内积的 Plancherel 测度。这里核查的是这些指定版本，不宣称搜索已排除所有后续文献。

在一般非正交基中，真正要插入的是 Gram 核的逆：
\[
\mathbf1=\iint da\,db\ |a\rangle (G^{-1})_{ab}\langle b|,
\qquad G_{ab}=\langle a|b\rangle.\tag{4}
\]
此式仅在测试空间上存在所需逆核时成立，零模需先处理。不能为了套用(2)而假定 G 对角。正性、共同域与收敛性仍是物理任务。

## 3. 共同角路径不包含双真空

直接使用 CLP(3.2)、(3.4)、(3.6)、(3.44)的标签，而不是从时域曲线猜 identity：
\[
b_\pm^2=\pm i/\kappa,\quad P_\pm=\pm\frac{i\alpha}{4\pi b_\pm},
\quad h_\pm=Q_\pm^2/4+P_\pm^2.
\]
沿共同角 \(\alpha=2\pi(1-\nu/\kappa)\) 精确得到
\[
h_+=\frac{1-i\nu}{2}+\frac{i(1+\nu^2)}{4\kappa},\qquad
h_-=\frac{1+i\nu}{2}-\frac{i(1+\nu^2)}{4\kappa},\qquad h_++h_-=1.\tag{5}
\]
这是标签代数恒等式，和原文消除 framing 的关系一致。ν=i 给(1,0)，ν=−i给(0,1)，任何有限κ都成立；ν=0更不是恒等表示。这里的h±是两个延拓因子的权重，**不是**可直接相加得到物理标量维数Δ的一对普通左右权重。

各因子分别取真空的正角分支要求
\[
\alpha^{\rm vac}_\pm=2\pi(1+b_\pm^2)=2\pi(1\pm i/\kappa),\tag{6}
\]
它们不在共同角路径上。负角反射分支不消除两者差别。真正的双真空需要独立两套标签及其配对/积分规则。共同实ν上的模方，在离开实轴后作共同ν的亚纯延拓，与对每个复ν重新取模方是不同操作。

所以此前 NEXT 中“Δ→0涉及ν=±i”只适用于场质量关系 \(m_{\rm field}^2=1+\nu^2\) 的端点提示，不足以识别 CLP 配对图的双恒等腿。这个区别改变了下一步：不能将约化谱端点通过与双真空缝合通过混为一项。

## 4. 约化谱的可证明分布极限

现在完全转到固定质量κ→∞之后的两点函数；半径1。只对已得闭式作亚纯延拓。令 \(\nu=i(1-\epsilon)\)、\(0<\epsilon<1\)，除以发散的零时值
\[
g_\nu(0)=\frac{1-\epsilon}{4\pi\sin\pi\epsilon}
\]
后得到正概率谱
\[
p_\epsilon(\omega)=\frac{\sin\pi\epsilon}
 {2(1-\epsilon)[\cosh\pi\omega-\cos\pi\epsilon]},\qquad
\chi_\epsilon(t)=\int_{\mathbb R}d\omega\,e^{-i\omega t}p_\epsilon(\omega)
=\frac{\sinh((1-\epsilon)t)}{(1-\epsilon)\sinh t}.\tag{7}
\]
在t=0取连续值1。傅里叶规范不含额外1/(2π)。可由上一轮闭式在0<Imν<1内解析延拓，或直接留数求和；谱指数尾在该开条带的紧子集上支配积分。

由二阶导数得
\[
v_\epsilon=\int\omega^2p_\epsilon(\omega)d\omega
=\frac{2\epsilon-\epsilon^2}{3},\qquad
\int_{|\omega|>D}p_\epsilon\le v_\epsilon/D^2.\tag{8}
\]
因此pε弱趋δ，且对于有界二阶导数的C²测试函数，偶性消去线性项，误差≤\(\|f''\|_\infty v_\epsilon/2\)。有界连续测试函数的弱收敛由集中性及总质量1推出。这不是未正规化无质量标量的有限二点函数，也不构造Hilbert空间中的恒等算符；归一化已经除去了红外发散。

### 极点与穿越留数

谱极点在 \(\omega=i(2n\pm\epsilon)\)，最近两点±iε夹住实轴。上点留数为\(1/[2\pi i(1-\epsilon)]\)，下点相反。从ε>0逼近0时未穿过实轴，但端点发生夹逼。若继续到ε=−δ，0<δ<1，则固定实轴积分给
\[
\int_{\mathbb R}e^{-i\omega t}p_{-\delta}(\omega)d\omega
=-\frac{\sinh((1-\delta)t)}{(1+\delta)\sinh t},\tag{9}
\]
它不是(7)的解析延拓。两颗跨越极点补项为
\[
\chi_{-\delta}(t)-\int_{\mathbb R}e^{-i\omega t}p_{-\delta}(\omega)d\omega
=\frac{2\cosh(\delta t)}{1+\delta}.\tag{10}
\]
在t=0也有有限非零差。该等式同时由(7)、(9)的双曲恒等式和两极点留数得到。这里只解决约化谱的两极点，未证明完整有限κ的双Gamma/double-sine轮廓可沿同一路径交换极限。

## 5. 有限质量截断与残余归一化

固定M>m*>0，D=2(M−m*)，A=cosh2πM。定义
\[
R_{M,\epsilon}(\omega)=\left(1-\frac{\omega^2}{4M^2}\right)
\frac{A-\cosh\pi\omega}{A-\cos\pi\epsilon},\quad
q_{M,\epsilon,\lambda}=\mathbf1_{[-D,D]}p_\epsilon R_{M,\epsilon}
\left(1+\frac{\lambda\omega^2}{4M^2}\right),\quad\lambda\ge0.\tag{11}
\]
这是声明过的谱模型；在此区间0≤R≤1。记 \(v_D=\int_{-D}^D\omega^2p_\epsilon d\omega\)、\(T_D=\int_{|\omega|>D}p_\epsilon d\omega\)。精确质量亏损为
\[
1-\int q_{M,\epsilon,0}
=T_D+\frac{v_D}{4M^2}
+\frac{\sin\pi\epsilon}{(1-\epsilon)(A-\cos\pi\epsilon)}
 \left(D-\frac{D^3}{12M^2}\right).\tag{12}
\]
证明：展开1−R；其第二项含coshπω−cosπε，恰与pε分母抵消，剩下二次多项式积分。用(8)给显式上界
\[
B_{M,\epsilon}=\frac{v_\epsilon}{D^2}+\frac{v_\epsilon}{4M^2}
+\frac{\sin\pi\epsilon}{(1-\epsilon)(A-\cos\pi\epsilon)}
 \left(D-\frac{D^3}{12M^2}\right)=O_{M,m_*}(\epsilon).\tag{13}
\]
于是对所有实t，
\[
\left|\widehat q_{M,\epsilon,\lambda}(t)-\chi_\epsilon(t)\right|
\le B_{M,\epsilon}+\frac{\lambda v_\epsilon}{4M^2}.\tag{14}
\]
同一右侧也控制两密度的L¹距离。固定λ、M、m*时，qε弱趋δ；乘回z(M)即z(M)δ。所有λ都过关，却在实主系列上保留不同的有限M首修正。因此**领先QFT匹配加这个正规化δ端点测试，仍不唯一确定n**。这没有声称这些nλ都能满足物理完整crossing。

对于固定的有界连续偶因子n（或在该有限区间上一致收敛到连续函数的nε族），端点测试只读出n(M,0)；若要把极限配成所选恒等谱规范，只能要求其对角值为1。离对角值需要别的输入。若nε在离对角处随ε发散，则不能引用这一结论。

## 6. 非一致时间极限

固定ε>0，χε(t)在|t|→∞趋0；固定t而ε→0则趋1。更具体地，固定τ≠0，
\[
\lim_{\epsilon\downarrow0}\chi_\epsilon(\tau/\epsilon)=e^{-|\tau|}.\tag{15}
\]
由(14)，同一结论适用于固定M、m*、λ的q。这是模型的红外时间尺度，不是混沌指数或有限G动力学。δ弱收敛不提供全时间趋1；(14)是相对于χε的一致绝对误差。

## 7. 数值与研究决策

106个诊断全部通过，分为实轴傅里叶/二阶矩16项、跨极点18项、有限质量亏损/界/加精度36项、有限κ标签代数36项（后者包含同一真空标签在不同ν循环下的重复核对，不能当作36个独立物理检验）。50位计算，有限质量样本65位复算，最大加精度差6.41×10⁻⁵²；不是区间误差认证。ε=.2,.05,.01,.002的全实轴谱；有限M=1,3，m*=.2，ε=.1,.02,.004，λ=0,2；未以高时间极限公式冒充高振荡求积。版本和SHA256保存在JSON。

下一步不继续调n来逼近想要的曲线。优先建立G=0全局SL(2,C)表示/partial-wave内积与现有反冲H的完整标签、反射和谱测度字典，或明确证明其并非同一核；这是未来复Virasoro Gram问题的可检验基准。若研究双真空，应保留独立α±，先在实b理论缝合，再逐个跟踪轮廓继续；不得将(3)的δ平方。有限κ物理四点、有限G及完整对偶均保持开放。
