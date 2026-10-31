# 阶段三：量子观察者的质量、钟和制备字典

日期：2026-10-29。主来源固定为 Kolchmeyer–Liu，*Chaos and the Emergence of the Cosmological Horizon*，[arXiv:2411.08090v1](https://arxiv.org/html/2411.08090v1)，已逐节读取 §4.1—5.2、§6 的相关原文。本稿以原文等式作输入，推导本项目可实际实现的波包映射和两项有限制备任务。不修改此前 dS₃ 反冲基线。

**判断：存在可用的量子钟模型及可算顶角。** 原文已经把观察者总质量量子化，并给出物理内积和场导致的质量跃迁。此前 ζ(p) 可以按明确 Jacobian 等距嵌入其真空质量包扇区；这首先是态空间和自由演化的映射，并未识别两套完整可观测量或证明 dS₄、dS₃ 与 dS₂ 模型等价。若要让每次作用后仍有 p≥0，需要明确加入质量下界投影。单右观察者真空包二点对外质量确实对角；相位敏感的现成任务是双观察者传播核。

## 1 原文对象及可直接查验的字典

设 L=1，故质量与 proper time 都采用宇宙半径单位。物质自由场的 principal-series 参数在本稿记为 ν>0，以免与观察者谱变量 s 混淆；物质质量平方为 1/4+ν²。模型取严格 G→0，允许观察者质量有限并量子化其轨迹，未包含有限 G 度规反作用。

| 对象 | 原文定位及公式 | 本项目必须保留的区别 |
|---|---|---|
| 动态质量 | (4.2)—(4.4)：S=∫[P Q̇−Q√(−gττ)]dτ；一阶约束 −Pt²+Pθ²/cosh²t+Q²=0 | Q 是总静质量；P 是其共轭变量，不是我们的正钟能量 p |
| 单观察者预 Hilbert 空间 | (4.7)—(4.8)：\(\mathcal H_{obs}=\int_{s>0}^{\oplus}ds\,\mathcal P_s\)，\(\langle s'n'|sn\rangle=\delta(s-s')\delta_{n'n}\) | \(\mathcal P_s\) 是 principal-series 单粒子表示；这是约束前构件 |
| 精确质量 Hamiltonian | (3.9)、(5.2)：\(m(s)=\sqrt{s^2+1/4}\)，\(H_R=m(\hat s_R)\) | s 只在大质量极限近似等于 m，不能直接用 s 替代传播质量 |
| 质量下界 | §4.1：Q>1/2，因此 s>0 | 这是模型为限制 principal series 的选择，不是一般量子引力定理 |
| 真空物理扇区 | (4.29)—(4.31)：\(\mathcal H_0=L^2(\mathbb R_+,\rho_{Pl}(s)ds)\)，\(\rho_{Pl}=s\tanh(\pi s)/(4\pi^2)\) | 同一态中左右质量相等；不是两个独立钟波函数的张量积 |
| 含物质扇区 | (4.36)、(4.42)、(4.50)：\(\mathcal H=\oplus_N\mathcal H_N\)；N≥1 的 CG 通道带独立 sL、sR | 物质和两观察者由约束纠缠，顶角改变质量及轨迹表示标签 |
| proper-time 作用 | (5.6)：\(\phi_R(\tau)=e^{-iH_R\tau}\phi_Re^{iH_R\tau}\) | 符号按原文；τ是观察者固有时，不是全局坐标 t |
| 一般局域场插入 | (5.4) 的时空三模重叠 | 不存在简单 \(\varphi(\hat t_R,\hat\theta_R)\) 位置算符替代；也不是固定质量上的乘法 |

经典作用变分给出 Q̇=0、Ṗ=−√(−gττ)，所以 −P 随 proper time 增长。质量有下界后，不能由这项经典正则关系直接宣布存在全域自伴正则时间算符；本稿仍把读时 POVM 与质量谱分开。

原文主体只规范掉 SO⁺(2,1)，使用复 Hilbert 空间。若进一步规范时间反演，(4.51) 把波函数变为其复共轭，不变态要求相应系数为实。下文的任意复 ζ 嵌入默认主体版本；后面的实 ± 双峰比较同时符合这一实态条件。

## 2 物理内积：不能用普通投影抹掉群体积

两观察者真空态的形式表达是

\[
|\Psi_f\rangle=V^{-1/2}\int_0^\infty ds\,f(s)
\sum_{n\in\mathbb Z}(-1)^n|sn\rangle_L|s,-n\rangle_R\otimes|0\rangle_m.
\tag{D1}
\]

其中 V 是非紧 SO⁺(2,1) 的形式群体积。该式不属于预 Hilbert 张量积里的普通正规化 singlet 子空间。原文 (4.25)—(4.30) 在位置表示中固定左点及残余 boost，积分覆盖六个轨道区间，再得到

\[
\langle\Psi_{f_2}|\Psi_{f_1}\rangle_{phys}
=\int_0^\infty ds\,\rho_{Pl}(s)f_2(s)^*f_1(s).
\tag{D2}
\]

图计算中 (4.33) 的 \(\delta(0)\sum_n/V\to\rho_{Pl}(s)\) 是上述 Faddeev–Popov 归一化的简写，不能把 δ(0) 和 V 分别当有限数。加入物质后必须使用 (4.42) 的物理内积和连通 CG 图；任意断开 singlet 张量不一定有有限物理范数。

## 3 质量跃迁的真实顶角

把原文 (5.4) 写成

\[
\langle a,n_a|\phi_R|b,n_b\rangle
=\int_{dS_2}dvol_X\,
\overline{\Psi_{n_a}^{a}(X)}\,\varphi(X)\,\Psi_{n_b}^{b}(X).
\tag{D3}
\]

\(\Psi_n^s=\sqrt{s\tanh(\pi s)/\pi}\,\psi_n^s\) 是 (4.9) 的时空正频归一化模；该因子与物理内积的 ρPl 不是同一个量。创建一个物质量子的顶角由 (D3) 中 \(\varphi\) 换成 \(\overline{\psi_n^\nu}\) 得到，角动量要求 \(n_b=n_a+n\)。其时间积分是 CG 顶角，允许连续 a、b；没有强加普通平直空间式的质量相加条件。

质量变化的演化因子精确为

\[
\langle a|\phi_R(\tau)|b\rangle
=e^{-i[m(a)-m(b)]\tau}\langle a|\phi_R|b\rangle.
\tag{D4}
\]

因此各次插入须保留中间质量段。原文 §5.2 将 :φⁿ: 和导数按各自时空重叠另行穿衣；“先做复合场再穿衣”不是自动等于穿衣场的普通幂。

## 4 ζ(p) 的一个精确等距嵌入及其边界

这节是我们的变量变换与建模选择。选定一个基准静质量 m₀≥1/2，令

\[
s_0=\sqrt{m_0^2-1/4},\quad
p(s)=2\pi[m(s)-m_0],\quad p'(s)=\frac{2\pi s}{m(s)}.
\tag{D5}
\]

对原来归一化的 \(\zeta\in L^2(\mathbb R_+,dp)\)，定义

\[
\boxed{f_\zeta(s)=\mathbf1_{s\ge s_0}
\sqrt{\frac{p'(s)}{\rho_{Pl}(s)}}\,\zeta(p(s)).}
\tag{D6}
\]

直接换元给出

\[
\|\Psi_{f_\zeta}\|^2=\int_0^\infty|\zeta(p)|^2dp=1,
\qquad H_R=m_0+\frac{p}{2\pi}\quad\text{在 }\mathcal H_0\text{ 的该像上}.
\tag{D7}
\]

恢复单位时 m=LMobs、τ=τphys/L，原来的模时间为 tmod=τ/(2π)。于是 \(e^{-iH_R\tau}=e^{-i2\pi m_0t_{mod}}e^{-ipt_{mod}}\)，自由相位与此前钟能量 p 的演化精确相容。形式经典钟坐标为 −P/(2π)。式 (D6) 在 s=0 的点值不重要，按加权 L² 等价类理解。

但 (D6) 本身**不**提供所有算符之间的交织映射。尤其 m₀>1/2 时，\(\phi_R\) 可以把初态送到 a<s₀，此时按 (D5) 解释的内能为负。若要明确采用每次插入均保持非负内能的模型，应指定

\[
Q_{m_0}=\mathbf1_{H_R\ge m_0}\mathbf1_{H_L\ge m_0},
\qquad \Phi_R=Q_{m_0}\phi_RQ_{m_0}.
\tag{D8}
\]

初态只支撑于 b≥s₀ 与在每个插入间保留 (D8) 不同。原文 (5.21)—(5.24) 已研究质量窗压缩；(D8) 是这里选择的半线版本。在右侧二点上，它把中间 a 积分下限从 0 改为 s₀。本文以已知有限矩/二次型作计算，不借阶段二的另一套 Fock 证明来断言本模型所有压缩场的本质自伴性。

这还没有将读时 POVM、有限资源优化或 Han–Hu 尾核与有限质量轨迹顶角全面等同。原 dS₄ 自由场谱与本节 dS₂ 场质量参数也须分别记录。

## 5 可计算任务一：单右观察者的 Gaussian 涂抹二点

定义八个 Gamma 函数的乘积

\[
\mathcal G_\nu(a,b)=
\prod_{\epsilon_\nu,\epsilon_a,\epsilon_b=\pm1}
\Gamma\!\left(\frac14+\frac i2
[\epsilon_\nu\nu+\epsilon_a a+\epsilon_b b]\right).
\tag{D9}
\]

对于实 ν、a、b，它也等于四个模平方的乘积，严格为正。原文 (6.17)、(6.25) 的未归一化 HH 二点谱是

\[
D_\nu(a,b)=\frac{a\tanh(\pi a)b\tanh(\pi b)}{64\pi^7}
\mathcal G_\nu(a,b).
\tag{D10}
\]

HH 包为 fHH(s)=sech(πs)，其物理范数并非 1：

\[
Z_{HH}=\int_0^\infty\rho_{Pl}(s)\operatorname{sech}^2(\pi s)ds
=\frac1{8\pi^4}.
\tag{D11}
\]

由 tanh·sech²=−(sech²)'/(2π) 分部积分即可得到。60 位 mpmath 独立求积给出 0.00128324778183554189864409790838983676772847523170478284，与 (D11) 的差约 1.52×10⁻⁶⁴。这是原文归一化的核查，不给 D 再除一次 ZHH 后重复归一化。

原文 (5.8) 在任意外态间明确含 \(f_2(b)^*f_1(b)\)，所以由 (D10) 可推出

\[
\langle\Psi_{f_2}|P_{R,b}\phi_RP_{R,a}\phi_R|\Psi_{f_1}\rangle
=f_2(b)^*f_1(b)\cosh^2(\pi b)D_\nu(a,b).
\tag{D12}
\]

对角外质量并非事先猜测：右场与左 Hamiltonian 对易，而真空态中两侧质量相等。此处 P_R,s 是关于 ds 的谱密度记号，不是可归一化单质量本征态的普通投影；式 (D12) 积分时不再额外乘两份 Plancherel 测度。该任务不能读出 f 的外质量相位。更一般地，若所考虑的无界右侧乘积定义域由 exp[iθ(HL)] 保持，且各因子与它强对易，则任意纯右有限阶矩在 f→exp[iθ(m(s))]f 下不变。这项结论来自左质量谱分解，并不只限于二点；后续右四点任务应研究幅度、质量段和投影，不能期待仅换外包相位就出现响应。

令 w(b)=ρPl(b)|f(b)|²，为归一化质量密度。定义便于数值实现的条件**谱核**

\[
\boxed{V_\nu(a|b)=\frac{\cosh^2(\pi b)D_\nu(a,b)}{\rho_{Pl}(b)}
=\frac{a\tanh(\pi a)\cosh^2(\pi b)}{16\pi^5}
\mathcal G_\nu(a,b).}
\tag{D13}
\]

V 不是条件概率密度，其 a 积分一般不等于 1。取 proper-time Gaussian \(g_\sigma(u)=e^{-u^2/(2\sigma^2)}/(\sqrt{2\pi}\sigma)\)，对每个 φ 单独涂抹。以 (D6)、(D8) 的归一化态和压缩选择，精确二点二次型为

\[
\boxed{C_{\zeta,\sigma}(\tau)=
\int_0^\infty dp\,|\zeta(p)|^2
\int_{s_0}^\infty da\,V_\nu(a|b(p))
 e^{-\sigma^2[m(a)-m(b)]^2}
 e^{i[m(a)-m(b)]\tau},}
\tag{D14}
\]

其中 b(p)=√[(m₀+p/(2π))²−1/4]。正号相位来自原文反向 Heisenberg 约定，直接把 (6.15) 的 τ₂=τ、τ₁=−τ 即可核验。σ、τ 是 L 单位的 proper time，不能误用此前模窗口的 σ。若省略 (D8)，只将 a 下限改为 0，另行标注模型。

输出同时记录 C 本身与 C(τ)/C(0)，后者不是事先把 V 当概率。可直接检查 C(0)>0、C(−τ)=C(τ)*、|C(τ)|≤C(0)。相同 |ζ|² 的不同相位给相同右二点，这是所选任务的盲区，不是所有观测量都无相干效应。

## 6 可计算任务二：双观察者传播与实符号相干

原文 (4.34)—(4.35) 已给出从右观察者发射、左观察者吸收一个物质量子的规范不变传播算符 \(\widehat G_\nu\)。其对真空包的矩为

\[
\langle\Psi_{f_2}|\widehat G_\nu|\Psi_{f_1}\rangle
=\iint_0^\infty da\,db\,f_2(a)^*K_\nu^{LR}(a,b)f_1(b),
\]

\[
\boxed{K_\nu^{LR}(a,b)
=\frac{a\sinh(\pi a)b\sinh(\pi b)}{64\pi^7}\mathcal G_\nu(a,b)
=\cosh(\pi a)\cosh(\pi b)D_\nu(a,b).}
\tag{D15}
\]

它不含 δ(a−b)。改用 ζ 的平坦 dp 内积后，核是

\[
\boxed{\mathcal K_\nu(p',p)
=\frac{K_\nu^{LR}(b(p'),b(p))}
{\sqrt{\rho_{Pl}(b(p'))p'(b(p'))\rho_{Pl}(b(p))p'(b(p))}}.}
\tag{D16}
\]

这里分母的 p′(s) 表示 (D5) 的导数，不是外变量 p′ 的新函数定义。实现中应命名 `dp_ds` 避免歧义。这个核改变两端质量，属于双观察者任务；不能把它称为单右观察者关联函数。它实对称、逐点为正不意味着算符正定，本轮不作该宣称。

为保留相同能量分布而只改变可见的相干，选择两条互不重叠的实非负包

\[
u_c(p)=\sqrt{\frac4{3w}}\cos^2\!\frac{\pi(p-c)}{2w}\,
\mathbf1_{|p-c|\le w},\quad c>w>0,
\qquad \zeta_\pm=\frac{u_{c_1}\pm u_{c_2}}{\sqrt2},
\quad c_2-c_1>2w.
\tag{D17}
\]

二者都正规化、为实，且逐点 |ζ₊|²=|ζ₋|²。它们的质量均值相同；作为原平坦 dp 钟波函数，还具有相同的 \(\int|\zeta'|^2dp=\pi^2/(3w^2)\) 与零时间均值，所以同一个协变时间 POVM 约定下读时 RMS 也相同。此处不把该 POVM 自动认作所有物理观察者测量。

记 \(L_{ij}=\iint u_{c_i}(p')\mathcal K_\nu(p',p)u_{c_j}(p)dp'dp\)。则

\[
\langle G\rangle_\pm=\tfrac12(L_{11}+L_{22}\pm2L_{12}),
\qquad \langle G\rangle_+-\langle G\rangle_-=2L_{12}>0.
\tag{D18}
\]

严格不等式仅用两包非负、非零及 (D16) 的逐点正性，**不需要也没有证明核正定**。这给出一个明确的同概率分布、不同实相干符号的物理矩差别。它不等于发现新的对偶；原传播核来自论文，波包任务是其一个可复算应用。

## 7 实现规格和验证边界

Γ8 应在 log 域计算：

\[
\log\mathcal G_\nu(a,b)=2\sum_{\epsilon_a,\epsilon_b=\pm1}
\Re\log\Gamma\!\left(\frac14+\frac i2[\nu+\epsilon_a a+\epsilon_b b]\right).
\tag{D19}
\]

同时采用稳定 logcosh/logsinh，避免大质量指数先溢出再相消。D、V、LR 的解析比值及对称性分别验证，不互相代替独立参考。核归一化已经由另一路原文核查确认；正式数值还需完整保存参数、依赖、代码哈希、求积阶次和频谱截断。

建议第一组参数：m₀={1/2,2,5}，ν={1/2,1}，c₁=3、c₂=9、w=1，σ={1/4,1/2}，τ={0,1/2,1,2}。质量包紧支撑，双观察者核只在有限矩形求积。右二点的 a 尾应单独用截断扩展和独立自适应求积检查；不要把两阶 Gauss 结果相近叫严格误差认证。对固定紧支撑 b、ν，Stirling 公式给 V(a|b)=O(a⁻¹e⁻²πa) 于 a→∞，Gaussian 再抑制质量差；没有声称这已给显式全区间尾常数。

优先验证：HH 范数、(D6) Jacobian 归一化、八 Gamma 直接高精度乘积与 (D19)、D↔V↔LR 常数、左右核对称性、实 ± 包同 |ζ|² 和同导数范数、(D18)、右二点复共轭与正谱上界。有限样本 Gram 矩阵正性检查若做，也不能据此宣称完整 LR 算符正定。

原文 §6.3 的四点含 6j 顶角与非循环交叉项；HH 二点谱对称不意味着其全部代数上是迹态。大型质量微正则极限的次序见 (5.22)—(5.24)、(6.18)—(6.22)，不能用有限质量波包计算越过这些条件。

## 8 本轮推进了什么

已有精确行动入口：在原文物理内积上构造正规化质量包，按明确是否压缩中间质量的选择计算右侧涂抹二点；在同一模型中用原文双观察者传播核检验外质量相干。尚未完成的是这些观测量与此前 dS₄ Han–Hu 全部代数的算符级识别、有限 G、dS₃ 基线的跨维度匹配，以及有限制备四点动力学。Narovlansky 微观标签 s₀ 未在此被等同为观察者质量。
正式计算补记：上述两项波包任务现已执行并保存到 [阶段报告](../notes/stage3-packets.md) 和 [验证输出](../results/stage3/packets/validation.json)，145 项检查通过。独立 [核审查](stage3-observer-kernel-audit.md) 还证明了 Gaussian 双涂抹使右二阶矩对外质量一致有界，因此任意归一化真空质量包均可用于该二点二次型；本轮生产只使用明确紧支撑包。这个一般可积性结论不等于完整物理 Hilbert 空间上的场算符有界性。