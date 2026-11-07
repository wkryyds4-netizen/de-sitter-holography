# 有序 Wick 后代、静态回归与完整中性记忆核

2026-11-07，本地阶段。前驱是 size-bridge 与 neutral-marked-series。正式推导见 [LaTeX](../manuscript/susskind-memory-wick-descendants/main.tex) 和 [PDF](../output/pdf/de-sitter-susskind-memory-wick-descendants.pdf)。本轮前一goal turn为progress：已封存rank/size桥和离散分布；当前沿真实后代推进，没有重复issue或旧main。

## 已得到的结果

所有固定次数的有序Gaussian-Wick后代、它们的chaos投影和Gaussian OU像，都有明确的共同J两能量图规则。扣除确定性能量回归后，无field桥的图成组抵消；余图给外部q-Hermite系数的有限带宽和指数衰减。因此固定次数残量有L²(ρ_q⊗ρ_q)谱密度，exact pinching在annealed L²下消失。旧uniform series尾再解除次数截断，得到完整有限实时间neutral B的协方差收敛构造。全联合memory的charged A及AB/BA仍待计算。

模型保持标准CAR、集合可重叠的complex Gaussian p-body、跨电荷同一原J、β=μ=0归一annealed全迹，偶N，λ_N=p²/N→λ>0，固定J_s，ℏ=1/string时间。C=J_s²/λ，q=e^(−4λ)∈(0,1)。每个field次数和外部能量次数先固定，再取double scaling；没有联合小λ、增长时间或有限G结论。

## 1. 通用primitive与chaos计数

设Z_j=W_j(H^j)、U_j=(jp−R)Z_j，T_{v,w}=U_v diamond Z_w，F_{r,v,w}=H^rT_{v,w}。diamond是coupling Gaussian Wick乘法，系数矩阵按原顺序相乘。整体v+w source block禁止内部pair；物理size mark只作用其中前v个slot。dagger反转顺序，把prefix mark变为suffix。

整个field有L_F=r+v+w个raw slots；完整两场图π中，若i_F条弦的两端均在这个field内，则其Gaussian chaos degree是d_F=L_F−2i_F。W_d投影只保留d_F=d，q^Ncal只乘q^{d_F}。外部energy槽不计入field。对V_j和eta_j的显式C因子另行保留；程序C=1时记录physical degree j，恢复矩统一乘C^{(j+k+a+b)/2}。

两个局部mark的incident chord sets为I1/I2，K(I)=|I|²−4cr(π|I)，K12=K(I1∩I2)。0/1/2个mark分别给1、λK1、λK12+λ²K1K2，再乘C^{L/2}q^cr。负系数不删，同站点项保留。

rank到size的替换先在U_v上用已证L²桥，再用右Wick creation界

\[
\|X_v\diamond Z_w\|_2\le\sqrt{(v+w)!/v!}\,M^{w/2}\|X_v\|_2
\]

传递；H的固定degree乘法界和chaos/OU收缩再覆盖primitive。没有使用任意矩阵Gaussian场相乘的L²连续性。

## 2. 回归系数和切图

H_d^{(q,C)}记monic q-Hermite，P_d=[d]_q!，范数平方C^dP_d。纯source有

\[
\operatorname{Reg}T_{v,w}=\kappa_vH_{v+w}^{(q,C)}(H),
\qquad\kappa_v=\lambda\frac{(v^2-4q\partial_q)P_v}{P_v}.
\]

证明按全桥图的Mahonian置换分解：prefix v内部、suffix w内部和shuffle，prefix内部逆序权重独立于所选外部端点集合。suffix mark同理。

更一般的F_d=W_d[H^rT_{v,w}]先固定field内部partial matching γ，禁止source内部pair，剩d开放端点。内部弦之间及内部弦与开放leg的cross总数cγ只由内部配置决定；开放leg之间cross为某个排列σ的inv。Kγ(σ)也只依赖γ和σ，不依赖外部端点绝对位置或其属于哪个energy段。因此

\[
g_{\gamma,d}=C^{|\gamma|}q^{c_\gamma}\lambda
\frac{\sum_{\sigma\in S_d}q^{\operatorname{inv}\sigma}K_\gamma(\sigma)}{P_d(q)}.
\]

按γ及原系数求和得Reg F_d=g_F,d H_d。d=0包含在内。q^Ncal在同一γ组乘q^d。与单个外部H_d配对给完全相同的置换和，所以这确是能量回归；不是只凭无桥图定义的替代物。

## 3. 无桥抵消与全外部次数的谱界

固定两个field内部图、所有外部partner位置；若无弦直接连接两个field，则I1∩I2为空、K12=0。各自开放端点排列不改变与完全外部弦的cross；两个排列和独立。故无桥组恰为各自能量回归乘积，扣除回归后整组消失，单图不必为零。内部partial matchings、chaos projection及q^Ncal系数都不破坏此结论。

外部能量换归一q-Hermite p_n=H_n/√(C^nP_n)。标准多block q-Wick组合消去所有energy block内部pair；删除的弦不属mark incident sets，也不改剩弦间cross，因此marked权重可以作为该消去的旁系数。此处是图恒等式，不是把有限N coupling Wick等同spectral Wick。

两个field raw长度总和≤D。余图至少一条field桥；设r条energy–energy桥，端点数给|m−n|≤D、r≥min(m,n)−D。每条field桥与全部r条energy桥交叉，产生q^r。固定field端点选择最多const_D(m+n+1)^D，剩energy桥置换和为P_r，且

\[
P_r/\sqrt{P_mP_n}\le(q;q)_\infty^{-1}.
\]

故中心化残量的归一系数满足

\[
|M_{mn}|\le K_D(m+n+1)^D
q^{\max(0,\min(m,n)-D)}1_{|m-n|\le D}.
\]

所有常数依赖固定field/C/λ/q，不依赖m,n。平方可和给L²(ρ_q⊗ρ_q)密度；紧支撑上的多项式稠密性使系数确定原谱测度。该论证直接实现Bikram–Mukherjee的HS策略所需接口，无须猜测marked fields存在标准q-Fock实现。

量化尾为2(2D+1)K_D² Σ_{n≥L}(2n+D+1)^{2D}q^{2(n−D)}，L≥D。可用相邻比作几何上界，(q;q)_∞^−1≤exp[q/(1−q)²]。这是对极限系数的所有外部degree界，不是finite-N的degree一致误差。

## 4. 精确pinching与实际静态部分

对一组实际残量R_j,N，定义正matrix measure

\[
\nu_{jk,N}(A,B)=E\operatorname{tr}
[R_{j,N}^\dagger P_H(B)R_{k,N}P_H(A)].
\]

E是input/column，E′是output/row。正性和普通q-Gaussian高矩支配给弱收敛和每轴紧支撑±2√(C/(1−q))。HS密度说明极限对角线Δ质量为零。有限N exact pinching含整个简并块，||P_H R_N||₂²=ν_RR,N(Δ)，Portmanteau对闭集Δ给其趋零。这里是annealed L²，不是almost sure，也没有交换无限时间平均。

## 5. eta与V的全部确定性回归

固定degree极限中的递推为

\[
r_j=4\lambda C\sum_{v<j}q^v T_{v,j-1-v},\qquad
\eta_{j+1}=H\eta_j-C[j]_q\eta_{j-1}-r_j,
\quad\eta_0=\eta_1=\eta_2=0.
\]

设g_j=Reg eta_j，g_0=g_1=g_2=0，则

\[
g_{j+1}(x)=xg_j(x)-C[j]_qg_{j-1}(x)
-4\lambda C\left(\sum_{v<j}q^v\kappa_v\right)H_{j-1}^{(q,C)}(x).
\]

展开g_j=Σ_d g_jd H_d。V_j=4λq^j U_j+(q^j−q^Ncal)eta_j的回归为

\[
v_j(x)=4\lambda q^j\kappa_jH_j(x)
+\sum_{d\le j-2}(q^j-q^d)g_{j,d}H_d(x).
\]

R_j=V_j−v_j(H)即无静态质量的动态输入。低阶g3=−4λ²CqH1；g4=−8λ²Cq(1+3q)/(1+q)H2−4λ²C²q。最后scalar恢复旧eta4净chaos0，必须在所有chaos操作后才统一减去。R1=4λq(U1−λH)，其谱16λ³Cq² μ_q，乘记忆场prefactor的平方4C²后恢复旧force谱。

## 6. 完整finite-time neutral核

c_j(s)=∫e^{isE}H_j(E)dρ/(C^jP_j)。有限辅助Jacobi条目收敛且一致有界，故其指数在紧s窗强收敛，c_j,m(s)收敛；α_m→1由连续性处理。旧有限参数表达保留mε_m C_m/p→2C和两侧e^(−ihs/2)。

Reg是极限joint Gram空间到energy闭包的正交投影，故||R_j||₂²≤lim||V_j,N||₂²，旧K_S(j+1)²X_S^j/√j!尾同样支配残量。先固定J减去v_j(H)，利用无对角质量传递exact pinching，再由L²尾和pinching收缩去J，得到

\[
\lim_NE\operatorname{tr}[B_N(s)^\dagger e^{iTH}B_N(u)e^{-iTH}]
=4C^2\sum_{j,k}\overline{c_j(s)}c_k(u)
\int e^{i(s-u)(E+E')/2+iT(E'-E)}d\nu_{jk}(E,E').
\]

每个ν_jk由有限图给所有M_mn，再以L² Hermite级数重构密度；双time级数由CS绝对可和。这是全有限时间的收敛构造，不是仅列Taylor系数，也不宣称已经把核化为一个简单特殊函数。

场的L²尾给谱TV尾，AC对TV极限封闭，所以全时间残量仍AC；无限和密度仅需L¹，不另断言L²。次数截断尾对全部实外部T一致；N→∞的谱极限针对固定T，不能升级为整条时间轴上一致。全charge转换使用旧中央窗和binomial尾。

## 验证与继续

descendant-checks.json首次运行通过：42完整full-Fock halfsize与Gaussian-Hermite/坐标配对比较，最大误差6.661338147750939e−15；12项V1/2旧矩校准；10项source回归；6项非空无桥中心化精确相消，包含OU eta4/5；59个两场多项式中19行含负系数，共31负项；9个eta chaos范数。有限CAR用浮点归一矩阵、保留全部坐标；formal多项式用Fraction/SymPy精确算术。没有实际数值pinching或任意degree的数值扫查，general theorem由切图/HS证明承担。

独立source-regression-checks.json首次运行成功约0.46秒：d1–7共5913置换、56个prefix/suffix整数多项式恒等式、80638个mark样本，其中12316负K；7个八槽b0布局共735候选、168合法配对、38固定伙伴分组，factorization和中心化均精确；120个单图中心项非零，60负项保留。q=1/3、λ=1/4是独立形式变量的有理代入，不宣称满足物理q=e^(−4λ)。科学输出不覆盖；recorder验证前驱和输入hash，并绑定正式稿、编译日志与最终渲染。

下一项真实任务是charged A的forward/reverse源和AA/AB/BA共同核。旧A0–B2非零cross是必须满足的校准；不得因neutral已完成而删cross。完整memory联合核、de Sitter体质量/位置字典、finite-G四点仍未完成。Issue #1及旧KS支线保持原状态；本轮不Git、不新查评论、不外部回复。
