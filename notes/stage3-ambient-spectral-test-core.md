# 环境位置测试的完整谱核心与强轨道逼近

2026-10-29。本文补足生成FP空间与最大位置态类比较所需的一项分析输入：对固定角向模式，全部连续谱双通道加有限束缚态能够在明确的加权位置范数中逼近紧支撑测试；相应产品的群平均像在完整环境FP截面空间中总性成立。这里连续谱使用全部s>0，尚不施加实际observer硬下限s_*。

采用dS₂、L=1、proper G=SO⁺(2,1)，dx=cosh(t)dt dθ及既有六轨道测度。物质K可为任意可分Hilbert空间上的强连续酉表示V。已有输入是[径向散射稿](stage3-radial-scattering-core.md)R2—R16及[共同域稿](stage3-positive-frequency-orbit-core.md)P14、P16—P18。本文不以运动学L²稠密代替FP稠密，也不从有界运动学投影直接推出分布上的谱投影连续。

## 1. 一个较弱、但仍充分的加权范数

对单角向模式写

\[
u(t,\theta)=\frac{e^{in\theta}}{\sqrt{2\pi\cosh t}}Q(t),\qquad
b_\alpha(u):=\|Q\|_2+M_\alpha(Q),\qquad
M_\alpha(Q)=\int_{\mathbb R}(1+|t|)^\alpha|Q(t)|dt.
\tag{E1}
\]

固定1/4<α<1即可，以下可以统一取α=1/2。bα是固定n径向空间上的范数；其表达式不依赖n。

令f=|u|。它是K=SO(2)不变位置函数，且其H₀径向变换恰为|Q|。由R9，任意0<α≤1有

\[
|e_{0,\pm}(t,s)|\le C\min\{1,[s(1+|t|)]^\alpha\}.
\]

因此f的径向ds谱密度q_f满足

\[
q_f(s)\le C M_\alpha(Q)^2\min(s^{2\alpha},1),\qquad
\int_0^\infty q_f(s)ds=\|Q\|_2^2.
\tag{E2}
\]

注意这里的s是取模后H₀的谱变量，不是原u的质量；即使u有正质量下限，f也不能据此截去零谱。

群Schur权K_P(s)=4π²/[s tanh(πs)]在零点为O(s⁻²)，在s≥1有界。于是

\[
\begin{aligned}
\int_0^\infty K_Pq_f^2ds
&\le C M_\alpha(Q)^4\int_0^1s^{4\alpha-2}ds
 +C M_\alpha(Q)^2\|Q\|_2^2\\
&\le C_\alpha b_\alpha(u)^4.
\end{aligned}
\tag{E3}
\]

阈值α>1/4正来自第一个积分。按R15—R16同时在L¹(ds)及L²(K_Pds)中逼近q_f，得到合法的球函数L²延拓，故

\[
\|A_{|u|}\|_{L^2(G)}\le C_\alpha b_\alpha(u)^2,
\qquad A_f(g)=\langle f,U_X(g)f\rangle.
\tag{E4}
\]

P14的非负Tonelli恒等式随即给单模式产品的强轨道界：

\[
\boxed{\|a_{u_Lu_R\otimes v}\|_{L^2(Z)}
\le C_\alpha b_\alpha(u_L)b_\alpha(u_R)\|v\|_K.}
\tag{E5}
\]

E4—E5的常数不依赖两个磁指标，因为这里应用的始终是取模后的H₀径向定理。物质向量仅需属于K，不要求其矩阵系数具有统一球函数界。产品的运动学L²范数也由bα乘积控制，故这些产品属于D_orb。

同一磁指标下的差仍是单模式函数。因此若u_{i,j}→u_i在bα中成立，产品差拆成两个项后，E5证明其在强轨道半范数中趋零。有限产品和同理；不需要无限磁指标的一致散射估计。

## 2. 固定磁指标的完整径向谱

令k=|n|，

\[
H_k=-\partial_t^2+(1/4-k^2)\operatorname{sech}^2t,
\qquad D(H_k)=H^2(\mathbb R).
\tag{E6}
\]

H₀的无束缚态、无零共振及双通道酉ds变换已在径向稿中证明，所用外部输入为Germain–Pusateri–Rousset [Theorem 2.3](https://ems.press/content/serial-article-files/15932)。该文的无束缚态假设不直接套到k≥1；这些扇区用以下因式分解补齐。

对k≥1设ℓ_k=k−1/2，A_k=∂_t+ℓ_k tanh t，定义域H¹。其闭算符伴随是A_k†=−∂_t+ℓ_k tanh t。直接计算给

\[
A_k^\dagger A_k=H_k+\ell_k^2,\qquad
A_kA_k^\dagger=H_{k-1}+\ell_k^2.
\tag{E7}
\]

归纳时H_{k−1}+ℓ_k²严格正：k=1时下界为1/4；k≥2时其最低本征值加ℓ_k²为2k−2。因此

\[
V_k=A_k^\dagger(H_{k-1}+\ell_k^2)^{-1/2}
\tag{E8}
\]

是L²到ker(A_k)⊥的等距满射。具体而言，A_kA_k†严格正使A_k†闭值域，而其值域的正交补为ker A_k；极分解给E8。E7及闭算子的函数演算给H_kV_k=V_kH_{k−1}的谱交织。另有

\[
\ker A_k=\operatorname{span}\{\operatorname{sech}^{\ell_k}t\},
\quad H_k\operatorname{sech}^{\ell_k}t=-\ell_k^2\operatorname{sech}^{\ell_k}t.
\]

迭代得到完整酉谱分解

\[
\boxed{H_k\simeq H_0\oplus\bigoplus_{\Delta=1}^k
\{-({\Delta}-1/2)^2\}.}
\tag{E9}
\]

每个负本征值单重。其归一化本征函数记b_{k,Δ}；先从sech^(Δ−1/2)t出发，再施加后续A_j†并除以相应正常数即可构造。它们及固定阶导数在两端指数下降，所以属于全部Mα空间。Casimir本征值为1/4−(Δ−1/2)²=Δ(1−Δ)。这里列出的是固定磁指标中的束缚态，不把它们误作额外的principal质量。

## 3. Jost归一化与零能非共振

取F_{k,+}(t,s)∼e^{ist}于t→+∞，F_{k,−}(t,s)=F_{k,+}(−t,s)。由[Ferrers超几何表示](https://dlmf.nist.gov/14.3#E1)及[连接公式](https://dlmf.nist.gov/15.8#E4)，

\[
F_{k,+}(t,s)=e^{ist}{}_2F_1\left(\tfrac12-k,\tfrac12+k;1-is;
\tfrac{1-\tanh t}{2}\right),
\]

\[
T_k(s)=\frac{\Gamma(1/2-k-is)\Gamma(1/2+k-is)}
 {\Gamma(1-is)\Gamma(-is)},\qquad
R_k(s)=\frac{-i(-1)^k}{\sinh\pi s}T_k(s).
\tag{E10}
\]

偶势的左右反射幅相同，R_k按TF_{k,+}=overline(F_{k,−})+R_kF_{k,−}固定。Gamma递推给

\[
T_k(s)=T_0(s)(-1)^k\prod_{j=1}^k
\frac{j-1/2-is}{j-1/2+is}.
\tag{E11}
\]

故|T_k|²=tanh²πs、|R_k|²=sech²πs。再用Γ(1/2−k)Γ(1/2+k)=(−1)^kπ，有

\[
T_k(s)=-i(-1)^k\pi s+O(s^2),\qquad R_k(0)=-1.
\]

沿用W_k=F_{k,−}F'_{k,+}−F'_{k,−}F_{k,+}和T_k=2is/W_k，得到

\[
W_k(0)=\frac{2(-1)^{k+1}}{\pi}\ne0.
\tag{E12}
\]

所以固定k无零能共振。这里不能只凭“半整数势参数”略去判断，也不能把H₀的势非负论证照搬到负势。

定义ds归一化入射模式

\[
e_{k,\pm}(t,s)=(2\pi)^{-1/2}T_k(s)F_{k,\pm}(t,s),\quad s>0.
\tag{E13}
\]

这些模式确实给E9连续部分的酉变换，而非仅是两个ODE解。为核对归一化，由Jost端点条件，

\[
F_{k,+}=\frac{A_k^\dagger F_{k-1,+}}{\ell_k-is},\qquad
F_{k,-}=\frac{A_k^\dagger F_{k-1,-}}{is-\ell_k}.
\]

结合E11，e_{k,+}和e_{k,−}分别是V_k作用于前一级广义谱模式后，乘以−(ℓ_k−is)/sqrt(ℓ_k²+s²)和+(ℓ_k−is)/sqrt(ℓ_k²+s²)。这两个乘数均为单位相位；从H₀酉变换归纳即得

\[
a_\pm(s)=\int\overline{e_{k,\pm}(t,s)}Q(t)dt,
\]

\[
Q=\sum_{\Delta=1}^k\langle b_{k,\Delta},Q\rangle b_{k,\Delta}
 +\sum_{\pm}\int_0^\infty e_{k,\pm}(t,s)a_\pm(s)ds,
\tag{E14}
\]

以及对应Parseval公式。E14首先在L²成立；以下对Cc∞输入建立更强的位置逼近。谱测度仍是ds，未加入群Plancherel权。

## 4. 好半轴上的均匀Jost余项

对t≥0写F_{k,+}=e^{ist}m_+(t,s)。Volterra方程为

\[
m_+(t,s)=1+\int_t^\infty D_s(u-t)V_k^{pot}(u)m_+(u,s)du,
\quad D_s(v)=\int_0^v e^{2isr}dr,
\]

其中V_k^{pot}=(1/4−k²)sech²，避免与E8的等距算符混用符号。对实s≥0，

\[
|\partial_s^jD_s(v)|\le\frac{2^j}{j+1}v^{j+1}.
\]

势及其矩指数可积。先用Volterra迭代估计j=0，再逐阶对参数微分：最高阶导数仍由相同Volterra算子控制，低阶源项由上述矩界控制。因为对每个固定p，∫_t^∞(u−t)^p|V_k^{pot}(u)|du≤C_{k,p}e^(−2t)，得到

\[
\sup_{s\ge0}|\partial_s^j(m_+(t,s)-1)|\le C_{k,j}e^{-2t},
\quad t\ge0.
\tag{E15}
\]

j≥1时减去1不影响导数。可先在有限s区间对Volterra级数微分，再用不含s的上界延拓；不在s=0使用含1/s的粗估计。m_−在t≤0满足同样结论，因势偶对称。

E10—E11在s=0具有可去奇点，在正实轴没有极点。Gamma比的[Stirling展开](https://dlmf.nist.gov/5.11#iii)及其解析导数控制表明：T_k、R_k在[0,∞)光滑，各阶s导数有界。高能区也可从E11的有限有理乘积和T₀的Gamma比单独验证。于是对t≥0，

\[
\sqrt{2\pi}\,e_{k,+}=T_ke^{ist}m_+,
\quad
\sqrt{2\pi}\,e_{k,-}=e^{-ist}\overline{m_+}+R_ke^{ist}m_+.
\tag{E16}
\]

负半轴交换±处理。因此e_{k,±}全t、实s≥0有界，固定阶s导数至多增长为C_{k,j}(1+|t|)^j；在任意紧t区间它们及全部s导数一致有界。这里常数可依赖固定k；后文从未交换无限k极限与这些常数。

## 5. 紧支撑输入的谱系数端点性质

令Q∈Cc∞(ℝ)。在其紧支撑上，F_{k,±}在s=0光滑，T_k(s)=O(s)。因此E14的谱系数满足

\[
a_\pm(s)=s\,\beta_\pm(s),\qquad \beta_\pm\in C^\infty([0,s_0]).
\tag{E17}
\]

这项零点来自非共振散射归一化，不能对任意半轴Fourier振幅预设。

另一方面，H_k^MQ仍光滑且支撑于同一紧集。对s≥1，由ODE及紧支撑分部积分，

\[
a_\pm(s)=s^{-2M}\int\overline{e_{k,\pm}(t,s)}(H_k^MQ)(t)dt.
\]

E16的紧t导数界允许任意固定阶s微分。因此对所有j,L，

\[
\sup_{s\ge1}s^L|\partial_s^ja_\pm(s)|<\infty.
\tag{E18}
\]

E17和E18既不要求Q已在连续谱子空间，也不遗漏其束缚投影；束缚部分由E14单独保留。

## 6. 低谱截断在Mα中的真实速率

固定χ∈C∞([0,∞))，χ=1于[0,1]，χ=0于[2,∞)，并令0<ε<1/2。定义低谱部分

\[
L_\epsilon Q(t)=\sum_\pm\int_0^\infty
 e_{k,\pm}(t,s)\chi(s/\epsilon)a_\pm(s)ds.
\tag{E19}
\]

在每个半轴，E16把它写成有限个普通振荡积分加Jost余项。每个普通积分的振幅具有sχ(s/ε)B(s)形式，B及其前两阶导数在固定小谱区间有界；含T_k的项还多一个s，但不需要利用这个额外零点。

设cε(s)=sχ(s/ε)B(s)。它满足

\[
c_\epsilon(0)=0,\quad \|c_\epsilon\|_1\le C\epsilon^2,
\quad |c_\epsilon'(0)|+\|c_\epsilon''\|_1\le C.
\tag{E20}
\]

第一种估计直接给积分O(ε²)。对t≠0分部积分两次，s=∞端项为零，s=0第一端项因cε(0)=0消失；第二次留下cε'(0)/(it)²，与∥cε''∥₁共同给O(|t|⁻²)。E15使Jost余项由Cε²e^(−2|t|)控制。故

\[
|L_\epsilon Q(t)|
\le C\min(\epsilon^2,|t|^{-2})+C\epsilon^2e^{-2|t|}
\le C'\frac{\epsilon^2}{(1+\epsilon|t|)^2}.
\tag{E21}
\]

最后一步对ε≤1使用e^(−2|t|)≤C(1+ε|t|)⁻²。常数可依赖Q,k,χ，与ε、t无关。这里的联动大时间区域已被控制，不是固定t后取ε→0。

于是对0≤α<1，

\[
M_\alpha(L_\epsilon Q)\le C_\alpha\epsilon^{1-\alpha},
\qquad \|L_\epsilon Q\|_2\le C\epsilon^{3/2}.
\tag{E22}
\]

L²界也直接来自Parseval与a±=O(s)。Mα界由变量r=ε|t|及∫₀∞(1+r)^(α−2)dr<∞给出。取α=1/2，低端bα误差为O(ε^1/2)。这解释了为什么原M₁条件过强，而本稿需要1/4<α<1。

## 7. 高谱截断与完整逼近

令R>1，定义高谱部分

\[
H_R^{tail}Q(t)=\sum_\pm\int_0^\infty
 e_{k,\pm}(t,s)[1-\chi(s/R)]a_\pm(s)ds.
\]

E18及T_k、R_k的导数界说明：上述每个半轴平面波振幅及其任意固定阶s导数的L¹范数为O(R⁻ᴸ)，任意L均可。振幅在s=0邻域为零，在∞快速下降，故对s任意多次分部积分没有边界项。E15控制其Jost余项。因此对任意L,N，

\[
|H_R^{tail}Q(t)|\le C_{L,N}R^{-L}(1+|t|)^{-N}
 +C_LR^{-L}e^{-2|t|}.
\tag{E23}
\]

选择N>α+1即得Mα(H_R^{tail}Q)→0；L²亦趋零。这里先固定Q,k，随后R→∞，不主张对所有磁指标一致。

保留全部有限束缚项，并对连续谱使用中间滤波，定义

\[
Q_{\epsilon,R}=\sum_{\Delta=1}^k
\langle b_{k,\Delta},Q\rangle b_{k,\Delta}
+\sum_\pm\int_0^\infty e_{k,\pm}(t,s)
[1-\chi(s/\epsilon)]\chi(s/R)a_\pm(s)ds.
\tag{E24}
\]

若R>2ε，其与Q的差在L²中恰为LεQ+H_R^{tail}Q；E21及E23使同一差也是所写通常位置函数。连续谱系数属于Cc∞((0,∞))。由E22—E23，

\[
\boxed{\|Q-Q_{\epsilon,R}\|_2+M_\alpha(Q-Q_{\epsilon,R})\longrightarrow0}
\tag{E25}
\]

于ε↓0、R↑∞成立，对任意0≤α<1。连续谱紧包本身在两端快速下降，束缚项指数下降，所以每个近似者的bα有限。

## 8. 可改用原KL正、负频双族，但不预设其正交

本段仅说明测试核心可换基，不用未验证的形式δ正交关系。KL [原(4.9)](https://arxiv.org/html/2411.08090v1#S4.E9)定义tildeΨ_n(t,θ)=(−1)^nΨ_n(−t,θ)。按共同域稿P7，固定k=|n|、h=k+1/2，可把双端系数写为

\[
A_{\epsilon\sigma}(s)=B_\sigma(s)
 e^{-i\epsilon\pi h/2}e^{-\epsilon\sigma\pi s/2},
\]

其中Bσ光滑且在s>0无零点。于t→+∞，以e^(±ist)为行、以原Ψ和tildeΨ为列的2×2系数矩阵行列式是

\[
(-1)^n[-2\sinh(\pi s)]B_+(s)B_-(s)\ne0.
\tag{E26}
\]

角向及径向酉归一化只给共同非零常数，不影响可逆性。两族在每个正s张成同一个二维ODE解空间，因此与E13入射双通道之间有光滑可逆矩阵，且该矩阵及逆矩阵在任意紧正谱区间有界。故E24中的任一Cc∞双通道包可精确重写为一个原KL正频Cc∞包与一个原KL负频Cc∞包之和。

E26只证明可逆换基。要把这些族解释成环境Hilbert空间的两个正交群不变扇区，还需独立核对其谱Gram归一化；本稿的环境总性不依赖这项正交性。这里的“正频”始终是原KL模式选择，未改成时间Fourier半轴。

## 9. 产品测试的强轨道逼近与环境FP总性

令S_all为以下产品的有限线性包：两个observer各取固定n的Cc∞((0,∞))连续谱双通道包，或固定n中任一束缚本征模式；物质取任意v∈K。也可依E26把连续双通道换成原KL±双族。E5保证S_all⊂D_orb。

对任意φ∈Cc∞(X²;K)，先在一个固定紧集上以有限物质向量、有限角向Fourier和及两个时间变量的紧支撑光滑张量和一致逼近φ。此处可先以有限秩正交投影一致逼近紧集上的紧像φ(C)，再用标量光滑张量逼近；所有时间cutoff保持在一个共同紧集。P16将这个一致逼近转为a_(误差)的L²收敛，运动学L²收敛亦成立。

每个有限张量因子只有固定磁指标及Cc∞径向函数。E25逐因子提供全谱近似，再用E5及产品差分解可得ψ_j∈S_all，使

\[
\|\phi-\psi_j\|_{L^2(X^2;K)}
+\|a_{\phi-\psi_j}\|_{L^2(Z)}\longrightarrow0.
\tag{E27}
\]

先选有限张量精度，再对其有限个径向因子选ε、R，最后对精度对角取列即可；从未要求对无限磁指标统一选谱cutoff。

记Bψ(z)=∫_G V(g⁻¹)ψ(gs(z))dg。Bochner三角不等式给∥Bψ∥₂≤∥aψ∥₂。故E27特别给Bψ_j→Bφ于FP范数。

最后证明环境总性。若k∈L²(Z;K)与全部B(S_all)正交，E27使其与全部B(Cc∞)正交。由P17的绝对配对，

\[
0=\langle k,B\phi\rangle_{L^2(Z;K)}
=\int_{X^2}\langle Rk(y),\phi(y)\rangle_Kdy
\quad\forall\phi\in C_c^\infty(X^2;K).
\]

P18说明Rk∈L¹loc并且R是注入的，因此先得Rk=0分布，再得k=0。于是

\[
\boxed{\overline{B(S_{all})}^{\,FP}=L^2(Z;K).}
\tag{E28}
\]

这是全连续双族及全部束缚测试共同生成的环境空间；尚未声称其中仅正频双observer部分已经等于全部环境空间。把某个谱约束定义为对互补测试族消失，并用扇区正交证明其等于相应生成闭包，是后续独立步骤。E28为该步骤提供的是非循环的完整测试族。

## 10. 硬质量下限与结论范围

E25不能直接改成在固定正点s_*附近删去窄带的同范数逼近。零端成功使用了a±(s)=O(s)；若在s_*处振幅不为零，宽ε窄带的普通Fourier主项尺度是ε·f(εt)e^(is_*t)，其Mα尺度一般为ε^(−α)，并不趋零。即使窄带的运动学L²范数趋零，也不能据此丢掉其强轨道贡献。

因此本稿先在完整s>0上证明环境测试定理。实际observer下限应在已识别的正频物理谱L²空间中，作为相应质量坐标的谱支持条件处理；绝不能把上段无效的环境切孔当作证明。若另行允许集中于硬边界的分布系数，需单独检验其FP有限性和谱约束含义，不能默认将其纳入或删去。

本文是既有一维散射、Darboux分解与轨道Fubini理论在当前位置FP问题中的具体应用。解析结论为E5、E9、E25、E27和E28；没有运行或以有限点数值替代这些证明，也未声称解决有限G、未调节位置核或完整de Sitter全息对偶。