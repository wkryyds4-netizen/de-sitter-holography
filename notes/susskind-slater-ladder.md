# 同一 J 的 Slater 极化与完整单模 ladder

2026-11-07。沿 Susskind 静态补丁主线，从上一轮完整 bilocal 四点转入物理态/极化。此轮先读 MSS、Xu、Rajgadia–Xu 及 Itakura/Rajeev 原文，再构造实际微观算符。正式推导见[LaTeX](../manuscript/susskind-slater-ladder/main.tex)，[PDF](../output/pdf/de-sitter-susskind-slater-ladder.pdf)。详细全word证明见[proof](susskind-slater-ladder-proof.md)，文献规范和首步H计数见[interface](susskind-slater-ladder-literature-interface.md)。

## 1. 为什么选这条路径

Rajgadia–Xu 2604.14387v3 已构造适应 KM 态的 dressed boson，bare b 为canonical、实际M则含占据数dressing。此一般策略是前驱结果。本模型是守物理电荷的complex SYK，不能直接移植Majorana空态和全部mixed chord规则。只换固定半填充纯态、保持旧U(N)-covariant singlet words，也不会改annealed矩：由固定charge空间不可约及Schur引理，有限N即有 E<O|W|O>=E tr_m W。无需再对态Haar平均。

本轮让算符明确依赖参考占据投影，选原H中只制造particle–hole的块。这是具体的新计算，不重复旧Rindler midpoint或旧低三species代数。

## 2. 定义、态与极限

标准CAR；H_N=sqrt(w_N) Σ J_AB C_A† C_B，J为unit-covariance Hermitian Gaussian，包含A/B重叠，w_N=4^p Js²/lambda_N (p!/N^p)²，lambda_N=p²/N。N偶数、m=N/2、p→∞、lambda_N→固定lambda>0，Js固定，hbar=1、string energy/time。所有word长度、ladder level先固定。

固定O占据m个site，E=O补集。态phi_O=E_J<O|·|O>，非旧full-Fock trace，也未声称O是H的stationary态。令

\[
V_N^\dagger=\sqrt{w_N}\sum_{A\subset E,B\subset O;\ |A|=|B|=p}J_{AB}C_A^\dagger C_B,
\quad L_N=\binom{N-m}p\binom mp,\quad
M_N^\dagger=V_N^\dagger/\sqrt{w_NL_N}.
\]

只用原同一J，M_N=(M_N†)†，严格湮灭O，phi(MM†)=1，phi(M†M)=0。每个算符净物理charge为零，但增加/减少p个particle–hole pair。

## 3. 完整固定word定理

从右往左读词，每creation高度+1、annihilation高度−1。unbalanced或non-Dyck词在有限N为零。2k长Dyck词的全forward Gaussian配对总贡献精确为

\[
F_N(w)=\prod_{\rm closing}h\prod_{\rm opening}
R_N(h),\quad
R_N(h)=\frac{\binom{N-m-hp}p\binom{m-hp}p}{L_N}.
\]

活跃区间重叠的标签必须在particle及hole两侧都不交叠；nested与crossed都计入，偶CAR单项在合法条件下为+1。closing有h个配对选择，opening避开hp个已用site，因此得到精确乘积。固定h时R_N(h)→q^h，q=e^(-4lambda)。

backward配对在有限N不全为零。首次遇到尚未creation的标签的annihilation时，其两侧p-set必须落在此前≤kp激发site中，独立标签条件比例≤rho_N(k)=binom(kp,p)²/[binom(N−m,p)binom(m,p)]→0。至多k!个配对给总误差≤k!rho_N(k)。例phi(MM†MM†)=1+1/L_N，新增检查保留该项。

故全部固定word收敛到

\[
M|n\rangle=\sqrt{nq^{n-1}}|n-1\rangle,
\quad M^\dagger|n\rangle=\sqrt{(n+1)q^n}|n+1\rangle.
\]

creation范数s_k,N=k!Π_(j=0)^(k−1)R_N(j)在有限N精确，趋k!q^[k(k−1)/2]。完整word还给实际规范化creation向量的annihilation leakage→0（固定level、annealed向量范数），不只是纯creation Gram拟合。

## 4. canonical化与Hamiltonian检验

极限有限粒子core上M=q^(n/2)b、M†=b†q^(n/2)，[b,b†]=1；逆b=q^(−n/2)M无界。更具体可取实际n_N=Q_E/p、q_N=e^(−4lambda_N)，定义b_N=q_N^(−n_N/2)M_N及其伴随。每一步的height是精确occupation本征值，故Dyck word的微观b矩精确等于q_N^(−Σ_closing(h−1))乘M矩。全部固定word因此趋canonical真空矩，任意固定非交换多项式P满足E||([b_N,b_N†]−I)P(b_N,b_N†)O||²→0。这是指定Slater态的固定激发core结果；全空间inverse不一致有界，未声称operator-norm或旧trace收敛。

M是compact weighted shift；ker M=真空一维，谱投影P0及(M†)^nP0M^j生成全部matrix units，故W*(M)=B(ell²)。同一微观M在旧normalized full-Fock迹下E tr(M†M)=4^(−p)精确，在新Slater态phi(MM†)=1。这是新适应代数，非旧tracial代数内部的直接拆分。另两副本Wick中跨副本标签必须由其余sites覆盖，给fixed scalar word方差≤(2k)!rho_N(2k−1)→0；annealed主定理不依赖此自平均附加结论。

完整H|O>产生k个particle–hole的权重严格为W_k=w_N binom(m,p)binom(p,k)binom(N−m,k)。总C_m=w_N binom(m,p)binom(N−m+p,p)→C=Js²/lambda。r=p−k服从超几何分布

\[
P_N(r)=\frac{\binom pr\binom{N-m}{p-r}}{\binom{N-m+p}p},
\quad E(r)_{\underline t}=\frac{(p)_{\underline t}^2}{(N-m+p)_{\underline t}}
\to(2\lambda)^t.
\]

阶乘矩由(p²/(N−m))^t控制，正指数矩一致有界，故r→Poisson(2lambda)。sigma_N²/C_m→e^(−2lambda)=sqrt(q)是捕获的**范数平方比例**；规范化overlap amplitude为e^(−lambda)=q^(1/4)。

在joint disorder/state Hilbert空间K_N=closure span{(M†)^j O}内，精确P_K H O=sigma M†O，余项范数平方=C_m−sigma²→C(1−sqrt(q))>0。因此不能把完整H直接认成单ladder的sigma(M+M†)。下一步保留缺失通道，计算同J joint(H,M,M†)律。

## 5. 检查、旁证及下一步

新driver code/check_slater_ladder.py首跑成功后封存：1,537,310个Gaussian label assignment；754,030个forward CAR/interval比较；46完整Dyck词；136个non-Dyck零词；25个非零backward有限修正；4个完整H首步枚举；6个超几何渐近诊断行；189个old mixed-Q线性四cumulant恒等式。参数/版本/输入哈希在results/susskind-slater-ladder/exact-checks.json。有限fixtures核身份/符号，不代替渐近证明；浮点TV诊断不作interval证书。

低旧三场Y=ah+bu+cx的κ4=−(1−q)A²−2(1−q²)AB−(1−q⁴)B²，A=a²+b²、B=c²，是旧mixed-Q律的标准后果，说明线性whitening不能得到quasifree Gaussian场。不把此旁证当本轮主成果。另存[KMS接口](susskind-bilocal-kms-polarization.md)，区分joint annealed Gibbs、exact thermofield与仅两点Gaussian dilation。

下一首任务：保留同J的r≥1首步缺额，构造joint(H,M,M†)首批interleaved四点及一般配对，求正表示与动力学闭包；再与已有真实bilocal共同矩接线。stationary态/真实flow、空间与内部生成元、meson质量及finite-G动力学仍开放。已有旧科学输出不重跑、不覆盖；本轮只本地，不Git、不issue回复。
