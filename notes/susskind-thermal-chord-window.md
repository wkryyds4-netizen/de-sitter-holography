# 热制备的大弦数窗口与 boost 单次作用

2026-11-08；正式证明见 manuscript/susskind-thermal-chord-window/main.tex，六页阅读版见 output/pdf/de-sitter-susskind-thermal-chord-window.pdf。本阶段为本地研究进展，未新增 Git 或 issue 回复。

这一轮从旧零弦种子缺陷出发，直接研究 Lin–Stanford 的 Euclidean 制备及其精确生成元。在固定 lambda 双标度、半填充的实际 charged 曲率单插入模块内，取 epsilon=4lambda、q=exp(-epsilon)、K=H/(2Js)=h/sqrt(epsilon)、r=q^d；实际通道 d=2。以下向量来自 tracial chord Hilbert space，尚未插入完整 selected rho。

## 1. 非空的可控联合窗口

端点热态 Psi=Zq^(-1/2) exp(-beta K_L/2)|d;0,0>。令 t=beta/(2sqrt(epsilon))、m=t²，正交端点基 e_n=|n,0>/sqrt([n]!)。奇偶变换后，q-Jacobi 热向量 x=exp(t Tq)e0 的系数为正，且逐项不超过普通振子 y=exp(t T1)e0。

标准 q-Wick moment 的每个 2j 配对最多 j(j−1)/2 个交叉。归一指数矩中 j~Poisson(2m)，于是 delta=2epsilon m²<1 时：

- 1−delta <= Zq/exp(2m) <=1。
- ||Psi−chi||² <=2(1−sqrt(1−delta))，chi_n=(-1)^n exp(-m/2)m^(n/2)/sqrt(n!)。
- |Psi_n|² <= Poisson(m)_n/(1−delta)。
- E_Psi(N/m−1)² <=1/[(1−delta)m]，|Psi_0|<=exp(-m/2)/sqrt(1−delta)。

充分窗口是 m→∞、epsilon m²→0，即 sqrt(epsilon)<<beta<<epsilon^(1/4)。其中弦数趋无穷，epsilon N趋零，Psi与原种子的重叠趋零。只称充分窗口，未证明最优。LS固定beta的大弦数鞍点并不自动提供这个联合极限的统一误差。

## 2. 同一Hilbert空间内的精确作用

总Hamiltonian弦数n处，c_n=q^((n+d)/2)，u_n=|n−1,1>/sqrt([n]!)。
Gram给 <e_n,u_n>=r、||u_n||²=r²+(1−r²)/[n]；f_n=(u_n−r e_n)/sqrt((1−r²)/[n])与e_n正交归一。

沿用LS的Jij=a_i†a_j−[n+d]，B=(JLL−JRR)/(2c)、P=i(JLR−JRL)/2、E=−(JLL+JRR)/(2c)，Bhat=B/sqrt(1−c²)、Phat=P/sqrt(1−c²)。取g_n=sqrt([n](1−r²))，精确作用是

\[
\widehat B e_n=\frac{[n](1-r^2)e_n-rg_n f_n}{2c_n\sqrt{1-c_n^2}},\quad
\widehat P e_n=-\frac{ig_n f_n}{2\sqrt{1-c_n^2}},
\]
\[
E e_n=\frac{[d]-(1-r)^2[n]/2}{c_n}e_n-\frac{rg_n}{2c_n}f_n.
\]

原Hamiltonian boost Btr=(K_L−K_R)/2会改变总弦数：其level-n分量是

\[
(B_{tr}\psi)_{e,n}=\frac{1-r}{2\sqrt\epsilon}
(\sqrt{[n]}\psi_{n-1}+\sqrt{[n+1]}\psi_{n+1}),\quad
(B_{tr}\psi)_{f,n}=-\frac{\sqrt{1-r^2}}{2\sqrt\epsilon}\psi_{n-1}.
\]

保留输出的最高升弦层，不将有限截断的最后一级丢弃。上述代数是已有LS representation的应用，n=0重现旧B=P=0的缺陷。

## 3. 本轮接通的限度

chi_(n−1)=−sqrt(n/m)chi_n给出符号。在该热窗口内

\[
\|(\widehat B+B_{tr})\Psi\|\to0,\qquad
\|(E-d)\Psi\|\to0,\qquad
\|(\widehat B+i\widehat P)\Psi\|\to0.
\]

Bhat、Btr均值趋零；Bhat、Phat、Btr二阶矩/方差趋d/2；< [Bhat,Phat] >趋id。证明用系数平方的共同包络 C_d(1+epsilon(n+1))²exp(2epsilon(n+d+1))，在Poisson bulk有界、尾经指数倾斜受控，再用态的范数差传图范数。并非仅由态范数接近推断无界算符接近。

同一方法给||(Tq/sqrt(m)+2)Psi||→0。因此“beta→0”不能理解成原制备的范数恢复。现在只有单次作用和弱commutator期望，不含所有后代上的commutator范数、全算符收敛、Poincare质量或完整selected矩阵元。下一步先构造移动表示的明确识别，再插入 rho=1−sqrt(q)(1+lambda)+sqrt(q)(h−sqrt(lambda)U)²；不能丢掉U源。

## 4. 计算与更正

4,152项有理精确检查：Gram symmetry 1224、adjoints 2448、endpoint Gram132、generator action144、endpoint norms144、Btr48、独立pairings/Jacobi moments12。q=1/4,4/9，d=1,2,3，n=0..7，vacuum moments到12阶。

热计算取epsilon=m^-3、beta=2/m，m=8,32,128,512；cutoff=122,226,530,1522，再增加64；Taylor阶数274,570,1451,3628。d=2的Bhat+Btr范数依次0.1717887、0.0875588、0.0440945、0.0220849；m512时Bhat²=0.9961203、Btr²=1.0000153、相反符号差范数1.9979448。最大相消指标90.48，未裁剪小分量。

第一次通用expm计算失败；把标量移位挪到指数外仍失败，故否定“只是标量移位导致”的初猜。分步/中间归一有所改善，但m128仍未过原精度，最终m512与正项级数差3.44e−9；只作诊断。主结果采用正项Taylor。

两个补充更正与初封存输出一同阅读：

1. Windows的numpy.longdouble实际nmant=52，初输出的“53/64-bit”描述不成立；不能将两种双精度求值路径说成提高精度。新增precision-checks.json用50位Decimal重算d2八个量，最大相对差3.1524e−11，发生在很小的coherent distance。旧输出原字节保留。
2. 原m512的Poisson logsf下溢使Taylor bound记录0，0不是误差证书。新增tail-bounds.json用正路径/总次数界及60位log-Chernoff覆盖全部行：scaled Taylor向量尾<=sqrt Pr(Pois(2m)>=ceil((L+1)/2))，最大8.406e−39；归一Jacobi cutoff尾<=sqrt[2Pr(Pois(2m)>=M+1)/(1−delta)]，最大1.680e−23。未重跑成功driver。

两种尾界只控制精确算术的截断，不认证所有浮点舍入。科学输出、环境和哈希均由local-evidence.json绑定；初失败脚本与诊断保留。下一轮不得重跑成功driver或recorder；纯helper可按新问题导入，必须注意其旧tail/longdouble元信息已被补充文件更正。
