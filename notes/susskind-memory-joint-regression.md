# 完整联合记忆核：电荷回归、精确去静态与有限时间重求和

2026-11-07。正式稿：[九页 PDF](../output/pdf/de-sitter-susskind-memory-joint-regression.pdf)、[LaTeX](../manuscript/susskind-memory-joint-regression/main.tex)。完整论证见[联合回归证明](susskind-memory-joint-regression-proof.md)及[独立电荷 tilt 证明](susskind-memory-charge-tilt-proof.md)；来源见[literature](../literature/susskind-memory-joint-regression-sources.md)。

本轮完成指定 complex DSSYK 模型中的完整有限时间记忆协方差，包含 charged A、neutral B 及 AB/BA 交叉项。表达式是可收敛计算的谱与级数构造。前轮缺少的内部收缩电荷背景、共同样本强归约与 exact pinching 已接通。体内粒子质量、空间位置及引力相互作用的解释继续作为下一问题。

## 1. 算符、态与极限

标准 CAR，重叠的 complex Gaussian p-body 相互作用，各电荷 sector 使用同一组原始耦合。态为 beta=mu=0 的归一 annealed full-Fock trace；偶 N、p 趋无穷，lambda_N=p^2/N 趋固定 lambda>0，固定 Js，hbar=1，使用 string time。先固定场次数与实时间；不包含与 lambda 趋零或增长时间的联合极限。记 C=Js^2/lambda、q=exp(-4lambda)、X_N=(Q-N/2)/p。

对 m>0，V psi=m^(-1/2) sum_i |i> c_i psi，K=I_N tensor H_(m-1)，Q_perp=1-VV^dagger，R=Q_perp KV，K_perp=Q_perp K Q_perp。实际记忆为 M_m(s)=R^dagger exp(isK_perp)R，实际涨落场

\[
 \mathsf F_{m,N}(s)=(1-\mathcal P_h)
 m e^{-ihs/2}M_m(s)e^{-ihs/2},\qquad
 \mathcal P_h Y=\sum_eP_eYP_e.
\]

m=0 场取零，pinching 保留完整简并块。前驱已给 F=A+B+o(L2) 及普通时间级数尾；本轮在该实际接口上完成极限。

## 2. 电荷与能量须一起回归

令 D Y=QY-sum c_i^dagger Yc_i，Ncal 为 coupling Gaussian number，

\[
 G_{n,N}=[p\mathcal N-D-\lambda_N\mathcal N(\mathcal N-1)]H_N^n.
\]

倾斜源有 exact tr exp(t X_N)=cosh(t/(2p))^N。每条 internal Gaussian chord 带来 1+p(exp(t/p)-1)z，极限为 exp(tz)，而 open singleton 精确相消。未标记 singleton 修正 cosh(t/(2p))^(-2p ell)=1+O(1/p) 必须在该精确相消后再丢弃。复 t 紧盘的局部一致收敛允许 Cauchy 求导；正范数在正负实 tilt 下的指数界控制全电荷尾。

标量 charge/energy 联合律为 gamma_lambda tensor rho_q，其中 X 方差 1/(4lambda)。这是 commuting 变量的联合谱律，不能推出任意动态场独立。两 whole-field 源旧有的 4i1i2 无桥关联正是这个共同守恒 X 的协方差。

写 monic q-Hermite 展开 x^n=sum_d a_nd H_d(x)，P_d=[d]_q!、m_d=q partial_q log P_d。显式联合回归为

\[
 g_n(X,H)=\sum_d a_{nd}
 [-4\lambda m_d-\lambda(n-d)(1-2X)]H_d(H).
\]

有限 N 的 g_n(X_N,H_N) 与 H_N 精确对易，但不称其为有限 N 条件期望。G2 例中除了能量多项式，还必须减去 -2lambda C+4lambda CX；只扣能量会多留 4lambda C^2 静态方差。

## 3. 同一耦合中的强归约与静态分离

设 Z_d=W_d(H^d)，

\[
 J_{d,N}=(dp-D)Z_{d,N}
 -\lambda[d(d-1)-4m_d]H_d(H_N),\quad
 R_{n,N}=G_{n,N}-g_n(X_N,H_N).
\]

切内部 partial matching 后只留下原有 C^i q^c_gamma。中心化源为 -4lambda(inv sigma-m_d)，连通权为 lambda[2b(b-1)-4cr(B)]。mark 不对 q^c_gamma 或 a_nd 再求导。RR、RJ、JJ 三种 same-copy Gram 一起一致，因而

\[
 \|R_{n,N}-\sum_d a_{nd}J_{d,N}\|_2\longrightarrow0.
\]

固定 Gaussian degree 的乘法界允许固定左右能量权重。纯源范数为

\[
 \|J_d\|_2^2=C^dP_d\{
 \lambda[2d(d-1)-4m_d]+16\lambda^2v_d\},
 \quad v_d=(q\partial_q)^2\log P_d.
\]

无 field bridge 的图由两个独立中心化 Mahonian 和成组抵消。对外部 Hermite degrees m,n，剩余桥给有限带宽及 q^min(m,n) 衰减，故每个固定场的双能量谱具有 L2(rho_q tensor rho_q) 密度。固定 mark 的普通矩界还给两坐标的 int E^(2m) dnu_dd <= K_d R_H^(2m)，R_H=2sqrt(C/(1-q))，明确其紧支撑。正矩阵测度控制 mixed density；全谱对角线无质量。闭集 Portmanteau 得到 ||P_H R_n,N||_2 趋零，不交换无限时间平均。

## 4. 完整有限时间公式

必须先在有限 N 使用 ordinary force K_n 的可和尾、截断、减去精确守恒 k_n=g_(n+2)-2Hg_(n+1)+H^2g_n，再传递 pinching。之后才在极限 Hilbert space 重求和，得到

\[
 c_d(s)=\frac{\int e^{isE}H_d(E)d\rho_q(E)}{C^dP_d},\qquad
 D_d(s;E,E')=-c_d''+i(E+E')c_d'+EE'c_d.
\]

中性残量此处记 N_j=V_j-v_j(H)，避免与 charged R_n 混淆。用场族 Xcal_a 属于 {J_d,N_j}，幅度 f_(J_d)=D_d、f_(N_j)=+2C c_j。谱约定 E 为输入/列、E' 为输出/行：

\[
 \nu_{ab,N}(A,B)=\mathbb E\operatorname{tr}
 [\mathcal X_{a,N}^{\dagger}P_H(B)\mathcal X_{b,N}P_H(A)].
\]

则

\[
 \lim_N\mathbb E\operatorname{tr}
 [\mathsf F_N(s)^\dagger e^{iTH_N}\mathsf F_N(u)e^{-iTH_N}]
 =\sum_{a,b}\int
 e^{i(s-u)(E+E')/2+iT(E'-E)}
 \overline{f_a(s;E,E')}f_b(u;E,E')d\nu_{ab}.
\]

所有 AA/AB/BA/BB 块保留；每个矩由有限共同耦合配对求和给出，能量依赖由收敛 Hermite 展开给出。谱 total variation 的绝对尾在 compact s,u 上一致，外部实 T 相位不改变尾界。原有限 N 谱收敛仍只对固定 T 宣称；无限场求和后只保证 L1 密度，不把固定次数 L2 界外推成统一 L2。

低阶恢复 A0=J2、A1=J3-3{H,J2}/2、B1=2C N1、B2=2C[N2-{H,N1}]；A0B1 为零，而 A0B2=d0(w+eta) mu_(q^2) 严格非零。旧命名 R 与 J3 的完整 same-copy 识别仍留作独立校准，不作本定理前提。

## 5. 独立检查与保存

三份科学 JSON 均首次成功后保留，不覆盖、不重跑旧 main。

- joint-regression：68 joint X/H projections；39 two-energy layouts、117 Gram 等式；30 无桥抵消，其中28原值非零；44 mixed layouts，32非零；16 sealed-neutral engine crosschecks；18 A0B2 非零校准、12 A0B1 零校准。总93277 accepted term diagrams，34370负权及5203非零charge-connected保留。最高10 raw slots，G至5、V至4。两组有理 q/lambda 是形式独立参数。
- charge-tilt：19502 完整 support/pairing/layout 案例、315192系数；21258非零、6929负项。系统 N2p1/N3p1/N4p2，所有有序p-set assignments；16 singleton、16 collision ratios；216 complex tilts、20 full-charge G2 tilts。70/110 digits最大差约3.01e-65与3.00e-67，有限参数偏差保存，不称区间证书。
- joint-series：90普通/Hermite有理恒等式、1782 centered-force coefficients、405 c/c'/c'' 数值；physical lambda=.125/.25/.5、Js=1，五个时间，d=0..8，Jacobi 48/72。最大差1.2800225526324438e-14；解析截断尾最大4.18e-33，不含浮点舍入。环境及输入SHA见原输出。

最终主稿已双遍 pdflatex 编译，9页95dpi逐页目检，无裁切、重叠、未定义引用或公式溢出。重建命令为 `powershell -NoProfile -ExecutionPolicy Bypass -File code/build_memory_joint_regression.ps1`。输入/输出SHA、父证据、编译日志和渲染记录见 results/susskind-memory-joint-regression/local-evidence.json；该文件由只读验证器在全部文档保存后 exclusive-create。

下一物理任务是从这个核恢复正确减背景的完整 bilocal 涨落，再匹配 canonical algebra/dressing、内部空间及质量生成元。此前已复现的 FLZ/'t Hooft 谱是文献目标，不能当成从本模型新提取的质量。memory time、中心 boost time 和介子内部变量须继续区分。有限 G 四点仍开放；长期 goal 保持 active。本轮仅本地，无 Git 或外部 issue 回复。
