# 完整双局域涨落的强重建

2026-11-07。阅读版：[七页 PDF](../output/pdf/de-sitter-susskind-bilocal-reconstruction.pdf)；[LaTeX](../manuscript/susskind-bilocal-reconstruction/main.tex)、[完整证明](susskind-bilocal-reconstruction-proof.md)、[独立高矩/插值核对](susskind-bilocal-multiplier-proof.md)、[来源](../literature/susskind-bilocal-reconstruction-sources.md)。

这轮把上一阶段的完整 memory covariance 接回原始 singlet bilocal。结果覆盖完整有限相对时间，并给实际有限 N 场到 memory 线性滤波的 O(p^(-1/3)) 强误差；带电谱幅度简化为一个 Hermite 系数，中性幅度由无重节点奇性的指数差商给出。体内空间、粒子质量和相互作用的字典尚待建立。

## 对象与精确有限起点

标准 CAR、重叠 complex Gaussian p-body Hamiltonian、所有 charge sector 共用原始 J。态为 beta=mu=0 的归一 annealed full-Fock trace；偶 N、p 趋无穷，lambda_N=p^2/N 趋固定 lambda>0，固定 Js、hbar=1、string 时间。C=Js^2/lambda，q=exp(-4lambda)。先固定时间窗，不交换增长时间或小 lambda 极限。

对 m>0，h=H_m、K=I_N tensor H_(m-1)，V psi=m^(-1/2) sum |i>c_i psi；P=VV^dagger、R=(1-P)KV、D_perp=(1-P)K(1-P)。F(t)=V^dagger exp(iKt)V 满足

\[
 F'=i\alpha hF-M*F,\quad F(0)=1,\quad
 M(s)=R^\dagger e^{iD_\perp s}R,\quad\alpha=1-p/m.
\]

目标场是 Bcal_N(t)=(1-P_h)S_c(t)，S_c(t)=m exp(-iht/2)F(t)exp(-iht/2)，没有除以 m。P_h 保留完整简并能量块；m<p（含m0）的动态场为零。

M0=P_h M、deltaM=M-M0。每个能量块配一个独立 bath 副本，R0 psi=(R P_e psi)_e、D0=directsum D_perp，自伴矩阵 K0=[[alpha h,R0^dagger],[R0,D0]] 的系统压缩即 F0。故 F/F0 均逐样本 contraction，且F0与h对易；简并块内部无需为标量。精确有序 resolvent 差式为

\[
 F-F_0=-F_0*\delta M*F
 =-F_0*\delta M*F_0+F_0*\delta M*F_0*\delta M*F.
\]

取 Q_N(s)=m deltaM(s)=exp(ihs/2) Fcal_N(s)exp(ihs/2)。前驱已证 Q_N 的 compact-time full-trace L2 一致界。

## 返回余项的定量控制

直接复用旧正力高矩：G=R^dagger R <= A_m/m，E tr_m A_m^k <=(2k-1)!!(N v_N)^k，Nv_N=O(p)。矩形 Schatten Holder 和 pinching 收缩性给 M/M0/deltaM 的 sector k范数 <=K_k p/m；本轮只用 k=2,4,8。

在 m>=N/4 上，Q_N 的二阶界给 ||deltaM||2=O(p^-2)，旧八阶粗界给 ||deltaM||8=O(p^-1)。标准非交换 Lp 插值于是给

\[
 \|\delta M\|_4\le\|\delta M\|_2^{1/3}
                         \|\delta M\|_8^{2/3}=O(p^{-4/3}).
\]

各自的精确 Volterra 方程、contraction 和固定 H 四矩给 ||F-exp(iht)||4、||F0-exp(iht)||4=O(p^-1)。一阶 exact 差式中替换两侧传播子，每项的集体二范数至多 N·O(p^-1)·O(p^-4/3)=O(p^-1/3)。bad active 区 p<=m<N/4 使用 ||Q_N||2,bad<=Kp exp(-N/16)，inactive 区为零。这不需要推断完整 m deltaM 的四阶一致界，也不需要新四标记定理。

因此对固定 S，

\[
 \sup_{|t|\le S}\|\mathcal B_N(t)-\widetilde{\mathcal B}_N(t)\|_2
 \le K_Sp^{-1/3}+K_Sp e^{-N/16},
\]

\[
 \widetilde{\mathcal B}_N(t)=
 -\int_0^t ds\int_{-(t-s)/2}^{(t-s)/2}du\,
 e^{iH_Nu}\mathcal F_N(s)e^{-iH_Nu},\qquad t\ge0.
\]

负时间按伴随。内积分来自 a+b=t-s、u=(a-b)/2，da=du，无额外2。F0的贡献精确pinch掉，而滤波后场已在P_h核中，最后的pinching可精确删除。两个场二阶范数均一致有界，故对应 covariance 的差同为 O(p^-1/3)，对所有外部实中心时间 T 一致。原 N 谱极限仍只对固定 T 断言。

## 全部谱幅度与协方差

用前驱的共同正矩阵谱 nu_ab，场族 {J_d(d>=2),N_j(j>=1)}；N_j 是 neutral residual。E为输入/列、E'为输出/行，e=(E+E')/2、omega=(E'-E)/2，chi_L=sin(omega L)/omega，chi_L(0)=L。

记 monic q-Hermite H_d 的范数 C^dP_d，P_d=[d]_q!，c_d(t)=int exp(itx)H_d(x)rho_q(dx)/(C^dP_d)。实际双局域的 limiting covariance 是

\[
 \mathcal K_{\mathcal B}(T;t,v)
 =\sum_{a,b}\int e^{iT(E'-E)}
 \overline{b_a(t;E,E')}b_b(v;E,E')\,d\nu_{ab}.
\]

charged 的 centered force 为 -g_d''-omega^2g_d，g_d=exp(-iet)c_d；d>=2 的两个初值均零。Green 积分因此精确化简为

\[
 b_{J_d}(t;E,E')=e^{-iet}c_d(t).
\]

neutral 的等价无奇性形式为

\[
 b_{N_j}(t;E,E')=\frac{2Ce^{-iet}}{C^jP_j}
 \int H_j(x)\,[E,E',x](e^{it\cdot})\,d\rho_q(x).
\]

二阶差商按连续延拓处理重复节点；其simplex表达为 (it)^2 int_(a,b>=0,a+b<=1) exp(it[aE+bE'+(1-a-b)x]) da db，故绝对值<=t^2/2，三重节点时为 -t^2exp(itE)/2。数值实现用小三角矩阵expm，避免近重节点直接分式抵消。

所有AA/AB/BA/BB保留。紧时间滤波范数<=S^2/2，继承前驱的Hilbert/谱TV可和尾；无限谱绝对连续、无对角静态质量，只称L1密度。低阶是

\[
 \mathcal B(t)=-\frac{t^2}{2}J_2
 -\frac{it^3}{6}\left[J_3-\frac32\{H,J_2\}+2CN_1\right]
 +O_{L^2}(t^4).
\]

恢复旧 -t^2Y_-/2，以及 covariance 首项 t^2v^2 K_(Y_-)(T)/4。这里展开的是已构造的极限场，不是只从有限Taylor项猜测全时间结果。

## 检查与保存

- finite-volterra-checks.json：4系统（两个非平凡CAR、m=p零记忆对照、带二维简并h的独立bath）；8Laplace、16时间积分、24contraction、12curvature。最大二阶resolvent误差2.69e-15；8/14阶时间积分交叉。t=.91的三个非平凡集体投影返回余项约.00100997/.000206864/.000184283，均保留。generic M0(.3)/M0(.8) commutator范数.0330555；错误逐特征向量pinch损失.258588，说明不能拆简并块。
- spectral-filter-checks.json：1260 charged、1440 neutral有理系数恒等式至degree14；physical lambda=.25/.5、Js=1，四个signed times、三组能量（含重合和近重合）、d1..6，共144field cases。spectral/time尺寸(40,24)/(64,36)，最大identity差2.3915e-15、refinement差2.5342e-15。另20scalar repeated-node cases对照3x3 expm/simplex积分。浮点数据不是区间证书，不作DS拟合。

两个新driver均首次成功、exclusive-create；没有运行旧main。最终七页PDF双遍pdflatex，95dpi page-1至page-7逐页目检通过，无undefined/overfull/裁切重叠。构建命令 `powershell -NoProfile -ExecutionPolicy Bypass -File code/build_bilocal_reconstruction.ps1`。参数、依赖版本、父证据与输入/输出hash由 results/susskind-bilocal-reconstruction/local-evidence.json 保存。

下一步转向这个完整bilocal的代数和物理生成元：先检验同一耦合的四点/交换子结构，明确canonical构造所需的态、极化与正规序。有限tr([A,B])=0只排除在同一tracial态下直接得到非零vacuum CCR期望，并不证明场彼此交换。随后才匹配已知Itakura/FLZ介子内部核与真正微观P^-。模型四点不能未经字典当作有限G引力相互作用。当前长期goal active，仍只本地保存。
