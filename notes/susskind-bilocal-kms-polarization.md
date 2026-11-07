# Bilocal 的 Gibbs 倾斜、能隙极化与两种 doubled representation

2026-11-07，本地研究附录。采用项目既定 common-J complex DSSYK、normalized annealed full-Fock trace τ_N、固定 λ>0、string H/time 单位。输入为完整 bilocal joint moments 及所有固定 Schatten 阶数的紧时间尾。

**本稿构造的是新态下的可控诊断。** 下文 Gibbs 态明确采用 annealed joint normalization，不能改写为 quenched/samplewise normalization，也不将它认作 Susskind 的 bulk Rindler 态。已封存 [Rindler packet](susskind-rindler-packet.md) 已给 midpoint 两点构造、等 β/4 四点的 33≠3 及随插入数变化的热夹心反例；这里不重复这些结果。

## 1. 可以传递到实际有限 N 的固定 β 态

令 τ_N(Y)=E_J tr_N(Y)，定义
\[
 Z_{\beta,N}=\tau_N(e^{-\beta H_N}),\qquad
 \phi_{\beta,N}^{\rm ann}(Y)=
       \frac{\tau_N(e^{-\beta H_N}Y)}{Z_{\beta,N}}.        \tag{1}
\]
它不同于
\[
 E_J\frac{\operatorname{tr}_N(e^{-\beta H_N}Y)}
                 {Z_{\beta,N}(J)},\qquad
 Z_{\beta,N}(J)=\operatorname{tr}_N(e^{-\beta H_N}).       \tag{2}
\]
(1) 相当于同时按 Z_{β,N}(J) 倾斜 disorder 分布。(2) 的极限需要 partition-function 集中或额外 multi-trace/reciprocal 控制，本稿不由单迹 joint moments 推出它。

不能仅凭每个固定 S_r 矩有界便交换 Gibbs 指数。这里有更强输入
\[
 \sup_N\|H_N^n\|_r\le K_r^n\sqrt{n!},
\]
故对固定 B、r，
\[
 \sup_{N,|\beta|\le B}
 \left\|e^{-\beta H_N}-\sum_{n=0}^M
                    \frac{(-\beta)^n}{n!}H_N^n\right\|_r
 \le\sum_{n>M}\frac{(BK_r)^n}{\sqrt{n!}}\to0.            \tag{3}
\]
这同时给所需统一指数矩。因 τ_N(H_N)=0，Jensen 给 Z_{β,N}≥1。旧能量矩极限和 (3) 给
\[
 Z_{\beta,N}\to Z_\beta=\int e^{-\beta E}d\rho_q(E),
 \qquad\phi_\beta(Y)=\frac{\tau(e^{-\beta H}Y)}{Z_\beta}. \tag{4}
\]
对任意固定 k 个 actual bilocal、其伴随或其已证局域场，先以 (3) 截断 Gibbs 因子，对有限能量插入使用共同矩极限，再由 Hölder 和更高固定 S_r 去截断。actual-to-proxy 的强差也以同样 Hölder 控制。因此全部固定 k 点，尤其完整四点，在 (1) 下收敛到 (4)，对 β 和各时间的固定紧集局部一致。

原 actual bilocal 仍满足 exact energy pinching 为零，所以 φ_{β,N}^{ann}(B_N)=0；对其极限也是如此。改变态不需要改写原 finite-N 算符。

## 2. KMS、谱权及详细平衡的方向

统一取 \(\alpha_t(X)=e^{itH}Xe^{-itH}\)。有限 N 的 trace cyclicity 直接给
\[
 F_{AB}(t)=\phi_\beta(A\alpha_t(B)),\qquad
 F_{AB}(t+i\beta)=\phi_\beta(\alpha_t(B)A).               \tag{5}
\]
固定复时间条带的能量指数由 (3) 控制，故该边界关系传给这里的 analytic moment domain。标准有限代数的 KMS 与 Hilbert–Schmidt 实现可核于 [Dereziński, §§8.3--8.6](https://arxiv.org/html/math-ph/0511030v2#S8.SS3)。

使用项目原谱约定：E 为 input/column，E′ 为 output/row，ω=E′−E。若
\[
 d\nu_{ab}(E,E')\quad\text{表示}\quad
 \tau[X_a^\dagger P_H(dE')X_bP_H(dE)],
\]
则 Gibbs 谱精确为
\[
 d\nu_{ab}^{\beta}(E,E')=
              Z_\beta^{-1}e^{-\beta E}d\nu_{ab}(E,E').   \tag{6}
\]
注意权重乘在 input 能量 E，不是任意选 E′ 或平均能量。

对自伴 channel family（非自伴 bilocal 可先取实、虚 quadratures），定义
\[
 C_{ab}^\beta(t)=\phi_\beta(X_a\alpha_t(X_b))
           =\int e^{it\omega}\,dS_{ab}^\beta(\omega).
\]
循环迹和 (6) 给
\[
 \boxed{dS_{ab}^\beta(-\omega)
              =e^{-\beta\omega}dS_{ba}^\beta(\omega).}  \tag{7}
\]
这里正 ω 是从 input 向较高 output 能量的矩阵元；改变 Fourier 符号时必须同时改变这句话。

对单个自伴 X，
\[
 \phi_\beta([X,\alpha_t(X)])
 =\int e^{it\omega}(1-e^{-\beta\omega})dS_X^\beta(\omega).
\]
β=0 时平均交换子为零，是 trace cyclicity；它不推出交换子算符为零。项目旧 J_2/W_1 交换子平方严格正已经提供实际反例。

## 3. 精确 thermofield/Liouville 实现保留所有四点

有限 N 取 joint Hilbert–Schmidt 空间
\[
 \mathcal K_N=L^2(d\mu_J;\operatorname{HS}_{\rm normalized}),
 \qquad\Omega_{\beta,N}=Z_{\beta,N}^{-1/2}e^{-\beta H_N/2}.
\]
定义左乘 \(\mathsf L(A)Y=AY\)、右乘 \(\mathsf R(A)Y=YA\)，后者作为 opposite algebra 的表示；左右相互对易。于是
\[
 \langle\Omega_\beta,\mathsf L(A_1)\cdots
             \mathsf L(A_k)\Omega_\beta\rangle
       =\phi_\beta(A_1\cdots A_k).                       \tag{8}
\]
这个 doubling 是精确的，保留非 Gaussian 四点。其时间生成元、模算符为
\[
 \mathscr L=\mathsf L(H)-\mathsf R(H),\qquad
 \Delta_\beta=e^{-\beta\mathscr L},\qquad
 \sigma_u=\alpha_{-\beta u}.                            \tag{9}
\]
\(\mathscr L\) 的本征值是能隙 E′−E，通常有正有负；它不是正能量的 bulk Hamiltonian。β=0 时 Δ=1。

能隙正负分解 \(1_{\mathbb R_\pm}(\mathscr L)X\Omega_\beta\) 是明确的二范数极化。它只是 operator-state 向量分解，不等于已经得到微观 annihilation/creation operators 或 CCR。此标准形式与 [Dereziński §8.5--8.6](https://arxiv.org/html/math-ph/0511030v2#S8.SS5) 一致。

## 4. 可直接用于实际四点的平滑正频率波包

为了将极化传到高范数和实际有限 N 四点，优先用平滑频带。取
\(f\in C_c^\infty((\omega_0,\omega_1))\)、0<ω_0<ω_1，并用约定
\(f(\omega)=\int k_f(t)e^{it\omega}dt\)。定义
\[
 X_{f,N}=f(\operatorname{ad}H_N)X_N
       =\int k_f(t)\alpha_t(X_N)dt.                     \tag{10}
\]
有限 N 精确等于
\(\sum_{E,E'}f(E'-E)P_{E'}X_NP_E\)。由于 k_f∈L¹，
\[
 \|X_{f,N}\|_r\le\|k_f\|_1\|X_N\|_r.                 \tag{11}
\]
这给 actual fields 的全固定高矩与四点传递：先截断 t 积分到紧区间使用已有时间极限，再由 k_f 的 L¹ 尾去截断；unitary conjugation 的范数与 t 无关。额外能量窗也可用有 Fourier-L¹ 表示的光滑函数从左右乘入。

由 (6) 可计算波包 Gram
\[
 \phi_\beta(X_{a,f}^\dagger X_{b,g})
       =\int\overline{f(\omega)}g(\omega)
                         dS_{ab}^{\beta}(\omega).        \tag{12}
\]
选有限 channel bank、频带及 f，先对其正 Gram 的非零谱子空间正交归一，再计算旧共同图定义的四点余项，便得到一个可操作的 mode basis。精确 sharp projector 在 S₂ 上收缩；本稿不把任意 sharp spectral multiplier 自动当作 S₄ 收缩。

## 5. 同一个二点核的辅助 canonical thermal Gaussian dilation

另一个构造只复现二点。对 β>0，在正频率定义正矩阵测度
\[
 dJ_\beta(\omega)=(1-e^{-\beta\omega})dS_\beta^+(\omega),
 \qquad n_\beta(\omega)=(e^{\beta\omega}-1)^{-1}.          \tag{13}
\]
以 dJ_β 为 one-particle 内积、h 为乘 ω，构造 bosonic doubled Fock 空间。在 band-limited domain 上可写
\[
 b_\beta(f)=a_1(\sqrt{1+n_\beta}\,f)
                 +a_2^\dagger(\sqrt{n_\beta}\,\bar f),   \tag{14}
\]
其中第二个 Fock 因子是 conjugate one-particle space。它满足辅助 CCR，真空中的 b_β 二点含 1+n_β 与 n_β。由
\[
 (1+n_\beta)dJ_\beta=dS_\beta^+,
 \qquad n_\beta dJ_\beta=e^{-\beta\omega}dS_\beta^+,
\]
其自伴场 \(b_\beta+b_\beta^\dagger\) 精确重现 (7) 的二点，且标准 Liouvillean 为 \(d\Gamma(h\oplus-\bar h)\)。这是 [Dereziński §9.3, Theorem 9.4](https://arxiv.org/html/math-ph/0511030v2#S9.SS3) 的 Araki–Woods 构造在本谱上的应用。ω→0 时须保留 \(\operatorname{Dom}\sqrt{n_\beta}\)；固定正频带避免这个域问题，不能无条件把 dJ 的逆平方根延伸到零频率或零 Gram 特征值。

更抽象地，对实 channel space 令 C_β(f,g)=φ_β(X_fX_g)，V=Re C_β、σ=2 Im C_β。C_β≥0 给 Gaussian Weyl 态存在所需的不确定性条件。其 CCR 是 \([\Phi(f),\Phi(g)]=i\sigma(f,g)I\)，但这是**辅助表示的定义**，不是由实际 X 的平均交换子反推出的算符恒等式。

β=0 时 σ=0，该二点 Gaussian completion 退化为 commuting Gaussian 模型；若另选非零 symplectic form，就改变了二点 ordered kernel。正频率分解仍能在 §3 的真实 Liouville 空间进行，但不能除以 (13) 中已为零的因子。

Gaussian dilation 的四点强制 Wick 配对。真实四点提供独立可算的检验量
\[
 \mathcal D_4(A,B,C,D)=\phi_\beta(ABCD)
 -\phi_\beta(AB)\phi_\beta(CD)
 -\phi_\beta(AC)\phi_\beta(BD)
 -\phi_\beta(AD)\phi_\beta(BC),                          \tag{15}
\]
对这里的中心化 fields 使用。旧归一 J_2 在 β=0 已有 D₄=q⁴−1≠0。故真实 thermofield 实现 (8) 与二点 Gaussian dilation (14) 是不同对象。由 (3) 的局部解析性，这个非零余项也在 β=0 的某个邻域持续；不声称任意 β、任意波包都非零。

## 6. Tracial state 与有限温度 Rindler 字典的精确障碍

设一个候选 *-map 同时保持态并 intertwine 微观 H 时间与某个非零倍数的 bulk boost 时间。tracial τ 的拉回仍 tracial。对其任意自伴动态 channel，正谱满足 S(−ω)=S(ω)。若同一态还对同一 flow 满足有限 β_R>0 KMS，则详细平衡要求
\[
 (1-e^{-\beta_R\omega})dS(\omega)=0,
\]
从而非零频率谱全部消失。因此已有非平凡动态谱不能通过这样的 state-preserving、time-intertwining *-map 被解释为 finite-β_R Wightman 态。非零常数时间缩放不解除此障碍。

这限定一个具体字典类型；不排除改变保留的代数、state-adapted observables、不同制备、不同 bulk/micro 时间关系或非均匀极限。特别是 (1) 的 Gibbs 倾斜已经改变微观态，不能用“选择 β=2π”代替 bulk boost 生成元和字典的推导。旧 midpoint 构造和插入数依赖热夹心问题见本稿开头引用。

## 7. 可继续的建设性接口

先用 (10)--(12) 在实际 common-J bilocal 家族上建立频带波包及全四点 Gram；它保留完整非 Gaussian 信息，能成为后续 state-adapted 算符的检验基准。辅助 (14) 提供相同二点的 canonical 对照，(15) 测量需要由 dressing/非线性结构承担的差异。

[Xu, arXiv:2403.09021v7](https://arxiv.org/html/2403.09021v7) §3.3 明确使用 Q-twisted Fock ladder realization，并保留 empty-chord state 的 tracial 性；不能把其中 ladder 的 vacuum polarization 与已有 tracial field 代数混同。[Rajgadia–Xu, arXiv:2604.14387v3](https://arxiv.org/html/2604.14387v3) 则通过 KM-state-adapted dressed operators 改变保留的极限可观测量，得到不同的态/代数描述。这提供下一项明确任务：为本项目的具体 complex-CAR/bilocal family 给出对应微观 dressed sequence，先验证其与原场的 mixed 二点和四点，再检查频带极化；不能仅因为两篇论文存在就宣布本项目已获得同一 dressing。

本稿只新增 (1) 的固定 β 共同矩传递、明确的谱权/平滑极化、真实 thermofield 与辅助 Gaussian dilation 的区分，以及指定 *-字典的 tracial 障碍。没有改变原 β=0 模型结论，也没有完成 de Sitter 体生成元或 canonical dressing。
