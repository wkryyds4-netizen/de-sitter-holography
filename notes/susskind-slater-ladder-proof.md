# Slater 极化 ladder 的全部固定 word 极限

2026-11-07，本地证明。N 偶，O 为 N/2 个 occupied sites，E=O^c；Ω_O 是相应 Slater determinant。p→∞、λ_N=p²/N→λ>0，q=e^{-4λ}∈(0,1)。以下首先取 disorder-annealed **Slater vector state**，不是原 full-Fock tracial state。最后另证固定 scalar words 的 disorder 集中。

在原 Hermitian complex Gaussian coupling 中，取 A⊂E、B⊂O、|A|=|B|=p 的一侧 block。其标准化系数 g_AB 是相互独立的 circular complex Gaussian，E g_AB=E g_AB²=0、E|g_AB|²=1；共轭 block 的系数为其复共轭。记
\[
 X_{AB}=C_A^\dagger C_B,\quad
 D_N=\binom{|E|}p\binom{|O|}p,\quad
 V^\dagger=\sqrt{w_N}\sum_{A,B}g_{AB}X_{AB},
\]
\[
 \sigma_N^2=w_ND_N,\qquad
 a_N^\dagger=V^\dagger/\sigma_N
       =D_N^{-1/2}\sum_{A,B}g_{AB}X_{AB},\qquad
 a_N=(a_N^\dagger)^\dagger.                              \tag{1}
\]
σ_N² 精确等于 E||V†Ω_O||²。耦合 g 的方差约定必须保留，不能混入实 Gaussian diagonal 的因子。

## 1. 相对 Slater 真空的 CAR 与符号

对 i∈E 令 f_i†=c_i†，对 i∈O 令 f_i†=c_i；这些 f 满足标准 CAR，且 f_iΩ_O=0。每个 X_AB 是 2p 个 f† 的纯 creation monomial，至多带一个由既定 C_A/C_B 次序确定的相位。因 2p 为偶数，若两个 labels 在 E、O 两侧分别不交，则其 monomials 以及互为 adjoint 的 monomials相互交换，交换符号为 +1。若两个 creation labels 有交集，creation 乘积为零。

保持任意原始 monomial 相位；一次 creation 与其 matching adjoint 的相位乘积为一。以下不通过丢弃 CAR 负号作 majorization，而是先证明 surviving forward words 的相位精确相消。

## 2. Gaussian pairing 中的 forward 与 backward 弦

对一个 operator word，按其作用于 Ω_O 的顺序，即从书面右端向左端，编号 t=1,...,L。creation 标 +，annihilation 标 −。令 h_t 为第 t 步之后的 creation-minus-annihilation 数，h_0=0。

若 L 不为 2m、总正负数不等，或某一 prefix 的 h_t<0，则 vacuum matrix element 精确为零。这可直接用 E-sites occupation：每次 a† 增加 p、每次 a 减少 p。

对 balanced Dyck path，complex Gaussian Wick 公式把期望展开为 m! 个正负端点双射。每条 Gaussian 弦赋予独立均匀 label (A,B)；不同弦的 labels 允许偶然相同，不能事先排除。每张图的归一系数恰为这些独立 labels 上的平均，且每个 raw CAR vacuum matrix element 的绝对值≤1。

称弦 forward，若其 creation 在 annihilation 之前；否则称 backward。**Backward pairing 并非有限 N 逐项为零。** 例如顺序 +−+− 的逆向 Gaussian 配对可以通过两条弦取同一 label 而非零，贡献量级 D_N^{-1}。

定义（大 N 时 rp<min(|E|,|O|)）
\[
 \rho_N(r)=\frac{\binom{rp}p}{\binom{|E|}p}
                  \frac{\binom{rp}p}{\binom{|O|}p}.      \tag{2}
\]
若 rp 超过某边的 site 数，将相应分子上截到该数，以保持≤1。

给定含 backward 弦的图，取第一个尚未打开其 matching creation 的 annihilation。此前非零路径只用已打开的 forward labels；目标 label 与它们独立。要使该 annihilation 非零，其 A 必须全含于此前已产生的≤mp 个 E particles，其 B 必须全含于此前已产生的≤mp 个 O holes。因此该图绝对贡献≤ρ_N(m)。对固定 m，
\[
 \rho_N(m)\le\left(\frac{K_m p}{N}\right)^{2p}
       \longrightarrow0.                               \tag{3}
\]
所有 backward 图总误差≤m!ρ_N(m)。这一步也覆盖极罕见的重复 Gaussian labels。

## 3. Forward 图的精确 finite-N 指示函数

固定一个全 forward pairing π。每条弦从 creation 到 annihilation 之间称为 active；两条弦的 active intervals 重叠，等价于它们 nested 或 crossed，disjoint intervals 不重叠。

**逐路径恒等式：** 对任意 label assignment，该图的 vacuum matrix element 为
\[
 \boxed{\prod_{\{u,v\}\in\mathcal E_\pi}
  1_{A_u\cap A_v=\varnothing}\,
  1_{B_u\cap B_v=\varnothing},}                          \tag{4}
\]
其中 \(\mathcal E_\pi\) 是 active-interval overlap graph 的边集。

证明用时间步归纳。此前路径非零时，状态是所有 active X_labels 的乘积作用于 Ω_O，且 active supports 两两不交。新 creation 若与任一 active label 相交，重复 f† 使该步为零；否则它加入这个乘积。一次 forward annihilation 的 matching label 仍 active；其 even monomial 可无负号地越过其他 disjoint active factors，随后 X†X 在该 support 上回到真空，精确删除该 factor。最后全部 labels 关闭，振幅为 +1。反之第一次违反不交条件的 creation 已使整个 word 为零。

已经关闭的 interval 与后续 interval 可以使用重叠甚至相同的 labels；(4) 不禁止它们。每个 monomial 的相位仍逐对相消。因此这里没有额外 q^{crossings} 的费米子交换因子：nested 和 crossed active overlaps 都受同一个 Pauli exclusion 条件。

## 4. Pair-Poisson 极限、triple overlaps 和重复 labels

也可不经Poisson直接给finite-N forward总式。对任一合法历史，h个active labels在E/O两侧的union大小恰为hp；新Gaussian label条件于过去仍独立均匀，故成功率精确为R_N(h)=binom(|E|−hp,p)binom(|O|−hp,p)/D_N。closing删除自己的active label，不增加概率因子；每次closing有h种matching选择。因此整个给定Dyck词的forward贡献精确为Π_close h Π_open R_N(h)。R_N(h)→q^h即给下文(8)。这是有限条件计数的直接证明，与下面pair-Poisson证明一致。所有归一式假定1≤p≤min(|E|,|O|)、J_s>0；固定k的规范化creation向量在kp≤min(|E|,|O|)时定义，DS极限中最终自动满足。

固定 m 个独立 p-subsets 在 E 上，全部 pair intersections
\(Y_{uv}^{E}=|A_u\cap A_v|\) 联合趋独立 Poisson(2λ)；O 上有独立的同一极限。可由联合 factorial moments 证明：为 K 次指定 pair intersections 选互异 sites 的主项为
\[
 (|E|)_K\prod_u\frac{(p)_{d_u}}{(|E|)_{d_u}}
       \to(2\lambda)^K,\qquad\sum_ud_u=2K.
\]
不同 pair 要求落在同一个 site 的修正至少含三标签。其总发生概率至多
\[
 2\binom m3\frac{p^3}{(N/2)^2}\to0.                    \tag{5}
\]
固定高阶 factorial moments 中的合并布局同样趋零。重复整块 labels 的概率≤binom(m,2)/D_N→0；可留在精确平均中，不影响联合分布极限。

(4) 是有界 0/1 函数，所以这里无需指数 overlap majorant，也不产生 degree² 指数损失。联合 Poisson 给
\[
 \mathbb E(4)\longrightarrow
 e^{-4\lambda|\mathcal E_\pi|}=q^{|\mathcal E_\pi|}.       \tag{6}
\]
本节只取每个 word 次数固定后 N→∞，不声明随 N 增长的 word 长度估计。

## 5. 所有 word 的闭式矩与 weighted-shift 识别

对给定 Dyck height path，每次 creation 之前有 h 个 active labels，新 interval 与它们各新增一条 overlap edge。因此
\[
 |\mathcal E_\pi|=\sum_{\mathrm{creation}}h_{\rm before}
       =\sum_{\mathrm{annihilation}}(h_{\rm before}-1),   \tag{7}
\]
与具体 forward matching 无关。每次 annihilation 可选择当前任一个 active creation，故 forward pairing 数为
\(\prod_{\mathrm{annihilation}}h_{\rm before}\)。结合 (3)、(6) 得全部 word 定理：
\[
 \boxed{\lim_NE\langle\Omega_O,w(a_N,a_N^\dagger)\Omega_O\rangle
 =\prod_{\mathrm{annihilation}}
       h_{\rm before}\,q^{h_{\rm before}-1}}             \tag{8}
\]
对 balanced Dyck path 成立，其他 word 为零。

在 \(\ell^2(\mathbb N_0)\) 上令 n 为 number operator、b 为普通 bosonic lowering，定义
\[
 a^\dagger=b^\dagger q^{n/2},\qquad a=q^{n/2}b.
\]
于是
\[
 a^\dagger|k\rangle=\sqrt{k+1}\,q^{k/2}|k+1\rangle,
 \qquad a|k\rangle=\sqrt k\,q^{(k-1)/2}|k-1\rangle.      \tag{9}
\]
逐 height step 相乘恰好给 (8)，所以这是完整 joint *-moment 极限，超出了纯 creation norms。q 固定小于一时 a 为有界 compact weighted shift，\(\|a\|^2=\max_{k\ge0}(k+1)q^k\)。它的交换关系为
\[
 aa^\dagger-q a^\dagger a=q^n,
\quad[a,a^\dagger]=(n+1)q^n-nq^{n-1},                   \tag{10}
\]
并非原始微观 CCR。普通 b 只在有限 occupation domain 上由 \(b=q^{-n/2}a\) 恢复；逆 dressing 在全空间上无界。

纯 creation calibration 为精确有限 N 式
\[
 \mu_{k,N}=E\|(a_N^\dagger)^k\Omega_O\|^2
 =k!\prod_{j=0}^{k-1}
 \frac{\binom{|E|-jp}p\binom{|O|-jp}p}
      {\binom{|E|}p\binom{|O|}p}
 \longrightarrow k!q^{k(k-1)/2}.                        \tag{11}
\]
此外，完整 word 定理真正补上 annihilation 的闭合：对每个固定 k≥1，
\[
 E\|a_N(a_N^\dagger)^k\Omega_O
       -kq^{k-1}(a_N^\dagger)^{k-1}\Omega_O\|^2\to0.     \tag{12}
\]
证明是展开差的三个 word，使用 (8)；极限交叉项为 μ_k，范数比 μ_k/μ_{k-1}=kq^{k-1}。有限 N 可以有 leakage，(12) 才说明它在固定 occupation 的 DS 极限消失。不同 k 的 creation vectors 在有限 N 已因 E-particle number=kp 严格正交。

### 5.1 实际有限 N 的 inverse dressing 与固定 core CCR

令 \(n_N=Q_E/p\)，其中 \(Q_E=\sum_{i\in E}c_i^\dagger c_i\)，并取 \(q_N=e^{-4\lambda_N}\)。有限 N 精确有 \([n_N,a_N^\dagger]=a_N^\dagger\)。定义实际算符
\[
 b_N=q_N^{-n_N/2}a_N,\qquad
 b_N^\dagger=a_N^\dagger q_N^{-n_N/2}.                   \tag{12a}
\]
它们精确互为伴随。在任何从 Ω_O 开始的 word 中，所有中间向量的 n_N 值都是相应 height h；即使 a_N 产生有限 N leakage，occupation number 仍相同。因此 creation 步的 dressing 因子为 q_N^{-h/2}，annihilation 步为 q_N^{-(h-1)/2}。对闭合 Dyck word，(7) 给整个矩阵元的精确关系
\[
 \langle\Omega_O,w(b_N,b_N^\dagger)\Omega_O\rangle
 =q_N^{-\sum_{\rm annih}(h-1)}
       \langle\Omega_O,w(a_N,a_N^\dagger)\Omega_O\rangle.
                                                               \tag{12b}
\]
故其所有固定 words 的极限是普通 bosonic vacuum law
\[
 \lim_NE\langle\Omega_O,w(b_N,b_N^\dagger)\Omega_O\rangle
           =\prod_{\rm annih}h,                         \tag{12c}
\]
非 Dyck words 仍为零。特别对任意固定 noncommutative polynomial P，展开平方范数后使用 (12c) 得
\[
 E\|([b_N,b_N^\dagger]-1)P(b_N,b_N^\dagger)\Omega_O\|^2
                         \longrightarrow0.             \tag{12d}
\]
这给真正的实际微观 sequence、指定新态与固定 height core 上的 asymptotic CCR。它不声称 finite-N 全矩阵算符满足 CCR；q_N^{-n_N/2} 在完整 Fock 空间上随 N 不一致有界，不能把该结论延拓到原 full-trace 代数或增长 occupation。下面 scalar 集中界在乘以固定 word 的 dressing 常数后同样适用。

## 6. Scalar disorder self-averaging：并非假设，亦可证明

annealed (8) 本身不需要样本自平均。若要 fixed typical-disorder scalar words，也有直接界。记
\(F_{w,N}(g)=\langle\Omega_O,w(a_N,a_N^\dagger)\Omega_O\rangle\)，word 长度 2m。

展开 E|F_w|²−|EF_w|² 的两个 vacuum replicas。所有 Gaussian pairing 完全位于各自 replica 内的项精确抵消。每个剩余图至少有一个跨 replica 的 Gaussian label α；它在每个 replica 内只出现一次。要使相应 vacuum matrix element 非零，α 的每个 E site 及每个 O site 必须也出现在至少一个其他 Gaussian label 中，否则该 site 在此 replica 内只被 creation 或 annihilation 一次，无法返回真空。

其他最多 2m−1 个 labels 各占 p sites。条件于它们，α 的 A/B 独立均匀，故该图的绝对贡献≤ρ_N(2m−1)。因此
\[
 \boxed{\operatorname{Var}F_{w,N}
              \le(2m)!\,\rho_N(2m-1)\longrightarrow0.} \tag{13}
\]
不使用 matrix-operator-norm 集中，也不假设不同 N 的耦合独立。沿严格增长的整数 N 双标度序列，(3) 的界可求和；若这些随机矩阵定义在共同概率空间上，Chebyshev 与第一 Borel--Cantelli 引理给每个 word 几乎处处收敛，再对可数多个 words 取交。没有声明增长 degree、uniform-in-word 或 H 动力学的样本极限。

## 7. 这个构造解决什么、尚未解决什么

由 (9) 可直接看见该 ladder 的极限 von Neumann algebra 是 B(ℓ²)：a 的 kernel 为真空一维，故 spectral projection P_0=1_{\{0\}}(a†a) 在 W*(a) 中，而 \((a†)^iP_0a^j\) 经非零常数归一给全部 matrix units。真空 vector state 为纯态。这个结论只属于新选择的 Slater-adapted ladder sector，不能移回原 tracial bilocal algebra。

(12) 证明 ladder 的自身 annihilation 闭合；它不证明原 full Hamiltonian H 保持这个 cyclic sector，也不证明该 weighted shift 就是 bulk boost generator。full H 的首步激发扇区及剩余权重由 [独立文献/微观接口](susskind-slater-ladder-literature-interface.md) 处理，应另行保留其 H-leakage 结论。下一步是这个新 ladder 与 H、原 singlet/bilocal 的 mixed words，而不是把 a+a† 直接替换为 H。

差别在有限 N 已可量化：每个 X_AB†X_AB 是“B 全 occupied、A 全 empty”的 projector，所以
\[
 E\operatorname{tr}_{\rm full}(a_N^\dagger a_N)=4^{-p},
 \qquad E\|a_N^\dagger\Omega_O\|^2=1.                  \tag{14}
\]
因此这个具体 sequence 在旧 annealed full-trace L² 中消失，在新的 Slater state 中却保留非零归一矩。它并不是把旧 tracial 字段经态保持 *-map 变成 canonical oscillator；它明确改变了极限保留的 observable/state interface。

本稿没有运行 scientific main 或修改任何封存输出；有限 CAR 检查应独立保存，验证 algebra 而不替代上述任意固定 word 证明。
