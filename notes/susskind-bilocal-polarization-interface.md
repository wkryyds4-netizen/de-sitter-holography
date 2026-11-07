# 从 tracial 集体场到 state-adapted bilocal

2026-11-07。沿同 J overlapping complex p-body CAR 模型，固定 lambda_N=p²/N→lambda>0，
q=e^(−4lambda)、C=J_s²/lambda。旧 h/u/x 与完整 bilocal 的态仍是 beta=mu=0
normalized annealed full-Fock trace。本稿单独引入参考 Slater 态时会显式标记；
不把这个新态的结果回写成原 trace 的 CCR。原文位置见
[canonical sources](../literature/susskind-bilocal-canonical-sources.md)。

## 1. 正规序中心项实际上记录什么

为避开不同论文的 Fourier/质量规范，先直接做有限 CAR。令 a,b 为一粒子模式，
r=1,...,N_c 为独立颜色，E_ab=sum_r c_(a,r)† c_(b,r)。对任意参考一粒子矩阵
gamma_ab=omega(c_a†c_b)，定义

\[
 M_{ab}=N_c^{-1}E_{ab}-\gamma_{ab}I.
\]

CAR 精确给

\[
 [M_{ab},M_{cd}]=\frac1{N_c}
 \{\delta_{bc}M_{ad}-\delta_{ad}M_{cb}
   +\delta_{bc}\gamma_{ad}-\delta_{ad}\gamma_{cb}\}.\tag{1}
\]

gamma 对角、占据为 n_a 时，最后一项是
N_c^(−1) delta_bc delta_ad(n_a−n_b)。它依赖占据差，并不是凭选取复系数
就出现的 universal scalar。若改用另一个参考 P 正规序、却保留原态，
M 的非零均值会补偿这个变化；例如 gamma=I/2 的 trace 下，人工选择 sea P
产生的中心项，被 mean(M)=I/2−P 的项精确抵消。

对 Slater sea P²=P，定义 epsilon=I−2P。纯一粒子 projector 的经典轨道给
Phi²=I；相对于 epsilon 的微小变化满足 {epsilon,deltaPhi}=0。
因而 tangent 是 filled/empty 两块之间的 particle–hole 块。
这说明 Itakura/Rajeev 最可迁移的建构次序：先定 occupation/最低权态，再定
非对角块和非线性 dressing。它不是先用 many-body adH 的正频率自行定义海。
本项目的 N 个模式没有现成的额外 N_c 颜色极限。

## 2. large-color boson contraction 的具体态条件

取 occupied 模式 i 和 empty 模式 a，令

\[
 b_{ai}^{(N_c)}=N_c^{-1/2}\sum_r c_{i,r}^\dagger c_{a,r}.
\]

在每色相同、颜色间乘积的对角 quasi-free 态中，精确二点为

\[
 \omega(bb^\dagger)=n_i(1-n_a),\qquad
 \omega(b^\dagger b)=n_a(1-n_i),\qquad
 \omega([b,b^\dagger])=n_i-n_a.\tag{2}
\]

当 kappa=n_i−n_a>0，除以 sqrt(kappa) 后才有单位中心项。
极限 occupation 是 n_a(1−n_i)/kappa；选定通道是 vacuum 模式的条件是
n_a(1−n_i)=0，而 filled/empty sea n_i=1,n_a=0 是自然的充分选择。
对于有限个模式、固定 word degree，颜色独立性使高于二阶的 fluctuation
cumulants 为 O(N_c^(1−k/2))，只留下有序二点配对。
所以这个标准 fluctuation 极限给具有 (2) 的 quasi-free CCR；这里是可直接
由颜色分组证明的模型构造，不是对原 DSSYK 增添一个已存在的颜色参数。

在 filled/empty sea 上 b annihilates reference vacuum；固定 K 个 particle–hole
激发子空间内，(1) 的 density correction 为 O(K/N_c)。这才对应 HP 展开
sqrt(N_c−A)/sqrt(N_c)→I 的低占据范围。在 n_a=n_i=1/2 的 trace 中，
中心项为零，颜色 fluctuation 极限则为 commuting Gaussian fluctuations。

## 3. 旧 tracial 代数内的直接障碍

在 normalized trace tau 中，任何 a 有 aa†、a†a∈L¹ 时，

\[
 \tau(aa^\dagger)=\tau(a^\dagger a),\qquad
 \|[a,a^\dagger]-\kappa I\|_r\ge|\kappa|,\quad r\ge1.\tag{3}
\]

第二式只用 trace 和 Hölder。这也适用于 N-dependent 线性组合、时间滤波、
以及有相应矩控制的极限；不允许在保持该态和矩收敛时得到非零 vacuum CCR。
GNS 中 ||aOmega||²=tau(a†a)=||a†Omega||²，所以非零 a 不能只湮灭 Omega。
用 adH 的正频谱投影取 a，仍是同一个 tracial 算符空间中的元素，不能改变 (3)。

旧 h/u/x 的 mixed-q 表示可以有辅助 creation/annihilation；但这些 ladder
不是原 tracial field algebra 的元素。它们湮灭空向量，而 field algebra 上
该向量 separating。这是代数范围的区别，不能仅靠把“正频”改名解决。
本轮另有低扇区 Q-Fock 检查；此处不重复其矩计算或旧 Rindler 四点结果。

## 4. 换成 Slater 纯态，为什么还需要改变 observables

令 S⊂{1,...,N}，|S|=m；|S> 是该 occupied pattern 的 Slater determinant。
对于依赖同一随机 J 的 neutral 算符 W(J)，假设其在物理 U(N) 旋转下协变：
W(U·J)=U_F W(J)U_F†。GUE 耦合分布不变，所以 E_J W_m(J) 与 wedge^m U
对易。该表示不可约，由 Schur 得到精确有限 N 恒等式

\[
 \boxed{\mathbb E_J\langle S|W(J)|S\rangle
       =\mathbb E_J\operatorname{tr}_m W(J).}\tag{4}
\]

所有旧 h/u/x 的 ordered words、协变的时间演化及完整 bilocal 都满足这个条件。
右侧是 fixed-charge trace；不擅自把它等同于有限 N 的 full-Fock charge 混合。
由于乘积仍协变，限制在这些旧 observables 上，左边已经是 tracial functional。
所以仅更换微观纯态的名字并不够。需要显式依赖 S、从而不属于这个固定-S
协变类的 observables。

同理，unnormalized Euclidean-prepared numerator
E<S|e^(−beta H/2)W e^(−beta H/2)|S> 等于 sector thermal numerator。
对每个样本先除以随机 norm 再取平均，则涉及 quenched ratio；(4) 不直接证明
其等于 annealed ratio，不在这里跳过该区别。

## 5. 在原复杂费米子模型中定义一个适配 S 的同 J 算符

设 S 为 filled、S^c 为 empty，定义纯 particle–hole creation 部分

\[
 V_S=\sqrt{w_N}\!
 \sum_{\substack{A\subset S^c,\ B\subset S\\|A|=|B|=p}}
 J_{AB}C_A^\dagger C_B,\qquad
 \sigma_S^2=w_N\binom{N-m}p\binom mp,\qquad
 M_S=V_S^\dagger/\sigma_S.\tag{5}
\]

这里使用原 Hamiltonian 的同一 J 子块，不引入独立 probe。
M_S 净物理电荷为零，且 M_S|S>=0 精确；annealed 后
omega_S(M_SM_S†)=1、omega_S(M_S†M_S)=0。
它突破 (3) 的原因是态和 observable family 都已改变。
事实上在旧 normalized full-Fock trace 中，标准 CAR 还给精确式

\[
 \mathbb E\tau_{\rm Fock}(M_S^\dagger M_S)
 =\mathbb E\tau_{\rm Fock}(M_SM_S^\dagger)=4^{-p}.\tag{5a}
\]

每个 disjoint A/B monomial 的迹范数平方是4^(−p)，Gaussian covariance 去掉交叉项。
所以该 family 在旧 tracial L² 中消失，却在 Slater 态中有非零真空二点；
不能将两个极限的 observable classes 不加说明地认同。

固定 k 时，可直接数互不交叠的 particle/holes 集合，得

\[
 \mathbb E\omega_S[M_S^k(M_S^\dagger)^k]
 =k!\prod_{j=0}^{k-1}
 \frac{\binom{N-m-jp}p\binom{m-jp}p}
      {\binom{N-m}p\binom mp}.\tag{6}
\]

在 m/N→1/2 下，右边趋 k! q^(k(k−1)/2)，例如四阶为2q而非canonical 2。
因此“存在真正湮灭参考态的算符”仍不等于“该算符本身已是 canonical boson”。

本轮 [完整单 ladder 证明](susskind-slater-ladder-proof.md) 进一步处理所有固定 words：
forward Gaussian pair 的 active intervals 只要共存，不论 nested 或 crossed，
两侧底层 p-set 都必须不交叠。合法项的 even-CAR 符号为正，pair-Poisson 给
q^(crossings+nestings)。在每次 closing height=h 处，配对选择与权重合并为
h q^(h−1)。这给单模 weighted shift

\[
 M|n\rangle=\sqrt n\,q^{(n-1)/2}|n-1\rangle,
 \qquad M=q^{\hat n/2}b,\quad[b,b^\dagger]=I.\tag{7}
\]

有一个重要的有限 N 余项：Gaussian matching 中 annihilation 先于其配对
creation 的项不一定逐项为零。取第一个未打开 label 的 annihilation，
它的 p-set 必须落入此前至多 kp 个已激发 sites；两侧独立集合给上界
binom(kp,p)^2/[binom(N−m,p)binom(m,p)]，固定 k 时消失。
不能将这些项在 finite N 直接删除。

## 6. Canonicalization 的定义域与实际 Hamiltonian 的第一道检验

(7) 的 weights 非零且趋零，M 是 bounded weighted shift。其 kernel 是真空一维，
故 P_0=1_{0}(M†M) 属于 W*(M)；M†^n P_0 M^k 给全部 matrix units。
于是单 ladder 的 vacuum GNS 表示满足 W*(M)=B(ell²(N_0))。
在有限粒子 core 可以定义

\[
 b=q^{-\hat n/2}M,\qquad b^\dagger=M^\dagger q^{-\hat n/2}.\tag{8}
\]

inverse dressing 无界；这是新代数中的 affiliated canonical ladder，
不是原 tracial h/u/x 中的 bounded linear combination。

还可把 (8) 写成实际有限 N 的 state-adapted family，而不只给抽象存在性。
令 n_(S,N)=Q_(S^c)/p；在从 |S> 出发的 M/M† words 中，它每次精确变化±1。
定义 b_(S,N)=q^(−n_(S,N)/2)M_S，其伴随为 M_S†q^(−n_(S,N)/2)。
固定 word 只经过有限高度；每次 closing height=h 的两个 dressing 合计
q^(−(h−1))，精确抵消极限 M word 中的 q^(h−1)。于是这些实际 b_(S,N)
在新 Slater 态的全部固定 words 趋 canonical vacuum moments。
这是固定激发 core 的 moment statement，不是全空间算符范数或旧 trace 的收敛；
全空间 inverse dressing 随 N 增长，不能据此保持原 uniform tracial estimates。

同时，单 ladder 并不足以容纳实际 H。作用 H|S> 后，若 k 个粒子从 S 到 S^c，
其 annealed 范数权重是

\[
 W_k=w_N\binom mp\binom pk\binom{N-m}k,
 \qquad\sum_{k=0}^pW_k=C_m.
\tag{9}
\]

令 r=p−k。归一分布为

\[
 \Pr(r)=\frac{\binom pr\binom{N-m}{p-r}}{\binom{N-m+p}p},
 \quad\mathbb E(r)_t=\frac{(p)_t^2}{(N-m+p)_t}
 \longrightarrow(2\lambda)^t.\tag{10}
\]

所以 r→Poisson(2lambda)，sigma_S²/C_m→exp(−2lambda)=sqrt(q)。
在 annealed GNS 里，K_S=closure span{(M_S†)^k Omega_S} 的不同 k 按 kp
个激发正交；omega_S(H) 的 disorder 平均为零，精确有

\[
 P_{K_S}H\Omega_S=\sigma_S M_S^\dagger\Omega_S,\qquad
 \|(1-P_{K_S})H\Omega_S\|^2=C_m-\sigma_S^2
 \longrightarrow C(1-\sqrt q)>0.\tag{11}
\]

这给最小的 same-J closure 检验，已在二点层面识别缺失方向。它排除把整个 H
直接认成单 ladder 的 sigma(M+M†)，但不排除更大的 joint(H,M) 适配代数。

## 7. 2026 KM 构造实际帮助到哪里

Rajgadia–Xu 的有用之处是提供相同的操作次序：保留 state-adapted strings，
再识别可见 dressed operator 与 bare canonical ladder 的区别。
本轮 (5)–(11) 是这个思路在 complex charge-conserving、原始同 J 模型内的具体适配。
它的单 M sector 已可独立证明；完整 mixed(H,M) 仍要从本模型计算。
原 Majorana 模型的 spin-charge、长度规范、probe coefficient family 与这里不同，
所以不能直接搬入其混合 chord kernel 或 Type-I 全代数结论。

下一项最有效的工作是保留 (11) 的额外 H 方向，计算其与 M、M† 的首批 ordered
mixed words及渐近 pairing 规则，再寻找闭合的多方向表示。
Itakura/FLZ 的纵向一粒子变量及质量算符仍是后续独立字典；这里的 excitation
number与many-body能差均没有被命名为介子质量。
