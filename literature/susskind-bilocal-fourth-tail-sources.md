# 完整 bilocal 四阶尾：矩阵 Gaussian 与 Clifford 工具的适用接口

2026-11-07，本地文献核查与推导备忘。接续
[低通道实际接口](../notes/susskind-bilocal-four-point-transfer.md)、
[普通传播子 rank 界](../notes/susskind-memory-wick-rank.md)、
[中性 marked series](../notes/susskind-memory-neutral-marked-series.md)。
本文件不改动已经封存的两点/低通道四点结果，也不宣称完整 bilocal 四点已完成。

模型仍为 standard CAR、允许重叠的 complex Gaussian p-body Hamiltonian、同一个原始 J、
归一化 annealed beta=mu=0 full-Fock trace，偶 N、p²/N=lambda_N→lambda>0，固定 Js，
hbar=1/string 时间。把独立复耦合拆成标准实 Gaussian 坐标后写
H(g)=sum_a g_a B_a，B_a 为确定性 Hermitian 矩阵，耦合方差已吸收入 B_a。
以下 j 是 coupling-chaos 次数，不能与 p 或 Clifford length 混同。

## 1. 这轮真正需要的范数

令 d=2^N、tau_d=Tr/d，定义

\[
 \|X\|_{r,\mathrm{ann}}
 =\bigl(\mathbb E\tau_d|X|^r\bigr)^{1/r}.
\]

这是 L^r(Gaussian;S_r(tau_d))，不是 L^r(Gaussian;S_2)。本轮足够目标是

\[
 \lim_{K\to\infty}\limsup_{N\to\infty}
 \sup_{|t|\le S}\|
 \widetilde{\mathcal B}_N(t)-\widetilde{\mathcal B}^{(K)}_N(t)
 \|_{4,\mathrm{ann}}=0.
\tag{1}
\]

现有 actual−filtered 的 L4 返回误差 O(p^(-1/7)) 不提供 (1)。
普通 Taylor 系列的一种足够输入是，对实际完整系数 F_{n,N} 有

\[
 \|F_{n,N}\|_{4,\mathrm{ann}}
 \le A(n+1)^b L^n\sqrt{n!},
\tag{2}
\]

常数在所需 N/charge 区域一致；再接原来的全电荷尾。
因为 sum_n (LS)^n(n+1)^b/sqrt(n!) 在每个有限 S 收敛。
对 Hermite/marked series 则使用该系列已有的系数界，不能把两种系数混换。
只有“每个固定 n 的 limsup_N 满足 (2)”还不能交换 limsup_N 与无穷求和：
须有统一尾，或增长截断加独立高次数余项估计。

## 2. 已核主来源与定位

1. De Huang–Joel A. Tropp，*From Poincaré inequalities to nonlinear matrix concentration*，
   Bernoulli **27** (2021), 1724–1744，
   [作者原文](https://tropp.caltech.edu/papers/HT21-Poincare-Inequalities.pdf)，
   [arXiv:2006.16561](https://arxiv.org/abs/2006.16561)。
   §2.4 Example 2.3 给 Gaussian carré-du-champ；§2.6 Assumption 2.6 给条件；
   Theorem 7.1/(7.1)，印刷页 1739，是本轮可用的 Schatten 矩界。
   Example 7.3 明确指出二次 matrix chaos 的进一步化简受到 partial-transpose 障碍。

2. Madhur Tulsiani–June Wu，*Simple Norm Bounds for Polynomial Random Matrices via Decoupling*，
   ITCS 2025, 91:1–91:21，
   [正式全文](https://drops.dagstuhl.de/storage/00lipics/lipics-vol325-itcs2025/html/LIPIcs.ITCS.2025.91/LIPIcs.ITCS.2025.91.html)。
   §5 Lemmas 12–13 把前项写成 Gaussian 导数和 row/column 形式；Claim 15 给导数块的递归。
   本笔记只取一步不等式及递归结构，下面另写显式四阶迭代，不依赖隐藏在 degree 记号里的常数。

3. Joel A. Tropp，*User-Friendly Tail Bounds for Sums of Random Matrices*，
   [作者原文](https://tropp.caltech.edu/papers/Tro11-User-Friendly-preprint.pdf)，
   §4.7/(4.17)，PDF 第 19 页，列出 Gaussian 非交换 Khintchine 矩界，
   最优常数归于 Buchholz, Theorem 5。这里只用线性系列；不把同一原 J 的乘积当独立系列。

4. A. S. Bandeira–K. Lucca–P. Nizić-Nikolac–R. van Handel，
   *Matrix Chaos Inequalities and Chaos of Combinatorial Type*，
   [arXiv:2412.18468v2](https://arxiv.org/html/2412.18468v2)。
   §2.2、§2.4.2 和 Appendix A.2 是有向 flattening/逐坐标迭代的主来源。
   Theorems 2.4–2.5 在所核版本是算符范数估计，带 log(d+m) 与依赖阶数的隐常数；
   不能直接作为本模型 fixed-S4 的全时间可和尾定理。

5. Eric A. Carlen–Jan Maas，*An analog of the 2-Wasserstein metric in non-commutative probability
   under which the Fermionic Fokker–Planck equation is gradient flow for the entropy*，
   [arXiv:1203.5377v1](https://arxiv.org/html/1203.5377v1)，§2.3/(19)–(24)、§2.4。
   这里使用 Clifford 梯度、grading、Leibniz 与 number operator 的规范。

6. Eric A. Carlen–Elliott H. Lieb，*Optimal Hypercontractivity for Fermi Fields and Related
   Non-Commutative Integration Inequalities*，
   [hep-th/9209053v1](https://arxiv.org/pdf/hep-th/9209053)，(1.3)–(1.6)、Theorem 4。
   它控制物理 Clifford semigroup；不是仅在 coupling Gaussian 坐标上作用、在矩阵上恒等的半群。

## 3. Gaussian Sobolev：一个能直接检查的充分条件

Huang–Tropp 的条件是可逆齐次 Markov 过程、trace Poincaré 常数 alpha，及使极限/导数合法的函数域。
对标准 Gaussian OU，alpha=1；矩阵多项式在域内。Hermitian F 的专门化为

\[
 \|F-\mathbb EF\|_{2s,\mathrm{ann}}
 \le \sqrt2\,s\,
 \left(\mathbb E\tau_d
 \left[\sum_a(\partial_aF)^2\right]^s\right)^{1/(2s)},
 \qquad s=1\ \text{或}\ s\ge3/2.
\tag{3}
\]

特别地，s=2 的四次幂常数是 64，不是一般的 L2→L4 常数。
以下是从 (3) 的 Hermitian dilation 得到的项目规范换算：对一般复矩阵 F，设

\[
 G_R=\sum_a(\partial_aF)(\partial_aF)^*,\qquad
 G_C=\sum_a(\partial_aF)^*(\partial_aF).
\]

\[
 \mathbb E\tau_d|F-\mathbb EF|^4
 \le32\,\mathbb E\tau_d(G_R^2+G_C^2).
\tag{4}
\]

证明只需把 F 放入 [[0,F],[F*,0]]；其四阶 trace 是 2Tr|F|^4，梯度平方是 diag(G_R,G_C)。
Tulsiani–Wu Lemma 13 写的是较松的常数；(4) 的 1/2 来自这个显式 dilation。
矩形/导数块同样用普通 Tr，最后统一除**原始** d，不另除增大的辅助矩阵维数。

因此有两条具体可查路径：

- 对每个 ordinary 系数 F_{n,N} 控制
  (E tau_d(G_R²+G_C²))^(1/4) 为 (2) 的增长；
- 或对完整有限时间原场 F_N(t) 及其截断余项控制该平方函数，直接得到 (1)。

两条都必须实际证明平方函数界。只知 E G_R 或 E F² 是标量，不够控制 E tau_d G_R²。

**先去均值，再 pinching。** 若未 pinching 的 F 保持 U(N) 等变且守电荷，Schur 引理使
EF 在每个电荷扇区是标量，所以 P_H(EF)=EF。逐样本 pinching 是 Schatten 收缩，故

\[
 \|(1-\mathcal P_H)F\|_{4,\mathrm{ann}}
 \le2\|F-\mathbb EF\|_{4,\mathrm{ann}}.
\tag{5}
\]

应把 (3)/(4) 用在未 pinching 的 F；无需对随 H 变化的精确谱投影求导。
全 Fock 的 EF 可以是 Q 的函数，不必是整个空间上的单一标量。

## 4. 线性 scalar covariance 能做什么；非线性为何另需证明

线性 Hermitian X=sum_a g_a A_a 的已知矩界是

\[
 \mathbb E\tau_d X^{2k}
 \le(2k-1)!!\,\tau_d\left(\sum_a A_a^2\right)^k.
\tag{6}
\]

若 sum A_a²=vI，右端为 (2k−1)!!v^k。这正是低生成元 B1 的合法输入。
对于多项式，下面的**自算反例**表明 U(d) 等变和 scalar covariance 本身不够：
z∈C^d 为标准复 Gaussian，Ezz*=I，取

\[
 X=\frac{zz^*-I}{\sqrt d}.
\]

X 是纯二阶 coupling-Wick chaos、Hermitian 且 U(d) 等变，并且

\[
 \mathbb EX=0,\qquad\mathbb EX^2=I,\qquad
 \mathbb E\tau_dX^4=d+2+\frac5d+\frac1{d^2}.
\tag{7}
\]

验证：r=||z||²~Gamma(d,1)，X 的特征值为 (r−1)/sqrt d 和 d−1 个 −1/sqrt d；
用 Er^k=d(d+1)…(d+k−1) 展开即可。
这不是本模型的反例，而是排除“任何协变固定 chaos 都能由二范数控制四范数”的论证。

更基本地，coupling OU tensor Id_matrix 在确定性 rank-one 投影 P 上不变，
||P||4/||P||2=d^(1/4)。所以把标量 Gaussian 超可积直接张量上矩阵恒等映射，
得不到归一化 annealed S2→S4 无维数定理。

## 5. 显式四阶 row/column 迭代：把待证问题落在确定性系数上

以下是基于线性矩界的短推导，不作为本模型已经满足的估计。
令 g^(1),…,g^(j) 为独立标准实 Gaussian 向量，A_{i_1…i_j} 为确定性 d×d 系数，

\[
 Y_j=\sum_{i_1,\ldots,i_j}A_{i_1\cdots i_j}
 g^{(1)}_{i_1}\cdots g^{(j)}_{i_j}.
\]

对 S⊂{1,…,j}，F_S(A) 把 S 内 coupling indices 和原矩阵行放在行端，
补集 coupling indices 和原矩阵列放在列端。定义

\[
 f_{4,S}(A)=\left(d^{-1}\operatorname{Tr}|F_S(A)|^4\right)^{1/4}.
\tag{8}
\]

辅助 coupling 指标使用 counting trace；只除原始 d。
对任意矩形线性系列 Y=sum_i g_i A_i，(6) 的 dilation 给

\[
 \mathbb E\operatorname{Tr}|Y|^4
 \le\frac32\left\{
 \operatorname{Tr}\left(\sum_i A_iA_i^*\right)^2+
 \operatorname{Tr}\left(\sum_i A_i^*A_i\right)^2\right\}.
\tag{9}
\]

每一步先条件于其余 Gaussian 副本，(9) 的两项分别把该指标变成行或列。
递归 j 次，保持所有辅助 counting traces，得到

\[
 \mathbb E\tau_d|Y_j|^4
 \le(3/2)^j\sum_{S\subset[j]}f_{4,S}(A)^4
 \le3^j\max_S f_{4,S}(A)^4.
\tag{10}
\]

这个常数对 j 显式，且未出现 log d。所有 2^j 个方向都要保留；
只计算 S=empty 或 S=[j] 的两种 covariance flattening，不足以使用 (10)。

### 从同一原 J 的 Wick chaos 到 (10)

对称系数 A 定义

\[
 X_j=\sum_I A_I:g_{i_1}\cdots g_{i_j}:.
\]

令 g=j^(-1/2)sum_r g^(r)。有限 Gaussian 条件期望恒等式为

\[
 \mathbb E[Y_j\mid g]=j^{-j/2}X_j.
\tag{11}
\]

可用条件协方差 Cov(g^(r),g^(s)|g)=(delta_rs−1/j)I 的生成函数逐项验证；
同坐标重复也保留 Wick contraction，因此不要求 square-free。
条件 Jensen 与 (10) 给

\[
 \|X_j\|_{4,\mathrm{ann}}
 \le j^{j/2}3^{j/4}\max_S f_{4,S}(A).
\tag{12}
\]

j^(j/2)≤e^(j/2)sqrt(j!)，所以若真正证明所有定向系数满足

\[
 \max_S f_{4,S}(A_{j,N})\le A_0(j+1)^bK^j,
\tag{13}
\]

则 (12) 给 (2) 型可和增长。这把困难转成**实际 marked/Wick 后代全部有向系数张量**的估计。
Gaussian 解耦和固定 S4 迭代本身只付指数与 sqrt(j!)；
尚未建立的是 (13)，尤其作用 rank/body 缺额之后的系数。
一般多项式应先按 coupling chaos 分层，再用三角不等式；不能只保留最高 chaos。

二范数等距 ||X_j||2²=j! sum_I tau_d(A_I* A_I) 只控制 Frobenius 量；
用未归一化 S4≤S2 生搬到 (8) 会留下 d^(1/4)，并不补足 (13)。

## 6. Clifford product rule 给出的精确 deficit 接口

用 2N 个 Majorana Q_i，{Q_i,Q_j}=2delta_ij，grading 为 gamma。
Carlen–Maas 的梯度是 nabla_i A=(Q_i A−gamma(A)Q_i)/2，
且 nabla_i(AB)=gamma(A)nabla_iB+(nabla_iA)B。
其 number operator N_Cl=sum nabla_i* nabla_i 的特征值是 Clifford length；
项目的 L=N_Cl/2。

下面是该规范下的有限代数推论。对 even A,B，

\[
 L(AB)=(LA)B+A(LB)+\sum_i(\nabla_iA)(\nabla_iB).
\tag{14}
\]

例如可从 L A=(1/4)sum_i(A−Q_i A Q_i) 展开；最后一项的正号来自 graded 导数。
对应正 carré-du-champ 为
Gamma_L(A,B)=(1/2)sum_i(nabla_iA)*nabla_iB。

coupling Gaussian number Ncal=sum_a(g_a partial_a−partial_a²) 满足

\[
 \mathcal N(AB)=(\mathcal NA)B+A(\mathcal NB)
 -2\sum_a(\partial_aA)(\partial_aB).
\]

因此对 even F 及 H=sum_a g_a B_a，令 D_size=pNcal−L，严格有

\[
 [D_{\rm size},H_{\rm left}]F
 =(pH-LH)F-2p\sum_a B_a\partial_aF
 -\sum_i(\nabla_iH)(\nabla_iF).
\tag{15}
\]

实时间 U(t)=exp(itH) 又有

\[
 \partial_aU(t)=i\int_0^tU(t-s)B_aU(s)\,ds,
 \quad
 \nabla_iU(t)=i\int_0^tU(t-s)(\nabla_iH)U(s)\,ds.
\tag{16}
\]

在有限 N、Gaussian Sobolev 共同域上，(15) 可以放进 commutator Duhamel。
但它**还不是统一四阶界**：右端带 p 的 Gaussian 项与 Clifford 项要先共同重组，
逐项取范数通常把所需抵消毁掉。old rank=D_rank 使用 R 而非 L；
(15) 不自动把 R 换成 L，也不证明后代的 rank/size L4 强桥。

对任一独立于 g 的 body/rank/size 算符 A，有一个可用的导数恒等式：

\[
 \partial_a[(p\mathcal N-A)F]
 =(p\mathcal N-A+p)\partial_aF.
\tag{17}
\]

(16)–(17) 明确了 (4) 应估计的行/列梯度。一次 Gaussian 导数后不能漏掉 +p。
旧 creation/annihilation 证明工作在 Gaussian L2 tensor HS；其交叉正交并不自动是矩阵序或 S4 估计。

Carlen–Lieb Theorem 4 的半群满足 P_t^Cl Q_A=e^(−t|A|)Q_A，
||P_t^Cl||_{L^r→L^s}≤1 当 e^(−2t)≤(r−1)/(s−1)。
若直接按 physical Clifford degree ≤2pj 使用 L2→L4，它产生 3^(pj)，
仍随 p 发散。因此这个正确的超可积定理本身不能代替 (13) 或梯度抵消。

## 7. 接实际 descendant 的最低证明清单

这不是新物理条件，而是把现有 L2 桥接成四点所需的数学步骤：

1. 针对实际 charged J_d、neutral N_j/其有限有序 Wick 后代，而非只 h,u,x，证明每个固定截断的共同矩与高矩。
   可用完整局部 source 公式；固定次数 exp(B degree²) 界可用于这一步，不能用于全时间尾。
2. 若 actual−proxy 已在 L2 趋零，则再证明它们有统一某个 L^r，r>4。
   插值得 ||actual−proxy||4→0，权重 theta=(r−4)/(2(r−2))；r=8 时 theta=1/3。
   不把 linear B1 的 scalar covariance 证明照抄给 nonlinear 后代。
3. 通过 (4) 的统一梯度余项，或 (12)–(13) 的全方向系数估计，完成 (1)。
   相邻 charge 的实际传递、能量左右权重和既有 binomial 外尾仍须保留。
4. 用四因子 Hölder 传极限：每个替换差乘其余三个 L4 范数。
   unitary center-time 共轭不改变范数；共同 J 的具体 mixed moments 仍由步骤 1 决定。

**当前判断：** 最直接的主文献路线是 Gaussian Sobolev 的 row/column 平方函数；
若能利用局部 CAR collision 把它写成正的低缺额量，就有机会得到真正 uniform-N 的完整尾。
显式 Wick 解耦/flattening 路线则把 degree 常数完全看清，可作为独立备选和固定后代检查工具。
截至本笔记，没有查到一个仅凭现有 L2 rank 界就自动完成该 S4 桥的定理；
也没有在这里补写未经证明的模型特定平方函数或 flattening 估计。

本文是原文定位和可直接验证的有限代数推论，不作这些标准不等式的原创性声明。
未执行科学 main、未覆盖旧证据、未改稿件/PDF、未进行 Git 或外部 issue 操作。
