# 实际 ordinary charged 场的有限 N、全次数 Schatten 尾

2026-11-07。本稿接 [operator draft](susskind-bilocal-fourth-tail-operator-draft.md)
与 [pure/repeated marked contraction proof](susskind-bilocal-fourth-tail-diagram-draft.md)。
模型、common-J、beta=mu=0、normalized annealed full-Fock trace、固定正紧区间内
lambda_N=p²/N、固定 string 实时间窗与前两稿完全相同。
本文先给有限 N 的 uniform degree bound，再交换 degree 截断与 double scaling。
所用 Clifford Riesz 定理已有 [原文逐项核查](../literature/susskind-bilocal-clifford-riesz-sources.md)。

## 1. 目标及三个已证输入

以下记 Ncal 为 coupling Gaussian number，L 为 physical half-Clifford length，D 为 body loss。
三个输入是：

1. 对每个固定 r，所有 sectors 的 ||H_m^n||_r ≤ K_r^n sqrt(n!)，以及
   Gaussian-chaos top projection ||W_n P||_r ≤ 2^n ||P||_r；任意 degree-n
   chaos multiplier 的范数至多 K(1+sqrt2)^n max|multiplier|。
2. 在每个固定内部 charge 区间，纯块 Z_n=W_n(H^n) 上任意固定次
   Delta_D=np−D、Delta_L=np−L 的组合具有
   K(n+1)^b K_r^n sqrt(n!) 的 finite-N 包络。特别需要次数至多2。
3. 实际 covariance residual A_v=p(Phi−C theta_v)Z_v 已有同型包络，
   且有完整 normalized CP path telescoping，而非 coefficientwise 点值近似。

下面证明存在固定 b 与每个固定 r 的常数，使

\[
 \boxed{\sup_N\|G_{n,N}\|_r
 \le K_r(n+1)^b L_r^n\sqrt{n!},\qquad
 G_{n,N}=[p\mathcal N-D-\lambda_N\mathcal N(\mathcal N-1)]H^n.}\tag{1}
\]

更高固定 r 自动覆盖非偶数 r，以下只需证明固定 r≥2。
常数可覆盖 double-scaling 序列有限个小 N 的初始项。

## 2. normalized CP 步骤的精确交换式

在 neutral sector families 上设未归一化 maps

\[
 T X=\sum_i c_i^\dagger Xc_i,\quad
 S X=\sum_i c_iXc_i^\dagger,\quad
 D=Q-T,\quad\bar D=N-Q-S,\quad L=(D+\bar D)/2.
\]

CAR 给 [Q,T]=T、[Q,S]=−S、[T,S]=2Q−N。令
t=Q^(−1)T、s=(N−Q)^(−1)S。在中央区，这些是不同 sector 的 normalized
Schatten contractions。直接计算得

\[
 [D,t]=t-t^2,\quad
 [D,s]=-s+\frac Q{N-Q}ts-\frac{2Q-N}{N-Q}I,\tag{2}
\]
\[
 [\bar D,s]=s-s^2,\quad
 [\bar D,t]=-t+\frac{N-Q}{Q}st+\frac{2Q-N}{Q}I.\tag{3}
\]

对 charge 函数 f，有

\[
 [D,f(Q)]=[f(Q)-f(Q-1)]T,\quad
 [\bar D,f(Q)]=[f(Q)-f(Q+1)]S.\tag{4}
\]

若 P 是长度 ell≤2p 的 normalized CP chain，(2)–(3) 逐项 Leibniz 展开给
在留足 p+2 个 charge 余量的固定中央区

\[
 \|[D/p,P]X\|_{r,m}+\|[L/p,P]X\|_{r,m}
 \le K\frac{\ell}{p}\sup_{|k-m|\le p+2}\|X_k\|_{r,k}.\tag{5}
\]

右边每项只出现常数长度的额外 CP 步骤；它们不产生维数因子。
这一步只是 CP contraction 与有界 scalar coefficients，没有使用 sector Clifford Riesz。

## 3. 实际 A_v 的 D/p 与 L/p 包络

先取 m∈[N/3,2N/3]、p≤N/24；路径及 (5) 所需邻居都留在一个固定更宽中央区。
在 vp≤N/8 时，operator draft 的精确 telescope 是至多2p项之和，形如

\[
 A_v=pC(Q)\sum_{h=1}^{2p} f_h(Q)P_h E_h,\quad
 E_h\in\{(vp-D)Z_v,(vp-\bar D)Z_v\}.\tag{6}
\]

其中 P_h 长度≤2p；f_h 是一个 1/(Q+shift) 或 1/(N−Q+shift)
与至多2p个 1−vp/(Q+shift)、1−vp/(N−Q+shift) 的乘积。
所有被使用的因子在 [0,1]，精确的 charge shifts 保留。
因此

\[
 |f_h|\le K/N,\quad
 |f_h(Q\pm1)-f_h(Q)|\le K(1+v\lambda_N)/N^2.\tag{7}
\]

此处只需粗界 |C(Q±1)−C(Q)|≤Kp/N、C≤C_*。
对 (6) 作用 D/p 或 L/p：

- 作用于 E_h 的项由 repeated-pure 输入控制，因为
  D E_h/p=v E_h−(vp−D)E_h/p；L 同理，barD=2L−D。
- 作用于 P_h 的 commutator 由 (5) 控制，ell/p≤2。
- 作用于 pC f_h 的 charge 差分由 (4)、(7) 控制。
  具体每项额外系数 ≤K N |Delta(C f_h)|
  ≤K[p/N+(1+v lambda_N)/N]。

前两类求和的总系数至多 K p²/N；第三类求和至多
K[p²/N+p(1+v lambda_N)/N]。均只有 v 的固定多项式增长。
故三个量 A_v、DA_v/p、LA_v/p 有同一个形式的 degree 包络。

当 vp>N/8，直接使用 A_v=p(Phi−Ctheta_v)Z_v 的 CP 粗界。
此时 p≤8lambda_N v，N=p²/lambda_N≤64lambda_N v²。
使用 D、L 的 coarse sector-neighbor norm O(N)，得到
||DA_v/p||+||LA_v/p||≤K N max_neighbor||Z_v||，仍为 polynomial(v) K_r^v sqrt(v!)。

外 charge 区不使用 inverse charge denominators；A_v 本身的 coarse bound 为
Kp K_r^v sqrt(v!)，其 D/p、L/p 在相邻 sectors 的 bound 为
KN K_r^v sqrt(v!)。二项外尾≤K exp(−cN) 吸收所有 N、p 多项式。
全过程对原始 A_v 估计，没有插入带跳跃的 sharp cutoff。
由此已组装出可直接输入 full-Clifford 定理的

\[
 \boxed{\|A_v\|_r+\|DA_v/p\|_r+\|LA_v/p\|_r
 \le K_r(v+1)^{b_0}L_r^v\sqrt{v!}.}\tag{8}
\]

## 4. full-Clifford 乘积界与一个 CP 有限范数引理

Ben Efraim–Lust-Piquard 的 Theorem 6.8 和 Lemma 5.6（见原文核查笔记）给
fixed-r、dimension-free 的 Clifford square-function 估计。其直接推论是：
对 neutral X,Y，D 或 L 的 product-rule cross 项 B 满足

\[
 \|B(X,Y)/p\|_r
 \le K_r\left[
 \|X\|_{2r}\|LX/p\|_{2r}
 \|Y\|_{2r}\|LY/p\|_{2r}
 \right]^{1/2}.\tag{9}
\]

这里所有范数是 normalized annealed **full-Fock** 范数；(8) 已先完成该组装。
从未把物理 Clifford OU 当作固定 charge 的 semigroup。

另一个无需维数因子的引理用于 ordinary powers。设 Phi 是任意 CP map，Phi†
相对于同一有限迹取伴随；它可逐 disorder 样本依赖耦合。对正 X,Y 和共轭指数 u,v，
分别对 X,Y 谱分解，再对非负权重 sum_i |(F_i)_{ab}|² 作标量 Hölder，得到

\[
 \operatorname{tr}Y\Phi(X)
 \le[\operatorname{tr}X^u\Phi^\dagger(I)]^{1/u}
      [\operatorname{tr}Y^v\Phi(I)]^{1/v}.\tag{10}
\]

disorder 积分可以一并纳入这个 Hölder。令 u=2r、v=2r/(2r−1)，再对迹用 Hölder
和 L^r duality，便有

\[
 \|\Phi(X)\|_r
 \le\|X\|_{4r}\|\Phi^\dagger(I)\|_2^{1/(2r)}
                  \|\Phi(I)\|_{2r-1}^{1-1/(2r)}.\tag{11}
\]

正负实虚部分分解把 (11) 推广到任意 X，只多至多4的常数。
因此无需把 auxiliary I_N tensor X 的 Schatten 范数错误地当作 X 的范数。

取 force CP map Gamma_H(X)=sum_i[c_i†,H]X[H,c_i]。
H 与 LH/p 的所有固定高范数已由 pure n=1 输入控制；(9) 的原 square-function
估计控制正反 Gram Gamma_H(I)、Gamma_H†(I) 的所有固定范数为 O(p)。
反向 force map 同样成立。于是 (11) 给

\[
 \|\Gamma_H(X)\|_r+\|\bar\Gamma_H(X)\|_r
 \le K_r p\|X\|_{4r}.\tag{12}
\]

## 5. ordinary H 幂及 eta 的微分包络

精确 body product rule 迭代为

\[
 D H^w=pwH^w-
 \sum_{a+b+c=w-2}H^a\Gamma_H(H^b)H^c.\tag{13}
\]

L 的公式则是 w 个单槽 LH 插入与正反 Gamma_H cross 项的平均，保持原次序。
对 (13) 的每项只用固定阶 Hölder，再用 (12) 和 ordinary power 输入；例如
||H^a Gamma_H(H^b) H^c||_r≤K_r p L_r^w sqrt(a!b!c!)。
阶数始终是 r 的固定倍数，与 w 无关。w² 个位置求和后

\[
 \|H^w\|_r+\|DH^w/p\|_r+\|LH^w/p\|_r
 \le K_r(w+1)^2L_r^w\sqrt{w!}.\tag{14}
\]

实际 eta_j=p[Z_j−Q_{j,m}(H)] 的受迫 Jacobi 递推已精确给出

\[
 p r_\ell=\sum_{v+w=\ell-1}A_v\diamond Z_w
          =\sum_{v+w=\ell-1}\mathcal W_{v+w}(A_vH^w),\qquad
 \eta_j=-\sum_{\ell<j}P_{j,\ell,m}(H)\,p r_\ell.\tag{15}
\]

先完全展开这个 recurrence；不在每次迭代把 Hölder 阶数加倍。
由 (8)、(9)、(14)，A_v H^w 及其 D/p、L/p 有
polynomial(v+w) L_r^(v+w) sqrt((v+w)!) 包络。
Gaussian top projection与 D、L 可交换，只增加2^(v+w)。
因此三个 source p r_ell、D r_ell、L r_ell 也都有该包络。

在中央区，P_(j,ell,m) 的系数是长度≤j的 Jacobi path sums。
所有 recurrence 权重 b_k,m=C_m sum_(v<k)theta_(v,m) 一致有界。
其一步 charge 差分的粗界为

\[
 |\Delta b_{k,m}|\le K\frac pN+K\frac{k^2}{N},\tag{16}
\]

其中第二项来自 2p 个 bounded factors 的差分；下一节给出更强版本。
故任意 polynomial coefficient f_m 的界与 NDelta f_m/p 的界都至多
polynomial(j) K^j。利用 (4) 可把 coefficient 先移到乘积外。
再对 H^a(p r_ell) 用一次固定阶的 (9) 与 Hölder。
次数和、Jacobi path 数及 sqrt(a!ell!)≤sqrt(j!) 保证

\[
 \boxed{\|\eta_j\|_r+\|D\eta_j/p\|_r+\|L\eta_j/p\|_r
 \le K_r(j+1)^{b_1}L_r^j\sqrt{j!}.}\tag{17}
\]

这里外 charge 的补全可直接用 Q_j 的全域 coarse bound，不能仅把中央区结论
未经解释送进 (9)。在全域 0≤sum_(v<k)theta_v≤k、C≤C_*，Jacobi 系数满足

\[
 |[x^{j-2k}]Q_j(x)|\le C_*^k\frac{j!}{2^k k!(j-2k)!}.
\]

乘上 ||H^(j−2k)||_r 后，用
sqrt(j!/[k!²(j−2k)!])≤3^(j/2)，得到
||Q_j(H)||_r≤polynomial(j)L_r^j sqrt(j!) 的所有 sector 粗界。
所以外区 eta_j 为 p 乘该包络，Deta_j/p、Leta_j/p 为 N 乘同型包络；
原始完整 fields 的这些外尾由二项指数小量吸收。
在 source products 使用 (9) 时，source 已通过 (8)、(14) 组装成 full-Fock；
在最外层 H^a source 的 (9) 中同样如此。只有 scalar coefficients 的
central/bad-charge估计按区域分开，没有对 sharp-cutoff fields 求导。

## 6. 逆 Jacobi 系数的 sharp charge 差分

唯一需要真正 O(1) charge-derivative 的地方是从 pure Z_j 返回 ordinary H^n。
令

\[
 H^n=\sum_{j=0}^n h_{nj,m}Q_{j,m}(H),\quad
 Q_{j,m}(H)=Z_j-\eta_j/p.\tag{18}
\]

从 xQ_j=Q_(j+1)+b_jQ_(j−1) 得到 h 是长度 n 的 nonnegative weighted
up/down paths；在中央区 b_j≤B，所以 sum_j|h_nj|≤K^n。

现在计算精确 covariance ratio：

\[
 \frac{C_m}{C_{m-1}}
 =\frac m{m-p}\frac{N-m+1}{N-m+p+1},\qquad
 \frac{C_m}{C_{m-1}}-1
 =\frac{p(N-2m+p+1)}{(m-p)(N-m+p+1)}.\tag{19}
\]

写 X_N=(m−N/2)/p，中央区因此有

\[
 N|\Delta C_m|\le K\lambda_N(1+|X_N|).\tag{20}
\]

theta_v,m 是2p个 g_(vp)(x)=(1−vp/x)_+ 的乘积，各分母 x≥cN。
g_b 的 Lipschitz 常数≤b/(cN)²，所有因子≤1；故**任意 v**都有

\[
 N|\Delta\theta_{v,m}|\le K\lambda_N v.\tag{21}
\]

这包含跨越 positive-part 截断点的情况，不需对大 v 失去 N 因子。
又 theta_v≤exp(−c lambda_N v)，故 sum_v theta_v≤A。
由 (20)–(21) 得

\[
 N|\Delta b_{j,m}|\le K(j+1)^2(1+|X_N|).
\]

逐条 path 用乘积差分、path 数≤2^n，得到

\[
 \boxed{\sum_j N|\Delta h_{nj,m}|
 \le K(n+1)^3L^n(1+|X_N|).}\tag{22}
\]

full-Fock 二项 charge 下 X_N 的每个固定阶矩一致有界。
所以 (22) 与邻居 sector 的纯 Z_j 包络可以直接进行 charge Hölder。

## 7. 完整 ordinary charged 组合及尾

保留所有 charge-coefficient 差分，对 (18) 精确作用 pNcal−D，得到

\[
 (p\mathcal N-D)H^n=
 \sum_j\left[
 h_{nj}(pj-D)Z_j-(\Delta h_{nj})T Z_j
 -h_{nj}\mathcal N\eta_j
 +h_{nj}D\eta_j/p+(\Delta h_{nj})T\eta_j/p
 \right].\tag{23}
\]

Delta h_m=h_m−h_(m−1)。没有把 Gaussian number 在 ordinary 内部 pairs 上
替成外部 degree n；(23) 正是这些 pairs 与 body source 抵消后的有限 N 恒等式。

五类项分别受控如下。

- 第一类：pure deficit 输入与 sum|h|≤K^n。
- 第二类：normalized T/m 收缩、(22)、X_N 的高矩和 sector pure bound。
- 第三类：eta_j 的 coupling degree≤j，number multiplier由 §1 的
  Chebyshev estimate 处理，至多 jK^j 倍 (17)。
- 第四类：直接使用 (17)。
- 第五类：只需 (16) 型粗差分 NDelta h/p≤polynomial(n)K^n，结合 (17)
  和相邻 sector 估计。若仅有 full-Fock 范数，也可用 full-trace T 的
  L^r norm≤N 与有界中央系数；邻居权重比在中央区有界。

这些论证在中央区给 polynomial(n)L_r^n sqrt(n!)。
外 charge 可对原始 G_n 直接使用 Ncal、D 的粗范数：
Ncal 与 Ncal(Ncal−1) 只需 n²K^n 的 chaos multiplier，D≤O(N)，
因此 outer contribution≤K(N+pn+lambda_N n²) e^(−cN/r)L_r^n sqrt(n!)。
该界对 N 一致且保持 degree 包络。最后，反项
lambda_N Ncal(Ncal−1)H^n 由同一 Chebyshev multiplier 独立控制。
这证明 (1)。

因此在任意固定紧时间窗，实际 common-J charged 时间场

\[
 G_N(s)=\sum_{n\ge0}\frac{(is)^n}{n!}G_{n,N}
\]

有真正的 finite-N 一致尾

\[
 \boxed{\sup_N\sup_{|s|\le S}
 \left\|G_N(s)-\sum_{n\le K}\frac{(is)^n}{n!}G_{n,N}\right\|_r
 \le K_{r,S}\sum_{n>K}\frac{(n+1)^bL_{r,S}^n}{\sqrt{n!}}
 \longrightarrow0.}\tag{24}
\]

任意固定 s 导数与左右能量幂通过固定高阶 Hölder 和次数平移保持 (24)。
pinching 与真实 unitary 相位均在固定 Schatten 阶数有界。
operator draft 已给实际 charged-force→Q_1G、neutral memory 的高范数强替换，
以及实际返回误差的高范数归约。故这里补上的是完整 relative-time bilocal
共同四点所缺的实际 ordinary 尾；可将既有固定截断的 common-J 多源矩、
回归和 exact-pinching结论先逐截断传递，再以 (24) 及中性同型尾解除截断。

本结论限于既定模型的 normalized annealed joint moments、固定 lambda>0 与
固定实时间紧窗；不据此增加体空间传播、引力相互作用或长时间的声明。

## 8. 独立代数检查的范围

新的 [sparse Fraction driver](../code/check_bilocal_coherence_channel.py) 首次成功输出
[coherence/CP checks](../results/susskind-bilocal-fourth-tail/coherence-channel-checks.json)：
N=2,...,6，p=1,2，a=0,1,2 的全部 admissible charges；340 coherence 恒等式、
265 normalized CP-path 恒等式、75 inactive zero cases、2,609 kernel entries 和
1,305 tagged-site bounds。132 个非常数输入与14个非平凡 CAR gauge sectors。
成功运行0.091秒，输出独占创建。它验证微观 algebra，不替代上面的 degree-uniform证明。

driver SHA256：6b599080c461f8feabcbc075ce2febb46ef7b0756711c186a41acb9d0472455d。
output SHA256：99b2e3d852c5e706853c703310ee8fa2494aa3ac611bef59e5fb4dd89f3e430a。

第二个独立 [charged-reassembly driver](../code/check_bilocal_charged_reassembly.py)
首次成功保存 [输出](../results/susskind-bilocal-fourth-tail/charged-reassembly-checks.json)，
0.528秒。N=2,...,5 的348个完整 neutral matrix units 验证未归一化 CAR maps、
四个 normalized CP commutators 和两种 charge-coefficient Leibniz rules。
另在 (N,p)=(2,1),(3,1),(3,2),(4,2) 的 Hermitian p-body fixtures 上核到次数6：
D/L ordinary-power product rules、pure Hermite number eigenvalues、inverse Jacobi
reconstruction 与 (23) 的全部五项。总计3,312个矩阵恒等式、9,949个逐 entry 比较、
84个 polynomial identities；16个 fixtures 的 charge-difference 项确实非零。

这里 H(g)=gA 只含一个标准实 Gaussian 坐标，number=g∂g−∂g²，所有运算为 Fraction。
它用于独立检查精确代数与符号，不代表完整 original covariance ensemble，亦不检验
渐近范数不等式。28条 per-degree records 都保存实际算出的零残差计数。
新 driver 只导入 sealed coherence driver 的 sparse CAR primitives，未执行其 main。

charged driver SHA256：abc20991b3da5a5ef432b993735bbb1eff026c2b20237b38c1205ed74e1ca283。
charged output SHA256：bb76a4b5b195ecb203612afbf9b285e9d91fd59bb5c57ace47018d65aa18a0b4。
