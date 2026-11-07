# 第四范数尾：OU 投影、完全正 covariance 路径与实际中性场

2026-11-07，新阶段工作稿。模型为 common-J overlapping complex p-body CAR，normalized annealed full-Fock trace，beta=mu=0，固定 string 时间窗，偶 N、lambda_N=p²/N 留在正紧区间并趋 lambda>0，Js>0 固定。这里所有 L^r 都是 Schatten 阶数，不是 Gaussian-number Sobolev 阶数。

本稿的完整新结果是 **实际 neutral memory 分量的所有固定 Schatten 阶数统一可和尾**，以及它到既有局域代理的强高范数替换。关键是把确定性 covariance map 写成逐步完全正 charge-transfer 通道，绕过未建立的非线性 rank-map 高范数不等式。本稿末节记录的 charged ordinary 缺口，随后由 [独立 charged 证明](susskind-bilocal-fourth-tail-charged-proof.md) 闭合；该证明额外控制 CP 路径的 D/L 导数、实际 eta 与 charge 系数差分，不以纯 Wick 的尾替代 ordinary exponential 的尾。

## 1. 不改变范数阶数的 OU/Chebyshev 工具

对矩阵值 Gaussian 多项式 P，令
\[
 (\mathsf O_tP)(g)=\mathbb E_{g'}P(tg+\sqrt{1-t^2}g'),\qquad -1\le t\le1.
\]
Gaussian 测度不变性与条件 Jensen 给每个 r≥1 的
\(\|\mathsf O_tP\|_r\le\|P\|_r\)。这里从头到尾是同一个 r，不是 L² 到 S_r 的 hypercontractivity。

若 coupling degree(P)≤n，则
\[
 F(t)=\mathsf O_tP=\sum_{k=0}^n t^k\mathcal W_kP
       =a_0+\sum_{j=1}^n a_jT_j(t),
\]
其中 T_j 是 Chebyshev 多项式。积分公式给
\(\|a_0\|_r\le\|P\|_r\)、\(\|a_j\|_r\le2\|P\|_r\)。T_n 的最高系数为 2^(n−1)，所以
\[
 \boxed{\|\mathcal W_nP\|_r\le2^n\|P\|_r.}\tag{1}
\]
T_j 的系数绝对值之和等于 \(|T_j(i)|\)，至多 \(R^j\)，\(R=1+\sqrt2\)。故任何确定性有界 chaos 乘子满足
\[
 \boxed{\left\|\sum_{k\le n}\theta_k\mathcal W_kP\right\|_r
 \le K R^n\max_{k\le n}|\theta_k|\,\|P\|_r.}\tag{2}
\]
常数 K 不依赖矩阵维数、n 或 Gaussian 坐标个数。逐 charge 的有界 theta 也可使用。此界避免固定-degree Vandermonde 插值在 degree 上的不必要坏常数。

由所有 sector 的 \(\mathbb E H_m^2=C_mI\)、\(C_m\le C_*\)，Gaussian 矩界给固定 r 下
\[
 \|H_m^n\|_r\le K_r^n\sqrt{n!},\qquad
 \|Z_{n,m}\|_r\le(2K_r)^n\sqrt{n!},\quad Z_n=\mathcal W_n(H^n).
\tag{3}
\]
第二式只使用 (1)，无需一般矩阵 chaos 的 L²→L^r 断言。

## 2. 本稿调用的纯标记界

[独立 tagged contraction 证明](susskind-bilocal-fourth-tail-diagram-draft.md)先保留 Gaussian pairing 和 CAR 符号，再用 ancilla row/column covariance 范数控制每张图。它给出任意固定 r、任意次数 v，在固定内部 charge 区间上的 sector 一致界
\[
 \|(vp-D)Z_v\|_{r,m}+\|(vp-\mathcal L)Z_v\|_{r,m}
 \le K_r(v+1)^2L_r^v\sqrt{v!}.\tag{4}
\]
其中 D 是 body-loss，\(\mathcal L\) 是 half-Clifford size。这里的区间可取 N/4≤m≤3N/4；固定数目 CAR 插入只使端点常数改变。它不是从全 charge 平均反推出 sector 上界。

定义 reverse loss
\[
 \bar D X=(N-Q)X-\sum_i c_iXc_i^\dagger=2\mathcal LX-DX.
\]
故 \((vp-\bar D)Z_v=2(vp-\mathcal L)Z_v-(vp-D)Z_v\) 也满足 (4)。整个证明只需这些纯 Gaussian-Wick 场；不预设 ordinary charged G_n 的有限 N 同型界。

## 3. covariance map 是精确的 down/up 完全正通道

令 \(\mathscr H_k=\wedge^k\mathbb C^N\)。两个相邻 charge 通道是
\[
 T_k:B(\mathscr H_{k-1})\to B(\mathscr H_k),\qquad
 T_kX=\frac1k\sum_i c_i^\dagger Xc_i,
\]
\[
 S_k:B(\mathscr H_{k+1})\to B(\mathscr H_k),\qquad
 S_kX=\frac1{N-k}\sum_i c_iXc_i^\dagger.
\tag{5}
\]
二者都 unital。使用
\(d_{k-1}/d_k=k/(N-k+1)\)
和其相邻版本，直接核得它们也保持各自 **归一 sector trace**。因此是两个不同有限矩阵代数之间的 bistochastic CP maps；对所有 Schatten-r，包含 disorder 平均后，
\[
 \|T_kX\|_{r,k}\le\|X\|_{r,k-1},\qquad
 \|S_kX\|_{r,k}\le\|X\|_{r,k+1}.\tag{6}
\]

设原 Gaussian covariance 为
\[
 \Phi_m(X)=w_N\sum_{I,J:|I|=|J|=p}
 C_I^\dagger C_J X C_J^\dagger C_I,
 \quad C_m=w_N\binom mp\binom{N-m+p}p.
\]
先连续 p 次 S 从 m 降到 m−p，再连续 p 次 T 回到 m，得到精确式
\[
 \boxed{\Psi_m:=\Phi_m/C_m
 =T_m\cdots T_{m-p+1}\,S_{m-p}\cdots S_{m-1}.}\tag{7}
\]
证明只需展开 Kraus products。每个无序 p-set 的 p! 个顺序在左右各出现一次，CAR 重排符号成对相消；分母
\((m)_p(N-m+p)_p=(p!)^2\binom mp\binom{N-m+p}p\)
恰好消去这两个 p!。因此 \(\Psi_m\) 对同一个样本的随机 X 也是确定性 Schatten contraction。

对任何 neutral full-Fock lift X，每步还有精确误差恒等式
\[
 T_kX_{k-1}-(1-b/k)X_k=\frac{[(b-D)X]_k}{k},
\]
\[
 S_kX_{k+1}-(1-b/(N-k))X_k
     =\frac{[(b-\bar D)X]_k}{N-k}.\tag{8}
\]
这些是完整算符恒等式；它们保留原 Gaussian 系数与 Kraus 结构，避免仅凭单个 monomial 的小点值推导高 Schatten 范数。

## 4. 实际 Phi 缺陷的全次数高范数界

定义
\[
 \theta_{v,m}=\mu_m(vp),\quad
 \mu_m(b)=\prod_{j=0}^{p-1}
 (1-b/(m-j))_+(1-b/(N-m+p-j))_+,
\]
\[
 A_{v,m}=p(\Phi_m-C_m\theta_{v,m})Z_{v,m}.\tag{9}
\]
取原 charge 的 good 区 N/3≤m≤2N/3，并取充分大 N 使 p≤N/12。路径 m→m−p→m 全部留在 N/4≤k≤3N/4，正好落入 (4) 的 sector 一致范围。

先设 vp≤N/8，令 (8) 中 b=vp。所有 2p 个 scalar reference 因子均在 [0,1]，其乘积正是 \(\mu_m(vp)\)。沿路径 telescoping；transported error 由 (6) 收缩，剩余 scalar factors 的模至多一，每个分母≥N/4。于是
\[
 \|\Psi_mZ_{v,m}-\mu_m(vp)Z_{v,m}\|_{r,m}
 \le\frac{8p}{N}
  \sup_{m-p\le k\le m}
  \max\{\|(vp-D)Z_v\|_{r,k},\|(vp-\bar D)Z_v\|_{r,k}\}.
\tag{10}
\]
乘 pC_m 后，\(p^2/N=\lambda_N\) 和 (4) 给所需一致界。

若 vp>N/8，不沿用可能变负的 scalar factors。此时
\(p<8\lambda_Nv\)，直接用 (6)--(7)、\(|\theta_v|\le1\) 和 (3)：
\[
 \|A_{v,m}\|_{r,m}\le2pC_m\|Z_{v,m}\|_{r,m}
 \le K_r(v+1)L_r^v\sqrt{v!}.
\]
两区合并即
\[
 \boxed{\|A_{v,m}\|_{r,m}
 \le K_r(v+1)^2L_r^v\sqrt{v!}
 \quad(N/3\le m\le2N/3),}\tag{11}
\]
对所有 degree v、所有充分大 N 同时成立。此步没有对 actual rank map 使用高范数估计。

bad charge 区使用同一 CP 粗界 \(\|A_{v,m}\|_r\le2pC_m\|Z_v\|_r\)。其二项质量≤2exp(−N/18)，所以全迹贡献受
\(K_rp e^{-N/(18r)}L_r^v\sqrt{v!}\)
控制。吸收到常数后 (11) 也给全 charge 的 degree-uniform 界。inactive m<p 的场直接置零。

## 5. actual Jacobi 误差的全次数范数

以下在各 sector 保留有限参数 \(C_m,\theta_{v,m}\)，简写 C,theta。令
\[
 a_j=\sum_{v<j}\theta_v,\quad
 Q_{j+1}(x)=xQ_j(x)-Ca_jQ_{j-1}(x),\quad
 \eta_j=p[Z_j-Q_j(H)].
\]
已有精确 Gaussian annihilation 恒等式为
\[
 p r_j=\sum_{v<j} A_v\diamond Z_{j-1-v},
 \qquad
 \eta_{j+1}=H\eta_j-Ca_j\eta_{j-1}-p r_j,
 \quad\eta_0=\eta_1=\eta_2=0.\tag{12}
\]
\(A_v\) 保持纯 coupling degree v。故
\[
 A_v\diamond Z_w=\mathcal W_{v+w}(A_vH^w),
\]
用 (1)、固定的 2r 阶 Hölder、(3)、(11) 即得
\[
 \|A_v\diamond Z_w\|_r
 \le2^{v+w}\|A_v\|_{2r}\|H^w\|_{2r}
 \le K_r(v+1)^2L_r^{v+w}\sqrt{(v+w)!}.\tag{13}
\]
因此 \(p r_j\) 有 polynomial(j)·L_r^j·sqrt(j!) 包络。

不能在 (12) 中每迭代一步就把 Hölder 阶数加倍。正确办法是先完全展开强迫递推：
\[
 \eta_j=-\sum_{\ell<j}P_{j,\ell}(H)\,p r_\ell.
\]
\(P_{j,\ell}\) 次数 d=j−1−ell，系数满足
\[
 |[x^{d-2s}]P_{j,\ell}(x)|
 \le(C_*A)^s\binom{d-s}s,\qquad 1\le a_j\le A.
\]
现在只使用一次 2r 阶 Hölder 控制每个 \(H^{d-2s}p r_\ell\)，而 (13) 只再用固定 4r 阶。阶乘合并
\(\sqrt{(d-2s)!}\sqrt{(\ell-1)!}\le\sqrt{j!}\)，有限系数和只按 j 指数增长。于是
\[
 \boxed{\|\eta_{j,m}\|_{r,m}\le K_r(j+1)^bL_r^j\sqrt{j!}}
 \tag{14}
\]
在 good 区一致，b 为某个固定整数（例如放宽取 b=6 足够，不优化）。

bad 区不使用 (11) 的精细路径界；直接从 \(\eta_j=p(Z_j-Q_j(H))\)、(1) 及 Q_j 的系数界得到
\(\|\eta_j\|_{r,m}\le pK_r(j+1)^bL_r^j\sqrt{j!}\)。\(a_j\le A\) 的几何和界也适用于这里：若某正部因子变零，theta_j 只会更小。二项尾再次吸收 p，因此 (14) 具有全 charge 版本。

## 6. 完全绕过 rank 的 actual neutral 时间级数

有限参数 Jacobi 系数满足全复平面的精确展开
\[
 e^{isH}=\sum_{j\ge0}c_{j,m}(s)Q_j(H),\qquad
 |c_{j,m}(s)|\le e^{2\sqrt{AC_*}S}(2\sqrt A S)^j/j!.
\tag{15}
\]
减去逐样本精确守恒函数
\(pC\sum_jc_j(s)\theta_jQ_j(H)\)，得到
\[
 \boxed{p(1-\mathcal P_h)\Phi_m(e^{isH})
 =(1-\mathcal P_h)\sum_{j\ge0}c_j(s)T_{j,m},}
\]
\[
 T_{j,m}=A_{j,m}+(C_m\theta_{j,m}-\Phi_m)\eta_{j,m}.\tag{16}
\]
因为 \(\Phi_m/C_m\) 是同一 sector 的 Schatten contraction，
\[
 \|(C\theta_j-\Phi)\eta_j\|_r\le2C_m\|\eta_j\|_r.
\]
结合 (11)、(14) 得真正有限 N、全次数的 actual 系数包络
\[
 \|T_{j,N}\|_r\le K_r(j+1)^bL_r^j\sqrt{j!}.\tag{17}
\]
于是 (16) 在任意固定实时间紧窗绝对一致可和，尾由
\[
 K_{r,S}\sum_{j>K}(j+1)^b L_{r,S}^j/\sqrt{j!}\longrightarrow0
\tag{18}
\]
支配，且对 N 一致。这不是先 N 后 j 的逐次数极限包络。

原实际 neutral memory 分量为
\[
 B_{m,N}(s)=\frac{m\epsilon_m}{p}(1-\mathcal P_h)e^{-ihs/2}
       \sum_jc_{j,m}(\alpha_ms)T_{j,m}e^{-ihs/2}.
\tag{19}
\]
\(m\epsilon_m/p\le2\)、真实相位是 unitary，\(1-\mathcal P_h\) 在 S_r 的范数≤2。因此 (18) 直接成为 **实际 B_N 的** 全 charge 高范数尾，不需假定 Gaussian-chaos 乘子或 rank multiplier 守恒。

## 7. 实际中性场的局域多点极限

前驱 neutral-marked-series 已证明实际 B_N 到 finite-parameter local-size/ordered-Wick 代理的紧时间 L² 强替换。新的纯标记界 (4)、(1)--(2) 和同样的展开递推给该局域代理所有固定 S_r 的可和尾。另一方面 (19) 已给实际 B_N 相同的高范数一致界。

所以旧 L² 强差通过固定更高阶范数插值升为每个固定有限 S_r 的紧时间强差。对四点可用 r=4、上界阶数 8。随后先对有限截断使用任意固定多源共同矩及其 exact pinching 传递，再由 (18) 和 Hölder 去掉截断。这闭合 actual neutral memory 分量的全有限时间共同四点；其 Volterra 滤波后的 neutral bilocal 分量同样成立。任何 fixed charged-generator 插入可使用其已有高矩；完整 charged ordinary 时间场仍需末节列出的另一个尾接口。

其中 exact pinching 仍需按正确顺序处理：固定后代的 L² 无对角谱给 pinching 残量趋零，新高矩和收缩性将其升为所需 S_r；不能单凭极限多源矩声称有限 N pinching 可去。

## 8. OU number remainder 和 actual charged 强接口也可升阶

沿 oriented-transfer 定义
\[
 E_m=[m(1-\alpha_m^{\mathcal N})-D]e^{isH}|_m,
 \quad G_N=[p\mathcal N-D-\lambda_N\mathcal N(\mathcal N-1)]e^{isH}.
\]
D 在差中精确抵消。good 中央区已有
\[
 |R_m(k)|\le a_m(k+1)^3,\quad
 a_m=4\lambda_+|m-N/2|/N+8\lambda_+^2/(3p).
\]
用 (2) 作用到 H^n，再用 (3)，以及二项 \(\|a_m\|_{L^r(b)}=O_r(p^{-1})\)，得到任意固定 r、时间导数和左右能量权重下
\[
 \sup_{|s|\le S}\|H^a\partial_s^j(E-G_N)H^b\|_r=O(p^{-1}).\tag{20}
\]
外区粗 number multiplier 为 O(p)(k+1)²，被二项尾吸收。其次数包络只是 polynomial(n)·L_r^n·sqrt(n!)，因此全时间求和合法。

还可控制 \(\alpha\to1\) 的强误差。定义随机 CP force map
\[
 \Gamma_H(X)=\sum_i[c_i^\dagger,H]X[H,c_i].
\]
body product rule 和 Gaussian Euler/Laplacian 给精确式
\[
 (p\mathcal N-D)H^n
 =\sum_{a+b+c=n-2}H^a[\Gamma_H(H^b)-2p\Phi(H^b)]H^c.\tag{21}
\]
其指数形式是带负号的两重 Duhamel 积分。对 unitary U，矩形 Hölder 给
\(\|\Gamma_H(U)\|_r\le\|\sum_i[H,c_i]^\dagger[H,c_i]\|_r=O_r(p)\)；\(\|\Phi(U)\|_\infty\le C_*\)。固定导数/能量插入用更高固定阶 Hölder 和各 force 的 Gaussian covariance 界仍为 O(p)。数目算符项由 (2) 控制。因此 G_N 的这些所有固定高范数有 O(p) 粗界。

与旧 weighted L²=O(1) 插值，S_r 范数为 \(O(p^{(r-2)/(r-1)})\)（r≥2，使用 2r 阶粗界）。中央 \(1-\alpha=O(1/p)\)，外区以更高范数和二项尾处理，故 actual charged force 到 \(\mathscr Q_1G_N\) 的强替换在 S_r 中成立，速率可取 \(O(p^{-1/(r-1)})\)，r=4 为 O(p^-1/3)。它没有自动证明 G_N 本身的统一第四范数或 ordinary 时间尾。

## 9. actual memory 到两分量的高范数归约

旧力任意固定高矩给 pLambda(s) 的全 charge 所有固定 S_r 一致界；其中 Lambda=R†exp(iKs)R。已有 pLambda→b_s(H) 的紧时间 L² 强极限，目标 b_s(H) 为有界 energy function。因此插值将该强极限升到每个固定 S_r。

第一次返回乘 m 后写为两个 pLambda 的乘积，前因子 m/p² 一致有界；选取 2r 阶 Hölder，替成保留 m、alpha 的精确守恒 reference。下一返回的矩形 Hölder 界为 mG³，旧力第 3r 阶矩给 O(p^-1)，外区只有多项式增长与指数尾。于是实际 memory 与普通 Lambda 项的差在 pinching 后于每个固定 S_r 趋零。

OU force 展开的四 auxiliary-Gaussian 项已有紧时间有界高范数（前驱甚至给 samplewise 有界的条件平均），其 energy regression 的旧 L² 强差同样升阶。六 auxiliary-Gaussian 因子的余项用固定 6r 阶 Hölder，乘 m epsilon³ 后中央为 O(p^-1)，外尾消失。因此旧实际 memory→charged A+neutral B 的强归约也在固定 S_r 中成立。新的困难集中在 charged ordinary 尾，而不是这些实际返回/两分量接口。

## 10. 精确 coherence Markov 表示（独立微观接口）

取 disjoint exclusive sets A,B，|A|=|B|=a；rest 有 L=N−2a 个 sites，共同 occupation R 有 M=m−a 个粒子。以 \(T_{AB}=C_A^\dagger C_B\) 的真实 CAR gauge 写
\(X=T_{AB}f(n_{\rm rest})\)。Phi 的一个 Kraus 对只能从 R 的共同粒子移走 p 个 J，再在移除后 rest 的空位中放回 p 个 I。不能碰 exclusive A/B，否则对应 bra 或 ket 振幅为零。与 T_AB 交换这些 rest CAR 字母后，相对 gauge 的全部权重非负。

令 P_a 是上述 p-removal/p-addition Markov chain，若 M≥p 则
\[
 \kappa_a=w_N\binom{m-a}p\binom{N-m-a+p}p=C_m\mu_m(a),
\]
\[
 \boxed{\Phi_m[T_{AB}f]=\kappa_a T_{AB}P_af.}\tag{22}
\]
若 M<p，左边与 kappa_a 均为零。具体 transition：对 |R|=|R'|=M、k=|R\setminus R'|，
\[
 P_a(R,R')=\frac{\binom{M-k}{p-k}}
                 {\binom Mp\binom{L-M+p}p},\qquad 0\le k\le p,
\]
其余为零。它是对称且 stochastic 的 kernel。

对只依赖 c 个 rest sites 的有界 f，coupling 原 occupation 与一步结果给
\[
 |P_af-f|\le2\|f\|_\infty\,c
       [p/M+p/(L-M+p)].\tag{23}
\]
若 f 为 0/1 pattern indicator，可删去前面的 2。一项 raw coefficient product 的 exclusive rank a 与 diagonal-site 数 c 满足
\(vp-a\le R_{\rm overlap}\)、\(c\le R_{\rm overlap}\)，所以 (23) 确实识别了 collision 导致的 O(p/N) 小量。

单凭这个 coefficientwise 点值小量不足以推出 summed Gaussian Schatten 界。上文真正闭合此缺口的是 (7)--(10) 的 **完整 CP 路径 telescoping**，它保留全部 Gaussian row/column 结构，不把点值 majorization 当作任意矩阵系数 tensor 的范数不等式。本节独立的有限 CAR/Markov/CP-path 检查另存新输出，不覆盖任何旧 driver。

## 11. 本稿留下、随后闭合的 charged ordinary 尾

纯 Wick 的 (4) 允许所有次数统一求和；ordinary \(G_n\) 的内部 Gaussian pair 可让一个 physical tag 只约束同一个独立 p-set，粗贡献变为 p。必须在有限 N 先把这些内部 tags 与 number/source counterterms 精确组合，才能恢复 degree-uniform 的 order-one 界。

已知逐 degree 的 \(\limsup_N\|G_n\|_4\) 可和包络只构造 limiting charged 四范数级数，不等于实际有限 N 的 uniform tail。上文 (20)--(21) 补上强接口，但本身未改变这个极限次序。本稿因此只将中性部分与实际传递闭合。

随后完成的 [charged finite-N tail 证明](susskind-bilocal-fourth-tail-charged-proof.md) 使用本稿的实际 A_v telescope，先证明 D A_v/p、L A_v/p 的同型包络，再用已核原文的 full-Clifford Riesz 定理处理乘积项，控制实际 eta_j 及 inverse-Jacobi charge 差分，最后精确重组 ordinary G_n。这一步才补上完整 actual finite-relative-time bilocal 共同四点需要的 charged uniform tail。
