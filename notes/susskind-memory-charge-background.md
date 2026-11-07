# Whole-field G2 的静态电荷背景

2026-11-07，本地独立核查。主线是 whole-field oriented transfer 的 Gaussian 零 chaos；不是 charged 动态谱的完成。

态为 beta=mu=0 的归一 full-Fock trace，标准 CAR，同一原耦合 J、各电荷使用同一个 w_N，hbar=1。偶 N，p→∞，lambda_N=p²/N 留在正紧区间；明确取 lambda_N→lambda 后才将极限写作 lambda。物理尺度为固定 Js，

\[
w_N=4^p\frac{J_s^2}{\lambda_N}\left(\frac{p!}{N^p}\right)^2,
\qquad C_m=w_N\binom mp\binom{N-m+p}{p},\quad C_0:=C_{N/2}.
\]

这里下标 0 是中央参考值，并非真空扇区；真空方差为零。

## 精确零 chaos

对任何电荷函数 f(Q)，在 m>0 扇区有

\[
(Df(Q))_m=m[f(m)-f(m-1)],\qquad
D X=QX-\sum_i c_i^\dagger Xc_i.
\]

因此，Gaussian 零 chaos 的 H² 是 C_Q，而 whole-field G2=[p Ncal-D]H² 的零 chaos 是 -D C_Q。若 G2 另外包含 -lambda_N Ncal(Ncal-1)H²，或者使用精确 finite-m transfer 乘子，零 chaos 仍相同。D 与 Gaussian 平均可交换；此处没有将 D 换成 rank 或 half-size。

精确二项比给出

\[
\frac{C_{m-1}}{C_m}=\frac{m-p}{m}
\frac{N-m+p+1}{N-m+1},\qquad
(DC_Q)_m=C_m\frac{p(N-2m+p+1)}{N-m+1}.
\]

该式在 m=p 成立，m<p 两边均零；m=0 按 D 的定义为零。

固定中央坐标 m=N/2+xp、x 属于紧集，有

\[
\frac{(G_2^{(0)})_m}{C_0}
=-2\lambda_N(1-2x)+O_X(p^{-1}).
\]

特别是 x=1/2 的极限为零，但有限 N 的 +1 项还在，不能用相对误差衡量该点。

## Full-binomial 的 L² 证明

令 M_N~Bin(N,1/2)，Z_N=(M_N-N/2)/sqrt(N)。精确地 E Z_N=0、E Z_N²=1/4。下面证明

\[
\left\|\frac{G_2^{(0)}}{C_0}
+2\lambda_N-4\sqrt{\lambda_N}Z_N\right\|_{L^2(b_N)}\longrightarrow0.
\]

首先，C_m 关于 m 在 (N+p)/2 附近达到最大值。对 m>=N/2 的上升段，

\[
\frac{C_{m+1}}{C_m}-1
=\frac{p(N-2m-1+p)}{(m+1-p)(N-m+p)}
\le \frac{8p^2}{N^2}\quad(N\ge4p).
\]

上升步数至多 p/2+1，故保守而足够的全局界为

\[
0\le\frac{C_m}{C_0}\le
\exp\!\left(\frac{16p^3}{N^2}\right)=1+O(p^{-1}).
\]

在 |Z_N|<=L 上，两个二项乘积的对数展开给 C_m/C0=1+O_L(1/p)，而精确差分的有理系数为

\[
\frac{p(N-2m+p+1)}{N-m+1}
=2\lambda_N-4\sqrt{\lambda_N}Z_N+O_L(p^{-1}).
\]

两式都是紧参数上一致的。为解除 L，先取广中央区 N/4<=m<=3N/4。精确分母 >=N/4 给

\[
\left|\frac{G_2^{(0)}}{C_0}\right|
\le K\left(8\sqrt{\lambda_N}|Z_N|+4\lambda_N+\frac{4p}{N}\right).
\]

binomial 的一致 sub-Gaussian 尾保证平方一致可积。广中央区外，用全局粗界

\[
\left|G_2^{(0)}/C_0\right|\le Kp(N+p+1),\qquad
\Pr(|M_N-N/2|>N/4)\le2e^{-N/8},
\]

其平方贡献随 N 消失。于是先 N→∞、再 L→∞ 即得上述 L² 近似；不需要从有限计算拟合速率。

由此直接得出，仍保留 finite lambda_N：

\[
\frac{\mathbb E G_2^{(0)}}{C_0}=-2\lambda_N+o(1),
\quad \frac{\mathbb E(G_2^{(0)})^2}{C_0^2}
=4(\lambda_N+\lambda_N^2)+o(1),
\quad \frac{\operatorname{Var}G_2^{(0)}}{C_0^2}=4\lambda_N+o(1).
\]

物理规范中的精确乘积

\[
C_0=\frac{J_s^2}{\lambda_N}
\prod_{j=0}^{p-1}\left(1-\frac{2j}{N}\right)
\left(1+\frac{2(p-j)}{N}\right)
=\frac{J_s^2}{\lambda_N}[1+O(p^{-1})]
\]

给 C0→C=Js²/lambda。于是 full trace 的均值、二阶矩、方差分别趋于 -2lambda C、4C²(lambda+lambda²)、4lambda C²。该场是电荷函数，逐样本守恒；不是新的动态谱。

## 独立的有限组合校验

binomial 移位恒等式给

\[
\mathbb E_{b_N}DC_Q=\mathbb E_{b_N}[(2M_N-N)C_{M_N}].
\]

另一种计数按内部收缩的两个 p-set 的交集 r=|I intersect J| 分组。每项固定 p 个 occupied 位点及 p-r 个 empty 位点。令

\[
a_r=w_N\binom Np\binom pr\binom{N-p}{p-r}2^{-(2p-r)}.
\]

则精确地 E C_M=sum_r a_r，E D C_Q=sum_r r a_r。新 driver 独立计算 binomial 和 overlap 两边；81 个系统、1643 个扇区通过精确 Fraction 检验。

## 两个内部 Gaussian 弦的 pair-source

此部分只检查 G2 零 chaos 的必要校准，不代替一般 whole-field 图定理。

取单点 c=[[0,1],[0,0]]、n=c†c、bar n=1-n、tau=Tr/2。对 z_i²=0，标记 -D 给

\[
L_i(n)=(1-z_i)n,\qquad L_i(\bar n)=\bar n+z_i n.
\]

两个 field 发生共同站点碰撞时，必须除以各自单点 trace，并在 nilpotent 环内求逆：

\[
R_{XY}=\frac{\tau[L_1(X)L_2(Y)]}
{\tau L_1(X)\,\tau L_2(Y)}.
\]

独立 2x2 计算得到

\[
R_{nn}=2,\quad R_{\bar n\bar n}=2-2z_1-2z_2+4z_1z_2,
\quad R_{n\bar n}=2z_2,\quad R_{\bar n n}=2z_1.
\]

四种跨 field 碰撞的 sum(R-1) 仅余 4z1z2。每个 field 内部 I/J 自碰撞的原局域词是 n，不是 n bar n，因此

\[
R_\mathrm{own}(z)=\frac{\tau L(n)}{\tau L(n)\tau L(\bar n)}
=2(1-z).
\]

用未标记 C² 归一后，pair-source 指数为

\[
\kappa=-2z_1-2z_2+4z_1z_2,
\quad [z_1z_2]e^{\lambda\kappa}=4\lambda+4\lambda^2.
\]

两个单标记均值为 -2lambda，connected 系数为 4lambda，完全恢复上面从实际 finite C_Q 得到的 full-trace 极限。

这两个 H² field 之间没有 Gaussian 配对桥，但共享站点仍产生静态电荷涨落。charged whole-field 路线不能直接套用 neutral Wick-source 的“无桥项全部被能量回归扣除”结论；须保留电荷背景并进一步构造正确回归。此处未宣布 AA/AB/BA 核或一般 pinching 已完成。

## 保存与数值范围

新 driver：`code/check_memory_charge_background.py`；新科学输出：`results/susskind-memory-oriented-transfer/charge-background-checks.json`，首次运行成功，拒绝覆盖。

另做 12 组完整 binomial 求和、72 个中央点，70/110 位独立重算；m<p 项严格为零，其余全部扇区求和，没有尾截断。计算保留 lambda_N，其中包含 N=4p²+2p 的漂移序列。p=128、lambda_N=1/4 的方差相对偏差仍约 1.7257%，lambda_N=1 时约 8.6505%；原始偏差完整保留。高精度一致性不是区间证书，有限点也不承担上述极限证明。
