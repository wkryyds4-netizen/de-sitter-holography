# Joint charge tilt：singleton 修正、共同矩与一致可积性

2026-11-07，新 joint-regression 阶段的独立推导。沿用已封存 oriented-transfer 的 fixed-degree 共同 J 图及局域 nilpotent 规则；本记录不修改旧证据，不以有限计算替代渐近证明。

标准 CAR、overlapping complex Gaussian p-body Hamiltonian、各电荷同一原 J 和共同方差 w_N、beta=mu=0 的归一 annealed full-Fock trace。偶 N，p→∞，lambda_N=p²/N 留在正紧集并趋于 lambda>0，Js 固定，C=Js²/lambda，q=exp(-4lambda)。这里 t 是无量纲 charge source，不是物理时间。所有 word 次数、X 的多项式次数、标记数在 N 极限前固定。

定义 X_N=(Q-N/2)/p。讨论一个或两个标记；D 标记作用于整个 field，允许 Gaussian 内部配对，源为 pd-D 或 G 源 pd-D-lambda_N d(d-1)。neutral L 标记仍只作用于 Gaussian-Wick 禁内部配对区中的连续子块，不扩大其作用区。

## 1. 精确 tilted occupation 和遗漏的未标记 singleton

对任意复 t，有限 N 的未标记 charge MGF 是

\[
M_N(t)=\operatorname{tr}e^{tX_N}=\cosh(t/(2p))^N.
\]

在任意固定复圆盘 |t|<=R，p 充分大时 cosh 的分母无零点；可使用精确 normalized product-state occupation

\[
\rho_t=\frac{e^{t/p}}{1+e^{t/p}}
=\frac{1+\tanh(t/(2p))}{2},\qquad
\tau_t(\operatorname{diag}(a,d))=(1-\rho_t)a+\rho_t d.
\]

这只是复解析的计算工具；只有实 t 时才是概率态。

固定 Gaussian pairing 有 ell 条弦、2ell 个 p-set labels。未标记且无 site 重叠时，每弦有 p 个 n singleton 与 p 个 bar n singleton。相对于 t=0 的乘积，不能在有限 N 直接删去

\[
S_{\ell,N}(t)=[4\rho_t(1-\rho_t)]^{p\ell}
=\cosh(t/(2p))^{-2p\ell}.
\]

它与 M_N 精确合成 cosh(t/(2p))^(N-2pell)。在 |t|<=R 上

\[
\log S_{\ell,N}(t)=-\frac{\ell t^2}{4p}+O_R(p^{-3}),\qquad
\log M_N(t)=\frac{t^2}{8\lambda_N}+O_R(p^{-2}).
\]

因此遗漏 S 的影响是 O_R(1/p)，不是严格零。它在先完成下面的精确 source 抵消后可以在极限中置为 1；若仍保留分别为 O(p) 的源项，则不能凭这个量级擅自丢弃它。

可用的显式复数界为：若 z=t/(2p)、|z|<1，则 cosh 的无穷乘积给

\[
|\log\cosh z|\le\frac{|z|^2}{2(1-|z|^2)},\qquad
|\log\cosh z-z^2/2|\le\frac{|z|^4}{4(1-|z|^2)}.
\]

这些界足够控制上面两个余项，不要求 t 为实数。

## 2. Source singleton 的精确有限 p 化简

对 z_a²=0，局域对称 D loss 满足

\[
e^{-z_aD}c^{(\dagger)}=e^{-z_a/2}c^{(\dagger)},\quad
e^{-z_aD}n=(1-z_a)n,\quad
e^{-z_aD}\bar n=\bar n+z_an.
\]

开放于该 D field 的弦在两个 labels 上各有 odd factor e^(-z_a/2)。总 2pd_a 个 singleton 恰抵消 e^(pd_a z_a)，对任意 rho_t 都成立。neutral L 子块的 incident 弦也只有该种 odd 缩放，因此其 source singleton 同样精确抵消，并不产生 charge-linear 项。

一个 D-field 内部弦的两个 singleton，在除以未标记 singleton 后给

\[
(1-z_a)\left(1+\frac{\rho_t}{1-\rho_t}z_a\right)
=1+(e^{t/p}-1)z_a.
\]

每个 label 有 p 个站点。因此 i_a 条内部弦的**精确**残因子为

\[
\exp_{\rm nil}\{i_a\,p(e^{t/p}-1)z_a\},
\]

其中 nilpotent 指数可按有限多项式理解。一个 field 的一阶系数就是 i_a p(e^(t/p)-1)，不是有限 p 下的 i_a t。对一个或两个 sources，

\[
|p(e^{t/p}-1)-t|\le\frac{|t|^2}{2p}e^{|t|/p}
\]

保证所有所需 nilpotent coefficients 在复 t 紧集上一致趋于 exp_nil[t sum_a i_a z_a]。G 源的 -lambda_N d_a(d_a-1) 是固定次数的额外常数 shift，不改变该论证。neutral L source 在这条式中取 i_a=0；不能把较大 neutral field 的其他内部 contractions 误算为 L 的内部弦。

## 3. 重叠项的支配与复 t 局部一致收敛

带源碰撞比值必须除以两个或多个**带源 singleton**。先作第 2 节的化简，再处理碰撞。固定 label 数时，其各 nilpotent coefficients 都是 rho_t 的有理函数；rho_t 在复 t 紧圆盘上一致趋于 1/2，分母远离零。因此每个 pair-site coefficient 相对 rho=1/2 的变化为 O_R(1/p)。局部因子为零时也如此：比较的是完整 signed 多项式，不能把局部 trace 为零的图先删掉。

令 R_tot 为全部 label-pair intersections 的总数。对固定 label 数、至多两个标记，精确单点展开有系数绝对值上界

\[
K_R(1+R_{\rm tot})^2 A_R^{R_{\rm tot}},
\]

其中已抵消的 pd 因子不再留在右端。一个 pair overlap Y 是 Hypergeometric(N,p,p)。对 a>=0，其 factorial-moment 展开给

\[
\mathbb E e^{aY}
\le\left[1+\frac pN(e^a-1)\right]^p
\le e^{\lambda_N(e^a-1)}.
\]

有限个 pair 的和再用 Holder，得到 R_tot 的任意固定指数矩一致有界。故固定 pair 比值的 O_R(1/p) 改变可以在求和/期望后控制。

任意三 labels 同站点的概率至多常数乘 p³/N²=O(1/p)。坏事件上的上述支配配 Cauchy--Schwarz 已足以给 o(1)；这里不借此声称整个 marked 图误差一定具有 O(1/p) 速率。旧 fixed-label pair-Poisson 极限与指数矩控制因此在整个复 t 紧圆盘上一致适用。物理 C 的有限 N normalization 也沿用旧极限；若 lambda_N-lambda 没有给定速率，不能为最终目标写统一 O(1/p) 收敛速率。

结论：对一个固定 pairing，省略其原有 C^ell 权重，one/two-source 的极限是

\[
\exp\!\left(\frac{t^2}{8\lambda}\right)
\exp_{\rm nil}\!\left\{
\lambda\left(\kappa_0+\sum_a\kappa_a z_a+\kappa_{12}z_1z_2\right)
+t\sum_a i_a z_a\right\}.
\]

对 G，kappa_a 使用已减 d_a(d_a-1) 的值；对 L，其 i_a=0，kappa 来自原有 neutral mark 区。有限个允许 pairing 求和保留该局部一致性。

## 4. 固定 X 多项式共同矩与真正的 UI

有限 N 的 F_N(t)=E tr[e^(tX_N) W_N] 是整函数；固定 word 的 disorder 平均是有限 Gaussian 矩。第 3 节给出复 t 局部一致收敛，故 Cauchy 导数公式允许任意固定 r：

\[
\partial_t^rF_N(0)=\mathbb E\operatorname{tr}[X_N^rW_N]
\longrightarrow\partial_t^rF(0).
\]

仅有实 t 的点态收敛不足以支持这一步。这里 W_N 可含外部 H powers、一个或两个 whole-field G sources，或与 neutral Wick-subblock L source 的混合；次数固定。

还可以给出不用符号矩替代的正范数 UI。令 P_N 是固定系数、固定次数的多项式 word，至多包含一个 G 或 neutral source，可有有限个左右 H 因子。P_N†P_N 的各项正好在 two-source 图范围内。对任意固定 a>0，复局部界特别给

\[
\sup_N\mathbb E\operatorname{tr}[e^{\pm aX_N}P_N^\dagger P_N]<\infty.
\]

X_N 与所有这些 charge-preserving operators 对易，而且 e^(a|X|)<=e^(aX)+e^(-aX)。因此对任意固定 r，

\[
\sup_N\mathbb E\operatorname{tr}
[|X_N|^r\mathbf1_{|X_N|>L}P_N^\dagger P_N]
\le K_{r,a}e^{-aL/2}\longrightarrow0.
\]

one/two-field cross 的绝对尾由相应加权 Hilbert--Schmidt Cauchy--Schwarz 控制。有限线性组合也成立。这覆盖本阶段所需的联合 Gram 与多项式电荷插入，不宣称任意增长次数、增长 t 或无限时间级数已在这里交换极限。

## 5. 为什么系数恰为 4lambda iX

令 X~Normal(0,1/(4lambda))。无 D 内部源时，charge MGF 因子与 energy/neutral 固定图分离，给相应共同矩中的 Gaussian charge 因子。含 D 内部源时，它们不再独立。

取 I(z)=i_1z_1+i_2z_2，则严格 Gaussian 恒等式为

\[
\mathbb E_Xe^{tX+4\lambda X I(z)}
=e^{t^2/(8\lambda)}
\exp_{\rm nil}\{tI(z)+4\lambda i_1i_2z_1z_2\}.
\]

线性 t i_a 和 DD connected 项 4lambda i_1i_2 同时被同一个随机变量重现。因此逐内部 matching 的 charge-linear 部分是 **4lambda i_a X**，不是 2lambda i_a X，也不是一个新的独立 Gaussian。

从 full kappa12 中扣去 4i1i2 后，其余桥弦几何与该 Gaussian charge 表示相乘。对一个 D 与一个 neutral L，i_L=0，不产生额外 charge covariance。这个表示确定了 joint regression 必须保留的电荷项；一般 energy regression、残量谱及 exact pinching 仍需另作证明，不能由本节的 mixed coefficients 单独推出。

## 有限验证接口

配套新 driver 仅写 `results/susskind-memory-joint-regression/charge-tilt-checks.json`，拒绝覆盖。有限 CAR 使用完整 charge 生成多项式而非抽样 charge；每一 u^m、每一 source 子集的系数与独立单站点乘积比较。这一恒等式同时支持有限阶 X 插入与任意有限 tilted trace。另记录 rational occupation 的精确 singleton 化简、复 t 的有限 p 修正，以及 G2 零-chaos 的全 charge 高精度倾斜和。实际计数、参数、依赖及源码哈希以科学 JSON 为准。

首次运行 22.34 秒通过，exclusive 保存科学输出。三个完整有限系统 (N,p)=(2,1),(3,1),(4,2) 的全部 p-set assignments，给 19502 个 pairing/layout 案例及 315192 个逐 charge/source 系数相等，其中 21258 个非零、6929 个负系数均保留；原始整矩阵最大绝对值 256。N=3 仅用于有限 CAR 恒等式，并未纳入偶 N 的双标度序列。18 个代表性 charge 多项式另保存 X^0 至 X^4 的精确矩。

另有 16 个 rational singleton 检验、16 个 tilted collision ratios、216 个复 t scalar 修正案例，以及 20 个 G2 全 charge 倾斜和。70/110 位复核的最大绝对差分别约 3.01e-65 与 3.00e-67。有限 N 与极限的实际偏差仍保存，例如 N=65536、p=128、lambda_N=1/4、t=-1/2 时，归一二阶倾斜矩的绝对偏差约 0.05849；不会用数值精度冒充渐近精度。全 charge 求和只略过数学上严格为零的 m<p 项，没有统计或 charge 尾截断。
