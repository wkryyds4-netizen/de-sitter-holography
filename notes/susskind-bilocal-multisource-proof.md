# 固定次数多源：有限 N 矩界、Poisson 极限与 W1 接口

2026-11-07，新 bilocal-four-point 阶段。模型为标准 CAR、同一原始复 Hermitian Gaussian 耦合的 overlapping p-body Hamiltonian，beta=mu=0 归一 annealed full-Fock 迹；hbar=1，string 单位。偶 N、p→∞，lambda_N=p²/N→lambda>0，lambda_N 始终留在正紧区间，Js 固定。本稿证明固定次数字段和任意固定有限源数 r 的矩界与逐图极限；新有限显式检查覆盖 r≤4。不主张完整实际有限时间 bilocal 的四点极限。

本稿沿用 [whole-field 源定义](susskind-memory-oriented-transfer-proof-draft.md) 与 [local-size 桥](susskind-memory-size-bridge.md)，把之前两源之后未用的有限多源界写成真实证明。所有旧文件、旧科学 JSON 和 mains 保持不动。旧 [J2/H 任意固定联合矩的 mixed-q 结果](susskind-singlet-multimark.md) 作为前驱，不重复算成新发现。

## 1. 字段、禁配对区与所证范围

令 Z_d=W_d(H^d) 为 coupling Gaussian 的纯 Wick 字段；它在 Gaussian 图中删除整个 d 槽区内部的配对。这里 W_d 不是 q-Hermite 谱多项式，也不是下面的新字段 W1。允许每个字段携带一个 mark：

- whole-field D mark 为 (pd−D)Z_d，其中 DY=QY−Σ_i c_i†Yc_i；
- local-size L mark 作用于某个连续 k 槽子区，源为 pk−L，整个所在 Wick 区仍禁止内部配对；L 是 half-Clifford-size map，而非 U(N) rank。

不同字段的 mark 区互不交叠。可插入任意固定个未标记 H 槽；也可在字段中加有限线性组合的未标记 q-Hermite/能量多项式反项，只须其系数在参数紧集上有界。展开反项后，按其真实 ordinary 或 Wick 类型分别施加禁配对规则，不替未标记 ordinary 能量块虚构内部禁配对。

固定总槽数 2ell 和 r 个 mark。每个允许的 Gaussian pairing pi 将槽两两匹配，给 2ell 个独立均匀 p-sets I_u；槽 u 的词为 C†_(I_u) C_(I_pi(u))。所有场与 adjoint 的原始次序保留。每个 source 只出现一次，故用
\[
 \mathfrak A_r=\mathbb R[z_1,\ldots,z_r]/(z_1^2,\ldots,z_r^2).
\]
记各 mark 区长度 k_a，sigma_a 为 D 或 L。有限源是 exp_nil[Σ_a z_a(pk_a−sigma_a)]，其全标记系数恰为指定逐字段乘积，不是一般 real-source 指数的渐近替代。

所需结论：每个固定 trace word 的所有 source 系数在有限 N 上一致有界，并逐图收敛；因此其有限线性组合及由重复字段构成的固定 Schatten 偶矩一致有界。下述界实际对每个固定 r 都成立；新科学检查及本阶段直接使用先限于 r≤4。

## 2. 精确局域源与 singleton 相消

按 site 重排 CAR 字母，保留完整费米子符号。对于单 site 的矩阵 X=[[a,b],[c,d]]，两个局域 loss 为
\[
 D_iX=\begin{pmatrix}0&b/2\\c/2&d-a\end{pmatrix},\qquad
 L_iX=\frac12\begin{pmatrix}a-d&b\\c&d-a\end{pmatrix}. \tag{1}
\]
前者是 {n,X}/2−c†Xc；对全局 neutral 词求和即上述 D。后者是 [X−c†Xc−cXc†]/2，求和即 L。先将一个 mark 区在该 site 上的完整词相乘，再作用 1−z_a sigma_(a,i)，最后按字段原次序取 Tr/2，得到 T_S(z)，S 是此 site 所属的 p-set 标签集合。即使零源局部 trace 为零，也保留整个有符号多项式。

记 singleton t_u(z)=T_{u}(z)。因 mark 区不含同一 Gaussian 弦的两个端点，singleton 在该区若出现，局部词只能为一个奇字母；(1) 对两种 loss 都给 sigma(c)=c/2、sigma(c†)=c†/2。于是
\[
 t_u(z)=\frac12\prod_{a:\,u\text{ 在 mark }a\text{ 恰有一个 incidence}}
                           (1-z_a/2).
\]
每个 k_a 槽区恰有 2k_a 个这样的标签 incidence。每个 I_u 含 p 个 site，因此 nilpotent 环中精确地
\[
 e_{\rm nil}^{\sum_a pk_az_a}\prod_{u=1}^{2\ell}(2t_u(z))^p=1. \tag{2}
\]
这同时消去所有混合 source 的 O(p) 项；不是只核对各一阶导数后忽略交叉项。若带内部 Gaussian 配对，singleton 规则需另外处理，本稿的纯 Wick 限定避免引入它们。

对至少两个标签相撞的 site，将其局部多项式除以所属标签各自的带源 singleton：
\[
 R_S(z)=\epsilon_S\,T_S(z)\Big/\prod_{u\in S}t_u(z). \tag{3}
\]
epsilon_S 是 CAR 重排符号。分母零阶是 2^(−|S|)，总能在有限 nilpotent 环中求逆；没有除以可能为零的 T_S(0)。重排符号按字母交换的 parity 分解到各碰撞 site；无交叠的基准符号为正。对双交叠，epsilon_{uv} 可直接由将标签 v 合并到 u 后的完整序列 (x,pi(x)) 的 inversion parity 计算。三重及更高交叠在有限 N 的恒等式中全部保留。

## 3. 任意固定 mark 数的有限 N 支配界

取系数范数 ||Σ_S a_S z_S||_A=Σ_S|a_S|，它对 nilpotent 乘法次可乘。固定槽数、r、mark 布局和 pairing 后，(3) 只有有限种非空局域类型；所有 ||R_S||_A 都被某个 B≥1 控制，且与 p,N 无关。对固定全部布局也可取最大值。

令
\[
 R_{\rm tot}=\sum_{u<v}|I_u\cap I_v|.
\]
含至少两个标签的 site 数至多 R_tot。由 (2)--(3)，除去未标记 singleton 2^(−2ell p) 后，完整源词的系数范数至多 B^(R_tot)。因此无需先删三重交叠即可证明支配。

对 Y=|I_u∩I_v|，精确 factorial moments 给
\[
 \mathbb E\binom Yk=\frac{(p)_k^2}{(N)_k k!}
 \le\frac{\lambda_N^k}{k!}.
\]
故对 a≥0，
\[
 \mathbb Ee^{aY}\le\exp[\lambda_+(e^a-1)].
\]
若 M=binom(2ell,2)>0，Hölder 不要求各 overlap 独立，直接得到
\[
 \mathbb E B^{R_{\rm tot}}
 \le\exp[\lambda_+(B^M-1)]. \tag{4}
\]
每条 Gaussian 弦的确定 prefactor 为
\[
 b_N=w_N\binom Np^2 2^{-2p}
 =\frac{J_s^2}{\lambda_N}\left[\frac{(N)_p}{N^p}\right]^2
 \le J_s^2/\lambda_-.
\]
有限 pairing 数、(4) 和此 prefactor 给每个 source 系数的有限 N 绝对一致界。将 (P†P)^k 展开时仅增加固定槽数与固定 mark 数，故得到各固定 Schatten 2k 矩界。反项和未标记 H 幂只增加有限项，结论保持。它不是从 annealed L² 推出 matrix Schatten 高矩的 hypercontractivity，也不把 rank mark 自动当作已证的 local-size mark。

## 4. 保留全部 source 的 pair-Poisson 极限

对固定标签数，各 Y_uv=|I_u∩I_v| 联合趋独立 Poisson(lambda)。一个直接的矩证明是：对固定多重 factorial moment，为每次所要求的碰撞选互异 site，得到主项
\[
 (N)_K\prod_u\frac{(p)_{d_u}}{(N)_{d_u}}
 \longrightarrow\lambda^K,\qquad \sum_ud_u=2K.
\]
任何不同 pair 要求合并到同一 site 的修正至少含三个标签，阶数趋零；固定阶只有有限种这种布局。(4) 在任意更大常数 B 上成立，保证 moment 确定性与所需一致可积。

更直接地，存在任意三重交叠的概率至多
\[
 \binom{2\ell}{3}\frac{p^3}{N^2}\longrightarrow0.
\]
其加权贡献也消失，因为 Cauchy--Schwarz 和 (4) 给
\(\mathbb E[B^{R_{tot}}1_{triple}]\le
(\mathbb EB^{2R_{tot}})^{1/2}\Pr(triple)^{1/2}\to0\)。所以三重 site 被删去发生在建立有限 N 支配之后。

在无 triple 事件上，精确源归一词为 ∏_(u<v) R_uv(z)^(Y_uv)。联合 Poisson 与一致可积允许逐 nilpotent 系数取期望；即使某些 R_uv(0)=0 或为负也无问题。有限维代数的 Poisson 生成函数给
\[
 \mathbb E\prod_{u<v}R_{uv}^{Y_{uv}}
 \longrightarrow\exp\left[\lambda\sum_{u<v}(R_{uv}(z)-1)\right].
\]
又 b_N→C exp(−lambda)，C=Js²/lambda。因此每个完整 diagram 的极限为
\[
 \boxed{C^\ell\exp[\lambda\kappa_\pi(z)],\qquad
 \kappa_\pi(z)=\sum_{u<v}(R_{uv}(z)-1)-\ell.} \tag{5}
\]
零源是既有 unmarked CAR 计数 κ_pi(0)=−4 cr(pi)，令 q=exp(−4lambda)。此处 (5) 首先作为完整 A_r 多项式指数使用，尚未截断任何三/四源项。Hermite 反项可在此后作有限线性组合；对同一标记字段减 cP 等价于将源乘 exp_nil(−cz)，不能把反项改成 disorder 独立的探针。

## 5. 在本限定 D/L 类中，高阶 κ 为什么为零

此性质来自额外局域论证，不来自“pair-Poisson 就是二阶源”。一个双标签碰撞有四个 incidence。对任意 mark 区：

- 零 incidence 时，局部词是 I，D_iI=L_iI=0，该源不出现；
- 一个 incidence 时，局部词是 c 或 c†，其因子 (1−z_a/2) 与带源 singleton 分母精确相消，包括其他源的全部系数；
- 所以一个 source 若在 R_uv 中存活，该 mark 区至少消耗两个 incidence。

mark 区互不相交，四个 incidence 因而至多支持两个存活 source。故每个 R_uv 以及 κ_pi 的三阶及更高 source 系数在这里均严格为零。这一证明对任意固定 r 成立，但不延伸到互相重叠的 mark 区、同一字段多次 loss 导数或别种局域超算符。

尽管如此，完整四源 moment 不为零。去掉 q^(cr) 后，全标记系数必须用
\[
 [z_1z_2z_3z_4]e^{\lambda\sum_{\varnothing\ne A}\kappa_Az_A}
 =\sum_{\mathcal P\in\operatorname{Part}([4])}
          \lambda^{|\mathcal P|}\prod_{A\in\mathcal P}\kappa_A. \tag{6}
\]
在本类中仅含 singleton/pair blocks 的 partitions 存活。这里的 κ 是固定 Gaussian diagram 的局部碰撞指数；不是 ensemble 四点的普通 cumulant。因此高阶 κ 消失不意味着这些算符是普通 Gaussian 联合场。

## 6. 单槽 L 字段 W1 与 J2 的新接口

定义有限 N 的 local-size 代表
\[
 W_{1,N}=(p-L-\lambda_N)H_N.
\]
若采用 transfer 稿的固定极限参数 lambda 反项，两种代表的差仅为 (lambda−lambda_N)H_N。所有 sector 的 E H_m²=C_m I、C_m 一致有界及 Gaussian 矩估计使 ||H_N||_r 在每个固定 r 上一致有界，故该差在每个固定 Lr 趋零；以下极限与此选择无关。
对一个单槽 L source，该槽的两个标签是其 Gaussian 弦的 own I/J。与其他标签碰撞时，该 source 在 ratio 中是 odd singleton 因子，精确消去；只有 own I/J 相撞会留下 source。

这一个 own collision 的局部矩阵是 n。若弦另一端无 mark，source ratio 为 2+z_W；若另一端也为单槽 L，则 ratio 为 2+z_W+z_{W'}+z_Wz_{W'}；若另一端处于 whole-D mark 中，ratio 对 D source 无依赖，仍为 2+z_W。最后一句直接用 D_i n=n，而 L_i n=n−I/2；不要求 D field 只有两个槽。

于是 κ_W=1，两个 W source 之间 κ_WW'=1 当且仅当它们是原 Gaussian pairing 的同一弦两端，其他情况为零；κ_WD=0，且第5节排除更高 source。乘 exp_nil(−lambda_N z_W) 后，单源项在极限精确消失。因此任何非零固定图中，每个 W 槽必须与另一个 W 槽直接 Gaussian 配对，每对贡献 lambda。特别地
\[
 \|W_1\|_2^2=C\lambda,\qquad \langle HW_1\rangle=0.
\]

现在接上前驱已经证明的 J2/H 端点交换相消。既然 surviving W 弦不接 J2 槽，它们只是额外的未标记 spectator 弦。交换一个连续两槽 J2 字段的两个端点不改变与任何这种 spectator 弦的总交叉数，所以旧相消逐个保持。剩余 J2 依旧两两形成长度2的双弦 bundle，H/W 各为长度1弦。令
\[
 w=4C^2\lambda+\frac{16C^2\lambda^2q}{1+q},\qquad
 \widehat H=H/\sqrt C,\quad \widehat J=J_2/\sqrt w,
 \quad\widehat W=W_1/\sqrt{C\lambda}.
\]
在前驱 J2 纯源代表的有限次数矩识别下，三字母固定 word 的 surviving 配对只连同种字母，交叉权为 q^(ell_a ell_b)，其中 ell_H=ell_W=1、ell_J=2。因此 HH、HW、WW 交叉均为 q，HJ、WJ 为 q²，JJ 为 q⁴。这个推论的新增内容是同一微观耦合下 W1 的接入；J2/H 的压缩本身是旧结果。

将这些代表换成实际 rank 字段或 exact-pinched 字段，需要另外的高范数强替换；本稿提供 local 字段所需的全部固定高矩，不将已知 L² 桥直接当成四点乘法连续性。完整 actual 有限时间 bilocal 四点还需要时间重建余项的相应高范数控制，未由此宣布完成。

## 7. 本轮独立科学检查

新 code/check_bilocal_four_sources.py 首次运行通过，17.19s；results/susskind-bilocal-four-point/finite-four-source-checks.json 使用 exclusive-create 封存。旧 mains 未运行。

- N2p1、N3p1 的显式完整 Fock CAR，用 integer matrices 直接施加全局 D 和 L，对照独立 Fraction 局部超算符；四个二槽 pure-Wick fields 的全部60合法pairing、4种D/L布局与指定有限p-set赋值，加四单槽LLLL的全部N⁴赋值及W1反项，共40032系数比较，全部16 masks保留。
- 新局部引擎以 anticommutator/creation-annihilation 定义实现，逐局域输出对照封存 helper；298个 diagram/type 案例给4768个κ系数比较，另以 set partitions 独立核对4768个指数系数。1490个三/四阶κ为零，216例完整四mark系数非零。
- 混合1/2/1/2槽 LDLD 布局独立核40个 W-D cross 为零、10个 W-W own-chord 规则。有限CAR中保留N3二槽的68负项和W1的992负项；所记录最终矩阵乘积的最大绝对整数为81，所有乘法先检查int64界。

精确参数、Python/NumPy版本和所有输入 SHA256 均在 JSON。输出 SHA256 为 f79e9ab2b810f74549a43fc97e409d47ebf7f8b3c72c57f3d91f1999de17c92c。本检查证实有限代数与指定图的系数，任意固定次数的有限N支配和Poisson极限依靠第2--5节证明。
