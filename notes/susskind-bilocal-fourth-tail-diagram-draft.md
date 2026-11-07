# 四范数尾的图方法：纯 Wick 的有限 N 界与 ordinary 缺口

2026-11-07，未封存工作稿。接 NEXT 中完整 bilocal 的第四范数尾任务。模型、态、单位沿用 common-J overlapping complex p-body CAR、beta=mu=0 归一 annealed full trace、string 时间；lambda_N=p²/N 留在正紧集并趋于 lambda>0，Js 固定。

本稿不改旧代码或科学 JSON，不运行旧 main。以下区分三个层级：纯 Wick/local-size 代表的有限 N 可和界；带内部 Gaussian 收缩的 ordinary charged 系数的逐次数极限界；尚未补齐的实际全时间第四范数接口。

## 1. 任意交叉配对的有界 contraction：保留符号

考虑一个固定 Gaussian pairing 的系数矩阵词。每个独立 Gaussian index alpha 在两个端点各出现一次，复 Hermitian 耦合对应 A_alpha 和 A_alpha†。假设对每个允许的中间 charge，相关受限系数族满足
\[
 \Big\|\sum_\alpha A_\alpha^\dagger A_\alpha\Big\|\le v_-,\qquad
 \Big\|\sum_\alpha A_\alpha A_\alpha^\dagger\Big\|\le v_+.
\]
从右向左扫描该矩阵词。首次见到某 index 时使用 column map
\[
 V\psi=\sum_\alpha A_\alpha\psi\otimes e_\alpha,
 \qquad \|V\|\le\sqrt{v_-};
\]
第二次出现时用 row map W(Σ_alpha psi_alpha tensor e_alpha)=Σ_alpha B_alpha psi_alpha，其范数由该端点的 row covariance 控制。未闭合 index 存在独立的 ancilla 因子中；交叉只要求把欲闭合的因子通过单位ary permutation 移到指定位置。其余系统算符、固定 site 的 c_i/c_i† 插入与这些 ancilla identity 作 tensor 乘积。

因此一个任意交叉 pairing 的完整矩阵 contraction，其 operator norm 至多为全部端点 row/column 范数之积，乘固定插入算符的范数。没有先对每个 CAR 配置取绝对值，也没有删掉负图。尤其若各弦的两个端点 covariance 都不超过 v_chord，则
\[
 \|\operatorname{Contraction}_\pi\|
 \le\prod_{\rm chords}v_{\rm chord}. \tag{1}
\]
这个有限矩阵不等式不依赖平面性；它正是避免 exp(const·degree²) 损失的关键。

## 2. 给 p-set 加有限 membership 限制的 covariance

系数 A_IJ=sqrt(w_N) C†_I C_J。在 charge m 上，未受限 row/column covariance 均为 C_m I，且 C_m≤C_*=(Js²/lambda_−)exp(2lambda_+)。

固定若干 site，要求 J 含指定集合 T_J、I 含指定集合 T_I，记大小 t_J、t_I；含交集时按两个独立标签的要求分别计数。对 occupation basis O、|O|=m，A†A 非零要求 J 属于 O；移除 J 后，I 可从 N−m+p 个 site 中选。令 K=T_J union (T_I intersect O)。若 T_J 不属于 O，column covariance 的该对角项为零；否则其精确值为
\[
 w_N\binom{m-|K|}{p-|K|}
       \binom{N-m+p-t_I}{p-t_I},
\]
其中非法二项系数置零。这个式子直接保留 T_I 中原 occupied site 也必须进入 J 的条件。放宽 K 到 T_J 后，受限和的对角项可用
\[
 C_m\frac{(p)_{t_J}}{(m)_{t_J}}
       \frac{(p)_{t_I}}{(N-m+p)_{t_I}}
\]
支配。若某个要求的 I site 原已 occupied 且未被 J 移除，该项为零，故不会破坏上界。AA† 同理交换 I/J 的角色。

在两侧中央区 N/4≤m≤3N/4，经过至多 O(r) 次固定 CAR 插入的中间 charge 仍距端点 O(N)；当 N 足够大时，固定 t=t_I+t_J≤2r 有
\[
 v_{\rm restricted}\le C_*\left(K_r\frac pN\right)^t. \tag{2}
\]
K_r 与 degree 无关。若同一 Gaussian 弦两个端点各有 membership 条件，应先取其并集作为这条弦允许的 index 集合，然后在两个端点使用同一受限族。这样每条独立 site/set 条件只算一次。

CAR sandwich 导数可展开为有限项 L A_IJ R，L/R 是固定 site 的 c 或 c† 的乘积、范数不超过1。每项保留其 membership 条件，再把固定插入作为 (1) 中的收缩算符。因此不需要给导数矩阵另猜一个 scalar covariance。

## 3. 纯 Wick mark 的 exact tagged-site 展开

对 neutral X,Y，body map 满足精确 product rule
\[
 D(XY)=(DX)Y+X(DY)-\sum_i[c_i^\dagger,X][Y,c_i]. \tag{3}
\]
每个 p-body 系数矩阵 A 有 DA=pA。故对 k 个槽的乘积
\[
 (kp-D)(A_1\cdots A_k)
 =\sum_{a<b,i} A_1\cdots[c_i^\dagger,A_a]\cdots
                         [A_b,c_i]\cdots A_k. \tag{4}
\]
省略处保持原始次序，不交换中间矩阵。[c_i†,A_a] 只在 i∈J_a 时非零，[A_b,c_i] 只在 i∈I_b 时非零；每个 commutator 再展开为两项固定 CAR sandwich。

对纯 coupling-Wick 字段 Z_k=W_k(H^k)，D 与 Wick 投影可交换。将 (4) 放入多字段期望后，所标记的两个槽 a,b 不可能被 Gaussian pairing 互相配对，因为同一 Wick 字段内部禁止配对。因此 J_a 与 I_b 是两个不同的独立 p-set 标签：每个 tag site 至少要求两个独立 set membership。

为明确 local half-size 的单槽项，记 reverse body map
\(\bar D X=(N-Q)X-\sum_i c_iXc_i^\dagger\)，则 L=(D+barD)/2。对齐次奇偶算符定义 graded derivative
\(\partial_aX=aX-(-1)^{|X|}Xa\)。直接 normal ordering 给
\[
 \bar D A_{IJ}=pA_{IJ}-\sum_i\partial_{c_i^\dagger}\partial_{c_i}A_{IJ},
 \qquad
 (p-L)A_{IJ}=\frac12\sum_i\partial_{c_i^\dagger}\partial_{c_i}A_{IJ}.
\]
双导数只在 i∈I∩J 时非零。因为 A_IJ 为偶，双导数展开为
\[
 c_i^\dagger c_i A-c_i^\dagger A c_i
       +c_i A c_i^\dagger-Ac_ic_i^\dagger,
\]
是四项范数至多1的固定 CAR sandwiches。例子 A=n_i 给 (1−L)n_i=I/2。

L 的 cross product rule 是 (3) 与 reverse 规则的平均，含
\(\tfrac12[c_i^\dagger,X][Y,c_i]\) 和
\(\tfrac12[c_i,X][Y,c_i^\dagger]\)。迭代后同样选择两个槽和一个 site；两种 orientation 分别要求 J_a/I_b 或 I_a/J_b 的 membership。加上上述 own I/J 项，每个 local-size mark 至多有 O(k²) 种槽选择与固定数目的 CAR 项，每项含两个独立 set 标签的条件。

本论证也适用于 S_v diamond Z_w 的 mark 位于前 v 个槽的情形：整个 v+w Wick 区禁止内部配对，mark 不会失去两个独立标签的条件。外部普通 H 幂可以保留，且不必禁止它们与 source 区配对。

## 4. 任意固定 mark 数的 finite-N degree 界

设 trace word 总槽数为 L，带 r 个上述纯 Wick D/local-L mark。固定一个 Gaussian pairing、每个 tag 的槽选择和 CAR sandwich 项，再按 r 个 tag 所用 site 的相等关系分组。若一组有 s 个不同 site，则至多 N^s 种 site 赋值；每个 site 至少约束两个不同 set 标签，合并条件总数 t≥2s。

将 (2) 代入 (1)，受限 covariance 提供 (K_r p/N)^t。因为 p/N≤1、t≤2r，这些因子与 N^s 合并后至多被只依赖 r 和 lambda_+ 的常数控制：N^s(p/N)^(2s)=lambda_N^s。site partitions、CAR 展开项及槽选择数合并为 K_r L^(2r)。

每个 tag 展开项至多插入两个固定 CAR 字母，所以中间 charge 偏移至多2r。总独立 membership 条件数至多2r；取例如 N≥64r，在原来的两侧中央区，所有 falling-factorial 分母均≥N/8，(2) 可取固定数值 K_r=8 的粗估计（不优化）。较小 N 只影响最终依赖 r 的常数，不影响 degree。

中央区之外仍使用同一个 tagged 展开，只忽略 membership 带来的小量。它最多留下 N^r，而全 charge 二项尾≤2exp(−N/8)；sup_N N^r exp(−N/8)<∞。少数 N 不满足中间 charge 的中央估计时，同一粗界也可吸收到固定 r 常数。由此得到真正的有限 N、任意 degree 估计
\[
 \boxed{|\mathbb E\operatorname{tr}({\rm marked\ word})|
 \le K_r(1+L)^{2r} C_*^{L/2}(L-1)!! .} \tag{5}
\]
这里 pairing 数由 (L−1)!! 支配，所有符号先保留在矩阵 contraction 内，最后才用 operator norm。若 L 为奇数，Gaussian 期望为零。

对一个 n 槽纯 source P_n 的第四矩，取 r=4、L=4n。由于
\[
 (4n-1)!!=\prod_{j=1}^n(4j-3)(4j-1)
 \le16^n(n!)^2,
\]
得到
\[
 \boxed{\sup_N\|P_n\|_4
       \le K(n+1)^2(2\sqrt{C_*})^n\sqrt{n!}.} \tag{6}
\]
更一般固定 Schatten 阶数用对应固定 r，仍是 polynomial(n)·K_r^n·sqrt(n!)。未标记能量反项的系数若按 degree 至多指数增长，也不破坏这种形式。

### 4.1 同一纯 Wick 块上的 repeated D/L deficit

这一步不依靠 nilpotent 源区互不交叠。对 parity 为 a 的齐次算符，定义 graded extension
\[
 D_gX=\tfrac12\{Q,X\}-(-1)^a\sum_i c_i^\dagger Xc_i,
 \quad\partial_bX=bX-(-1)^aXb.
\]
在 neutral 偶场上 D_g=D。对有 r_c 个 creation、r_a 个 annihilation 的 normal-ordered 系数，D_g 的本征值是 (r_c+r_a)/2；这可由 creation/annihilation 两个 Euler 计数直接验证。任意齐次 X,Y 满足精确 graded product rule
\[
 D_g(XY)=(D_gX)Y+X(D_gY)
       +(-1)^{|X|}\sum_i(\partial_{c_i^\dagger}X)(\partial_{c_i}Y).
\]
另定义 L_g=D_g−(1/2)Σ_i partial_(c_i†)partial_(c_i)；在 neutral 场上就是 L。此处 graded extension 只是有限 CAR 代数工具，不主张它在奇场上是普通 CP semigroup。

从 n 个原始 p-body 因子出发，每次 tag contraction 删去两个字母。若某项已删 t 对，剩余因子的 normal-ordered 总 degree 为 np−t。因此再作用 Delta_D=np−D_g，产生 t 倍该项和新的两因子 contraction；再作用 Delta_L=np−L_g，还允许同一因子上的 own I/J contraction。所有系数只依赖既有 contraction 数，而不含新的 p 或 N。

归纳可知：同一纯 Wick n 槽块上任意固定顺序的 k 次 Delta_D/Delta_L，是至多 k 个 tag 的有限和，项数至多 K_k(1+n)^(2k)。每个新的 contraction 仍要求两个不同底层 p-set 标签：cross-factor 情形由纯 Wick 禁内部 Gaussian pair 保证，own-factor 情形来自 I/J 两个独立标签；重复选到已经删去的同一字母只给零。旧 tag 的倍数项无需新 site。

故 (1)--(5) 继续适用，但 r 应取期望中所有 deficit 的总次数。例如
\[
 \boxed{\sup_N\|\Delta_{\sigma_1}\cdots\Delta_{\sigma_k}Z_n\|_s
 \le K_{k,s}(n+1)^{2k}L_{k,s}^n\sqrt{n!}}
\]
对任意固定偶数 s 成立。中间 parity/charge 变化仍只由固定数量的 CAR 导数控制。特别 (np−D)^2 Z_n、(np−D)(np−L)Z_n 均有真正有限 N 的可和 degree 包络。它不直接处理 Gaussian 内部 pairs，也不自动把 D 交换过带 charge 系数的 covariance CP 链。

## 5. 这给出的实际可和结果：local proxy 系列

纯 S_j=(jp−L)Z_j、H^a(S_v diamond Z_w) 直接满足 (5)--(6)。中性局域 proxy 的受迫 Jacobi 递推仅含这些 primitives；有限参数 a_j 一致有界，theta_j 几何衰减。因此将递推展开成未标记 H 幂和 primitives，确定性系数的总绝对值至多按 j 指数增长。有限 degree 的 bounded Gaussian-chaos multiplier 可用另稿的 OU/Chebyshev 估计控制，额外代价同样只是 K^j，不需假设 L²→S4 hypercontractivity。

所以局域 eta_hat_j、V_j 的第四范数具有
\[
 \sup_N\|V^{\rm loc}_{j,N}\|_4
 \le K(j+1)^b L^j\sqrt{j!}, \tag{7}
\]
其中 b 固定。旧有限 Jacobi 系数满足
\(|c_{j,m}(s)|\le e^{2\sqrt{AM}S}(2\sqrt A S)^j/j!\)。于是 local proxy 的时间级数在每个固定实时间紧窗上具有一致可和的 S4 尾。exact pinching、真实 unitary 相位及有界 charge prefactor 不破坏此结论。

这一步给的是 finite-N **local proxy** 的可和第四范数界。将完整 actual rank/covariance-filter 中性场替成这个 proxy，仍需把旧 L² 强归约或其尾提升到第四范数；不能只因每个固定 j 的局域矩已受控便宣布 actual 全级数已替换。

## 6. Ordinary charged 系数：逐 degree 的 sharper 极限界

对 G_n=[p Ncal−D−lambda_N Ncal(Ncal−1)]H^n 必须保留同一 field 的内部 Gaussian pairs。给定 pairing 时，若 field 内有 i 条完整弦，则 d=n−2i，number source 是 exp(pd z)，不是 exp(np z)。内部弦的两个 singleton 比值为 1−z 与 1+z，乘积在 z²=0 中精确为1；open 弦的 odd singleton 仍与 pd 源抵消。这个论证同时适用于四个互不交叠的 ordinary fields。

固定四源时，任何 pair-site 只有四个 incidence，局域 2×2 ratio 及其带源 singleton 分母的有限逆只有常数种类型。因此全部非零 source 系数有一个与总槽数无关的常数上界 B4。指数 kappa 是 O(L²) 个 pair ratios 的和，额外 −d(d−1) 反项仍为 O(L²)。用完整 set-partition 指数公式，不论是否先证明高阶 kappa 为零，四源系数均被 K(1+L)^8 支配，乘未标记 q^(crossings)≤1。

故固定 n 后取 DS 极限，可严格得到
\[
 \limsup_N\|G_n\|_4\le K(n+1)^2 L^n\sqrt{n!}. \tag{8}
\]
对 K_n=G_(n+2)−{H,G_(n+1)}+HG_nH，外部能量槽仍为有限普通槽，同样给 polynomial(n)·L^(n+2)·sqrt(n!) 的包络。因此已经确定的 limiting charged series 在 S4 中整函数式可和。

然而 (8) 的 sup/limsup 顺序不同于 (6)。有限 N 的 (4) 若选中了同一 ordinary field 内的一条 Gaussian 弦两端，则 J_a=I_b 是**同一个**独立 p-set 标签；tag 只需一个 membership，N·(p/N)=p，不再有 lambda 尺度。p Ncal 与 Gaussian 内部收缩必须在这个层面重新精确相消，不能沿用纯 Wick 的两个独立标签论证。

低阶例子揭示剩余项：G_2 的零 chaos 等于 −D C_Q，是有界中心涨落但非零的 charge 函数。它来自相邻 charge covariance 的差，不能在 tagged contraction 前把内部 pair 一律丢掉。可能的下一步是把这些单标签内部 tags 重组为明确的 covariance 差/charge 多项式，再让其余 tags 用 (2)；目前本稿没有证明这个 ordinary 重组在所有 degree 的统一常数。

## 7. 不能交换尾极限：明确剩余条件

逐次数 (8)、整个实际场的统一 L² 与 L² 阶乘尾，仍不足以推出 actual S4 尾。一个标量逻辑反例是取 k_N→∞、p_N=(k_N!)²、概率为 p_N^(−2) 的事件 A_N，令
\[
 f_N(t)=p_N\mathbf1_{A_N}\frac{t^{k_N}}{\sqrt{k_N!}}.
\]
其 L² 范数在任意固定时间紧窗一致有界，且具有统一的 factorial L² degree 尾；每个固定 Taylor 系数最终为零。所有固定高范数仍满足 O_S(p_N) 的粗界。但 ||f_N(1)||4=sqrt(k_N!) 发散。因此“逐 n 的 limiting 可和包络”不能替换所需的 lim_K limsup_N actual remainder。

现有以 B^(R_overlap) 取绝对值的 finite-N 误差，会产生 exp(c·degree²)；它既不证明时间尾可和，也通常不足以在需要的 degree≈log p/log log p 截断上压住 O(p) 高范数损失。新的纯 Wick contraction 界解决了其中一类代表，但完整 bilocal 仍需：

1. ordinary charged 内部 pair 的精确 tagged cancellation，或另一条实际 charged S4 尾界；
2. actual neutral rank/filter 到已可和 local proxy 的 S4 强桥与次数尾；
3. 固定实际后代的共同高点强接口，再接已经证明的 bilocal S4 返回误差。

本稿没有把 limiting S4 系列或 local proxy 的有限 N 系列等同于尚未受控的完整 actual bilocal。所有新的科学检验、正式 LaTeX 与 PDF 应在上述证明被复核并确定本阶段范围后单独生成。

## 8. D 与 normalized CP 步骤的精确交换式

为下一步处理 actual force/eta，记 neutral sector families 上的未归一化 maps
T X=Σc†Xc、S X=ΣcXc†，Q 为 charge 左乘。直接 CAR 给
\[
 [Q,T]=T,\quad[Q,S]=-S,\quad[T,S]=2Q-N,\quad D=Q-T.
\]
在远离 charge 端点处定义 t=Q^(−1)T、s=(N−Q)^(−1)S；它们分别是 normalized-sector CP contractions。对任意确定性 charge 函数 f，
\[
 [D,f(Q)]=[f(Q)-f(Q-1)]T.
\]
因此有精确恒等式
\[
 \boxed{[D,t]=t-t^2,\qquad
 [D,s]=-s+\frac Q{N-Q}ts-\frac{2Q-N}{N-Q}I.}
\]
在两侧中央区这些系数均有界。对长度不超过2p的 CP 链展开 commutator 后，[D/p,chain] 的 sup-sector Schatten 范数至多 O(1) 倍输入的 sup-sector 范数；只需留足 p+O(1) 个 charge 的中央余量。

这说明应对 covariance residual 的 **CP telescoping 表达** 逐项作用 D/p，结合 4.1 的 repeated pure bounds 与 charge 有限差分；直接粗估 [D,Phi−Ctheta] 会丢失其抵消并留下 O(p)。本段尚未完成 D 对 eta 递推中全部普通 H 乘积及系数的统一估计，因此不把 charged ordinary 尾写成已证。

## 9. 后续闭合的交叉引用

上文第6--8节保留该次证明的边界。后续 [charged-tail proof](susskind-bilocal-fourth-tail-charged-proof.md) 已用精确 CP 导数、full-Clifford Riesz、展开后的 eta 递推及 inverse-Jacobi charge 差分证明 ordinary charged 的 finite-N 全次数高范数尾；[独立审计](susskind-bilocal-charged-tail-audit.md) 给出交叉核查。其与实际中性场、同样本后代桥、exact pinching 和 bilocal 返回的完整接线见 [full-four-point transfer](susskind-bilocal-full-four-point-transfer.md)。这些后续证明解决了先前明确留下的尾接口，没有把逐次数 limiting bound 当作 finite-N uniform tail。
