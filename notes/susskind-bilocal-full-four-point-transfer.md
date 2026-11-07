# 完整有限时间 bilocal 四点：最后的接线审计

2026-11-07，本地附录。沿用 common original Gaussian couplings、standard CAR、overlapping complex p-body Hamiltonian、β=μ=0 normalized annealed full-Fock trace、固定 J_s 与 string 实时间、λ_N=p²/N→λ>0 的正紧区间规范。结论针对固定有限时间窗的 annealed ordered moments。

**审计结论：** 新 neutral/charged all-fixed-S_r 尾，加已有同样本 L² 桥、pinching 无对角质量及 bilocal 强返回，足以构造实际完整有限时间四点。多源定理并不直接覆盖带内部 Gaussian pair 的 ordinary charged source；该部分通过已有同样本强归约到所有 J_d，再由新高范数界升阶。一般后代必须使用完整图规则，不能压缩成低阶 H/W_1/J_2 三 species 模型。

## 1. 多源定理确实覆盖所有需要的中性后代

[多源证明](susskind-bilocal-multisource-proof.md) §1--5 允许一个 mark 作用在某个 pure Gaussian-Wick 块内的连续子区，并允许任意额外 ordinary H 槽。中性后代的通用项恰为
\[
 F_{a,v,w;d}=\mathcal W_d\!\left[H^a
       \bigl((vp-L)Z_v\diamond Z_w\bigr)\right].          \tag{1}
\]
整个 v+w 槽的 Wick 块禁止内部配对；L mark 只在其中前 v 槽，伴随后成为后 v 槽。H^a 槽在这个块以外，可以相互配对，也可以与块内槽配对。因此每个允许图的 marked subblock 都不含完整 Gaussian 弦，正好满足多源 singleton 相消的前提。

对 composite field 整体施加 \(\mathcal W_d\)，在 Gaussian 配对展开中只保留其内部弦数 \(i=(a+v+w-d)/2\) 的图；\(q^{\mathcal N}\) 只乘 \(q^d\)。这些筛选不改变原 v+w Wick 块的禁配对规则，也不把 larger composite 的 open set 误用作 mark incident set。任意有限个 (1)、其伴随、pure charged J_d 和能量多项式相乘时，各 mark 区仍互不交叠。

所以固定 joint word 的所有矩由原多源公式给出：逐项展开确定性反项与后代，在每张允许 pairing π 上按原始次序计算
\[
 C^{\ell}[z_1\cdots z_k]
 \exp\{\lambda\kappa_\pi(z)\},\qquad
 \kappa_\pi=\sum_{u<v}(R_{uv}(z)-1)-\ell.                \tag{2}
\]
每个 primitive 只有一个 mark，k 为该项的实际 marked-field 数；无 mark 的反项取相应较低 source 系数。chaos/OU 筛选及 scalar 系数随后按线性相加。三/四源 κ 消失并不使完整四源 moment 消失，须保留指数的全部 set-partition 项。

有限 N 正性及循环迹传给这个共同矩泛函 τ。它确定 H、所有 J_d、所有中性残量 N_j 的共同 ordered polynomial moments；二能量谱 ν_ab 本身不能确定这些四点。

## 2. Ordinary charged 内部 pairs 用同样本强桥处理

[joint-regression proof](susskind-memory-joint-regression-proof.md) §4 已证明
\[
 R_{n,N}=G_{n,N}-g_n(X_N,H_N),\qquad
 \left\|R_{n,N}-\sum_d a_{nd}J_{d,N}\right\|_2\to0,       \tag{3}
\]
其中 X_N=(Q−N/2)/p，g_n 是其明确的联合 charge/energy 回归多项式，a_nd 是固定极限 q-Jacobi 逆系数。其证明包括 R_n--J_d 的 mixed common-copy Gram 项；不是仅比较两个边缘律。固定左右 H 幂的加权版本也已给出。

新的 [charged-tail proof](susskind-bilocal-fourth-tail-charged-proof.md) 给 G_n 的每个固定 S_r 一致界。g_n(X_N,H_N) 的所有固定高范数由 sector H 矩界和 binomial X_N 固定矩控制；J_d 的高范数由 pure marked bound 控制。因此 (3) 两侧差在每个固定更高范数一致有界。对任意固定 r>2，取 R>r，用
\[
 \|Y\|_r\le\|Y\|_2^\theta\|Y\|_R^{1-\theta},
 \quad 1/r=\theta/2+(1-\theta)/R,                         \tag{4}
\]
将 (3) 及其固定能量加权版本升为 S_r 强桥。r=4 可用 R=8。

这正是处理 ordinary charged 内部 pairs 所需的步骤。无需给纯 Wick 多源定理增加一个未证的 ordinary singleton 假设。尤其不需要假定任意 ordinary 四源 κ 只有二阶。

## 3. 中性局域桥和 exact pinching 均可升阶

中性 [marked-series](susskind-memory-neutral-marked-series.md) 与 [wick-descendants](susskind-memory-wick-descendants.md) 已给：actual B_N 到 rank proxy 的紧时间 L² 桥；固定 rank mark 到 local L mark 的同样本 L² 桥，经右 Wick creation、固定 H 乘法、chaos/OU 投影传给 (1)；两边 L² 可和尾允许解除次数截断。

新 [operator-tail proof](susskind-bilocal-fourth-tail-operator-draft.md) 给 actual B_N 与 local proxy 的所有固定高范数可和尾。故这里可直接对 actual 与 local 两端使用 (4)，不需要先证明一般 rank-map 的高范数定理，也不需要 rank proxy 的高范数界。

把 local V_j 的能量回归写为 v_j(H)，令
\[
 N_j=V_j-v_j(H).
\]
v_j 的明确递推见 wick-descendants §5。其 β=C[j]_q 一致有界，强迫项系数含 \(\sum_v q^v\kappa_v\)、\(\kappa_v=O(v^2)\)，因此普通多项式系数总和至多按 j 指数增长。结合 H 幂矩界，v_j(H_N) 继承同型 \(\operatorname{poly}(j)L_r^j\sqrt{j!}\) 包络；减回归不会破坏时间可和尾。

旧两能量谱定理已证明每个固定 J_d、N_j、R_n 的残量无对角质量，从而其 exact pinching 在 L² 趋零。pinching 对所有 S_R 收缩，结合刚得到的固定高矩与 (4)，它在每个固定 S_r 也趋零。g_n(X_N,H_N)、v_j(H_N) 始终是同样本精确守恒函数，去 pinching 前先减这些反项。

这里不交换无限时间平均与 N 极限；顺序是固定 degree、传递 exact pinching、最后用新高范数尾去 degree 截断。

## 4. 明确的共同极限场构造

令
\[
 R_n^{\infty}=\sum_d a_{nd}J_d,\qquad
 L_n=R_{n+2}^{\infty}-\{H,R_{n+1}^{\infty}\}+HR_n^{\infty}H.
\]
以 (2) 定义所有有限 joint words 的 τ。H 的指数和以下级数可由多项式截断，在继承自有限 N 的 Schatten 范数/Hölder 不等式下完成；无需猜测一般后代的标准 q-Fock 表示。

记 c_j(s) 为前驱 monic q-Hermite 指数系数。完整 memory 极限是
\[
 \mathcal F(s)=e^{-iHs/2}
 \left[\sum_{n\ge0}\frac{(is)^n}{n!}L_n
                 +2C\sum_{j\ge1}c_j(s)N_j\right]e^{-iHs/2}. \tag{5}
\]
新有限 N 尾、fixed-degree 桥及 Hölder 给 (5) 在每个固定 S_r、每个固定实时间窗绝对可和。其共同多点由 (2) 逐截断确定。

对于 δ≥0，定义
\[
 \mathcal B(\delta)=-\int_0^\delta ds
       \int_{-(\delta-s)/2}^{(\delta-s)/2}du\,
       e^{iHu}\mathcal F(s)e^{-iHu},\qquad
 \mathcal B(-\delta)=\mathcal B(\delta)^\dagger.           \tag{6}
\]
这是所有 fixed S_r 中的积分，积分算符范数至多 δ²/2。charged 部分也可写成
\[
 \mathcal B_A(\delta)=e^{-iH\delta/2}
              \sum_{d\ge2}c_d(\delta)J_d\,e^{-iH\delta/2}. \tag{7}
\]
此时 Hermite 重求和在 S_r 合法：\(\sum_d a_{nd}\le K^n\)，pure J_d 的阶乘界给双重绝对和 \(\sum_n\operatorname{poly}(n)L_{r,S}^n/\sqrt{n!}<\infty\)。因此旧 charged 滤波简化不再只限于二范数解释。

neutral 部分保留
\[
 \mathcal B_B(\delta)=-2C\sum_{j\ge1}\int_0^\delta ds\,c_j(s)
 \int_{-(\delta-s)/2}^{(\delta-s)/2}
 e^{iH(u-s/2)}N_j e^{-iH(u+s/2)}du,                       \tag{8}
\]
等价地使用前驱无重节点奇性的 divided-difference amplitude。

## 5. 从实际 finite N 对象传递四点

令 \(\mathcal F_N^{(K)}\) 为 (5) 的 finite-N local representative 截断，保留共同 H_N、J_d,N、N_j,N，先处理 actual two-component coefficients，再用上述固定-degree 强桥。operator-tail §§8--9 与 charged-tail 的能量加权结论已把实际返回和 two-component 误差升至固定 S_r。因此
\[
 \lim_{K\to\infty}\limsup_N\sup_{|s|\le S}
       \|\mathcal F_N(s)-\mathcal F_N^{(K)}(s)\|_4=0.     \tag{9}
\]
同样成立于每个固定更高阶数。有限参数 C_m、θ_j,m、c_j,m 与 α_m 的替换先在固定 degree 上进行；中央收敛和 binomial 外尾控制其差。

把 (9) 送进线性滤波，误差至多 S²/2 倍。旧 [bilocal transfer](susskind-bilocal-four-point-transfer.md) §7 已给 actual bilocal 到 actual-memory 滤波的 S_4 强差 O(p^{-1/7})。由此 actual bilocal 有统一四范数，并有对应的 local 截断尾。

固定 K 后，H_N 的时间指数用普通 Taylor 多项式逼近。固定更高阶 H 幂矩界给阶乘尾，对固定 K 及紧时间窗一致；于是 (2) 的 polynomial joint moments 传给带所有真实 energy phases 的固定截断。随后去 K。

具体地，设 \(\mathcal B_N(T,\delta)=e^{iH_NT}\mathcal B_N(\delta)e^{-iH_NT}\)，\(\epsilon_i\in\{1,\dagger\}\)。则
\[
 \boxed{\lim_N\mathbb E\operatorname{tr}
 \prod_{i=1}^4\mathcal B_N(T_i,\delta_i)^{\epsilon_i}
 =\tau\!\left(\prod_{i=1}^4
       [e^{iHT_i}\mathcal B(\delta_i)e^{-iHT_i}]^{\epsilon_i}\right).} \tag{10}
\]
乘积按给定原始顺序，包含 OTOC。右端由 (2)、(5)--(8) 和绝对可和截断定义。误差的四点 telescoping 使用
\[
 |\tau(Y_1Y_2Y_3Y_4)-\tau(Z_1Z_2Z_3Z_4)|
 \le\sum_i\|Y_i-Z_i\|_4
       \prod_{k<i}\|Y_k\|_4\prod_{k>i}\|Z_k\|_4.        \tag{11}
\]
固定参数下可在紧 T_i、δ_i 窗局部一致；这里没有全实轴 T 一致的 N 极限，也没有增长时间结论。

完成的是指定微观模型的完整有限时间共同四点构造和 actual-operator 传递。它不是仅低阶 Taylor 诊断，也不是由二点谱猜测四点；尚不增加 de Sitter 体位置/质量字典、canonical dressing 或有限 G 引力散射的结论。
