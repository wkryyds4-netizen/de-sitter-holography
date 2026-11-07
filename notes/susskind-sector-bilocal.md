# 固定电荷实际双局域与选态制备

2026-11-08。完整证明见 [LaTeX](../manuscript/susskind-sector-bilocal/main.tex) 与 [PDF](../output/pdf/de-sitter-susskind-sector-bilocal.pdf)。沿用静态 overlapping complex Gaussian p-body CAR 模型、同一耦合、归一 annealed sector Schatten 范数、string 能量/时间和 hbar=1。lambda_N=p²/N→lambda>0。所有范数阶与时间窗固定；统一界覆盖 |m−N/2|≤L0 p，极限矩结论限 m=N/2。

## 这一步接通什么

前一阶段已经控制总记忆到普通传播项的返回。本阶段先得到实际总记忆的所有固定有限 Schatten 界，再证明半填充下的四辅助场条件回归，最后将实际双局域重建误差改进到 O(p^-1)。由已有轨道密度界，有限个协变双局域组成的完整块可放在已控制的固定多项式制备 P、Q 之间。任意 selected 因子穿插在非多项式块内部仍不在此结论中。

## 1. 有限 alpha 的精确接口

令 h=H_m，alpha=1−p/m，epsilon=1−alpha²，F=V†exp(isK)V，D 为跨相邻电荷族的 body loss。则

\[
E_m=m(F-\alpha^{\mathcal N}e^{ish})
=[m(1-\alpha^{\mathcal N})-D]e^{isH}|_m,
\qquad \mathsf G_m=[p\mathcal N-D-\lambda_N\mathcal N(\mathcal N-1)]e^{isH}|_m.
\]

微分算符 Qalpha[X]=−X''+i alpha{h,X'}+alpha²hXh 满足 Lambda=Qalpha[F]，以及

\[
m(\Lambda-\Lambda_{OU})=\mathscr Q_\alpha E_m,
\qquad \Lambda_{OU}=\epsilon E_u[u e^{is(\alpha h+\sqrt\epsilon u)}u].
\]

外侧和指数内部的 u 是同一辅助 Gaussian 拷贝。标量 chaos 余项 r_m(k)=m(1−alpha^k)−pk+lambda_N k(k−1) 满足全整数 k 的三次多项式界 a_m(k+1)^3，a_m=O(1/p)。旧 OU/Chebyshev 乘子和新的 ordinary 时间尾允许固定导数与外侧能量幂，故 Qalpha E−Q1 G=O_r(1/p)。不从无权收敛直接求导。

## 2. 先建立实际总记忆界

偶数辅助场展开为 LambdaOU=epsilon Phi(U)−epsilon² I_s^alpha+R_s，R_s 的 r 范数至多 epsilon³|s|⁴||u||_(6r)^6/24。四 u 的三种配对全部保留。交叉配对按 B_i X B_j 与 Y B_i Z B_j 分解，行列平方和均为 C_m²I，因此 ||I_s^alpha||∞≤3C_m²s²/2。

带电 A_N 由上节受控，中性 B_N 由上一阶段实际 Jacobi 尾受控；m epsilon²=O(1)、m epsilon³=O(1/p)，再加旧返回误差 O_r(p^-1/2)，即得实际 Fscr_N 的所有固定有限范数一致有界。此步不依赖尚未证明的条件回归。

## 3. 半填充的条件回归

central-ordinary 的完整 Fourier 支配取所有 marks 为空：每个允许 pairing 的归一极限为 q^cross，q=exp(−4lambda)。独立颜色只删去不兼容配对，不改变单图占据迹，因此直接在 m=N/2 得标准多色 q-Gaussian 矩，未从 full-Fock 平均反推中央扇区。

BKS Theorems 2.11、2.13 的第二量子化条件期望给 E_(h,v)W(h,u)=E_hW(h,u)。对两个完整独立辅助块取 replica，将 T_N=E_u[u f1(h)u f2(h)u f3(h)u] 的 TT、Tt、tt 三项极限抵消，得到 T_N−t(h)→0 于 L2。

I_s 的 total-degree 时间级数有可和层界 S^(n+2) sqrt((2n+7)!!) C_*^((n+4)/2)/(2n!)。限制到极限能谱再 clip 的连续守恒函数有界，故 L2 收敛插值至所有固定 r。有限 alpha 的相位 telescope 仅需 5r 矩，误差 O(1/p)。exact pinching 消去守恒函数，于是 Fscr_N−A_N−B_N→0。此条件回归不提供速率。

## 4. 共同局域极限与实际传播

物理能量规范的 J_d=C^(d/2)j_d，H^n=sum_d a_nd^(q,C) Q_d^(q,C)(H)。设 R_n=sum_d a_nd J_d，L_n=R_(n+2)−{H,R_(n+1)}+H R_n H。则

\[
\mathcal F_\infty(s)=e^{-iHs/2}\left[
\sum_{n\ge0}\frac{(is)^n}{n!}L_n+
2C\sum_{j\ge1}c_j^{(q,C)}(s)N_j\right]e^{-iHs/2}.
\]

先固定截断再取 N→∞，随后以实际/local 尾去截断。整个 charged/neutral 塔保留同一耦合的混合矩。

M0=P_H M 的静态 bath dilation 给收缩 F0，且保留整块简并子空间。精确 Volterra 恒等式为 F−F0=−F0*(M−M0)*F。G=R†R=O_r(1/p) 给 F、F0 到 exp(iht) 的误差 O_r(1/p)；实际 Q=m(M−M0) 已在 2r 范数有界。因此 Hölder 直接给

\[
\left\|\mathcal B_N(t)+\int_0^t ds
\int_{-(t-s)/2}^{(t-s)/2}du\,
e^{ihu}\mathcal F_N(s)e^{-ihu}\right\|_{r,m}=O(p^{-1}),\quad t\ge0.
\]

负时间由伴随定义。积分负号、面积、du=dv 均由实际时间 simplex 导出。F0 在简并块内可为非标量，不可用选定能量本征基的逐对角替换。

## 5. 制备与证据

固定多项式制备的 orbit density 在 L2 一致有界；整个协变块误差由高阶 Hölder/telescope 压入 L2，然后用密度及 polarization 转成 P†(bilocal word)Q 的矩阵元。任意固定次序、伴随和紧外部中心时间可用，包括固定阶 OTOC。归一化还要求极限制备范数非零。

单 selected 态仍含 W 边界：sigma1(Y)=[1−sqrt(q)(1+lambda)]tau(Y)+sqrt(q)tau((H/sqrt(C)−W)Y(H/sqrt(C)−W))。含中性场时不能改用 charged-only 能量密度。

新输出 interface-checks.json 首跑通过后封存：124 配对源、1640 多色词、400 整块 replica、105 有理 number 余项、10 CAR/OU 行、6 实际 bilocal 时间和18范数不等式。科学输出 SHA256 为 0442b1715b4d80fb55703aaf6f474773e98f7708765e2bb4cd6e2ef671ea80f3。最大 charged/OU 相对差 3.75e-15 以下，Volterra 残差 6.94e-15 以下；这些浮点值为有限代数诊断，不是 DS 速率拟合或区间认证。C_m=1 指 CAR fixture；独立简并 bath fixture 不是 SYK ensemble。

下一步把已控制的实际双局域共同极限接到静态观察者的具体谱解释，先核已有 Susskind/DSSYK 空间、质量、singlet 字典的适用态及能量规范。任意非多项式 selected 交错和有限 G 动力学仍开放。
