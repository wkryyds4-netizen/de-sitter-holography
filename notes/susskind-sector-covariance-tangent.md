# 固定电荷协方差的一阶展开

2026-11-08。正式稿见 [LaTeX](../manuscript/susskind-sector-covariance-tangent/main.tex) 与 [PDF](../output/pdf/de-sitter-susskind-sector-covariance-tangent.pdf)。

沿用同一静态 Gaussian 耦合、standard CAR、归一 annealed sector Schatten 范数。先在完整中性块族上作用 D、bar D、L=(D+bar D)/2，再限制电荷；固定能量规范、hbar=1。lambda_N=p²/N 在正紧集，|m−N/2|≤L0 p，v、r、L0 固定。

此前 covariance 的精确路径 Psi=Phi/C 为 t^p s^p，外到内分母为 m,...,m−p+1; N−m+p,...,N−m+1。令 a_j=1−b/d_j，Delta_j 为 b−D 或 b−bar D。先对整条路径 telescope，再对每一前缀作用的 Delta_j X telescope，得到一阶项系数 product_(h≠j) a_h/d_j，余项为 i<j 的有序 Delta_i Delta_j、外接原 CP 前缀。整个有限式不除以 a_j，因此 theta=0 也成立。

0≤b≤d_min 时各参考因子在 [0,1]，前缀保持归一迹且收缩；余项范数至多 binom(2p,2)/d_min² 乘路径上的二次 deficit 最大范数。取 b=vp，一阶两类系数为 theta*zeta_T/p、theta*zeta_S/p，zeta_T 与 zeta_S 各为 p 个倒数之和再乘 p，各与 2lambda_N 相差 O(p^-1)。继承既有 pure repeated-tag 的一、二阶界，得到

\[
\|p(\Phi-C_m\theta_{v,m})Z_v
-4\lambda_N C_m\theta_{v,m}(vp-L)Z_v\|_{r,m}
=O_{v,r,L_0}(p^{-1}).
\]

保留全部有限参数，不把 lambda_N、C_m、theta 的极限收敛率擅自设为 p^-1。此处常数不对 v 一致。有限检查仅检验路径代数与 S2 余项界，一般强估计依赖正文证明及旧 pure-tag 输入。

已保存的 tangent-checks.json：96 个中性 matrix units，1824 个精确路径等式、1824 个精确归一 S2 余项不等式，608 个参考因子为零的案例，60 个有理标量行。首空乘积类型错误及脚本保留；成功计算不重跑。输出仍位于 susskind-sector-neutral-transfer 目录，因为它原是这一后续阶段的第一个输入；本短稿只封存 tangent 子阶段。

下一步：实际中性系数 T_j=A_j+(Ctheta_j−Phi)eta_j 的第二项含低 Gaussian chaos 后代。不能从 A_j 的直接桥推出 Phi 在这些复合场上的筛选作用。优先复核旧 neutral-marked-series 在 compact central window 内、尚未 binomial 求和时的 rank slope/graph tail，再接 size-bridge 的逐 sector 固定次数估计和已控制高矩。之后才是 exact pinching、完整中性时间场、实际 memory 返回和 selected bilocal。
