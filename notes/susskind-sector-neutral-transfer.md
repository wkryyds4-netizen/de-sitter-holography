# 固定电荷的实际中性记忆场

2026-11-08，本地阶段。正式推导见 [LaTeX](../manuscript/susskind-sector-neutral-transfer/main.tex) 与 [八页 PDF](../output/pdf/de-sitter-susskind-sector-neutral-transfer.pdf)。沿用 standard CAR、overlapping complex Gaussian p-body、同一静态耦合、归一 annealed sector Schatten、lambda_N=p²/N→固定正lambda、固定能量规范、hbar=1。强替换在 |m−N/2|≤L0 p；选态和下述中心矩结论在 m=N/2。H 是 string 能量，h=H/√C 才是无量纲能量。

## 1. 接上真正缺失的 composite filter

令 Dfrak=p Ncal−Rcal，Psi=Phi/C=mu(Rcal)，Theta=mu(p Ncal)，E_j=Z_j−Q_j(H)，eta_j=pE_j。既有 rank/chaos 的联合谱分解和 mu 的 12p/N Lipschitz 界给

\[
\|(\Psi-\Theta)\eta_j\|_2
\le12\lambda_N\|\mathfrak D E_j\|_2
\le12\lambda_N(\|E_j\|_2\|\mathfrak D^2E_j\|_2)^{1/2}
=O_j(p^{-1/2}).
\]

第一输入 E_j=O(1/p) 可直接用已有实际 eta 包络。第二输入来自旧 neutral-marked-series 的 fixed-degree graph 证明：中央窗的 Racah 六矩证书在 a/√N、r/√N、(m−N/2)/√N 紧集有效，覆盖固定 k 的 r≤kp；Gaussian sym 收缩给 Dfrak²Z_j 有界，一步四矩给 double commutator，递推得 Dfrak²H^n 与 Dfrak²Q_j 有界。上述全在同一固定 sector 的归一 HS 范数中成立。

这比重新拼接整个 rank-proxy 时间级数更直接，但仍使用已有 rank 图范数输入；不能称为完全不用 rank 的证明。

## 2. 固定次数与实际全时间级数

上一阶段 A_v→4lambda_N C_m theta_v S_v，S_v=(vp−L)Z_v，固定 v 的所有固定 S_r 误差 O(1/p)。把它代入实际受迫 Jacobi recurrence，top-chaos projection 与固定次数 Holder 给 eta_j−eta_j^S=O_j,r(1/p)。

T_j=A_j+(Ctheta_j−Phi)eta_j，V_j^S=4lambda_N theta_j S_j+(theta_j−Theta)eta_j^S，精确差为

\[
C^{-1}T_j-V_j^S=C^{-1}A_j-4\lambda_N\theta_jS_j
+(\theta_j-\Theta)(\eta_j-\eta_j^S)-(\Psi-\Theta)\eta_j.
\]

故固定 j 的 L2 差 O(p^-1/2)，高矩插值给所有固定有限 S_r 收敛，不把同一速率扩大到所有 r。实际与 local 的系数均有逐 sector 的 polynomial(j)L_r^j√j! 包络，Jacobi c_j(t) 有 K_S L_S^j/j!，所以先 fixed J 后去 J，得到 compact real time 的实际 B_N 与 local proxy 强差在每个 S_r 趋零。保留有限 alpha_m=1−p/m、epsilon_m=1−alpha_m² 和 prefactor m epsilon_m C_m/p。

## 3. 中央矩、exact pinching 与选态

每个 local 后代为 W_d[H^a(S_v diamond Z_w)]；L mark 落在纯 Wick 子块内，外部 H 与整体 chaos filter 在 higher-wick 源定理范围内。其 all-G source 与旧 tracial source 的逐图恒等式直接给 central 固定混合矩，不能由 fulltrace 平均反推。

继承旧回归 v_j(H)，N_j=V_j−v_j(H)。无 field bridge 的图按组抵消；normalized energy Hermite 系数有有限带宽和 polynomial×q^min 衰减，给 L²(rho×rho) 密度、对角质量零。对实际有限 N centered local 场的正两能量测度，用 central 高矩、弱收敛及 Portmanteau 得 ||P_H N_j,N||2→0，再插值至所有固定 S_r。投影包括完整简并块，严格先 finiteN exactpinch 后 DS，不交换长时间平均。

于是完整中性极限场为

\[
B_\infty(t)=2C e^{-iHt/2}\sum_{j\ge1}c_j^{(q,C)}(t)N_j e^{-iHt/2}.
\]

任意固定有限 S_r 级数尾由 finiteN pinched V_j 的包络传入极限；无需预设整个回归映射的高范数收缩。可与已控制 charged time fields 形成固定长度协变词，再由有限 N orbit density 传到 P†(covariant word)Q 的固定多项式制备。

单 selected 态继续为

\[
\sigma_1(\mathcal B)=[1-s(1+\lambda)]\tau(\mathcal B)
+s\tau((h-W)\mathcal B(h-W)),\quad
s=\sqrt q,\quad W=(S_1-\lambda H)/\sqrt C.
\]

中性词不能用 D-only 的 1−s+s h²。全态不默认 stationary，所有外部相位照原序保留。重新得到旧初始斜率校准：

\[
\sigma_1(B_\infty'(0)^\dagger B_\infty'(0))
=64\lambda^3 C^3q^2[1+\lambda\sqrt q(1+q)].
\]

后面的正修正正是漏掉 W 边界项时会丢失的项。这是旧斜率结果在完整场连接中的校准，不冒称新斜率发现。

## 4. 检查与下一步

新 interface-checks.json 首次成功：三有限 CAR 系统保留全部 Gaussian 坐标，15 个 degree 行核 composite filter 和 graph Cauchy–Schwarz，6 行 filter 非零；直接 Kraus 与 Casimir 谱作用最大误差9.65e−16，浮点容差检查。两形式有理参数点核16个旧两能量图比较、40个回归零式、16个 selected density 等式、2个 D-only 错用的非零差。参数 q=gamma^4、lambda 独立，不当作物理 q=e^-4lambda 的数值例。没有用有限N数值pinching替代谱论证明。科学SHA 6f10a97e7b014b0e267dbdcb2c5686bd758849d9d535bf34ca840050edc69635；成功driver不重跑。

下一首任务：原始 total memory 到 charged+neutral 的逐 sector 小误差，特别是 collective-return 守恒背景及投影；再将 Volterra 返回接原 bilocal。已有 force 高矩是上界，不直接证明小误差。多个 nonpolynomial selected 因子与 pinched 块任意交错亦不能由 fixed-preparation 论证自动得到。完整 selected bilocal、空间/质量体字典与finite-G仍开放。
