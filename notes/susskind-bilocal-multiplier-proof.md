# Bilocal 首返替换：正力矩与 Schatten 插值的定量证明

2026-11-07。本稿给出 [bilocal 重建证明](susskind-bilocal-reconstruction-proof.md) 的独立常数核算，只调用已封存的 [return-expansion 力高矩](susskind-memory-return-expansion.md) 和 [实际 memory 二阶界](susskind-memory-uniform-bound.md)。此新笔记初稿的四源 finite-proxy 路线被下述直接估计取代；未创建或运行四源 driver，旧证据不改写。

沿用标准 CAR、共同原始 Gaussian 耦合的 overlapping complex p-body 模型、beta=mu=0 归一 annealed full-Fock 迹、hbar=1/string 实时间。固定 Js 和时间窗 S，令 p→∞，lambda_N=p²/N 位于正紧区间 [lambda_−,lambda_+]。耦合方差为
\[
 w_N=4^p\frac{J_s^2}{\lambda_N}\left(\frac{p!}{N^p}\right)^2.
\]
不取增长时间或联合 lambda→0。

## 1. 归一化与实际 memory 输入

记 d_m=binom(N,m)、tr_m=Tr/d_m、b_(N,m)=2^(−N)binom(N,m)。对 charge 集合 I，
\[
 \|Y\|_{r,I}^r=\sum_{m\in I}b_{N,m}\,
          \mathbb E\operatorname{tr}_m|Y_m|^r. \tag{1}
\]
截取后不重新归一化；省略 I 表示全迹。r 是 Schatten 阶数，L∞ 是逐样本、逐 charge 的 operator norm 本质上界。Hölder 与插值均在该有限迹直和乘 disorder 概率的代数上使用，不需要各因子独立。

在 m>0 上令 h=H_m、K=1_N tensor H_(m−1)，V psi=m^(−1/2)Σ_i |i〉c_i psi，V†V=1。定义
\[
 Q_\perp=1-VV^\dagger,\quad R=Q_\perp KV,\quad
 D_\perp=Q_\perp KQ_\perp|_{Q_\perp},\quad
 a=V^\dagger KV=(1-p/m)h,
\]
\[
 F(t)=V^\dagger e^{iKt}V,\quad
 M(s)=R^\dagger e^{iD_\perp s}R,\quad G=R^\dagger R.
\]
精确 pinching P_hY=Σ_e P_eYP_e 保留整个简并能量块。设
\[
 M_0=\mathcal P_hM,\quad\delta M=(1-\mathcal P_h)M,
 \quad Q_N(s)=m\delta M(s).
\]
Q_perp 是 dilation 投影，Q_N(s) 是放大的 memory 场，均不表示 charge 数算符。m<p 时 H_m=H_(m−1)=0，memory 严格为零；m=0 也置零。既有实际场满足
\[
 \mathcal F_N(s)=(1-\mathcal P_h)
       [m e^{-ihs/2}M(s)e^{-ihs/2}],\qquad
 Q_N(s)=e^{ihs/2}\mathcal F_N(s)e^{ihs/2},
\]
\[
 L_S:=\sup_N\sup_{|s|\le S}\|Q_N(s)\|_2<\infty. \tag{2}
\]
这里相位与 pinching 精确相容；(2) 是同一实际场的界。

## 2. 正力高矩与全部 charge 的常数

令 F_i=[H,c_i]、A=Σ_i F_i†F_i。前驱证明
\[
 G=\mathcal A_m/m-(p/m)^2h^2,\quad 0\le G\le\mathcal A_m/m,
\]
\[
 \mathbb E\operatorname{tr}_m\mathcal A_m^k
 \le(2k-1)!!(Nv_N)^k,\quad
 v_N\le\frac{4J_s^2}{\lambda_N}\frac pN e^{2\lambda_N}. \tag{3}
\]
既有 Gaussian 分部积分证明对矩形 force 作 Hermitian dilation 并停在二阶矩，再用退火 Schatten Minkowski；不假定 F_i 独立。

定义
\[
 A_*=\frac{4J_s^2}{\lambda_-}e^{2\lambda_+},\qquad
 C_*=\frac{J_s^2}{\lambda_-}e^{2\lambda_+}=A_*/4.
\]
于是 Nv_N≤A_*p。正算符本征值单调性给
\[
 \mathbb E\operatorname{tr}_mG^k
 \le(2k-1)!!(A_*p/m)^k. \tag{4}
\]
这里没有声称矩阵幂是算符单调的。矩形 Hölder 的每个因子使用同一源维数 d_m 归一，包括 R† 的 bath 迹也除 d_m，而非除 bath 维数。逐样本有
\[
 \|R^\dagger e^{iD_\perp s}R\|_{r,m}
 \le\|R\|_{2r,m}^2=(\operatorname{tr}_mG^r)^{1/r}. \tag{5}
\]
另有 sector 方差 C_m=w_N binom(m,p)binom(N−m+p,p)≤C_*：直接用 m(N−m+p)≤(N+p)²/4 及 (1+p/N)^(2p)≤exp(2lambda_N)。Gaussian Hamiltonian 四矩给
\[
 \|h\|_4\le3^{1/4}\sqrt{C_*}. \tag{6}
\]

## 3. 广中央区的 2/8→4 插值与低 charge 尾

只分 good 集 G_N={m≥N/4} 和 bad active 集 B_N={p≤m<N/4}；上 charge 尾无需剔除。充分大 N 时 p<N/4。设
\[
 B_r=4\lambda_+A_*[(2r-1)!!]^{1/r},\qquad r=4,8.
\]
由 (4)--(5) 及 pinching 的 Schatten 收缩性，对所有实 s 有
\[
 \|M(s)\|_{r,\mathcal G_N},\ \|M_0(s)\|_{r,\mathcal G_N}
 \le B_r/p,\qquad
 \|\delta M(s)\|_{r,\mathcal G_N}\le2B_r/p. \tag{7}
\]
由 (2) 则得到更小的二范数
\[
 \sup_{|s|\le S}\|\delta M(s)\|_{2,\mathcal G_N}
 \le4L_S/N\le4\lambda_+L_Sp^{-2}. \tag{8}
\]
非交换插值 1/4=(1/3)/2+(2/3)/8 在 (1) 的未重新归一化范数中仍成立。定义
\[
 E_S=(4\lambda_+L_S)^{1/3}(2B_8)^{2/3},
\]
即得
\[
 \boxed{\sup_{|s|\le S}\|\delta M(s)\|_{4,\mathcal G_N}
                \le E_Sp^{-4/3}.} \tag{9}
\]
本式没有声称放大场 Q_N=m deltaM 的全时间四范数一致有界。

bad 区使用 1−P_h 的二范数正交投影性质，因子为 1 而非 2。二项下尾 Σ_(m<N/4)b_(N,m)≤exp(−N/8) 和 (4) 给
\[
 \begin{aligned}
 \|Q_N(s)\|_{2,\mathcal B_N}^2
 &\le\sum_{\mathcal B_N}b_{N,m}m^2
            \mathbb E\operatorname{tr}_mG^2\\
 &\le3A_*^2p^2e^{-N/8}.
 \end{aligned}
\]
故
\[
 \boxed{\sup_s\|Q_N(s)\|_{2,\mathcal B_N}
             \le\sqrt3 A_*p e^{-N/16}.} \tag{10}
\]

## 4. 两端传播子的四范数误差

以有序 convolution (A*B)(t)=int_0^t A(t−s)B(s)ds 记，F′=iaF−M*F。令 F0 解换 M 为 M0 的同一初值问题，F0(0)=1。其收缩性有显式自伴 dilation：
\[
 R_0\psi=(RP_e\psi)_e,\quad D_0=\bigoplus_eD_\perp,\quad
 K_0=\begin{pmatrix}a&R_0^\dagger\\R_0&D_0\end{pmatrix}.
\]
因为 R0†exp(iD0s)R0=M0(s)，exp(iK0t) 的系统压缩就是 F0(t)。各系统能量块及对应 bath 副本分别不变，故
\[
 \|F(t)\|_\infty,\|F_0(t)\|_\infty\le1,
 \qquad [F_0(t),h]=0. \tag{11}
\]
简并块内部不要求 F0 为标量。令 Ua(t)=exp(iat)、Uh(t)=exp(iht)，则精确地
\[
 F-U_a=-U_a*M*F,\qquad F_0-U_a=-U_a*M_0*F_0.
\]
两重时间单纯形体积≤S²/2，(7)、(11) 给 F−Ua 与 F0−Ua 各自的 good 四范数≤S²B4/(2p)。又 a−h=−(p/m)h，由 (6) 给 Ua−Uh 的 good 四范数≤4lambda_+S 3^(1/4)sqrt(C_*)/p。设
\[
 D_S=S^2B_4/2+4\lambda_+S\,3^{1/4}\sqrt{C_*}.
\]
因而两个差分别满足
\[
 \boxed{\sup_{|t|\le S}\|F-U_h\|_{4,\mathcal G_N}\le D_S/p,
 \quad\sup_{|t|\le S}\|F_0-U_h\|_{4,\mathcal G_N}\le D_S/p.} \tag{12}
\]

## 5. 精确差式与最终全迹误差

保持乘法次序的 Volterra resolvent 恒等式是
\[
 m(F-F_0)=-F_0*Q_N*F.
\]
定义 Ecal_N(t)=m(F−F0)(t)+(Uh*Q_N*Uh)(t)。good 区精确展开为
\[
 \mathcal E_N=-m(F_0-U_h)*\delta M*F
             -mU_h*\delta M*(F-U_h). \tag{13}
\]
Hölder、(9)、(11)--(12) 及 m≤N≤p²/lambda_− 使每个 integrand 的二范数≤(D_SE_S/lambda_−)p^(−1/3)。两个项乘时间单纯形体积 S²/2，得到 good 误差 S²D_SE_S p^(−1/3)/lambda_−。

bad 区直接比较两个完整 integrand F0 Q_N F 与 Uh Q_N Uh，每个由收缩性控制；差≤2||Q_N||_(2,bad)。不在这里先提取 m≤N，避免多项式损失。结合 (10) 得
\[
 \boxed{\sup_{|t|\le S}\|\mathcal E_N(t)\|_2
 \le\frac{S^2D_SE_S}{\lambda_-}p^{-1/3}
      +S^2\sqrt3 A_*p e^{-N/16}.} \tag{14}
\]
先对正时间证明；负时间由 F(−t)=F(t)†、F0(−t)=F0(t)†、Q_N(−s)=Q_N(s)† 得到，也可用两重有向积分，其绝对体积不变。

左右包裹 exp(−iht/2) 并施加 1−P_h 不增加二范数，且 F0 的 centered 贡献由 (11) 精确为零。因此 (14) 同样控制实际 bilocal 与单次 memory 滤波的差。相位、Jacobian 和共同谱滤波见主体 [重建证明](susskind-bilocal-reconstruction-proof.md)。

本结果只用旧实际二阶强界、正力 2/4/8 阶矩和显式 dilation；不需要四标记极限、rank-to-size 四范数桥或完整 memory 四范数可和级数。这里没有用有限矩阵数值检查替代渐近证明。
