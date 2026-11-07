# 输出缺额的奇偶矩与实际二次 Gram

2026-11-06，本地。[五页 PDF](../output/pdf/de-sitter-susskind-memory-parity-gram.pdf)、[LaTeX](../manuscript/susskind-memory-parity-gram/main.tex)、[来源](../literature/susskind-memory-parity-gram-sources.md)。不 commit/push/ZIP。

本轮求出上一阶段输出缺额 S_N 的实际同耦合等时权重：用既有 Racah 差分方程求带奇偶符号的0/1/2阶矩，不假定完整Poisson分布。它与H²和输入缺额方向的Gram已闭合；与A的cross仍待求，不能据当前正残量宣布新的动态或体粒子。

## 精确带符号矩

固定 N=2m，标准CAR、maximally mixed charge-m态、归一迹。D=binom(N,m)，d_r=binom(N,r)²−binom(N,r−1)²，x_r=r(N+1−r)。独立单位方差的rank a,b Hermitian输入给旧正乘积律π_ab(r)=E||P_r(U_aV_b)||²，Σπ1。记 S=a+b−r，X=x_(a+b)−x_r，f=χ_a(b)，其中χ为旧纯rank4F3。

\[
E_\pi(-1)^S=f,\quad E_\pi[(-1)^S X]=-2abf,
\]
\[
E_\pi[(-1)^S X^2]=(4a^2b^2-2ax_b)f
+a^2(m+1-a)[\chi_{a-1}(b)-f].
\]

a0时最后项取0。公式包含全部允许output ranks，不由有限差分拟合。第二式/第三式可为负，因为这是signed moment。

证明入口：τ=star transpose star^-1在rank r上为(−1)^r；normalized Hermitian矩阵基完备性给Στ(Tμ)XTμ=(−1)^mDτ(X)，即Σ_r(−1)^r d_r χ_r(j)=(−1)^(m+j)D。Racah差分L=B_j(1−shift+)+D_j(1−shift−)，其中

\[
B_j=\frac{(j-m)(j-N-1)^2}{2(2j-N-1)},\quad
D_j=\frac{j^2(j-m-1)}{2(2j-N-1)}.
\]

Lχ_a=x_aχ_a，端点B_m=D_0=0；K=N(N+1)/2，2(B+D)=K−x，Lx=N(x−K/2)，εLε=K−x−L。把signed transform代入已有π三重和，令Xop乘x_a，得到signed x_r一阶=(Xop+L_a)χ，二阶=[(Xop+L_a)²+N(K/2−Xop−L_a)]χ。消去χ(a+1)后得到上面的稳定第二矩，O(N³)项已合并。

## DS 与奇偶选择

令a/p→α,b/p→β,p²/N→λ>0，θ=2λαβ。Racah第k项可一致由Kλ^k/k!控制，固定k趋(−4λαβ)^k/k!，所以χ→e^(−2θ)，相邻χ差→0。第二矩修正的a²(m+1−a)/N²有界。故

\[
(E(-1)^S,E(-1)^SX/N,E(-1)^SX^2/N^2)
\to e^{-2\theta}(1,-\theta,\theta^2-\theta).
\]

X=S[N+1−2(a+b)+S]，a+b≤2p，旧正矩界让前二阶signed极限可用S代替X/N。结合旧unsigned矩：

\[
E[1+(-1)^S]\to1+e^{-2\theta},\quad
E[S(1+(-1)^S)]\to\theta(1-e^{-2\theta}),
\]
\[
E[S^2(1+(-1)^S)]\to\theta^2(1+e^{-2\theta})+\theta(1-e^{-2\theta}).
\]

τ直接给P_r(V_bU_a)=(−1)^S P_r(U_aV_b)，从而sym强度=4·1_(S even)π。单缺额S1严格不进入对称乘积。仍未证明全Poisson律。

## 真正同耦合 Gram

完整p-body complex Gaussian H=Σh_a，CaI=E h_a²，C=ΣCa，ωa=Ca/C，bar s_N=Σ(p−a)ωa。定义

\[
B_N=W_2H^2=H^2-CI,\quad
L_N=W_2\{H,D_N-\bar s_NH\},\quad D_N=\sum_a(p-a)h_a,
\]
\[
\mathcal S_N=\sum_{a\le b,r}(a+b-r)W_2Z^+_{ab;r}.
\]

L_N在本笔记只表示输入缺额二次方向，不是旧memory tangent或Liouvillian。不同输入pair及输出rank正交，a=b的W2通过独立副本decoupling处理，给任意对称filter的精确式：

\[
\left\langle\sum fW_2Z^+,\sum gW_2Z^+\right\rangle
=C^2\sum_{a,b}\omega_a\omega_b\sum_r[1+(-1)^S]\pi_{ab}(r)fg.
\]

输入缺额的两个权重变量独立趋Poisλ，bar s_N→λ；固定输入deficits时θ=2λ，q=e^(−4λ)。factorial input尾和第三output矩控制所有二次/混合矩。以(B,L,S)排序，Gram/C²极限为

\[
\begin{pmatrix}
1+q&0&\theta(1-q)\\
0&2\lambda(1+q)&0\\
\theta(1-q)&0&\theta^2(1+q)+\theta(1-q)
\end{pmatrix}.
\]

减B投影μ=θ(1−q)/(1+q)，S_perp=S_N−μB_N。物理C→Js²/λ：

\[
\|\mathcal S_\perp\|_2^2\to J_s^4\left[\frac{2(1-q)}\lambda+\frac{16q}{1+q}\right]>0,
\quad \langle S_\perp,L_N\rangle\to0.
\]

再取小λ范数趋16Js4。旧F bridge给J_HF−4λqC L_N→0；线性误差的Gaussian四阶矩界和Hölder将其推到乘积，因此S_perp–JhF cross也趋0。这里只排除B/L两个方向，不包含A或完整commutant，不是时间谱/动力学结论。

## 验证和继续

- signed-moments.json：N2..12偶数、全部139input rank pairs，417 signed矩及139差分式；符号Gram投影和smallλ16核对。
- car-gram.json：N4p1m2/N4p2m2/N6p1m3，分别136/666/666全部Gaussian对称坐标对，完整3×3Gram由整数CAR/Fraction计算，独立重现π奇偶公式。
- limit-checks.json：24组λ=.1/.25/.5,p16..128,deficits(0,0)/(1,2)，exact Fraction公式后80/120位求值。有限变量是X/N，不是S；偏差保留，最后一组signed差约(.004323,−.001067,−.002125)，非区间认证。

初版(−1)**负指数在π=0通道把Fraction和式转成float导致断言失败；改为整数奇偶分支，保留初版与说明，无公式/容差调整。五页双遍编译首次通过，95dpi最终五页逐页核版完成。

下一步专攻S_N–A的下charge相干cross和A加权尾，不能用已有S正范数替代完整旧span投影。随后恢复W4后代low chaos/raw Z2谱、中央charge窗、finite-memory和体空间/质量字典。整体goal仍active。
