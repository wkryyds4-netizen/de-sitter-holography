# 相邻电荷的相干乘积与 A 交叉项

2026-11-06，本地。[六页PDF](../output/pdf/de-sitter-susskind-memory-coherent-cross.pdf)、[LaTeX](../manuscript/susskind-memory-coherent-cross/main.tex)、[来源](../literature/susskind-memory-coherent-cross-sources.md)。不commit/push/ZIP。

本轮求出了上一阶段缺失的交叉项。半填充、先固定正λ的DS极限，输出缺额S_N扣掉H²投影后，与同样扣H²的A的内积，等于S自身的残量范数平方。下一步应计算A−S的范数，不能再把两者当作正交方向。

## 约定与相干变换

标准CAR、完整重叠complex Gaussian p-body系综；固定charge m、归一迹、annealed HS；h=H_m、k=H_(m−1)使用同一耦合，C=Eh²，O=W2A，A=Rem(k²)+(2p−m)h²，B=W2h²。W2是普通Gaussian二次chaos。

一般m的共同rank范围R=min(m−1,N−m)，D+=binom(N,m)、D−=binom(N,m−1)，d_r=binom(N,r)²−binom(N,r−1)²。charged(t+1,t)维数e_t=binom(N,t)[binom(N,t+1)−binom(N,t−1)]，ell_r²=m(m−r)(N−m+1−r)/(N−m+1)。

将同一primitive rank-a tensor在两charge分别归一为EU±²=I，跨charge covariance Psi_a(X)=EU+ X U−的特征值：

\[
\zeta_a(t)=\ell_a T_a(t)/m,\quad
T_a(t)={}_4F_3[-a,a-N-1,-t,t-N;1-m,m-N,-N-1;1].
\]

它由旧embedded charged binomial特征值和上charge rank权重的三角分解固定。精确有限求和恒等式见稿(5)，不由小N外推。

定义独立输入a,b的有序相干权重M_ab(r)=Etr_m[P_r(U+V+)†Rem(P_r(U−V−))]，mu_ab=m M_ab/(ell_a ell_b)。composition与covariance reshuffling给

\[
T_a(t)T_b(t)=\sum_r\mu_{ab}(r)T_r(t),\quad
\mu_{ab}(r)=\frac{d_r\ell_r^2}{m^2D_+D_-}\sum_t e_tT_a(t)T_b(t)T_r(t).
\]

保留所有multiplicity及相位。sum mu=1，推导不依赖mu正性；本轮有限表中最小值均0，不能称已观察到负权重。

## 精确奇偶矩与大项抵消

半填充N=2m；x_r=r(N+1−r)，S=a+b−r，X=x_(a+b)−x_r。chi_ab为旧neutral Racah核，xi_ab只将它的一个下参数−m换成1−m：

\[
\xi_{ab}={}_4F_3[-a,a-N-1,-b,b-N-1;1-m,-m,-N-1;1].
\]

Cf=m²(N+1)/N，Ct=m(m−1)(N+1)/(N−1)。精确

\[
\sum\mu=1,\quad\sum X\mu=-2ab+x_ax_b/C_t,
\]
\[
\sum(-1)^S\mu=\xi,\quad
\sum(-1)^SX\mu=-2ab\xi+m(\xi-\chi).
\]

证明使用KS已有Racah差分式，参数(α,β,γ,δ)=(-m,-m−1,-N−2,0)，端点先消因子。B_a=(a+1−m)(a−N−1)²/[2(2a−N−1)]，D_a=a²(a−m−2)/[2(2a−N−1)]，特征值t(N−t)。signed completeness为Σ_r(-1)^r d_r ell_r² T_r(t)=(-1)^(m−1+t)mD−(m−t)。

关键contiguous简化：xi=T_a(b)−b[T_a(b−1)−T_a(b)]/(N+1−2b)。若c_j为chi第j项，则xi第j项=mc_j/(m−j)，所以x_b[T_a(b−1)−T_a(b)]/(N+1−2b)=−m(xi−chi)。b0时整个差项定义为0。未证明完整输出分布的Poisson极限。

真实p-body输入的上下方差和ell因子相乘恰为g=(m−p)²/m，与a,b无关。故same-J精确cross为

\[
\langle\sum fW_2Z^+,O\rangle
=C^2\sum_{a,b}\omega_a\omega_b\sum_r f(1+(-1)^S)
[g\mu_{ab}(r)+(2p-m)\pi_{ab}(r)].
\]

a=b decoupling系数已包含；mu在r=m补0。eta=p²/m，Delta=xi−chi，uX=−2ab+x_ax_b/Cf。每个有序输入对的B、X/N cross（除CaCb）简化为

\[
b_{ab}=\eta(1+\chi)+g\Delta,
\]
\[
z_{ab}=N^{-1}\{\eta(u_X-2ab\chi)+g[x_ax_b/(m(m-1)(N+1))+(m-2ab)\Delta]\}.
\]

所有O(m)抵消先在精确式完成。

## 双标度的实际交叉项

固定λ=p²/N>0，q=e^-4λ，theta=2λ，先DS再考虑smallλ。终止级数有Kλ^j/j!统一支配，m(xi−chi)→−4λq；且b、z对全部a,b≤p一致有界。input-deficit factorial tail给完整rank和的极限，中心输入缺额L_N=W2{h,D_N−bar s h}与O的cross→0。

X/N仍不是S。用X=S[N+1−2(a+b)+S]和上轮positive Gram得||S_N−X_N/N||_m≤(4p/N)||S_N||_m=O(N^-1/2)。旧full-Fock W2A范数有界，而半填充charge的binomial权重~N^-1/2；正性仅给||O||²_m=O(N^1/2)。这个弱界已足够使proxy cross误差O(N^-1/4)→0。**没有把full-Fock范数当成固定半填充范数。**

所以

\[
\langle B,O\rangle/C^2\to\theta(1-q),\qquad
\langle L_N,O\rangle/C^2\to0,
\]
\[
\langle S_N,O\rangle/C^2\to\theta^2(1+q)+\theta(1-q)
=\lim\|S_N\|^2/C^2.
\]

令c=theta(1−q)/(1+q)、Sperp=S_N−cB、Y=O−cB，则cross(Sperp,Y)→Vs=lim||Sperp||²=Js4[2(1−q)/λ+16q/(1+q)]。因此Sperp与O−S_N正交于极限。尚需Delta_N=||O−S_N||²_m；若它有有限极限Delta*，扣B/O后的S残量为Vs Delta*/(Vs+Delta*)。这只是条件公式，未求Delta*，也未完成完整旧span或动态谱。

## 验证与保存

4JSON全部第一次执行通过：transform-checks（1925三角协方差、139正交、139差分、556矩），charged-map-checks（10完整矩形map、306152 complex entries），car-cross（26227全部Gaussian对称坐标对；4×4Gram），limit-checks（384contiguous、2符号、3独立exact全和、12组80/120位完整rank求和）。有限(4,1,2)/(6,1,3)/(6,2,3)扣B/A后S残量均正，仅作有限证据。N6p2计算68.3秒后正常终态。

高精度为X/N cross；λ.5,p128的(B、proxy、L)绝对偏差约(1.95318e−5,−.0548285,3.32957e−5)。无rank cutoff，保存相位抵消因子，非区间证书。

LaTeX两次失败为mathscr宏缺包、随后mathrsfs包本机未安装；改标准mathcal，初始稿/两份失败日志保留，公式未改。最终双遍6页，无undefined/Overfull，95dpi六页逐一检查通过。输入/输出、旧阶段文件与HEAD由local-evidence核对。

下一步先求固定半填充Delta_N及A的正weighted tail；不重复已完成cross，也不将旧full-Fock Wd直接代入Delta*。随后处理g_A取常数、中央charge窗、W4低chaos/rawZ2，最后finite-memory及体空间/质量字典。goal active。
