# 从全次数 rank 界到完整普通传播子

2026-11-07，本地。[八页PDF](../output/pdf/de-sitter-susskind-memory-wick-rank.pdf)、[LaTeX](../manuscript/susskind-memory-wick-rank/main.tex)、[来源](../literature/susskind-memory-wick-rank-sources.md)。本轮不commit/push/ZIP。

## 本轮实际推进

上一轮已把第一次返回归入守恒背景；剩余实际场是 X0=(1−P_h)m exp(−ihs/2)Λ(s)exp(−ihs/2)。本轮尚未证明 X0 有界，但已取得**完整普通传播子 exp(isH) 的 rank–Gaussian-chaos 缺额统一界**，包含所有低 chaos 收缩。这比仅对每个固定 Taylor 次数、或只对最高 chaos 估计更进一步。

规范沿用标准CAR、含指标重叠的complex Gaussian p-body H、共同J、β=μ=0 annealed normalized full trace、ℏ=1/string实时间。偶N，λ_N=p²/N→λ>0在正紧集；先固定实时间窗。下面有限参数的中央界要求N≥8、p≤N/8、N/4≤m≤3N/4，DS序列最终满足。

## 正 rank 链与全次数缺额

固定charge，C=Etr_m H²，Φ(X)=E_(独立H') H'XH'是确定性协方差map。中性算符rank r的特征值κ_r=w binom(m−r,p)binom(N−m+p−r,p)，μ_r=κ_r/C。令H1…Hk为独立副本乘积，π_k(r)=C^−k E||P_r(H1…Hk)||²。既有Racah角色χ_r(s)给

\[
\sum_r\pi_k(r)\chi_r(s)=\mu_s^k.
\]

正性来自真实乘积的U(N)不变sesquilinear协方差及Choi分解；乘积不必Hermitian。一步转移P(r,t)=Σ_a ω_a π_ra(t)是真概率，rank三角性给r≤kp和|t−r|≤p。这里rank为表示论标签，不是矩阵秩。

设x_r=r(N+1−r)、e=m(N−m)(N+1)/N、u=μ1、v=μ2，X_j=x_(r_j)。前两个Racah角色给精确条件均值ux+c和条件方差(v−u²)x²+Bx+V1。全部有限N有理系数见稿第2–3页。v≤u²，使方差中的大二次项可以直接舍去，避免高次数大项相减。

中央区间有|B|≤102400p、V1≤512p²N、c≤8pN、1−u≤16p/N，故VarXk≤410112 k²p²N。kp≤N/8时对x_(kp)−Xk转换，大kp时直接D²≤k²p²，得到对**所有k**成立的

\[
E_{\pi_k}(kp-r)^2
\le2\cdot10^6(\lambda_N k^2+\lambda_N^2k^4).
\]

常数粗但完全显式，不是拟合。固定k半填充时ED→λk²、VarD→λk²、ED²→λk²+λ²k⁴；只证明这些矩，没有声称完整Poisson分布。

## 从Wick级数到普通指数

对最高coupling chaos Zk=W_k(H^k)，输入Gaussian张量的Sym为正交投影，给

\[
\|w(r)Z_k\|_2^2\le k! C^k E_{\pi_k}|w(r)|^2.
\]

只有一个k!，因此缺额加权Wick指数的平方范数由e^(C|s|²)乘四次多项式控制。这是最初获得的结果；随后没有停在Wick指数，继续处理实际普通指数。

在Gaussian L²⊗HS的U(N)协变子空间中，令Ncal数coupling chaos，Rcal=ΣrP_r，Dfrak=pNcal−Rcal≥0。写H=ΣgαBα、ΣBα²=CI，左乘H等于creation Ccal加其adjoint。degree-k协变系数输出rank≤kp。

一步转移缺额δ=p+r−t∈[0,2p]。同一条件矩给全k界

\[
\sup_{r\le kp}E_r\delta^2
\le10^6[\lambda_N(k+1)+\lambda_N^2(k+1)^2].
\]

关键的cross消去是确定性事实：Σα L_(Bα)P_t L_(Bα)为U(N)等变正算符，在不同输入rank间严格正交，同rank标量为CP(r,t)。先作用未Sym creation张量，再用Sym收缩，得到[Dfrak,Ccal]块范数界；annihilation用负adjoint。于是

\[
\|[\mathfrak D,H_{\rm left}]\psi\|^2
\le4\cdot10^6 C(\lambda_N+\lambda_N^2)
\|(\mathcal N+1)^{3/2}\psi\|^2.
\]

对F(s)=exp(isH)，任意ℓ阶Duhamel有序积分项W=U0 Bα1 U1…Bαℓ Uℓ，逐样本ΣαtrW†W=C^ℓ。ℓ!个排列与单纯形体积抵消，给Σ||∂^ℓF||²≤(Cs²)^ℓ。Gaussian Hermite等距于是给

\[
\sum_k\zeta^k\|W_kF(s)\|_2^2\le e^{Cs^2(\zeta-1)},\quad\zeta\ge1,
\]

以及||(Ncal+1)^(3/2)F||²≤1+7z+6z²+z³，z=Cs²。这控制完整chaos分布。

由于Dfrak I=0，有限N的Gaussian Sobolev共同域上可作commutator Duhamel，得到实时间的主要结果：

\[
\boxed{\|\mathfrak D e^{isH}\|_2^2
\le4\cdot10^6 C(\lambda_N+\lambda_N^2)s^2
(1+7z+6z^2+z^3).}
\]

全charge外尾由Dfrak≤pNcal给2p²exp(−N/8)(z²+z)，z=M s²、M=sup C_m。因此完整β=μ=0迹也在固定实时间窗统一有界。这里没有把复时间的unitarity界用于普通指数；最高chaos Wick级数的entire结论是另一个范围。

## 实际covariance map与仍缺的电荷接口

正部延拓μ(x)=∏_(j<p)(1−x/(m−j))_+(1−x/(N−m+p−j))_+，中央区间全局Lipschitz≤12p/N。保留有限N的μ(kp)，得到

\[
\|p[\Phi/C-\mu(p\mathcal N)]e^{isH}\|_2
\le12\lambda_N\|\mathfrak D e^{isH}\|_2=O_S(1).
\]

右边μ(pNcal)是coupling-chaos乘子，**没有证明它逐样本与H对易**。因此该结论不自动给pinched X0的统一范数；剩余应接实际相邻charge的coherent transfer，不能停留于普通neutral covariance map。

另一方面，取独立Gaussian线性charged probes Mξ=Σξ_i c_i，F_m(s)=V†exp(iKs)V=m^−1 Eξ Mξ†exp(iH_-s)Mξ。两个F插入精确等于两个独立颜色的四点。对B_a(t1,t2)=exp(iat1)F(t2−t1)exp(−iat2)，混合端点导数为exp(iat1)Λ(t2−t1)exp(−iat2)，a=(1−p/m)h。它保留共同H，不假定Fi独立。

用正二能量测度ρ(dx,dy)=Etr[P_h(dx)V†P_K(dy)V]，定义精确确定能量回归ℓ_s(x)=Eρ[(y−(1−p/m)x)²exp(isy)|x]，则

\[
\|X_0(s)\|_2^2\le\sum_m b_m m^2\mathcal V_{N,m}(s),\qquad
\mathcal V=E\tr\Lambda^\dagger\Lambda-\int|\ell_s|^2d\mu\ge0.
\]

这是充分上界而非等式；差为样本依赖守恒部分的剩余方差。中央一致V=O(p^−4)加尾界足够；旧strong仅给Σb_mV=o(p^−2)，不误称sector一致。下一步需要将新all-chaos控制接到这个实际余量，或直接控制去能量块的charged场。

## 文献边界与证据

BNR有限迹(3.11)–(3.12)保留多重交叠；其后领先双标度简化不含本题所需全部site-one修正。一个标记位点同时命中两独立H集合的概率(p/N)²=λ_N/N，两probe同site概率1/N，都与connected阶有关。已有终式neutral probe不能直接代替charged单c_i。任何候选核必须恢复旧等时W=4Js⁴[1/λ+4q/(1+q)]，不能丢后面的共同J相干项。

三科学JSON首次运行通过且禁止覆盖：252精确Racah/Markov等式、1800有理界、4一般符号式和6个固定k符号极限；3独立CAR covariance channel、9完整最高chaos tensors、3charged probe finite-time系统；18完整Hermite系数案例（含H⁴低chaos）、43一步分布。80/120位五组exact-formula核对不作舍入区间认证，k100例仅核计算而非固定k极限。

最终八页双遍编译、95dpi逐页检查；首轮六页稿的literal qquad排版错误源/PDF保留diagnostics。旧PDF skill路径失效，定位新版26.1105.12137后marker成功一次。所有输入/输出和历史文件哈希绑定local-evidence。整体目标active，研究尚未完成。
