# 相邻电荷的完整时间传播与 Gaussian OU 近似

2026-11-07，本地。[八页PDF](../output/pdf/de-sitter-susskind-memory-charged-transfer.pdf)、[LaTeX](../manuscript/susskind-memory-charged-transfer/main.tex)、[来源](../literature/susskind-memory-charged-transfer-sources.md)。本轮不commit/push/ZIP。

本轮将上一阶段的普通传播子rank–chaos界接到实际charge-removal核。标准CAR、共同重叠complex p-body Gaussian耦合、β=μ=0、ℏ=1、string实时间；λ_N=p²/N在正紧区间，Js固定。定义T_m X=Σci†Xci/m，F_m(s)=T_m exp(isH_{m−1})。证明

\[
 \sup_{|s|\le S}\sum_{m=1}^N b_m
 \|F_m(s)-\mathcal O_{\alpha_m}(e^{isH_m})\|_{2,m}^2
 =O_S(N^{-2}),\qquad \alpha_m=(1-p/m)_+.
\]

O_a是原Gaussian耦合空间的OU算符，O_a f(J)=E_{J'}f(aJ+√(1−a²)J')。这是保留原样本J的条件平滑；不是把两个电荷独立化，也不是先前的能量q-Mehler。结果覆盖完整指数全部chaos，未作时间截断。

## 两个方向的精确系数

记n=N−m+1，α=(m−p)/m，C=C_m，A=T†，Xcal的rank本征值为x_r=r(N+1−r)。AT=1−Xcal_-/(mn)、TA=1−Xcal_+/(mn)是已知EGUE(1)谱的直接应用，零模按各自定义域保留。

对共同独立副本(u,v)和确定rank-r输入，正向缺陷精确为

\[
 E\|T(vX_r)-\alpha uTX_r\|_+^2
 =\frac{Cp(m-p)x_r}{m^3(m-1)n}\|X_r\|_-^2.
\]

反向须保留V=v−αA(u)=Xcal_-v/(mn)。令D_u(Y)=A(uY)−A(u)A(Y)，则

\[
 vA(Y)-\alpha A(uY)=VA(Y)-\alpha D_u(Y),\qquad
 E\langle VA(Y),D_u(Y)\rangle_-=0.
\]

后项norm²= Cp x_r/[mn²(n+p−1)]·||Y_r||²；源variance ν0=Cp²(m−p)N(N+1)/[m³n(n+p−1)]。所以原通道的annihilation方向不能照抄forward界。正交是协方差恒等式，不是假设独立。

## 极分解与完整时间界

取T=U S_-=S_+U，S_r=√[(m−r)(n−r)/(mn)]，U为partial isometry且UI=I。低chaos (k+1)p≤N/8时仅经过公共rank，S_t≥1/2，不在边界零模上取逆。

按每个r→t比较：δ=p+r−t，|S_t−S_{r+p}|≤32δ/N；α²S_r²−S_{r+p}²有稿(22)的精确式。配合旧一步Eδ²界与新forward Gram，得到polar creation块平方≤10¹² C B_m(k+1)³/N²，其中B_m=λ_N+λ_N²+λ_N(m−n)²/N。低degree的annihilation用公共子空间的adjoint；高degree只用U收缩与Gaussian number指数尾。

Duhamel于是给U exp(isH_-)−exp(isH_+)的norm²≤R_m(S)，稿(30)列有限参数常数及2^{-floor(N/8p)}尾。由T=S_+U，实际F到S_+exp(isH_+)的误差不增。

全部合法rank都有

\[
 |S_r-\alpha^k|\le\frac{4(kp-r)}N+
 \frac{8kp|m-n|}{N^2}+\frac{528k^2p^2}{N^2}.
\]

因此上一阶段的Dfrak=pNcal−Rcal普通指数界、number二/四阶矩，给S_+exp到α^{Ncal}exp的O(√B_m/N)范数误差。结合polar界得到主结论。全charge不需再固定L：Σb_m(m−n)²=N+1，M≤Js²e^{2p/N}/λ_N，宽中央窗外误差平方≤8e^{-N/8}。m≤p的F=I；m=p的OU参考一般不等于I，但在外尾中处理。

## 接下来算什么

本轮误差是无导数L²界。真实力核Λ=−F''+i{αh,F'}+α²hFh需要导数及两侧能量加权；不能直接对O(N^-1)误差求导。

OU参考本身的力组合有精确表达：H'=αh+√(1−α²)u，

\[
 \Lambda_{OU}(s)=(1-\alpha^2)E_u[u e^{isH'}u].
\]

同一个u在两侧及指数内部，不能换成独立Φ sandwich。下一步建立带导数/能量权重的传递界，保留m放大后可能order-one的缺陷，接精确能量回归和实际X0协方差。OU参考一般不与h(J)逐样本对易；本轮未宣布X0紧性或体介子字典完成。

## 核验

三新driver全部首次exit0，科学JSON拒绝覆盖。六个完整CAR系统覆盖全部Gaussian坐标、全部矩阵输入、非方形operator空间和零模；验证forward/reverse完整Gram、complex cross0和polar maps。另有3个符号恒等式及230个exact rational源方差案例。没有从这些小系统拟合统一速率，速率由稿内解析证明给出。

八页最终PDF双遍编译，无undefined/overfull，95dpi全部逐页检查通过。第一次mathscr未加载导致编译失败，源与日志保留diagnostics，改成mathcal后通过。旧阶段证据未覆盖。
