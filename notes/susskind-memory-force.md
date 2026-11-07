# 从带电force到平均记忆核

2026-11-06，本地阶段。[LaTeX](../manuscript/susskind-memory-force/main.tex)、[PDF](../output/pdf/de-sitter-susskind-memory-force.pdf)、[来源](../literature/susskind-memory-force-sources.md)。沿用Alan Day与Chengheng的研究记录署名。本轮不commit/push/ZIP。

## 结果与规范

标准CAR、含指标重叠的complex Gaussian p-body H=Hs，w=4^p(Js²N/p²)(p!/N^p)²，全归一β=μ=0迹、ℏ=1。先N,p→∞且p²/N→λ>0，再考虑λ→0；q=exp(−4λ)、σ²=Js²/λ，string时间固定在紧区间。

固定charge m，Vψ=Σ|i〉c_iψ/√m，h=H_m，K=1⊗H_(m−1)，P=VV†、Q=1−P、R_ch=QKV、K⊥=QKQ。上一轮M(s)=R_ch†exp(iK⊥s)R_ch；定义能量居中Mhat(s)=exp(−ihs/2)M(s)exp(−ihs/2)，m0取0。

本轮得到标量annealed mean：

\[
p\,E\operatorname{tr}_0\widehat M(s)\longrightarrow4J_s^2G_q(s),
\]

固定紧时间一致。G是文献已有单位归一q-Mehler探针核，probe参数r=q；后续小λ为sech²(2Js s)。这里求出的是平均记忆，未把它当成算符值极限或涨落协方差。既有Mehler公式不属于原创结果。

## 精确renewal与导数接口

B_n=V†K^nV、B0=1，b1=B1，b_n=B_n−Σ_(j=1)^(n−1)b_j B_(n−j)。按首次回到P拆分，b_(j+2)=R_ch†K⊥^jR_ch。乘法次序不可改。该构造来自标准first-return/非交换renewal方法。

r=p/m、J_h X={h,X}/2、A=Σ[H,c_i]†[H,c_i]、M3=Σc_i†ad_H³c_i、T3=(M3+M3†)/2，则

\[
b_3=T_3/m+(1+2r)J_h A/m-r^2(1+r)h^3.
\]

G=M(0)=A/m−r²h²≥0，居中线性系数Z1=b3−J_h G满足

\[
mZ_1=T_3+2rJ_hA-mr^3h^3.
\]

由B2=A/m+(1−2r)h²、B3=T3/m+3J_hA/m+(1−3r)h³直接代入得到。三组精确有理复矩阵核renewal至六阶及两条一阶式。R_ch与旧三阶singlet R是不同对象。

## 全阶带电收缩

D_i=[H,c_i]，F_N(s)=(2p/N)ΣEtr0 D_i†D_i(s)。其双能量权重为非负的矩阵元模方，支持任意固定双能量矩分析。

两个D端点彼此Wick配对时，求导site i有p个选择。剔除与本弦另一集合或其它弦集合相撞的标签，平均只损失O(1)。孤立标签在原H字中给n_i或1−n_i，trace=1/2；两导数去掉该单site占据因子，剩余字保持原次序，trace精确变为2倍。带标签弦去掉一个site只改变集合长度p→p−1，固定λ下其余交叠规则不变，故是一条被强制配对的H型matter chord，与每条H弦交叉权q。

如果两个D端点分别配外部不同H，则同一site必须落入两条独立弦的相关集合交集，期望选择数O(p²/N)=O(1)，比同弦p少一阶。固定字长下，无放回抽样的阶乘矩控制交叠数的指数矩，额外2^overlap局部trace权不会改变该计数阶数。三重交叠概率O(p³/N²)消失；同弦有效计数p+O(1)，异常及跨弦计数加权后仍O(1)。这保留同J收缩，而非假设D是独立probe。

未归一系数2pσ²，乘2p/N得到4λσ²=4Js²。所有固定双能量矩趋4Js² μ_q。正测度与更高偶矩给紧性/一致可积性，候选q-Mehler测度紧支撑、矩确定，故弱收敛。二阶频率矩统一界给固定紧时间一致收敛。有限CAR字检查仅核局部代数，不替代这一全阶论证。

## 正交记忆到force的误差

固定m写D=KV−Vh=R_ch−rVh、g=tr_m G、c=tr_m h²。精确确定性界：

\[
|\operatorname{tr}_m e^{-ihs}R_{ch}^\dagger(e^{iK_\perp s}-e^{iKs})R_{ch}|
\le s^2\operatorname{tr}_m G^2/2,
\]

\[
|\operatorname{tr}_m e^{-ihs}(R_{ch}^\dagger e^{iKs}R_{ch}-D^\dagger e^{iKs}D)|
\le2r\sqrt{gc}+r^2c.
\]

第一条用block-diagonal K0与K间两次Duhamel，Q→P→Q至少两次跃迁；外面再有两个R_ch，Schatten-4 Hölder给trG²，积分面积s²/2。第二条直接展开并用HS Cauchy–Schwarz。

中央窗|m−N/2|≤N^(2/3)，r=O(1/p)、Eg=O(1/p)。旧singlet结果Etr0(A−a(Qf))²=O(1)，且a_m=C_m p(N+1−p)/(N−m+1)=O(p)，给中央窗平均EA²=O(p²)。由0≤G≤A/m得EtrG²=O(1/p²)。因此乘p后两误差为O(s²/p)、O(p^−1/2)。中央窗1/m对2/N相对差O(N^−1/3)。

尾控制：T(X)=Σc†Xc；A≤2T(H_-²)+2m h²、T(X)²≤mT(X²)，Gaussian Etr_m h⁴≤3C_m²，得Etr_mG²≤24(C_(m−1)²+C_m²)。C_m≤(Js²/λ_N)exp(2λ_N)及binomial尾≤2exp(−2N^(1/3))压过多项式。m<p的H块为0。由此接完整全迹，而不凭半填充点值跳过charge尾。

## Mehler独立检验及原始传播

已有μ_r=ρ_q(θ)ρ_q(φ)K_r(θ,φ)dθdφ、Eθ=2σcosθ/√(1−q)。直接拆q乘积n0项得

\[
\frac{(E-E')^2}{1-r}\,d\mu_r\xrightarrow{r\uparrow1}2\sigma^2d\mu_q.
\]

用|1−e^{i(θ+φ)}|²|1−e^{i(θ−φ)}|²=4(cosθ−cosφ)²消去零点；剩余q乘积直接变K_q。对角线上左密度恒0，故只作a.e.密度与weak measure结论。|1−r e^{iα}|²≥r|1−e^{iα}|²及剩余乘积的正下界给支配，无额外δ。此式交叉验证归一，不代替上面的同J微观证明。

f_N(t)=(2/N)ΣEtr0 c_i†c_i(t)，f_N(0)=1，Gaussian符号对称使平均f_N′(0)=0。trace Hilbert空间的ad_H自伴性给−p f_N″=F_N，故

\[
p[1-f_N(t)]\longrightarrow4J_s^2\int_0^{|t|}(|t|-s)G_q(s)\,ds.
\]

后取小λ得log cosh(2Js t)，恢复已知large-interaction-order传播展开。MS规范interaction order q_MS=2p，区别于本项目chord q。此为固定紧时间首个1/p项，不推出t~p的指数重求和。

## 保存的验证与下一步

bridge.json：3精确系统，18个renewal阶数、6个一阶/居中式；3系统共12个时间点核两比较界。首轮exact helper用int截掉半整数导致失败，保留初版与错误，改用精确dyadic Rational后通过。

kernel-checks.json：130 tagged CAR字、17非零；25精确弦多项式；9组65位离对角Mehler极限。首轮120随机字非零覆盖不足，保留失败，追加10个定向非零字，不降阈值。四λ曲线由Jacobi及+12 refinement；解析vacuum Taylor截断<1e−13；λ=.25/.1另用128/192点Mehler积分。舍入与积分误差非区间认证，有限λ偏差完整保留。

下一步明确减守恒部分、放大归一化后的Mhat(s)算符协方差，保留两个memory参数及独立外部中心时间。零阶mG_residual=Y，一阶有上述精确Z1接口；完整有限memory时间的涨落尚未得到。此对象才可继续与体介子内部核比较。普通平均核不直接确定空间平移或质量。完整goal active。
