# 原始三阶 singlet 的一次 Gaussian 分量

2026-11-06，本地阶段；[六页PDF](../output/pdf/de-sitter-susskind-cubic-linear.pdf)、[LaTeX](../manuscript/susskind-cubic-linear/main.tex)。本轮不commit/push/ZIP。

本轮补全了此前W₃投影之外的算符映射。W₁T₃并不总是守恒；有限N反例之外，扣除逐电荷sector的H投影后，其残量F在DS中有正范数64Js⁶q²。它相对R三阶通道的等时强度小λ为2λ，提供明确的下一谱计算对象。

## 精确有限N公式

标准CAR，H=Hs，L=adH，A=Σ(Lc)†Lc，M₃=T₃+(3/2)LA。Gaussian covariance E J_IJ J_KL=w δ_ILδ_JK，w是Hamiltonian的方差系数，不是λ=p²/N。以下每个binomial越界取0。

H_m=Σν H_m^(ν)，ν≤min(p,m,N−m)，为已有U(N)不可约分解。

C_m=w binom(m,p)binom(N−m+p,p)；
κν(m)=w binom(m−ν,p)binom(N−m+p−ν,p)；
jc=w binom(m−1,p)binom(N−m+p,p)；
j−=w binom(m−ν,p)binom(N−m+p−ν+1,p)；
j+=w binom(m−ν−1,p)binom(N−m+p−ν,p)；
ellν=ν(N−ν−p+1)/(N−2ν+1)。

W₁T₃|m=Σν tν(m) H_m^(ν)，其中

tν=−(m+2p)C_(m−1)+(m−3p)C_m+3p jc
 +(m−p)κν(m−1)−mκν(m)
 +3[ellν j−+(p−ellν)j+]。

这些是有限N表达式，不是DS拟合。来源中的κ和rank分解已存在；本轮计算是其在raw T₃一次chaos的应用。

推导从L³c的四个monomial中每次收缩两个H。跨电荷协方差J(Y)=E K_(m−1)Y K_m，对traceless charged(r,s) tensor的系数为w binom(m−s,p)binom(N−m+p−r,p)，由disjoint最高权正规序计数得到。令Hν_m=binom(m−ν,p−ν)Oν，向量cOν和Oνc只含charged下/上两种rank；Hilbert–Schmidt投影给ellν及p−ellν。全部系数实，所以W₁M₃=W₁T₃且W₁LA=0。

N5p2m3、w1时(t0,t1,t2)=(-96,-167/2,-88)。saved validation给一个显式整数Hamiltonian，tr₃[H,W₁T₃]†[H,W₁T₃]=8933634/5，精确非零。前两个小系统恰好对易，不能只用它们推出守恒。

## 逐sector能量投影和DS残量

已有rank权重

ων=dν Λν(N,m,m−p)/[binom(N,m)Λ0(N,m,p)]，
dν=binom(N,ν)²−binom(N,ν−1)²，
Λν(N,m,r)=binom(m−ν,r)binom(N−m+r−ν,r)。

定义tbar(m)=Σωνtν，F=W₁T₃−tbar(Qf)H。则〈HF〉m=0且

〈F²〉m=C_m³ Var_ω(tν/C_m)。

tbar(Qf)H每个realization守恒；F只是正交于这一个方向，不是已经投影掉完整commutant。

半填充先DS λ=p²/N固定，s=p−ν的精确相邻权重比趋λ/(s+1)，且对所有s可被32λ/(s+1)控制，因此s→Poisson(λ)，全部多项式矩收敛。t_(p−s)/C−t_p/C→8λq s；相邻系数差由rational factors统一有界，故可用factorial tail取variance。得到

Var_ω(tν/C_m)→64λ³q²，q=e^(−4λ)。

完整β=μchem=0 trace的电荷binomial在m/N=1/2集中，相关估计在附近密度窗口一致；C_m最大在(N+p)/2附近、Cmax/C_(N/2)=exp[O(p³/N²)]→1。远端电荷指数尾压过系数的多项式增长。因C→Js²/λ，

〈F²〉_DS=64Js⁶q²。

相对前轮R的范数，比例64λ³q²/V~2λ。一次与三次chaos等时正交，但同J实时演化不保留chaos，不能把这些谱未经计算独立相加。也不能只由ν/p→1断言F为d1/q探针。

完整精确算符账本：

M₃=R+aS+bZ+tbar(Qf)H+F+(3/2)LA。

它不意味着raw M₃的完整两点谱已经求完。

## 验证

- cubic_linear_probe.py：探索性有限Gaussian极化，N3p1/N4p2提示守恒，N5p2m3明确非零；保留全部浮点拟合结果，后续以精确通道公式取代拟合判断。
- cubic_linear_channels.py：N3p1,N4p2,N5p2,N6p2,N6p3的42个sector/rank，直接Gaussian整数CAR极化和公式完全相等，最大保存整数entry57600。
- check_cubic_linear.py：31个charged协方差公式、199个精确rank概率归一及相邻比、3个symbolic投影恒等式、精确非守恒witness。
- cubic_linear_limit.py：λ=.1,.25,.5，p16/32/64/128，所有p+1个rank无截断，80/110位复算；输出60位variance一致，记录max_relative_precision_change为0.0指保存精度一致，不是严格零误差。实际finite-p偏差完整保留：p128三组约−1.723%、−.805%、+.744%。

下一步优先算F的动态矩与strength density，保留H–F同J cross contractions；结合精确rank协方差而不是继续重复有限N拟合。待此谱明确再做X(d2)/R(d3)/F混合高点与相互作用。
