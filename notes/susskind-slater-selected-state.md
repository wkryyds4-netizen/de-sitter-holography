# 单次Slater选态激发：完整协变密度与实际曲率动力学

2026-11-08。沿用同原J、半填充固定Slater O、annealed phi_N、standard overlapping complex CAR、DS p²/N→固定lambda>0。h=H_m/√C_m，W=√lambda u，j=C_m^-1(2p-D-kappa)(H²-C_Q)|m，x=j/√v；D/L先fullFock作用后restrict，C_m固定目标sector。q=e^-4lambda，gamma=e^-lambda，s=gamma²=√q，kappa=2lambda(1-q)/(1+q)，v=4lambda+16lambda²q/(1+q)。h时间无量纲、hbar1。M†为同J selected occupied→empty块，phi(MM†)=1。

## 1. 密度扩展至包括曲率的整个协变代数

令D=W*(h,u,x)，tau为其已证faithful trace，sigma1(B)=phi(M B M†)。本轮证明对全部B∈D：

sigma1(B)=tau(rho1 B)，rho1=s(h-W)²+r_P I，r_P=1-s(1+lambda)>0。

tau rho1=1，r_PI≤rho1≤[r_P+4s(1+lambda)/(1-q)]I；上界是实际norm。旧结果仅在h,u子代数成立；曲率的加入现在有独立全词证明。两端可添任意固定协变dressing A：A M†制备的密度是A rho1 A†，归一需除tau(A rho1 A†)。不把这个公式用于任意多次interleaved M，也不替代off-diagonal phi(j(M†)²)。

## 2. 一般证明的边界图机制

词M B M†两端放在0,L，B只含H/W/J。若两边界互配，outerF包围全部内部G；没有innerF nesting，所有marked coefficients不变，贡献tau(B)。否则恰有两条F，分别触到0或L；它们都不可能成为nesting的inneredge，故nu=ordinary crossing。

每个内部J的alpha也等同allG：若某incidentF的伙伴是0且另一伙伴同在左，指定的较大leftendpoint属于另一条G；若伙伴是L且另一伙伴同在右，指定较小rightendpoint也属于G；伙伴分居两侧本来无G指示。两条F同时触到一个J时partners为0,L，仍为split情形。K_JJ/WJ与其他源不变。唯一改变是内部W若连接边界，K_W由1变0，减去lambda中心化后留下-lambda。

比较allG词chi B chi，chi=h-W：边界chi连接内部H/J槽时权1，连接内部W时权-lambda，因为H-W交叉为0而W-W为lambda。与selected两F的内部source完全相同，另乘gamma²=s。chi两边界互配权tau chi²=1+lambda，扣掉这部分，补上真正selected直接互配的tauB：sigma1(B)=tauB+s[tau(chi B chi)−(1+lambda)tauB]。trace cyclic给rho1。所有比较是有限图和有限source derivatives；一般性不靠数值拟合。normal covariant representation使polynomial恒等式延至D。

## 3. 不同观测代数看到的平稳性

对B∈D_hx=W*(h,x)，mixed-q配对中每个含单W的项为0；WW两端相邻，其唯一配对不cross，tau(WWB)=lambda tauB。所以

sigma1(B)=tau(f(h)B)，f(h)=1-s+s h²。

f(h)严格正且与h对易；因此该restriction在alpha_t=Ad e^(ith)下stationary，且f有界。进一步B只含x时，两个H也相邻且独立于x种类，tau(h²B)=tauB，故sigma1|W*(x)=tau|W*(x)。整个x的谱测度仍是单位方差q^4-Hermite律，不只是前几个矩。另一方面，整个D的rho1通常不与h对易，精确

tau([h,rho1]†[h,rho1])=2q lambda(1+lambda)(1+q)(1-q²)>0。

推导可在h,u两色qGaussian空间旋转：chi/√(1+lambda)=z，h=chi/(1+lambda)+√[lambda/(1+lambda)]y，y,z单位正交。tau([y,z²]†[y,z²])=2(1+q)(1-q²)，由三种zz配对得出。完整state不stationary，而h,x restrictionstationary。这里平稳性指特定态相对于h-flow及其观测代数，不指定dS温度或Gibbs态。

## 4. 实际曲率的全固定时间关联

前阶段a_N=A0,N/(C_m√v_N)=(1-P_H)x_N的精确identity和polynomial-density桥，适用于M外置、内部整个covariant块。因此每个固定实t：

lim_N phi_N[M_N a_N(t)a_N M_N†]=C1(t)=tau[f(h)e^(ith)x e^(-ith)x]。

这是先实际finiteN exact pinching、后DS的极限；不是以finitewindow替换exact pinching。M制备归一为1。令r=q²、mu_q为单位方差q-Hermite能量measure，K_r(E,E')=Σr^n H_n(E)H_n(E')/[n]_q!（已有Mehler核），则

C1(t)=∫f(E)e^{it(E−E')}K_r(E,E') mu_q(dE)mu_q(dE')。

f>0、K>0且总质量1，故为正spectral measure的特征函数。simultaneous E,E'反号给C1实偶；不要求C1逐点非负。曲率等时谱相同并不决定这个two-energy measure，因为选态改变能量权f(E)。令C0为f=1的迹态核，有

−C1''(0)=sigma1([h,x]†[h,x])=2(1-q²)+s(1+q)(1-q²)²。

Mehler条件矩E[E'|E]=rE、E[E'²|E]=1−r²+r²E²，加sigma1(h²)=1+s(1+q)，即得此式。短时衰减率严格大于trace的2(1-q²)。

## 5. 可复查的物理参数计算与严格负号

q-Jacobi矩阵T offdiag√[n]_q，D_r=diag(r^n)。C0=<0|e^(itT)D_r e^(-itT)|0>，利用T²|0>=|0>+√(1+q)|2>得到C1=C0+s√(1+q)<2|e^(itT)D_r e^(-itT)|0>。与直接q-Pochhammer角度谱积分比较。q=.0625,.25,.64，lambda=-logq/4，t=0,.5,1,2,4,6；L192/256与GL(192,128products)/(288,192products)分别核对，最大两表示绝对差7.33e−14。L误差的exact-arithmetic共同界4[1+s√(1+q)]e^z z^(L−2)/(L−2)!，z=2|t|/√(1-q)，要求z²≤L(L−1)，本轮所有设置满足。浮点roundoff、GL离散不是interval认证；保存cutoff相对/绝对差、phase condition、product尾。t0界None表示exact0。

q=1/4（lambda=log2/2）,t=2时独立有理证书给：
0.0081500723940 < C0(2) < 0.0081500738176；
−0.0044147491323 < C1(2) < −0.0044147477086。

权威端点在sign-certificate.json中为Fraction，以上是向外放宽的十进制数。证书由monic Hermite递推精确计算Delta=E−E'的0..40阶偶矩，cos展开到40阶；|Delta|≤8/√3<5，所以余项≤10^42/42!≈7.1174e−10。正单位质量允许积分后保留此界。符号比较完全为有理数，不含float/GL/eigensolver。由此确认同一实际选态曲率在该极限时间有负关联，而trace核正。这里不是负概率，不是有限G混沌结论。

## 6. 已封存证据与下一步

exact-checks：两有理点，全部1093个长度≤6 covariant词共2186density比较；2953边界pairing fullsource指数；1022 reduced-hx词；26纯j矩；短时comm和full-density nonstationarity精确校准。curvature-dynamics：三物理q/18时点/两维数/两独立积分；图已保存。sign-certificate独立认证正负号。三个driver首次完成均passed，成功后不重跑。

下一步把可计算态接口用于高bilocal后代；优先核哪些已有rank/charge与fulltrace尾可经轨道density转入single-selected covariant块。若需要更复杂态，可先测试双selected制备的density，但不要自动沿用本轮仅两边界的几何理由。完整空间/质量字典和finiteG仍为长期目标，旧fulltrace完整bilocal四点保持已完成。
