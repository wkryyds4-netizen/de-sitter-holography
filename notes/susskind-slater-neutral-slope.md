# 实际中性斜率、协变子代数的放大表示与有限分辨率去相干

2026-11-08。本轮从slater-covariant继续，固定半填充m=N/2、同J overlapping complex CAR、λ_N=p²/N→λ>0；φ_N=E〈O|·|O〉。物理H_N方差C_m=w_ND_N，h_N=H_N/√C_m，无量纲时间共轭h_N。新stage的t_N是归一实际中性斜率，不与时间变量混淆；旧笔记中的h=H_m是有量纲Hamiltonian。

## 1. 从已有EGUE covariance到实际斜率

Φ_m(X)=E_{H'}H'_m X H'_m（独立H'只用于定义确定性映射，不是新matter），Φ_m/C_m=μ_N(R)。沿用Kota公式

μ_N(r)=binom(m−r,p)binom(N−m+p−r,p)/[binom(m,p)binom(N−m+p,p)]。

H的归一rank权重为ω_N(r)=[binom(N,r)²−binom(N,r−1)²] binom(m−r,m−p)binom(N−p−r,m−p)/[binom(N,m)D_N]，r=0,...,p。令μbar=Σωμ、κbar=C_m μbar，α=1−p/m、ε=1−α²。实际未pinch一次系数T_N=mεα[Φ_m(H_N)−κbar H_N]，实际exact-pinched B1=(1−P_H)T_N。归一

t_N=T_N/[8λ_N C_m^(3/2)q_N√λ_N]=f_N(R)h_N，
f_N(r)=mεα[μ_N(r)−μbar]/[8λ_Nq_N√λ_N]，q_N=e^(−4λ_N)。

注意T_N前减的是精确annealed regression；在最终1−P_H或1−D_T下，该H项都精确消失。

令d=p−r。相邻权重比精确为a²(1−b²)/(1−a²)·(m−r+1)/(d+1)·(N−p−r+1)/(N−m−r+1)，a=r/(N−r+1)、b=(r−1)/(N−r+2)。固定d极限λ/(d+1)；p≤N/8时≤8λ_+/(d+1)，故所有d矩统一可积，d→Poisson(λ)。这是rank缺额，区别于纯Slater particle-hole hop缺额Poisson(2λ)。

固定d，p[μ_N(p−d)−μ_N(p)]→4λq d；乘积差有界0≤p差≤6λ_Nd。结合mεα/p→2和精确回归，得Σω|f_N(r)−(p−r−λ_N)/√λ_N|²→0。差场是Hermitian linearGaussian equivariant，scalar covariance等于此平方差。因此沿上一轮Wick row/column界，t_N可在所有固定selected mixed words和紧时间关联中强替换为u_N；这不是只在trace上的边缘矩比较。选中纯p-hop为rank p，精确φ_N(t_N M_N†)=sqrt(L_N/D_N) f_N(p)→−√λγ。

## 2. 全协变表示是标准qFock的放大

使用已证H/U/P正Gram，定义K_n=H_n∩ker a_H∩ker a_U。有限grade内，K_n是ran a_H†+ran a_U†的正交补。由qCCR，把任意两个ordinary创建词的内积正规序，得到标准两色qFock Gram乘primitive内积。不同primitive grade的词即使总grade相同，也因剩余annihilator杀primitive而正交。逐grade分解补空间穷尽全H，故

H≃F_q(R²)⊗K，K=⊕K_n，h=s_H⊗I、u=s_U⊗I。

这是上一轮指定态GNS之外的新信息：当前完整Slater表示对D=C*(h,u)的限制本身就是标准表示的放大。因此D''≅Γ_q(R²)⊗I。它是diffuse有限因子，无非零compact；结合A+=D+K(H)得到A+/K(H)≅C*_q(R²)。这里终于有了full quotient识别所需的faithful amplification证明，不能把此结论倒写成上一轮已经证明。

所有有限grade向量ψ在该分解下可写Σ_j P_j(h,u)Ω⊗κ_j（κ正交），故其D态有有界正密度ρψ=Σ_j P_jP_j†；σψ(B)=τ(ρψ B)。一般向量给normal L1密度，未声称一律bounded。

## 3. 一个selected激发的精确密度

c=−√λγ，v=γh+cu，r_P=1−s(1+λ)>0。ζ=M†Ω−vΩ同时被a_H、a_U杀掉，norm²=r_P，是另一个primitive copy。于是对每个B∈D（并normal延拓至D''）

φ(M B M†)=τ([v²+r_PI]B)，ρ1=v²+r_PI。

τ(ρ1)=1，r_PI≤ρ1≤[r_P+4s(1+λ)/(1−q)]I，后者恰为norm。r_P不是上一轮u相对hM的残量r_u=1−λ/(e^(2λ)−1)：这里投影方向相反、子空间不同。r_P→0而r_u→1/2并不矛盾。

## 4. 有限分辨率能量平均可接回实际场

D_T(X)=(2T)^−1∫_{−T}^T e^{ith}Xe^{−ith}dt。qMehler两能量核K_q和μ_q给

V_q(T)=τ[(D_Tu)²]=∫sinc²[T(E−E′)] K_q(E,E′)dμ_q(E)dμ_q(E′)。

对固定q，K_q连续有界、ρ_q是平方根边缘密度，故T V_q(T)→κ_u=π∫K_q(E,E)ρ_q(E)²dE>0。又V_q(T)=T^−1∫_0^(2T)(1−t/(2T)) C_u(t)dt。单selected态误差≤||ρ1||V_q(T)，每个有限gradeψ则≤||ρψ||V_q(T)。D_Tu一致有界，先在dense finitegrade证明，再延拓得完整Slater H上的strong D_Tu→0。

有限N定义同样的D_{N,T}，使用dimensionless h_N时间，并设b_{N,T}=(1−D_{N,T})t_N。这正是归一实际一次系数的有限时间平均去除，精确回归H项自动取消。对固定T，上一轮紧时间全word定理及本轮线性Gaussian桥给b_{N,T}→u−D_Tu。继而对任意固定selected词W，lim_{T→∞}lim_{N→∞}φ_N W(h_N,M_N,b_{N,T})=φ W(h,M,u)。同样允许固定外部实时插入。

顺序不能换成finiteN exact pinching。固定Jacobi截断L的T→∞会出现正plateau；不是continuum残余。数值报告该plateau、L比较及exact-arithmetic误差界4exp(x)x^L/L!，x=2TE_*；不enclose浮点roundoff。

## 5. 精确pinching已经能传到哪里，缺什么

仅协变h,u族：在central sector，正两能量测度ν_N由t_N定义，其矩/time极限为K_q μ_q²。这个极限没有E=E′质量。Portmanteau直接给||P_H t_N||_(2,m)²=ν_N(diagonal)→0，不交换时间极限。Gaussian高矩与pinching Schatten收缩再给每个固定S_r→0；Schur和tracial Hölder便使actual exact-pinched slope可用于所有纯协变固定words/time。

含M的finiteN情况仍缺一种uniform selected-state绝对连续性控制。若能对Y_N=(P_Ht_N)²证明E〈O|P_N†Y_NP_N|O〉≤C_P(E tr_m Y_N^r)^(1/r)，C_P独立N，则旧S_r小量足够。此处只列足够条件，不宣称对所有covariant正Y_N成立。极限D上的bounded density不能自动推广给N-dependent精确谱投影。下一步应直接核这个特殊Y_N或用实际算符结构给替代估计，然后J2/fullbilocal；不能把有限分辨率定理当作该缺口的完成。

## 6. 证据

rank-checks：20精确rank系统、12全rank DS行、80/120位比较，保存全部p+1概率与系数，有限误差完整保留。amplification：1022精确density恒等式、8个primitive向量、32个tower内积；有理γ,d独立测试点不在physical曲线。dephasing：3physicalq、36有限时间cutoff行、512/768 Jacobi、两套对角积分、独立时间积分差1.82e−13。无限维结论靠上述证明，非finite核维数外推。三个成功driver不重跑。仅本地保存，无Git/issue操作。
