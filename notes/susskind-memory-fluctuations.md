# 能量分辨的记忆背景与首对涨落系数

2026-11-06，本地。[六页PDF](../output/pdf/de-sitter-susskind-memory-fluctuations.pdf)、[LaTeX](../manuscript/susskind-memory-fluctuations/main.tex)、[来源](../literature/susskind-memory-fluctuations-sources.md)。本轮无commit/push/ZIP。

## 规范与本轮结果

标准CAR，完整overlap complex Gaussian p-body H=Hs，w=4^p(Js²N/p²)(p!/N^p)²，全β=μ=0归一迹、ℏ=1。先N,p→∞且p²/N→λ>0，q=e^(−4λ)、σ²=Js²/λ，之后才取小λ。memory时间s、u与外部中心时间T分开。

上一轮完成的是p〈Mhat(s)〉→4Js²G_q(s)。本轮补出其能量分辨背景，并求mMhat首两个系数在明确守恒扣除后的完整双能量谱矩阵。对完整有限s,u的涨落尚未证明收敛；本轮的两系数极限先在有限N求导再取DS，不交换这两个操作。

## 1. 平均值之外仍有随能量变化的守恒背景

设ν_q是方差σ²的q-Hermite谱，μ_r(dxdy)=ν_q(dx)P_r(x,dy)是已有Mehler测度。BKS已给P_r H_n=r^nH_n、P_r自伴、P_rP_t=P_(rt)。这里r是能量转移核参数，不是实际memory时间。

对任意有界连续f，上一轮memory-force两个误差界加插f(H)只多||f||∞，故

\[
p\langle f(H)^\dagger\widehat M(s)\rangle
\to4J_s^2\int\bar f(x)g_s(x)\,d\nu_q(x),\qquad
g_s(x)=e^{-isx}P_q[e^{is\cdot}](x).
\]

正force双能量测度的弱收敛给此式；charge窗/尾沿用原界。这是完整固定s的能量测试回归关系，不是Mhat强算符收敛证明，也不是有限样本commutant条件期望的完整构造。

由P_q²=P_(q²)精确得到

\[
V_{bg}(s)=\int|g_s-G_q(s)|^2d\nu_q=G_{q^2}(s)-G_q(s)^2\ge0.
\]

条件前两矩给g_s=1+is(q−1)x−s²[(1−q)²x²+(1−q²)σ²]/2+O(s³)，所以Vbg=σ²(1−q)²s²+O(s⁴)，固定λ时足够小非零s严格正。

令X_N=pMhat(s)，以g_s−G_q作有界能量test并用HS Cauchy–Schwarz：

\[
\liminf\|X_N-\langle X_N\rangle1\|_2^2\ge16J_s^4V_{bg}(s).
\]

中央窗2m/N→1、||Mhat||₂²≤〈G²〉=O(p^−2)及指数charge尾进一步给

\[
\liminf\frac{\|Q_f\widehat M-\langle Q_f\widehat M\rangle1\|_2^2}{p^2}
\ge\frac{4J_s^4}{\lambda^2}V_{bg}(s).
\]

因此在m∼p²的集体尺度只扣一个scalar mean会发散；必须处理能量函数背景。尚未证明只扣领先(4Js²/p)g_s(H)就能得到紧的涨落，次领先守恒项也可能保留。

## 2. 精确扣除零阶与一次守恒多项式

固定m、h=H_m、r=p/m、J_hB={h,B}/2。上一轮G=A/m−r²h²、mZ1=T3+2rJ_hA−mr³h³，Mhat(s)=G+isZ1+O(s²)。

沿旧规范C_m=〈h²〉_m、a_m=〈A〉_m、c_m=〈(A−a_m)(h²−C_m)〉/〈(h²−C_m)²〉，Y=A−a_m−c_m(h²−C_m)。W_j为普通coupling chaos，

- a=−(1−q)²/(1+q²)，b=2λ(1−q)(1+q³)/[(1+q²)(1+q+q²)]；
- R=W3T3−aW3(J_hA)−bW3h³；
- F=W1T3−tbar h，B_S=W1(J_hA)−sbar h，B_Z=W1h³−zbar h；bar均为相应linear-chaos到h的投影系数；
- P_m=(ac_m+b)h³+[tbar−(ac_m+b)zbar]h。

旧精确T3−P_m=R+aJ_hY+F−aB_S−bB_Z。定义

\[
C_{0,m}=a_m+c_m(h^2-C_m)-mr^2h^2,
\]

\[
C_{1,m}=P_m+2r[a_mh+c_m(h^3-C_mh)]-mr^3h^3.
\]

二者逐样本与h对易，且

\[
X_{0,N}=mG-C_{0,m}=Y,\quad
X_{1,N}=mZ_1-C_{1,m}=R+(a+2r)J_hY+F-aB_S-bB_Z.
\]

用旧rank tail估计B_S/B_Z=O_L²(1/p)、固定λ下J_hY有界二阶矩，极限一次系数为L=R+aJ_HY+F。L只是新记忆系数的符号，不是旧Liouvillian；本文Liouvillian显写[H,·]。有限N不宣称扣完全部commutant。

## 3. 首对系数的正谱矩阵

旧F/Y/R交叉双能量测度均为0，对角依次VF μ_q、W μ_(q²)、VR μ_(q³)。令e=(x+y)/2，则

\[
d\nu_{(Y,L)}=
\begin{pmatrix}
W d\mu_{q^2}&aeW d\mu_{q^2}\\
aeW d\mu_{q^2}&V_Fd\mu_q+V_Rd\mu_{q^3}+a^2e^2W d\mu_{q^2}
\end{pmatrix}.
\]

W=4Js⁴[1/λ+4q/(1+q)]、VF=64Js⁶q²；VR=σ⁶(λA1+λ²A2)/D²，D=(1+q²)(1+q+q²)，A1=4(1−q²)(1+q²)(1+q+q²)²，A2=16q²(1+q)(1+q²)(1+q+q²)(2q²−q+2)。相对ν_q⊗ν_q的密度矩阵行列式WK_(q²)[VF K_q+VR K_(q³)]≥0。

普通中心时间交叉〈Y†L(T)〉=0：变换(x,y)→(−y,−x)保留频差但翻转e。能量分辨cross不为零；再乘e积分得aWσ²(1+q²)/2。

J2(T)=∫e²exp[iT(y−x)]dμ_(q²)，则CL(T)=VF G_q(T)+VR G_(q³)(T)+a²WJ2(T)。首阶memory多项式Xi(s)=Y+isL的协方差为WG_(q²)(T)+su CL(T)。对真实memory只确定Taylor系数(0,0),(1,0),(0,1),(1,1)，不能把这个多项式当完整有限s,u解。

## 4. 居中自然去掉时间后代

raw M3扣P_m仍含(3/2)[H,Y]；memory tangent L没有这一项。于是

\[
N_L=V_R+V_F+a^2W\sigma^2(1+q^2)/2,
\]

\[
D_L=2\sigma^2[V_R(1-q^3)+V_F(1-q)]+a^2W\sigma^4(1+q)(1-q^4)/2.
\]

小λ：λW→4Js⁴、λNL→32Js⁶、DL/NL→24Js²；F范数占比2λ+o(λ)，能量后代占比8λ³+o(λ³)。归一L中心时间核趋sech⁶(2JsT)，Y趋sech⁴。这是既有q核在实际memory系数上的应用，不是重新发明核，也不把d标签当体质量。

## 5. 验证、失败与下一步

checks.json：3有限系统各5个有理复矩阵恒等式，通过直接V†K^nV的renewal与另一套Gaussian polarization的cubic chaos交叉核对；28个条件Hermite精确多项式；三个Laurent系数4/32/768；谱矩阵行列式恒等式。四λ Jacobi/+12，λ.25/.1独立128/192节点Mehler、直接conditional variance、ordinary cross0及energy cross非零均核对。

Js1、s=.5的Vbg依λ=.25/.1/.025/.01分别为.124133785950749、.047414703352544、.010685786898544、.004165911045489。一次系数中R的等时比例分别为.81169764、.88637262、.95765673、.98132544。有限λ偏差不当浮点误差。

首跑独立对照失败：Mehler代码的e² descendant遗漏K_(q²)，λ.25的归一t0错误为.9910814。初版与诊断保留，补上该因子后通过，不改容差。解析Jacobi截断界与舍入refinement分开；后者非区间认证。

下一步：先求二次居中记忆Z2及其G²首返项，明确第四singlet与守恒多项式，计算与Y/L的mixed谱；再控制两个有限memory插入的余项和涨落尺度，而非反复做单核扫描。已有能量回归不能替代该两插入证明。继续接体介子内部核与空间字典，完整goal active。
