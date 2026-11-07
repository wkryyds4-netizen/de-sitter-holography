# 三阶一次残量的动态矩与保留谱权重

2026-11-06，本地。[LaTeX](../manuscript/susskind-cubic-dynamics/main.tex)、[PDF](../output/pdf/de-sitter-susskind-cubic-dynamics.pdf)。无 Git 提交/推送/ZIP。

## 物理推进

沿 Susskind singlet 主线，上一轮得到 F=W₁T₃−tbar(Qf)H 的 DS 正范数，本轮进一步证明它在极限中仍有动态谱权重。

标准 CAR、完整重叠 complex p-body Gaussian 系综、β=μchem=0、Hs 的 string 时间、ℏ=1。先固定 λ=p²/N>0 取 DS，再小 λ。q=exp(−4λ)，C=Js²/λ：

- 〈F²〉DS=64Js⁶q²，沿用前轮。
- 本轮 〈FL²F〉DS=128Js⁸q²(1−q)/λ。
- 归一第二频率矩 μ₂,F=2Js²(1−q)/λ，随后小 λ→8Js²。

已有 X、R 的小 λ 曲率为16Js²、24Js²，因此 F 的全阶 q 探针等价值得直接检验；这轮没有把一个曲率当成完整核或体质量。

## 精确有限电荷公式

以下 p≤min(m,N−m)，其他电荷只在 full trace 的指数小尾控制中使用。D=binom(N,m)，C_r=binom(m,r)binom(N−m+r,r)，d_a=binom(N,a)²−binom(N,a−1)²。

ω_a^(r)=d_a binom(m−a,m−r)binom(N−r−a,m−r)/(D C_r)。实际 C_m=w C_p，v_a=C_m ω_a^(p)。这是已有嵌入系综权重，w不是λ。

H=ΣH_a，各 Gaussian 不可约通道独立。纯 rank 协方差 E H_a Y_b H_a=v_a g_ab Y_b，由已知完整 r-body covariance 的三角式唯一确定：

Σ_(a≤r) ω_a^(r) g_ab=binom(m−b,r)binom(N−m+r−b,r)/C_r。

显式解 g_ab=₄F₃(−a,a−N−1,−b,b−N−1;−m,m−N,−N−1;1)，终止于 min(a,b)。它是已有 Racah 多项式；本轮给当前参数匹配和协方差归一。g对称，g_0b=1；归一协方差是 unital trace-preserving CP map，故 |g_ab|≤1。

证明的中间恒等式：z_k(a)=(−a)_k(a−N−1)_k，Σω_a^(r)z_k(a)=(−r)_k(m−N)_k(−N−1)_k/(m−N−r)_k。把权重写成 very-well-poised 形式，乘 z 后移 a=k+j，使用 Rogers–Dougall 终止求和；再用 Pfaff–Saalschutz，把 g 的 ω 平均还原成已知 binomial 协方差。详见稿第2页。这些标准恒等式不作原创声明。

取上一轮精确 t_a，f_a=t_a−Σωt，F=Σf_aH_a，V_F=Σv_af_a²。完整 Gaussian Wick 给

B20=Etr F²H²=C_m V_F+Σv_av_b f_af_b g_ab；
B11=Etr FHFH=Σv_av_b f_a² g_ab；
M2=2(B20−B11)=Σv_av_b(f_a−f_b)²(1−g_ab)≥0。

B20的第二项就是同J交叉收缩。N5p2m3、w1时 C18，V180，B203420，B111800，M23240，额外交叉180。这个 ensemble 平均不同于上轮单个确定Hamiltonian的witness。

## DS 与谱尾

固定 rank deficits s=p−a,u=p−b，Racah 第 k 项相邻比趋−4λ/(k+1)，全范围模可由16λ/(k+1)支配，所以 g_(p−s,p−u)→q。先证明纯rank而不是只用平均κ。

前轮已证 s→Poisson(λ)、c_(p−s)−c_p→8λq s，其中 c=t/C_m，并有统一 factorial tail 和线性系数界。结合 |g|≤1，可以对带系数平方的双rank和取极限，M2/C_m⁴→128λ³q²(1−q)。分别得到 B20/(C_m V_F)→1、B11/(C_m V_F)→q；同J交叉项在完成极限与中心化后才趋0。

full trace：先取 |m−N/2|≤N^(2/3)，该窗内密度、C_m/C_(N/2)和上式一致收敛。补集 binomial 指数尾压过 |t_a|≤K N Cmax 的多项式增长，Cmax/Ccenter→1。不能用太宽的 N^(3/4) 窗口直接声称 C_m 一致趋C；当前正文已采用2/3。

为了排除“极少高频权重逃到无穷，留下假曲率”，对 Hermitian Gaussian Y 且 E Y²=s²I 用 Gaussian 分部积分和 Young 不等式得 Etr Y^(2k)≤(2k−1)!! s^(2k)。F也是Gaussian且 E F²=V_F I，故归一频率测度满足 μ4≤16 sqrt(315) C_m²。不需要 H 与 F 独立。

因此每个 DS 弱子列极限均继承正的第二矩，不能是δ0。Paley–Zygmund还给 P(|ω|≥sqrt(C(1−q)))≥(1−q)²/[16sqrt315]，对固定λ成立。这一保守下界没有排除额外静态atom，也未证明完整谱唯一。

## 证据与继续位置

1. check_cubic_dynamics.py：1515个精确有理三角反解与Racah比较，928个factorial moment，5个直接有限CAR动态矩。最大Gaussian整数矩阵元48、范数和9216。
2. check_cubic_dynamics_moments.py：3个系统全部B20/B11的三配对分别核查，6组高精度加速和与精确双rank和比较。它针对同J交叉项，而不仅核最终相减。
3. cubic_linear_dynamics.py：λ=.1,.25,.5，p16/32/64/128，全rank和、全终止级数，80/120位，保存65位归一矩一致。记录相消系数；非舍入区间认证。数值约定有限sector C_m=Js²/λ，Js1。
4. p128曲率相对DS仍偏−.63674%、−.45805%、−.25000%。保留有限p偏差，不用稳定精度冒充到达极限。

下一步直接计算〈F H^u F H^v〉更高阶并核全阶收缩。当前只证明纯rank两弦 crossing，尚不自动证明多弦因子化；同J F–H pairing 在高阶的消去仍须论证。随后求F与R/S/Z实时cross，完成raw M3谱，再研究X/R/F混合高点及体平移字典。
