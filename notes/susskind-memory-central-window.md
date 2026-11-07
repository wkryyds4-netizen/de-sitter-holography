# 中央电荷窗口：秩尾与守恒反项

2026-11-06。标准CAR、完整overlapping complex Gaussian p-body、同一样本各charge使用共同耦合方差w_N、固定charge归一化迹与annealed HS范数、ħ=1。λ_N=p²/N在正紧集，m=N/2+xp=N/2+z√N；x或z在固定紧集。详细正式推导见 `manuscript/susskind-memory-central-window/main.tex`。

## 六阶矩的一致扩展

将旧Racah三角递推直接改为M=N/2+z√N，而不是沿用bounded shift的证书。ε=N^(-1/2)，a=A/ε、b=B/ε；缩放系数的分母变成(1+2zε−2kε²)(1−2zε−2kε²)(1+(1−k)ε²)。在QQ[A,B,z]展开至ε12，第一至第六阶的21个负幂全部相消，常数仍为Touchard Tj(2AB)。抽掉已知ε幂后分母首项非零，精确有理函数在参数紧集解析，故这是紧集一致估计。

X=x_(a+b)−x_r，x_j=j(N+1−j)，S=a+b−r≥0。X≥S(N+1−4p)给中央窗内sup Eπ S6≤K，进而π(S>R)≤KR^-6。这是独立等向rank输入的正乘积律；实际memory还须接相邻charge相干范数及输入尾，不能直接宣布完整Gram已扩展。

四完整有限分布共28个0至6阶矩精确一致；40组中央窗口有限参数以Fraction记录，并核M↔N−M对称。p128、λ1/4、z2、a=b=p时六阶proxy为22.4617786805，极限24.046875，保留实际偏差。

## 能量尺度与滤波器有不同的有限修正

共同w_N下，精确C_m=w_N binom(m,p)binom(N−m+p,p)。令C0=C_(N/2)，固定x时

p(C_m/C0−1)=4λ²x(1−x)+O(1/p)。

由两个binomial比值的逐因子对数展开得证；精确反射是C_m=C_(N+p−m)，即x↔1−x。q_m=exp[−p²N/(m(N−m))]是DS参数，log q_m=−4λ−16λ³x²/p²+O(p^-4)。因此谱形的变化更小，仍不能在pC_m反项中换成C0。30个精确比值及反射检验已保存。

从旧两charge精确五字滤波器直接展开，固定u,v,S时

gH/C_m=2p(1−q)²+K0+4λx(1−q)²+8λ(q²−q)(u+v)+8λq²S+O(1/p)，

K0=λ(1−10q+9q²+16λq−32λq²)。新增4λx(1−q)²与rank无关，Gaussian居中后乘C_m(h_m²−C_mI)，属于守恒背景。gA四个rational prefactor至1/p仍为(3/2,−1,−1,1/2)+(λ/p)(8,−5/2,−5/2,3)，固定deficit领先γ=(3−4q+q²)/2未变。

用独立有限五字实现核60个精确gH/gA恒等式；30组中央减半填充差以80/120位复核。最后一组p256、x1、(u,v,S)=(2,1,3)为0.39841612746439，极限0.39957640089373。精确Fraction和有限差值保留，没有拟合消除偏差。

## 下一步

以本轮六阶矩和修正后的守恒背景为输入，推实际charged norm transfer、A的input-rank尾及marked overlap的中央窗一致性，然后重组Gram与投影。固定L窗口之外还需可积控制才能做全charge迹。之后是两插入finite-memory余项与体介子字典。本轮三JSON、六页双遍编译/95dpi核版仅本地保存。
