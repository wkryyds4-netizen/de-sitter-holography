# 中央窗口内真实A范数、输入尾与core强展开

2026-11-06。固定标准CAR、完整overlapping complex Gaussian p-body、同J相邻charge、annealed normalized HS、ħ=1。m=N/2+xp，|x|有界，λ_N=p²/N在正紧集；共同w_N使C_m有界。正式公式见 `manuscript/susskind-memory-central-norm/main.tex`。

## 一般charge的关键变化

n=N−m+1，d=m−p，x_a=a(N+1−a)。已知精确传递为R_a=C_a^−/C_a=d²n/[m(mn−x_a)]，||Rem X_r||²=(m/n)(mn−x_r)||X_r||²。它们在中央窗分别有界和≤m²。

χ是上charge neutral Racah核，ξ的下参数为(1−m,m−N,−N−1)。带电T的两个分量给ξ=T_a(b)−b[T_a(b−1)−T_a(b)]/(N+1−2b)，不使用半填充乘积奇偶选择规则。逐项有ξ_j=m c_j/(m−j)，而χ^-_j=m(n−j)c_j/[n(m−j)]，故χ^-=χ+(n−m)(ξ−χ)/n。中央窗n−m可为O(p)，不能换成半填充值1。

令g=d²/m，h0=2p−m，η=p²/m，P=g²x_ax_b/[(mn−x_a)(mn−x_b)]，cf=(m−1)n(N+1)/N，D=ξ−χ，则实际同J pair norm为

f_ab=η²(1+χ)+P(mn/cf−1−χ)+[(g²−P)(n−m)/n+2h0g]D≥0。

||W2A||²/C_m²=Σω_aω_b f_ab，B cross的pair系数为η(1+χ)+gD。四个非半填充完整CAR系统31466个Gaussian坐标对核对(B,O)Gram，1331个有理pair另核contiguous、norm重组与正性。

## 一致相消和输入尾

中央χ summand的log二阶比半填充多4λ²x²j/p²。将其与精确f的有理系数组合，p²、p项Poisson生成函数求和均0，有限series项4λ²q，unsigned项4λ(1+λ)。对u,v,j≤ceil(log p) Taylor余项受控，窗外j有factorial envelope。因此f_(p−u,p−v)→4λ+4λ²(1+q)，在对数输入窗及中央x紧集一致。

输入ω_(ν−1)/ω_ν的精确式在中央窗≤K/(u+1)，固定u趋λ/(u+1)，得到Poisson(λ)输入deficit及factorial尾。窗外粗界f≤Kp4由相邻charge norm transfer和Cauchy–Schwarz给出。任何固定多项式输入权重都可保留，特别Σ(1+u+v)²||W2 A_ab;r||²≤KC_m²。

定义D_N=Σ(p−a)h_a，L_in=W2{h,D_N−bar u h}；L_in不是早先名为L的三次memory系数。得到Gram(B,L_in,O)/C_m²的极限

[[1+q,0,2λ(1−q)],[0,2λ(1+q),0],[2λ(1−q),0,4λ+4λ²(1+q)]]。

扣B的Y范数为4λ+16λ²q/(1+q)。本轮未把含S_N的半填充四方向Gram直接推广。

## 真实滤波器的强替换

上一轮六阶output矩配合传递界给Σ_(S>R)||W2A_ab;r||²≤KC_m²p4/R6。gA/C_m−γ≤K(1+u+v+S)/p，γ=(3−4q+q²)/2。输入/输出块正交给误差平方≤KC_m4[p^-2+R²/p²+p4/R6]。取R=p^(3/4)，得到ΣgA W2A_ab;r=C_mγO+O_L²(p^-1/4)，中央紧窗一致。

一般rank比例的gH有限前因子已另符号核对。以factorial输入尾和六阶输出矩控制该部分，完整第二chaos core为

W2Ecal=C_m[2p(1−q)²+Kbar+4λx(1−q)²]B+C_mγO+8λC_m(q²−q)L_in+8λC_mq²S_N+o_L²(1)，

Kbar=λ(1−10q+9q²)−16λ²q²。增长反项使用实际C_m、λ_N、q_N、x_N；总误差不套用gA单项速率。

36组全p+1秩和完整terminating series的80/120位检查通过，九个Gram元素及加权A范数在65存储位一致。最后λ.5,p128,x2的A范数为3.07899090506715，极限3.13533528323661；加权范数15.6459435799004，极限15.6766764161831，实际偏差与相消条件数保留，不将精度一致当DS收敛。

下一步先补S_N相关Gram，特别S–O的一般charge公式，避免半填充parity捷径；随后energy insertions、全charge尾、finite-memory及体介子字典。六页PDF真实双遍编译与逐页核版，全部仅本地保存。
