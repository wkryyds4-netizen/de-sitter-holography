# 同一 J 的全部固定 word 与有限时间构造

2026-11-07。本轮将 joint-fourth 的两弦规则扩为任意固定阶，并给实际有限N统一时间尾。完整证明见 [LaTeX](../manuscript/susskind-slater-joint-words/main.tex)、[九页PDF](../output/pdf/de-sitter-susskind-slater-joint-words.pdf)，[原文核查](../literature/susskind-slater-joint-word-sources.md)。

## 对象

沿用标准CAR、半填充m=a=N/2、固定Slater O、原overlapping p-body Hermitian Gaussian J。φ_N=E_J〈O|·|O〉非stationary。L=C(m,p)²、D=C(m,p)C(m+p,p)，M†=原纯O→E块/√L、h=H/√(wD)，q=e^(−4λ)、γ=q^(1/4)，p²/N→固定λ>0，ℏ=1/string能量，字长固定后取极限。h无量纲，真实H时间乘√(wD)→J_s/√λ。

## 全word规则

每张Gaussian pairing中HH边标G，其余边S由M端点固定方向；++/−− edge不相容，creation先作用（从右读）为forward。

\[
\lim_N\varphi_N(w)=\sum_{\pi\;\mathrm{compatible,forward}}
\gamma^{b(\pi)}q^{\nu(\pi)},\qquad
\nu=\#\mathrm{crossings}+\#\{\mathrm{nested\ pairs:\ inner\ S}\}.
\]

b为HM混合边数，奇数总阶为零。G完全处于另一边内的nesting不惩罚；SS任何active overlap惩罚q。这个定理现在由一般k的解析证明给出，不再只作为两弦猜想。

## 证明核心：局部多项式和归一化缺陷

每条边右端X=P†Q、左端X†。G的P/Q任意；forward S要求P⊂E,Q⊂O；backward S相反。每site membership S⊂{P_e,Q_e}，按实际局部CAR次序算初始空/占据矩阵元，再乘(−1)^(Σ_cross t_e t_f)，t_e是site属于该edge两集合的个数。得到f_o(S)∈{−1,0,1}，F_o=Σ f_o(S)z_S。原始配对精确等于[Π z_i^p]F_0^m F_1^m。四个thin crossings的4p²为偶，BNR精确重排符号给剩余site符号；高重交集尚未丢弃。

除L^k后，设多重membership的个数n_(o,S)，d_i=Σ_(o,S含i)n_(o,S)，singleton=p−d_i，U_o=kp−Σ_(correct i in o)d_i+Σ_S n_(o,S)。精确每项为

\[
\frac{(m)_{U_0}(m)_{U_1}}{(m)_p^{2k}}
\prod_i(p)_{d_i}\prod_{o,S}\frac{f_o(S)^{n_{o,S}}}{n_{o,S}!}.
\]

绝对值被C_(k,λ)Π α_(S,N)^n/n!支配，α=m(p/m)^|S|。总activity有界，|S|≥3 activity为O(p³/N²)→0；总defect数>R由指数级数尾控制。因此可以先固定defect数取极限，再求和，包括所有signed cancellations。Backward S强迫一个正确singleton数为零、d_i=p，至少p个defect，受同一尾控制趋零，覆盖包含G的图。

双membership总系数：自身G:1、自身S:0、disjoint:2、crossed:0、nested innerG:2、nested innerS:0。Falling-factorial背景e^(−2λk(k−1))与双cluster指数相乘，得raw/L^k→e^(2λk_G)q^ν。再除(D/L)^(n_H/2)且n_H=2k_G+b，即得γ^b。正文有逐局部情形和系数证明。

## 实际有限时间，而非仅Taylor系数

将h,M写成同一组独立实Gaussian x_α的矩阵系数A_α,B_α。精确ΣA†A=ΣAA†=I；ΣB†B及ΣBB†在E占据x的基向量上分别取C(x,p)²/L、C(m−x,p)²/L（顺序按M或M†互换），均≤1。

每张Wick图从右读，首次标签用column map打开noise tensor slot，第二次用row map关闭；crossing只需置换打开的noise slots，范数1。故任意有限N、任意2r字的每图算符范数≤1，|φ_N(w)|≤(2r−1)!!。字v长d时||vΩ||_annealed≤F_d=√((2d−1)!!)。这避免了固定k cluster估计随k²增长而无法直接求时间尾的问题。

全word极限正性给GNS；任意polynomial P degree d有||h^nPΩ||≤C_P F_(n+d)，所以全polynomial core为entire analytic vectors，Nelson定理给h本质自伴。任意固定polynomial插入的有序指数关联，设T=Σ|z_j|、d=总插入度、C=系数ℓ¹范数乘积，总指数次数≥R的尾≤CΣ_(n≥R) T^n F_(n+d)/n!，在ρ_R=T√(2(R+d)+1)/(R+1)<1时≤C T^R F_(R+d)/(R!(1−ρ_R))。先N极限再R极限，紧复时间窗一致。插入的graph-norm收敛由同一增加固定degree的界保证。覆盖真实有限时间模型关联，但不覆盖随N增长的时间或宣称quenched收敛。

## 六阶校准与记录

φh⁶=5+6q+3q²+q³，φM³M†³=6q³。g=h−γ(M+M†)、s=√q，

\[
\frac{\varphi(g^6)}{(1-s)^3}=5+6s+14q+10sq+13q^2+7sq^2+6q^3.
\]

在γ=1/2、2/3两个有理点，degree≤2 Gram秩7、degree≤3秩15，正pivots已保存。这只是两点精确核对，不据此断言一般15维或二叉树表示。

- site-checks.json：1,461,120全CAR标签/局部式比较，1,215精确系数提取，1,814双cluster全图表核对；N4p1与N6p1分别有8、72个负贡献，保留并通过。
- word-checks.json：729六字的精确γ多项式、g六矩、两有理点Gram。
- time-bound-checks.json：348精确finite covariance计数、9个70位精度尾诊断；后者非interval证书。

旧stage/driver未重跑。新阶段本地保存，不Git/不issue。下一目标：把Gram与左乘法转为可高效计算的全层transfer/波函数表示，计算state preparation和谱核，然后明确stationary观察者态及bilocal/体生成元映射。当前全joint dynamics限于annealed Slater模型；完整de Sitter字典和finite-G仍未完成。
