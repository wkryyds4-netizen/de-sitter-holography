# 固定电荷梯度与实际 ordinary 全次数尾

2026-11-08。本轮接 central-ordinary 的逐阶强桥，使用相同 standard CAR、同一 Hermitian Gaussian p-body disorder、半填充及 p²/N 留在固定正紧区间的双标度。严格陈述与完整推导见 [LaTeX](../manuscript/susskind-sector-gradient-tail/main.tex)，阅读版见 [PDF](../output/pdf/de-sitter-susskind-sector-gradient-tail.pdf)。这里是模型内的有限 N 估计，不改变 finite-G 或 bulk 字典的状态。

## 结果及适用范围

对任意固定有限 Schatten 阶 s、固定 L0，在 |m−N/2|≤L0 p 的中央窗口，

G_(n,N;m)=C_m^(−n/2)[p Ncal−D−lambda_N Ncal(Ncal−1)]H^n|m

具有 K_s(n+1)^b L_s^n sqrt(n!) 的 N、n 一致包络。D 先作用全中性块族再限制扇区，C_m 是目标扇区固定常数，不移入 D。Gaussian number 不能替成普通外部次数 n。时间用 h=H_m/sqrt(C_m)、hbar=1 的无量纲单位。

因此实际 ordinary 时间场的尾由 sum_(n>K)(n+1)^b(L_s S)^n/sqrt(n!) 控制，在每个固定紧时间窗一致趋零。与上一阶段 m=N/2 的逐阶 same-copy 强桥相接后，可解除 ordinary 展开的次数截断；逐样本 exact pinching 及已控制的固定多项式制备也可接入。多个字段限于已定义的协变观测块和固定次数制备，不暗中扩展到所有非多项式交错制备。

## 一般估计的关键

1. 相邻电荷梯度 Y_i=[c_i,X] 同时保留 row/column 平方函数和 X†。归一 t/s 通道梯度的主项，平方函数收缩因子分别为 sqrt((m−1)/m) 与 sqrt((N−m)/(N−m+1))。中央区余项每步≤8/sqrt(N)倍邻域算符范数，所以 ell 步代价线性为 ell/sqrt(N)，不出现 K^ell。
2. 在分次 CAR 上，partial 与 D_g、L_g 的交换子均为 partial/2。于是梯度可先移到原始 Wick 块上，所有 np−T 缺额同时变为 np−1/2−T_g。normal-ordered 单项式的导数只删掉一个匹配字母，保持其原 p-set membership；不必在反复标记后追踪 sandwich 因子的表面导数。
3. 2r 次梯度矩含 r 个被成对求和的 gradient index。按 site 相等关系分组，A 个带 pure tag 位点各至少有两个不同底层 set 标签，B 个 gradient-only 位点至少一个，且 B≤r。保留 CAR 符号的受限 covariance/ancilla 收缩给 N^(A+B)(p/N)^(2A+B)=lambda^A p^B。开 2r 次根后，梯度代价为 sqrt(p)。k 次 repeated deficits 的包络为 sqrt(p)(n+1)^(2k+1)L^n sqrt(n!)。
4. 梯度与 H 幂相乘不能直接把 I_N tensor H 的辅助 trace 因子丢掉。对矩形 CP 映射使用 normalized-trace adjoint Phi^sharp=(d_source/d_target)sum V†(·)V，谱分解加两次 Holder 得由正反 Gram 控制的固定阶范数估计。由此得到 g_s(XY)≤K_s[g_(4s)(X)||Y||_(4s)+||X||_(4s)g_(4s)(Y)]，相邻维数比在中央区有界。
5. 对实际 A_v=p(Phi−Ctheta_v)Z_v，用旧精确 CP telescope 加上述梯度界。vp≤N/8 时系数总量 p²/N 有界；vp>N/8 时 p=O(v)、N=O(v²)，粗界仍只有次数多项式代价。旧 A_v、DA_v/p、LA_v/p 的证明原本已逐扇区，无须重做 fulltrace Riesz。
6. 先完全展开受迫 Jacobi recurrence，再作固定次数 Holder。A_vH^w 和 H^a(pr_ell) 只有两层乘积，范数阶数不随 degree 迭代翻倍。物理梯度与 Gaussian top-chaos 投影可交换，OU/Chebyshev 对行列数组仍适用，因此 eta_j 及其 D/p、L/p、gradient/sqrt(p) 都有阶乘包络。
7. 逆 Jacobi 系数差分保留 1+|X_N|，X_N=(m−N/2)/p。正是这个因子限制最终定理到 O(p) 中央窗口；没有声称整个 [N/3,2N/3] 区间都有同样常数。旧五项 ordinary 重组恒等式接回真实 H^n，得到统一尾。

## 验证与证据

上一轮成功输出 results/sector-cp-gradient/finite-car-checks.json 保持原字节：1032 梯度恒等式与216 unital/normalized-trace 恒等式精确通过，168 PSD与324 Schatten检查为浮点容差诊断。

新 [graded checker](../code/check_sector_graded_gradient.py) 首次运行成功，遍历 N=1..4 的340个matrix units、340个normal monomials：5008个分次导数交换子、840个至多三次重复缺额移位、340个Euler关系及2504个带符号删除关系精确通过，其中1252个导数非零。采用乘2的 D_g/L_g 矩阵以整数计算；没有冒充完整 Gaussian ensemble 数值模拟。

一般尾定理依赖正文的代数和范数证明；有限检查核符号与归一，不代替次数一致性证明。两个成功科学输出及所有父阶段文件不覆盖、不重跑。原来的全 Fock Riesz 证明保留；本轮替换的是它不能直接用于固定电荷的那一环。

## 下一步

ordinary charged 级数的实际中央尾已接通。下一首任务转为完整 selected bilocal 的真实 neutral/rank/filter 与 memory-return 逐扇区审核，确定旧 fulltrace 接口哪些可直接保留、哪些必须另证；随后才是空间/质量生成元和 finite-G。低频 Git 约定继续，本阶段仅本地保存。
