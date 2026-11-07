# 三阶两点谱闭合：带电协方差提升与真实时间后代

2026-11-06，本地阶段。正式 [LaTeX](../manuscript/susskind-filtered-spectrum/main.tex) 与 [七页 PDF](../output/pdf/de-sitter-susskind-filtered-spectrum.pdf)。沿 Susskind 静态补丁 → complex DSSYK∞ → 导数 singlet 主线；本轮无 Git 提交、推送或 ZIP。

## 主要结果

上一轮留下的能量分辨交叉测度 ν_FY、ν_FR，在固定正 λ 的双标度极限中均为零。结合已有 ν_RY=0，F/Y/R 的两点谱矩阵对角化。这里说的是两点，不是有限 N 统计独立，也没有推出混合高点独立。

标准 CAR、含重叠指标的 complex Gaussian p-body H_s，退火归一 trace、β=μchem=0、ℏ=1；L=[H,·]、J_H={H,·}/2。先 N,p→∞、λ=p²/N>0 固定，q=e^(−4λ)、σ²=J_s²/λ。下列小 λ 极限随后取，时间固定。

## 真正新增的接口：charged lift

在 sector m，h=H/√C=Σhν，y=F/√〈F²〉，s=p−ν。前轮精确权重比为

R₁(s)=ω_(p−s)^(p−1)/ω_(p−s)^p
=s(N−2p+s+1)(N−m+p)/[(m−p+1)p²]。

定义 yhat=√λ Σ[R₁(p−ν)−1]hν，则 ||y−yhat||₂→0。这次使用算符 L² 桥接，而不只是协方差总变差。因为 s→Pois(λ)，真实滤波趋(s−λ)/√λ，R₁→s/λ；factorial tail 和多项式界给 L² 支配。ω 加权 R₁ 和严格为1。

一条 yhat–h Wick 边等于 √λ(K_(p−1)−K_p)。不能仅在 sector m 替换后就对导数 mark 使用普通 chord 规则，因为 mark 也经过 m−1。

设 V^(p−1) 是在参考 m 归一到方差1的独立 embedded ensemble，定义

Vtilde=g(Q_f)V^(p−1)，g(n)=(n−p+1)/(m−p+1)。

g 把它提升为真正的 p-body 算符。rank ν 在 body r、charge n 的因子为 binom(n−ν,r−ν)，从而

g(n) binom(n−ν,p−1−ν)/binom(m−ν,p−1−ν)
=binom(n−ν,p−ν)/binom(m−ν,p−ν)。

所以任意允许的 m/n charged covariance 恰好等于 R₁(p−ν) 乘 p-body covariance。不是忽略邻扇区的近似。

## 导数 mark 的精确变换

D B=Q_fB−Σc†Bc。对任意中性 B 和 f(Q)：

D[f(Q)B]_m=f(m−1)(DB)_m+m[f(m)−f(m−1)]B_m。

若 k=2或3 个因子的 product B 中恰有一个 p−1-body，其余是 p-body，η=(m−p)/(m−p+1)，则精确

[(kp−D)(gB)]_m=η[(kp−1−D)B]_m+[(k−1)p/(m−p+1)]B_m。

第一项是自然低体阶 product 的 contraction-count mark；第二项 O(p/N)=O(1/p)，1−η=O(1/N)。T₃/R 用 total三槽 mark、12/23双槽 mark 的有限线性组合写出，故包括13非相邻 mark。若低体阶边落在 vertex 外，mark不变；落在内由上式处理。

固定 diagram 的集合大小为 p 或 p−1。已有 BNR generic-length 图规则和本项目局部 determinant mark 中，二重交叠均值 r_i r_j/N→λ，三重交叠 O(p³/N²)→0。无放回抽样指数矩控制复 mark 邻域，一致收敛允许固定 mark 导数过极限。每个带mark图在这次单体阶移动下具有同一极限，且有界，所以 additive O(1/p) 确实消失。中央 charge 窗 N^(2/3) 一致，窗外指数尾压过多项式系数。

## 任意固定两弧矩到谱测度

对 B=Y_N或R_N、固定u,v，误差先在算符层用

|〈(y−yhat)H^u B H^v〉|≤||y−yhat||₂ sqrt〈B H^(2u) B H^(2v)〉→0。

第二因子由已有完整对角核控制。中心化 Y/W₃R 先精确删除块内 Wick 配对，再作 signed edge 替换。每张剩余图仅含一个 yhat，故恰是一条 p−1 lifted edge 减去 p edge；两项有相同极限，逐图相消。有限张图相加给所有固定 u,v 的 FY/FR 矩为0。

有限谱矩阵正性给 spectral Cauchy–Schwarz，三个对角 Mehler 测度给紧性与高阶一致可积。任意 cross 子列极限支撑在同一个紧能量方形，所有多项式矩为0，故测度为0。这才得到完整双能量 cross，而不只是时间对称性或低阶拟合。

有限系统的〈F J_HY_N〉=〈FG_N〉通常非零；p〈FJ_HY_N〉/〈F²〉→(1−2λ−c)/2，c=2λ(1−q)/(1+q)。数值保留这条1/p泄漏，用来确认我们没有把有限系统误当精确解耦。

## 扣除明确守恒多项式后的完整三阶谱

沿用上一轮 a,b,c_N、bar t、bar z：P_N=(a c_N+b)H³+[bar t−(a c_N+b)bar z]H，U_N=M₃−P_N。

精确 U_N=R+aJ_HY_N+F+(3/2)LY_N−aB_S−bB_Z；最后两个线性项 L²消失。

设 μ_r 是单位归一 q-Mehler 双能量测度，Σ=E+E′、Δ=E−E′，则

dν_(U,U†)=V_F dμ_q + V_R dμ_(q³) + W[aΣ/2+(3/2)Δ]² dμ_(q²)。

V_F=64J_s⁶q²，W=4J_s⁴[1/λ+4q/(1+q)]。

V_R=σ⁶(λA₁+λ²A₂)/D²，D=(1+q²)(1+q+q²)，
A₁=4(1−q²)(1+q²)(1+q+q²)²，
A₂=16q²(1+q)(1+q²)(1+q+q²)(2q²−q+2)。

这是正测度，M₃非厄米故必须用U†。归一范数

N_U=V_F+V_R+Wσ²[a²(1+q²)/2+(9/2)(1−q²)]。

ΣΔ项在固定频率积分中由(E,E′)→(−E′,−E)相消，但双能量测度中仍保留。固定正λ无零频atom。未扣除的原始M₃还含P_N静态部分；本轮没有声称投影掉有限N全部commutant。

## 小λ的物理形状

V_R~32J_s⁶/λ，W~4J_s⁴/λ，V_F~64J_s⁶，a~−8λ²。能量和项的范数O(λ²J_s⁶)可忽略；二阶时间后代贡献144J_s⁶/λ，因此N_U~176J_s⁶/λ。

由已有C_d→sech^(2d)(2J_st)及更高频率矩控制，规范化U两点趋

[47 sech⁶(2J_st)−36 sech⁴(2J_st)]/11。

等时R占2/11，(3/2)LY占9/11；第一正零点J_st≈.2639534223，归一μ₂=552J_s²/11。负时间相关不违背谱正性：S_U=[47(ω/J_s)²+128]S₂/880≥0，S₂=πω(ω²+16J_s²)/[96J_s⁴sinh(πω/(4J_s))]，S_U(0)=16/(165J_s)。这是对原始导数插入结构的修正，不能把纯R的sech⁶直接贴给M₃。

## 计算与保存

- CAR Casimir代替rank基底显式实现R₁：C_ad(B)=1/2Σ[Eij,[Eji,B]]，本征值ν(N+1−ν)。三个有限系统8个charged blocks、288772 entrywise等式；15个双/三product marked identities。整数交叉相乘并限制绝对值<2^51，无蒙特卡洛。
- 12组全rank，λ=.1/.25/.5，p16/32/64/128，80/120位L²误差及finite泄漏；不是区间认证。
- 五λ Jacobi cutoff35/55/142/309/582，+16核对。给包含J²初态的Taylor截断界<1e−12；独立128/192节点Mehler在λ.25/.1核时间、范数、μ₂，最高时间差3.57e−14。λ.005仍与极限有最大绝对差.01307；过零处不使用相对误差假装高精度。
- 九个独立弦矩多项式、精确Laurent系数176/8832、四个频率的50/80位Fourier核对。保留两time kernels相消条件数。
- 初始PDF有1.11525pt溢出，缩短一处句子后七页重编通过，初始log保留。新谱checker首断言误把expand作用于2+q而非r(2+q)，结构比较失败；保留旧脚本，改为展开差精确比0后通过，公式和容差未改。

下一步不再重复单个两点谱：先算真正混合Y/R四点与F混合高点，检查不同singlet的互交叉权和非配对项。以X R X R、X X R R及允许的三点为入口，先文献核查再扩marked计算。随后与已有体平移/质量字典结合；完整de Sitter全息和有限G动力学仍开放。
