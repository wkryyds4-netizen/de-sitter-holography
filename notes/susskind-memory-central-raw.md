# 中央电荷窗口的原始二次记忆系数

2026-11-06。主稿 `manuscript/susskind-memory-central-raw/main.tex`，阅读版 `output/pdf/de-sitter-susskind-memory-central-raw.pdf`。标准 CAR、完整 overlapping complex Gaussian p-body、hbar=1、同 J 相邻电荷、归一最大混合 sector trace。偶 N，m=N/2+xp，|x|≤L 固定，λ_N=p²/N 在正紧集，C_m→C=Js²/λ。先此极限，再小 λ。

## 核心结果

有限 N 的五个 raw quartic filters 重新展开后，第二 chaos 表只多

    Δ(W2P)=2 λ_N x ell_P C_m B，B=h²−C_m I，
    ell_(E,M,S,C,Z)=(2,2,2q,2q²,0)。

所以五行的新增 B 系数是 4λx(1,1,q,q²,0)；O、JhF、output-deficit S_N 的领先系数均不变。实际 core 的额外 B 系数为 4λx(1−q)²C_m。这里 S_N 按每个普通二次乘积块的 S=a+b−r 加权，不是 2p−r。

这不是把 O(1/p) 的 central trace 误差乘 p。计算直接在一般 charge 的精确 rational filters 上合并，再展开；κ_(tp−d)/C_m 的 1/p 项为 exp(−4λt)[4λd−4λ²t²]/p，与 x 无关。中央 A 范数、六阶输出矩和 factorial input tails 给同样的强替换；H-filter 减 p ell 后以 K(1+u+v+S) 控制。A 部分取输出窗口 p^(3/4)，总表余项只声明 oL²。

## F 接口

精确 cubic W1 channel t_a 的 rank-dependent 部分（除 C_m）为 κ_a E_a/C_m，E_a 的 rational expression 见稿式(8)。中央固定 u=p−a 的差满足 (t_(p−u)−t_p)/C_m→8λqu；全 rank 导数界给 ≤Ku。rank recurrence 的 factorial 尾统一，Poisson 均值仍 λ。有限权重中心化后

    F=8λ_N q_N C_m(D_in−bar u h)+oL²，
    ||F||²→64λ³ C³ q²=64Js⁶q²，
    JhF=4λ_N q_N C_m L_in+oL²。

线性 Gaussian 误差的第四矩≤3variance²，因此最后乘 h 有明确控制。一般 R_j 比值趋 (u)_下降j/λ^j，支持旧 signed covariance 的 total-variation 整图桥接。R1 proxy 是算符 L² 近似，不是只核对范数。

旧 upper/lower charge lift 是 general m 恒等式，不需重新发明。中央 DbV=pV−2D_in V−E_in V，E_in 系数 u(u−1+N−2m)/(m−p+1) 为 O(1/p)；d=m−p+1，lift correction 为 p/d=O(1/p)。因而 F–Y±、F–R、F–Q4 的有符号边相消可移到中央窗口。其余 endpoint certificates 用上一轮中央 trace theorem 转移。带能量乘子的替换需 factorial weighted errors 与固定阶 Gaussian/marked 高矩，不能只引用裸 L²。

## 原始 U2 和谱

令 Qcal=C0−alpha E0−beta M0−gamma S0−delta Z0，ell_Q、k_Q 等由上表同线性组合得到。旧 half counterterm 的公式按实际 sector 的 C_m、c_m、bar t、bar z、e_Q、K_A、P_m 计算，并在 B 括号中加 2λ_N x ell_Q，得到本轮 C2,m（主稿式(21)）。每一项逐样本与 h 对易。K_A²/m 的有限守恒贡献保留；p 增长项不可提前换 C 极限。

扣除此项后，原始 U2=mZ2−C2,m 的强展开及完整 (Y,L,U2) 正双能量谱，与旧 half raw-spectrum 相同。新增的是中央适用范围及 sector counterterm，而不是新发现 Mehler 核。五个正谱来源为 Q4、T、JhR、JhF 和 H2(E,E')Y；H2=bJ e²+bD ω²+Cκ。显式动态系数见新稿及原 raw-spectrum。旧有序小λ形状 (9sech⁸−2sech⁶)/7 因此继承，不重新数值扫描。

## 执行与继续

- 接续已有 `symbolic-descendants.json`，核对其 source hashes，没有重跑或覆盖。
- 新 `central-filter-checks.json`：123 个 finite rational coefficient references；36 组 λ=.1/.25/.5、p16/32/64/128、x−2/0/2 的全部 p+1 rank，80/120 位，保存65位一致（不表示实际舍入误差为零）。保存 variance、proxy L²、HF/FF signed covariance L1、pBZ/pBS 误差及 weighted tail。
- λ=.5,x2 的 variance 相对极限差从 p16 的 .77769 降到 p128 的 .07308；proxy L²误差平方从 .01803 降到 .0003023。有限尺寸差没有隐藏。
- 新 `central-assembly.json` 核七份 prior/current certificates 的输入 hash、central Q shift 和正谱矩阵 determinant；没有重跑旧 CAR 或 endpoint scans。

下一步全 charge 求和：控制固定中央窗外的**已扣守恒项字段**，与 binomial charge law 合并。须先明确端点 sector 的反项定义，处理 rank 缺失及 c_m 分母退化；不能凭固定 x 一致收敛覆盖整个和。之后 finite-memory 两插入余项，再接 Susskind 体介子的空间与质量字典。仅本地保存；无 commit/push/ZIP。
