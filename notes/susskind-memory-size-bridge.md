# 秩标记、局域费米子大小与固定次数的缺额分布

2026-11-07，本地研究阶段。正式推导见 [十页 LaTeX 稿](../manuscript/susskind-memory-size-bridge/main.tex) 与 [PDF](../output/pdf/de-sitter-susskind-memory-size-bridge.pdf)。这是 neutral-marked-series 的下一步：给其中的最高 Gaussian chaos 标记场提供可计算的共同耦合矩，不改写已经完成的时间级数尾界。

## 模型和对象

标准 CAR、允许集合重叠的 complex Gaussian p-body Hamiltonian；所有电荷扇区使用同一原始耦合 J。态为 β=μ=0 的归一化 annealed full-Fock trace。偶 N，p²/N=λ_N→λ>0，固定 J_s，ℏ=1，string 时间。所有次数 j,k,a,b 先固定，再取 double scaling。

记 C_N=E tr H²→C=J_s²/λ，q=exp(−4λ)，Z_j=W_j(H^j)，U_j=(jp−R)Z_j。W_j 是随机耦合的普通 Gaussian chaos 投影；R 是 U(N) 表示的秩标签，不是矩阵的秩，也不是 q-Hermite 次数。

## 1. 有限代数桥与实际误差

在中性 full-Fock 算符上定义

\[
DX=QX-\sum_i c_i^\dagger Xc_i,\qquad
\bar DX=(N-Q)X-\sum_i c_iXc_i^\dagger.
\]

先复合映射，再限制电荷 m。Casimir 恒等式为

\[
\mathcal C_x=(N-m+1)D+m\bar D-D\bar D=R(N+1-R).
\]

L=(D+\bar D)/2 是 Majorana Clifford length 的一半：LΓ_A=|A|Γ_A/2。设 S_j=(jp−L)Z_j，d=jp，A_d=d−D，B_d=d−\bar D，则精确有

\[
(N+1-2d)(U_j-S_j)=A_dB_dZ_j-(d-R)^2Z_j
-\frac{N-2m+1}{2}(B_d-A_d)Z_j.
\]

乘法顺序固定为 A_dB_d。清分母后的恒等式覆盖有限共振 N+1−2d=0。原始词中的 singleton 站点被局部缺额标记消去；一次标记不使其他 singleton 产生新活动。固定次数的 overlap 指数矩和旧秩四矩给中央紧窗 |m−N/2|≤Xp 上 ||U_j−S_j||_{2,m}=O_{j,X}(p^−1)，精确半填充改善到 O_j(N^−1)。固定次数的 H 左右乘法用 Gaussian degree 乘法界控制，因此可在两场、固定能量插入矩中使用此桥。

## 2. 非 Abelian 大小代数与全迹秩尾

Tarakemeh–Xu 的辅助 SU(2) 满足 C_x=N(N+2)/4−J⃗²。其大小 s=N/2−J 正是这里的 R；L=N/2−J_x。R 与 L 在 full trace 下可交换，0≤R≤L；不能把 spin 标量 J 当成 L。

令 T=(D−\bar D)/2，δ=L−R。精确 ladder 为

\[
[T,L]=T,\quad T+T^\dagger=(2Q-N)\cdot,\quad
T^\dagger T=\delta(N+1-L-R).
\]

\[
\|T^rX\|_2^2=
\left\langle (\delta)_r^\downarrow
\prod_{v=0}^{r-1}(N+1-L-R+v)\right\rangle_X.
\]

局部表 I→0、2n−I→I、c†→c†/2、c→−c/2。中性词中后两类的总贡献相消；重复作用收缩 r 个不同的原始 overlap 站点，带 r! 重数。固定 j,r 的系数由 const·(1+K)^{2r}c_j^K 支配，故 ||T^rZ_j||₂ 有统一界。对 jp≤N/4，ladder 各因子≥N/2，再用 Stirling 展开得所有固定 k 的 ⟨δ^k⟩_{Z_j}=O_{j,k}(N^−1)。特别是全迹 ||U_j−S_j||₂=O_j(N^−1/2)。Q−N/2 是不能对任意叠加态直接套固定支撑近似的反例；共同 J 的 Wick 结构在证明中实质使用。

## 3. 两标记的显式弦图规则

对两个长度 j,k 的连续 coupling-Wick blocks，中间隔普通 H^a,H^b。总长度 2n=j+k+a+b。对所有配对 π 求和，但禁止同一 Wick block 内配对。令 I_1,I_2 是接触两个 block 的弦集合，I_1∩I_2 是桥弦；定义

\[
K(I)=|I|^2-4\operatorname{cr}(\pi|_I),\qquad
K_1=K(I_1),\ K_2=K(I_2),\ K_{12}=K(I_1\cap I_2).
\]

\[
C^{-n}\lim E\operatorname{tr}(U_jH^aZ_kH^b)
=\sum_{\pi\ \mathrm{allowed}}q^{\mathrm{cr}(\pi)}\lambda K_1,
\]

\[
C^{-n}\lim E\operatorname{tr}(U_jH^aU_kH^b)
=\sum_{\pi\ \mathrm{allowed}}q^{\mathrm{cr}(\pi)}
(\lambda K_{12}+\lambda^2K_1K_2).
\]

奇总长度为零。full trace 与中央紧窗均成立。单图的 K 可负，不可去掉符号；同一站点的 λK_{12} 不能省略。

物理 Clifford semigroup 为 P_t Γ_A=e^{−t|A|}Γ_A。对 block 使用 e^{tjp}P_{t/2}，源导数是 S_j。局部 double incidence 的规则为 n→n+(e^t−1)I/2 和 1−n→1−n+(e^t−1)I/2。先保留费米子重排符号，再取局部 trace，不能用可能为零的无标记 trace 作除数。BNR 的 pair-Poisson 极限给归一化指数 λ[Σ_{u<v}(R_uv(z)−1)−n]；−n 扣除单弦方差。对两个平方为零的源，指数化为 q^{cr} exp[λ(K_1z_1+K_2z_2+K_{12}z_1z_2)]。两集合碰撞只有四个 incidence，三阶及以上 connected 局部标记为零；这并不自动证明任意多个秩标记的 L^p 替换。

## 4. 闭式矩与完整离散分布

设 P_j(q)=[j]_q!，𝒟_j=j²−4q∂_q。等时有

\[
\lim\langle Z_j,U_j\rangle=C^j\lambda\mathcal D_jP_j,
\qquad
\lim\|U_j\|_2^2=C^j(\lambda\mathcal D_j+\lambda^2\mathcal D_j^2)P_j.
\]

j=1,2 的归一范数分别 λ+λ²、4λ+16λ²；j=3 为 λ(9+10q+2q²−3q³)+λ²(81+50q+2q²+9q³)。

实际 size 源函数在固定 j、紧实源窗收敛到

\[
G_j(t)=e^{\lambda j^2(e^t-1)}[j]_{e^{-4\lambda e^t}}!.
\]

它尚未归一化。其整数权重为

\[
w_{j,d}=e^{-\lambda j^2}\frac{\lambda^d}{d!}
\sum_{\sigma\in S_j}(j^2-4\operatorname{inv}\sigma)^d,
\quad \sum_dw_{j,d}=P_j(q).
\]

概率是 w_{j,d}/P_j(q)，约定 0^0=1。反转置换把底数 a 配成 2j−a，所有整数次幂的配对和非负。R、L 的整数可交换谱给两个加权分布的 sum-norm TV 距离≤2C_N^{−j}||(L−R)Z_j||₂²。因此实际 rank 缺额也在 TV 下收敛到该分布；所有固定多项式矩用 general-r ladder 收敛。这里没有断言 finite-N rank 正指数矩收敛。

j=1 为 Poisson(λ)；j=2 为 q/(1+q) 的零点质量，加 1/(1+q) 的 Poisson(4λ)。一般 j：E d=λ(j²−4E_q inv)，Var d=E d+16λ²Var_q inv。图用 λ=1/4、j=1,2,4,6 的解析极限概率；显示 d≤28，最大显示尾约1.3753e−8，保存的更大截断有≤1e−42 的精确有理尾界，非有限 N 模拟。

## 5. 两能量测度与接下来的静态问题

固定有限 U_j 列表的两能量矩确定唯一正矩阵测度。|K_1|≤j²、|K_2|≤k²、|K_{12}|≤min(j,k)² 给普通 q-Gaussian 矩支配；极限每个能量坐标支撑于 [−R_H,R_H]，R_H=2√(C/(1−q))。正性、高固定矩的一致可积性和谱 Cauchy–Schwarz 给有界连续测试函数收敛，包括固定实时间的两点相位。

U_j 的确定性能量回归为 a_j^reg H_j^{(q,C)}(H)，a_j^reg=λ𝒟_jP_j/P_j；H_j 是 monic spectral q-Hermite。有限扇区 Q_{j,m} 的中央系数收敛与旧全迹尾界提供转换。回归本身不能排除额外随机对角谱质量。

j=1 可完整分解：ν_11=λ²E²ρ_q(dE)δ_E(dE′)+λCμ_q。扣 λH 后仅留 λCμ_q，无静态原子；F=8λCq(U_1−λH)+o(L²) 恢复旧 ν_F=64λ³C³q²μ_q。

## 实际验证和继续位置

- finite-bridge.json：96 整数词/扇区、384 恒等式含6共振；N3p1 九坐标 j1–4 的全部1428最高chaos扇区系数、5712恒等式；168偶 Clifford strings。96例有限桥误差均非零。
- marked-chord-checks.json：160有限CAR词、640比较，覆盖两个标记所有子集；54矩多项式。
- geometry-checks.json：34布局、5716合法配对、25884精确比较；1085高阶connected项为零；保留298负单标记和34负双标记系数。153 all-bridge配对独立核等时式；16个U1能量例独立核q-Hermite/Mehler。
- rank-law-checks.json：96 neutral matrix units及精确联合投影；1542整数正性检查；18概率例、54生成函数例；截断19–227，精确有理尾≤1e−42。80/120位最大scaled差约3.6833e−78。

四个科学输出均首次运行通过，未重跑或覆盖历史数据。版本、输入输出SHA、最终图/PDF/日志/十页渲染及前驱科学文件在 local-evidence.json 绑定。初次编译的 mathscr 缺宏失败、核版时发现的求和下标排字问题均保留诊断副本，最终双遍编译与十页目视检查通过。

下一步对 forced recurrence 的有序 Gaussian-Wick 后代 eta-hat 应用局部规则，保留 composite blocks 的禁配对；由已有阶乘尾组装全时间 neutral kernel，单独识别零频质量后再传递 exact pinching。随后把 charged defect 和 reverse source 放进同一 J 图，计算 AA/AB/BA/BB。完整有限时间联合核、体质量/位置字典、finite-G 四点动力学仍开放。

Issue #1 的已有独立算术核查保留；9890旧端点重现失败未删，21143单元重积分未做。本阶段未重新查看新评论、未回复issue、未commit/push/ZIP。
