# 实际有限时间记忆场的统一二阶界

2026-11-07，本地阶段稿。正式推导见 `manuscript/susskind-memory-uniform-bound/main.tex`，九页 PDF 为 `output/pdf/de-sitter-susskind-memory-uniform-bound.pdf`。本轮推进的是实际有限时间场，不再只限于它的有限阶 Taylor 系数。

## 对象与结论

标准 CAR、同一原始 J 的重叠 complex p-body Gaussian 模型，β=μ=0 全 Fock 归一迹，ℏ=1/string 实时间；偶 N，λ_N=p²/N 在固定正紧区间，J_s 与时间窗 S 固定。所有范数默认 annealed normalized L²；带下标 r 的辅助范数明确是 Gaussian number Sobolev 范数，非 Schatten-r。

沿用上一阶段的 T、F=T exp(isH_-)，以及 exact dilation 的 R、K、D=QKQ。能量投影 P_h 保留全部简并块。定义

\[
X_{0,m}(s)=(1-P_h)m e^{-ihs/2}R^\dagger e^{isK}Re^{-ihs/2},\qquad
\mathsf F_m(s)=(1-P_h)m e^{-ihs/2}R^\dagger e^{isD}Re^{-ihs/2}.
\]

本轮证明全 charge 的 `sup_N sup_|s|≤S ||X0(s)||₂ < ∞`。结合 return-background 阶段已有 `d_N(S)=sup||Fscr-X0||₂→0`，实际场也统一有界，且

\[
\sup_{|s_1|,|s_2|\le S,\ t\in\mathbb R}
|\mathcal C_{\mathsf F}(s_1,s_2;t)-\mathcal C_{X_0}(s_1,s_2;t)|\to0.
\]

这是两个有限 N 两点函数的渐近等价；其共同极限的存在和值仍待求。二阶矩界本身也不是 path-space tightness 定理。旧 Return 定理的 λ_N 收敛假设由正紧区间的子列论证衔接。

## 推导链条

1. **加权 charged transfer。** 低 chaos polar 缺陷加上高次数吸收，给 G=U H_-−H U 从 Sobolev r+5 到 r 的 `O(√B_m/N)` 界，`B_m=λ_N+λ_N²+λ_N(m-n)²/N`，`n=N-m+1`。配合 number-weight 传播与 `[p Ncal−Rcal,H]` 界，直接对初值 `h_-^j` 用 Duhamel；得到任意固定时间导数、两侧 h 幂和 Sobolev 权重下的 charged OU 误差同阶。没有对已有无权误差求导。
2. **精确 OU force。** `Λ=-F''+iα{h,F'}+α²hFh`，所以 `||Λ−Λ_OU||₂≤K_S√B_m/N`。乘 m 后仍为有界项，不能扔掉。`Λ_OU=ε E_u[u exp(is(αh+√εu))u]`，`ε=1−α²`；外侧和指数内用同一个 u。
3. **只在独立 Gaussian 副本强度上展开。** 保持 `U_t=exp(iαht)` 精确，偶 Dyson 展开为 `Λ_OU=ε Φ(U_s)−ε² I_s(h)+E_s`。四 u 项保留全部三种 Wick 配对，`||I_s||∞≤3 C²s²/2`；余项 `||E_s||₂≤ε³|s|⁴√10395 C³/24`。用单位单纯形的 s² 参数化兼容正负 s。中央区 `mε≤2p`，`mε²≤16λ_N`，`mε³≤128λ_N²/p`。因此最后仅中性项的 p 放大需要新控制。
4. **coupling chaos 到逐样本能量函数的定量桥。** 令 `Z_j=W_j(H^j)`，`a_j=Σ_{v<j} μ(vp)`。Gaussian annihilation 精确给 `A Z_j=C a_j Z_{j-1}+r_j`，余量界只带一个 `(j−1)!`。此前正 rank 链的全次数二矩控制 r_j。保留有限 μ(vp)，用 `μ(jp)≤μ(p)^j` 得 `1≤a_j≤A=(1−exp(−3λ_-))⁻¹`。
5. **辅助 Jacobi 递推与全复平面展开。** `Q_{j+1}=xQ_j−Ca_jQ_{j−1}`、`E_j=Z_j−Q_j(H)`，E_j 的 coupling degree ≤j−2。递推给 `||E_j||₂≤(p/N)L K_*^j√(j!) C^(j/2)`，`K_*=4+A`。有界辅助 Jacobi 算符给展开系数的 `1/j!` 界。逐项 Jacobi 方程证明 `Σ c_j(s)Q_j(z)=exp(isz)` 对所有复 z 成立；另用上述范数界保证实际 Gaussian H 的 annealed L² 绝对收敛。辅助谱测度从未当作实际有限 N DOS，也不需给实际 H 截谱。
6. **精确守恒参考。** 对确定性 `|θ_j|≤1`，`g_{θ,s}(H)=Σ θ_j c_j(s)Q_j(H)` 逐样本与 H 对易，且 `||θ(Ncal)e^{isH}−g||₂≤2(p/N)R_S`。取 θ_j=μ(jp)，结合已证明 rank–chaos 缺额，得到

\[
p\|(1-P_h)\Phi(e^{isH})\|_2\le C\lambda_N(12D_S+2R_S).
\]

θ 不能依赖原样本 J。我们没有宣称 Gaussian chaos 乘子本身守恒。

7. **全 charge 与两点函数。** 中央区四项上界依次为 `K√B_m`、`2Mλ_N(12D_S+2R_S)`、`24λ_N M²S²`、`O(1/p)`。`Σ b_m(m-n)²=N+1` 控制中央求和；外部 charge 用固定导数的粗界与 `N² exp(−N/8)`。最后 Cauchy–Schwarz 将旧 d_N→0 提升为协方差替换，对所有外部实时间 t 一致。

## 核查

- `jacobi-recurrence-checks.json`：三个完整 CAR 系统 N3p1m2、N4p1m2、N4p2m3，最高 degree 5、4、3，共12个完整 Hermite 案例；独立 covariance-word 重构 annihilation，并核单 factorial 的正 tensor 界与低 chaos 多项式项。
- `dyson-force-checks.json`：Pauli/2 四坐标 Hermitian 2×2 Gaussian，C=1，h=.4σ_x+.7σ_z；3个完整四点 word，9组 ε=.04/.16/.36、s=−.7/.7/1.4；Gaussian quadrature 10/14 与 time quadrature 12/18 交叉。最大 Gaussian 阶差约4.61e−15。这是非对易恒等式测试，不是 DS 模拟，也非区间积分认证。
- `jacobi-entire-checks.json`：105个精确有理几何不等式，3套参数在辅助谱之外 x=±6 的6项检查；65多项式、时间级数150阶、80/120位对照。误差小于1.1e−56，精度差小于5.2e−80。统一结论来自解析证明。
- 三个科学 driver 首次执行均 exit 0；不重跑覆盖。最终九页双遍 pdflatex，无 undefined/overfull/fatal 警告，95 dpi 九页逐页核版。首次视觉检查发现全 charge 式的转义文字及积分微分排版，修正前源/log/第8页保留在 diagnostics；另补 λ 紧性子列衔接。

## 继续位置

下一步求 charged-transfer 缺陷、中性 covariance-filter 残量、四 u 项的共同 J 联合极限及交叉项。统一界允许有界贡献，不证明每个分量都非零或都有极限。可先从精确三分量定义及零/一次时间校准进入，再处理整个时间级数。已有 q-Gaussian Wick 条件期望可提供领先守恒背景，但不能代替有限 N collective-scale 残量。

issue #1 初核独立保留，9890旧端点回归失败未改变；没有全量重积分、外部回复或 Git 发布。本轮所有材料只在本地。
