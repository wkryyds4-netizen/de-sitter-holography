# 带电记忆场：精确转移、完整局部源与电荷背景

2026-11-07。本轮接上已封存的完整 neutral B 核，处理实际 charged A 的全部 Gaussian chaos。正式稿在 `manuscript/susskind-memory-oriented-transfer/main.tex`，阅读版在 `output/pdf/de-sitter-susskind-memory-oriented-transfer.pdf`。详细证明见本轮新建的 [proof draft](susskind-memory-oriented-transfer-proof-draft.md) 和 [charge background](susskind-memory-charge-background.md)。

标准 CAR、overlapping complex p-body、跨电荷同一原 J、β=μ=0 归一 annealed full trace，偶 N、λ_N=p²/N→λ>0，固定 J_s、ℏ=1/string 实时间；固定次数及固定时间窗先于双标度极限。C=J_s²/λ，q=e^(−4λ)。以下是模型内推导，不是 de Sitter 体字典或有限 G 引力动力学。

## 本轮完成

精确 full-Fock loss map 为 D X=QX−Σc_i†Xc_i，因此相邻 charge 传递严格写为 F_m=e^(ish)−[D e^(isH)]_m/m。令 α_m=(1−p/m)_+，则

\[
E_m=m(F_m-\alpha_m^{\mathcal N}e^{ish})
=[m(1-\alpha_m^{\mathcal N})-D]e^{isH}|_m.
\]

这里 Ncal 是 coupling Gaussian number。全 k Taylor 积分余界配旧 Gaussian Sobolev 传播，证明任意固定导数/左右能量/Sobolev 权重下，紧实时间窗内

\[
E-G_N=O_{L^2}(1/p),\qquad
G_N=[p\mathcal N-D-\lambda_N\mathcal N(\mathcal N-1)]e^{isH}.
\]

不能提前 λ_N→λ 后仍声称相同速率。m=0 与 m<p 为零，m=p 保留 OU 的 0^0=1。

完整普通幂系数 G_n 包含所有内部 Gaussian 配对。Polar defect telescoping 给 ||G_n||₂≤K(n+1)³ L^n√n!；force 系数 K_n=G_(n+2)−{H,G_(n+1)}+HG_nH 的界为 K′(n+3)^4L^(n+2)√n!。因此所有固定实时间窗的 Taylor 尾可和，能接精确 force/pinching 接口。未假设 body loss 等于 rank，也未删掉 charge-dependent Jacobi 系数的导数。

## 完整图规则

整个 field 的内部弦数 i、只有一个端点在 field 的弦集 O、d=|O|，决定其 coupling chaos。对 nilpotent z²=0，使用 exp[z(pd−D)]。对称局部 D 把 c/c† 乘 e^(−z/2)、n 变 e^(−z)n、1−n 变 1−e^(−z)n。open singleton 与 exp(pdz) 精确相消，internal chord 的两个 singleton 因子乘积亦为1。每个 pair collision 必须除以两个**带源 singleton**的乘积。

闭式为

\[
\kappa_a=d_a(d_a-1)-4\operatorname{cr}(O_a)-2i_a,
\quad
\kappa_{12}^{DD}=4i_1i_2+2b(b-1)-4\operatorname{cr}(B),
\]

其中 B=O₁∩O₂ 是直接桥弦集，b=|B|。G 的 one-source 只再减 d(d−1)，故变成 −4cr(O)−2i。两 G 场的图权重是 C^ℓ q^cr[λκ12+λ²κ′1κ′2]。mixed D/L 要另取 I_L 为neutral marked subblock的incident集，B_DL=O_D∩I_L，connected coefficient为 |B_DL|(|B_DL|−1)−4cr(B_DL)。I_L可能含较大neutral field的内部弦，不能换成整个field的open集。neutral 整个 Wick-source 禁配对区与其较小 mark 区不可混淆。

局部分类证明：internal own chord 贡献 −2；两个 open chords 留下原 pure-source 图；其他 single-source spectator 碰撞总和相消。双标记额外留下 internal–internal 的 4i₁i₂。这正是从 pure Wick source 扩展到普通 H^n 必须增加的部分。

## 电荷背景与旧谱校准

无 field-to-field Gaussian bridge 时仍可有 κ12=4i₁i₂。最小图是两个 H² 各自内部收缩，其源指数为 λ(−2z₁−2z₂+4z₁z₂)。独立有限差分给

\[
\mathcal W_0G_2=-DC_Q,
\qquad (DC_Q)_m=C_m\frac{p(N-2m+p+1)}{N-m+1}.
\]

对 X=(Q−N/2)/p，其全迹极限为 −2λC+4λCX，Var X=1/(4λ)。均值 −2λC、二阶矩 4C²(λ+λ²)、方差4λC²；与 pair-source 完全一致。它是 charge 函数，exact pinching 逐样本消去。

G₂总范数/C²为8λ+4λ²+16λ²q；仅扣确定性能量回归后仍为8λ+16λ²q/(1+q)。再扣静态charge方差4λ，才得到旧 Y₋ 方差 w/C²=4λ+16λ²q/(1+q)。该校准说明不能把全部无桥相关项直接当作已被平均能量回归消去；这里只校准norm，并未从norm替代全部两能量谱。

## 证据与继续位置

本轮五个科学 JSON 均首次运行后 exclusive 保存，不重跑 main。四类独立输入包括：完整 finite CAR/Hermite/OU、pure oriented mark、whole-field source、finite charge variance；第五类将closed geometry与局部算法及G₂旧方差分层比较。参数、依赖、精确计数及哈希以 `results/susskind-memory-oriented-transfer/local-evidence.json` 为准。有限 charge p128 的方差相对偏差仍可达8.65%，未隐藏；高精度复核不是区间认证。

下一步按 proof draft 第8节构造联合 X/H 回归，先核所有多项式共同 Gram，再传递 exact pinching。候选系数写成普通幂的 q-Hermite 展开 H^n=Σ_d a_nd H_d，i=(n−d)/2：

\[
g_n(X,H)=\sum_d a_{nd}
\left[-4\lambda q\partial_q\log[d]_q!-2\lambda i+4\lambda iX\right]H_d(H).
\]

目前保存为下一步构造入口；未宣布完整 AA/AB/BA、charged static subtraction、体字典或 finite-G 已完成。沿用旧非零 A₀–B₂ cross 作为必须恢复的校准。长期 goal 保持 active；本轮无 commit/push/ZIP/外部 issue 回复。
