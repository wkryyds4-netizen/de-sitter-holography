# 原始三阶 singlet 的归约与交叉谱

2026-11-06，本地研究阶段。正式推导见 [LaTeX](../manuscript/susskind-cubic-cross/main.tex) 与 [七页 PDF](../output/pdf/de-sitter-susskind-cubic-cross.pdf)。本轮无 commit、push 或 ZIP。

## 这一步解决什么

沿用 Susskind 静态补丁 → complex DSSYK∞ → 导数 singlet 主线。上一轮已求 F 的完整 q 核，但原始 M₃ 还含 R、S、Z 及时间后代，不能把各个对角谱直接相加。本轮求出 S/Z 的精确归约、主要有限 p 余项，并消去 R–Y 的完整交叉谱。尚未求完 F–R 和能量分辨的 F–Y，因此 raw M₃ 的最终谱仍有具体待算项。

标准 CAR，包含重叠指标的 complex Gaussian p-body H=H_s，退火归一 trace、β=μchem=0，ℏ=1。L=[H,·]，J_H B={H,B}/2。先 N,p→∞、λ=p²/N>0 固定，q=e^(−4λ)、σ²=J_s²/λ；随后才考虑小 λ。

定义 A=Σ(Lc_i)†Lc_i，M₃=Σc_i†L³c_i，T₃=(M₃+M₃†)/2；T/S/Z=W₃(T₃,J_H A,H³)，R=T−aS−bZ，a=−(1−q)²/(1+q²)，b=2λ(1−q)(1+q³)/[(1+q²)(1+q+q²)]。W_j 为普通 Gaussian coupling chaos。

## 精确一次通道

在 charge m，H=ΣνHν，〈Hν²〉=Cων：

- C=w binom(m,p)binom(N−m+p,p)，C_−=C_(m−1)。w 是协方差，不是 λ。
- κν=w binom(m−ν,p)binom(N−m+p−ν,p)。
- W₁H³=ΣzνHν，zν=2C+κν。
- W₁J_H A=ΣsνHν，sν=CK₀+κνJν。
- K₀=2p−m+mC_−/C+p²/m。
- Jν=p[p(N−m+1)−ν(N−ν+1)]/[(m−ν)(N−m−ν+1)]。

边界可去奇点使用 sν=[(3p−m)zν(m)+(m−p)zν(m−1)−tν(m)]/3；tν 是前轮已有 W₁T₃ 系数。它来自 T₃=(3p−D)H³−3J_H A 和 D[zν(Q)Hν]_m=[mzν(m)−(m−p)zν(m−1)]Hν,m。

令 bar z=Σωz、bar s=Σωs，B_Z=Σ(z−bar z)Hν，B_S=Σ(s−bar s)Hν，F=Σ(t−bar t)Hν。于是 Z=H³−bar z(Q)H−B_Z，有限 N 下通常并不守恒。

## 双标度余项与尾控制

得到 ensemble HS L² 极限：pB_Z−F/2→0，pB_S−(1−2λ)F/2→0。

证明用 s=p−ν→Pois(λ)，κν/C 的精确乘积与 Jν 有理式。固定 s：

- p[(κ_(p−s)−κ_p)/C]→4λq s。
- J_p→−2λ，p[J_(p−s)−J_p]→4λs。
- p(s_(p−s)−s_p)/C→4λq(1−2λ)s。
- 旧 t 通道差 (t_(p−s)−t_p)/C→8λq s。

在有限权重下先中心化，再过极限。中央电荷窗 |m−N/2|≤N^(2/3) 内分母≥cN，κ/C≤1，rank 导数满足 |∂ν(κ/C)|、|∂νJ|≤Kλ/p，|J|≤Kλ，故 p 倍差被 Kλs 控制。旧 t 的 Lipschitz 界与 factorial rank tail 给完整 L² 支配；窗外用精确无奇点式的多项式界及指数电荷尾。不是只对有限 deficit 点态展开就直接交换求和。

## F–Z 的有限 p 交叉项

由 Z 的守恒多项式分解与〈FZ〉=0，所有有限系统/实 t 均有

〈F(t)Z〉=〈FB_Z〉−〈F(t)B_Z〉。

所以 unscaled cross 在 DS 中均匀趋零；固定时间的首项为

2p〈F(t)Z〉→V[1−C_F(t)]，V=64J_s⁶q²，C_F 是上一轮 r=q Mehler 核。

动态矩单独用 ||LB||₂²≤12C||B||₂² 控制线性 Gaussian 余项，而非对点态时间极限求导。得到 2p〈FL²Z〉→−128J_s⁸q²(1−q)/λ。

有限例 N5,p2,m3,w1：三两弧矩依次为 T:(4590,5400,4590)，S:(−6570,−5760,−6570)，Z:(−4680,−3600,−4680)，交叉曲率为 −1620、−1620、−2160。等时 coupling-chaos 正交不能替代动态正交。

## S 与原始 M₃ 的精确归约

a₀=EA=(2p−m)C+mC_−，O=A−a₀I，Y_N=O−c_N(H²−CI)，其中 c_N=〈O(H²−C)〉/〈(H²−C)²〉。精确

c_N=(bar s−a₀)/(C+bar κ)→c=2λ(1−q)/(1+q)。

G_N=W₁J_HY_N=B_S−c_NB_Z，S=c_NZ+J_HY_N−G_N，且 pG_N→(1−2λ−c)F/2。

令 P_N=(a c_N+b)H³+[bar t−(a c_N+b)bar z]H，它逐电荷守恒。则

M₃−P_N=R+aJ_HY_N+F+(3/2)LY_N−aB_S−bB_Z。

这是有限 N 精确恒等式，L A=L Y_N。DS 中 Y=√W_d X，W_d=4J_s⁴[1/λ+4q/(1+q)]。能量矩阵元中的 Y 系数为 a(E+E′)/2+(3/2)(E−E′)。M₃ 非厄米，组合谱须正确放置 dagger。

## R–S/R–Z 全阶相消与 R–Y 谱

复用 cubic-kernel 已穷尽的连续三槽端点模板，仅更换双线性 mark 权重。R 系数 P/D 与 −λB/D，S 对应 mark12/23 的各 1/2，Z 无 mark。177 模板、6234 配对中，RS 的 λ、λ² 系数及 RZ 的 λ 系数全部作为 q 多项式恒等于 0。37404 个局部 minor 与已有 determinant 实现逐项一致。

解析覆盖任意阶的原因：外部 spectator 弦在三槽内部置换下的 crossing 数不变，可提出；余图由 k=0..3 桥数、外端点 A/B 颜色及两弧 split 完全分类。这不是低阶样本外推；但本轮确实复用旧模板生成器，不声称做了第二套独立枚举。

因而任意固定 u,v：〈RH^u S H^v〉DS=〈RH^u Z H^v〉DS=0。将 S=cZ+J_HY−G 代入，G 的固定多项式矩贡献由线性 Gaussian 矩界及 Hölder 消失，故 (E+E′)dν_RY=0。

R/Y 对角双能量测度为 r=q³/q² 的绝对连续 Mehler 测度。谱 Cauchy–Schwarz 保证交叉测度也绝对连续，不能支撑在 E+E′=0 的零面积直线上。因此 dν_RY=0。对任意子列极限均成立；紧性、一致可积与紧支撑矩确定性完成极限论证。

## 实际验证与下一步

五份 JSON 保存于 results/susskind-cubic-cross：三 finite cross 系统、32 个 S 一次通道；15 组全 rank 80/120 位余项比较；177 个精确模板和37404 minor；另以 rank 双线性式独立复现 Z cross 曲率，三个确定性矩阵共12个有理恒等式验证。λ=.1,p256 的 Z/S/G scaled L²平方相对 F norm 分别6.537e−6、8.038e−6、7.476e−6；λ=.5,p256 的 S 为1.936e−5。精度比较不是区间认证。

新增 checker 初次对未展开 SymPy 表达式使用结构相等而停止；诊断 residual 精确展开为零。保留初始脚本及失败记录，改为 exact expand 后全通过，没有放宽容差或改物理公式。原有限 T/S 六配对结果尚无第二种完整独立枚举，只对 Z 曲率提供 rank 公式交叉核对。

下一步集中求 ν_FR 和能量分辨 ν_FY；普通时间 F–Y 的对称性零不消掉 F–J_HY。优先把 rank filter 的协方差桥接接到 marked chord，保留同 J 正常序相消。完成 raw M₃ 谱后再做混合高点和体空间/质量字典。
