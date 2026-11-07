# 相邻电荷打开普通乘积禁戒通道

2026-11-06，本地。[四页 PDF](../output/pdf/de-sitter-susskind-memory-charge-lift.pdf)、[LaTeX](../manuscript/susskind-memory-charge-lift/main.tex)、[来源](../literature/susskind-memory-charge-lift-sources.md)。不 commit/push/ZIP。

本轮解决 NEXT 的第一个分岔：上一轮 N5p2m3 四个 map 在每块都共线，不能推广为一般的 H² scalar rank filter。半填充的普通乘法有精确奇偶选择规则；charge removal 与实际记忆收缩能进入被这个乘法消掉的通道。后续改用两个电荷上的相干 reduced amplitudes。

## 选择规则与一般反例

标准 CAR，固定 charge m、归一 HS 迹，完整重叠 complex Gaussian p-body ensemble。N=2m 时定义线性反自同构 τ(X)=star X^T star^{-1}，star 为有符号 complement 映射。τ(Eij)=δij−Eij，因此与 adjoint U(N) 作用可交换；每个中性 rank 仅出现一次，在最高权 traceless monomial 上求符号给 τP_r=(−1)^rP_r。

所以 P_r{X_a,Y_b}=0 当 a+b−r 为奇数。尤其同 rank 平方只有偶 output ranks。这里不意味着 tensor product 不含该 irrep，也不排除另一个 equivariant map 进入它。

取任意 1≤p<m=N/2：

\[
H=\prod_{a=1}^p(n_{2a-1}-n_{2a}),\quad
t_a=n_{2a-1}+n_{2a},\quad b_a=(n_{2a-1}-n_{2a})^2.
\]

每次 partial contraction 都由被遗漏 pair 的正负项相消，因此输入为纯 rank p。H²=∏b_a 为偶，而直接移除占据模式给

\[
A=\sum_a t_a\prod_{b\ne a}b_b,\qquad
A_{odd}=\sum_a(t_a-1)\prod_{b\ne a}b_b\ne0.
\]

恰一个 pair 双占/空占、其余选中 pair 单占时为 ±1，其他配置为0。于是

\[
\|H^2\|_m^2=\frac{2^p\binom{N-2p}{m-p}}{\binom Nm},\quad
\|A_{odd}\|_m^2=\frac{p2^p\binom{N-2p}{m-p-1}}{\binom Nm},\quad
\frac{\|A_{odd}\|_m^2}{\|H^2\|_m^2}=\frac{p(m-p)}{m-p+1}.
\]

这族沿 DS 序列也存在，但它是确定性输入的范数，不能当作 Gaussian 平均强度或新动态通道。N4..12 半填充全部 p<m 的15系统穷尽配置核 removal、Casimir eigenvalue 与两个范数式。

## 真正 Q_Ecal 的精确见证

沿用 core Ecal=Σ(ad_H²c)†(ad_H²c)+[H,[H,A]]/4，Q=ΔEcal/2，W2Ecal=Q−2EEcal。取 N6,p2,m3,w1，h=(n1−n2)(n3−n4)|3、k 同一多项式在 charge2。输入 rank2（Casimir20）。定义 D=n1+n2+n3+n4−2I，则

\[
P_1h^2=0,\quad P_1A=D,\quad P_1Q_{Ecal}=99D,
\quad \|D\|_3^2=2/5,\quad \|P_1Q\|_3^2=19602/5.
\]

P1 不受 scalar mean 扣除影响，因此同时是 W2 的反例。纯输入 rank 上 F∝h，P1JhF=0，也不能补回此通道。独立对所有225实Gaussian坐标做原始 commutator core 的四次 polarization，精确重现 explicit covariance Q；raw word identity 也相等。全部4个 full matrix Casimir projectors 核 τ 的符号，不只看 diagonal 子空间。

不改变旧 N5 共线性/Gram 结果。这里的 diagonal realization 与 D 对易，所以不把它当动态见证；也尚未求本系统完整 block Gram/multiplicity。

## 可用的一般接口与正界

分别保留 Z⁺_ab;r=P_r(h_a⊙h_b)、Z⁻_ab;r=P_r Rem(k_a⊙k_b)，a<b 时 ⊙ 为反对易子，a=b 时为平方。

精确 A_ab;r=Z⁻+(2p−m)Z⁺，不能除以可能严格为0的 Z⁺。

primitive rank tensor 的跨电荷归一取 h_a=e_a(m)T_a^(m)，e_a(m)=binom(m−a,p−a)。令 K^(m)_ab;r,ρ 为统一 intertwiner 基中的 sym product coefficient，则

\[
\mathsf A_{ab;r,\rho}=(m-r)e_a(m-1)e_b(m-1)K^{(m-1)}_{ab;r,\rho}
+(2p-m)e_a(m)e_b(m)K^{(m)}_{ab;r,\rho}.
\]

这是前驱 Rem(T_r^-)= (m−r)T_r^+ 在 reduced-product 上的应用。ρ、相位及两个电荷的 cross 必须保留。Kota–Vyas(34)/(35)提供现有 recoupling 框架，本轮未解出所有 K。

对共同 rank 范围内任意 X_r：

\[
\|Rem X_r\|_m^2=\ell_r^2\|X_r\|_{m-1}^2,\qquad
\ell_r^2=\frac{m(m-r)(N-m+1-r)}{N-m+1}\le m^2.
\]

由 disjoint primitive monomial 范数 binom(N−2r,m−r)/binom(N,m) 在相邻电荷的比直接得出；equivariance/Schur 推广到整个 rank 空间。共同范围外由 intertwiner 消掉。8个直接 CAR 检查全部相等。

故对 a<b，EE||Z⁻||²≤4ell_r² Ca⁻Cb⁻ π_ab^(m−1)(r)，a=b 的 W2 版本系数为2(Ca⁻)²。它将上轮正 π 接到下电荷乘法，但没有求 EE〈Z⁺,Z⁻〉，不能忽略相干抵消。对整个 Q 尚需含 Ψ 插入的更多幅；不宣称两个 Z 已张成 Q。

## 保存与继续

final-checks.json 为最终验证。checks.json 是添加解析计数核对前的第一次通过记录，preliminary-checker.py 的 SHA 与该记录原始代码 SHA 相同；未覆盖。一次临时 python -c 从根目录导入失败，改在 code 目录执行后通过，无科学计算失败或容差调整。

四页双遍编译与95dpi逐页视觉检查通过，无未定义引用或公式溢出。下一步从含 multiplicity 的相邻电荷 coherent amplitudes 推 Q_Ecal，不再假设 scalar filter；之后才计算 DS 范数/cross、raw Z2 与有限 memory，并接体空间/质量字典。
