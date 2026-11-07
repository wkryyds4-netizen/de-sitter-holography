# 记忆 core 的一般双电荷公式

2026-11-06，本地。[六页 PDF](../output/pdf/de-sitter-susskind-memory-two-charge.pdf)、[LaTeX](../manuscript/susskind-memory-two-charge/main.tex)、[来源](../literature/susskind-memory-two-charge-sources.md)。本轮无 commit/push/ZIP。

本轮把上一阶段所需的一般 Q_Ecal 振幅真正求出：复杂的 Ψ 插入可以由一个普适 transfer identity 化为两个相邻 charge 的乘积。两个乘法 map 保留各自所有 multiplicity，并非假设块内只有一个表示论方向。随后先合并大项，取得 H² 部分的固定缺额极限和 L² 近似；完整 A 权重及 output-deficit cross 仍待计算。

## 一般 charged transfer

标准 CAR、完整重叠 complex Gaussian p-body ensemble，E J_IJ J_KL=wδ_ILδ_JK，ℏ1。固定 charge maximally mixed state，归一迹 tr_m，L²=Etr X†X。先有限N，不取稀薄极限。以下 m>p≥1，d=m−p，n=N−m+1。

沿用 C=w binom(m,p)binom(N−m+p,p)、C_-=w binom(m−1,p)binom(N−m+p+1,p)、jc=w binom(m−1,p)binom(N−m+p,p)，以及

\[
\kappa_r=w\binom{m-r}p\binom{N-m+p-r}p,
\quad e_r=\frac wp\binom{m-r-1}{p-1}\binom{N-m+p-r}{p-1}.
\]

越界组合数取0；内部 e_r=pκ_r/[(m−r)(n−r)]。令 Ψ(X)=E uXv，其中 u,v 为相同新耦合在两个charge的Hamiltonian。对任何矩形 Y_i:H_(m−1)→H_m：

\[
P_r\sum_i\Psi(Y_i)c_i
=\kappa_r P_r\sum_iY_ic_i
-e_rP_r Rem\left(\sum_i c_iY_i\right).
\]

证明：primitive T_r 的 V_i=T_r^+c†、U_i=c†T_r^- 只含 charged (r,r−1)/(r+1,r) 两型；V−U 是低型。用未归一HS，A_-/A_+=(m−r)/(n−r)，||V−U||²=r(N+1−2r)A_+/(n−r)，〈V,V−U〉=rA_+。低型投影为(n−r)(V−U)/(N+1−2r)。两型旧协方差特征值 L_r/U_r 满足 L−U=(N+1−2r)e_r、U+(n−r)e_r=κ_r。得到 Ψ(V)=κV−(n−r)eU。再取 S(Y)=ΣYc、T(Y)=ΣcY 的HS adjoint，因 Ψ自伴，推出任意Y的上述式。r0及边界用存活型或binomial延拓。

对于真正 p-body 输入 ranks：

\[
\Psi(h_ac^\dagger)=\kappa_a h_ac^\dagger-\frac{p\kappa_a}{d}c^\dagger k_a,
\]
\[
\Psi(c^\dagger k_a)=-de_a h_ac^\dagger+[\kappa_a+(m-n)e_a]c^\dagger k_a,
\quad\kappa_a(m-1)=\kappa_a+(m-n-p)e_a.
\]

## 完整 Q 的两个振幅

Z⁺_ab;r=P_r(h_a⊙h_b)，Z⁻_ab;r=P_rRem(k_a⊙k_b)，不同rank用反对易子，同rank用平方；x_a=a(N+1−a)。直接将transfer应用到旧五字收缩，得到

\[
Q_{Ecal}=\sum_{a\le b,r}(g^+Z^++g^-Z^-),
\]

\[
g^+=(8p-2m)C_-+\frac{19}2djc+(12p-9m)C
+m(\kappa_a+\kappa_b)+\frac{3p-m}2\kappa_r-\frac{11}4dn(e_a+e_b),
\]
\[
g^-=3C_-+\frac52C-4jc-\frac{4m+p}{4d}(\kappa_a+\kappa_b)+\frac12\kappa_r
-(m-n+p)(e_a+e_b)
+\left[m+3n-p-\frac{5(2mn-x_a-x_b)}{4d}\right]e_r.
\]

关键二次项 P_rΣΨ(h_ac†k_b)c=dκ_rP_r h_ah_b−e_r(m−a)(n−a)P_rRem(k_ak_b)/d。其他项只需线性恒等式。五个单收缩各自的系数表由 memory_two_charge_filter.py 保存；symbolic-filter.json 独立合并并证明与以上简式相等。

W2Ecal=Σ(g⁺W2Z⁺+g⁻W2Z⁻)，因为 E Q=2E Ecal。m=p 边界 k0、Ecal=p h4、Q=pD40，另行处理。不会除以为0的 Z⁺，上一轮 odd-output 反例完整保留。

## 双标度中合并大项

用 A_ab;r=Z⁻+(2p−m)Z⁺，记 g_A=g⁻、g_H=g⁺+(m−2p)g⁻。固定半填充 m=N/2，N=p²/λ，λ>0；u=p−a、v=p−b、S=a+b−r 固定，q=e^(−4λ)。

\[
g_A/C=(3-4q+q^2)/2+O(p^{-1}),
\]
\[
g_H/C=2p(1-q)^2+K_0+8\lambda(q^2-q)(u+v)+8\lambda q^2S+O(p^{-1}),
\]
\[
K_0=\lambda[1-10q+9q^2+16\lambda q-32\lambda q^2].
\]

O(m)项已在有限式合并后严格相消；剩余O(p)系数独立于rank，对应守恒 W2H²=H²−C。

证明使用 g_H/C=B0+Baκa/C+Bbκb/C+Brκr/C 的精确有理系数。固定deficits下为 (2p+λ,−2p−5λ,−2p−5λ,2p+9λ)+O1/p。binomial product 的 log 在 j=tp−s 时为 −4λt+(4λs−4λ²t²)/p+O1/p²。所有固定窗口余项由log Taylor界控制，不是数值拟合。

## H² 部分的 L² 结论

全比例窗口 a=αp,b=βp,r=ρp、αβ∈[0,1],ρ∈[0,2]，四个有理系数一致展开为

- B0=2p+λ+O1/p；
- Ba=−2p+λ(6−11α)+O1/p，Bb同β；
- Br=2p+λ(5α+5β+6ρ−13)+O1/p。

缩放后的有理分母在紧集远离0；κ相邻rank差/C≤Kλ/p。因此 |g_H/C−2p(1−q)²|≤Kλ(1+u+v+S) 一致。

不同input pairs在Gaussian二次chaos正交，不同output ranks在HS正交。记 Ca=E h_a²，ωa=Ca/C，则

\[
\left\|C\sum f_{ab;r}W_2Z^+_{ab;r}\right\|_2^2
\le4C^4\sum_{a,b}\omega_a\omega_b\sum_r\pi_{ab}(r)|f_{ab;r}|^2.
\]

旧input factorial tail及S三阶矩给右边线性滤波平方的一致可积。先截u,v，再用S的O1/K平方尾；输入窗外由factorial尾和条件S²一致界控制。结合固定窗口收敛，物理 C→Js²/λ 下：

\[
\mathcal H_N:=\sum g_HW_2Z^+
=C[2p(1-q)^2+K_0]W_2h^2
+8\lambda C(q^2-q)W_2\{h,D_N\}
+8\lambda Cq^2\mathcal S_N+o_{L^2}(1),
\]

D_N=Σ(p−a)h_a，S_N=Σ(a+b−r)W2Z⁺。这只证明 H² 部分、精确半填充。A部分系数虽然固定deficit收敛，尚需A本身的加权尾，不能对两电荷各自O(m)范数用粗界后声称误差消失。全charge trace也须补中央窗。

已有 cubic-linear 工作将中心化 D_N 接到 F；后续可使用 F≈8λq C(D_N−λh)，但应在乘积层保留相应L²控制。新问题聚焦 S_N 与 A/JhF 的 norm/cross。

## 验证与下一步

五系统 (4,1,2)/(5,2,3)/(6,2,3)/(6,2,4)/(6,3,4)：522个精确rank单收缩/总和恒等式；另对全部rectangular matrix units验证87个完整linear maps，共5,470,368 complex entries，最大清分母整数60。3个符号投影式通过。后者是普适transfer的有限全map核对；前者只是选取纯rank realization的quadratic核对，不冒充全quadratic tensor穷尽。

5份JSON均passed：filter-checks、transfer-checks、symbolic-filter、limit-checks、uniform-prefactors。24组有限λ=.1/.25/.5、p16..128、80/120位核对保留实际偏差；λ.5,p128,u2v1S3时gA差.024866，gH有限部差−.021443，非区间认证。全部首次计算通过，无失败或阈值调整。六页双遍编译及最终逐页95dpi检查通过。

下一步：半填充 P_r(V_bU_a)=(−1)^S P_r(U_aV_b)，所以 sym强度=4·1_(S even)π。需求奇偶加权S矩，而前三个普通矩不够；同时处理下电荷和cross Gram，证明 g_A 可在算符和中取常数。之后恢复 W4 后代low chaos、求raw Z2完整谱，再继续finite-memory和体空间/质量字典。整体goal未完成。
