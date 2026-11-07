# 中央电荷窗口：产品传递与输出缺额

2026-11-06。阶段稿见 `manuscript/susskind-memory-central-cross/main.tex`，六页阅读版见 `output/pdf/de-sitter-susskind-memory-central-cross.pdf`。这一轮只本地保存，不 commit、push 或 ZIP。

## 这次接通的部分

沿用完整重叠 complex Gaussian p-body 系综、同一组 J、归一最大混合固定电荷迹、普通第二 Gaussian chaos W2、hbar=1。N 为偶数，p²/N→λ>0，m=N/2+xp，x 在固定紧集内，λ 在正紧集内。C_m 保持实际扇区方差，不能在放大反项里换成 C_half。

已有半填充输出缺额范数，现在可以通过相邻扇区等距映射推广到整个中央紧窗口。B=W2H²，L_in=W2{h,D−bar u h}，S_m=Σ(a+b−r)W2Z_ab;r。这里 L_in 是输入秩缺额，不是此前三阶记忆场 L。

设 q=exp(−4λ)，b=2λ(1−q)，s=4λ²(1+q)+2λ(1−q)，则 G(B,L_in,S)/C_m² 的极限为

    [1+q, 0, b]
    [0, 2λ(1+q), 0]
    [b, 0, s]

这是本轮已完成的解析结论。它利用已知半填充 Gram 作参考，不在一般 charge 中假定乘积奇偶恒等式。

## 传递证明

n=N−m+1，x_r=r(N+1−r)，ell²=m(mn−x_r)/n。Rem X=Σc†Xc 在 rank r 上把下扇区归一 HS 范数乘 ell；I_m=Σell_r^−1 Rem P_r 因此是公共 rank≤2p 算符空间上的等距映射。沿中央路径，充分大的 p 保证两侧都包含这些 rank。逆向路径用其 HS 伴随。

每个扇区都嵌入同一个 p-body Gaussian 耦合空间，输入 rank a,b 的标号一致。I_m 只作用输出矩阵，保留 a,b,r，与 W2、S 以及输入权重交换。无需把它解释成正量子信道。

R_m=C_(m−1)/C_m=(m−p)(n+p)/(mn)，t_r=R_m ell_r。令 Bhat_ab;r=W2Z_ab;r/C_m，则真实同 J 的 A=Rem(k²)+(2p−m)h² 给精确式

    t_r I_m Bhat_(m−1) − (m−2p) Bhat_m = O_ab;r/C_m.

中央窗口、r≤2p 上 t_r≥cN，|t_r−(m−r)|≤K，所以 |t_r−(m−2p)|≤K(1+u+v+S)。上一轮实际 O 的有界范数和正 product 矩给一步误差 K/N。共 O(p) 步，归一 B 误差 O(1/p)。标量符号检查还给

    t_(rho p)−(m−rho p) → 2λ[x(2−rho)−1],
    p²(R_m−1) → 4λ²(2x−1).

S 的处理：统一截断 S≤R 后先 telescope，再乘 S，低部分误差 KR/p。两端的正六阶矩给 ||S 1_(S>R) Bhat||≤K/R²。只算两个端点，不沿每一步重复累加尾项。因此

    ||U_(m←half)(S_half/C_half)−S_m/C_m|| ≤ K(R/p+R^−2).

取 R=p^(1/3)，得到 O(p^−2/3)。这是扇区间传递速率，不是整个有限 p Gram 到其极限的收敛速率。输入权重版本使用已有 Σ(1+u+v)²||O_ab;r||²/C_m²≤K；bar u_m→λ 一致，因此 L_in 也可传递。三方向 Gram 随之成立。

## 剩下的交叉项：检验通过，解析一步仍待补

完整四方向 (B,L_in,S,O) 只剩 S–O 的一般 charge 极限。本轮记录了可用的候选公式，并且没有把有限系统核对当作普适证明。

χ 是 neutral Racah，ξ 下参数为 (1−m,m−N,−N−1)，D=ξ−χ；g=(m−p)²/m，η=p²/m，Cf=m(N−m)(N+1)/N，Ct=(m−1)(N−m)(N+1)/(N−1)。候选 coherent exchanged Casimir moment 为 g[(x_a+x_b)ξ−mD]，direct 为 g(x_a+x_b−x_ax_b/Ct)。相应 X/N–O 的稳定 pair 式为

    z_ab={η[−2ab(1+χ)+x_ax_b/Cf]
          +g[x_ax_b/(m(m−1)(N+1))+(m−2ab)D]}/N.

候选中央极限是 s，与 half 一致。X=x_(a+b)−x_r=S[N+1−2(a+b)+S]，已证正矩给 ||S−X/N||/C_m=O(1/p)，已有 O/C_m 有界，所以代理转实际 cross 没有新的尾缺口。

关键尚缺：从 Casimir product rule 及 [E_ij,c_k]=−δ_ik c_j，把 exchanged gradient contraction 解析化成 g x_b[T_a(b−1)−T_a(b)]/(N+1−2b)。后者等于 −gmD 的 terminating-series 恒等式已经核实，但这并不自动完成前面的矩形算符收缩。不能用半填充 alternating completeness 代替这一步。

如果补齐后 c=s，则完整四 Gram 的 Schur complement 给 S 扣 B/L_in/O 后的方差 vΔ/(v+Δ)>0，其中 v=2λ(1−q)+16λ²q/(1+q)，Δ=2λ(1+q)。本阶段只将这写成条件结论。

## 验证记录

- `transport.json`：标量符号式；N6p1m3 完整 rank 0/1/2 等距矩阵，维数1/35/189；1296个有序 Gaussian 坐标对。sqrt 数值等距差约2.22e−16，传递恒等式最大差3.47e−17。统一尾/telescope 证明在稿中，不把小系统当作极限证书。
- `charged-casimir.json`：15个纯 rank pairs、60条有理恒等式、66531个有序坐标对，涵盖四个系统。没有 half parity 替代。
- `car-exploration.json`：四套完整4×4 CAR Gram，共31466个对称坐标对，candidate X–O差都精确为0。保留 exploratory 状态。
- `cross-limit.json`：624条精确 contiguous 恒等式、四组有限完整rank核对；λ=.1/.25/.5、p16/32/64/128、x−2/0/2，共36组80/120位完整终止和。输出存65有效位，精度变化记录0仅表示这65位相同，不代表任意精度或区间认证。第三列仍是候选cross的函数求值。

λ=.5,p128 的候选cross在 x−2/0/2 分别约1.97682736、1.94517151、1.91541182，目标2；有限尺寸误差没有隐去。

双遍编译使用 `powershell -NoProfile -ExecutionPolicy Bypass -File code/build_memory_central_cross.ps1`，渲染 `pdftoppm -r 95 -png ...`，最终六页逐页查看。首版 LaTeX 行尾转义遗漏导致 align 多 label 和溢出，首稿/失败日志保留 diagnostics；修正后六页无引用警告与公式溢出。输入、输出、渲染和前阶段文件哈希由 local-evidence 保存。

下一步先补一般 charged gradient 收缩，再把四 Gram 接中央能量插入与 raw U2 谱。全电荷尾、有限 memory 两插入余项及体介子空间/质量字典继续在后。整体目标仍 active。
