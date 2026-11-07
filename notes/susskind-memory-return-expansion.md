# 有限记忆时间：两次返回的受控归约

2026-11-06。六页正式推导见 `manuscript/susskind-memory-return-expansion/main.tex`，已双遍编译为 `output/pdf/de-sitter-susskind-memory-return-expansion.pdf`。

## 结果与物理位置

沿用全重叠 complex Gaussian p-body H、标准 CAR、共同 J、β=μ=0 全退火迹、ħ=1/string time。偶 N，λ_N=p²/N 在正紧集，固定 Js、固定有限记忆时间窗 |s|,|u|≤S。

上一轮实际对象 Fscr_m(s)=(1−P_h)[m exp(−ihs/2) M_m(s) exp(−ihs/2)] 继续保持；P_h 是逐样本完整能量块 pinching。令 Λ(s)=R†exp(iKs)R、U(s)=exp(ias)，a=(1−p/m)h。保留 J0=Λ、J1=Λ*U*Λ、J2=Λ*U*Λ*U*Λ 后：

- 集体场的退火 Hilbert–Schmidt 误差为 O_S(p^-2)。
- 两插入协方差的绝对误差为 O_S(p^-1)，对所有实外部中心时间 T 一致。
- 协方差可进一步写为 n+l≤2 的六个实际共同 J、逐样本 pinched 产品加 O(p^-1)。完整九项代理保留正性。

这是有限时间的归约：指数保持完整，未截为 s 的有限次多项式。尚待计算六个产品的极限，证明有限 s 的紧性，随后连接 Susskind 体介子空间/质量字典。

## 有限 Volterra 恒等式

分块 K=[[a,R†],[R,D]]，D=QKQ|Q。消去 P 块的两次跃迁带 i²=−1，得到

Λ=M−M*U*Λ，M=Λ+M*U*Λ。

迭代固定 j 次后 M=Σ_(n=0)^j Jn+Ej，Ej=M*(U*Λ)^{*(j+1)}。未假设无限级数收敛。使用同一源维数 d_m 归一矩形 Schatten 范数，2j+4 个 R/R† 因子加 2j+2 维时间单纯形给

||Ej(s)||₂ ≤ |s|^(2j+2)/(2j+2)! · [tr_m G^(2j+4)]^(1/2)，G=R†R。

相似地 ||Jn||₂≤|s|^(2n)/(2n)! [tr_m G^(2n+2)]^(1/2)。负时间由 adjoint 与反向单纯形保留。标量 K=σx、V=e1 的独立检查是 M=1、J0=cos s、J1=s sin s/2、J2=(s sin s−s²cos s)/8。

## 用 Gaussian square functions 补齐高矩

改用比早先设想的高阶 marked-chord 枚举更直接的路线。Fi=[H,ci]，A=ΣFi†Fi；单项 [cI†cJ,ci] 只有 i∈I 才非零。计数完整 I∩J 重叠，得到本稿 (16)–(18) 三个占据态对角协方差。其统一上界为

E Fi†Fi、E FiFi† ≤ v_N I，v_N=4(Js²/λ_N)(p/N)(1+p/N)^(2p−1)≤4(Js²/λ_N)(p/N)exp(2λ_N)。

矩形 Fi 的 Hermitian dilation 是 Gaussian matrix。分部积分、在 X 本征基中的加权算术几何平均给

E Tr X^(2k)≤(2k−1)v_N E Tr X^(2k−2)。

停在二阶矩，避免用增大的目标空间维数归一。再用退火 Schatten k 的 Minkowski，得到所有 sector 的

E tr_m A_m^k≤(2k−1)!! (N v_N)^k=O_k(p^k)。

不同 i 的力无需独立。由 G≤A/m 的本征值单调性（不是声称矩阵幂算符单调），分 m≥N/4 和 m<N/4；active m≥p，低尾概率≤exp(−N/8)。因此

B_k,N=Σ b_m E m² tr_m G^k≤(2k−1)!!(N v_N)^k[(N/4)^(2−k)+p^(2−k)exp(−N/8)]=O(p^(4−k))。

本阶段只用 k=2,4,6,8。无需上一轮增长窗的额外推广，也不交换无限时间与 N 极限。

## 为什么保留两次返回

现有粗界 ||Fscr||₂=O(p)。只保留一次返回时，场误差 O(1/p) 乘 O(p) 仅给协方差 O(1)；不能把场误差趋零直接当作协方差误差趋零。两次返回场误差 O(1/p²)，Cauchy–Schwarz 才给 O(1/p)。这不预设实际场紧性。

## 可计算的普通传播公式

F(s)=V†exp(iKs)V，则 Λ=−F''+i{a,F'}+aFa。保留精确 a。源能量 Ea、下电荷能量 εα 给

Wα_ab=(1/m)Σi(cαa^i)*cαb^i(εα−aa)(εα−ab)，Λ_ab=Σα Wα_ab exp(iεαs)。

J1/J2 对应三/五个能量节点的有序路径。时间积分用 confluent divided difference；数值用对角 iξ、超对角1的双对角矩阵指数右上角，避免重合能量相减除零。

L=J0+J1+J2 的逐样本代理协方差：

(m²/d_m) Σ_(Ea≠Eb) exp[iT(Ea−Eb)+i(s−u)(Ea+Eb)/2] conj(L_ab(s)) L_ab(u)。

简并能量整个块去掉。下电荷与源 H、W、pinching 同属一个样本，不能用独立探针因子化替换。

## 本地检查

`code/memory_return_expansion.py` 首次运行 exit0；`finite-checks.json`：864精确整数 CAR 协方差对角项（同时检查所有非对角零）、118132有理计数界、4实际共同 J 系统×4正负时间×3返回阶数=48余项界、导数/adjoint/协方差相位、单独谱路径及5标量检查。N3p1动态 pinched 协方差为数值零，另外三个 p2 系统非零；不以平凡对照代替动态覆盖。

传播测试用带种子的复整数 Hermitian couplings，再除以全 H 算符范数；“Gaussian integer”指复整数，并非 Gaussian 分布抽样。这些是有限恒等式检查，非 DS 拟合。complex128 与经验舍入 allowance 如 JSON/稿内所列；渐近阶由解析不等式得出。

首次排版有3处长行溢出，日志保留于 diagnostics/first-layout.log。修订后双遍编译无 undefined/overfull，六页95dpi已逐页检查。既有科学输出保持不动；本轮无 commit/push/ZIP。
