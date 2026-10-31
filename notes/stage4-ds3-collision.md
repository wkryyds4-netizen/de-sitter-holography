# dS₃ 交叉核的可去奇点与受控重质量因子

2026-10-31。Alan Day / Chengheng，内部自查。[五页PDF](../output/pdf/de-sitter-ds3-collision.pdf)、[LaTeX](../manuscript/ds3-collision/main.tex)、[来源](../literature/stage4-ds3-collision-sources.md)。接续[同species Bose](stage4-ds3-bose.md)，不改写旧科学结果。

## 模型与目标范围

继续ℓ=1、G=0、proper SO⁺(3,1)、BD标量观察者与KG物质、固定观察者下界s>0及半周期创建块。研究的是已有物理U_00在完整Racah变换中的矩阵元，不换态、不加新的内部滤波。完整Bose自由场极限仍是目标，本轮先修复可去奇点并证明一个关键Omega因子的重质量控制。

## 三个实质量面与合并级数

LPRS(3.38)的每个Omega为Γ(z)A(z)+Γ(−z)B(z)，z=H−h′+h₁+h₃−1，H=h或1−h，h′可取shadow。正实质量的碰撞面是a+c=b+d、d=a+c+b、b=a+c+d。原Euler乘积积分在端点为x^−1/2(1−x)^−1/2乘至多对数因子，局部可积；这些面不是有限实质量的物理极点。

固定H、h_i，h′₀=H+h₁+h₃−1，h′=h′₀−z。四个上参数A=(H+h₁+h₂−1,H+h₃+h₄−1,H+h₃−h₄,H+h₁−h₂)，ReA_j=1/2。Gamma约分后得到Omega=C(z)Σ(GA_n−GB_n)/z，C与两组Gamma比的完整公式见PDF§2。每项原来n^−2的首项逐项消去，合并后为2C(z)n^−3+O(n^−4)。固定参数的局部一致收敛允许在z=0取导。

精确碰撞处不用挪动质量，直接为

\[
\Omega_0=\sum_{n\ge0}t_n q_n,\quad
t_n=\frac{\prod_j(A_j)_n}{(2h'_0)_n(2H)_n(n!)^2},
\]
\[
q_n=\psi(n+2h'_0)+\psi(n+2H)+2\psi(n+1)-\sum_j\psi(n+A_j).
\]

这与独立的A′(0)−B′(0)−2γA(0)一致。非零z用GB·expm1(logGA−logGB)/z稳定合并，保持Gamma相位。所有这些均是已知LPRS公式的代数重组，未作原创特殊函数声明。

## 一个严格的重质量窗口

取c=M、b=M+u、a=M+u+v、d=M+v，H=(1−ib)/2。假设|u|+p≤1、|v|+q≤1及所有观察者质量≥s。A_j=1/2+iα_j，此时|α_j|≤1/2，m=min(b,d)>0。

由|n+A_j|²≤n²+n+1/2和(x²−x+1/2)(x+1)≤x³（x≥1），逐项得到

\[
|t_n|\le\frac1{(n+1)^2}\prod_{j=1}^n\frac{j^2}{j^2+m^2}
\le\frac1{(n+1)^2(1+m^2)}\quad(n\ge1).
\]

用digamma差级数，|q_n|≤B_M=2+asinh b+asinh d+4log2 Σ|A_j−1|。于是

\[
\Omega_0=q_0+R_M,\qquad
|R_M|\le\frac{(\zeta(2)-1)B_M}{1+m^2}=O(\log M/M^2).
\]

q₀=ψ(1−ib)+ψ(1+id)+2ψ(1)−Σψ(A_j)，故Omega₀=2logM+2ψ(1)−Σψ(A_j)+O(M^−1)。这只是一项Omega，不是完整6j或Bose范数。

求和n=0,…,N−1的严格绝对级数尾为

\[
\left|\Omega_0-\sum_{n<N}t_nq_n\right|
\le\frac{B_M}{N}\prod_{j=1}^N\frac{j^2}{j^2+m^2}.
\]

它是精确算术中的遗漏级数界，不包括mpmath舍入误差。M32,N128样本约2.61703e−39；更高M时尾可以小于现有算术误差很多，不能拿这个数字当总误差。

## 邻近面与对数尺度

d=M+v+r，z=−ir/2。若A_j及A_j−z的虚部模均≤1/2，两组4F3各自具有相同类型的比值界。其n0近似S=Γ(z)a₀(z)+Γ(−z)b₀(z)保留完整Gamma因子，显式a₀,b₀见PDF§4。定义m_A=min(|M+v+r|,|M+u−r/2|)、m_B=min(|M+v+r/2|,|M+u|)，则

\[
|\Omega-S|\le(\zeta(2)-1)
\left(\frac{|\Gamma(z)a_0|}{1+m_A^2}+\frac{|\Gamma(-z)b_0|}{1+m_B^2}\right).
\]

各自的N项尾用P_N(m)/N再乘对应Gamma系数。r非零但趋零时这个直接界含1/|r|；r=0用上节加权级数。

在转移窗口严格内部，对每个固定λ，取r=λ/logM。Gamma比给a₀∼(−iM)^z、b₀∼(iM)^−z，结合Γ(±z)=±1/z+O(1)，得到点态极限Omega/(2logM)→sin(λ/2)/(λ/2)，λ=0取1。λ≠0的级数误差O(logM/M²)，λ0由前节控制。尚未证明对λ一致，更不是全U的分布极限，不能在物理积分内直接换成δ核。

M1024仍不在一个很好的领先近似区：λ0的scaledOmega≈1.31857+.08411i，离1的绝对差.32949；λ1的scaledOmega≈1.27753+.09388i，离sinc(.5)差.33222。保留慢收敛，没有只展示漂亮样本。

## 数值与修正记录

[初始扫描](../results/stage4/ds3-collision/study.json)保留M2,4,8,16与r0,.1,.4，共12点双精度层次（50/80十进制位）。完整U峰上参考约.96379723、1.15801422、1.36981862、1.58731944；另一个M32高精度参考约1.80680785。点值不代表已完成物理质量积分。

[诊断](../results/stage4/ds3-collision/validation.json)三面、Hermiticity、近面与独立Euler共9项通过；Euler相对差约3.04e−26。初始prototype的digamma Bernoulli系数先用Python float做1/k，导致一处约2.81517e−27共同偏差，单纯升dps无法揭示。原代码和输出原字节保留，[修正模块](../code/ds3_racah_merged.py)改全mp系数，与独立参数导数相差约9.74e−58。

完整M32核在110位算术下，(N,K)=(512,24)相对错2.42e24；(1024,32)错5.78e6；(2048,40)、(4096,48)、(8192,56)依次约1.74e−13、5.39e−19、8.60e−25。全部失败设置保留，不能据最后一个点宣传统一稳定算法。说明可去极点修复、舍入控制与大参数渐近余项是三件不同的事。

[重质量界诊断](../results/stage4/ds3-collision/asymptotic.json)含两个转移组、M2至64共12点，均位于解析界内；[受控直接和](../results/stage4/ds3-collision/controlled.json)含M16至1024、λ0/1，70/100位比较及独立的解析尾。[额外验证程序](../code/validate_ds3_controlled.py)正在与完整hyper/参数导数做100位比较，session 47841仍活跃，controlled-validation.json尚未生成，不能计作通过。数值库尚未区间认证。

下一步要把其余Omega、Gamma前因子、shadow两项及物理归一一起控制，再对真实创建幅作质量积分和谱尾估计。不能把单因子的点态sinc极限当完整Bose自由极限。全Fock、有限G和微观全息字典仍开放。
