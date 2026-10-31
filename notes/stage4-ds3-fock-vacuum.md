# 受调节 Fock 场、真空返回与完整四点

2026-10-31。正式推导见[九页阶段稿](../output/pdf/de-sitter-ds3-fock-vacuum.pdf)及[LaTeX](../manuscript/ds3-fock-vacuum/main.tex)。这是指定模型内的自查推导，使用已有 CG、群平均、Bose 对称化和解析向量定理，不作原创性或完整全息对偶声明。

## 两种调节必须分开

取 ℓ=1、G=0、proper SO⁺(3,1)，两个可区分 BD 标量观察者，标签下界 s>0；同种自由 KG 物质 p>0，耦合为1。观察者静能量 h(a)=√(1+a²)，不是标签 a。

半周期创建块的乘子为 exp(πΔh/2)exp(−itΔh)。创建加 Hilbert 伴随得到自伴场，但这不是将原来的完整 Hermitian 场同向平移到虚时间。另取原场的实 Gaussian 时间涂抹，创建乘子为 exp(−σ²Δh²/2)exp(−itΔh)，σ>0；湮灭仍为其伴随。本轮同时保留这两个不同 observable，不能用前者替代后者的热性问题。

## 全粒子数的构造与弱下降

单观察者空间为 ∫a≥s da P_a。沿旧局部顶角规范，ρ(b)=2b²/π、w_a(b)=ρ(b)K₃(a,p,b)，标量通道 CG synthesis J_a 从 ∫db P_b 等距嵌入 P_a⊗P_p。实际局部顶角确定

\[
(T_j(t)\xi)(a)=J_a\left[ E_{abp}m_j(a,b;t)\xi(b)/\sqrt{w_a(b)}\right].
\]

这里是 L² 谱综合，不把固定标签的分布 CG 列当可平方求和向量。半周期率 R_h=|E|²exp(πΔh)/w=(g(a)/g(b))R̂_p，g(x)=x exp(π(h(x)−x))/(1+exp(−2πx))，R̂_p≤k_p(a−b)，k_p(u)=1/[8π(coshπu+coshπp)]。

令 B_s=2exp(π/(√(1+s²)+s))，则 ∥T_h∥²≤𝓑_h=B_s(1/(2π²)+1/(2sπ³))。Gaussian 率满足 R_σ=exp(−σ²Δh²−πΔh)R_h≤Λ_σ R_h，Λ_σ=exp(π²/(4σ²))，故 𝓑_σ=Λ_σ𝓑_h。界不在 s↓0 或 σ↓0 一致。

在 W_N=H_O⊗SymᴺP_p 的完整主系列分解中保留全部允许自旋与内部树路径；Bose 投影作用于完整重数空间。标量 spectator 只选出总表示的标量部分，不删除内部非零自旋。W_N 的 SO(3) 不变部分可写为 ∫dc M_N,c。N=0 单独是 c≥s 上的一维重数，不增加独立总谱标签。

种子取 f_L(c)|0〉⊗v，f∈Cc∞((s,∞))、v∈W_N^K。球函数 φ_c(r)=sin(cr)/(c sinh r)，统一界 r/sinh r。先积分左质量，紧光滑 f 的 Fourier 衰减与右侧统一球函数界、Haar 的 sinh²r 合并给 C_n r(1+r)^(−n)，n>2 可积。因此不需要对两个固定标签球函数的绝对 Haar 积分作错误假设。得到

\[
Q_N(f\otimes v,k\otimes w)=\int dc\,\rho^{-1}\bar f k\langle v,w\rangle,
\quad A_N(f\otimes v)=\rho^{-1}fv,
\quad H_N=\int_{c\ge s}^{\oplus}\rho(c)dc\,M_{N,c}.
\]

紧 c 支持 F 选 f=1、v=ρF 证明像稠密；固定 f 的有界估计允许有限角模、紧实际观察者质量逼近。这个结论针对明确的群平均完备化，不断言未指定的最大位置分布空间。

创建块 C_j,N=√(N+1)(1⊗P_Sym)(T_j⊗1) 为有界等变映射，∥C_j,N∥≤√((N+1)𝓑_j)。用共同紧群截取，先有限模式强极限，后去群截取；种子空间的绝对群平均给零范数代表独立及反向块的伴随。完整内部重数保留，不能把 scalar 子块当完整酉重耦合。

在 D_fin=⊕alg H_N 上，Φ_j=C_j+C_j* 对称且保域。若 ψ 粒子数≤r，则

\[
\|\Phi_j(t)^k\psi\|\le(2\sqrt{\mathcal B_j})^k\sqrt{(r+k)!/r!}\,\|\psi\|.
\]

比值检验给整解析向量；Nelson 定理给本质自伴。不同实时间的有限乘积也在此共同域存在。各扇区块有界，完整场并非有界算符。未确定代数的 factor type。

## 真空返回不能遗漏

归一真空制备满足 ∫ρ|f|²=1。返回函数 q_j,c(τ)=∫s∞da R_j(a|c)exp(−iτ(h_a−h_c))，且 C_j,0(t₁)*C_j,0(t₂)=q_j,c(t₂−t₁)，对 conserved spectator c 对角。

四场真空期望只有 01010 和 01210 两条粒子数路径：前者是 ∫ρ|f|²q_c(t₂−t₁)q_c(t₄−t₃)，后者为两创建态的 Bose 内积，已含 direct/exchange。前者共享同一个 c，不能改成两个分别质量平均的二点乘积。

\[
D_{\Phi,j,f}(t)=D_{B,j,f}(t)+4\int dc\,\rho|f|^2(\Im q_{j,c}(t))^2.
\]

这是零粒子和两粒子扇区正交的精确恒等式，不是近似 Wick 分解。单步 |a−c|>D 的解析尾为 B_s exp(−πD)/(2π)×[1/π+(D/π+1/π²)/s]，Gaussian 乘 Λ_σ；不含求积和舍入误差。

## 重质量极限及热性修正

固定 s,p,σ>0，条件 c=M→∞，紧实时间：q_j,M=L_j+iL_j′/M+O(M⁻²)。半周期谱为偶函数 k_h=k_p，L_h(t)=sin(pt)/(4πsinhπp sinh t)；Gaussian 为 k_σ(u)=exp(−σ²u²−πu)k_p(u)，L_σ 为其 exp(−itu) Fourier 变换。

Gaussian 满足 k_σ(−u)=exp(2πu)k_σ(u)、L_σ(t+2πi)=L_σ(−t)。因 τ=t₂−t₁，通常 W(t)=L_σ(−t) 满足下条带 KMS。这是极限热性，未证明有限 M 的精确热性。半周期谱在非零频率不满足这个热详细平衡。

旧任意紧转移核定理适用于一般四时间的紧 C² 相位族。Gaussian 每个 emission 相对半周期乘子的模≤√Λ_σ，两创建链幅≤Λ_σ，平方范数及 Hilbert 交叉尾界乘 Λ_σ²。先 M 后截取 K，恢复三个标准收缩：

\[
G_{4,j}^{\infty}=L_j(t_2-t_1)L_j(t_4-t_3)
+L_j(t_3-t_2)L_j(t_4-t_1)+L_j(t_3-t_1)L_j(t_4-t_2).
\]

也覆盖归一质量概率逃向无穷的制备；不覆盖固定不动的初态包或增长时间。两种调节均有 D_B→0，但完整场 D_Φ→4(Im L_j(t))²。半周期极限为零，Gaussian 一般非零。这纠正把两创建次序差解释为完整热场不同时刻对易的潜在误读。

半周期返回另有 q_h,M=L_h+iL_h′/M−(t+iπ)L_h′/(2M²)+O(M⁻³)。这不确定总 D_Φ 的领先速率，因为 D_B 尚无优化全局速率。全 N 存在也不等于全部关联函数的 Gaussian 极限定理。

## 计算、独立复核与失败

[结果目录](../results/stage4/ds3-fock-vacuum/)保存参数、依赖和哈希。半周期 p=.6、s=.2，M=.7/4/16/64/256，t=0/.5/1/2，60位，每组 D24/32，共20组。五个原始 E/w 核对相对差最高约2.01e−58。t=1 的真空贡献由 M=.7 的2.47376e−5降至 M=256 的1.61887e−9。

辅助模型有3个观察者能级、2个 oscillator modes、非对易复矩阵及 cap4；验证4个√N界、4个四场路径和3个正交分解，不冒充 dS 完整自旋空间验证。原 driver 已保存全部求积，之后因 np.bool_ 的 JSON 序列化失败；另 recovery 复用 raw 并显式转换标量，36项通过。原脚本及失败不覆盖。

原 independent audit 的153项仅失败于 SymPy 结构 ==；另以 expand(lhs−rhs)==0 验证同一二阶恒等式，154项通过，不改数据或阈值。独立96/160阶 Gauss 对60位结果最大差4.34e−17，加阶最大差8.25e−17；不是总误差认证。

Gaussian p=.6、s=.2、σ=.8，M4/16/64/256及同四时间，共16组双截取；16返回、4 KMS、4详细平衡共24项通过。最大 KMS 残差1.09e−62。t=1 的极限真空贡献 .00411327832，M256 为 .0040762936。另115项独立审计通过，96/160阶 Gauss 对参考最大差1.31e−16。这115项保存在附属证据中，未回写本次已编译 PDF 的数值段。

无新 U 核求值、无完整四场全谱数值求积；四点极限依赖解析链，不能用检查项数量替代证明。

## 下一研究入口

回到有限 G 插入和微观字典比较，先读[tqft-recoil-dictionary](../literature/tqft-recoil-dictionary.md)、[microscopic-match](../literature/stage4-microscopic-match.md)、[October TQFT update](../literature/october-tqft-update.md)。保留真实 Gaussian 调节、真空返回和有限质量制备；先核有限 G 的态、测度、轮廓和两点归一，再判断有无能比较的四点 observable。G=0 完整 Fock 构造本身不确定有限 G 动力学、type II 或微观全息对偶。
