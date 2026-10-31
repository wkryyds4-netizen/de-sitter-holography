# FP 物理内积与 multiplicity 直积分：已核校准及下一步

日期：2026-10-29。本稿记录本轮末的短路线判断，不启动新计算。来源为 [Kolchmeyer–Liu 2411.08090v1](https://arxiv.org/html/2411.08090v1) 的 (4.22)、(4.31)—(4.33)、(4.42)、(4.44)—(4.50) 与附录 E。顶角归一化采用[独立重建](stage3-vertex-analytic-audit.md)，不采用已发现不相容的 literal A.7。

## 1. N=0 与 N=1 的测度校准

令

\[
\rho_{\rm Pl}(d)=\frac{d\tanh\pi d}{4\pi^2},\qquad
w(d)=\frac{\pi d\tanh\pi d}{4}=\pi^3\rho_{\rm Pl}(d).
\]

E.22 的 raw principal CG 满足 C†C=[1/w(d)]δ(d−d′)，因此

\[
\widehat C^{d}_{a,\nu}=\sqrt{w(d)}\,C^{d}_{a,\nu}
\]

是采用 dd 谱测度的归一化 CG。E.22 对 S、A 两通道给相同自正交常数；通道完备性须连同 E.1.3 使用。

N=0 的物理空间已由原文确定为 L²(ρ_Pl(d)dd)。N=1 至少在原 (4.44)—(4.45) 明写的 S 通道上，校准可精确完成：设 raw CG 前的任意联合质量系数为 f(d,a)，归一化 CG 前的 multiplicity 系数为 g(d,a)=f(d,a)/√w(d)。原 (4.45) 给

\[
\|\Psi_f\|_{\rm FP}^2=\frac1{\pi^3}\int dd\,da\,|f(d,a)|^2
=\int dd\,\rho_{\rm Pl}(d)\int da\,|g(d,a)|^2.
\tag{M1}
\]

任意联合系数先由有限可分离测试函数取得，再按该范数完备化。这与“保留左 observer 的连续标签 d，物理内积取同标签 multiplicity 的加权直积分”完全吻合。它不是在一个既有运动学 Hilbert 空间里取单个测度为零的 trivial fiber。这里也没有仅凭 S 通道算例就宣称所有 N 的 FP 酉识别已完成。

## 2. 全 N 酉桥假设下的条件下降命题

固定 ν>0、Gaussian 窗 σ>0 及左右 observer 下界 s*>0。以受右质量截断的运动学空间

\[
W_N=H_{R,\ge s_*}\otimes\operatorname{Sym}^N H_\nu
\]

作 proper SO⁺(2,1) 表示分解。其总 principal 部分若统一写成

\[
W_N^{P}\simeq\int_0^\infty{}^{\oplus}dd\,igl(H_d\otimes\mathcal M_{N,d}\bigr),
\]

必须把所有内部 principal/discrete 表示与多重度纳入 M_{N,d}。总 discrete 表示不与左 observer 的 principal 表示配对；这不意味着可以删掉产生总 principal 的内部 discrete 路径。

**条件命题。** 假定对每个 N 已证明：

1. 上述分解是以与 dd 等价的 principal 测度给出的酉分解，归一化CG及内部测度均明确。
2. 与原 (4.42) FP 内积一致的酉映射存在：
   U_N:H_{N,Q}^{phys}→∫_{d≥s*}^{⊕}ρ_Pl(d)dd M_{N,d}。
3. 在一个完整物理稠密测试域上，原 (5.4) 的创建/湮灭模式积分，与下面的运动学等变创建映射及其伴随经 U_N 传送后的矩阵元一致；相关积分、CG求和及规约次序已获控制。

由已重建的顶角 bubble 与Gaussian包络，Tσ:H_R→H_R⊗Hν有界且等变，||Tσ||≤√𝔐。第N层创建映射为

\[
C_N=\sqrt{N+1}\,(1_R\otimes P_{\mathrm{Sym}^{N+1}})(T_\sigma\otimes1),\qquad
\|C_N\|\le\sqrt{(N+1)\mathfrak M}.
\tag{M2}
\]

适当的张量重排在式中省略。等变性使它分解成同一总标签 d 上的 multiplicity 映射 C_N(d)，其范数几乎处处满足同一界。因ρ_Pl(d)>0于d>0，换到物理测度不会增加新的非零测度例外集。故

\[
\|C_N^{phys}\|\le\sqrt{(N+1)\mathfrak M},\qquad
\|(C_N^{phys})^\dagger\|\le\sqrt{(N+1)\mathfrak M}.
\tag{M3}
\]

创建与湮灭在固定N态上落到正交的N±1扇区，因此候选对称场满足 ||Φσψ_N||²≤(2N+1)𝔐||ψ_N||²。若三项假设成立，有限粒子数代数直和便提供完整稠密对称域候选，而不只是真空顶角通道。此处没有证明全N酉桥，故不据此宣称原场已闭合或本质自伴。

## 3. 下一步最小任务：N=2 的完整树与 FP 范数

先做一个完整可检查的两粒子版本，再考虑归纳：

- 分解 H_R⊗Hν⊗Hν 时保留中间 principal 的 S/A，以及 D⁺/D⁻ 路径；用附录 E.1/E.2 的正交、完备系数统一归一化，明确每条连续谱测度及离散求和。
- 在有限通道和紧谱测试系数上，直接证明原 (4.42) 的 FP 二次型等于 ∫ρ_Pl(d)||g_d||²dd；不把形式 (4.33) 的发散量抵消当作这一步的全部证明。
- Bose 对称化作为与群作用对易的正交投影处理。它会混合耦合树，不能把投影后的树直接当成仍互相正交的基；应在 multiplicity 空间保留对应投影。
- 核对原 (5.4) 的局部模式创建与 (M2) 的 fiber 算符在这些测试态上的作用相同，再扩到完备化。

原 (4.50) 给的是抽象 L² 直和描述，把内部连续标签及多重度压入 j，并未在公式中展示这座具体酉桥。不能跳过归一化、测度与FP等距识别。

先完成 proper SO⁺(2,1) 空间，再施加与显式模式约定一致的空间反射约束；D⁺/D⁻会被反射交换。时间反演的反线性实结构另行处理。原脚注33与 E.16 的反射文字张力已在[解析审查](stage3-vertex-analytic-audit.md)单列，不用该脚注把所有A通道当作反射偶态。