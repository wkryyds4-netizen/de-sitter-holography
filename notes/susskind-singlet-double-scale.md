# 首个导数 singlet 的双标度谱矩

2026-11-05，本地研究阶段。沿 Susskind 主线，研究对象仍是同一样本的 A=Σ[H,cᵢ]†[H,cᵢ]，现在从有限矩阵推进到一般 N,p 的系数公式和固定 λ 的解析极限。

[七页PDF](../output/pdf/de-sitter-susskind-singlet-double-scale.pdf)；[LaTeX](../manuscript/susskind-singlet-double-scale/main.tex)；[文献记录](../literature/susskind-singlet-double-scale-sources.md)；[验证](../results/susskind-singlet-double-scale/validation.json)。

## 得到的物理量

标准 CAR，H=H_s，完整 Hermitian complex Gaussian 系综，标准 CAR 系数方差 w=4ᵖ(J_s²/λ)(p!/Nᵖ)²。全态为β=μ=0，tr=2⁻ᴺTr，λ=p²/N。先 N,p→∞ 且 λ>0 固定。

取 O=A−E A，其中 E A=f(Q) 是确定的电荷函数。扣除它不改变[H,A]。记 q=e⁻⁴λ，本轮得到

\[
V_0^{\rm DS}=\mathbb E\operatorname{tr}O^2
=4\mathcal J_s^4(\lambda^{-1}+1+q),
\]

\[
M_2^{\rm DS}=\mathbb E\operatorname{tr}[H,O]^\dagger[H,O]
=\frac{8\mathcal J_s^6}{\lambda^2}
\left[1-q^2+4\lambda(q-q^2)\right].
\]

两者在固定λ>0都有限且正。其比值是有限系统归一关联函数初始曲率的双标度极限：

\[
\Omega_{\rm DS}^2=\frac{2\mathcal J_s^2}{\lambda}
\frac{1-q^2+4\lambda(q-q^2)}{1+\lambda(1+q)}
\xrightarrow{\lambda\to0^+}16\mathcal J_s^2.
\]

这把候选 singlet 在弦尺度上的初始演化变成了可检查的定量结果。下一步是求频谱形状并接体传播核。O仍可能含样本相关守恒部分；此比值不是粒子质量或晚时衰减率。这里也未交换时间极限与双标度极限。

## 为什么能算一般 N,p

正常序算符的body rank算符为 D(X)=QX−Σc†Xc。H²中每个左右块收缩降低一个body rank，因此给这些收缩乘(1+h)后，对h求导恰好给 A=(2p−D)H²。这一步使用普通费米 Wick 正常序。

一个六次Hamiltonian字在Gaussian耦合配对π之后，只需六个独立的p元创建指标集。自由无限温迹是一个6p阶Wick行列式；按模式重新分组后，变成最多6阶的局部行列式乘积。列块配对产生全局符号(−1)³ᵖ，不能漏掉。

对四种 K†K 字（AHHA、AHAH、HAHA、HAAH）分别插入两组标记，得到

\[
M_2/w^3=(-1)^{3p}[z_1^p\cdots z_6^p hg]
\sum_{\pi\in\mathcal P_\times,\,\mathcal M}
\epsilon_\mathcal M\det(I+ZT_\mathcal M(h,g)P_\pi)^N.
\]

P×是两组三个位置之间的6种配对；可只留这些配对的原因是上一阶段已证明K的一次Gaussian分量为零。普通系综中心化A的Gram则剩两种跨配对。

把T乘2后可全程用整数系数，最终抽取项统一除2⁶ᵖ。因而 M₂/w³=2⁻⁶ᵖΣⱼ binom(N,j)cⱼ(p)，j≤6p。每个p只算一次，之后任意N直接代入；不再构造2ᴺ维矩阵。

特别地，p=2的精确结果为

\[
\frac{M_2(N,2)}{w^3}=
\frac{N^2(N-3)(N-2)(N-1)(N+1)}{512}
(N^4+4N^3-53N^2+120N-16).
\]

## 双标度推导

六个随机p元集的两两交叠数趋于独立Poisson(λ)。三重交叠概率≤20p³/N²→0。单个交叠的超几何阶乘矩给指数矩界，进而控制带符号/复权的总和。界在标记变量h,g的复邻域一致，因此可用Cauchy公式交换标记求导与极限，而不只凭交叠数依分布收敛。

无标记的指数逐个恢复BNR已有q弦权重。24个带标记项按交叉数0、1、2、3合并，对λ和λ²的系数依次为(8,0)、(0,32)、(−8,−32)、(0,0)，给出上面的M₂。四位置的两项同法给V₀。

## 验证和计算边界

- 720个带标记字、2880个有理数对照与直接CAR矩阵完全一致，含非零负项；保存了全部指标集。
- p=1,2,3的完整15配对与只留6种跨配对的整数数组完全相同。p=1全零；p=1至5的N≤p+1控制均零。
- N=5,p=2另用上一阶段独立全协方差张量方法得到2700；N=4,p=2复现120。这些不是Monte Carlo估计。
- 保存p=1至5的精确数组与因式分解。在λ=1/4，M₂/极限值从p=2的.691904升至p=5的.878486；λ=1时p=5仍只有.628144。有限p偏差不可忽略，极限结论来自证明。
- 初版交叠指数审计遗漏零局部子式对应的−1，断言失败。修正后24个标记图与15个已有六阶弦权重通过；有限系数算法和数据未受影响。原代码和失败说明保留。

所有求和为完整有限代数，没有谱截断或随机误差。Python3.14.5、NumPy2.4.6、SymPy1.14.0。七页PDF已双遍编译、全页1200px检查；初版长文件名溢出已修复，日志保留。仅本地保存，无新Git提交/推送或ZIP。

## 下一步

先分辨O中剩余的零频成分，再求第四频率矩M₄=E tr((ad_H)²O)†((ad_H)²O)。它是八阶耦合平均，四次Gaussian多项式可能有二次分量，不能照搬本轮“只留跨配对”。优先查已有带标记弦／复合探针结果，随后建立相应精确生成式或直接能量核；目标仍是完整singlet谱及Rindler体传播映射。
