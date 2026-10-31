# 群平均文献更新与本项目的测度字典

2026-10-29。新增主来源为 Donald Marolf、Xuyang Yu，[On the Positivity of de Sitter Group Averaging，2610.24073v1](https://arxiv.org/html/2610.24073v1)，提交于2026-10-16。本轮精读§3.1—3.3，并核对引言给出的适用范围；未把其余维度及轻场部分当作已逐式复现。[本地原文](../sources/marolf-yu-2610.24073v1.html)的SHA-256为 `4b55e4b58532dc8d7483f00dd7193f362bb326e958f9137bb24387567c0eb3cc`。

## 直接相关的已有结果

原§3.2—3.3以主系列的分布正交和完整CG分解建立dS₂主/离散系列的群平均正性。S、A反射奇偶在§3.3明确区分；这是此前KL脚注33张力的独立文献依据。矩阵元收敛的测试域及表示积分换序须保留，正性不意味着群平均在运动学范数下有界。以下只使用其(3.19)、(3.24)、(3.35)的归一化及§3.3的分解输入。

| 名称 | 定义 | 与本项目的对应 |
|---|---|---|
| MY群Schur系数 | K_P(d)=4π²/[d tanh(πd)] | 1/ρ_Pl(d) |
| MY的ρ_P与KL的PP连续CG权 | w(d)=πd tanh(πd)/4 | π³ρ_Pl(d) |
| KL单observer运动学谱测度 | dd，基矢δ归一 | 不是ρ_Pl(d)dd |
| normalized PP合成系数 | Ĉ=√w C | 连续输出采用裸dd |

MY采用KL附录E的基底；本项目的物理模式采用KL正文相位。涉及磁指标的公式必须先按KL(3.11)、E.24—E.25转换，不能直接删除MY(3.25)—(3.26)中的相位。K_P、w与ρ_Pl的正实关系不依赖这个基底相位。

## 这改变了哪一步研究

对全主系列的两observer加物质模型，群平均正性与CG约化已有明确前驱，不再把它们列为本项目待首创的命题。我们需要完成的是这些结果与**本项目具体的**KL位置FP内积、质量投影和局部模式创建之间的识别。

[谱系数推导](../notes/stage3-ga-spectral-bridge.md)把seed收缩b和不变幅度g区分：g=K_P b，故∫K_P|b|²=∫ρ_Pl|g|²。不能把群平均seed的权K_P直接当成不变态的物理权。这个倒数关系是现成Schur公式的应用，没有新的正性定理之意。

[两粒子通道表](stage3-two-particle-cg-map.md)另核principal-discrete顶角的权v_Δ(d)=d tanh(πd)/(Δ−1/2)。它与w不同；若位置FP桥成立，裸外PD顶角系数的权为ρ_Pl/v_Δ=(Δ−1/2)/(4π²)，而非PP情形的1/π³。

当前推进条件仍是[FP内积稿](../notes/stage3-fp-two-particle-inner-product.md)的覆盖性及[局部创建下降](../notes/stage3-creation-descent-audit.md)。新文献不自动证明这两项。此前已冻结的原文和归一化疑点记录保持不变；本稿补充后续依据。
