# Gaussian 与硬边界速率：前驱及证据范围

2026-10-30。重读[KL 2411.08090v1](https://arxiv.org/html/2411.08090v1)的B.50—B.57和之后的轮廓说明，主源沿用[已缓存原文](../sources/kolchmeyer-liu-2411.08090v1.html)。原文已有奇数Gamma极点抵消及指数2，不将F协议的同一指数称为原创。

本轮添加的具体步骤是：完整外腿与V₋相乘后的横向分离、固定proper Gaussian的条带支配及正系数，以及指数外部密度/逐次移动硬阈值下的角点尾推导。每次阈值与原文只在外部施加microcanonical窗不是同一制备。

Gamma的极点/残数从[DLMF 5.2(i)](https://dlmf.nist.gov/5.2)在线核读，条带渐近复核[DLMF 5.11](https://dlmf.nist.gov/5.11)，后者已有[缓存](../sources/dlmf-5.11.html)。前者仅在线读取，未尝试新增缓存。极点与Fourier跳跃的标准数学结论不构成新发现；新的模型系数与支配推导由本项目承担核查责任。

定向搜索了“de Sitter quantum observer OTOC mass threshold power law finite preparation”和“OTOC hard energy cutoff algebraic late time tails principal value observer”。结果含既有观察者代数工作与其他机制的代数尾，未据搜索摘要认定重叠或新颖性，也未完成系统查新。这里不作首创性声明。

Q的结论目前限定π₀(r)=β exp(−βr)。其|T|⁻³来自真实积分区域随χ变化的角点，在减掉D_edge后仍存在；未添加有限G相互作用，也没有证明所有制备都如此。F实际曲线已有经验求积；Q只有局部角点诊断，完整Q曲线尚待复算。

见[推导](../notes/stage3-preparation-asymptotics.md)、[LaTeX](../manuscript/preparation.tex)及[数值输出](../results/stage3/preparation/diagnostics.json)。目前是作者自查证明草稿，没有独立审稿认证。
