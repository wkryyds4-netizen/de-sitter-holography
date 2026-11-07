# 固定扇区记忆返回：来源与本轮实际使用范围

2026-11-08。先回读已有推导与来源，再采用精确反演路线。

- [Berkooz–Narovlansky–Raj, 2006.13983v2](https://arxiv.org/html/2006.13983v2)：complex SYK、固定电荷背景；本轮重开原文。其独立 probe 公式不替代实际 force 的共同J收缩。本轮强返回证明不调用独立probe四点。
- [Bożejko–Kümmerer–Speicher, funct-an/9604010](https://arxiv.org/abs/funct-an/9604010)：q-Gaussian条件期望背景；本轮重开摘要和来源入口。旧return-background已核的Mehler接口保持归属文献；本轮通过有限参数Jacobi函数完成守恒比较，不新增依赖其定理的强桥步骤。
- [Montoya-Castillo–Reichman, 1603.01903v1](https://arxiv.org/pdf/1603.01903)：重开原PDF，核§III式(11)–(13)，以及引言对Shi–Geva辅助核方法的归属。Dyson展开去掉投影传播是已有方法。本轮有限矩阵Volterra符号和有序卷积另行写出并检查。猜测的v2 HTML链接404，随后查摘要确认仅v1，原PDF读取成功；未将失败链接作为核读依据。

本地直接证明输入：

1. [operator draft](../notes/susskind-bilocal-fourth-tail-operator-draft.md) 式(2)、(21)：有限degree OU/Chebyshev multiplier与body/Gaussian product rule。
2. [sector-gradient-tail](../manuscript/susskind-sector-gradient-tail/main.tex) 主定理、(O)、(T)：中央O(p)窗口的实际ordinary factorial包络。归一C_m因子在loss后作用，转换回string时间只改变常数。
3. [sector-neutral-transfer](../manuscript/susskind-sector-neutral-transfer/main.tex) §4：精确Jacobi指数展开、实际T_j尾和有限参数theta。这里只需要未pinching的exact decomposition及actual bounds；不以极限谱回归替代有限N误差。
4. [memory-return-expansion](../manuscript/susskind-memory-return-expansion/main.tex) §§2–3：有限Volterra等式、两侧force covariance及逐sector Gaussian Gram高矩。将固定r的Holder阶写成6r，直接控制E1，不需要binomial charge汇总。
5. [memory-return-background](../manuscript/susskind-memory-return-background/main.tex) §§4–5：DF与R关系、统一源维数矩形范数、有限alpha守恒比较结构和小时间校准。旧双tag强桥以fulltrace写出，本轮不直接限制它。
6. [polynomial-density](../manuscript/susskind-slater-polynomial-density/main.tex)：已控制固定Gaussian-degree制备的有限N轨道密度矩。只用于误差在固定制备core上的传递。

本轮新增的是这些接口的组合：exact inverse→force defect紧时间O(1)→finite-N conserved proxy→sector all-r return suppression。标准逆式消去、CP收缩、Holder、Volterra与Mehler结构不作原创性声明。仍未处理X0的两分量强极限及完整selected bilocal。
