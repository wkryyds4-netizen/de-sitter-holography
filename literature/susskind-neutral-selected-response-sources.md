# 中性选态响应：前驱与公式来源

2026-11-08，先读前驱再构造新计算。

- Bożejko–Kümmerer–Speicher，[q-Gaussian processes: non-commutative and classical aspects](https://arxiv.org/html/funct-an/9604010)。Theorems 2.11、2.13给第二量子化和条件期望；4.4、4.5、4.6给有序Markov过程的经典版本及q-Mehler转移。此处三段能量协方差(q,q,q²)属于该已有结构，不作为新Markov定理。转移参数与形变参数在本算例恰好同为q，正文区分了二者。
- Effros–Popa，[Feynman diagrams and Wick products associated with q-Fock space](https://pmc.ncbi.nlm.nih.gov/articles/PMC166362/)，Corollary 2.1、Theorem 3.2。已有交叉权q-Wick规则。当前四U中，扣去最前两U互配后，两剩余拓扑给相同三块权和系数1+q。
- 原始q-Mehler乘积继续来自BINT [1811.02584v2 §3.2](https://arxiv.org/html/1811.02584v2#S3.SS2)。数值角度实现沿用已核公式，程序独立写入新pure helper；不重复旧driver。

本地实际接口：slater-neutral-slope的rank归一；slater-orbit-density的单selected实际exact-pinching桥；slater-selected-state的完整rho；sector-bilocal的有限协变块与actual重建；weighted-mehler-limit的全谱主控、带权展开和固定整数d结论。本轮使用其d=1，新增U的实际解释来自已有slope桥，并非仅凭整数参数猜测。

本轮工作是该模型选态密度的三能量精确响应、端点安全的联合尾估计、紧时间中性记忆与odd能量—斜率响应。尚未穷尽相似四场表达的文献，不作普适公式优先权主张。空间/质量生成元与finite-G不是这些q-Gaussian结果的直接推论。
