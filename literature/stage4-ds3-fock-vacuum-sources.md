# Fock 与真空返回：前驱和引用范围

核读日期2026-10-31。正文见[LaTeX](../manuscript/ds3-fock-vacuum/main.tex)。以下为本轮已核读原文的用途，不宣称这些工具或标准 Wick 公式为原创。本轮未增加第三方全文缓存。

1. Derkachov–Spiridonov，[1711.07073v2](https://arxiv.org/html/1711.07073v2)，§3.1 的(22)、(29)、(30)及§4：完整 SL(2,C) CG 正交、完备和 reflection 代表。构造保留内部非零自旋，标量 synthesis 仅为一个正交子通道。具体 E/√w 归一仍由项目既有局部顶角确定，不能仅从协变性猜定。
2. [On the Positivity of de Sitter Group Averaging, 2610.24073v1](https://arxiv.org/html/2610.24073v1)，§§5.1–5.2、footnote29：tempered 张量分解、奇数 d 的主系列范围，以及两观察者无物质时必须作质量包。此文此前已读，本轮核对其适用条件。我们的种子绝对收敛证明明确先积分左质量，不将固定标签两球函数当绝对可积。
3. Kolchmeyer–Liu，[Chaos and the Emergence of the Cosmological Horizon, 2411.08090v1](https://arxiv.org/html/2411.08090v1)，§§4–6：观察者局部场、完整 Wick 路径的前驱。其 dS₂ 构造不直接证明这里 dS₃ 的全自旋、域或测度。两创建配对缺少真空返回这一点必须在具体模型中补齐。
4. Giulini–Marolf，[A Uniqueness Theorem for Constraint Quantization, gr-qc/9902045v2](https://arxiv.org/html/gr-qc/9902045v2)，§II：RAQ 测试域、null 商和伴随要求。这里只用明确受调节有界块及共同群截取证明弱下降，不声称直接满足其所有唯一性假设。
5. Nelson，*Analytic vectors*, Ann. Math.70 (1959),572–615。原 DOI 页面本轮未成功获取；判据用已有 Neeb 缓存核查：*On Analytic Vectors for Unitary Representations of Infinite Dimensional Lie Groups*, AIF61(5)(2011),1839–1874，Theorem4.8。缓存为 sources/neeb-2011-aif2660.pdf，本机 pdftotext 后核读：对称算符拥有稠密解析向量即本质自伴。此处用阶乘范数界验证假设，不把有限矩阵检查当替代证明。

依赖项目旧推导：[局部顶角](../notes/stage4-ds3-local-vertex.md)、[群配对](../notes/stage4-ds3-group-pairing.md)、[两物质](../notes/stage4-ds3-two-matter.md)、[Bose](../notes/stage4-ds3-bose.md)、[任意紧转移及全谱尾](../notes/stage4-ds3-compact-global.md)。本轮四时间极限在旧核紧 C² 条件及 Hilbert 尾上延伸；Gaussian 乘子有界但不允许移除 σ。

原创性尚未作穷尽查新；本稿定位为将已有工具连接到指定局部模型的可检查研究记录。有限 G 与微观字典需继续回读既有 TQFT/反冲资料。
