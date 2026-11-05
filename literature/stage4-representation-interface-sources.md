# 谱表示接口：本轮文献核查

2026-11-05。局部原文阅读与适用性检查；不是穷尽查新，不宣称标准定理原创。

| 来源 | 实际读取位置 | 可直接使用的内容与限制 |
|---|---|---|
| [Kim–Scarinci 2112.13329v2](https://arxiv.org/abs/2112.13329v2) | 在线版本页；本机 PDF 提取的 §6.3 Def.6.3、Prop.6.6(SB1/SB4)，§6.5 Prop.6.15 及证明 | 正形变参数 t 的 Barnes 参数是 h=it；SB1 的“实 h 有理时零极点简并”不是 t=1 排除条件。Prop.6.15 假设完整 Weyl 关系，并在平直分支证明中直接使用 Stone–von Neumann。读取这些段落不等于独立验证全部证明。 |
| [Kim 2112.13962v2](https://arxiv.org/html/2112.13962v2) | 在线 §7 Props.7.3–7.4 及两命题之间的证明来源说明 | 明确引用 KS 的参数解析延拓；不作为另一份独立五边形证明。 |
| [Dereziński math-ph/0511030v2](https://arxiv.org/pdf/math-ph/0511030v2) | 在线 PDF 文本 §3.1 正规性定义，§3.2 Theorem 4，印刷第12页；相邻 Schrödinger 表示定义 | 有限维非退化辛空间的正规 CCR 表示等于 Schrödinger 表示张量重数空间。不要求不可约；不可约对应一维重数。这里的“正规”指 Weyl 单参数群强连续，不是正规算符。本轮直接引用该定理，不另造分类证明。截图调用失败，未声称原文页面图像核版或本机全文下载。 |

一般参数产品由 KS 的 compact 产品及共轭配对代入得到；强连续性用一致收敛的产品尾和支配收敛，有限个有界算符乘积用标准强拓扑连续性。两条推论依赖此前指定算符的正残差输入；本轮没有重新认证该数值证书。

## 勘误与版本检索范围

检索式：`"2112.13329" erratum correction pentagon`、`"2112.13962" pentagon correction`、`"2112.13329" "irrational"`。

在线 arXiv 页仍列 KS v2（2024-06-21）、Kim v2（2023-07-05）。另检查 [KS 出版社页面](https://link.springer.com/article/10.1007/s00220-024-05012-8) 和 [Kim 作者论文页](https://sites.google.com/site/hyunkyukimmath/papers)，页面文本没有匹配的 correction 条目。没有定位到相关勘误，不代表证明它不存在；本轮没有得到出版 PDF 与 arXiv 的逐字差异比较。

## 本轮避免的重复

此前 flip-domain 已处理形式/强对易、单因子域和相位奇点，weak-continuation 已处理移动极点弱积分；本轮不重做这些推导。此前 residual-certificate 已认证固定参数的连续残差，不重跑其求积。本轮新增的是其适用范围推论：正参数开邻域、正规 Weyl 表示的重数等价，以及稠密域不能改变最终有界词恒等式。

参见[本轮 LaTeX](../manuscript/representation-interface/main.tex)。定量邻域宽度、非正规或改谱完成、实际有限 G 物理态仍未确定。
