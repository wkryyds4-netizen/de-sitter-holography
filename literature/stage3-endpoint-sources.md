# Euler端点方法的来源与归属

2026-10-30。承接[上一轮来源图](stage3-crossed-large-mass-sources.md)。

- [KL 2411.08090v1](https://arxiv.org/html/2411.08090v1)：本轮重新读取B.21—B.24、B.29—B.33、B.41—B.42。[缓存](../sources/kolchmeyer-liu-2411.08090v1.html)不变。超几何参数为原式；新工作是对其端点保留两支并控制单因子余项。
- [DLMF 16.5.2](https://dlmf.nist.gov/16.5#E2)：Euler型降阶积分，需Re b₀>Re a₀>0；p=q+1情况亦需注意自变量条件。这里的z=1通过额外绝对可积证明取得，不只引用开域公式。
- [DLMF 15.6.1](https://dlmf.nist.gov/15.6#E1)：Gauss Euler积分。DLMF使用Olver正则化函数；转换成普通₂F₁须乘Γ(g)。
- [DLMF 15.8.4](https://dlmf.nist.gov/15.8#E4)：z与1−z的连接公式，两个端点支均须保留。已有[15.8缓存](../sources/dlmf-15.8.html)。另外在线核对15.10.21的等价普通函数写法。
- [JKMS 2209.02131v3](https://arxiv.org/html/2209.02131v3)：本轮额外读取A.37—A.39的Barnes表示与极点分隔规则。它属于JT Wilson函数，仍不能直接替换dS₂核；本轮选择从dS₂的原₄F₃使用Euler变换。

16.5与15.6在线公式已读，本地Invoke-WebRequest均因DNS解析失败，未形成缓存。不要把来源链接或手抄公式说成原文归档；manifest记录失败。15.10仅在线核对，未尝试下载。

Fourier端点论证在[第10节](../manuscript/endpoint.tex)自行给出：平滑cutoff、Abel边界值、dyadic分部积分、复δ圆周估计与最大模原理。它应用标准方法，没有做足以支持原创性声明的完整重叠检索。

本轮结论不否认KL的完整四支渐近；单项遗漏的非衰减支可能在组装时消去，数值确实提示抵消。必须先做代数合并及可去极点处统一余项，才可讨论完整分布结论。
