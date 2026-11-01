# 热截断与矩阵元来源边界

2026-11-01，延续[上一阶段具体字典](stage4-paired-kernel-sources.md)。本轮未声称找到新的完整配对五边形证明。

- [Kim–Scarinci 2112.13329v2](https://arxiv.org/abs/2112.13329v2)：沿用已核§6.3产品/极点与Prop6.15。hbar1产品相位和已知平直对照来自该定义，不算新函数。
- [Garoufalidis–Kashaev 1706.08132v2](https://arxiv.org/abs/1706.08132v2)：沿用Appendix B已有有界Fourier框架。本轮的热截断不自动保持五边形。
- [NIST DLMF §7.8](https://dlmf.nist.gov/7.8)：在线核读Mills比定义及Gaussian尾界。本轮只用更弱的初等积分上界，并在稿中直接证明，不依赖特殊函数数值。
- [mpmath 1.3.0 Contexts](https://mpmath.readthedocs.io/en/1.3.0/contexts.html)：在线核读interval arithmetic的包含语义和实验性提醒，版本与本机依赖匹配。最终常数另用不调用超越函数库的精确Fraction区间独立认证。

已有Sobolev挖孔、热半群谱不等式、Gaussian积分和Taylor/交错级数余项属于标准方法；本轮工作是指定相位/测试下的明确常数与实际矩阵元复现，不作一般数学原创性声明。第一次打开mpmath原域名403，随后readthedocs固定1.3.0文档成功；不把失败访问列为已核依据。新在线两页只记链接与阅读范围，没有下载或分发第三方全文。

证明应按稿第2页产品梯度→第3页挖孔和加权积分→第4页有理常数检查，不能仅从audit passed推断证明成立。求积准确性与热截断界是两个不同问题，当前只有后者完成。
