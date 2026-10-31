# 实 principal 探针、轮廓与谱矩：来源边界

2026-10-31。配套[阶段稿](../manuscript/principal-eikonal/main.tex)。这是对已有 eikonal 模型的推导与计算，不作首创性声明。

- [CSTY v1](https://arxiv.org/html/2607.14042v1)：核读(2.25)、§2.4.1原动量积分(2.32)–(2.36)、§2.4.2矩展开(2.37)–(2.40)与附录G。Bessel波函数、固定平均质量表示、shock与两个sector相加减断开项都是既有输入；本轮保留真实principal标签，推四时间Gaussian平均。负轴必须保留侧向轮廓。原文脚注12明确要求相反动量旋转使p+p−保持实数，本轮直接用±π/4旋转复算，不能将负耦合发散的正实Laplace积分作为原公式。
- [DLMF 10.32.9](https://dlmf.nist.gov/10.32.E9)：K的cosh积分，Re z>0；实principal给cos(νr)，权重可能变号。[10.43.19](https://dlmf.nist.gov/10.43.E19)：Mellin变换的Gamma乘积，条件|Re(order)|<Re(exponent)；本轮exponent=(n+2)/2−iω、order=iν满足。Gaussian谱正性来自标准自相关/Plancherel，不能包装成全量子引力正性构造。
- [MNX v1](https://arxiv.org/html/2607.13137v1)：已核§3.2、§4、§5.1–5.4及§5.5结论、§6，相关中间公式以本机提取和原网页互核，未复现全部附录。§5.23–5.28已包含一般Bessel矩及质量规约下展开；§5.4横向微分作用量展示L=1耦合。graviton质量规约的双缩放与不带该规约的tree leading G不同，也不等于CSTY观察者反冲。不能借此替换我们的轮廓或假定external dressing已完成。
- [Bonifacio–Hinterbichler–Law–Lochab v1 PDF](https://arxiv.org/pdf/2610.28719v1)：HTML有效作用量处截断，实际改读PDF第8–15页，覆盖§§4–6及附录A；此前§§2–3已读，附录B未完整读。D≥4、测试观察者范围；Z_edge是一圈W_edge的倒数。c3=−c2复现指定二次谱；第一性bulk-edge分裂、量子化、边界条件/角点歧义及有应力张量的物理观察者耦合仍未给出。本轮不直接外推到dS3 massive observer。

本轮最接近的前驱就是CSTY的未涂抹矩展开与MNX的Bessel矩。新增推导的定位：真实四时间Gaussian自相关给正谱平均，完整lateral核以四rapidity积分保留相位；数值仅在明确模型内。未完成足以支持优先权声明的系统查新。原文缓存只在本机，版本/哈希见本阶段source-record.json；没有修改旧sources/manifest或分发第三方全文。
