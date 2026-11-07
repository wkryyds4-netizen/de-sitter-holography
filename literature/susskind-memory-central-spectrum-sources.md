# 中央谱扩展的文献与既有证书

2026-11-06 本轮实际在线打开Berkooz–Narovlansky–Raj，*Complex Sachdev–Ye–Kitaev model in the double scaling limit*，arXiv:2006.13983v2，§3.3及式(3.11)–(3.15)：https://arxiv.org/html/2006.13983v2 。原文处理有向fermion word、保留交叉符号、一般index-set lengths、Poisson pair overlaps与高重叠消失；这些是既有输入，不能声称由本项目创立。CAR归一常数沿本项目已匹配规范。

回读项目主稿：

- `susskind-memory-descendants` 的half trace bridge、lower-chaos表，确认后者仍含p增长counterterms，不能用O(1/p)迹误差直接升级raw U2。
- `susskind-memory-oriented-kernel` 的实际双向commutator mark、local minors、overlap double-mark bound、128 reverse endpoint templates、Mehler matrix及spectral moment argument。
- `susskind-memory-mixed-spectrum` 的upper/lower charge lift与filtered F接口，留待中央rank-filter检查；未把half F bridge直接当中央结果。
- `susskind-memory-raw-spectrum` 的实际raw U2 restoration，明确下一步需要保留所有动态后代。
- 最近central-gram的实际有限N投影与strong core展开，作为本轮R2谱的算符输入。

本轮新增：central occupation ratio的显式log余项与uniform O(1/p) graph transfer；一般charge Casimir/force桥接中的N−2m修正；由既有kernels/certificates接central oriented measures与实际projected W2core谱。原有五份kernel/endpoint科学输出及其输入hash由 `code/memory_central_spectrum.py` 验证后复用，没有重复运行旧计算，也没有重新声称既有Mehler核的原创性。
