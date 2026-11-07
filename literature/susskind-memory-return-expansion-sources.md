# 有限返回展开的来源与使用范围

2026-11-06，先读原文再推导。

- Montoya-Castillo、Reichman：[1603.01903v1](https://arxiv.org/html/1603.01903)，实际核 §III (11)–(26)，包含正/反 Q 的自洽辅助核与导数改写。消去投影传播不是原创方法；该文将路线归于 Shi–Geva。这里对有限 Hilbert-space 分块重新推导符号，不机械抄 Liouville 投影的最后一项。
- Widder、Zimmer、Schilling：[2503.20457](https://arxiv.org/html/2503.20457)，实际核引言和 §2 开头的有限秩/秩一及 Volterra 背景。我们的有限矩阵没有无界生成元定义域问题，不声称直接调用其所有定理。
- 本项目 `susskind-bilocal-channel` 的 V/K/Q/R/G 精确字典，`susskind-memory-force` 的普通传播接口，`susskind-memory-charge-sum` 的实际 pinching 场定义继续使用；不改变旧系数谱。

本轮自己写出的有限迭代、Schatten 余项、单个力 square functions 的组合计数和 Gaussian 分部积分，都是在明确模型中给出可检查推导。未作优先权断言。此处可核的新进展是：将现有方法落实到共同耦合的集体尺度，并证明只需两次返回即可使绝对协方差误差消失。完整有限记忆极限与体字典尚待求。
