# 中央电荷窗口的文献核对

2026-11-06，实际在线读取以下原文相关节。

- Tall–Tomsovic，https://arxiv.org/html/2604.14522v2 ，§2.1式(2.5)–(2.13)、§6式(6.3)–(6.7)。使用固定sector归一化、共同coupling variance与reduced moments的区别。全文公式有转换/排版不一致的地方，本轮中央修正从项目精确binomial公式直接展开，不照搬§6鞍点表达式。
- Förste–Kruse–Natu，https://arxiv.org/html/2512.07715v2 ，§3.1–3.2、式(3.1)–(3.5)及charge直和讨论。复用canonical/grand canonical的sector结构作为背景；未据此推定marked memory算符的中央窗一致性，也未复用未经规范对齐的能量系数。
- Berkooz–Narovlansky–Raj，https://arxiv.org/abs/2006.13983 ，本轮检索核其固定charge主题；此前已有原文笔记，本轮不记为重新全文阅读。
- 项目 `results/susskind-memory-two-charge/symbolic-filter.json` 和 `results/susskind-memory-filter-limit/filter-certificate.json` 提供精确raw滤波器，本轮直接读取并哈希。
- `code/memory_sixth_moment.py` 提供Racah矩的递推框架；当前新代码显式改变charge scaling，不把bounded shift证书直接推广。

本轮增加的是中央窗一致六阶矩证书、共同方差的p倍有限漂移及raw滤波器的有限守恒修正；这些是模型内的计算结果，未作原创优先权检索完成的声明。
