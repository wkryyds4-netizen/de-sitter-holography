# raw memory 谱的来源与复用范围

2026-11-06。来源先于新组装：本轮使用既有 Mehler 测度和项目已经保存的 raw contraction 表，没有把这些核重新命名为新结果。

- J. Tall、S. Tomsovic，https://arxiv.org/html/2604.14522v2 ，第4节、式(4.7)–(4.11)：strength density 与能量本征基两点函数的 Fourier 接口。本轮实际在线读取。这里只复用两能量核及频率加权的解释，不据此声称已构造 de Sitter 对偶。
- M. Berkooz、V. Narovlansky、H. Raj，https://arxiv.org/abs/2006.13983 ：complex SYK 的 chord/probe 背景。原文在前阶段核查；本轮沿用既有 probe 核，不声称重新完整阅读。
- `results/susskind-memory-descendants/symbolic-filters.json`：raw 五字的 W2 表，当前符号代码直接读取，并保存输入哈希。
- `manuscript/susskind-memory-quartic/main.tex`：Q4 投影系数及 raw renewal 分离。
- `manuscript/susskind-memory-oriented-kernel/main.tex`：输出缺额与双向 force、X 的两部分谱。
- `manuscript/susskind-memory-mixed-spectrum/main.tex`：Y/F/R/T/Q4 的完整两能量正交性。本轮的谱加法使用这一接口，不仅使用等时 chaos 正交。

本轮组装的内容：明确逐样本守恒反项，恢复动态后代，得到 U2 的谱及 Y/L 交叉谱，并取有序小 λ 极限。未作优先权判断；更高记忆时间阶与体质量解释仍需继续研究。
