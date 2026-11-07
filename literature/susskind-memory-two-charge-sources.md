# 双电荷记忆公式的前驱与使用范围

2026-11-06实际核读 [Kota–Vyas 1504.00759](https://arxiv.org/pdf/1504.00759)：

- §3，PDF页14–15，(34)/(35)及相邻文字：reduced product须保留multiplicityρ，乘积用U(N) Racah重耦合。这里没有把中性算符空间无重数误用为两个输入tensor product无重数。
- §7，PDF页45，(113)：相邻charge的乘积reduced amplitudes在同一个ρ上相干相乘，并有相位/Racah因子。它支持后续cross Gram应保留符号的组织方式。
- 本轮读取原文提取文本，未声称视觉逐页审查。原文dilute asymptotics不是当前半填充DS，不直接挪用其极限。

局部前驱：cubic-linear笔记和check_cubic_linear.py的charged(creation,annihilation)协方差特征值、primitive embedding；memory-quadratic-contraction五个精确单收缩；memory-charge-lift的primitive Rem/Add和半填充选择规则；memory-rank-blocks的正π、前三矩和尾；cubic-linear的input factorial tail。

本轮是以上工具对实际Q_Ecal的应用：通过两个charged component和HS adjoint推普适transfer；给一般有限N的两个显式滤波幅；由精确合并和旧尾界推H²部分L²近似。没有宣称首次发现charged representations、Wigner–Eckart或整个de Sitter字典。
