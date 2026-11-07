# F完整核：现有弦结果与本轮接口

2026-11-06，实际阅读优先。这个阶段没有新的核公式；任务是证明同J的F满足已有核的适用条件。

1. [Tall–Tomsovic 2604.14522v2](https://arxiv.org/html/2604.14522v2)：本轮实读§2.1 (2.17)–(2.21)的固定阶多弦、§4开头关于独立probe的使用、§4.1 (4.1)–(4.18)的generic t-body crossing与strength density、§4.2的不同probe标签。实际使用的是独立embedded ensembles的DS图规则及Mehler核。不能把§4开头O=αH+R的独立化叙述直接套给当前rank依赖滤波：本轮给出HF covariance仍非零，另证总变差桥接。
2. [Berkooz–Narovlansky–Raj 2006.13983v2](https://arxiv.org/html/2006.13983v2)：本轮复核§3.2 generic p_j index-set lengths、正常序trace/交叠和§8.1 (8.6)/(8.7)的不同neutral probe参数。它支持body长度可不同、各p−j固定偏移在DS有同样crossing的使用；不是只引用一个单probe二阶矩后推广。
3. [Kota 2208.11312v2](https://arxiv.org/html/2208.11312)：沿用前两轮已实读§IV (28)–(31)不可约covariance和维数权重，不重新推导嵌入系综的基础谱分解。当前精确shift ratio由这些binomial权重直接相除。
4. [Mukhametzhanov 2303.03474v3](https://arxiv.org/abs/2303.03474v3)：本轮核版本/标题，具体小λ公式沿用前轮已读§4及本项目任意d的antisymmetric后代证明。当前d1，−logq=4λ，Jsrc=2Js。没有声称本轮重新阅读该PDF全部内容。

本轮数学接口：rank滤波协方差以p/p−1/p−2系综的有限有符号组合逼近；rank权重L1与Gaussian/Hölder图界保证任意固定图可替换；逐边系数和消去HF配对。核、矩确定性方法、Charlier正交多项式、Hölder/Gaussian矩界均为已有工具。不作全领域优先权断言。

当λ>.5时FF组合含负系数；只在diagram multilinearity中使用，未从负权构造概率系综。核心极限是固定word degree和固定正λ，之后再λ→0，不包含N依赖的晚时间或scrambling窗口的统一界。
