# 全电荷求和的文献与本地输入

2026-11-06。本轮在线打开并核读 Berkooz–Narovlansky–Raj，*Complex Sachdev–Ye–Kitaev model in the double scaling limit*，arXiv:2006.13983v2，§4.3.1及§5开头、式(4.16)、(5.1)–(5.4)：https://arxiv.org/html/2006.13983v2 。原文已有grand-canonical/fixed-charge分解、缩放chemical potential以及partition function factorization。我们的对象是已减守恒背景memory的双能量matrix；不把原有partition-function因子化称作本项目发现，也不直接用它替代marked字段一致界。

本轮使用共同coupling规范的精确C_m，完成平方直接得到更紧的全sector方差界。Chernoff界从binomial mgf现场推导；finite-dimensional pinching等于时间平均由谱展开直接证明。没有为这些标准工具主张原创性。

本地实际读取：

- `manuscript/susskind-memory-central-raw/main.tex`：真实C0/C1/C2以及(Y,L,U2)central谱，不改counterterm参数。
- `manuscript/susskind-memory-central-norm/main.tex`、`code/memory_central_norm_asymptotic.py`：Racah summand阶数、x²项、p²/p相消、factorial tails、whole-rank粗界。
- `code/memory_central_window.py`、`results/susskind-memory-central-window/central-moments.json`：固定长度六阶递推的rational denominators及QQ[A,B,z]多项式相消；新增长窗论证不只引用旧“compact z”措辞。
- `code/memory_central_trace.py`：occupation ratio、global upper bound与log remainder，重估x=O√logp后误差。
- `manuscript/susskind-memory-force/main.tex`、`manuscript/susskind-memory-fluctuations/main.tex`：早期平均memory的尾界和C0/C1定义。原平均核的尾界不能自动证明本轮更大的集体涨落tail。
- `code/cubic_linear_channels.py`及`code/cubic_linear_dynamics.py`：一般charge t_a、positive covariance eigenvalues与rank weights。
- `code/singlet_wick_moments.py`：Gaussian整数CAR基，用于新的全sector分母验证，不运行旧main覆盖结果。

新增内容为：全sector反项定义、统一polynomial norm界、增长窗余项审核及fulltrace UI/谱转移；由此把精确pinching减背景接到既有Taylor谱。有限memory covariance及体空间/质量字典尚未完成。
