# 交叉字典：原文、规范与适用范围

2026-10-31。接续两物质记录；本轮重新读取本机版本化全文，并在线复核 LPRS/KKS 书目。没有新增或改写原文缓存。

1. Liu, Perlmutter, Rosenhaus, Simmons-Duffin, *d-dimensional SYK, AdS Loops, and 6j Symbols*, [1808.00612v3](https://arxiv.org/abs/1808.00612v3)。本机 [全文](../sources/liu-et-al-1808.00612v3.html)：(3.4)–(3.5) 定义 s 通道共轭与 1↔3 后 t 通道的 quotient 配对；(3.14)–(3.17) 给两 shadow 项、block 与 K；(3.20) 显示 1/4 规范固定因子；(3.30)、(3.35)–(3.38) 固定 t 外部顺序 (3,2,1,4) 和解析 6j。沿用旧 global_racah.py，不将这些已知公式作为新结果。
2. Karateev, Kravchuk, Simmons-Duffin, *Harmonic Analysis and Mean Field Theory*, [1809.05111v1](https://arxiv.org/abs/1809.05111v1)。本机 [全文](../sources/karateev-kravchuk-simmons-duffin-1809.05111v1.html)：§2.3 (11)–(15) 固定 quotient pairing、Faddeev–Popov 因子 2^{-d} 及 vol SO(1)=1；§2.7 (36)、(39)、(49) 给 shadow 积分、delta 范数及两 block 展开；(74) 给 μ=1/[2^d N vol SO(d)]。d=2 标量 shadow-square N=π²/b²，因此 μ=b²/(8π³)，原始 CPW 的正标签 db 范数为 π⁴/b²。这里的 μ_KKS 不是观察者时空谱密度 μ_obs=b tanh(πb)/π。
3. Derkachov, Spiridonov, *The 6j-symbols for the SL(2,C) group*, [1711.07073v2](https://arxiv.org/abs/1711.07073v2)。本机 [全文](../sources/derkachov-spiridonov-1711.07073v2.html)：§3 的完整 CG 正交/反射与 §4 (32)–(33) 的两耦合次序变换。本项目正标签代表与全部整数自旋属于完整谱。DS 到物理有向权重的约定沿用旧 global-racah 笔记 J=−m/2、ν=−σ；这里仅算 J=0 元素，不能据此略掉全变换的自旋和。

本轮新增的是将上述已知 harmonic analysis 与本项目局部 BD/KG 模型的 000=1 通道、物理逆范数及创建振幅对接，并设计独立 Euclidean 面积积分。规范常数由三点和四点范数推得，没有数值拟合。解析尾界是本轮自查推导；数值不是区间认证。全自旋完备性引用上述理论，不声称用几个标量样本检验了完整酉性。
