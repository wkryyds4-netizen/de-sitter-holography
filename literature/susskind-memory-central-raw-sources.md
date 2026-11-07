# 中央 raw memory 的来源与复用范围

2026-11-06。本轮重新在线核读 Berkooz–Narovlansky–Raj，*Complex Sachdev–Ye–Kitaev model in the double scaling limit*，arXiv:2006.13983v2，§3.3，式(3.11)–(3.15)：https://arxiv.org/html/2006.13983v2 。原文保留有向 word 的 fermion signs，允许不同 index-set sizes，给 Poisson pair overlaps 及 triple overlaps 在极限中的消失。这提供旧弦规则输入；本项目的 central raw counterterm 不是该文现成结论。

本轮逐项回读并复用：

- `code/memory_descendant_filters.py` 与 `results/susskind-memory-descendants/symbolic-filters.json`：一般 m 的精确五 raw words 与 contraction filters。新中央展开来自这些有限 N 公式，不来自已取极限的 half 表。
- `code/cubic_linear_channels.py`、`code/cubic_linear_limit.py`、`manuscript/susskind-cubic-cross/main.tex`：t_a 精确式、rank recurrence、s_a/z_a、中心化和 cubic conserved polynomial。
- `code/cubic_covariance_bridge.py`、`code/singlet_filtered_bridge.py`：R_j 精确概率比、全图 TV 界、下 charge lift。
- `manuscript/susskind-memory-mixed-spectrum/main.tex`、`code/memory_mixed_lift.py`：一般 charge 的 upper commutator lift 与 k-slot normal mark lift。旧代码中有限 tensor 核对的参数和证明范围保留，未将旧 finite checks 说成本轮新增。
- `manuscript/susskind-memory-central-norm/main.tex`、`manuscript/susskind-memory-central-spectrum/main.tex`：中央 A 范数/尾、strong core、quantitative trace transfer、一般 Casimir 中 E_in 的 N−2m 项。
- `manuscript/susskind-memory-raw-spectrum/main.tex`、`results/susskind-memory-raw-spectrum/symbolic-spectrum.json`：raw identity、renewal counterterm、动态后代及正谱分解。

本轮新增范围：中央 low-chaos table 的 2λx ell B 修正；F 的中央 slope/factorial 控制与 signed lift 接口；由这些接口组装中央 raw U2。Mehler 核、旧 endpoint 相消、half raw kernel 均为既有输入，不计为新原创。完整文献重叠与 Susskind 体字典仍沿原项目主线推进。
