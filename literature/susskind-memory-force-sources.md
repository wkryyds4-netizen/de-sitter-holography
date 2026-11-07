# memory-force阶段来源核查

2026-11-06。先核已知方法，再处理本项目同J接口；不将现成公式改写成原创。

1. J. Bourgain, F. A. Grünbaum, L. Velázquez, J. Wilkening, *Quantum recurrence of a subspace and operator-valued Schur functions*, [1302.7286v1](https://arxiv.org/pdf/1302.7286)。已读PDF正文6–7页的首次返回振幅、非交换renewal与生成函数逆。本文用的是其标准P/Q字分解，应用于self-adjoint K的代数幂；不直接把离散监测概率解释搬到Hamiltonian时间。
2. M. Berkooz, V. Narovlansky, H. Raj, *Complex Sachdev-Ye-Kitaev model in the double scaling limit*, [2006.13983v2全文](https://arxiv.org/html/2006.13983v2)。已核§4交叠Poisson及(4.8)、§7 (7.1)–(7.3)独立Gaussian probe与大小∼√N、(7.9)–(7.15)charge规则。q-Hermite谱、Mehler probe kernel是既有结果。本轮不把size1 fermion代进终式后放大误差，而先对D_i=[H,c_i]进行显式tagged contraction，保留同J外部配对。
3. J. Maldacena, D. Stanford, *Comments on the Sachdev-Ye-Kitaev model*, [1604.07818v1全文](https://arxiv.org/html/1604.07818)。已核§2.4 (2.14)、(2.17)–(2.19)及(3.161)：large-interaction-order传播的1/q展开与双侧续接。这里固定时间log cosh只恢复其已知结构；interaction order为2p、能标2Js，本项目q=e^(−4λ)是另一参数。

本地前驱：[bilocal-channel](../notes/susskind-bilocal-channel.md)的V、M、G与charge界；[singlet-kernel](../notes/susskind-singlet-kernel.md)的中心化A二阶统一界；[linear-kernel](../notes/susskind-linear-kernel.md)的r=q核与small-λ规范。旧Y外部中心时间核不等于M(s)的正交记忆核。本阶段第一次把指定平均记忆对象接到已有q核；不作优先权声明。

读取方式：原文在线核查，未新增原文PDF缓存；链接与版本明确。有限计算只核代数及数值实现，全阶结论依据稿中的标记计数、测度论与误差界。
