# 中央四Gram的来源与路线变化

本轮实际在线回读Tall–Tomsovic arXiv:2604.14522v2 §2.1的fixed-charge primitive decomposition及covariance式(2.11)–(2.19)：https://arxiv.org/html/2604.14522v2 。此前已有的embedded-ensemble协方差框架继续作为输入，不将其重新包装为新理论。

直接引用的项目内已证输入：

- `manuscript/susskind-memory-residual-norm/main.tex` 第2节：neutral exchanged Casimir=(x_a+x_b)χ，基于self-adjoint trace covariance的明确证明；不是上一轮尚缺的charged exchanged式。
- `manuscript/susskind-memory-central-window/main.tex`：中央正输出六阶矩。
- `manuscript/susskind-memory-central-norm/main.tex`：实际O范数、输入加权尾、(B,L_in,O)Gram、完整core强展开。
- `manuscript/susskind-memory-central-cross/main.tex`：等距I_m、一步产品差O(1/N)、S传递O(p^−2/3)、(B,L_in,S)Gram。

本轮关键工作：对真实同J算符O=tV−h0U作极化，利用exact normalized weight ratio与neutral Racah contiguous关系证明F_(m−1)−F_m=O(1/(Np))，以此完成所需cross极限。所有新验证脚本都不调用候选charged first Casimir公式。

强于当前目标的有限N候选g[(x_a+x_b)ξ−m(ξ−χ)]仍缺general rectangular gradient的完整解析证明，保留旧文件的待证标记。改走极化并不是默认该式成立，而是独立完成了本阶段目标。下一步不再把该有限N命题作为energy插入研究的前置条件；若动态计算真正需要它，再回到原文recoupling结构核实。
