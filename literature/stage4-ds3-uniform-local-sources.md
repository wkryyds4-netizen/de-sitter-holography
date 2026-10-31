# 跨零合并与局部积分：原文依据

2026-10-31。保留已有全文缓存原字节。

- [Liu et al. 1808.00612v3](https://arxiv.org/html/1808.00612v3)，(3.35)–(3.38)：沿用上一轮已逐项核对的两shadow及两4F3核；本轮将其中Gamma可去极点合并后估计，不把原公式称为原创结果。
- [NIST DLMF §5.7.6](https://dlmf.nist.gov/5.7)：本轮核读digamma差级数，结合首项加积分比较推得σ>0下的显式|ψ(σ+it)|上界，用于正则Gamma系数的导数常数。
- [NIST DLMF §5.11](https://dlmf.nist.gov/5.11)：核对logGamma与digamma渐近的sector及余项条件。C¹一致展开由两者分别控制，不靠对未知点态余项形式求导。大参数离开负实切口，小参数保持紧集。
- [Groenevelt math/0306424v1](https://arxiv.org/html/math/0306424v1)，Proposition4.7证明及(4.4)：本轮核读可去项、Riemann–Lebesgue与Dirichlet积分的作用。采用标准测试积分论证，但本项目的物理核一致误差另行推导；没有移用Wilson正测度下的酉性来证明物理标量变换完整性。

局部dr分布结论和真实物理谱积分不同。后续必须核对实际质量权重、时间相位、非紧谱尾及极限交换，不能由本文局部定理直接宣布Bose自由极限完成。
