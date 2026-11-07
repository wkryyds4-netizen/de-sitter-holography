# finite-time rank/chaos 与 charged 四点接口：来源核查

2026-11-07。本轮先核已有四点与协方差，再构造全次数估计。

- Kota, *Bivariate moments of the two-point correlation function for embedded Gaussian unitary ensemble with k-body interactions*, [2208.11312v2](https://arxiv.org/html/2208.11312v2)，§IV(28)–(33)附近的不可约张量协方差、维数与两种耦合。用途是既有Φ/κ/rank分解，不把多trace谱密度关联当本项目单迹内算符协方差，也不搬其dilute渐近。
- [DLMF18.26.3](https://dlmf.nist.gov/18.26.E3)：终止Racah定义；本项目χ参数沿既有rank-product笔记。该文献不直接给本轮uniform-k缺额界，后者从真实正卷积和条件方差推导。
- Berkooz–Narovlansky–Raj, [2006.13983v2](https://arxiv.org/html/2006.13983v2)，本轮重读§3.2(3.11)–(3.12)明确任意N/p_i的精确site trace与方向交替条件；§3.3(3.13)之后才忽略三重交叠。§7(7.1)–(7.3)是独立Gaussian probe，§8终式为neutral probes。线性Mξ的精确系数平均可使用该思想，但charged size-one与1/N精度须回有限迹，不可只展开领先tilde q。λ_N/N的标记双命中尺度是本轮计数，不称BNR已给完整该阶余项。
- Aoun–Banna–Youssef, *Matrix Poincaré inequalities and concentration*, [1910.13797v3](https://arxiv.org/html/1910.13797)，Theorem1.2及§5 Prop5.5/其证明：Gaussian矩阵Poincaré、OU semigroup与梯度form。读后作为备用路线：对D=mΛ−B_s(h)，U(N)协变使ED为scalar，pinching contraction再接梯度平方和。尚未证明减谱背景后的该梯度和O(1)，未把Poincaré本身说成答案。最终普通指数定理直接写出有限Gaussian Hermite/creation代数与Duhamel证明，不靠未经核读的额外定理。

协作文献检查还排除了直接套用：[Bulycheva1706.07411](https://arxiv.org/html/1706.07411)的fundamental complex四点/共形ladder不提供固定λ的所需一致误差；[Pluma–Speicher1905.12999](https://arxiv.org/html/1905.12999)讨论多trace经典cumulants；[Tall–Tomsovic2604.14522v2](https://arxiv.org/html/2604.14522v2)§4.1 probe为neutral EGUE(t)。这些均非本轮证明依赖。

为Gaussian Fock规范另外检索了Janson *Gaussian Hilbert Spaces*及Nualart教材官方目录，只见到目录/简介，未将其当作已读定理证明。稿内所用a†/a和Hermite导数等距都是有限Gaussian代数，可直接验证；本地完整系数检查覆盖H⁴所有chaos。

本地前驱：memory-rank-blocks的两输入正卷积与Racah反演；memory-return-expansion的高矩与charge-tail；memory-return-background的实际场→X0；memory-fluctuations的共同J等时W。新范围为全k统一rank链、一步缺额、完整普通指数rank–chaos缺额及实际charged回归上界。不作原创性/优先权声明；没有完成X0有限时间谱或体介子重建。
