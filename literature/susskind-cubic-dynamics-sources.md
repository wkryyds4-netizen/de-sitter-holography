# 三阶一次残量动力学：来源与实际使用

2026-11-06，先复读现有公式，再确定计算；不缓存第三方全文，不公开过程记录。

- Kota，[2208.11312v2](https://arxiv.org/html/2208.11312)，§IV (28)–(31)、(37)–(39)：实际核对 U(N) 不可约 covariance 的两个通道、维数与Λ。用于确定ω及纯rank map的三角方程。文章的谱密度two-point不是本项目F的transition-strength correlator。
- Tall–Tomsovic，[2604.14522v2](https://arxiv.org/html/2604.14522v2)，§2.1 (2.8)–(2.13)，另核(2.19)：实际核对正常序计算和不可约tensor上的完整embedded covariance。本文没有直接引用其独立probe结论替代same-J F。更高点仍需控制纯rank产品分解。
- [DLMF 18.26.3](https://dlmf.nist.gov/18.26.E3)：Racah的终止₄F₃定义。本轮参数α=−m−1，β=m−N−1，γ=−N−2，δ=0；无需套用其一般正参数正交测度。
- [DLMF 16.4.9](https://dlmf.nist.gov/16.4.E9)：Rogers–Dougall very-well-poised ₅F₄的终止求和，用于ω的阶乘矩。实际用有理有限和形式，避免逐项代入Gamma极点。移a=k+j后 A′=−N−1+2k,b′=−N−1+k,c′=m−N+k,d′=−(r−k)。
- [DLMF 16.4.3](https://dlmf.nist.gov/16.4.E3)：Pfaff–Saalschutz终止平衡₃F₂和，把Racah的ω平均还原为已知κ/C。参数n=r,a=−b,b_src=b−N−1,c=−m,d=m−N−r。

新写入的是这些标准结果在当前F滤波上的应用：精确same-J动态矩、DS极限及频率权重控制。Gaussian矩界、Young、Cauchy–Schwarz和Paley–Zygmund是标准工具；稿中给了所需Gaussian矩界的简短直接证明。没有把已有协方差、特殊函数恒等式或不等式改名为原创。

检索还定位Kota nucl-th/0410069；本轮使用并核读的是上述2023文章中的明确公式，没有把未重读的2004全文作为新增证明依据。前轮t_a与Poisson尾控制引用本地cubic-linear原稿，科学文件原字节保留。
