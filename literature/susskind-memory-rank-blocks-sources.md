# rank 乘法与尾控制：前驱核读

2026-11-06实际核读：

1. [Kota–Vyas 1504.00759](https://arxiv.org/pdf/1504.00759)，§3，页14–15的(34)/(35)：不可约张量的耦合乘积含重数标签；其reduced matrix element通过U(N) Racah recoupling表达。本文后面的dilute limit与本项目半填充DS不同，本轮没有照搬其渐近。工具可读完整PDF提取文本，截图及HTML fallback cache miss，未声明视觉核读。
2. [Kota 2208.11312v2](https://arxiv.org/html/2208.11312)，§IV协方差两种耦合、Λν与dν；复用其矩阵协方差结构，不把谱密度two-point当真实memory two-point。
3. [DLMF 18.26.3](https://dlmf.nist.gov/18.26.E3)，终止Racah多项式定义。χ_a(s)的参数识别沿用本项目cubic-dynamics阶段，α_DLMF=−m−1，β_DLMF=m−N−1，γ_DLMF=−N−2，δ_DLMF=0。矩递推的probe阶数不得越过R=min(m,N−m)。

本轮在已有公式上的应用：Casimir projectors保留multiplicity的完整二次tensor Gram；协方差复合给正output-rank概率及Racah三重和；由前三个Casimir矩推固定λ输出deficit的矩极限和尾界。没有声称新发明Racah线性化、Wigner–Eckart、Gaussian decoupling或Markov不等式。

本地前驱：notes/susskind-memory-quadratic-contraction.md 的精确Q_Ecal，notes/susskind-cubic-dynamics.md 的纯rank χ及变换，notes/susskind-cubic-linear.md 的ω尾控制。本轮未证明memory filter的一般reduced coefficients或完整Poisson定律；输出π仅为independent isotropic products，same-J对应的是已写出的sym上界。
