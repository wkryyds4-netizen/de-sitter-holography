# dS₃ 群配对、线性配对态与两点：来源核查

2026-10-31。为[本轮推导](../notes/stage4-ds3-group-pairing.md)保存证据；不是系统性的首创性查新。

| 原文及固定版本 | 实际核查位置 | 本轮用途与边界 |
|---|---|---|
| [Yang–Zhang–Zheng, 2505.08647v1](https://arxiv.org/html/2505.08647v1) | §2.2 (2.23)–(2.27)、Appendix A (A.1)–(A.4) | 直接给径向 sinh²t dt 规范的 dS₃ 密度 2s²/π；本轮正弦变换是归一核对，不是新密度。其 HH 权重与大质量球面鞍点匹配，不等于本轮精确线性 gluing 判据 |
| [Kolchmeyer–Liu, 2411.08090v1](https://arxiv.org/html/2411.08090v1) | (4.20)、(6.1)–(6.10)、(6.13)–(6.21) | dS₂ singlet、左右插入、sech(πs) 态及真实质量时间相位的直接前驱。跨维需重算群密度和 CG 范数 |
| [Marolf–Morrison, 0810.5163v1](https://arxiv.org/html/0810.5163v1) | §2.3、boost 矩阵元及收敛讨论 | KAK、固定 K 型径向衰减。双粒子固定质量不能当绝对收敛积分 |
| [Marolf–Yu, 2610.24073v1](https://arxiv.org/html/2610.24073v1) | §§5.1–5.2、两粒子质量涂抹相关脚注29 | 高维 tempered 群平均、正性与紧质量测试。不是本轮新定理 |
| [Derkachov–Spiridonov, 1711.07073v2](https://arxiv.org/html/1711.07073v2) | 标量三点 intertwiner、反射及 CG 正交 | 泛型主系列标量不变三线性形式的唯一性；不直接移入其原始 CG 测度 |

YYZ 缓存使用 v1 标题 **Comments on the de Sitter Double Cone**。后续发表版本更名不回写这个版本。新缓存为 [原文](../sources/yang-zhang-zheng-2505.08647v1.html)，其余沿用既有字节。下载首次遇沙箱 DNS 失败，随后通过审批下载公开原文。

YYZ (2.23) 写球函数范数为 δ(s−s′)/ρ(s)，(2.25) 在 d=2（即 dS₃）给 ρ=2s²/π。Appendix A 用 c-function 固定常数。本轮从 φ=sin(sr)/(s sinh r) 得同一结果；时空谱 μ=s tanh(πs)/π 则来自另一内积，不可混同。YYZ (2.26) 的态系数采用不同规范，不能只比较符号名 ρ_HH。

本轮新的项目级进展，是把此前已固定的 BD 时空模式、局部插入与本轮物理群范数组装，得到指定 G=0 模型的单物质正核及旧 recoil ansatz 的 n_loc。没有把标准 Plancherel、tanh 积分、三点唯一性或 KL 的 sech 态称为原创。完整物理构造与发表层面的优先权还需另行核查。

只浏览了 Chen–Xu 2026 年协变观察者代数文章的摘要，其参考系模型使用群 L²；未用于本轮公式或证明，不登记为已精读来源。
