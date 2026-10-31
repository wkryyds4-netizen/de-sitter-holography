# 任意紧转移与物理尾控制：来源核查

2026-10-31。本轮按原文公式与已有项目规范核对；无新全文缓存，不作穷尽查新或首创声明。

- [NIST DLMF 4.36](https://dlmf.nist.gov/4.36)，4.36.1 的 sinh 乘积与4.36.2的 cosh 乘积。后者给半整数 Pochhammer 任意实虚部的正乘积上界；这是标准恒等式的应用。
- [NIST DLMF 5.11](https://dlmf.nist.gov/5.11)，log-Gamma、digamma 与 shifted ratio 渐近。用于固定紧移位和固定扇形的余项；不将点态大O形式求导当作一致导数证明。
- [Derkachov–Spiridonov, arXiv:1711.07073v2](https://arxiv.org/html/1711.07073v2)，§3.1正交/完备，式(22)、(29)、(30)；§4式(33)、(34)为重耦合基变换。完整空间包含自旋与等价表示识别；本项目正标签与球周期归一在 two-matter/crossing-dictionary 中已换算。只由完整酉变换推出 scalar compression 的范数≤1，不声称 scalar 块自身酉或逐点核≤1。
- [Liu–Perlmutter–Rosenhaus–Simmons-Duffin, arXiv:1808.00612v3](https://arxiv.org/html/1808.00612v3)，式(3.35)–(3.38)的完整交叉公式。新计算保留全部 shadow/伴随项和物理周期相位；公式本身不是我们的原创。
- [Groenevelt, arXiv:math/0306424v1](https://arxiv.org/html/math/0306424v1)，Proposition4.4 的 Wilson 对偶及 Proposition4.7 的可去极点/Dirichlet 机制为前驱；没有把 Wilson 的独立 Hilbert 测度直接移植为物理质量测度。

直接项目依赖：[任意紧物理权重](../notes/stage4-ds3-physical-local.md)、[跨零核](../notes/stage4-ds3-uniform-local.md)、[横向相位](../notes/stage4-ds3-transverse-phase.md)、[完整伴随](../notes/stage4-ds3-peak-limit.md)、[真实范数与非交叉](../notes/stage4-ds3-two-matter.md)、[Bose配对](../notes/stage4-ds3-bose.md)。

本轮推导：将小窗口系数界替换为 cosh 常数，逐项核其余紧移位论证；以实际两链平滑截取与质量一致 Hilbert 尾拼接迭代极限。自查范围仅为固定 s>0、G=0、半周期调节的两次创建扇区。增长时间、全 Fock、有限 G 和微观全息字典未完成。
