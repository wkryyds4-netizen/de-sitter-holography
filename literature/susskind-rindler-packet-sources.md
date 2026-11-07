# Rindler 波包与多点轮廓：本轮原文对应

核查日：2026-11-06。所列为实际打开并阅读的原文位置；没有以检索摘要代替关键公式。已有结果归原作者；本轮只作具体参数应用与推论。新来源在线核查，未另存原文 PDF。

| 来源 | 已核位置 | 用途与边界 |
|---|---|---|
| Miyashita–Sekino–Susskind, [2506.18054v3](https://arxiv.org/html/2506.18054v3) | §5.3 (5.22)、§6 (6.1)–(6.3)、§6.2 | singlet tower/meson 目标，Rindler flat limit 与 T_U=1/(2πρ)；这些段落没有当前 X 的径向重建。 |
| Rahman–Susskind, [2401.08555v3](https://arxiv.org/html/2401.08555v3) | §4.3 (4.42)–(4.44)、§5 (5.45)–(5.47) | cord 核参数 J0、T_cord=J0/π，vT_B 与 T_B 区分。当前 J0=2J_s 是按两点核的参数匹配，不自行推出同一个物理 proper clock。 |
| Crispino–Higuchi–Matsas, [0710.5373](https://arxiv.org/pdf/0710.5373) | PDF §II.C–E、PDF第10页 (2.92)；HTML §II.4/§III.1.2 | K_iν 模式、KG归一、Minkowski/Rindler热态。原式是四维，有横向 (2π)^{-1} 因子；本轮二维归一由 KG 与 KL 再核。HTML全局公式编号与PDF节内编号不同。 |
| [DLMF 10.43(v)](https://dlmf.nist.gov/10.43.v) | (10.43.30)–(10.43.32) | KL 正逆变换及点态条件。本轮用下一来源的 L² 版本，不假定自动满足 DLMF 所列每条点态条件。 |
| Yakubovich, [Uncertainty principles for the Kontorovich–Lebedev transform](https://doi.org/10.3846/1392-6292.2008.13.289-302) | pp.289–290 (1.1)–(1.3)、§4 (4.1)，原PDF已读 | forward 测度 dy、范数 y dy；取其输入为我们的 f(y/m)/y，得到 dρ/ρ 与2/π²的精确因子。 |
| Lin–Stanford, [2307.15725](https://arxiv.org/html/2307.15725) | §4.2 (158)–(164)、§4.3 (167)–(170)、§4.4、§4.5 (181)–(183)，§4.6开头 | bra/ket不同fake-circle映射、timefold、双对易子和boost生成元。§4.6完整推导与附录G尚待下一轮精读/代入，不能标成已经计算。 |
| Maldacena–Shenker–Stanford, [1503.01409](https://arxiv.org/pdf/1503.01409) | Appendix A，PDF pp.16–19，(31)–(35) | 左右楔时间方向、β/4移动与Minkowski ordering、用外波函数折叠散射。文中 d>2 的因子化论证不能直接移植到当前1+1维。HTML编号比本PDF相应编号大1，记录时以各自版本为准。 |

## 本轮真正推导的内容

把已知 sech⁴ 核做 Beta/Fourier 变换；将其正谱送入已知自由Rindler/KL公式，得到任意质量的Gaussian调节波包；算出一致ε误差、全空间范数及质量/径向伸缩自由度。然后对该具体自由场做普通等热间隔Wick收缩，得到η=0时33的极限，与微观第四矩3比较。这是候选字典的具体诊断，不当作新普适定理。

## 下一轮先用的文献公式

Lin–Stanford §4.3 (170) 与 §4.6/附录G，是固定时间首阶λ四点的直接前驱。需要先匹配λ_src=4λ、Δ=2、bra/ket位置与same-species三个通道；原文两种探针记号不自动保留同种的额外收缩。Rindler波包外态和fake-circle位置的关系仍需构建。当前不能从 Tricomi 的实Laplace权直接认定unitary eikonal相位或λ=G的物理量纲。
