# 固定电荷双局域：前驱与接口

2026-11-08。本阶段使用已有构造再补逐扇区强接口，不将 q-Gaussian 条件期望、OU 公式或 Volterra 恒等式列为新发现。

- Bożejko–Kümmerer–Speicher, *q-Gaussian processes: non-commutative and classical aspects*, [funct-an/9604010](https://arxiv.org/abs/funct-an/9604010)，原文 Theorems 2.11、2.13：正交投影的第二量子化、保迹条件期望及 C* 版本。PDF 正文相应位于第12–14页附近；2.10 是 factoriality，不能作为条件期望定理号。用于两个 whole auxiliary blocks 的 replica，不意味着任意交错经典独立。
- Berkooz–Narovlansky–Raj, *Complex Sachdev-Ye-Kitaev model in the double scaling limit*, [2006.13983v2](https://arxiv.org/html/2006.13983v2)，固定电荷 complex SYK 的背景。独立 probe 弦规则不自动替代本项目同一 J 的 body-loss 与 neutral 后代。
- [central-ordinary](../manuscript/susskind-central-ordinary/main.tex) §§2–5：完整 Fourier 支配和普通内部 Gaussian 收缩。取所有 marks 为空即可逐 pairing 证明半填充多色矩；颜色只删配对。这比引用旧 two-component 中一句 compact-central 判断更直接。
- [sector-gradient-tail](../manuscript/susskind-sector-gradient-tail/main.tex)：实际 ordinary 系数所有固定范数的阶乘尾，允许固定导数和 exterior H。
- [sector-neutral-transfer](../manuscript/susskind-sector-neutral-transfer/main.tex)：中性 actual/local 全时间强替换、central pinching、W 选态边界。
- [sector-memory-return](../manuscript/susskind-sector-memory-return/main.tex)：finite Jacobi 守恒参考、实际 Fscr−X0 的 O_r(p^-1/2) 返回，以及 force 高矩。
- [higher-wick](../manuscript/susskind-slater-higher-wick/main.tex) 与 [polynomial-density](../manuscript/susskind-slater-polynomial-density/main.tex)：完整 marked source 塔、固定 Gaussian degree 的 selected 制备密度及 whole covariant block 传递。
- [memory-two-component](../manuscript/susskind-memory-two-component/main.tex) §§2–4 与 [operator draft](../notes/susskind-bilocal-fourth-tail-operator-draft.md) §8：有限 alpha OU identity、四 u 展开、有限时间回归尾。本阶段另证 half-filled 三色桥，不由 full trace 限制得到。
- [bilocal reconstruction proof](../notes/susskind-bilocal-reconstruction-proof.md) §§1–6：static memory 自伴 dilation、精确 Volterra 比较及中心化时间积分。这里改用实际总记忆的所有固定范数界，得到 O_r(1/p) 传播误差。

新增检查与这些解析输入互相补充，有限有理/Fock 比较和小矩阵 CAR 输出不能单独证明双标度极限。最强完成范围为半填充实际双局域的固定有限协变词在已控制多项式制备之间的收敛；尚无完整 de Sitter 空间/质量或 finite-G 字典。
