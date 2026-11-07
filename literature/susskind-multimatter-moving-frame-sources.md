# 多物质移动框架：新增前驱与读区

2026-11-08。在线核查后先读原文，再写计算。本轮发现此前项目遗漏的重要后续文献，保留旧sealed稿的历史状态，以新记录更正继续方向。

- Lin–Stanford [2307.15725v1](https://arxiv.org/html/2307.15725v1)：本轮复读§2.3 (31)–(41)、§3.3.1 (72)–(104)、§3.3.2 (105)–(110)，以及App.F (354)–(355)的措辞。裸concatenation coproduct不是普通Hilbert乘积上的primitive加法；ordered direct Gram及outer action沿用原文。(72)的r1按邻近定义理解为rV。旧猜测证据不自行升格。
- van der Heijden–Verlinde–Xu，*Quantum Symmetry and Geometry in Double-Scaled SYK*，[2511.08743v2](https://arxiv.org/html/2511.08743v2)。[arXiv版本页](https://arxiv.org/abs/2511.08743)显示v2为2025-11-17修订。本轮读§2.3 (30)–(31)、§3.1 (40)–(57)、§3.2 (58)–(65)、§3.3 (66)–(79)及(82)–(83)、§3.4 (87)–(101)；相关diagram解释及(85)也已读，未将所有diagram展开代码当作逐图核验。核心输入为(55)/(92)的factorization、(62)/(69)的coproduct和(82)/(83)的swapping channel transform。其q_HVX²=q_project，d_i=Delta_i。文中(57)的aL/R是在另一双振子实现中的算符，不能直接同名替代本项目LS aL/R。
- 该文footnote16明确没有直接等同其coproduct与LS原concatenation，故我们的连接仍用direct Gram及outer lowering证明；不把该脚注的未解对应当成已知恒等式。HTML (41)的K/H记号与(40)/(44)存在表面因子不一致，能量intertwiner (67)/(75)/(77)也不直接拿来做数值归一化；本轮只使用读清的表示结构，所有具体frame norm独立推导。此处不是全论文勘误结论。
- (82)–(83)通过Askey–Wilson变换导出fusion-channel swapping权重，足以更新旧two-matter稿“仅LS有限级guess”的来源状态。该论文的能量态归一化、所引Jafferis变换与本项目same-species Gram的全量逐式对应留作下一步；当前只采用旧独立证明的k0,k1 physical norms。
- 已有项目输入：notes/susskind-two-matter.md及其41项sealed检查；notes/susskind-chord-module-frame.md及有限词/Poisson尾证明。无旧driver重跑。

论文提供已有量子群factorization与fusion，不把它们称为本项目发现。本阶段新增连接是共同位移下实际LS生成元的ordered additive finite-word极限及正强对易候选的有限q搬运。未读§4的Schwarzian细节、未重建完整6j核，不据此作质量/引力相互作用结论。
