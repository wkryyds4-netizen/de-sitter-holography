# 双局域四点：已有弦公式、共同耦合标记与物理解释

2026-11-07。先从项目的 singlet-multimark、mixed-cubic、memory-size-bridge、memory-wick-descendants 与 memory-to-bilocal 来源记录定位，再在线核读下列原文。本页没有运行旧 main，也不把两点谱等价提升成四点定理。

当前对象沿用标准 CAR、共同 complex Gaussian p-body 耦合 J、β=μ=0 的归一 annealed full-Fock trace、ℏ=1/string 时间；先固定次数和时间，再取 p²/N→λ>0，q=e^(−4λ)、C=J_s²/λ。关心 J₂ 及线性 size 残量

\[
W_1=(p-\mathcal L)H-\lambda H,
\qquad \mathcal L\Gamma_A=\tfrac12|A|\Gamma_A.
\]

此处 W₁ 是场名，不能与 coupling-Gaussian 第一 chaos 的投影符号混用；有限 N 回归若采用 λ_N 或精确系数，须保留到相应误差估计完成。

## 可直接复用的原文

1. **BNR，complex DSSYK**：[Berkooz–Narovlansky–Raj, 2006.13983v2](https://arxiv.org/html/2006.13983v2)。本轮读 §3.2–3.3 (3.11)–(3.15)、§7 (7.1)–(7.3)、§8.1 (8.1)–(8.7)。前者给一般 index-set 长度的有向 CAR 局部迹及 pair-Poisson 极限；后者明确 probe 系数独立于 H，且不同 flavor 系数独立。中性四点分 crossed/uncrossed channel；μ=0 时 H–probe 交叉权为 exp(−4pp_M/N)，probe–probe 的互交叉因子来自 (8.4) 的 exp(−4p_M1p_M2/N)。这些公式可以作为已证共同耦合归约之后的答案，不能直接替代归约。

2. **BINT，通用多点与 R-matrix**：[Berkooz–Isachenkov–Narovlansky–Torrents, 1811.02584v2](https://arxiv.org/html/1811.02584v2)。本轮核 §4.1 (4.9)、§4.2 (4.10)–(4.13) 及 §5 (5.6)。已有 uncrossed 积分、crossed transfer series、量子群 6j/R-matrix；§4.2 开头及脚注 8 明确其展示式省略互 probe 因子，而同种算符完整四点不能省略它。若已经证实某个单位归一场的交叉矩阵为 HH:q、HX:q²、XX:q⁴，应恢复 q⁴ 并合并三个 Gaussian 收缩 channel；不必重新推导 R-matrix。

3. **mixed-q 的正迹实现**：[Junge–Zeng, *Ultraproduct methods for mixed q-Gaussian algebras*, 1505.07852v2](https://arxiv.org/html/1505.07852v2)，§3.1 Theorem 3.1 与 (3.2)。对称 crossing 矩阵 Q 的条目在 [−1,1]，矩由同色配对及交叉 Q 因子给出，形成正的 tracial 表示。这里“矩阵”指固定的交叉参数矩阵；它不是本项目能量变量上的正矩阵谱 ν_ab。拥有后者只给共同两点 Gram 空间，不会自动给前者的 q-Wick 多点律，更不是 canonical boson 代数。

4. **复合块的 q-Wick 乘积**：[Effros–Popa, *Feynman diagrams and Wick products associated with q-Fock space*, math/0303045](https://arxiv.org/pdf/math/0303045)，Theorem 3.2（PDF 第 8–9 页）、Theorem 3.3/(11)（第 10 页）。前者给多块谱 q-Wick 的期望：禁止块内连线，保留块间所有完整配对及交叉数；后者给允许剩余 singleton 的算符乘积。它解释为什么“复合块”不等于一个新的独立 Gaussian 字母。有限 N 的 coupling-Gaussian Wick ordering 和这个谱 q-Wick ordering 仍须分别建立接口。

5. **小形变及长时间四点**：[Mukhametzhanov, *Large p SYK from chord diagrams*, 2303.03474v3](https://arxiv.org/html/2303.03474v3)，本轮核 §5.1 (5.1)–(5.10)、§5.2 (5.11)–(5.17)。原文从已有精确 6j 公式取小形变/平均能量缩放，区分欧氏时间与长 Lorentzian 时间的相位抵消。它适合在有效 probe 已识别后分析 OTOC，不提供当前 W₁、J₂ 的共同 J 标记四点。当前固定 λ、固定时间问题应先留在 BINT/BNR 精确框架，不把该长时间近似当作本轮证明输入。

6. **canonical dressing 的态与极化**：[Itakura, hep-th/9604032v2](https://arxiv.org/pdf/hep-th/9604032)，本轮复读 (2.15)–(2.18)，dressing 接口沿用 [memory-to-bilocal 来源](susskind-memory-to-bilocal-sources.md) 对 (4.18)–(4.29) 的核读。原文 bilocal 中心项由 Fock 真空正规序和动量符号分区产生；Holstein–Primakoff 以 \(\sqrt{N_c-A}B\)、\(A\sim B^\dagger B\) 处理非线性代数，并有物理子空间限制。应先确定这些结构，再比较介子生成元；本项目无限温迹与原文真空态不相同。

## 先处理项目已有的 J₂ 前驱

[singlet-multimark](../notes/susskind-singlet-multimark.md) **已经声称**旧
\(X=(O-cB)/\sqrt{W_d}\) 与 H 的任意固定阶 mixed-q 归约，并记录
\(\tau(X^4)=2+q^4\)。其 \(W_d=C^2[4\lambda+16\lambda^2q/(1+q)]\) 与当前 J₂ 二点规范相符。它是本轮必须核对的历史结果，不能把重新计算同一四点称为新预测。

接入需要核当前 J₂ 与旧 Y₋ 的实际共同耦合定义、charge 背景和有限 counterterm；只有方差相同还不够。若两者仅有 L² 强等价，四点转移仍需适当的高矩一致可积性或乘积估计。旧多标记证明对所有同 site 标记系数及交换轨道的处理，也应由本轮独立四源引擎核对；如果出现不符，保存旧输出，在新的 journal/proof 中说明具体修正，而不静默改写旧证书。

对于 W₁，本轮所读原文没有给出这个特定 Clifford-size/common-J 残量的现成四点闭式。这个判断仅限上述针对性核查，不是不存在文献或原创性声明。可以直接继承局部 CAR、Gaussian 收缩、Poisson 交叠与 q-Wick 组织；待做的是它们在实际四个外部插入上的组合。

## 四源计算不能提前丢掉的项

BNR 的“三重交叠消失”是三个独立底层 index-set 的交叠结论。经过 Gaussian Wick 配对，同一集合会在多个外部插入重复出现；两个底层集合的一个交叠 site 可以同时影响四个外部 field。因此四源生成函数应先保留局部 ratio 的全部 nilpotent 系数，包括可能的 z₁z₂z₃z₄ 项，再做中心化和图求和。两源的 λK₁₂ 并不规定四源函数只能是二次式。

同样，局部未标记迹为零的布局，加入 Clifford/有向 source 后未必仍为零；符号和局部矩阵必须在求导前保留。这些要求来自当前对象的有限 CAR 结构；不能用独立 probe 终式或把 body size 当一个固定整数来绕开。

推荐将新四点和三个不同基准分开比较：实际 common-J 有向 source；已存 J₂/旧 Y₋ 的双弦归约；假设独立 mixed-q 字母的参考值。最后一项是诊断基准，不是新计算的定义。

## 四点之后最直接的物理问题

先判断重建出的 matter 场是否能用一组简单 mixed-q 生成元描述，还是必须保留复合占据结构和非线性 dressing。对两个 Hermitian、已减背景场 A、B，tracial 恒等式给

\[
\tau([A,B]^\dagger[A,B])
=2\{\tau(A^2B^2)-\operatorname{Re}\tau(ABAB)\}\ge0.
\]

在单位范数、互相正交的 mixed-q **参考模型**中，右边为 \(2(1-Q_{AB})\)；同一单位场的普通四阶 cumulant 为 \(Q_{AA}-1\)。因此同时计算顺序四点、交替四点、交换子范数和二次乘积 Gram 矩阵，比只报一个四阶 cumulant 更能判断代数。非零交换子不受 \(\tau([A,B])=0\) 排除；偏离参考模型也不直接意味着体引力相互作用。

若参考归约成立，就调用已有 BINT 多点公式，继续寻找适当态/极化下的 canonical 变量；若不成立，就用四点中保留的共同 J 结构构造更大的复合场代数，再检验 Itakura 型占据依赖修正。这是本轮四点诊断的物理用途。微观纵向变量、实际 P⁻ 与体内部质量核仍须独立建立；当前中心时间或 memory 频率不被重命名为介子质量。

本页使用 arXiv 原文在线核读；未新增整篇缓存、未作宽泛优先权检索。四点最终数值或定理由本阶段的新 proof/results 给出，本页不预写其结论。
