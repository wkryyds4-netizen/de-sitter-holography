# Susskind 主线的态、极化与可迁移的 dressed ladder

2026-11-07。先检索项目既有来源，再在线核读下列作者原文。此文件只增加新阅读记录及下一任务；未修改既有封存文件，未运行科学 main。本轮检索没有找到一篇已完成 MSS singlet 到介子 canonical 算符映射的后续论文；这只是本轮定向检索的结果，不是穷尽性断言。

## 1. 接续位置：已有结论不重复计作进展

已读 [rindler-packet 笔记](../notes/susskind-rindler-packet.md) 及其 [来源](susskind-rindler-packet-sources.md)，也读了 mainline、bulk-target、bilocal-channel 来源。以下均是旧成果：Rindler boost 的逆温度为 \(2\pi\)；midpoint 两点能配准而不固定质量；普通等热间隔四点给 \(33\ne3\)；有限矩阵的 \(A_n\) 随插入数变化；Lin–Stanford fake-circle 必须区分 bra/ket。此次不重做这些计算。

新阅读任务是确定：原方案究竟指定了哪一个态与哪一个代数，以及最近的纯态研究是否提供一种能在 complex、守电荷、同一 \(J\) 模型中实现的极化。

## 2. MSS 原文真正给出的输入

主文献是 Miyashita–Sekino–Susskind，[*Holograms and Standard Models*, 2607.05678v2](https://arxiv.org/abs/2607.05678v2)。摘要页确认 v2 于 2026-10-19 修订，当前最新可核版本仍为 v2；[HTML](https://arxiv.org/html/2607.05678v2) 与[实际 PDF](https://arxiv.org/pdf/2607.05678)页首版本一致。下表页码先列正文印刷页，再列 PDF 物理页，避免错位。

| 原文位置 | 明确输入 | 尚未由该处提供 |
|---|---|---|
| §2，(2.4)，p.5 / PDF 6 | \(H=i\partial_t\)，static/boost 时间 | 独立空间平移与光前 \(P^-\) |
| §2.1，(2.5)；§3.1，(3.4)，p.9 / PDF 10 | Boltzmann 态趋于规范化迹；它不是 Gibbons–Hawking 温度 | 保持乘法的真空态映射 |
| §§3.2–3.5，(3.14)–(3.17)、(3.27)–(3.29) | 同时旋转费米子与 \(J\) 的 ensemble 对称；讨论单样本自平均 | 微观时间局域的 Gauss 约束 |
| §§4.1–4.2，(4.1)–(4.18) | 候选 \(U(N)\) QCD₂ 参数字典与两个极限的区分 | 由微观 \(H\) 推导介子内部核 |
| §4.3，(4.19)–(4.26)，pp.23–25 / PDF 24–26 | 导数 singlet 族、CP 分型、守恒 \(Q/H\) | 定量 DSSYK 介子谱与 canonical dressing |

§2.1明确以尚未证明的对偶作假设；§4.3明确把定量谱比较留作检验。原文的导数简写仍须经本项目 CAR、相位与 \(p\) 约定转换，不能直接复制其 \(M_1=H\) 为当前归一化恒等式。

前驱 [MSS 2506.18054v3](https://arxiv.org/html/2506.18054v3) 本轮复读 §3、§4.1–4.2、§5.2–5.3、§6.1–6.2：§6 的体比较使用 Minkowski 真空在 Rindler 楔的热限制，(6.1)–(6.2) 给 boost 温度 \(1/(2\pi)\)。§5.2 的 fundamental 字典替代更早 adjoint 版本。§6 的储熵层位置须结合 [Susskind 2511.10907v1](https://arxiv.org/html/2511.10907v1) §§1–4 的修正；本轮沿用旧项目已记录的修正，不再展开。

**本轮判读。** 微观迹与体楔 KMS 态是两个明确写出的输入；在所核段落中，没有把二者通过一个保持同一算符乘法及真实时间流的映射连接起来。这个缺口不能仅用温度重命名补齐；旧 midpoint 四点反例已经说明只调两点不足。

### “nonlocal gauge theory”应怎样使用

2607 原文 §§3–4 给出的是从 hologram 对称到候选体 gauge theory 的论证。检索其 HTML/PDF 及前驱的 nonlocal/non-local，并核读相关章节，本轮未找到一个另行定义的“nonlocal gauge theory”作用量或微观 gauge fixing。不能据此给我们的双时间算符添加尚未构造的 Wilson line，或将 ensemble 的 \(U(N)\) 协变性当成已经对微观 Hilbert 空间做过 singlet 投影。真正的 gauge-invariant 体 bilocal 需要在指定体规范理论中另立定义。

作为清楚的对照，新核读 Zhang–Peng，[*Gauging the complex SYK model*, 2502.18595v2](https://arxiv.org/html/2502.18595v2)，§2 (2.1)–(2.19)、§4.4 (4.26)。此文确实加入 \(U(1)\) gauge field 与 charge-\(k\) Wilson line；holonomy 积分给固定电荷扇区的 partition function。它研究另一种具体 gauging，不是 MSS 的全部 \(U(N)\) 体对偶，也没有把固定电荷的无限温迹自动变成极化真空。v2 为 2025-10-09；本轮最初请求 v3 失败后，已按元数据改读 v2。

## 3. Xu：应区分场代数与辅助 ladder

[Xu, *Von Neumann Algebras in Double-Scaled SYK*, 2403.09021v7](https://arxiv.org/html/2403.09021v7)，2026-06-25。本轮复读引言、§3.2 末及 §3.3 (77)–(110)，特别 (79)–(84)、(86)–(87)。在原文参数 \(\rho=\max(|q|,|r|,|r_V|)<1\) 下，空 chord 向量对场代数 \(\mathcal A=\operatorname{vN}(H,M)\) 是 faithful trace 的 cyclic/separating 向量，模算符为 1。\(Q\)-twisted Fock ladder 实现同一 Hilbert 空间，但这不等于单独的 creation/annihilation 算符都属于该场代数。本文 v7 特别补强旧版本 factoriality 证明的无限 chord 尾步骤；本项目不重复声称重证其 factor 类型。

这给一个直接适用于选路线的标准推论：若 \(a\in\mathcal A\) 且 \(a\Omega=0\)，faithfulness 给 \(a=0\)。所以想得到非零真空 annihilator，必须明确改变可观测代数、态或二者的微观实现；对旧 Gram 矩阵作平方根归一化没有完成这件事。此推论是既有 faithful-GNS 事实的应用，不作为新定理。

## 4. 2026 年真正提供极化的前驱

Rajgadia–Xu，[*Emergent States and Algebras from the Double-Scaling limit of Pure States in SYK*, 2604.14387v3](https://arxiv.org/abs/2604.14387v3)，2026-09-29。读取 [HTML](https://arxiv.org/html/2604.14387v3) §§3–4 的相关定义/构造和 §5.1；并用[实际 PDF](https://arxiv.org/pdf/2604.14387)核对公式。**两种格式编号不同**：

| 内容 | HTML 编号 | PDF 编号与印刷页 |
|---|---|---|
| 用 \(2N\) Majorana 配成 \(c_i\)，KM 空态 | (34)–(36) | (3.2)–(3.4)，p.12 |
| spin-flip \(Q/p\) 与适应态的 \(M,M^\dagger\) | (37)–(43) | (3.5)–(3.11)，pp.12–13 |
| \(H\mapsto a+a^\dagger\)，\(b\lvert0\rangle=0\)，\([b,b^\dagger]=1\) | (70)–(75) | (3.38)–(3.43)，pp.20–21 |
| dressed 表示与全部固定多项式矩 | (85)–(86) | (3.53)–(3.54)，p.22 |
| generic 与适应态的代数区别 | §§4.1–4.2 | §§4.1–4.2，pp.24–29 |

可借用的核心不是一个新的两点核，而是原文已经构造的算符—态对：

\[
M^\dagger\longmapsto b^\dagger q^{\Delta\hat n_{\rm tot}/2},\qquad
M\longmapsto q^{\Delta\hat n_{\rm tot}/2}b,
\qquad \hat n_{\rm tot}=\hat n_H+\Delta\hat n_M.
\]

Generic probes 在纯 KM 态中的双标度矩仍给 tracial 代数；纳入与参考态相关的微观 \(M,M^\dagger\) 后，原文才得到纯态的 Type \(I_\infty\) 表示。§4 强调区别在系数相关性；§5.1 的二点也保留插入与态制备的位置。这不是只选择一个纯态、保持原 probe 不变就能达到的效果。

**直接移植为何不成立。** 原文 spin-flip \(Q\) 不是我们的守恒物理电荷。我们的 \(H=\sum J_{AB}C_A^\dagger C_B\)、\(p\ge1\) 作用于物理空态为零，中性双时间 \(\sum c_i^\dagger(t_1)c_i(t_2)\) 也湮灭空态。照搬“取 complex 空态”会得到平凡通道。新的工作要在有非平凡动力学的参考态上重新构造适应算符。

## 5. 本轮完成的微观极化及下一项物理任务

本轮检验了 **固定半填充 Slater 参考态与同一 \(J\) 的粒子—空穴极化**，并完成其纯 ladder 的固定 word 双标度极限。下面记录文献启发与已完成结果之间的接口；与原 bilocal 的完整共同矩及体字典仍是下一任务。

令 \(O\subset\{1,\ldots,N\}\)、\(m=|O|\sim N/2\)，\(\lvert O\rangle\) 为相应占据态。这里 \(J\) 已吸收物理方差，约定 \(\mathbf E J_{AB}J_{CD}=w_N\delta_{AD}\delta_{BC}\)，与[微观接口](../notes/susskind-slater-ladder-literature-interface.md)一致。在原 Hamiltonian 中保留完全从 occupied 到 empty 的块：

\[
R_O^\dagger=\sum_{\substack{A\subset O^c,\ B\subset O\\|A|=|B|=p}}
J_{AB}C_A^\dagger C_B,\qquad R_O=(R_O^\dagger)^\dagger.
\]

标准 CAR 立即给 \(R_O\lvert O\rangle=0\)，而 \(R_O^\dagger\lvert O\rangle\) 非零。这里每个算符仍守物理 \(Q\)，参考占据数定义了额外结构。若原 Gaussian 方差为 \(w_N\)，令

\[
L=\binom{N-m}{p}\binom mp,\qquad B_O^\dagger=R_O^\dagger/\sqrt{w_NL},
\]

则独立标签的两点计数给 \(\mathbf E\langle O|B_OB_O^\dagger|O\rangle=1\)。本轮已由有限 CAR/Gaussian Wick 计数建立完全 creation words 的精确公式

\[
\mathbf E\langle O|B_O^n(B_O^\dagger)^n|O\rangle
=n!\prod_{r=0}^{n-1}
\frac{\binom{N-m-rp}{p}\binom{m-rp}{p}}{L}.
\tag{P1}
\]

连续激发必须在粒子与空穴两侧同时避开已用标签，Gaussian 配对给 \(n!\)；完整推导见[接口 §3](../notes/susskind-slater-ladder-literature-interface.md)。在 \(m/N\to1/2\)、\(p^2/N\to\lambda>0\)、固定 \(n\) 时右边趋于 \(n!e^{-2\lambda n(n-1)}\)，即 dressed oscillator \(b^\dagger q^{\hat n/2}\)、\(q=e^{-4\lambda}\) 的 creation 范数。本轮进一步完成[纯 ladder 证明 §5，(8)、(12)](../notes/susskind-slater-ladder-proof.md)：全部固定 Dyck-word 矩取同一极限，实际 annihilation 离开固定 creation level 的范数趋零。有限 \(N\) 的 backward-pair 修正仍存在；(P1) 和上述纯 ladder 极限均不等于带完整 \(H\) 的交错矩，也不证明 canonical 介子存在。

下一轮的顺序应是：

1. 以本轮已证明的纯 ladder 作为输入，构造同一 \(J\) 的 \(B_O,H\) 混合词，保留共享耦合的交叉协方差。完整 \(H\Omega_O\) 的纯激发权重仅趋于 \(\sqrt q\)，不能用 \(B_O+B_O^\dagger\) 替代完整 \(H\)。
2. 把 return-deficit 扇区一并纳入联合状态空间，再核 \(\langle O|B_O H^rB_O^\dagger H^s|O\rangle\) 等交错矩。首步 return 权重趋于 \({\rm Poisson}(2\lambda)\)；下一步要判定其在多步动力学中如何组合，而不是重复纯 ladder 检验。
3. 求 \(B_O,B_O^\dagger\) 与现有完整 bilocal 及其 \(W_1,J_2\) 分量的共同矩，确认适应态信息到底进入何处。独立随机的新 probe 只能作辅助对照，不能替代这一接口。
4. 核参考态是否随真实 \(H\) 平移保持、需要什么制备；只有给出 stationary/KMS 的真实态与流后，才把此结构接到 Rindler 粒子解释。这里还没有建立 \(P^\pm\)、介子质量或有限 \(G\) 散射。

有一个标准代数检查应先保留：单独对固定 \(m\) 的参考态做 \(U(N)\) Haar 平均，会由 \(\wedge^m\mathbb C^N\) 不可约性得到

\[
\int dU\,U|O\rangle\langle O|U^\dagger=P_m/\binom Nm.
\]

它在中性扇区恢复迹态；不能一边平均掉参考信息、一边假定真空 annihilation 仍被同一固定算符保留。若恢复旋转协变性，应把参考态、\(J\) 与适应算符作为相关数据一起处理。此处是标准 Schur 引理的应用；不是声称现已完成规范不变的体构造。

## 6. 阅读边界与继续位置

新近原始文献的实质输入来自 2604.14387v3 的态适应构造与 2502.18595v2 的真实 charge gauging；2403.09021v7 是 2026 年修订的代数前驱。检索词包括 MSS 题名、DSSYK meson、nonlocal/non-local gauge、pure states、gauging complex SYK。检索返回的第三方解读和镜像没有用作论证来源。

本轮不新建 PDF 原文缓存；PDF 在线版本与 HTML 都已核。MSS 后续定量介子谱仍需持续查新，但当前最短的探索路径是上面的可证伪微观极化接口，而非再次扫描既有两点温度或重新复现已知 ’t Hooft 谱。
