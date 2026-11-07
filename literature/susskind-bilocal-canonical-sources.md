# Canonical bilocal 的参考态、极化与适用范围

2026-11-07。先读既有 `susskind-bulk-target`、`susskind-bilocal-interface`、
`susskind-memory-to-bilocal-sources` 与 `susskind-bilocal-channel-sources`，
再在线核对以下原文。没有修改旧来源缓存、科学输出或 mains。
本轮的有限 CAR 推导与本模型适配见
[polarization interface](../notes/susskind-bilocal-polarization-interface.md)。

## 本轮实际核读的主要来源

1. **Kazunori Itakura**，*Boson Expansion Methods in (1+1)-dimensional Light-Front QCD*，
   [hep-th/9604032v2 PDF](https://arxiv.org/pdf/hep-th/9604032)，首页版本日期1996-08-11。
   本轮直接读 (2.15)–(2.18)、(3.1)–(3.4)、(4.18)–(4.29)、Appendix (A.1)–(A.6)；
   §5 的介子投影沿用已核笔记。
   (2.15) 的中心项依赖 Fock 真空正规序和动量符号，(2.17) 指定湮灭真空的三个块。
   HP 表示以 canonical B 及 sqrt(N_c−A) 的算符核实现原双线性；有限颜色有物理子空间。
   N_c→∞ 的近真空展开保持激发占据相对 N_c 小。这里 N_c 是颜色数，不能改名为本项目
   complex 模式数 N。原光锥质量规范 M²=2P⁺P⁻ 仍与项目采用的 M²=P⁺P⁻ 分开。

2. **S. G. Rajeev**，*Quantum Hadrondynamics in Two Dimensions*，
   [hep-th/9401115v1 HTML](https://arxiv.org/html/hep-th/9401115)，1994-01-24；
   [PDF](https://arxiv.org/pdf/hep-th/9401115)。核读 §4.1、§4.2、§4.4、§5，
   尤其 (83)–(95)、(106)–(108)、(133)–(142)。
   参考分裂 H=H_-⊕H_+ 给 epsilon，Grassmannian 约束是 (epsilon+M)²=1；
   线性化只保留相对于参考 sea 的非对角块。量子表示还要指定最低权态、
   色单态与固定 baryon sector，1/N_c 是其量子参数。
   本项目的有限 CAR 统一公式重新推导，不直接移植原文 M 的符号与因子。

3. **Jiuci Xu**，*Von Neumann Algebras in Double-Scaled SYK*，
   [2403.09021v7](https://arxiv.org/html/2403.09021v7)，2026-06-25。
   本轮核读引言的 algebra/state 区分，并沿已有 §3.3 笔记定位 Q-Fock。
   空 chord 向量在指定 field algebra 上给 trace 且 separating；把 chord number
   或独立 ladder 加入会改变所讨论的代数。不能由“Fock vacuum”这个名称推断它
   已是物理 canonical annihilation operators 的真空。

4. **Harshit Rajgadia–Jiuci Xu**，*Emergent States and Algebras from the Double-Scaling
   limit of Pure States in SYK*，
   [2604.14387v3](https://arxiv.org/html/2604.14387v3)，2026-09-29。
   本轮直接核 (33)–(48)、§3.2 (69)–(86)、§4.1–4.2 (87)–(110)、(123)。
   原模型是 Majorana SYK；适配 KM 自旋态的 strings 湮灭该态，通用 probes 则可只见到 trace。
   特别 (75) 中 bare matter b 为 canonical，实际适配算符却是
   M=q^(Delta n_tot/2)b；(85) 的次序不能改变。加入这些算符后原文构造不同的 GNS 代数。
   它提供了寻找 state-adapted observables 的具体范例，未证明本项目同 J、净电荷零
   particle–hole family 的 mixed(H,M) 律。该文 q=e^(−lambda) 与本项目 q=e^(−4lambda)
   使用不同 microscopic length 规范。

## 延用而不重复计算的靶标

[FLZ 0905.2280v2](https://arxiv.org/abs/0905.2280) 的 a=0 内部质量方程及谱，
保持 [bulk-target](../notes/susskind-bulk-target.md) 已核规范；本轮没有重算该谱。
参考态/极化、canonicalization、实际空间生成元与这个质量核是不同接口。
本轮也不把已有 `susskind-rindler-packet` 的二点质量歧义或四点非 Gaussian 检验
重新登记为新发现。

## 本轮对来源的使用方式

引用原文的是建构条件与已有表示。新 interface 中的 finite-CAR commutator、
tracial 下界、U(N)-Schur 恒等式、选定 Slater 态的同 J block 及其完整 H 首步
泄漏量均另行给出推导。Slater 单 ladder 的全 word 证明由本轮
`notes/susskind-slater-ladder-proof.md` 保存；不借用原文 Majorana chord rules
替代该推导。文献核查是针对这个接口的阅读，不作优先权或穷尽检索声明。
