# Oriented transfer：本轮原文与继承接口

核读日期：2026-11-07。继续 Susskind static-patch/common-coupling complex DSSYK 主线。本轮的具体问题是实际相邻电荷传播缺陷，而不是另开一个独立 probe 模型。

1. [Berkooz–Narovlansky–Raj, Complex Sachdev-Ye-Kitaev model in the double scaling limit, arXiv:2006.13983v2](https://arxiv.org/html/2006.13983v2)。再次核读 §3.2 的 (3.10)–(3.12)、§3.3 的 (3.13)–(3.15) 和 (4.3)。可直接继承的是有限 oriented CAR 的符号/局部迹，以及固定有限个 p-set 的 pair-Poisson 组织。原文 CAR 为2，本项目为1；原文能量方差不可直接移植到 finite driver 的整数规范。新版带源 singleton 分母、whole-field D mark 和内部弦电荷背景是本项目在该框架内的推导，不归为原文已有公式。

2. [前驱 charged transfer](../notes/susskind-memory-charged-transfer.md)：共同原 J 的 T_m、polar map、forward/reverse defects；独立 J′ 仅在 Gaussian OU 条件积分表示中出现。此次 full-Fock D 恒等式自动保留两方向，不把 H_m 与 H_(m−1) 独立化。

3. [前驱 uniform bound](../notes/susskind-memory-uniform-bound.md)：Gaussian Sobolev 左右乘法/实时间传播、number-weighted polar defect、rank-filter 估计。全时间 E→G 的新误差是纯 number multiplier；真正的时间导数在该估计内处理，未从无导数误差直接求导。

4. [前驱 size bridge](../notes/susskind-memory-size-bridge.md) 与 [neutral descendants](../notes/susskind-memory-wick-descendants.md)：pure Wick source 的 half-size local mark、source 禁内部配对、whole-field chaos bookkeeping、切图回归及HS判据。rank/size桥不能自动变成 body-loss桥。neutral 无桥相消也不能直接覆盖 charged ordinary field 的内部收缩。

5. [前驱 two-component](../notes/susskind-memory-two-component.md)：实际 memory=A+B+o(L²)，A₀→Y₋、A₁→R+aJ_HY₋，以及严格非零 ν_(A₀,B₂)=d(w+η)μ_(q²)。本轮 G₂在扣 energy 与charge两层背景后恢复旧w，下一步仍须恢复全部mixed谱。

本轮未进行新宽泛文献扫查；先复用已核原文和既有严格工具，随后建立所缺的局部源接口。未将现有 q-Hermite、Poisson 或 Gaussian Sobolev 工具包装为原创数学定理。issue #1属于旧KS支线，已完成的条件算术复核状态不变，未替代当前记忆核任务。
