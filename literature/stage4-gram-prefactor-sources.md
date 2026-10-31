# Gram 与横向能量推导：原文核查

核读日期：2026-10-30。本记录承接[WKB来源](stage4-wkb-overlap-sources.md)及[TQFT字典](tqft-recoil-dictionary.md)。新增计算见[中文推导](../notes/stage4-gram-prefactor.md)。以下区分原文公式、我们从公式做出的代数推导和未获得的输入。

## CLP

J. Chandra, D. Liska, B. Post，*3D de Sitter from complex Virasoro TQFT: Observers and knots*，[v1原文](https://arxiv.org/html/2610.30247v1)。[本地缓存](../sources/chandra-liska-post-2610.30247v1.html)，SHA-256：`fe8580884e81ba9817e292bbca7c70b7c7d494e89efef62695e94bcf07801d0f`。

- (3.43)—(3.47)：完整四面体图为双手征6j除以四个顶角C的平方根，再取实Euclidean模平方。它还不是已证明等同于CSTY规范归一关联函数的对象。
- (3.57)：顶角领先指数抵消；没有给出所有order-one因子抵消。
- (3.64)后讨论：手征一圈因子含(det G)^(−1/4)，并有顶角C的一圈修正。
- (3.35)：物理两点含N₁₂;W、两个Z_obs和|C₀|²；不能只拿裸6j作物理四点归一化。
- (3.37)：固定探针数据下 N₁₂;W 的领先项是 2π²/[κ(α₁α₂)^(3/2)] exp(−πν)|Γ(1+iν)|²，加O(κ⁻²)。这不自动给联合增长探针极限的一致误差。
- 脚注12：外部iR变换说明只适用于约化作用量；内部6j轮廓不等于物理外部能量轮廓。
- 附录B：q₁²−4q₀q₂=16 det G用于核对面指标排列，包括复角样本。

精确Gram因式分解及其两种极限是本轮代数计算。把几何因子用作退化构型的整个6j一圈项是额外假设，不能归给原文已经证明。

## CSTY

Y. Chen, D. Stanford, H. Tang, Z. Yang，*Negative shocks versus static patch holography*，[v1原文](https://arxiv.org/html/2607.14042v1)。[本地缓存](../sources/chen-stanford-tang-yang-2607.14042v1.html)，SHA-256：`9e352212021abf7135852f4b5dc8cd933f47334bd72e78f8def9d84d34fed5e2`。

- (2.3)：固定平均质量的混合表示保留四段能量并做delta约束，仍有三个能量方向。
- (2.73)：反极点时间选择。与(2.3)合用产生本轮ω相位和u、v边缘积分；Jacobian为4。
- 脚注19：明确需要对正交于ω的能量积分。此前仅比较对称截面的动作量不足以确定振幅。
- (2.54)、(B.5)—(B.7)：纯反冲四面体经典作用量及对偶平面四面体角公式，是本轮等质量Hessian的输入。
- §2.6：重背景、小折角、纯反冲约化函数。前一轮相对合流窗口仍只针对这个约化函数。

本轮Hessian由原六角作用量导出，没有当作新的有限G动力学。其局部高斯积分要进入物理关联函数，还须证明适当轮廓及一致余项。

## Meng–Yang

Y. Meng, T. Yang，*Asymptotics of complex b-6j symbols*，[v1原文](https://arxiv.org/html/2606.03653v1)。本轮在线核读Theorem 1.3，未缓存；来源主清单仍为26份文件。

定理针对固定hyperideal双曲二面角及规定arg b范围，给含Gram行列式的渐近和O(b²)余项。未将其解释成退化球面构型、联合增长质量、外部复杂能量路径上的一致定理。这些适用条件正是后续使用该一圈公式时要补足的部分。
