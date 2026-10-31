# 全 CG 创建顶角与物理定义域：A.7 归一化审查

日期：2026-10-29。来源为 Kolchmeyer–Liu，[arXiv:2411.08090v1](https://arxiv.org/html/2411.08090v1)，重点核对 (4.8)—(4.9)、(5.4)、(5.8)—(5.10)、(6.17)、(6.23)—(6.26)、(A.7)、(A.14)、(D.8)、(E.14)、(E.22)、(E.24)。本稿是新的源公式核验记录，不更改此前代码、结果或[闭平方域笔记](stage3-closed-domain-next.md)。

**结论：原版本 A.7 的 literal 读法与 D.8 所需的输入交换对称性不相容；它经 E.22 收缩所得的标量核，也与正文核相差非恒定正因子。PDF 印刷页 76 已视觉确认与 HTML 一致。当前不能用 A.7 闭合“正文谱权 → 全 CG 顶角 → 物理最小域”的证明。** 下文的假设修正只用来定位问题，不当作原文或作者确认的勘误。

## 1. 源公式与取向

观察者基底按 (4.8) 用 ds 归一化；(4.9) 给出 Ψ⁽ˢ⁾ₙ=√[s tanh(πs)/π] ψ⁽ˢ⁾ₙ。因此创建一个物质主系列粒子 ν>0，右观察者由 b 变为 a 的矩阵元，应从 (5.4) 的创建项与 D.8 的复共轭取得：

\[
T^{b n_b}_{a n_a;\nu n}
=\int d^2x\sqrt{-g}\,(\Psi^{a}_{n_a})^*(\psi^\nu_n)^*\Psi^b_{n_b}
=\sqrt{\frac{\pi}{\nu\tanh\pi\nu}}\,
\overline{\mathcal V^{b}_{\nu,a}}\,
\bigl(\mathcal C_{PP}^{P,S}\bigr)^{b n_b}_{a n_a|\nu n}.
\tag{I1}
\]

这里 b 是 CG 的输出标签、物理创建操作的输入观察者标签；两者方向不同。不能仅凭图的箭头交换 a,b。E.24 的纯相位改写不改变 E.22 的谱正交常数：

\[
\sum_{n_a,n}\overline{C^{b'n_b'}_{a n_a|\nu n}}
 C^{b n_b}_{a n_a|\nu n}
=\frac{4}{\pi b\tanh\pi b}\,\delta(b-b')\delta_{n_b n_b'}.
\tag{I2}
\]

这些式子理解为正标签紧集上的测试函数分布恒等式，随后才讨论完备化。

## 2. 不依赖 bubble 归一化的对称性检验

D.8 左边的两个未共轭波函数可交换。E.14 的 S 通道 CG 在两个输入表示同时交换时对称；E.24 的两输入相位乘积也对称。因此 D.8 所使用的顶角系数必须满足

\[
\mathcal V^{s_1}_{s_2,s_3}=\mathcal V^{s_1}_{s_3,s_2}.
\tag{I3}
\]

但 A.7 的已刊表达式包含 **Γ(1/2+i s₂) 的倒数**与**Γ(1/2+i s₃) 的乘积**；其余四个四分之一实部 Gamma、平方根与指数因子均对 s₂↔s₃ 对称。逐字采用该公式得到

\[
\frac{\mathcal V^{s_1}_{s_2,s_3}\big|_{\rm literal}}
     {\mathcal V^{s_1}_{s_3,s_2}\big|_{\rm literal}}
=\left[\frac{\Gamma(1/2+i s_3)}{\Gamma(1/2+i s_2)}\right]^2.
\tag{I4}
\]

一般 s₂≠s₃ 时该比值不为 1，其模为 cosh(πs₂)/cosh(πs₃)。这甚至不是可由单位模相位消除的差异。独立审查者另行核查 D.8、E.14、E.19、E.24 后确认了这一交换对称问题。

PDF 与 HTML 采用相同排版，故不能把疑点归于 HTML 转换。证据保存在 [原 PDF](../results/stage3/cg-domain-source-audit/kolchmeyer-liu-2411.08090v1.pdf) 和 [印刷页 76 截图](../results/stage3/cg-domain-source-audit/v1-page76.png)。截图由 Poppler 以 160 dpi 渲染并经视觉检查。

## 3. Literal A.7 的模平方与 kinematic bubble

定义正实函数

\[
\Gamma_8(s_1,s_2,s_3)=
\prod_{\epsilon_1,\epsilon_2,\epsilon_3=\pm1}
\Gamma\!\left(\frac14+\frac i2(\epsilon_1s_1+\epsilon_2s_2+\epsilon_3s_3)\right).
\]

由 |Γ(1/2+is)|²=π/cosh(πs) 与 1+e²ˣ=2eˣ cosh x，可直接化简 A.7：

\[
|\mathcal V^{s_1}_{s_2,s_3}|^2_{\rm literal}
=\frac{\cosh^2\pi s_1}{64\pi^3\cosh^2\pi s_3}
\left[\prod_{j=1}^3s_j\tanh\pi s_j\right]\Gamma_8(s_1,s_2,s_3).
\tag{I5}
\]

为使常数可复查，指数、平方根与 cosh 的模平方合成 cosh(πs₁)/[64π² cosh(πs₂)cosh(πs₃)]；三个半整数 Gamma 的模平方合成 cosh(πs₁)cosh(πs₂)/[π cosh(πs₃)]。两者相乘即 (I5) 前因子。

将 (I1) 作为由 literal A.7 定义的候选核并用 (I2) 收缩，得到

\[
K_{\rm literal}(a|b)
=\frac{4|\mathcal V^b_{\nu,a}|^2_{\rm literal}}
 {\nu\tanh\pi\nu\,b\tanh\pi b}
=\frac{a\tanh\pi a\,\cosh^2\pi b}{16\pi^3\cosh^2\pi a}
 \Gamma_8(\nu,a,b).
\tag{I6}
\]

而此前由 (6.17)/(6.25) 使用的正文核为

\[
V_\nu(a|b)=\frac{a\tanh\pi a\,\cosh^2\pi b}{16\pi^5}\Gamma_8(\nu,a,b),
\qquad
\boxed{K_{\rm literal}(a|b)=\frac{\pi^2}{\cosh^2\pi a}V_\nu(a|b).}
\tag{I7}
\]

这个差异依赖最终质量 a；它不能用一次整体归一化去掉。主文核也可由 (6.25) 的 HH 二点谱权除以初态 HH 密度取得。

**仅用于诊断的反事实检验：**如果把 A.7 最后一项 Γ(1/2+i s₃) 从分子移到分母，则 (I3) 恢复，(I5) 去掉 π²/cosh²(πs₃)，(I6) 恰等于 Vν。但没有原始积分的独立计算或作者勘误，这不能升级为已确认的修正公式。

## 4. 可保留的数学小步，以及物理解释的边界

固定 σ>0、观察者谱下界 s*>0，令

\[
h_\sigma(a,b)=e^{-\sigma^2[m(a)-m(b)]^2/2},\quad
m(s)=\sqrt{s^2+1/4},\quad
\mathfrak M=\sup_{b\ge s_*}\int_{s_*}^\infty V_\nu(a|b)|h_\sigma(a,b)|^2da<\infty.
\]

对 compact-mass、有限 n 的运动学观察者测试态 u，把 (I1) 乘 hσ 并积分 b 定义 Tσ,literal。先使用 CG 的谱正交性，再作非负 a 积分，得到

\[
\|T_{\sigma,\rm literal}u\|^2
=\int_{s_*}^\infty db\sum_{n_b}|u(b,n_b)|^2
 \int_{s_*}^\infty da\,|h_\sigma(a,b)|^2K_{\rm literal}(a|b)
\le\pi^2\mathfrak M\|u\|^2.
\tag{I8}
\]

因而 literal A.7 的 CG 核本身定义了一个有界的运动学创建映射，范数≤π√𝔐。结论只涉及按该公式定义的映射。由于 (I3)—(I7) 冲突，当前不能说这个映射就是原 (5.4) 的物理场创建块。

即便日后从独立三模积分确定了正确归一化，(I8) 型运动学估计仍不能直接当作完整物理 Hilbert 空间的算符域证明。物理 singlet 是分布式向量，原文采用 Faddeev–Popov 内积；它们不是运动学 Hilbert 空间的普通闭子空间。需要另外在物理稠密测试域上证明：

1. 与原 (5.4) 相同的创建/湮灭核对所有 relevant principal、discrete 与 multiplicity 通道有定义，并下降到规约等价类。
2. 创建范数估计和湮灭的伴随关系在 FP 物理内积中成立，尤其是控制规约极限、谱积分、CG 求和之间的交换。
3. 该测试域确实在完整 H_Q 稠密，而不仅在真空顶角生成的通道稠密。

若这些条件得到证明，才可补上[闭平方域笔记](stage3-closed-domain-next.md)中全域最小对称算符的存在条件；还不能因此直接宣称四次方域或本质自伴性。

## 5. 对既有计算的影响与版本核查

只读检查确认：`quantized_observer.py` 的正谱权直接实现正文 (6.17)/(6.25)；`observer_crossing.py` 与 `observer_crossing_regular.py` 直接实现 (A.14)、(B.21)—(B.24)。它们没有由 A.7 重建顶角。因此本次发现不要求默改这些代码，也不改变它们以已采用核为对象的数值核验与分析估计。

另一方面，当前源公式存在内部不一致，故既有“显式核的数学性质”不能被扩大解释为独立确认了 A.7、D.8 与全部物理场矩阵元的一致性。本次暂停的是这条新尝试的全 CG 归一化/物理域识别证明。

2026-10-29 访问 [arXiv 摘要与 submission history](https://arxiv.org/abs/2411.08090)，页面只列 2024-11-12 的 v1。搜索 `"2411.08090" erratum correction vertex`、`"Chaos and the Emergence of the Cosmological Horizon" "A.7"`、`"Kolchmeyer" "Liu" "Cosmological Horizon" corrigendum` 未找到相关勘误；[MIT CTP 刊目](https://physics.mit.edu/ctp-papers/2411-08090/)没有提供额外版本信息。这是本次检索未找到，不是证明不存在勘误。没有联系作者。

下一步最小可检验任务是从 (A.5) 的模式函数与 (E.24) 的 CG 独立计算一个非零三模积分，并比较相位和模；仅检验交换对称性已足够否定 literal A.7 与 D.8 同时成立，但不足以证明唯一正确修正。源归一化澄清后再推进 FP 下降，避免把计算投入到不相容的起点。