# 从 KM 文献到守电荷 Slater ladder 的微观接口

2026-11-07。本稿给出本轮标准 CAR、Gaussian Wick 与组合计数的解析推演；独立有限 CAR 检查和完整交错 word 审计由同阶段其他记录承担。本稿未调用任何科学 main，不修改旧封存输出，不作优先权声明。

## 1. 输入与文献接口

采用 \(N\) 个 complex fermion 的标准 CAR、\(H=H_s=\sum_{A,B}J_{AB}C_A^\dagger C_B\)，\(|A|=|B|=p\ge1\)，\(J_{BA}=\overline{J_{AB}}\)，协方差
\[
\mathbf E J_{AB}J_{CD}=w_N\delta_{AD}\delta_{BC}.
\]
完整 Hamiltonian 仍为旧同一样本 \(J\)，没有新增独立 Hamiltonian。本文把方差 \(w_N\) 纳入 \(C_m,C_R\)。算符范数以下均明确指参考态与 Gaussian disorder 的 annealed 向量范数，非旧 full-Fock trace，也非逐样本归一。

选固定占据集合 \(O\)，\(m=|O|\)，\(a=N-m\)，\(p\le\min(m,a)\)；\(\Omega_O=|O\rangle\) 为其 Slater 态。先保持 \(N,m,p\) 有限，再取 \(p^2/N\to\lambda>0\)、\(m/N\to1/2\)，每个 word 长度固定。

文献给出的可借路线见[本轮来源](../literature/susskind-bilocal-polarization-sources.md)。Rajgadia–Xu [2604.14387v3](https://arxiv.org/html/2604.14387v3) 在与 KM 态相关的算符族上得到
\[
M^\dagger\mapsto b^\dagger q^{\Delta\hat n_{\mathrm{tot}}/2},
\qquad M\mapsto q^{\Delta\hat n_{\mathrm{tot}}/2}b.
\]
其 PDF (3.53)–(3.54)，HTML (85)–(86)。原文 \(H\) 不守构造 KM 真空所用的 spin-flip 数；我们不能把物理 complex 空态直接当成有非平凡中性动力学的 KM 态。本稿用半填充 Slater 态引入粒子—空穴极化。

先保留一个有限 \(N\) 的精确对照：Gaussian \(J\) ensemble 对 \(U(N)\) 协变，而 \(\mathcal H_m=\wedge^m\mathbb C^N\) 不可约，故 Schur 引理给
\[
\mathbf E\langle O|f(H)|O\rangle=\mathbf E\,\tau_m(f(H))
\tag{S0}
\]
对每个有限多项式 \(f\) 成立，其中 \(\tau_m=\binom Nm^{-1}\operatorname{Tr}_{\mathcal H_m}\)。证明是 \(\mathbf Ef(H)|_{\mathcal H_m}\) 与 \(U(N)\) 的表示对易，因而为 scalar。选择 Slater 态本身不会改变 \(H\)-only annealed 谱；极化信息必须由与 \(O\) 相关的联合 probe 保留。(S0) 是标准不变性与 Schur 引理的应用。

## 2. 原 Hamiltonian 的纯激发块

定义
\[
R_O^\dagger=
\sum_{\substack{A\subset O^c,\ B\subset O\\|A|=|B|=p}}
J_{AB}C_A^\dagger C_B,\qquad
R_O=(R_O^\dagger)^\dagger.
\tag{S1}
\]
它将 \(p\) 个 occupied fermion 移入 empty sites，仍守物理 \(Q\)。\(R_O\Omega_O=0\)，而 \(R_O^\dagger\Omega_O\) 的各标签对应互不相同的占据态。令
\[
L=\binom ap\binom mp,\quad C_R=w_NL,\quad
B_O^\dagger=\frac{R_O^\dagger}{\sqrt{C_R}},
\quad B_O=(B_O^\dagger)^\dagger.
\tag{S2}
\]
故
\[
\mathbf E\|B_O^\dagger\Omega_O\|^2=1,\qquad
\mathbf E\|B_O\Omega_O\|^2=0.
\tag{S3}
\]
这已与同一 full-Fock trace 下两个范数相等的情况不同；改变的是态及保留的态适应微观算符。

## 3. 完全 creation word 的计数

对固定 \(n\ge0\)，取越界 binomial 为零。范数的计数推演为
\[
s_{n,N}:=\mathbf E\| (B_O^\dagger)^n\Omega_O\|^2
=n!\prod_{r=0}^{n-1}
\frac{\binom{a-rp}{p}\binom{m-rp}{p}}{L}.
\tag{S4}
\]
其有限代数理由如下：

1. 连续 creation 中所有 \(A\) 必须两两不交，所有 \(B\) 也必须两两不交；否则重复创建或重复湮灭使该项为零。已选 \(r\) 个标签后，下一标签数为 \(\binom{a-rp}{p}\binom{m-rp}{p}\)。
2. 非零项的 \(2p\)-fermion 算符支撑不交，皆为偶算符，所以换序无负号；作用于 \(\Omega_O\) 的向量范数为 1。
3. Raising family 中没有一对互为 Hermitian conjugate 的标签；其系数是独立 proper complex Gaussian。bra 与 ket 的 \(n!\) 个 Wick 匹配贡献相同。乘上 \(w_N^n\) 后与 (S2) 的归一相消。

由于不同 \(n\) 具有不同粒子—空穴激发数 \(np\)，creation 向量在有限 \(N\) 已正交。对 \(n\) 固定，
\[
\log\frac{\binom{a-rp}{p}}{\binom ap}
=-\frac{rp^2}{a}+o(1),
\qquad
\log\frac{\binom{m-rp}{p}}{\binom mp}
=-\frac{rp^2}{m}+o(1).
\]
因此
\[
s_{n,N}\longrightarrow n!e^{-2\lambda n(n-1)}
=n!q^{n(n-1)/2},\qquad q=e^{-4\lambda}.
\tag{S5}
\]
有限 \(N\) 末端因 Pauli 排斥而截断；(S5) 没有声称对增长的 \(n\) 一致。

单个 canonical oscillator 的 \(D^\dagger=b^\dagger q^{\hat n/2}\) 满足同一 creation 范数。仅由 (S4)–(S5) 可确定 creation span 上的 raising 作用及其压缩 adjoint；完整的纯 ladder 结论还需控制实际 annihilation。这个步骤已由[本轮解析证明 §5，(8)、(12)](susskind-slater-ladder-proof.md)完成：全部固定 Dyck-word 矩趋于 dressed oscillator 的真空矩，且 annihilation 离开固定 creation level 的范数趋零。有限 \(N\) 并不精确闭合，backward Wick 配对可给非零修正；双标度结论依赖固定 word 长度。§6另给 scalar self-averaging，不能混同为对全部算符的统一范数控制。这里完成的是纯 \(B_O,B_O^\dagger\) ladder；加入同一 \(J\) 的完整 \(H\) 后如何闭合仍待建立。

## 4. \(H\Omega_O\) 中保留了多少纯激发

令 \(\Pi_k^{(O)}\) 投影到相对 \(O\) 有 \(k\) 个新粒子及 \(k\) 个空穴的占据态子空间。完整 \(H\) 先从 \(O\) 湮灭 \(p\) 个 fermion，再把 \(p\) 个 fermion 创建到刚移出的 \(p\) 个位置或原来的 \(a\) 个 empty sites。因此
\[
C_m=\mathbf E\|H\Omega_O\|^2
=w_N\binom mp\binom{a+p}{p},
\tag{S6}
\]
\[
\mathbf E\|\Pi_k^{(O)}H\Omega_O\|^2
=w_N\binom mp\binom ak\binom p{p-k}.
\tag{S7}
\]
这与项目已用的固定 charge covariance 公式相同；这里的新增用途是把其作用于一个 Slater 参考态的激发数分开。恒等式
\[
\Pi_p^{(O)}H\Omega_O=R_O^\dagger\Omega_O
\tag{S8}
\]
在每个样本严格成立。

特别是纯激发**范数平方权重**
\[
\gamma_N:=\frac{C_R}{C_m}
=\frac{\binom ap}{\binom{a+p}{p}}
=\prod_{j=0}^{p-1}\frac{a-j}{a+p-j}
\longrightarrow e^{-2\lambda}=\sqrt q.
\tag{S9}
\]
如果定义 annealed 单位向量
\[
|\mathsf B\rangle=\frac{R_O^\dagger\Omega_O}{\sqrt{C_R}},
\qquad
|\mathsf H\rangle=\frac{H\Omega_O}{\sqrt{C_m}},
\]
则 overlap **振幅**是
\[
\langle\mathsf B,\mathsf H\rangle_{\rm ann}
=\sqrt{\gamma_N}\longrightarrow e^{-\lambda}=q^{1/4}.
\tag{S10}
\]
所以 \(\sqrt q\) 是权重，不是两个单位向量的 overlap 振幅。正交余项的范数平方趋于 \(1-\sqrt q\)，固定正 \(\lambda\) 时不能忽略。

## 5. Return deficit 是一个精确超几何分布

由 (S7) 得权重
\[
\mathbb P_N(k)
=\frac{\binom ak\binom p{p-k}}{\binom{a+p}{p}},
\qquad r=p-k,
\]
\[
\mathbb P_N(r)
=\frac{\binom p r\binom a{p-r}}{\binom{a+p}{p}}.
\tag{S11}
\]
这是从 \(a+p\) 个位置中抽 \(p\) 个、其中 \(p\) 个标记为“刚被移出位置”的超几何分布。它是 annealed 的正交范数权重，不是宣称每个单样本具有这些精确测量概率。

其下降阶乘矩为
\[
\mathbf E_N(r)_j=\frac{(p)_j^2}{(a+p)_j}
\longrightarrow(2\lambda)^j,
\tag{S12}
\]
而对 \(t\ge0\)，
\[
\mathbf E_N(1+t)^r
=\sum_{j=0}^p\frac{(p)_j^2}{(a+p)_j}\frac{t^j}{j!}
\le\exp\!\left(\frac{p^2t}{a+1}\right).
\tag{S13}
\]
有限正 \(\lambda\) 的该界控制指数矩，因此 (S12) 给
\[
r\ \Rightarrow\ {\rm Poisson}(2\lambda).
\tag{S14}
\]
固定 \(r\) 的扇区在 \(N\to\infty\) 前彼此正交。尽管 \(k/p=1-r/p\to1\)，不能仅凭这一归一化计数收敛就把所有 return 通道当成同一个向量；(S9) 已量化纯 \(r=0\) 通道丢失的有限权重。

## 6. 下一步要证明的物理连接

本接口给出了有明确 microscopic support、真空 annihilation 条件和可核 counting 的新算符族。接下来的决定性问题是完整同一 \(J\) 的混合矩：\(B_O,H\)、return-resolved \(H_k=\Pi_k^{(O)}H\Omega_O\) 的状态制备，以及它们与已重建 bilocal 的交叉关联。此处 \(H_k\) 只是向量标签，并非已定义的闭合算符。

这里的 \(\Omega_O\) 一般不是 \(H\) 的本征态，也没有证明是 stationary/KMS 态；\(B_O\) 不是无参考数据的 \(U(N)\) singlet。Haar 平均参考态而保持 probe 固定会恢复扇区迹。要讨论协变性，需要同时处理参考极化与算符。完整 MSS 的 singlet/meson 时间、空间及质量字典仍需从这些共同矩建立，不能将 \(p\)-激发、chord number 或 return deficit 直接命名为介子质量。

本轮优先级是证明/否定 ladder 与真实 Hamiltonian 的共同代数，再决定是否保留或如何补足这组适应变量。体光前 \(P^\pm\)、Rindler KMS 制备与规范 dressing 是下一层问题。
