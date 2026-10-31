# 创建算符下降审计：运动学有界映射、FP 缺口与条件自伴结论

后续更新：本稿的运动学构造继续作为输入；后续[全N生成空间](stage3-all-n-generated-space.md)及[受调节弱下降](stage3-regulated-local-descent.md)已在明确H_gen和共同群cutoff弱定义下实施原条件自伴步骤。未调节直接位置积分与最大分布空间的结论仍保持开放。

日期：2026-10-29。模型固定为 dS₂、L=1、G→0、自由实标量主系列 ν>0；σ>0 为 proper-time Gaussian 宽度，实际左右 observer 谱下界 s*>0。本文只新增此稿，不修改旧笔记、绑定计算或归档。

**审计结论：运动学 Tσ 的有界性及 SO⁺(2,1) 等变性可以从原 CG 合成定理、独立重建顶角和 Gamma 包络直接证明，无须借 FP 真空范数。下文补齐这一证明。尚不能把它直接称为原 (5.4) 在全部物理扇区上的创建块：完整 FP 等距、全部通道上的局部弱矩阵元兼容仍是独立前提。若这两步对全 N 成立，则有限粒子数域上的本质自伴性可由标准解析向量定理推出。**

本稿不声称重新证明了全部 CG 分解定理，也不把重建顶角称为作者确认勘误。

## 1. 原式、输入及非循环依赖

直接核对 [Kolchmeyer–Liu 2411.08090v1](https://arxiv.org/html/2411.08090v1) 的 (3.14)、(4.8)—(4.10)、(4.26)—(4.28)、(4.42)、(5.4)—(5.6)、(D.8)、(E.12)、(E.22)、(E.24) 及 E.1.3。[本地 HTML](../sources/kolchmeyer-liu-2411.08090v1.html) SHA256 为 `352ef29aa89c376203271bf984bd904fb471ff313a7d420bc6a0ebc875d78065`。

三个原式是

\[
\langle a,n_a|\phi_R|b,n_b\rangle
=\int d^2x\sqrt{-g}\,\overline{\Psi^a_{n_a}(x)}\,\varphi(x)\Psi^b_{n_b}(x),
\qquad
\varphi=\sum_n\bigl(\psi_n^\nu a_n+\overline{\psi_n^\nu}a_n^\dagger\bigr),
\tag{C1}
\]
\[
\phi_R(\tau)=e^{-iH_R\tau}\phi_Re^{iH_R\tau},\qquad
\langle a|\phi_R(\tau)|b\rangle
=e^{-i[m(a)-m(b)]\tau}\langle a|\phi_R|b\rangle,
\quad m(s)=\sqrt{s^2+1/4}.
\tag{C2}
\]

观察者 (4.8) 使用裸 ds，物质固定 ν 的单粒子基使用 Kronecker 内积。令 r(s)=s tanh(πs)，w(s)=πr(s)/4=π³ρ_Pl(s)。由 (4.9) 和 (D.8) 的共轭，创建模式核是

\[
t^{b n_b}_{a n_a;\nu n}
:=\int d^2x\sqrt{-g}\,\overline{\Psi^a_{n_a}}\,\overline{\psi_n^\nu}\,\Psi^b_{n_b}
=\sqrt{\frac{\pi}{r(\nu)}}\,
\overline{\mathcal V^b_{\nu,a,\mathrm{rec}}}\,
C^{b n_b,S}_{a n_a|\nu n}.
\tag{C3}
\]

这里 C 是原 E.24 的 untilded 系数。E.12 规定 C 本身是把输出表示嵌入两个输入表示的合成系数，因此 (C3) 中 C 不应再取共轭；反向的 CG 分解才使用其共轭。b 是创建前的实际右 observer 标签，同时是这一局部 CG 的总标签。a 是创建后的实际右 observer 标签。

独立输入按以下方向使用：

- [顶角解析重建](stage3-vertex-analytic-audit.md) 的 A9 给全部固定模式的三线性矩阵元，依赖原 E.1 完整两重度及明确反射/半相位约定；不依赖本稿或 FP 全扇区算符。
- 原 E.12、E.22、E.1.3 提供 CG 合成、正交及完备定理。下文将它当作原文表示论输入，不用真空二点重新标定 CG。
- [正核独立审查](../literature/stage3-observer-kernel-audit.md) §5 的 Gamma 包络给统一积分界。这是正函数的估计，不以全场算符存在为前提。
- [FP 路线](stage3-fp-multiplicity-bridge.md) 和[闭场域稿](stage3-closed-domain-next.md)仅列待完成目标，不能相互引用来证明全 FP 下降。

新发现的 [Marolf–Yu 2610.24073v1](https://arxiv.org/html/2610.24073v1) §3.3、(3.34)—(3.35) 独立说明同一 CG 合成/分解共轭方向，并使用 ρ_P=w；该文的 ρ_P 不是本项目的 ρ_Pl。它也明确 S 反射偶、A 反射奇。这提供额外来源支持，不自动替代本模型的 FP/局部核兼容证明。

## 2. Gaussian 相位与质量投影

以中心 τ₀ 的实正规化窗

\[
f_{\sigma,\tau_0}(\tau)=\frac1{\sqrt{2\pi}\sigma}
 e^{-(\tau-\tau_0)^2/(2\sigma^2)}
\]

涂抹 (C2)，得到

\[
h_{\sigma,\tau_0}(a,b)
=e^{-\sigma^2[m(a)-m(b)]^2/2}
 e^{-i[m(a)-m(b)]\tau_0}.
\tag{C4}
\]

每次插入的右投影 q_R=1[ŝ_R≥s*] 在输入 b 与输出 a 上都施加。下文直接把 a、b 积分限制于 [s*,∞)。这不是只限制外部包。左投影在后续物理分解中限制总标签 d≥s*；它与纯右创建对易。

(C4) 是先对固定模式的弱矩阵元作 Fourier 积分的结果。不能预设未涂抹全场已存在，再把它未经证明地写成有界算符值 Bochner 积分。零中心时 h 为实对称函数；非零中心时创建与湮灭必须按伴随取共轭相位，不能仍当作同一实核。

关键标签区别：N>0 时，右 observer 质量 b 与右 observer 加全部物质的总表示标签 d 通常不同。右质量 cutoff 只限制实际 observer 边；它不限制任意内部 CG principal 标签，更不删除内部 discrete 通道。

## 3. 用归一化 CG 合成严谨定义 Tσ

设

\[
H_R=\int_{s_*}^{\infty}{}^\oplus db\,H_b,
\qquad H_\nu=\ell^2(\mathbb Z),
\qquad \widehat C^{b n_b,S}_{a n_a|\nu n}
=\sqrt{w(b)}\,C^{b n_b,S}_{a n_a|\nu n}.
\]

原 E.22 在物理相位基中仍给

\[
\sum_{n_a,n}\overline{C^{b'n_b',S}_{a n_a|\nu n}}
 C^{b n_b,S}_{a n_a|\nu n}
=\frac{\delta(b-b')\delta_{n_b n_b'}}{w(b)}.
\tag{C5}
\]

它是谱分布正交式，不是固定 b 的 CG 列属于 ℓ² 的声明。不能在固定标签上写 δ(0) 再平方。

对固定 a，用原 CG 定理定义 S 通道合成等距映射

\[
J_a:\int_0^\infty{}^\oplus db\,H_b\longrightarrow H_a\otimes H_\nu,
\qquad
(J_a v)_{n_a,n}
=\int_0^\infty db\sum_{n_b}\widehat C^{b n_b,S}_{a n_a|\nu n}v(b,n_b).
\tag{C6}
\]

式中的积分先作谱合成理解；紧谱有限磁指标测试向量上的矩阵元由原系数给出，再按等距性延拓。完整 CG 分解还包括 A 与 D⁺/D⁻；J_a 只选其中 S 子空间，所以无需声称其像等于全部 H_a⊗Hν。

定义标量

\[
F_{\sigma,\tau_0}(a,b)
=\sqrt{\frac{\pi}{r(\nu)}}\,
\frac{\overline{\mathcal V^b_{\nu,a,\mathrm{rec}}}}{\sqrt{w(b)}}
 h_{\sigma,\tau_0}(a,b),
\qquad
(T_{\sigma,\tau_0}u)(a)=J_a\bigl[F_{\sigma,\tau_0}(a,\cdot)u\bigr].
\tag{C7}
\]

u 在 b<s* 处零延拓。CG 系数和 F 的可测性给可测合成场；先取紧 a 截止，再用下面的非负范数估计去截止。顶角 A9 的模平方直接给

\[
|F_{\sigma,\tau_0}(a,b)|^2
=\frac{4|\mathcal V^b_{\nu,a,\mathrm{rec}}|^2}{r(\nu)r(b)}
 e^{-\sigma^2[m(a)-m(b)]^2}
=V_\nu(a|b)e^{-\sigma^2[m(a)-m(b)]^2}.
\tag{C8}
\]

因此由 J_a 等距和非负 Tonelli 定理，

\[
\|T_{\sigma,\tau_0}u\|^2
=\int_{s_*}^{\infty}da\int_{s_*}^{\infty}db\,
 V_\nu(a|b)e^{-\sigma^2[m(a)-m(b)]^2}\|u(b)\|_{H_b}^2
=\int_{s_*}^{\infty}db\,A_\sigma(b)\|u(b)\|_{H_b}^2,
\tag{C9}
\]
\[
A_\sigma(b)=\int_{s_*}^{\infty}da\,
V_\nu(a|b)e^{-\sigma^2[m(a)-m(b)]^2},\qquad
\mathfrak M=\operatorname*{ess\,sup}_{b\ge s_*}A_\sigma(b)<\infty.
\]

统一界独立来自 Gamma 包络

\[
V_\nu(a|b)\le C_\nu\frac{e^{-2\pi(a-b)_+}}{1+|a-b|},
\]

及 |m(a)−m(b)|≥|a−b|/2（|a−b|≥1）。把积分分为 |a−b|≤1 与其补集即可获得与 b 无关的上界。这里没有用任何 FP 范数假设。

由 (C9) 及极化恒等式，T 唯一延拓为有界映射，且

\[
T_{\sigma,\tau_0}^*T_{\sigma,\tau_0}=A_\sigma(\hat s_R),
\qquad \|T_{\sigma,\tau_0}\|^2=\operatorname*{ess\,sup}A_\sigma\le\mathfrak M.
\tag{C10}
\]

若用 H_R 作自变量，右边准确写法是 Aσ(√(H_R²−1/4))；把 Aσ(H_R) 写作简记必须说明变量更换。

这补齐旧稿“由 bubble 可知 T 有界”的压缩论证。其成立依据是可用的 CG 谱合成定理及所列重建约定；不是把纯形式 CG 分布相乘视为已完成分析证明。

## 4. 等变性及 Schur 引理能给到哪里

J_a 根据 E.12、E.24 是 SO⁺(2,1) 的 intertwiner。F(a,b)、谱 cutoff 与 proper-time 相位只依赖表示标签，均与相应群作用对易，所以

\[
T_{\sigma,\tau_0}U_R(g)
=(U_R(g)\otimes U_\nu(g))T_{\sigma,\tau_0}.
\tag{C11}
\]

可先在谱合成测试向量上验证，再用 T 有界延拓到全部 H_R。无须把任意有限磁指标向量在有限群变换后的无穷展开塞进原时空积分；这避免了一个未证明的交换次序。

普通运动学空间上，取输出紧质量、有限 n_a,n 和输入紧质量、有限 n_b，(C7) 的弱矩阵元逐项正是 (C3)—(C4)。固定有限模式的三模积分在正标签紧集上有一致可积时间尾，故这里的有限求和、质量测试与积分交换可以控制。这已识别了运动学模式弱核。它仍不允许把任意物理 singlet 当成 H_L⊗H_R⊗Fock 中的普通向量。

Schur 的合法顺序是：先有有界等变映射，再有完整酉总表示分解，最后应用分解空间的交换代数定理。SO⁺(2,1) 的相应主系列以 d>0 不重复参数化；在明确的酉分解中，等变映射于同一 d 上形如 I_Hd⊗c(d)。它一般不是标量，因为 multiplicity 非平凡。

尤其 (C10) 的 Aσ(ŝ_R) 在加了旁观物质后不等于 Aσ(d)。不能把右 observer 的质量乘法算符，未经 CG 变换认同总 Casimir 的标量函数。

## 5. 完整运动学创建块及 Bose 投影

定义

\[
W_N=H_R\otimes\operatorname{Sym}^N H_\nu,
\qquad
C_N=\sqrt{N+1}\,(1_R\otimes P_{\operatorname{Sym}^{N+1}})
(T_{\sigma,\tau_0}\otimes1),
\tag{C12}
\]

其中把 T 新生成的物质因子与 N 个原因子排在一起。标准归一化对称张量约定下，分量式为

\[
(C_N\psi)(a,n_a;n_1,\ldots,n_{N+1})
=\frac1{\sqrt{N+1}}\sum_{j=1}^{N+1}\int db\sum_{n_b}
 h_{\sigma,\tau_0}(a,b)t^{b n_b}_{a n_a;\nu n_j}
 \psi(b,n_b;n_1,\ldots,\widehat n_j,\ldots,n_{N+1}).
\tag{C13}
\]

这是原实自由场创建作用的标准 Fock 组合因子；并非额外乘一个 (N+1)!。由于 P_sym 是正交投影且与群作用对易，

\[
\|C_N\|\le\sqrt{(N+1)\mathfrak M}.
\tag{C14}
\]

投影后通常只有这个上界，不能无证写 C_N* C_N=(N+1)Aσ(ŝ_R)。新增局部顶角只含 S，不意味着输入/输出完整总表示树只含 S；旁观粒子再耦合仍可涉及内部 A、D⁺、D⁻。

Bose 对称化在 FP 空间的收缩性还有一条独立证明：(4.26)—(4.28) 的 gauge-fixed 测度为正，P_sym 在每个 observer 坐标处只作用于物质 Hilbert 向量。因此在已经定义的 FP 测试波函数上逐点有 ‖P_sym Ψ(x_L,x_R)‖≤‖Ψ(x_L,x_R)‖，积分得 FP 收缩性；投影又保持群不变性。此事实不要求先完成全 CG/FP 酉桥。

反之，T 改变 observer 变量，不是这个逐点物质投影。上述点态不等式不能原封不动搬给 T 来绕过下降证明。

## 6. 从 C_N 到原 (5.4) 的物理算符：准确充分条件

设已构造完整运动学酉分解 G_N，主系列部分为

\[
G_N:W_N^P\simeq\int_0^\infty{}^\oplus dd\,(H_d\otimes\mathcal M_{N,d}),
\]

包含全部内部路径、连续测度、离散求和及 Bose 投影。总 discrete 不与左 principal 配对，但内部 discrete 不能删除。由 (C11)—(C14)，主系列之间的映射分解成

\[
G_{N+1}C_NG_N^{-1}
=\int_0^\infty{}^\oplus dd\,(1_{H_d}\otimes c_N(d)),
\qquad \|c_N(d)\|\le\sqrt{(N+1)\mathfrak M}\quad\text{几乎处处}.
\tag{C15}
\]

现在必须另证两项。

**FP 等距及完整性。** 把 g(d)∈M_N,d 与左 observer 的对偶 principal 配对，形成不变分布 R_Ng。采用原 singlet 的 (−1)^n 对偶相位，积分系数为相对于 dd 的 g；不要额外向这个不变系数乘一次 ρ_Pl。须在完整稠密测试空间证明

\[
Q_N^{FP}(R_Nh,R_Ng)
=\int_{d\ge s_*}\rho_{\rm Pl}(d)\langle h(d),g(d)\rangle\,dd,
\tag{C16}
\]

并证明其完成化给本模型的完整 N 粒子物理空间。已核 N=0 及 N=1 显式 S 通道不能替代 N=2 全通道证明。

**局部弱矩阵元兼容。** 对上述完整测试 g,h，由原 (5.4)、(3.14)、时间窗及逐次 Q 直接形成创建弱二次型 B_N^loc，证明

\[
B_N^{loc}(R_{N+1}h,R_Ng)
=\int_{d\ge s_*}\rho_{\rm Pl}(d)
 \langle h(d),c_N(d)g(d)\rangle\,dd.
\tag{C17}
\]

左边须是从原模式积分及 FP 规则独立定义、控制极限后得到的二次型；不能把右边用作左边的定义，然后宣称已识别原场。只需先对稠密测试 g,h 证明 (C17)：右边的统一界随即由 Riesz 定理给输出向量的存在和唯一性，并把创建块延拓到全部 H_N^phys。

在 (C16)—(C17) 成立后，设 U_N 为相应 FP 酉识别，便有

\[
C_N^{phys}=U_{N+1}^{-1}\left(\int^\oplus c_N(d)\,\rho_{\rm Pl}(d)dd\right)U_N,
\quad
\|C_N^{phys}\|\le\sqrt{(N+1)\mathfrak M}.
\tag{C18}
\]

ρ_Pl>0 对 d>0 成立，所以 dd 的零测集仍是物理测度的零测集；换权不增大 fiber 算符的本质上确界。这一步不需要在测度为零的单个 d 上定义物理态。

湮灭与 (C_N^phys)* 的识别也要在相同弱配对中核对。固定有限运动学矩阵元上，实场、实窗和自伴投影给复共轭配对；(C17) 及其共轭在完成化中给同一伴随。不能仅从运动学伴随关系，未经 (C16)—(C17)，跳到 FP 伴随关系。

## 7. 实际缺口及可执行的证明顺序

1. **(C16) 仍须独立证明。** 原 (4.33) 的 δ(0)Σ_n/V 替换或抽象 (4.50) 不能替代全通道规约积分。新的 group averaging 正性/Schur 来源有帮助，但 GA 谱种子系数与不变分布系数需明确换权；GA 与指定 FP 积分的等同也须证明。
2. **有限磁指标 cutoff 不保持 singlet。** 可以用它核查普通运动学矩阵元；不能把截断后的非不变态直接套用只为不变函数建立的 FP 轨道公式。应在不变 CG 系数上做谱/通道截止，或独立证明所用规约极限的收敛与截面无关性。
3. **Tonelli 的范围有限。** 它直接适用于 (C9) 的非负范数积分；不直接允许交换原复三模积分、无穷角动量和、全部 CG 积分及群规约。可用先合成后取 Hilbert 极限、受控弱配对来避免过强的绝对收敛要求。
4. **域的像通常不保留有限通道和紧谱。** C_N 有界只给其像属于整个 N+1 扇区。Gaussian 在质量差上衰减，不会把紧输入支撑变成紧输出支撑；Bose 再耦合也可生成无限通道。
5. **反射与时间反演分开。** 本稿先在 proper SO⁺(2,1) 复 Hilbert 空间工作。空间反射为酉约束，在一致模式约定下可另施加；D⁺/D⁻成对交换。反线性时间反演的“实子空间”不能直接当作复 Hilbert 子空间套用同一自伴定理。

建议顺序是先证 N=2 的完整 (C16)，再独立检查 N=0→1、N=1→2 的 (C17)，保留 A 与内部 discrete 输入测试。仅证明到 H₂ 并不包含 C₂:H₂→H₃ 的物理下降，更不足以宣布全 N 场已定义。

## 8. 全 N 下降成立后的条件本质自伴命题

本节仅为条件命题，不能反过来证明 (C16)—(C17)。假设对每个 N 已有同一模型的完整物理空间、上述创建块及其物理伴随，并满足 (C18)。定义

\[
\mathcal H^{phys}=\bigoplus_{N\ge0}\mathcal H_N^{phys},\qquad
\mathcal D_{fin}=\bigoplus_{N\ge0}^{alg}\mathcal H_N^{phys},
\quad \Phi=C+C^*\text{ 于 }\mathcal D_{fin}.
\]

D_fin 包含每个扇区的所有 Hilbert 向量；它不是只取有限 CG 通道和光滑紧谱系数。因每个块有界，D_fin 稠密、不变，Φ 在其上对称。

若 ψ 的粒子数支撑不超过 N，则

\[
\|C\psi\|\le\sqrt{(N+1)\mathfrak M}\|\psi\|,
\qquad
\|C^*\psi\|\le\sqrt{N\mathfrak M}\|\psi\|,
\]
\[
\|\Phi^k\psi\|\le
(2\sqrt{\mathfrak M})^k
\sqrt{\frac{(N+k)!}{N!}}\,\|\psi\|.
\tag{C19}
\]

所以任意 t≥0，级数 Σ_k t^k‖Φ^kψ‖/k! 的控制项相邻比为

\[
\frac{2t\sqrt{\mathfrak M}\sqrt{N+k+1}}{k+1}\longrightarrow0.
\]

D_fin 中每个向量均为整解析向量。由 Nelson 解析向量定理，Φ|D_fin 本质自伴。定理来源为 [Nelson, Analytic Vectors (1959)](https://doi.org/10.2307/1970331)；已缓存的 [Neeb 2011](../sources/neeb-2011-aif2660.pdf) Theorem 4.8 给所用准确表述。这是标准算符论推论，不主张新的自伴判据。

若进一步已证明每个完整 H_N 中的有限通道紧谱测试空间 D_CG,N 稠密，并且其局部场作用与上述块相同，那么固定扇区的界

\[
\|\Phi\xi_N\|^2
\le(2N+1)\mathfrak M\|\xi_N\|^2
\]

使任意 ξ_N 可由 D_CG,N 在图范数中逼近。有限扇区求和后，D_CG=⊕alg D_CG,N 是 D_fin 的图范数稠密子域，因而也是自伴闭包的核心。这个结论不要求 D_CG 在 Φ 下不变，但要求它对每个完整扇区确实稠密且作用已匹配。

当前状态是：(C3)—(C15) 的运动学构造已在明确源输入下闭合；(C16)—(C17) 的完整物理识别仍由独立 FP/GA 研究承担；本节自伴性保持条件性质。没有把 N=2 进展、真空四点或数值检查升级成全 N 定理。
