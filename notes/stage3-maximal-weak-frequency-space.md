# 弱正频条件下的最大有限 FP 空间

2026-10-29。模型保持 dS₂、L=1、proper SO⁺(2,1)、G→0，自由主系列物质 ν>0。固定物质粒子数 N，K_N=Sym^N P_ν。位置空间 Y=X×X、六轨道截面 Z、正测度 dμ 与前稿相同。本文不把集中于退化轨道上的奇异分布纳入 L²(Z;K_N)。

**问题与结论。** 原 KL (4.40) 用正频观察者模式写非正规 singlet，(4.42) 用位置 FP 范数选态，但没有指定其无穷展开的分布测试拓扑。本文给一个独立于群平均生成闭包的弱谱条件，并证明：在这个明确条件下，全部有限 FP 位置态恰等于先前构造的 H_gen。没有先把“正频”定义为 H_gen，也不声称任意未经指定的形式求和自动满足此条件。

证明输入：[完整环境谱测试核心](stage3-ambient-spectral-test-core.md)、[位置 FP–GA 定理](stage3-fp-two-particle-inner-product.md) F4—F7、[真实共同域](stage3-positive-frequency-orbit-core.md) P16—P18、[全 N 生成空间](stage3-all-n-generated-space.md)及[实际产品幅度与稠密性](stage3-product-seed-ga-density.md)。原单观察者正频子空间来自 [KL 2411.08090v1](https://arxiv.org/html/2411.08090v1) (4.9)—(4.16)、脚注25和39。以下是这些输入的应用；未作首创性声明。

## 1. 每个有限 FP 态的非紧测试配对是确定的

暂不加 observer 质量下限。记环境规约空间

\[
\mathcal F_N=L^2(Z,d\mu;K_N),\qquad
(Rk)(g s(z))=V_N(g)k(z).
\tag{W1}
\]

既有 P18 给 R:F_N→L¹_loc(Y;K_N) 连续、注入。因此每个 k 确实定义位置分布，不能对任意 L¹ 函数直接取切片来代替这一构造。

对强轨道域中的测试 φ，记

\[
B\phi(z)=\int_G V_N(g^{-1})\phi(g s(z))dg,
\qquad a_\phi(z)=\int_G\|\phi(g s(z))\|dg.
\]

内积反线性于第一项。coarea、Bochner–Fubini 及 Cauchy–Schwarz 给

\[
T_k(\phi):=\int_Y\langle Rk(y),\phi(y)\rangle dy
=\langle k,B\phi\rangle_{\mathcal F_N},\qquad
|T_k(\phi)|\le\|k\|_2\|a_\phi\|_2.
\tag{W2}
\]

绝对可积性由 ∫_Y‖Rk‖‖φ‖=∫_Z‖k‖a_φ 保证。它对所有下述谱波包、束缚态产品适用。这个非紧测试配对不是任意选出的分布延拓：取任意光滑位置 cutoff 0≤χ_j≤1、χ_j→1，则 a_{(1−χ_j)φ}→0 在 L²(Z) 中成立，先逐轨道支配收敛，再由 a_φ² 支配。故 T_k(φ)=lim_j T_k(χ_jφ)，与 cutoff 选择无关。

## 2. 完整观察者测试，而不只是正频测试

在 H_X=L²(X,dx) 中令 H_+ 为原 KL 正频模式的闭子空间，H_⊥=H_+^⊥。两者均在 proper 群作用下不变；这是单观察者运动学表示中的正交分解，与所要证明的物理空间无关。

取 S_+ 为所有有限磁指标、f∈Cc∞((0,∞)) 的原正频包 ∫f(s)Ψ_n^s ds 的有限线性包。取 S_⊥ 为负频连续谱包与全部固定磁指标离散径向束缚态的有限线性包。完整环境稿说明：每个固定 n 的连续谱重数为二；KL 正负频模式覆盖它，另外恰有 |n| 个束缚态。正负频波包的正规化与交叉正交由[LaTeX 稿附录 A](../manuscript/normalization.tex)以自伴径向时间演化和出射波包极限独立核对，不依赖形式 δ/PV 交换。S_+、S_⊥ 是两个运动学子空间中的具体谱测试类。

对 ε,δ∈{+,⊥} 定义

\[
\mathcal E_{\epsilon\delta}
=\operatorname{span}\{u_L(x_L)u_R(x_R)v:
 u_L\in\mathcal S_\epsilon,\ u_R\in\mathcal S_\delta,\ v\in K_N\}.
\tag{W3}
\]

完整环境稿使用 1/4<α<1 的位置矩估计，证明每个 E 都包含于 D_orb；而且它们合在一起的 B 像在整个 F_N 中稠密。此处关键是强轨道范数下的谱逼近，不能用普通 L² 运动学稠密取代。

不同 E sector 的运动学交叉矩阵系数对每个 g 都严格为零，因为至少一个 observer 因子属于相互正交、群不变的子空间。因此 FP–GA 定理给

\[
\langle B\phi,B\psi\rangle_{\mathcal F_N}=0
\quad (\phi\in E_{\epsilon\delta},\ \psi\in E_{\epsilon'\delta'},\
 (\epsilon,\delta)\ne(\epsilon',\delta')).
\tag{W4}
\]

定义 F_{εδ}=closure B E_{εδ}。由总性和 W4，得到实际位置 FP 空间中的正交分解

\[
\boxed{\mathcal F_N=\mathcal F_{++}\oplus\mathcal F_{+\perp}
\oplus\mathcal F_{\perp+}\oplus\mathcal F_{\perp\perp}.}
\tag{W5}
\]

这也构造出范数≤1的物理正交投影 Π_{++}。在完整谱测试类上，B(P_+⊗P_+⊗1)φ=Π_{++}Bφ。这里没有把运动学投影取到群正则表示的 trivial fiber，也没有断言 P_+ 保持任意 D_orb 元素。

## 3. 独立定义弱正频，并识别最大空间

**弱正频条件：** k∈F_N 满足

\[
T_k(\phi)=0\quad\text{对全部 }
\phi\in E_{+\perp}+E_{\perp+}+E_{\perp\perp}.
\tag{W6}
\]

这是用单观察者的完整补谱检验位置分布；每个配对由 W2 的绝对积分定义，不涉及 H_gen。允许所有 v∈K_N，也包括物质内部离散路径，不只局部创建可达的通道。

由 W2、W5，W6 当且仅当 k∈F_{++}。另一方面，[全 N 稿](stage3-all-n-generated-space.md)的实际产品幅度、全磁指标延拓和密度证明可令正 observer 支撑遍历 (0,∞) 中的紧集。所有估计只在各紧集使用，故无需统一控制 s→0，即得到

\[
\mathcal F_{++}=\mathcal H_N^{gen,0}
\simeq\mathcal H_N^{spec,0}.
\tag{W7}
\]

上标0表示没有正的 observer 下限，不表示物质真空。对 N≥1，

\[
\mathcal H_N^{spec,0}
=\int_{d>0}^{\oplus}\rho_{Pl}(d)dd
  \int_{a>0}^{\oplus}da\,\mathcal M_{N,d,a},\qquad
\rho_{Pl}(d)=\frac{d\tanh\pi d}{4\pi^2}.
\tag{W8}
\]

对 N=0，仅有 g(d)∈L²(ρ_Pl dd)，实际右 observer 标签 a=d，不能另造独立 da。上述扩展只是态空间定理，不宣称 Gaussian 创建块在撤掉下限后仍有此前的统一界。

## 4. 在物理谱识别之后施加真实 observer 下限

固定 s_*>0。令 S_{+,<} 是 f∈Cc∞((0,s_*)) 的正频包。除了 W6，再要求

\[
T_k(u_Lu_Rv)=0
\tag{W9}
\]

对 u_L∈S_{+,<},u_R∈S_+，以及 u_L∈S_+,u_R∈S_{+,<} 的全部产品成立。这是左右实际 observer 的谱支撑条件；并不对物质或 CG 内部谱施加 s_*。

已由 W7 把 k 写为 J_N^0 g。取左磁指标0、右总磁指标0的原粒子基 e_i(a)，其 CG 系数为 u_i(a,d)。实际产品幅度仍是

\[
g_{i;f,k_R}(d,a)=K_P(d)f(d)k_R(a)u_i(a,d),
\qquad K_P=\rho_{Pl}^{-1}.
\tag{W10}
\]

对 W9 测试，ρ_Pl K_P=1，故配对是 ∫f(d)k_R(a)〈g(d,a),u_i(a,d)〉da dd。先取 f 支撑于 (0,s_*)、k_R 任意紧正谱包。用既有产品密度证明的局部 L¹ 界、乘积测试总性和固定 a 的完整 CG 逆变换，推出 g=0 在 d<s_*。反过来取 k_R 支撑于 (0,s_*)，同理推出 g=0 在 a<s_*。这里须对可数 i 取共同零集，再在紧 A×D 上截断逆变换；不能把单个 d 的广义 CG 列贸然叫作可数正规基。

反向充分性也须覆盖 W9 的全部磁指标测试，而非只有上述小域。完整有限磁指标幅度是 g_test=K_PΣ_m(−1)^m f_m(d)v^P_{−m}(a,d)；左包质量支撑直接成为 d 支撑，右实际质量因子 k_R(a) 被逐 a 的 CG 变换保留。因此任一左/右低谱产品的物理幅度位于对应低谱区域，与高谱 g 严格正交。

因此 W9 恰等价于 g 支撑在 a,d≥s_*。边界 a=s_* 或 d=s_* 在这一 L² 测度下为零，不会新增 δ 支撑的有限范数态。N=0 时由 g_{f,k_R}=K_P f k_R 直接得到 d≥s_*，同样没有独立 a 条件。

由固定下限版本的全 N 定理，此支撑子空间的 J_N^0 像恰是原来构造的 H_N^gen(s_*)。所以如果定义

\[
\mathcal H_N^{max,w}(s_*):=
\{k\in\mathcal F_N:\text{满足W6和W9}\},
\]

便有带实际位置映射的结论

\[
\boxed{\mathcal H_N^{max,w}(s_*)
=\mathcal H_N^{gen}(s_*)=J_N\mathcal H_N^{spec}(s_*).}
\tag{W11}
\]

W6、W9是连续测试泛函的核，故左侧本身闭。硬阈 s_* 是在完整物理 L² 谱识别以后施加的；没有声称把环境紧支撑测试在 s_* 附近切出一个窄洞会在强轨道范数中趋近原测试。

## 5. 这个结论怎样解释原 (4.40)

原 (4.40) 指定正频模式以及物质 Hilbert 值，(4.42) 指定有限 FP 范数，脚注31还排除某些物质和不收敛的例子。原文没有规定一般非正规系数 A 的模式和在何种测试拓扑中收敛。因此本文给出的是一个明确的最大弱谱实现，而不是替原文证明所有可能的形式求和约定都相同。

若某个 (4.40) 展开在下述意义成立：对完整谱测试 E 的 W2 配对可逐项/受控极限求得，且不存在补频和下限外系数，那么它自动满足 W6、W9，W11 即保证无额外有限 FP 扇区。反过来 H_gen 中每个元素都满足这些弱谱约束，且通过明确 J_N 与完整 multiplicity 系数表示。这比只宣布“按生成闭包定义物理空间”多证明了一个覆盖方向。

另有一个容易检查的充分极限条件：若 k_j∈H_gen(s_*)、sup_j‖k_j‖_FP<∞，且 Rk_j 在局部分布意义下收敛到 F，那么 Hilbert 弱紧性给子列 k_j⇀k；R 连续到 D' 保证 F=Rk。H_gen 是闭线性子空间，故弱闭，k∈H_gen。注入性保证 k 唯一，且 ‖k‖≤liminf‖k_j‖（对实现下极限的子列）。任意弱聚点都因 R 的注入性等于同一 k，故整列弱收敛到 k。若另外范数收敛到 ‖k‖，则弱收敛提升为强收敛。

这里只把“局部分布收敛＋统一 FP 界”列为充分条件；不把未带任何范数控制的形式模式和说成已获相同结论，也没有声称所有有趣的分布态必然具有有限 FP 范数。

## 6. 对后续局部场的作用与边界

[Gaussian 弱场](stage3-regulated-local-descent.md)已在固定 s_*>0 的 H_gen 上构造，具有明确的共同群 cutoff、逐次 observer 质量投影和本质自伴性。W11允许将它称为“W6、W9所定义的最大有限 FP 空间上的该受调节弱场”，而无需假设多出未知的有限 FP 正频扇区。

但这不证明未调节位置乘法/时空积分等于该弱作用；也不证明所有形式展开、其他 α-vacuum 选择、反线性时间反演规约或有限 G 模型具有相同空间。旧四点与这个弱场的具体有序乘积匹配仍需逐式完成。完整 de Sitter 全息对偶没有由本定理得到。
