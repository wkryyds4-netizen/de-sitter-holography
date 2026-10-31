# 两粒子完整 CG 树：测度、离散路径与 Bose 投影

日期：2026-10-29。本文只闭合运动学表示分解的字典，供 [FP multiplicity 路线](../notes/stage3-fp-multiplicity-bridge.md) 的两粒子计算使用。最后写出的 FP 范数是待直接核验的目标，不是本稿已经证明的物理等距式。没有启动新数值扫描，也不修改既有论文缓存。

## 1. 来源、范围与记号

主要来源是 Kolchmeyer–Liu [2411.08090v1 附录 E](https://arxiv.org/html/2411.08090v1#A5)，具体使用 E.3–E.30、E.1.3、E.32–E.39；相位见 A.5–A.6 及 E.24–E.25。本地 [HTML](../sources/kolchmeyer-liu-2411.08090v1.html) SHA-256 为：

~~~text
352ef29aa89c376203271bf984bd904fb471ff313a7d420bc6a0ebc875d78065
~~~

独立数学交叉来源是 Groenevelt [math/0501511v1](https://arxiv.org/html/math/0501511v1)，§3.2 的 Proposition 3.6–3.7 与 §3.3：其酉谱分解独立支持 principal×discrete 的连续重数一与同号离散系列；§3.3 末段述及 principal×principal 的连续重数二及离散系列，并明确不展开该乘积的计算；因此它只交叉支持分解结构，不独立验证 KL 的 PP 系数与归一化。PD 分解对应取 k=Δ∈ℤ>0、ε=0，输出离散标签覆盖全部正整数。该文采用不同基和 Gamma 权，本稿的具体常数全部来自 KL，未直接移植其基相位或 Wilson 核。其显式 Racah 核讨论至少两条离散外腿，不能当作本题三个 principal 外腿的现成核。

模型是 dS₂、半经典引力极限 G→0，固定物质质量参数 ν>0；m(s)=√(s²+1/4)，曲率半径为 1。先处理 proper SO⁺(2,1)，磁指标为整数。用 P_s 表示主系列，Δ_s=1/2+is，s>0；D_Δ^± 表示 Δ=1,2,… 的正、负离散系列。这里不改变群的覆盖，也不引入半整数磁指标或 complementary 系列。

右观察者空间采用

\[
H_R=\int_0^\infty{}^\oplus da\,P_a,\qquad
W_2=H_R\otimes\operatorname{Sym}^2P_\nu.
\tag{C1}
\]

a 是右观察者谱标签，d 是右观察者与两物质合成的总 principal 标签，q 是两物质的内部 principal 标签。连续谱一律用裸 Lebesgue 测度，将权吸收入 CG；离散谱用计数测度。观察者下界 a≥s_* 可在最后限制，不能据此删去内部 q 或 Δ。

附录 E 的 tilded CG 是“耦合态在乘积基中的展开系数”；其复共轭才用于将乘积态分析到谱系数。下文先全部使用该基，最后统一转成正文的 untilded 模式相位。

## 2. 四类顶角及各自归一化

| 输入→输出 | raw 系数 | raw 完备关系中的权 | 裸谱测度的 normalized 系数 |
|---|---|---|---|
| P×P→P，x=S,A | \(\widetilde C_x\) | \(w(s)=\pi s\tanh(\pi s)/4\) | \(\widehat C_x=\sqrt{w(s)}\,\widetilde C_x\) |
| P×P→D_Δ^± | \(\widetilde D^\pm\) | 1 | \(\widehat D^\pm=\widetilde D^\pm\) |
| P×D_Δ^±→P | \(\widetilde B^\pm\) | \(v_\Delta(s)=s\tanh(\pi s)/(\Delta-\frac12)\) | \(\widehat B^\pm=\sqrt{v_\Delta(s)}\,\widetilde B^\pm\) |
| P×D_Δ^±→D_κ^± | \(\widetilde E^\pm\) | 1 | \(\widehat E^\pm=\widetilde E^\pm\) |

特别是 v_Δ 不等于 w。把两者都写成 √w 会给内部离散路径引入错误常数。只有

\[
w(s)=\pi^3\rho_{\rm Pl}(s),\qquad
\rho_{\rm Pl}(s)=\frac{s\tanh\pi s}{4\pi^2}.
\tag{C2}
\]

### 2.1 P×P→P：精确定义、两重数与完备关系

令 θ_ij=θ_i−θ_j，并按 E.15–E.16 定义

\[
I_S^{\delta_1\delta_2\delta_3}
=\left|\sin\frac{\theta_{31}}2\right|^{-\delta_3-\delta_1+\delta_2}
 \left|\sin\frac{\theta_{32}}2\right|^{-\delta_3-\delta_2+\delta_1}
 \left|\sin\frac{\theta_{21}}2\right|^{-\delta_2-\delta_1+\delta_3},
\]
\[
I_A=I_S\,\operatorname{sgn}\!\left[
\sin\frac{\theta_{31}}2\sin\frac{\theta_{32}}2\sin\frac{\theta_{21}}2
\right].
\tag{C3}
\]

E.17–E.20 等价于下面的完整 Fourier 定义，固定了 2π 常数：

\[
(\widetilde C_x)_{s_1,n_1|s_2,n_2}^{s,n}
=\frac1{(2\pi)^3}\int_{[0,2\pi]^3}d^3\theta\,
 I_x^{1-\Delta_{s_1},\,1-\Delta_{s_2},\,\Delta_s}
 e^{-in_1\theta_1-in_2\theta_2+in\theta_3}.
\tag{C4}
\]

共同平移使其只在 n=n₁+n₂ 时非零。这里实标签给各两点奇性实部 −1/2，三点重合的二维相对坐标奇性总次 −3/2，故 Fourier 积分局部可积；不需要将三点重合误当非可积极点。需要推广参数时再声明解析延拓。

E.21–E.22 的自正交常数为

\[
\sum_{n_1+n_2=n}
 \overline{(\widetilde C_x)^{s',n}_{s_1,n_1|s_2,n_2}}\,
 (\widetilde C_x)^{s,n}_{s_1,n_1|s_2,n_2}
 =\frac{\delta(s-s')}{w(s)},\qquad x=S,A.
\tag{C5}
\]

这里不能只由两个自范数推断 S/A 交叉正交。可联合 E.1.3 的完整谱分解使用；也可在 Fourier 规范下补核：同时翻转所有磁指标时 C_S 不变、C_A 变号，在总磁指标 n=0 上交叉内积相消；该重数间的交织内积由表示作用传到其余 n。所有连续谱正交式都先在谱波包上理解。

E.14 给

\[
(\widetilde C_S)_{s_2,n_2|s_1,n_1}^{s,n}
=(\widetilde C_S)_{s_1,n_1|s_2,n_2}^{s,n},\qquad
(\widetilde C_A)_{s_2,n_2|s_1,n_1}^{s,n}
=-(\widetilde C_A)_{s_1,n_1|s_2,n_2}^{s,n}.
\tag{C6}
\]

只有输入确实是相同粒子时，这个交换才是物理 Bose 交换。

### 2.2 P×P→D⁺/D⁻：磁指标、权与交换奇偶

用 k≥Δ 表示正整数。原 E.26–E.27 的记号是

\[
|k\rangle_\Delta^+
=\sum_{n_1+n_2=k}
(\widetilde D^+)^{\Delta,k}_{s_1,n_1|s_2,n_2}|n_1,n_2\rangle,
\]
\[
|-k\rangle_\Delta^-
=\sum_{n_1+n_2=-k}
(\widetilde D^-)^{\Delta,k}_{s_1,n_1|s_2,n_2}|n_1,n_2\rangle.
\tag{C7}
\]

负离散的 CG 上标仍是正 k；不可将其上标 k 偷换成 −k。精确 3F2 见下文附录 A，转录自 E.30 后的 HTML 组 A5.EGx54、A5.EGx55。

E.28–E.30 给

\[
\sum_{n_1+n_2=\epsilon k}
\overline{(\widetilde D^\epsilon)^{\Delta',k}_{s_1,n_1|s_2,n_2}}
(\widetilde D^\epsilon)^{\Delta,k}_{s_1,n_1|s_2,n_2}
=\delta_{\Delta\Delta'},\qquad \epsilon=\pm,
\]
\[
(\widetilde D^\epsilon)^{\Delta,k}_{s_2,n_2|s_1,n_1}
=(-1)^\Delta
(\widetilde D^\epsilon)^{\Delta,k}_{s_1,n_1|s_2,n_2}.
\tag{C8}
\]

因此该离散顶角已经归一化，不再乘 √w 或 √v。

为把完整性写清楚，固定总磁指标 N，记 i=(n₁,n₂)、j=(n₁′,n₂′)，两对之和均为 N。E.1.3 的三个情形合为

\[
\sum_{x=S,A}\int_0^\infty ds\,
 \overline{\widehat C_x^{s,N}(j)}\widehat C_x^{s,N}(i)
+\mathbf1_{N>0}\sum_{\Delta=1}^{N}
 \overline{\widetilde D_+^{\Delta,N}(j)}
 \widetilde D_+^{\Delta,N}(i)
+\mathbf1_{N<0}\sum_{\Delta=1}^{-N}
 \overline{\widetilde D_-^{\Delta,-N}(j)}
 \widetilde D_-^{\Delta,-N}(i)
=\delta_{ij}.
\tag{C9}
\]

原位置分别为 A5.EGx56、A5.EGx57、A5.EGx58。N=0 没有离散项，不能将这一磁切片的事实推广成整棵树没有离散表示。

### 2.3 P×D→P 和同号 D：另一种连续谱权

正离散支使用 E.32–E.33：

\[
(\widetilde B^+)^{s,n}_{a,m|\Delta,k}
=(-1)^k(\widetilde D^-)^{\Delta,k}_{a,m|-s,-n},
\quad k\ge\Delta,\quad m+k=n.
\tag{C10}
\]

右端的 −s 是源公式的 principal 对偶/等价参数，必须保留；它不是负观察者质量。直接替换成 +s 需要另写等价算符及相位。右端在特殊参数上的表观 Gamma 或 3F2 奇点须对完整表达式取解析极限，不能逐因子截断。

同号离散输出记

\[
(\widetilde E^+)^{\kappa,n}_{a,m|\Delta,k}
:=(\widetilde{\mathcal C}_{P,+}^{+})^{\kappa,n}_{a,m|\Delta,k},
\qquad n=m+k\ge\kappa\ge1,
\tag{C11}
\]

其精确公式是原 E.35，见附录 B。E.37 对 n≤0 只有连续项；E.38–E.39 对 n≥1 还含所有 κ=1,…,n 的离散项。合写为

\[
\int_0^\infty ds\,v_\Delta(s)\,
 \overline{(\widetilde B^+)^{s,n}_{a,n-k'|\Delta,k'}}
 (\widetilde B^+)^{s,n}_{a,n-k|\Delta,k}
+\mathbf1_{n\ge1}\sum_{\kappa=1}^{n}
 \overline{(\widetilde E^+)^{\kappa,n}_{a,n-k'|\Delta,k'}}
 (\widetilde E^+)^{\kappa,n}_{a,n-k|\Delta,k}
=\delta_{kk'}.
\tag{C12}
\]

相应谱正交写为 \(\widehat B^\dagger\widehat B=\delta(s-s')\)、\(\widehat E^\dagger\widehat E=\delta_{\kappa\kappa'}\)，不同谱支正交。式 C12 与完整 irreducible 分解一起固定了此规范；不能遗漏离散项后把连续部分单独当成满空间酉变换。

负离散情况原脚注 77 只说可由复共轭得到，未展开其基相位。本稿补选一个可逐生成元检查的**线性反射规范**，不声称已逐相位匹配脚注的复共轭写法：

\[
R_P|n\rangle_s=|-n\rangle_s,\qquad
R_D|k\rangle_\Delta^+=|-k\rangle_\Delta^-,
\]
\[
R:L_0\mapsto-L_0,\quad L_1\leftrightarrow L_{-1}.
\tag{C13}
\]

由 E.3–E.11，principal 有
\(R_PL_1R_P|n\rangle=-i(-n-1+\Delta_s)|n+1\rangle
=L_{-1}|n\rangle\)；离散的两组根式也逐项对应。故以下定义是合法交织规范：

\[
(\widetilde B^-)^{s,-n}_{a,-m|\Delta,-k}
:=(\widetilde B^+)^{s,n}_{a,m|\Delta,k},
\]
\[
(\widetilde E^-)^{\kappa,-n}_{a,-m|\Delta,-k}
:=(\widetilde E^+)^{\kappa,n}_{a,m|\Delta,k}.
\tag{C14}
\]

C14 的 D 输入/输出磁指标使用有符号记法，C7 的 \(\widetilde D^\pm\) 上标仍用正 k。负支权仍为 v_Δ(s) 与单位计数权；其离散输出在 n≤−1 时 κ≤−n。C13 已由独立数学核查逐生成元复核。

## 3. 完整两阶段树，然后才选总 principal

先在未 Bose 化的张量积上做完整分解：

\[
P_\nu\otimes P_\nu
\simeq\int_0^\infty{}^\oplus dq\,(P_q^S\oplus P_q^A)
\ \oplus\ \bigoplus_{\Delta\ge1}(D_\Delta^+\oplus D_\Delta^-).
\tag{C15}
\]

每条内部 principal 路径再分解

\[
P_a\otimes P_q
\simeq\int_0^\infty{}^\oplus dd\,(P_d^S\oplus P_d^A)
\ \oplus\ \bigoplus_{\kappa\ge1}(D_\kappa^+\oplus D_\kappa^-).
\tag{C16}
\]

每条内部 discrete 路径再分解

\[
P_a\otimes D_\Delta^\epsilon
\simeq\int_0^\infty{}^\oplus dd\,P_d
\ \oplus\ \bigoplus_{\kappa\ge1}D_\kappa^\epsilon.
\tag{C17}
\]

连续顶角分别采用 √w(q)、√w(d)、√v_Δ(d)，所有离散顶角用单位权。C15–C17 包含最终总离散支，是完整运动学树。

下一步与左 principal 观察者配成群不变态时，最终总离散不能与左 principal 配对，因此只保留 C16/C17 的总 principal 部分。C17 的内部 D_Δ^ε→总 P_d 路径依然存在；丢弃它们会破坏两粒子空间的完整性。

## 4. 合成 normalized 核与相位

把三条外腿按 (observer, matter 1, matter 2) 排列，磁指标为 (m,n₁,n₂)，总 n=m+n₁+n₂。总 principal 的两种合成核是

\[
T_{q,\beta,\alpha}^{d,n}(m,n_1,n_2)
=\sum_{k\in\mathbb Z}
(\widehat C_\alpha)^{d,n}_{a,m|q,k}
(\widehat C_\beta)^{q,k}_{\nu,n_1|\nu,n_2},
\quad \alpha,\beta\in\{S,A\},
\tag{C18}
\]
\[
T_{\Delta,\epsilon}^{d,n}(m,n_1,n_2)
=\sum_{k\ge\Delta}
(\widehat B^\epsilon)^{d,n}_{a,m|\Delta,\epsilon k}
(\widetilde D^\epsilon)^{\Delta,k}_{\nu,n_1|\nu,n_2}.
\tag{C19}
\]

磁守恒使 C18 中只有 k=n₁+n₂ 的项；C19 只有 εk=n₁+n₂ 的项。C19 中 B 的 discrete 输入是有符号指标，D 的输出上标 k 是正指标。

对固定 a，合成核在谱波包意义下满足

\[
\sum_{m,n_1,n_2}\overline{T_{\lambda'}^{d',n'}(m,n_1,n_2)}
T_\lambda^{d,n}(m,n_1,n_2)
=\delta(d-d')\delta_{nn'}\delta_{\lambda\lambda'}.
\tag{C20}
\]

λ=(q,β,α) 时最后一项为 δ(q−q′)δ_ββ′δ_αα′；λ=(Δ,ε) 时是 Kronecker δ；连续/离散内部路径的交叉项为零。证明顺序是先对两物质磁指标用 C5/C8/C9 的正交关系，再对外层用相应顶角谱正交。完整乘积空间的恒等分解还需要 C16/C17 的最终离散支；仅 C18/C19 给总 principal 投影。

正文 untilded 模式与这里的 E 基之间按 E.24–E.25 变换。取

\[
p_{s,n}=e^{i\phi_{n,s}/2},\qquad
e^{i\phi_{n,s}}=
\frac{i\cosh(\pi s)\Gamma(\frac12+is)
\Gamma(\frac12-is+n)\Gamma(\frac12-is-n)}
{\pi\Gamma(\frac12-is)}.
\tag{C21}
\]

半相位分支与原模式 A.5 及 (3.11) 一致，取 p_s,n=p_s,−n 和 p_s,0=e^{iπ/4}。则每个 PP→P 顶角乘 p_out/(p_in,1 p_in,2)。离散基可选 p_D=1，并在相邻两顶角使用同一选择。两顶角内部相位完全相消，整棵树仅余

\[
T_{\rm untilded}
=\frac{p_{d,n}}{p_{a,m}p_{\nu,n_1}p_{\nu,n_2}}\,
T_{\rm tilded}.
\tag{C22}
\]

C22 对内部 principal 和 discrete 路径相同。对每个顶角独立重选离散基相位却不补偿相邻顶角，会改变树系数，不是无害记号替换。C14 选定的外层负支与附录 A 原 D⁻ 内顶角可一致组合；此时空间反射交换整棵 D⁺/D⁻ 树的相位需由内顶角追踪，不预设它恒为 +1。

## 5. Bose 投影在物质优先树中对角

定义 \(\tau_{12}\) 只交换两物质外腿，\(P_B=(1+\tau_{12})/2\)。它与 proper 群作用对易。由 C6/C8：

\[
P_BT_{q,\beta,\alpha}=\mathbf1_{\beta=S}T_{q,\beta,\alpha},
\qquad
P_BT_{\Delta,\epsilon}
=\mathbf1_{\Delta\ {\rm even}}T_{\Delta,\epsilon}.
\tag{C23}
\]

因此内部 principal 只剩 β=S，内部 discrete 只剩 Δ=2,4,…，正负支均保留。外层 α=S/A **都保留**：外顶角交换的是 observer 和已耦合的物质对，不是两个相同物质。局域三模顶角只激发某个 S 通道，也不能用来删去一般两粒子态的 outer A。

固定 d,a 后的 Bose multiplicity 空间因此是

\[
\mathcal M_{2;d,a}
\simeq L^2(\mathbb R_+,dq;\mathbb C^2)
\oplus\ell^2(2\mathbb Z_{>0}\times\{+,-\}),
\]
\[
\mathcal M_{2;d}
=\int_0^\infty{}^\oplus da\,\mathcal M_{2;d,a}.
\tag{C24}
\]

未 Bose 化时第一项为 \(\mathbb C^4\)，第二项有全部正整数 Δ。

采用标准对称张量与 Fock 的单位映射：对 u_n₁n₂=u_n₂n₁，

\[
u\ \longmapsto\
\frac1{\sqrt{2!}}\sum_{n_1,n_2}u_{n_1n_2}
a^\dagger_{n_1}a^\dagger_{n_2}|0\rangle.
\tag{C25}
\]

其范数为 \(\sum_{n_1,n_2}|u_{n_1n_2}|^2\)。C23 中已经位于对称子空间的 normalized 树无需再乘任意 √2。后续从原 (4.38)–(4.42) 识别场态时必须显式核对这一约定，不能凭粒子记号猜测阶乘。

## 6. 与 observer-first 树的关系

observer-first 树先分解 P_a⊗P_ν₁：内部 principal 有 S/A 两重数，内部 D_Δ^± 有所有 Δ≥1；再接 P_ν₂。内部 principal→总 principal 又有 S/A 两重数；内部 D^±→总 principal 用 C10/C14 的顶角并显式交换张量因子到 P⊗D 顺序。相同方式保留最终总 discrete，先构成完整酉分解。

记两个完整分析变换为 U_mat、U_obs，则

\[
\mathcal R=U_{\rm mat}U_{\rm obs}^\dagger
\tag{C26}
\]

是酉 Racah 重耦合。固定 a、总 principal d 时，令 λ 为物质优先标签，η 为 observer-first 标签，均包含连续四重数和离散双支；其作用是

\[
(\mathcal R f)_\lambda
=\sum_{u,v=S,A}\int_0^\infty db\,
 \mathcal R_{\lambda;(b,u,v)}f_{b,u,v}
+\sum_{\Delta=1}^\infty\sum_{\epsilon=\pm}
 \mathcal R_{\lambda;(\Delta,\epsilon)}f_{\Delta,\epsilon}.
\tag{C27}
\]

该核可定义为两套合成磁核的重叠，但应先用总谱及内部谱波包规避 δ(0)，再解读其分布核。这里给出精确的酉变换定义与完整求和域，未声称算出了闭式三 principal Racah 矩阵。

Bose 投影在 observer-first multiplicity 上为

\[
P_B^{\rm obs}=\mathcal R^\dagger P_B^{\rm mat}\mathcal R.
\tag{C28}
\]

它一般混合连续/离散路径及其多重度。因此不能把 observer-first 树逐项投影后仍当成原来的正交基，也不能只留下“每顶角 S”来实现 Bose 对称化。原 A.14 用于特定四点顶角的单个 S 型 6j，并不等于 C27 的完整重数矩阵。

## 7. 给 FP 推导的目标及严格边界

采用 C18/C19 的 normalized 顶角及 C25 的 Fock 约定，先取 d,a,q 的紧支撑光滑谱波包（远离 0）、有限离散 Δ 支、有限通道组合。针对这些测试系数，待从原 (4.42) 直接核验的是

\[
\boxed{
\|\Psi_{g,h}\|_{\rm FP}^{\,2}\ \stackrel{?}{=}\
\int_0^\infty \rho_{\rm Pl}(d)\,dd\int_0^\infty da
\left[
 \sum_{\alpha=S,A}\int_0^\infty dq\,|g_\alpha(d,a,q)|^2
 +\sum_{\substack{\Delta=2,4,\ldots\\\epsilon=\pm}}
 |h_\epsilon(d,a,\Delta)|^2
\right].
}
\tag{C29}
\]

这是 normalized 运动学树经 N=0 与 N=1 原 (4.44)–(4.45) 显式 S 通道校准提出的 FP 目标。未在此处证明所有两粒子通道的群规约、整体常数、交叉项消失、积分极限交换及零空间识别。也没有把形式无限群体积相消当作证明。

若坚持以 raw 顶角乘积展开，并明确在同样裸 dd da dq 测度下写系数 F_α、F_ε，则顶角变换本身给

\[
g_\alpha(d,a,q)
=\frac{F_\alpha(d,a,q)}{\sqrt{w(d)w(q)}},\qquad
h_\epsilon(d,a,\Delta)
=\frac{F_\epsilon(d,a,\Delta)}{\sqrt{v_\Delta(d)}}.
\tag{C30}
\]

这只是同一广义核展开的系数重标度。若改为带权谱积分，必须同时重定义系数，不可混用 C30。目标 C29 中不再额外乘 w(q)、v_Δ(d) 或离散维数因子。

完成 C29 后仍须独立处理原 (5.4) 局部创建算符的矩阵元：检查它在以上所有测试通道上与运动学等变创建映射一致，并控制谱积分、磁求和、Bose 投影和 FP 规约的次序。C29 本身即使成立，也不单独证明全 N 的场闭包或本质自伴性。

空间反射应在 proper 群物理空间构成后另做投影；它交换 D⁺/D⁻，并依赖已选相位。时间反演的反线性结构另行处理。两者均不与 C23 的 Bose 交换混同。

## 8. 源文中须保留的疑点与本稿验收范围

- 原 A.7 的顶角归一化疑点及脚注 33 与显式 E.16 反射奇性的张力，见 [解析审查](../notes/stage3-vertex-analytic-audit.md)。本稿保留显式 E.15/E.16，不利用脚注把 A 通道删去；原缓存不改，也不称作者确认勘误。
- HTML (4.43) 左侧仍含 n₁，而右侧把 n₁,n₂ 求和并留下 n。这里仅记录指标标签需澄清；本轮未以 PDF 对此单项复核，不能据此宣称确定的物理错误。
- 脚注 77 的负离散复共轭规范未展开；C13/C14 是本稿显式补选且已逐生成元核对的规范。反射作用在树重数空间的具体相位没有全部算出。
- C15–C24 是附录 E 规范下的完整运动学分解和 Bose 筛选；C26–C28 给完整酉重耦合及投影的定义。没有提供三 principal 全重数 Racah 核的闭式特殊函数表达式。
- C29 是提交给 FP 二次型计算的目标，不是本稿结果。若某个通道失败，应保存该通道、调节方案与剩余常数，不以 N=1 或单个 S 通道的成功替代它。

## 附录 A. PP→D 的原始相位公式

以下转录原 E.30 后的两组公式并用 k 代替其正整数 n。设 k≥Δ、\((z)_r=\Gamma(z+r)/\Gamma(z)\)，定义实标签下取正根的

\[
A_{\Delta k}(s_1,s_2)
=\sqrt{\frac{(2\Delta-1)
\Gamma(\Delta+i(s_1-s_2))\Gamma(\Delta-i(s_1-s_2))}
{\Gamma(k-\Delta+1)\Gamma(k+\Delta)}}.
\tag{A1}
\]

当 n₁+n₂=k 时，

\[
\begin{aligned}
(\widetilde D^+)^{\Delta,k}_{s_1,n_1|s_2,n_2}
={}&
\frac{(-i)^{\Delta+k}\cosh(\pi s_1)
\Gamma(\frac12-n_1-is_1)\operatorname{csch}(\pi(s_1+s_2))
(\Delta+is_1-is_2)_{k-\Delta}}
{\Gamma(\frac12-n_1+is_2+\Delta)
\Gamma(1-\Delta-i(s_1+s_2))}
A_{\Delta k}(s_1,s_2)\\
&\times{}_3F_2\!\left(
\begin{matrix}\frac12-n_1-is_1,\ \Delta-k,\ \Delta-is_1+is_2\\
1-k-is_1+is_2,\ \frac12-n_1+is_2+\Delta
\end{matrix};1\right).
\end{aligned}
\tag{A2}
\]

当 n₁+n₂=−k 时，

\[
\begin{aligned}
(\widetilde D^-)^{\Delta,k}_{s_1,n_1|s_2,n_2}
={}&
\frac{i^{\Delta+k}\cosh(\pi s_2)\operatorname{csch}(\pi(s_1+s_2))
\Gamma(\frac12-k-n_1-is_2)
(1-k+is_1-is_2)_{k-\Delta}}
{\Gamma(1-\Delta-is_1-is_2)
\Gamma(\frac12-k-n_1+is_1+\Delta)}
A_{\Delta k}(s_1,s_2)\\
&\times{}_3F_2\!\left(
\begin{matrix}\frac12-k-n_1-is_2,\ \Delta-k,\ \Delta+is_1-is_2\\
1-k+is_1-is_2,\ \frac12-k-n_1+is_1+\Delta
\end{matrix};1\right).
\end{aligned}
\tag{A3}
\]

不满足相应磁守恒时系数为零。3F2 的 Δ−k 使级数终止；等质量等退化点仍须对完整乘积取连续/解析极限。C10 使用 A3 的解析参数 −s，不能只靠原 s₁,s₂>0 情形下各单因子的分支直觉评价它。

## 附录 B. PD⁺→D⁺ 的原始公式

为了连最终离散支的规范也明确，原 E.35 在本稿变量中为：δ=1/2+ia，n=m+k≥κ≥1、k≥Δ，

\[
\begin{aligned}
(\widetilde E^+)^{\kappa,n}_{a,m|\Delta,k}
={}&(-i)^m\Gamma(1-m-\delta)\Gamma(1-n-\delta+\Delta)
(\kappa+\delta-\Delta)_{n-\kappa}\\
&\times
\sqrt{\frac{\Gamma(2\Delta)\Gamma(\kappa-\delta+\Delta)
\Gamma(\kappa+\delta+\Delta-1)}
{\Gamma(2\kappa-1)\Gamma(1-\kappa-\delta+\Delta)
\Gamma(-\kappa+\delta+\Delta)}}\\
&\times
\sqrt{\frac{\Gamma(2\kappa)(2\Delta)_{k-\Delta}}
{(n-\kappa)!\,(k-\Delta)!\,\Gamma(n+\kappa)}}\\
&\times{}_3\widetilde F_2\!\left(
\begin{matrix}\kappa-n,\ 1-m-\delta,\ \kappa-\delta+\Delta\\
-m+\kappa+\Delta,\ 1-n-\delta+\Delta
\end{matrix};1\right).
\end{aligned}
\tag{B1}
\]

这里 \({}_3\widetilde F_2={}_3F_2/[\Gamma(b_1)\Gamma(b_2)]\) 是 regularized 函数；平方根沿原 CG 的连续规范使用。负离散输出由 C14 定义。C12 的离散求和采用这一单位权规范，而非另加 Plancherel 因子。
