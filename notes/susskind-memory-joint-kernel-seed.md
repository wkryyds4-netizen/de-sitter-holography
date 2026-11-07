# 联合涨落核的继续种子：尚待独立核查

2026-11-07。接 `susskind-memory-uniform-bound`。这是下一步计算的工作假设和校准表，不是已完成的有限 N 强极限定理；不纳入九页稿的已证结论。来自文献协作者的只读分析，root核了下述 q 背景公式在二/三次幂的代数系数，尚未编程独立认证或完成三色迹证明。

## 首先尝试去掉 four-u 项

目前 `m ε²=O(1)`，所以若能证明 `I_s(H)` 在紧实时间窗的 annealed L² 领先极限是逐样本能量函数，就无需精确到 O(1/p)：pinching 后 four-u 项及其与两个有界分量的 cross 全部消失。这会把全核计算降为 charged-transfer defect 与 neutral filter residual 的两分量问题。

可操作路径：以独立 u、v 副本将 `||E_u(word)−ι_s(H)||₂²` 写成三色 neutral Gaussian 混合矩；检查共享 H 的共同 J 收缩，应用已读 BKS 条件期望识别 q 背景；再用矩界/时间展开对紧 s 作统一强逼近。须证明 conditional average 的 L² 回归，不能只算其 annealed 均值。先固定 λ_N→λ 的中央紧 charge 窗，再推广 binomial 全迹。若该步骤不成立，回到完整三分量联合核。

## q-Gaussian 背景校准（已有 Wick 结构的推论）

在真正的 q-Gaussian 代数中，h、u 正交，方差 C，`v=αh+√εu`，`α²+ε=1`。W_n 为方差 C 的 monic q-Hermite；`P_r W_n=r^n W_n`，`D_q W_n=[n]_q W_{n−1}`，`N_q W_n=nW_n`。由 BKS Prop2.9/Thm2.11 的左右乘法与条件期望可推出候选校准式

\[
E_{W^*(h)}[uW_n(v)u]
=Cq^n\alpha^n W_n(h)
+\epsilon C^2[n]_q[n-1]_q\alpha^{n-2}W_{n-2}(h).
\]

后续应独立核所有低阶词及指数函数的延拓，再采用其推论

\[
\Lambda^{q}_{OU}(s)=\epsilon C P_{q\alpha}e_s(h)
+\epsilon^2 C^2P_\alpha D_q^2e_s(h),\qquad e_s(x)=e^{isx}.
\]

从 ε=0 的二阶系数抽取 I_s 的 q 背景，得到

\[
\iota_s(h)=\frac C2P_q[(N_q-x\partial_x)e_s](h)
-C^2D_q^2 e_s(h).
\]

此处 P_q 作用于括号内全部内容，x∂_x 是通常多项式微分，N_q 是 q-Hermite number operator，不是有限 N coupling chaos number。低阶校准为

\[
\iota_0=\iota'_0=0,\quad
\iota''_0=(2+q)C^2I,\quad
\iota'''_0=iC^2(1+4q+3q^2+q^3)h.
\]

取 q=exp(−4λ)、C=J_s²/λ 只限于本项目规定的 DS 极限。这些能量函数不能直接替成有限 N 的 I_s。

## 两个真实涨落分量的低阶待核接口

定义 `A_N=(1−P_h)m e^(−ihs/2)(Λ−Λ_OU)e^(−ihs/2)`，
`B_N=mε(1−P_h)e^(−ihs/2)Φ(e^(iαhs))e^(−ihs/2)`。
确切的有限 N 恒等式为

\[
B_N(0)=0,\quad B_N'(0)=iB_1,\quad B_N''(0)=-B_2,
\]
\[
B_1=m\epsilon\alpha(1-P_h)\Phi(h),\qquad
B_2=m\epsilon(1-P_h)[\alpha^2\Phi(h^2)-\alpha\{h,\Phi(h)\}].
\]

协作者提出 `p[μ(kp−d)−μ(kp)]→4λq^k d`（k=1,2）接旧 rank 尾与 charge 谱，或可得 B1→旧线性场 F，B2→`4λCq²(Y−+Y+)+2(q−1)J_HF`。这些记号须先回读 oriented/central/raw 谱的原定义和反项，再独立推导、数值校准；此处不把候选式记作已证。

如果候选正确，A0–B2 的 cross 一般不为零，因此不能从零/一次系数的正交推出完整 A(s)–B(u) 正交。下一份新证书宜先检验以上有限系数、rank slope、旧联合谱的系数与符号，再扩全部时间。

文献入口见 `literature/susskind-memory-uniform-bound-sources.md`；本记录只作可恢复的下一步定位，尚无新计算结果。
