# 从完整记忆核返回双局域涨落：来源与下一计算

2026-11-07，本地来源入口。承接 [完整联合记忆核](../notes/susskind-memory-joint-regression.md)；本页给出下一计算的精确有限维起点及待证极限，不另行宣布已完成 bilocal 重建。沿用共同 Gaussian 耦合、标准 CAR、complex p-body 模型、β=μ=0 全归一 annealed trace、ℏ=1 和 string 时间；偶 N，p²/N→固定 λ>0，J_s 固定。以下时间限于固定紧窗。

## 先完成记忆到实际 bilocal 的返回

[bilocal-channel](../notes/susskind-bilocal-channel.md) 已有固定电荷 m 的精确压缩表示。记 h=H_m，K=1_N⊗H_(m−1)，Vψ=m^(−1/2)Σ_i|i〉c_iψ，P=VV†，R=(1−P)KV，K_⊥=(1−P)K(1−P)。为避免与记忆涨落混名，压缩传播子记为

\[
F_{\rm comp}(t)=V^\dagger e^{iKt}V,\qquad
F_{\rm comp}'=i\alpha hF_{\rm comp}-M*F_{\rm comp},
\quad M(s)=R^\dagger e^{iK_\perp s}R,\quad\alpha=1-p/m.
\]

这里 \((A*B)(t)=\int_0^t A(t-s)B(s)ds\)，乘法次序保留。中心双局域算符及本轮目标为

\[
S_c(t)=m e^{-iht/2}F_{\rm comp}(t)e^{-iht/2},\qquad
\mathcal B_N(t)=(1-\mathcal P_h)S_c(t),
\]

\(\mathcal P_h\) 是逐样本完整能量块 pinching。目标是这个已减静态背景的场；旧结果 \(S/m\to1\) 没有给出它。m<p 的平凡块沿用旧约定。

取 \(M_0=\mathcal P_hM\)、\(\delta M=M-M_0\)，令 \(F_0\) 解相同 Volterra 方程、以 \(M_0\) 代替 M，且 \(F_0(0)=1\)。有限矩阵 Laplace resolvent 恒等式给出精确有序差式

\[
F_{\rm comp}-F_0=-F_0*\delta M*F_{\rm comp}
=-F_0*\delta M*F_0
+F_0*\delta M*F_0*\delta M*F_{\rm comp}.
\]

\(F_0\) 能量块对角，因而其对 \(\mathcal B_N\) 的贡献严格消失。即便能量简并，也不必假设 \(F_0\) 在每块内为标量。已完成的记忆场满足精确式

\[
\mathcal F_N(s)=(1-\mathcal P_h)m e^{-ihs/2}M(s)e^{-ihs/2},
\qquad
\delta M(s)=m^{-1}e^{ihs/2}\mathcal F_N(s)e^{ihs/2}.
\]

**下一待证命题**：在上述范数及时间窗中，将线性项两侧 \(F_0\) 换成 \(e^{iht}\)，并证明 m 放大后的第二项趋零。若此强替换成立，实际 bilocal 的极限由以下线性滤波给出（δ≥0）：

\[
\mathcal B(\delta)=
-\int_0^\delta ds\int_{-(\delta-s)/2}^{(\delta-s)/2}du\,
e^{iHu}\mathcal F(s)e^{-iHu}.
\tag{1}
\]

这里的 H 及场可先理解为已构造的联合矩/谱表示；(1) 尚不是新的微观强极限定理。推导中令两侧传播时长为 a,b，a+b=δ−s，u=(a−b)/2，故内积分没有额外 Jacobian 系数。

旧 [return-background](../notes/susskind-memory-return-background.md) 已给力高矩、Schatten Hölder 和 pΛ 的强背景近似，可作为余项估计的输入。中央电荷区的尺度是 \(\alpha-1=O(p^{-1})\)、背景记忆 \(O(p^{-1})\)、中心记忆 \(O(m^{-1})\)；真正需要证明的是这些量相乘后的范数控制及全 charge 尾。完整两点核只提供 L² 信息，不能单独控制两个随机 \(\delta M\) 的乘积。可检验的充分路线是固定次数高矩加可和时间尾，再用 L⁴/Hölder；也可寻找不要求完整 L⁴ 极限的直接余项界。

若 (1) 成立，定义
\(C_{\mathcal F}(T;s,t)=\tau[\mathcal F(s)^\dagger e^{iHT}\mathcal F(t)e^{-iHT}]\)，则

\[
\begin{split}
C_{\mathcal B}(T;\delta,\epsilon)
={}&\int_0^\delta ds\int_0^\epsilon dt
\int_{-(\delta-s)/2}^{(\delta-s)/2}du
\int_{-(\epsilon-t)/2}^{(\epsilon-t)/2}dv\,
C_{\mathcal F}(T+v-u;s,t).
\end{split}
\tag{2}
\]

双能量表示中，E 为输入、E′ 为输出，Δ=E′−E；内积分为实偶滤波

\[
\chi_L(\Delta)=\frac{2\sin(\Delta L/2)}{\Delta},\qquad
\chi_L(0)=L.
\]

应将它作用于已完成的整个 AA/AB/BA/BB 谱矩阵，保留记忆场原有的 \(e^{-is(E+E')/2}\) 相位和非零交叉项。正谱矩阵经该线性滤波仍正。必要校准是 \(\mathcal B(0)=0\)、\(\mathcal B(-\delta)=\mathcal B(\delta)^\dagger\)，以及 \(\mathcal B(\delta)=-\delta^2Y_-/2+o(\delta^2)\)，从而两点首项为 \(\delta^2\epsilon^2C_{Y_-}(T)/4\)。这些检查应在新输出中进行，不覆盖旧证书。

## 原文已经给出的下一物理接口

1. **Itakura，hep-th/9604032v2**，[*Boson Expansion Methods in (1+1)-dimensional Light-Front QCD*](https://arxiv.org/pdf/hep-th/9604032)。本轮核读 (2.15)–(2.18)、(4.18)–(4.29)、(5.1)–(5.6)：Fock 真空正规序及动量符号分区给 bilocal 代数的中心项；Holstein–Primakoff 表示以 \(\sqrt{N_c-A}B\) 和其伴随实现非线性 dressing，\(A\sim B^\dagger B\)。随后二次 Hamiltonian 给 ’t Hooft 核，三次项见 (4.29)。介子投影为
   \[
   b_r^\dagger(P)=\sqrt P\int_0^1dx\,\varphi_r(x)B^\dagger(Px,P(1-x)),
   \]
   内积为 L²(dx)。这里用 \(N_c\) 区分原文颜色数与本项目费米子数。原文质量约定是 \(M^2=2P^+P^-\)，项目用 \(M^2=P^+P^-\)，Hamiltonian 前因子需同步转换；有限颜色还须保留 boson 物理子空间。

2. **Cavicchi–Di Vecchia–Pesando，hep-th/9306091v1**，[*The Master Field of QCD₂ and the ’t Hooft Equation*](https://arxiv.org/html/hep-th/9306091)。本轮读 (26)–(35)：先取 master field 周围的二次 collective action，再用总动量和相对动量表示涨落；(31) 积去相对光前能量，相当于同一光前时间，(35) 得内部方程。原文脚注 4 明确指出不同光前时间的该 bilocal 不具有剩余规范不变性。可借用其构造顺序，不能直接把本项目双时间当成其等光前时间双空间坐标。

3. **FLZ 0905.2280v2** 的内部核及质量谱复用 [bulk-target 来源记录](susskind-bulk-target-sources.md) 和 [已存结果](../notes/susskind-bulk-target.md)，本轮不重算。项目规范下 \(\mathsf M^2=(\tau/\pi)K_0\)、\(K_0\varphi_r=2\pi^2\nu_r\varphi_r\)，所以 \(M_r^2=2\pi\tau\nu_r\)。这些是既有候选体理论的靶标。

## 怎样接上，而不跳过物理构造

先完成 (1) 的余项证明及 (2) 的实际 bilocal 协方差。协方差定义正的 Gram/GNS 度量；按其平方根作归一化并不等于 Itakura 的 canonical dressing。当前有限矩阵归一迹严格满足 τ([A,B])=0，不能在保持这个态和相应矩收敛时直接读成真空期望非零的 CCR。下一步需明确态、极化、正规序及物理子空间，再用 [bilocal-interface](../notes/susskind-bilocal-interface.md) 保存的四项 CAR 对易子检验代数。

其后需要从微观结构得到纵向变量 (P,x) 和实际 \(P^-\)，才可将二次内部算符与已有 \(K_0\) 比较，并继续相互作用。这里 x 是单介子内部夸克纵向分数；Δ 是中心时间能量差，s 的共轭量来自相邻 charge 的传播。三者尚未建立映射，不能把记忆频率命名为质量或任选 \(P^-\) 拟合体谱。

本页只新增下一步来源与计算入口；原文通过 arXiv 在线核读，未声称读取本机不存在的论文缓存。有限 Volterra 差式是标准 resolvent 应用，(1) 的模型极限及微观光前字典仍是待做工作。
