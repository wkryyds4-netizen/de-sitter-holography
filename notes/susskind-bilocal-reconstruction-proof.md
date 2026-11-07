# 从完整 memory 核到实际 bilocal：强返回定理

2026-11-07。接续 [来源与计算入口](../literature/susskind-memory-to-bilocal-sources.md) 及 [joint-regression 阶段](susskind-memory-joint-regression.md)。本稿只新增 memory 到实际 centered bilocal 的重建；旧文件和科学输出不改写。

采用标准 CAR、共同原始 Gaussian 耦合的 overlapping complex \(p\)-body 模型，\(\beta=\mu=0\) 归一 annealed full-Fock 迹，\(\hbar=1\)、string 实时间。偶 \(N\)、\(p\to\infty\)、\(\lambda_N=p^2/N\to\lambda>0\)，\(J_s>0\) 固定。以下常数允许依赖固定时间窗 \(S\)、\(J_s\) 和 \(\lambda_N\) 所在正紧区间，不依赖 \(N,p,m\)。不取增长时间或同时 \(\lambda\to0\) 的极限。

全迹范数为
\[
 \|Y\|_r=\left(\sum_m b_{N,m}\mathbb E\operatorname{tr}_m|Y_m|^r\right)^{1/r},
 \qquad b_{N,m}=2^{-N}\binom Nm.
\]
这里 \(r\) 是 **Schatten 阶数**，不是 Gaussian-number Sobolev 阶数。对 charge 集合 \(I\)，\(\|Y\|_{r,I}\) 指把求和限于 \(m\in I\)；保留原来的二项权重，不重新归一化。

## 1. 精确对象及唯一需要的新强极限

在 \(m>0\) 上令
\[
 h=H_m,\quad K=I_N\otimes H_{m-1},\quad
 V\psi=m^{-1/2}\sum_i|i\rangle c_i\psi,\quad V^\dagger V=1,
\]
\[
 Q_\perp=1-VV^\dagger,\quad
 R=Q_\perp KV,\quad D_\perp=Q_\perp KQ_\perp|_{Q_\perp},\quad
 a=V^\dagger KV=(1-p/m)h.
\]
为避免与已用的 body-loss \(D\) 混淆，这里把 bath 压缩写成 \(D_\perp\)。定义
\[
 F_m(t)=V^\dagger e^{iKt}V,\qquad
 M_m(s)=R^\dagger e^{iD_\perp s}R.
\]
\(F_m\) 是收缩的矩阵传播函数。以有序 Volterra convolution
\[
 (A*B)(t)=\int_0^t A(t-s)B(s)\,ds
\]
记，精确消去 bath 得
\[
 F'=iaF-M*F,\qquad F(0)=1.\tag{1}
\]
实际 centered bilocal 与已完成的 memory 场分别是
\[
 \mathcal B_{m,N}(t)=(1-\mathcal P_h)
     [m e^{-iht/2}F_m(t)e^{-iht/2}],
\]
\[
 \mathcal F_{m,N}(s)=(1-\mathcal P_h)
     [m e^{-ihs/2}M_m(s)e^{-ihs/2}].\tag{2}
\]
\(\mathcal P_hY=\sum_eP_eYP_e\) 保留完整简并能量块。\(m<p\) 时 \(H_m=H_{m-1}=0\)、\(R=M=0\)，两个 centered 场均为零；\(m=0\) 同样置零。

主结论是，任意固定 \(S\) 上
\[
 \boxed{\sup_{|t|\le S}\|\mathcal B_N(t)-\widetilde{\mathcal B}_N(t)\|_2
       =O_S(p^{-1/3}),}\tag{3}
\]
其中对 \(t\ge0\)
\[
 \boxed{\widetilde{\mathcal B}_N(t)=
 -\int_0^t ds\int_{-(t-s)/2}^{(t-s)/2}du\,
 e^{iH_Nu}\mathcal F_N(s)e^{-iH_Nu}.}\tag{4}
\]
负时间由伴随给出，亦可把 (4) 理解为两重有向积分。此处输入的是实际有限 \(N\) memory 场，不是先替成假定的极限算符。

## 2. 静态 memory 的显式收缩 dilation

令
\[
 M_0(s)=\mathcal P_hM(s),\qquad \delta M(s)=M(s)-M_0(s),
\]
并令 \(F_0\) 解 (1) 中将 \(M\) 换成 \(M_0\) 的方程。关键事实是 \(F_0\) 也逐样本为收缩；不能只凭一般 Volterra 方程断言这一点。

写 \(h=\sum_e eP_e\)。取 bath 直和空间及耦合
\[
 \mathscr K_0=\bigoplus_e\mathscr K_\perp,\qquad
 D_0=\bigoplus_eD_\perp,\qquad
 R_0\psi=(RP_e\psi)_e.
\]
则
\[
 R_0^\dagger e^{iD_0s}R_0
   =\sum_eP_eR^\dagger e^{iD_\perp s}RP_e=M_0(s).
\]
有限维自伴块 Hamiltonian
\[
 K_0=\begin{pmatrix}a&R_0^\dagger\\ R_0&D_0\end{pmatrix}
\]
的系统压缩 \(V_0^\dagger e^{iK_0t}V_0\) 满足定义 \(F_0\) 的同一 Volterra 初值问题，因而就是 \(F_0(t)\)。所以
\[
 \|F(t)\|_\infty\le1,\quad \|F_0(t)\|_\infty\le1,
 \quad F_0(-t)=F_0(t)^\dagger.\tag{5}
\]
\(K_0\) 按“系统能量块 \(P_e\)+其 bath 副本”分解，故
\[
 [F_0(t),h]=0.\tag{6}
\]
简并块内部不要求 \(F_0\) 为标量。\(M_0\) 是时间变量上的 positive-definite 矩阵核；这不是说每个非零时间的 \(M_0(s)\) 都是正算符。

由有限矩阵 Laplace resolvent 恒等式或 Volterra 唯一性，保持乘法次序有
\[
 F-F_0=-F_0*\delta M*F.\tag{7}
\]
这是一条精确恒等式，不依赖无限返回级数的收敛。因为 \(m\) 是 charge 块内的标量，定义
\[
 Q_N(s)=m\delta M(s)
       =e^{ihs/2}\mathcal F_N(s)e^{ihs/2},\tag{8}
\]
即得
\[
 m(F-F_0)=-F_0*Q_N*F.\tag{9}
\]
既有实际 memory 二阶界给
\[
 B_S:=\sup_N\sup_{|s|\le S}\|Q_N(s)\|_2
    =\sup_N\sup_{|s|\le S}\|\mathcal F_N(s)\|_2<\infty.\tag{10}
\]
因活跃区 \(m\ge p\)，(7)--(10) 已经推出
\[
 \sup_{|t|\le S}\|F-F_0\|_2\le\frac{S^2B_S}{2p},
 \qquad
 \sup_{|t|\le S}\|m(F-F_0)\|_2\le\frac{S^2B_S}{2}.\tag{11}
\]

## 3. 旧力高矩给出 memory 的高阶粗界

这一步只使用已封存的 [return-expansion](susskind-memory-return-expansion.md) 的 Gaussian square-function 估计。令
\[
 F_i=[H,c_i],\qquad A_m=\sum_iF_i^\dagger F_i|_m,
 \qquad G_m=R^\dagger R.
\]
精确地 \(KV-Vh\) 的平方为 \(A_m/m\)，而 \(R=Q_\perp(KV-Vh)\)，故
\[
 0\le G_m\le A_m/m.
\]
前驱已证所有固定整数 \(k\) 的 sector 一致界
\[
 \mathbb E\operatorname{tr}_m A_m^k\le C_kp^k,
 \qquad
 \mathbb E\operatorname{tr}_m G_m^k\le C_k(p/m)^k.\tag{12}
\]
后一式使用正算符的本征值单调性与 trace 幂，不声称矩阵幂算符单调。本稿只需 \(k=2,4,8\)。Hamiltonian 的固定 Schatten 矩也在所有 charge 上一致有界，来自前驱的 \(C_m\le M\) 和 Gaussian 矩估计。

矩形 Schatten Hölder 全部以源维数 \(d_m\) 归一，给
\[
 \|R^\dagger e^{iD_\perp s}R\|_{k,m}
 \le\|R\|_{2k,m}^2
 = (\operatorname{tr}_mG_m^k)^{1/k}.
\]
因此 annealed 后
\[
 \|M(s)\|_{k,m}\le C_k^{1/k}p/m.\tag{13}
\]
逐样本 pinching 是 Schatten-\(k\) 收缩，故 \(M_0\) 有同样界，\(\delta M=(1-\mathcal P_h)M\) 的右端至多乘 2。这里不需要对一般 matrix-valued Gaussian 多项式作 \(L^2\to S_4\) hypercontractivity。

## 4. 广中央插值与完整 charge 尾

取 good charge 集
\[
 I_N=\{m\ge N/4\}.
\]
高 charge 无需另外去掉，因为 (12)--(13) 对它们仍成立。对充分大 \(N\)，\(p<N/4\)。由 (13)，在 \(k=4,8\) 上
\[
 \sup_s(\|M(s)\|_{k,I_N}+\|M_0(s)\|_{k,I_N}
                  +\|\delta M(s)\|_{k,I_N})
 \le K_k\frac pN.\tag{14}
\]
另一方面，精确的 (8) 与强 memory 二阶界 (10) 给出更小的中心涨落二范数：
\[
 \sup_{|s|\le S}\|\delta M(s)\|_{2,I_N}
 \le\frac{4B_S}{N}.\tag{15}
\]
在有限迹直和与 disorder 概率的乘积代数上应用非交换 \(L^p\) 插值。由
\(1/4=(1/3)/2+(2/3)/8\)，得到
\[
 \boxed{\sup_{|s|\le S}\|\delta M(s)\|_{4,I_N}
 \le\left(\frac{4B_S}{N}\right)^{1/3}
      \left(\frac{K_8p}{N}\right)^{2/3}
 =O_S(p^{-4/3}).}\tag{16}
\]
即中心化在二阶矩上提供的额外 \(p^{-1}\) 小量，通过旧的粗八阶矩产生一个严格足够的四阶改进。没有推断 \(Q_N=m\delta M\) 本身的全时间四范数一致有界。

对 bad active charge 集 \(J_N=\{p\le m<N/4\}\)，二项下尾满足
\(\sum_{J_N}b_{N,m}\le e^{-N/8}\)。用 (12) 及 \(m\) 与 pinching 的可交换性，
\[
 \sup_s\|m\delta M(s)\|_{2,J_N}^2
 \le4\sum_{J_N}b_{N,m}m^2
             \mathbb E\operatorname{tr}_mG_m^2
 \le Kp^2e^{-N/8}.\tag{17}
\]
\(m<p\) 的场严格为零，所以没有遗漏低 charge 端点或 \(m=0\)。

## 5. 两个收缩传播子均接近自由能量相位

记 \(U_a(t)=e^{iat}\)、\(U_h(t)=e^{iht}\)。由 (1) 和其静态版本，
\[
 F=U_a-U_a*M*F,\qquad
 F_0=U_a-U_a*M_0*F_0.\tag{18}
\]
对 good 区使用 (5)、(14) 和时间单纯形体积，
\[
 \sup_{|t|\le S}\big(\|F-U_a\|_{4,I_N}
                       +\|F_0-U_a\|_{4,I_N}\big)
 \le K_S p/N.
\]
又因为 \(a-h=-(p/m)h\)，实轴 Duhamel 与固定 \(H\) 四矩给
\[
 \|U_a(t)-U_h(t)\|_{4,I_N}
 \le|t|\|(p/m)h\|_{4,I_N}\le K_Sp/N.
\]
所以
\[
 \boxed{\sup_{|t|\le S}
 [\|F(t)-U_h(t)\|_{4,I_N}
 +\|F_0(t)-U_h(t)\|_{4,I_N}]=O_S(p^{-1}).}\tag{19}
\]
在所有 charge 和样本上，两个差的算符范数都至多 2。上述四范数估计直接来自收缩传播及低强度 memory，并不是把 \(L^2\) 误差不经论证提升为 \(L^4\)。

## 6. 强首返替换及速率

由 (7) 精确展开差
\[
 m(F-F_0)+mU_h*\delta M*U_h
 =-m(F_0-U_h)*\delta M*F
  -mU_h*\delta M*(F-U_h).\tag{20}
\]
在 good 区，\(m\le N\)，Schatten Hölder 及 (16)、(19) 使每个 integrand 至多
\[
 N\,O_S(p^{-1})\,O_S(p^{-4/3})=O_S(p^{-1/3}).
\]
精确的 \(p,N\) 表达是 \(O_S(p^{5/3}/N)\)；用 \(N\asymp p^2\) 才得到最后的次数。在 bad active 区使用差传播子的算符范数 \(\le2\)、其余传播子范数 \(\le1\)，直接由 (17) 控制为
\(K_Sp e^{-N/16}\)。两重时间积分的体积至多 \(S^2/2\)，故
\[
 \sup_{|t|\le S}\|m(F-F_0)+U_h*Q_N*U_h\|_2
 \le K_S[p^{-1/3}+p e^{-N/16}].\tag{21}
\]

这也解释初始设想的 nonlinear 第二返回为何确实消失。再用一次 (7) 得
\[
 F-F_0=-F_0*\delta M*F_0
       +F_0*\delta M*F_0*\delta M*F.
\]
good 区该余项乘 \(m\) 后由两个 (16) 控制为
\(N\,O(p^{-8/3})=O(p^{-2/3})\)，四重时间单纯形只给固定体积。bad 区 (13) 的四矩与 \(p/m\le1\) 给 \(O(N e^{-N/16})\)。不过证明 (21) 本身不需要先走这个二次展开。

由 (6)，\(F_0\) 的 centered bilocal 贡献精确为零。对 (21) 两侧包裹 \(e^{-iht/2}\) 并施加 \(1-\mathcal P_h\)（二范数正交投影），即可得到 (3)。具体地，在 convolution 中记左右时间为 \(a,b\ge0\)、中间 memory 时间为 \(s\)，则 \(a+b+s=t\)。由 (8)
\[
 e^{-iht/2}U_h(a)Q_N(s)U_h(b)e^{-iht/2}
 =e^{ih(a-b)/2}\mathcal F_N(s)e^{-ih(a-b)/2}.
\]
令 \(u=(a-b)/2=a-(t-s)/2\)，有 \(du=da\)，内积分无额外 Jacobian。又因
\[
 \mathcal P_h[e^{ihu}\mathcal F_N(s)e^{-ihu}]=0,
\]
外层 pinching 可在这一步 **精确删除**。于是得到 (4) 的相位、积分域和负号。

(11) 还立即给实际 bilocal 的统一界
\[
 \sup_{|t|\le S}\|\mathcal B_N(t)\|_2\le S^2B_S/2.
\]
负时间的所有界由 \(F(-t)=F(t)^\dagger\)、\(M(-s)=M(s)^\dagger\) 及 pinching 保伴随得到。

## 7. 从强返回到实际 bilocal 两点函数

记既有 memory covariance 为
\[
 C_{\mathcal F,N}(T;s,t)=\mathbb E\operatorname{tr}
 [\mathcal F_N(s)^\dagger e^{iH_NT}\mathcal F_N(t)e^{-iH_NT}].
\]
对 \(\delta,\epsilon\ge0\)，(4) 的 covariance 精确是
\[
 \int_0^\delta ds\int_0^\epsilon dt
 \int_{-(\delta-s)/2}^{(\delta-s)/2}du
 \int_{-(\epsilon-t)/2}^{(\epsilon-t)/2}dv\,
 C_{\mathcal F,N}(T+v-u;s,t).\tag{22}
\]
由 (3) 和统一二阶界，(22) 与实际 \(\mathcal B_N\) covariance 的差为 \(O_S(p^{-1/3})\)，对外部所有实 \(T\) 一致。然后使用已证 memory covariance 的固定时间极限及有界支配，得到实际 bilocal 的 limiting covariance。此处“强近似误差对所有 \(T\) 一致”不等于“有限 \(N\) 谱本身在全实轴 \(T\) 一致收敛”。

对固定 \(T\)，memory 时参落在紧集，四重积分均有限。\(\delta\) 的变化仅改变 (4) 的三角积分域，统一场范数给局部 Lipschitz 控制；因此亦可在 \(\delta,\epsilon\) 的固定紧窗上取一致极限。负时参继续用伴随关系。

## 8. 共同正矩阵谱的滤波与 charged 简化

沿前驱谱约定，\(E\) 是 input/column、\(E'\) 是 output/row，设
\[
 e=(E+E')/2,\qquad \Delta=E'-E,\qquad \omega=\Delta/2,
\]
\[
 \chi_L(\Delta)=\frac{2\sin(\Delta L/2)}{\Delta}
               =\frac{\sin(\omega L)}\omega,\qquad \chi_L(0)=L.
\]
memory 的共同场族是 \(\mathcal X_a\in\{J_d\ (d\ge2),N_j\ (j\ge1)\}\)，谱为前驱的正矩阵 \(\nu_{ab}\)。其 amplitude 为
\[
 f_{J_d}(s)=D_d(s;E,E')=-c_d''+i(E+E')c_d'+EE'c_d,
 \qquad f_{N_j}(s)=2C c_j(s),
\]
且已有共同 centering phase \(e^{-ies}\)。因此 bilocal amplitude 为
\[
 b_a(\delta;E,E')=-\int_0^\delta
       \chi_{\delta-s}(\Delta)e^{-ies}f_a(s;E,E')\,ds.\tag{23}
\]
实际 limiting bilocal covariance 是
\[
 \boxed{C_{\mathcal B}(T;\delta,\epsilon)
 =\sum_{a,b}\int e^{iT\Delta}
 \overline{b_a(\delta;E,E')}b_b(\epsilon;E,E')\,d\nu_{ab}(E,E').}\tag{24}
\]
所有 AA、AB、BA、BB 都保留。\(|\chi_L|\le|L|\)，所以前驱的绝对可和场尾经过紧时间滤波仍绝对可和；固定能量支撑、total variation 收敛和无对角质量都保留。对角点的 \(\chi_L(0)\) 是连续延拓，不是移除零频质量的另一个假设。

charged 部分可完全化简。令 \(g_d(s)=e^{-ies}c_d(s)\)，则
\[
 e^{-ies}D_d(s)=-g_d''(s)-\omega^2g_d(s).
\]
\(\chi_t\) 是 \(\partial_t^2+\omega^2\) 的零初值 retarded Green 函数，因此
\[
 b_{J_d}(\delta)=g_d(\delta)-\cos(\omega\delta)g_d(0)
                         -\chi_\delta(\Delta)g_d'(0).
\]
对 \(d\ge2\)，\(c_d(0)=c_d'(0)=0\)，所以
\[
 \boxed{b_{J_d}(\delta;E,E')=e^{-ie\delta}c_d(\delta).}\tag{25}
\]
该恒等式也适用于 \(\omega=0\) 的连续极限。

neutral amplitude 保留一维积分，或者写成无重节点奇性的差商：
\[
 \boxed{b_{N_j}(\delta;E,E')=
 \frac{2C e^{-ie\delta}}{C^jP_j}
 \int H_j^{(q,C)}(x)\,
       [\exp(i\delta\,\cdot)](E,E',x)\,d\rho_q(x).}\tag{26}
\]
方括号是二阶 divided difference，包括重节点的 Hermite 连续延拓。其符号可由三个节点全零时
\([e^{i\delta\cdot}](0,0,0)=-\delta^2/2\)
核对；不要在近重合节点直接用相消严重的三个分式。

## 9. 低时间校准和本轮证明范围

由 (4)，\(\mathcal B(0)=0\)、\(\mathcal B(-\delta)=\mathcal B(\delta)^\dagger\)。前驱给 \(\mathcal F(0)=Y_-=J_2\)，且时间场在 limiting Hilbert 空间连续。故
\[
 \mathcal B(\delta)=-\frac{\delta^2}{2}Y_-+o_{L^2}(\delta^2),
\]
\[
 C_{\mathcal B}(T;\delta,\epsilon)
 =\frac{\delta^2\epsilon^2}{4}C_{Y_-}(T)
   +o(\delta^2\epsilon^2)
\]
在两正小时间都趋零时成立；这里的小时间极限是在已构造的固定 \(\lambda\) 理论内。正矩阵谱还保证等时范数非负及
\(C_{\mathcal B}(T;\delta,\epsilon)^*
 =C_{\mathcal B}(-T;\epsilon,\delta)\)。

本稿闭合的是实际 finite-\(N\) centered bilocal 到已知 memory 谱的强返回及完整两点重建。证明只用旧实际 memory 的统一二阶界、旧力的 \(2/4/8\) 阶矩、显式收缩 dilation 和 Schatten 插值。它不要求全时间 \(m\delta M\) 的四阶有界，也不需要新证明四标记图的任意高阶极限。此前提出的固定次数 \(S_4\) proxy 路线可作为备选，但不构成本定理的前置条件。

新的代数/数值滤波检查由本阶段单独输出记录；本稿没有运行或覆盖任何旧 main。将此 bilocal 正 Gram 度量进一步解释为 canonical meson dressing，仍须指定态、极化、正规序、物理子空间以及实际光前 Hamiltonian；(24) 本身不提供这些尚未建立的物理字典。
