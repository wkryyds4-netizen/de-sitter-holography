# 双局域低通道的共同代数与实时间四点

2026-11-07。承接完整 bilocal 两点重建，研究实际 charged 曲率 A0 和 neutral 一次系数 B1 的共同非交换结构。模型为 standard CAR、overlapping complex Hermitian Gaussian p-body、同一原始 J、beta=mu=0 归一 annealed full-Fock trace；偶 N，p²/N→固定 lambda>0，Js 固定，hbar=1/string 时间。所有极限先固定次数、固定中心时间。正式稿为 [LaTeX](../manuscript/susskind-bilocal-four-point/main.tex) 与 [九页 PDF](../output/pdf/de-sitter-susskind-bilocal-four-point.pdf)。

## 先辨认旧结果

[原文核读](../literature/susskind-bilocal-four-point-sources.md)覆盖 BNR local CAR/Poisson、BINT crossed/uncrossed kernel、Junge–Zeng mixed-q 正迹、Effros–Popa 多块 Wick 和 Itakura canonical dressing。旧 singlet-multimark 已给 J2 对应场与 H 的全部固定混合矩及四点，不能把 tau J2^4=(2+q^4)w² 重新说成新发现。

同一样本精确恒等式是 A=sum[H,c_i]†[H,c_i]=(2p−D)H²，O=(2p−D)Z2，J2−Y_old=kappa2(C−C_Q)。有界 sector 方差 C_Q→C 使其在每个固定 Schatten 范数下消失。实际 A0=(1−P_H)mR†R=(1−P_H)J2 完全精确；旧无对角质量加高矩将 pinching 差从 L² 升为所有固定 Lr。完整证明见 [实际接口](susskind-bilocal-four-point-transfer.md)。

## 本轮接入 W1 的结果

令 C=Js²/lambda、q=exp(−4lambda)，L 是 half-Clifford-size，W1=(p−L−lambda_N)H，Z2=Gaussian-Wick2(H²)，

\[
J_2=(2p-D)Z_2-\frac{2\lambda(1-q)}{1+q}(H^2-C),
\quad w=C^2\left(4\lambda+\frac{16\lambda^2q}{1+q}\right).
\]

归一三字母 h=H/sqrt(C)、u=W1/sqrt(C lambda)、x=J2/sqrt(w) 的任意固定 word 由同种字母配对与交叉矩阵给出：

\[
 Q=\begin{pmatrix}q&q&q^2\\q&q&q^2\\q^2&q^2&q^4\end{pmatrix},
 \qquad (\ell_h,\ell_u,\ell_x)=(1,1,2).
\]

新增的是 actual common-J 的 W1 接入旧 J2/H 定理，不是另创 mixed-q 数学框架。源推导先完整保留每个 nilpotent mask、带源 singleton 分母及有符号局部迹；精确 singleton 相消后，B^R_overlap 与 hypergeometric 指数矩控制所有固定高矩。先支配、后删 triple，再取 pair-Poisson。

在互不交叠 pure-Wick D/L mark 类，pair collision 只有四个 incidence；一个 source 至少消耗两个 incidence 才能在归一后存活。因此三阶以上的局部指数系数为零。这有额外局部证明，并非“没有 triple 交叠”自动意味着高阶 source 消失；完整四源 moment 仍有非零项。

W1 只有 own creator/annihilator collision 留源。单源 kappa=1 被 −lambda_N H 消去；WW connected kappa=1 当且仅当同一 Gaussian 弦两端，WD 为零。所以 W 必须直接配 W，每对 C lambda；它们作为 spectator 不破坏旧 J2 双弦交换相消。见 [多源完整证明](susskind-bilocal-multisource-proof.md)。

## 对应实际量，不只参考模型

actual neutral slope B1=m epsilon alpha(1−P_H)(Phi h−kappa_bar h) 已有 L² 桥到 8lambda Cq W1。本轮利用其未 pinching 部分是 U(N)-equivariant Hermitian 线性 Gaussian 场，Schur 给 E X_m²=v_m I，Gaussian 分部积分给 E tr X_m^(2k)≤(2k−1)!! v_m^k。旧 broad-central 方差界与 binomial 外尾闭合所有固定高矩，插值遂给每个固定 r：

\[
\|B_{1,N}-8\lambda CqW_{1,N}\|_r\to0,
\qquad\|A_{0,N}-J_{2,N}\|_r\to0.
\]

不能把这个一次 Gaussian 论证套到任意高阶 rank 后代。对固定生成元，H 高矩的 sqrt(s) 界使真实轴 exponential Taylor 的 Sr 余项紧时间一致趋零，再用 Hölder 传递任意固定实时 word。

新非交换诊断为

\[
\tau(x^2u^2)=1,\quad\tau(xuxu)=q^2,\quad
\operatorname{Gram}(xu,ux)=\begin{pmatrix}1&q^2\\q^2&1\end{pmatrix},
\quad\|[x,u]\|_2^2=2(1-q^2)>0.
\]

commutator 零迹但非零范数，故不是 scalar。二点正交并不意味着可交换；该结论也不将当前迹态改成真空 CCR。后小 lambda 的归一交换子平方范数为 16lambda+O(lambda²)。

## 实时间 mixed OTOC

F_xu(T)=tau[x(T)u x(T)u]，D=2(1−F)。uncrossed 两相邻同色对不与 H 弦交叉，剩余传播子相消，故 tau[x(T)²u²]=1。使用已有 BINT 正级数，长度 (2,1)，HX=q²、HU=q，互 probe 因子 q² 必须恢复；只有一个 species-compatible crossed channel。

\[
F_{xu}(T)=q^2\{1-2C(1-q^2)(1-q)T^2+O(T^4)\}.
\]

F 的非负性来自正级数；D≤2 不是任意 Hermitian pair 的迹恒等式自动给出的上界。新检查独立枚举 0/2/4/6 阶，共 130 arc、1376 colored pairings/每 species 次序。Js=1，lambda=.25/.1/.05，各22时点，L=61/130/245、L+16、swap 共198求值；t=1 的 D 约 1.99594313/1.97575004/1.92315023。最大细化差1.89e−15、swap差2.23e−16、small-time Taylor差1.45e−15；解析Jacobi截断界<5.44e−12，浮点检查非区间认证。绝对/相对误差、tail及相消条件均存原JSON，不能用小绝对差省去相对误差。

## 强返回与下一步

旧 contractive dilation/Volterra 证明直接推广：good charge 中 deltaM2=O(p^-2)、deltaM_(4r)=O(p^-1)、F/F0−U 的2r范数O(p^-1)。插值权1/(2r−1)，在 exact return 中乘m≤N，得到

\[
\sup_{|t|\le S}\|\mathcal B_N(t)-\widetilde{\mathcal B}_N(t)\|_r
=O_S(p^{-1/(2r-1)})+O_S(pe^{-N/(8r)}),\quad r\ge2.
\]

因此四阶返回误差已是O(p^-1/7)。下一首任务是完整滤波时间级数的 L4 可和尾：lim_K limsup_N sup_|t|≤S ||Btilde_N−Btilde_N^(K)||4=0，并识别固定截断所含更高后代的共同四点。当前任意固定次数矩界不提供degree统一可和常数；旧L²尾不能直接替换成L4尾。低通道的实时四点完成，完整finite-relative-time bilocal四点尚未完成。

随后再研究适用态/极化、canonical dressing与微观空间及内部质量生成元。中心能差、memory频率与quark momentum fraction仍是不同对象，不能自行选择P^-来拟合目标bulk谱；本轮OTOC没有被解释为finite-G散射振幅。

## 保存

四份新科学JSON均首次成功后exclusive-create且未覆盖：40032有限CAR/source系数、4768指数比较；53cyclic representatives覆盖233word/1763diagram；33精确矩阵恒等式/3030entry；66时间点/198kernel。源码、原文定位、参数、依赖版本及输入hash保留在 results/susskind-bilocal-four-point。双遍编译脚本 code/build_bilocal_four_point.ps1；95dpi逐页核版。只读recorder在最终核版后绑定科学输出、前驱、PDF与render，不执行旧main。本轮不commit/push/ZIP/外部回复；issue #1此前初核、9890旧端点失败和21143单元未独立复积分状态不变。
