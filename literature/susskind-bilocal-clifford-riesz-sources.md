# Clifford Riesz 界与 charged 乘积项的闭合接口

2026-11-07。本文件补充已经被计算哈希绑定的
[fourth-tail sources](susskind-bilocal-fourth-tail-sources.md)，不修改该文件。
目的仅是核查一个实际需要的数学工具：从 X 与 LX/p 的 Schatten 界控制
body/Clifford 乘积规则的 cross 项。这里 L 是物理 Majorana half-length，
不是 coupling Gaussian number，也不是表示论 rank。

## 原文已找到：不是待猜测的不等式

L. Ben Efraim–F. Lust-Piquard，*Poincaré type inequalities on the discrete cube and in the CAR algebra*，
Probability Theory and Related Fields **141** (2008), 569–602，
[期刊信息](https://link.springer.com/article/10.1007/s00440-007-0094-x)，
[已核全文 math/0702233v1](https://arxiv.org/html/math/0702233v1)。
具体用 §6.1 的 D'_i、N' 和 symmetrized gradient；§6.4 Theorem 6.8(1),(3) 的 beta=1/2 情形；
Lemma 5.6/(16) 的收缩半群插值。Theorem 6.1 是反方向的界，不能把它误当所需结论。

对固定 2≤s<infinity，Theorem 6.8 给

\[
 \left\|\left[\sum_i (D'_iX)^*D'_iX+(D'_iX)(D'_iX)^*\right]^{1/2}\right\|_s
 \le 2H_s\|(N')^{1/2}X\|_s. \tag{1}
\]

H_s 是文中 Hilbert-transform/UMD 常数，只依赖 s，可对矩阵维数一致选择。
对应的先驱是 F. Lust-Piquard，*Riesz transforms associated with the number operator on the Walsh
system and the fermions*，JFA **155** (1998), 263–285，
[DOI 10.1006/jfan.1997.3217](https://doi.org/10.1006/jfan.1997.3217)。
1998 期刊检索内容核到摘要；本文件具体公式依赖已经读到全文的 2007 版本，而不声称读过被拒绝访问的 1998 PDF。

维数一致性的另一个主来源是 G. Pisier–Q. Xu，*Non-commutative martingale inequalities*，
[math/9704209v1](https://arxiv.org/html/math/9704209v1)，Introduction (0.1)、(0.3)–(0.4)。
其非交换 martingale 常数只依赖范数指数，并说明 Schatten S_s 的 UMD 性。
有限维 Schatten 空间嵌入同一个 S_s；改变归一化 trace 只同时缩放两边，不改变这些常数。
这里不需要知道 H_s 的最优数值。

## 与项目规范完全一致的映射

项目有 2N 个 Majorana gamma_i，{gamma_i,gamma_j}=2delta_ij，parity automorphism 为 alpha。
定义

\[
 \nabla_iX=\frac12(\gamma_iX-\alpha(X)\gamma_i).
\]

对有序 Clifford 单项 gamma_A，若 i∈A，则 nabla_i gamma_A=gamma_i gamma_A，否则为零。
这与原文 D'_i 的定义一致，因此

\[
 N'=\sum_i\nabla_i^*\nabla_i=2L,\qquad
 L\gamma_A=\frac{|A|}{2}\gamma_A. \tag{2}
\]

原文用 n 个 Clifford 生成元的较大矩阵表示。取 n=2N 时，这只是项目 full-Fock Clifford
代数的 trace-preserving amplification；归一化 Schatten 范数和 (2) 均保持不变。

以下定义列、行与合并平方函数：

\[
 G_C(X)=\left[\sum_i(\nabla_iX)^*(\nabla_iX)\right]^{1/2},\quad
 G_R(X)=\left[\sum_i(\nabla_iX)(\nabla_iX)^*\right]^{1/2},
\]
\[
 G(X)=[G_C(X)^2+G_R(X)^2]^{1/2}.
\]

正矩阵序给 max(||G_C||_s,||G_R||_s)≤||G||_s。
由 (1)–(2) 得

\[
 \boxed{\|G(X)\|_s\le2\sqrt2 H_s\|L^{1/2}X\|_s.} \tag{3}
\]

这控制实际有限 Clifford 矩阵，适用于任意 X，不需要先取 double scaling。
它与先前不成立的“coupling OU tensor matrix identity 的 S2→S4 超可积”是不同命题。

## 无需额外高范数损失的能量插值

Clifford OU exp(−tL) 为归一化 trace 上的收缩半群。上述 Lemma 5.6，beta=1/2，给

\[
 \|L^{1/2}X\|_s\le4\|X\|_s^{1/2}\|LX\|_s^{1/2}.
\]

合并 (3)：

\[
 \boxed{\|G(X)\|_s^2\le128H_s^2\|X\|_s\|LX\|_s.} \tag{4}
\]

这里所有范数阶数都是同一个 s；不需要为这一步把 s 加倍。
对 common-J 随机矩阵，先逐样本使用 (4)，再对 Gaussian expectation 使用 Hölder，仍有

\[
 \|G(X)\|_{s,\mathrm{ann}}^2
 \le128H_s^2\|X\|_{s,\mathrm{ann}}\|LX\|_{s,\mathrm{ann}}. \tag{5}
\]

例如 E(a b)^(s/2)≤(Ea^s)^(1/2)(Eb^s)^(1/2)，其中 a=||X(J)||_s、b=||LX(J)||_s。
有限 N 的多项式均在所需算符域内；对完整时间场使用时还需原有 Sobolev/截断极限的域控制。

## 直接应用到乘积 cross 项

对 even X,Y，half-length product rule 为

\[
 L(XY)=(LX)Y+X(LY)+B_L(X,Y),\qquad
 B_L(X,Y)=\sum_i(\nabla_iX)(\nabla_iY). \tag{6}
\]

把第一组导数视为一行，第二组视为一列，矩形 Schatten Hölder 给

\[
 \|B_L(X,Y)\|_{r,\mathrm{ann}}
 \le\|G_R(X)\|_{2r,\mathrm{ann}}\|G_C(Y)\|_{2r,\mathrm{ann}}.
\]

因而对 p>0、固定 r≥2，

\[
 \boxed{\left\|\frac{B_L(X,Y)}p\right\|_{r,\mathrm{ann}}
 \le128H_{2r}^2
 \left[
 \|X\|_{2r,\mathrm{ann}}\left\|\frac{LX}p\right\|_{2r,\mathrm{ann}}
 \|Y\|_{2r,\mathrm{ann}}\left\|\frac{LY}p\right\|_{2r,\mathrm{ann}}
 \right]^{1/2}.} \tag{7}
\]

body-loss D 的 cross 项也直接受控。令
gamma_(2i−1)=c_i+c_i^dagger、gamma_(2i)=−i(c_i−c_i^dagger)，则对 even X

\[
 [c_i^\dagger,X]=\nabla_{2i-1}X-i\nabla_{2i}X,
 \quad [c_i,X]=\nabla_{2i-1}X+i\nabla_{2i}X.
\]

每个 2→1 坐标块的算符范数为 sqrt2，因此

\[
 B_D(X,Y)=\sum_i[c_i^\dagger,X][Y,c_i]
\]

满足与 (7) 相同的右端、常数改为 256H_(2r)²。
这正对应 D(XY)=(DX)Y+X(DY)−B_D(X,Y)；reverse-body 情形同理。
grading、导数顺序和原矩阵乘法顺序均保留。

由 (7) 得到两个具体充分条件：

- 若 X、LX/p、Y、LY/p 在 L^(2r)_ann 一致有界，则 cross/p 在 L^r_ann 一致有界。
- 若另外 X→0，而其余三项保持有界，则 cross/p→0，界中至少得到 ||X||^(1/2) 的消失。

可取 X=A_v、Y=H^w 接 CP-path 表示，再沿 eta 的受迫递推估计。
这里仅提供解析工具；A_v 的 L/p 界与 degree 包络应由该路径证明明确建立。
常数 H_(2r) 只依赖固定范数阶，故这一步不新增随 v,w 指数平方增长。

## 使用范围中的一个实质区别

(3)–(7) 是 **full Clifford normalized trace** 的定理，可以直接转为 annealed full-Fock。
不能不加证明地逐固定 charge 使用：物理 Clifford 梯度和 OU 会连接相邻 charge。
应先把 central-sector CP-path 界与 binomial 外尾组装成全迹的 X、LX/p 界，再用 (7)。
若 X 中含 C_m theta_(v,m) 等逐 sector 系数，L 作用于这些 Q 函数的有限差分必须保留；
不能把它们当全 Fock 常数移过 L。突兀的 charge cutoff 同样可能产生额外导数项。

本次只查原文和写数学接口，没有重跑科学代码、没有改动已封存输出或旧来源文件。
