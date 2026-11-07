# Charged transfer：精确 body 表达、全次数尾与 whole-field 局域源

2026-11-07。由 `susskind-memory-oriented-transfer-seed.md` 继续。本文是本轮新证明草稿；不改写旧阶段证据，不把尚未完成的 charged 静态分离或 AA/AB/BA 核写成已完成。全程采用同一原始 Gaussian 耦合、standard CAR、归一 annealed full-Fock 迹，\(\lambda_N=p^2/N\in[\lambda_-,\lambda_+]\)、\(\lambda_->0\)，先固定所有次数与实时间紧窗。时间为 string time，\(\hbar=1\)。

已有输入来自 `susskind-memory-uniform-bound`、`susskind-memory-charged-transfer`、`susskind-memory-wick-rank`：Gaussian Sobolev 乘法与传播界、polar transfer 的 number-weighted defect 界、rank-filter 界。以下新的全时间论证只调用这些已写明的界；局域 source 的独立有限 CAR 检查由本轮其他新 driver 另行记录。

## 1. 精确接口与记号

令 \(h=H_m\)、\(h_-=H_{m-1}\)，\(\mathcal N\) 为**耦合 Gaussian chaos number**。对 neutral full-Fock lift，

\[
 DX=QX-\sum_i c_i^\dagger Xc_i,
 \qquad (DX)_m=mX_m-\sum_i c_i^\dagger X_{m-1}c_i.
\]

因此在 \(m>0\) 上，若 \(T_mX=m^{-1}\sum_i c_i^\dagger Xc_i\)，则

\[
 F_m(s)=T_m e^{ish_-}=e^{ish}-m^{-1}[D(e^{isH})]_m.
\]

Gaussian OU 与 sector restriction 可交换。取 active sectors \(m\ge p\) 的 \(\alpha_m=1-p/m\)，定义

\[
 E_m(s)=m[F_m(s)-\alpha_m^{\mathcal N}e^{ish}]
       =[m(1-\alpha_m^{\mathcal N})-D]e^{isH}\big|_m,
\]
\[
 G_{m,N}(s)=[p\mathcal N-D-\lambda_N\mathcal N(\mathcal N-1)]e^{isH}\big|_m.
\]

在 \(m<p\) 上相关 Hamiltonian 均为零，以上两个场均置零；\(m=0\) 同样置零。\(D\) 是 body/loss map，不能换成 output rank \(\mathcal R\)。所有涉及 \(D\) 的复合先在 full Fock 做，再 restrict。

使用 Sobolev 范数

\[
 \|X\|_r=\|(\mathcal N+1)^{r/2}X\|_2,
 \qquad \|X\|_2^2=\sum_m2^{-N}\binom Nm\mathbb E\operatorname{tr}_m X_m^\dagger X_m.
\]

## 2. 全电荷、固定导数及左右能量加权的强替换

对每个整数 \(k\ge0\)、\(u=p/m\in[0,1]\)，Taylor 积分余项给出**不限制 \(k\)** 的界

\[
0\le m[1-(1-u)^k]-kp+\binom{k}{2}\frac{p^2}{m}
\le\binom{k}{3}\frac{p^3}{m^2}.
\]

约定 \(\binom{k}{3}=0\) 于 \(k<3\)。令

\[
 R_m(k)=m(1-\alpha_m^k)-kp+\lambda_Nk(k-1).
\]

则 \(E_m-G_{m,N}=R_m(\mathcal N)e^{ish}\)，其中 \(D\) **精确相消**。在广中央区间 \(N/4\le m\le3N/4\)，

\[
 |R_m(k)|\le a_m(k+1)^3,
 \qquad
 a_m=\frac{4\lambda_+|m-N/2|}{N}+
            \frac{8\lambda_+^2}{3p}.
\]

故 \(\|R_m(\mathcal N)\phi\|_r\le a_m\|\phi\|_{r+6}\)。已有 Gaussian Sobolev 传播界以及左右 \(H\) 乘法界适用于
\(\partial_s^j e^{ish}=(ih)^je^{ish}\)，给出任意固定 \(r,j,a,b,S\) 的

\[
 \sup_{|s|\le S}
 \|h^a\partial_s^j(E_m-G_{m,N})h^b\|_{r,m}\le K a_m.
\]

这里并非从无导数误差裸求导：先对实际导数作相同 number-multiplier 估计。二项分布的精确恒等式
\(\sum_m b_m(m-N/2)^2=N/4\) 给出 \(\sum_m b_m a_m^2=O(p^{-2})\)。

广中央外的 active sectors 用 Bernoulli 不等式取得
\[
 |R_m(k)|\le (p+\lambda_+)(k+1)^2.
\]
全 sector 有 \(C_m\le M\)，故相同 Sobolev 界仍一致；二项外尾不超过 \(2e^{-N/8}\)。因此

\[
 \boxed{
 \sup_{|s|\le S}
 \|H^a\partial_s^j(E-G_N)H^b\|_r=O_{r,j,a,b,S}(p^{-1}).
 }
\]

若提前将 \(\lambda_N\) 换成 \(\lambda\)，另加
\(O(|\lambda_N-\lambda|)\) 误差，不能无条件保留 \(O(p^{-1})\) 速率。

## 3. 全次数 Taylor 尾：全部低 chaos 均保留

定义精确多项式系数

\[
 E_{n,m}=m[T_mh_-^n-\alpha_m^{\mathcal N}h^n],
 \qquad
 G_{n,m}=[p\mathcal N-D-\lambda_N\mathcal N(\mathcal N-1)]H^n\big|_m.
\]

它们包含 \(0,\ldots,n\) 中所有实际存在的 Gaussian chaoses；不是最高 chaos 的替代。设 \(L=2\sqrt M\)，并沿旧 polar transfer 记号令

\[
 B_m=\lambda_N+\lambda_N^2+
        \lambda_N(2m-N-1)^2/N,
 \qquad T_m=S_+\mathcal U.
\]

旧 number-weighted polar defect 在 degree 不超过 \(a\) 时给出

\[
 \|G_{\rm pol}\psi_{\le a}\|_2
 \le K\frac{\sqrt{MB_m}}N(a+1)^{5/2}\|\psi\|_2,
 \quad G_{\rm pol}=\mathcal U H_{-,L}-H_{+,L}\mathcal U.
\]

Telescoping 展开 \(\mathcal Uh_-^n-h^n\) 的第 \(a\) 项中，前面的 \(h_-^a\) 范数至多 \(L^a\sqrt{a!}\)，defect 后的 degree 至多 \(a+1\)，剩余乘法给
\(L^{n-1-a}\sqrt{n!/(a+1)!}\)。因此

\[
 \|\mathcal Uh_-^n-h^n\|_2
 \le K\frac{\sqrt{MB_m}}N L^{n-1}\sqrt{n!}
              \sum_{a=0}^{n-1}(a+1)^2.
\]

乘 \(m\le N\)，并使用 \(S_+\) 的 contraction 性。其余 filter 部分由旧三项式

\[
 \|(S_+-\alpha_m^{\mathcal N})\phi\|_2
 \le\frac4N\|\mathfrak D\phi\|_2+
 \frac{8p|2m-N-1|}{N^2}\|\phi\|_2^{\rm Sob}+
 \frac{528p^2}{N^2}\|\phi\|_4^{\rm Sob},
 \quad \mathfrak D=p\mathcal N-\mathcal R,
\]

以及
\(\|\mathfrak D h^n\|_2\le K n(n+1)L^n\sqrt{n!}\)
控制。degree 不超过 \(n\) 使后两项分别只增加 \(n+1\) 与 \((n+1)^2\)。又有
\(p|2m-N-1|/N\le\sqrt{B_m}\)、\(\lambda_N\le\sqrt{B_m}\)，所以

\[
 \|E_{n,m}\|_2+\|G_{n,m}\|_2
 \le K\sqrt{B_m}(n+1)^3L^n\sqrt{n!}
\]

在广中央一致成立；number remainder 被同型界吸收。全二项平均满足
\(\sum_m b_m B_m=2\lambda_N+\lambda_N^2+\lambda_N/N\)。外中央用
\(\|E_{n,m}\|_2\le2N L^n\sqrt{n!}\) 及上述粗 number remainder，指数尾吸收所有 \(N,p\) 因子。因此在 full trace 上

\[
 \boxed{\|E_n\|_2+\|G_n\|_2\le K(n+1)^3L^n\sqrt{n!}.}
\]

固定 Sobolev 阶和左右能量次数只再增加 \(n\) 的固定幂。于是

\[
 G_N(s)=\sum_{n\ge0}\frac{(is)^n}{n!}G_{n,N}
\]

及任意所需的固定时间导数、左右能量加权，在每个实时间紧窗绝对一致可和。\(\sqrt{n!}\) 包络配 \(1/n!\) 足以控制次数尾。

## 4. 实际 charged force 的完整级数

精确微分接口是

\[
 \mathscr Q_\alpha f=-f''+i\alpha\{h,f'\}+\alpha^2hfh,
 \quad m(\Lambda-\Lambda_{\rm OU})=\mathscr Q_\alpha E.
\]

上一节的加权强替换与中央 \(|1-\alpha_m|=O(p^{-1})\)、外部二项尾，给出
\(\mathscr Q_{\alpha_m}E_m-\mathscr Q_1G_{m,N}\to0\) 的全迹强替换。可直接定义

\[
 K_{n,N}=G_{n+2,N}-\{H,G_{n+1,N}\}+HG_{n,N}H,
\]
\[
 \mathscr Q_1G_N(s)=\sum_{n\ge0}\frac{(is)^n}{n!}K_{n,N},
 \qquad
 \|K_{n,N}\|_2\le K'(n+3)^4L^{n+2}\sqrt{n!}.
\]

因此实际 \(A\) 可由此级数双侧乘 \(e^{-ihs/2}\)，再作原来的精确 pinching。\(G_0=G_1=0\)，所以 \(K_0=G_2\)。这条构造避免了尚未证明的
\((p\mathcal N-D)(Z_j-Q_j)\to0\)，也避免把依赖 charge 的 Jacobi 系数误当成 \(D\) 的常数。

## 5. Whole-field nilpotent source：先消 singleton，再取 pair-Poisson 极限

本节处理固定原始次数的矩。每个完整 field interval 只放一个 body mark，分别引入可交换 \(z_a\)，且 \(z_a^2=0\)。\(D\) 与 Gaussian number 可交换，因此

\[
 e^{z_a(p\mathcal N-D)}H^{n_a}
       =e^{z_ap\mathcal N}e^{-z_aD}H^{n_a}.
\]

固定 Gaussian pairing 后，该 field 的内部 pairs 数为 \(i_a\)，开放 degree 是
\(d_a=n_a-2i_a\)；number 因子正是 \(e^{z_apd_a}\)。外部 energy slots 不计入 \(n_a\)。

### 5.1 对称 loss 与 graded 局域实现

在 neutral full-Fock 空间，

\[
 DX=\sum_i\left(\tfrac12\{n_i,X\}-c_i^\dagger Xc_i\right).
\]

它与非对称局域表 \(n_iX-c_i^\dagger Xc_i\) 的差为
\(\tfrac12[Q,X]=0\)。对全局偶的 graded word，局域 sandwich 穿过其余 sites 带
\((-1)^{\operatorname{parity}(x_i)}\)；在上述对称形式下得到

\[
 D_iI=0,\quad D_ic_i=\tfrac12c_i,\quad
 D_ic_i^\dagger=\tfrac12c_i^\dagger,\quad
 D_in_i=n_i,\quad D_i(1-n_i)=-n_i.
\]

故局域源可直接实现为

\[
 c,c^\dagger\mapsto e^{-z_a/2}(c,c^\dagger),\quad
 n\mapsto e^{-z_a}n,\quad
 1-n\mapsto1-e^{-z_a}n.
\]

这些 maps 保持每个局域 field 的奇偶；neutral 条件也保持。非对称 creator-count 表适合纯 Wick 子块，但 whole-field normalization 应采用此对称实现，或者严格保留全局 charge-commutator 相消。

### 5.2 每条弦的 singleton 相消

一条 Gaussian chord 携带两个独立 \(p\)-sets。若它在某 field 中只有一个端点，则其两个 sets 在 singleton site 上各给 \(e^{-z_a/2}\)。全部 open chords 的这些因子，恰与 \(e^{z_apd_a}\) 抵消。

若该 chord 两端都在同一 field 内，其两个 sets 的局域词分别为 \(n\) 和 \(1-n\)。除以未加源的 singleton trace 后，两因子为

\[
 e^{-z_a},\qquad 2-e^{-z_a}.
\]

它们各出现 \(p\) 次，而

\[
 e^{-z_a}(2-e^{-z_a})=1\pmod {z_a^2}.
\]

所以完整 singleton product 在 nilpotent 环中**精确等于 1**，包括不同 sources 的 mixed 系数。这里不能先取一个真实 \(t\) 后丢弃 \(p t^2\) 项；使用每 field 一次 mark 的 square-zero 环正是所需范围。

### 5.3 Pair source 必须除 source-dependent singleton factors

对 creator sets \(u,v\) 的一个二重碰撞，将同一 site 在每个 field interval 中的矩阵词先按原顺序相乘，再施加该 field 的局域 loss map；最后按全词顺序取 \(\operatorname{Tr}/2\)。其归一化比 \(R_{uv}(z)\) 的分母，是这两个 sets 各自未合并时的 **source-dependent** singleton factors，而不是只取零源的常数分母。

Jordan/Koszul 符号必须在 trace 前计算为合并与未合并配置的相对符号。局域源保持 parity，故该符号不随 mark 子集改变；零的 unmarked numerator 不允许反推符号。分母的零源常数非零，其 nilpotent 逆是有限代数式。

若全图有 \(\ell\) 条 Gaussian chords，标准 pair-Poisson 极限给

\[
 \exp\left\{\lambda\left[\sum_{u<v}(R_{uv}(z)-1)-\ell\right]\right\}.
\]

零源指数为 \(-4\lambda\operatorname{cr}(\pi)\)。\(-\ell\) 是原来的 \(C_N^\ell\) variance normalization。由于 singleton 大项已经精确消掉，剩下每个固定 nilpotent 系数只带有限的局域权重；三重交叠概率 \(O(p^3/N^2)\) 和旧 overlap 指数矩足以控制极限。

写括号为
\(\kappa_0+\kappa_1z_1+\kappa_2z_2+\kappa_{12}z_1z_2\)。\(G_n\) 还需减掉 number 项 \(\lambda_Nd(d-1)\)。因此 one-/two-\(G\) 图权重分别为

\[
 C^\ell q^{\operatorname{cr}(\pi)}
       \lambda[\kappa_1-d_1(d_1-1)],
\]
\[
 C^\ell q^{\operatorname{cr}(\pi)}
 \left\{\lambda\kappa_{12}
 +\lambda^2[\kappa_1-d_1(d_1-1)]
             [\kappa_2-d_2(d_2-1)]\right\}.
\]

这里 \(d_a\) 依赖该图的内部 Gaussian pairs，不能用原始 slot 数 \(n_a\) 代替。

## 6. G2 的零 chaos 背景：有限恒等式与两个独立极限计算

令

\[
 C_m=w_N\binom mp\binom{N-m+p}p,
 \qquad
 w_N=4^p\frac{J_s^2}{\lambda_N}\left(\frac{p!}{N^p}\right)^2.
\]

由于 \(H^2=Z_2+C_Q\)，

\[
 G_2=(2p-D)Z_2-2\lambda_NZ_2-DC_Q,
 \qquad \WW_0G_2=-DC_Q.
\]

精确相邻比和 body 导数为

\[
 \frac{C_{m-1}}{C_m}
 =\frac{m-p}{m}\frac{N-m+p+1}{N-m+1},
\]
\[
 (DC_Q)_m=m(C_m-C_{m-1})
 =C_m\frac{p(N-2m+p+1)}{N-m+1}.
\]

对 \(m=N/2+xp\) 的紧窗，\(C_m\to C\)，所以

\[
 \WW_0G_2\longrightarrow-2\lambda C(1-2x).
\]

在全二项平均中，\(X_N=(m-N/2)/p\) 收敛到方差 \(1/(4\lambda)\) 的中心 Gaussian。广中央区间内上述精确式被 \(K(1+|X_N|)\) 控制，外区间只有多项式增长而二项尾指数小；因此前二矩以及相应一致可积性成立：

\[
 \lim\mathbb E\operatorname{tr}\WW_0G_2=-2\lambda C,
\]
\[
 \lim\|\WW_0G_2\|_2^2=4C^2(\lambda+\lambda^2),
 \qquad
 \lim\operatorname{Var}(\WW_0G_2)=4\lambda C^2.
\]

同一结果可完全从 nilpotent 局域图得到。一个内部 \(H^2\) 弦的 own creator/annihilator overlap 给

\[
 R_{\rm own}(z)
 =\frac{\tau[e^{-z}n]}{\tau[e^{-z}n]\tau[1-e^{-z}n]}
 =2(1-z).
\]

两个 \(H^2\) fields 各作内部配对时，跨两内部弦的四个 set 碰撞用
\(R_{XY}=\tau[L_1(X)L_2(Y)]/[\tau L_1(X)\tau L_2(Y)]\)，
\(X,Y\in\{n,\bar n\}\)、\(L_a(n)=(1-z_a)n\)、\(L_a(\bar n)=\bar n+z_an\)，给

\[
\begin{array}{c|c}
(X,Y)&R_{XY}\\\hline
(n,n)&2\\
(\bar n,\bar n)&2-2z_1-2z_2+4z_1z_2\\
(n,\bar n)&2z_2\\
(\bar n,n)&2z_1
\end{array}
\]

总 source 指数是
\(\lambda[-2z_1-2z_2+4z_1z_2]\)，其 mixed 系数正是
\(4\lambda+4\lambda^2\)。该图没有任何 Gaussian 弦直接连接两个 fields，却仍有非零 connected 静态涨落。

## 7. Whole-field 闭式系数：内部弦、开放弦与 spectator 分类

以下 \(O_a\) 专指在 field \(a\) 中**恰有一个端点**的 Gaussian 弦集；\(i_a\) 是两端均在该 field 内的弦数，\(d_a=|O_a|\)。令 \(B=O_1\cap O_2\)，\(b=|B|\)。局域分类给出

\[
 \boxed{\kappa_a=d_a(d_a-1)-4\operatorname{cr}(O_a)-2i_a,}
\]
\[
 \boxed{\kappa_{12}=2b(b-1)-4\operatorname{cr}(B)+4i_1i_2.}
\]

它们是任意允许 pairing 的逐图公式，不限于无 field bridge 的情形。源变量定义仍是第 5 节的 \(p\mathcal N-D\)，所以应用于 \(G\) 时须继续减 \(d_a(d_a-1)\)：

\[
 \kappa_a^G=-4\operatorname{cr}(O_a)-2i_a.
\]

单标记分类如下。

- 同一 internal chord 的两个 sets 碰撞贡献 \(-2\)；同一 open chord 的 own 碰撞，源的 odd 缩放与分母抵消，贡献零。
- 两条 open incident chords 给已核 pure-oriented 权重；全部相加为 \(d(d-1)-4\operatorname{cr}(O)\)。
- 两条 internal chords 若相交，合并 site 的完整四字词为零，loss map 不会复活零词；若不相交，两个非零碰撞的 linear coefficients 为 \(+2,-2\)，相消。
- internal--open 若相交，field 内三字词为零；若不相交，两个非零项同样以 \(+2,-2\) 相消。此交叉仍进入普通 \(q^{\operatorname{cr}(\pi)}\)，但不进入 \(\kappa_a\)。
- 与 field 外 spectator chord 的碰撞：open--spectator 逐项只有局域 odd 缩放，与 singleton 分母完全消去；internal--spectator 的四个 set 碰撞总 linear coefficient 为零。因而单标记闭式不依赖第二 field 的内部配对或其余 spectator pairing。

双标记分类中，除了 internal--internal 和两条共同桥弦，其余情形至少有一个 marked field 只见局域 odd singleton，其 source 在比值中约掉。每对两 field 各自的 internal chords，由第 6 节四项表给 mixed coefficient 4；共同桥弦给已核 pure-DD 权重 \(2b(b-1)-4\operatorname{cr}(B)\)。因此得到上面的闭式。特别在没有直接 field bridge 时仍有 \(\kappa_{12}=4i_1i_2\)，而非零。

这项分类说明：Gaussian 内部 contractions 的交叉因子 \(q^{c_\gamma}\) 不受单 body source 微分。不能另加一个未经推导的 \(\partial_q(q^{c_\gamma})\)。

## 8. 联合 charge source 与下一轮回归入口

令 \(X_N=(Q-N/2)/p\)，在 full trace 插入 \(e^{tX_N}\)。它的未标记归一化满足

\[
 \operatorname{tr}e^{tX_N}
 =\cosh(t/(2p))^N\longrightarrow e^{t^2/(8\lambda)}.
\]

单 site 的倾斜 occupation 为
\(\pi_t=(1+\tanh(t/(2p)))/2\)。open-field 的 odd 缩放仍精确约掉。一个 internal chord 的两 singleton source factors 相乘为

\[
 1+z\frac{2\pi_t-1}{1-\pi_t},
 \qquad
 \left[1+z\frac{2\pi_t-1}{1-\pi_t}\right]^p\longrightarrow e^{tz}.
\]

所以每个 internal chord 增加 joint exponent \(tz\)。倾斜使固定 pair-overlap 比值只变 \(O(1/p)\)；固定 \(t\) 紧窗上的旧 overlap 指数矩给出相同的支配。对两个 body sources，固定图的 joint source 因而由

\[
 e^{t^2/(8\lambda)}
 \exp\{\lambda(\kappa_0+\kappa_1z_1+\kappa_2z_2+
                       \kappa_{12}z_1z_2)
             +t(i_1z_1+i_2z_2)\}
\]

描述（省略原有 \(C^\ell\)）。这识别了 \(4\lambda iX\) 的 charge-linear 部分：\(\operatorname{Var}X=1/(4\lambda)\)，其交叉方差正是 \(4\lambda i_1i_2\)。这里的 \(t\) 是 charge source，不是物理时间。

对下一轮可保存一个明确的回归入口。设

\[
 x^n=\sum_d a_{n,d}H_d^{(q,C)}(x),\qquad
 a_{0,0}=1,\quad
 a_{n+1,d}=a_{n,d-1}+C[d+1]_q a_{n,d+1}.
\]

固定内部 matching \(\gamma\) 后，\(i=(n-d)/2\)，其普通权重为 \(C^i q^{c_\gamma}\)。开放端点置换的 Mahonian 平均给
\(\mathbb E_q\operatorname{inv}=q\partial_q\log P_d\)。因此候选 joint charge--energy regression 为

\[
 g_n(X,H)=\sum_d a_{n,d}
 \left[-4\lambda q\partial_q\log P_d
       -\lambda(n-d)(1-2X)\right]H_d^{(q,C)}(H).
\]

它在 \(n=2\) 恢复
\(-4\lambda q H_2/(1+q)-2\lambda C+4\lambda CX\)。由于闭式 source 不含 internal-cross 微分，\(a_{n,d}\) 外没有额外 \(q\)-导数。

更进一步的有限次数残量候选是
\(\sum_d a_{n,d}J_d\)，其中
\[
 J_d=(dp-D)Z_d-
 \lambda[d(d-1)-4q\partial_q\log P_d]H_d^{(q,C)}(H).
\]
这由 source 闭式提示，但本轮仍不把其 joint spectral/strong Gram 识别以及全时间 pinching 写成已完成。下一轮须将 charge source 的混合矩、开放端点 grouping 与无桥抵消整理成完整证明，再处理 AA/AB/BA。

本稿可靠交付范围：全加权 \(E\to G\) 强替换、包含所有低 chaos 的完整 Taylor/force 级数及统一尾、固定次数 whole-field nilpotent source 的 singleton 归一方式与闭式系数、G2 charge 背景及前二矩。旧数据与稿件保持；本轮新独立计算的参数、哈希和计数另见后续 evidence。
