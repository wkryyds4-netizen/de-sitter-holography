# Charged memory：联合电荷回归、共同 J 强归约与全时间谱

2026-11-07。本稿接续已封存的 [whole-field 源证明](susskind-memory-oriented-transfer-proof-draft.md)，将其中第 8 节的候选补成证明。这里的 J 在“共同 J”中指同一组原始 Gaussian 耦合；下文带下标的 \(J_d\) 是新定义的中心化算符场，两者不可混淆。旧文件、旧科学输出不改写。

采用 standard CAR、overlapping complex \(p\)-body Hermitian Gaussian Hamiltonian；态为 \(\beta=\mu=0\) 的归一 annealed full-Fock 迹。令偶 \(N\to\infty\)、\(p\to\infty\)、\(\lambda_N=p^2/N\to\lambda>0\)，固定 \(J_s\) 与 string time，\(\hbar=1\)。所有场次数先固定；有限时间结论在任意固定实时间紧窗成立，不包含 \(\lambda\to0\) 或随 \(N\) 增长的时间。

本稿调用的已证输入是：whole-field 单/双源逐图公式、全部低 chaos 的 ordinary Taylor 尾、固定 Gaussian degree 的左右 \(H\) 乘法界，以及 [neutral 后代阶段](susskind-memory-wick-descendants.md) 的中性残量谱与级数尾。以下没有将一个未经控制的 \(D\) 作用在 spectral/Wick 近似误差上。

## 1. 定义与待分离的守恒背景

记
\[
 C=J_s^2/\lambda,\qquad q=e^{-4\lambda},\qquad
 X_N=(Q-N/2)/p,\qquad \mathcal N=\text{coupling Gaussian chaos number}.
\]
\(D\) 是 full-Fock neutral 算符的 body/loss map，
\[
 DY=QY-\sum_i c_i^\dagger Yc_i,
\]
不是 output rank。定义
\[
 G_{n,N}=[p\mathcal N-D-\lambda_N\mathcal N(\mathcal N-1)]H_N^n.
\]
用 \(H_d^{(q,C)}\) 表示 monic \(q\)-Hermite 多项式，必要时简写为 \(H_d\)：
\[
 H_{d+1}(x)=xH_d(x)-C[d]_qH_{d-1}(x),\qquad
 \|H_d\|_{\rho_q}^2=C^dP_d(q),\qquad P_d(q)=[d]_q!.
\]
另设
\[
 m_d=q\partial_q\log P_d(q),\qquad m_0=m_1=0,
\]
\[
 x^n=\sum_d a_{n,d}H_d(x),\qquad
 a_{0,0}=1,\quad
 a_{n+1,d}=a_{n,d-1}+C[d+1]_qa_{n,d+1}.
\]
系数非负；只在 \(0\le d\le n\)、\(d\equiv n\pmod2\) 非零。

结论中的联合回归为
\[
 \boxed{g_n(X,H)=\sum_d a_{n,d}
 [-4\lambda m_d-\lambda(n-d)(1-2X)]H_d(H).}\tag{1}
\]
有限 \(N\) 反项使用同一个固定极限多项式 \(g_n(X_N,H_N)\)。它与每个样本的 \(H_N\) **精确对易**，因为 \([Q,H_N]=0\)。这里不声称它是有限 \(N\) 的精确条件期望。

## 2. 联合 charge tilt 及一致可积性

插入 \(e^{tX_N}\)，其中 \(t\) 是电荷源，不是物理时间。未加任何场时
\[
 \operatorname{tr}e^{tX_N}=\cosh(t/(2p))^N
 \longrightarrow e^{t^2/(8\lambda)}.
\]
倾斜的单 site 占据率为
\[
 \pi_t=\frac{1+\tanh(t/(2p))}{2}.
\]
一个 internal Gaussian chord 的两个 singleton 分别为 \(n\)、\(1-n\)。在该 field 的 square-zero 源 \(z^2=0\) 下，除以各自零源迹后的乘积为
\[
 1+z\frac{2\pi_t-1}{1-\pi_t}
 =1+z(e^{t/p}-1).
\]
它出现 \(p\) 次，因此
\[
 [1+z(e^{t/p}-1)]^p
 =1+p(e^{t/p}-1)z\longrightarrow e^{tz}.
\]
open chord 的 odd singleton 缩放仍与 \(e^{pdz}\) 精确抵消。故每条 field-internal chord 增加 \(tz\)，而不改变单源对内部交叉的处理。

下面说明这不是仅在典型电荷区间的形式计算。固定全部原始 slot 数，令 \(|t|\le T\)。相对 \(t=0\)，未标记 singleton 的归一因子是
\[
 [4\pi_t(1-\pi_t)]^{p\ell}\longrightarrow1,
\]
其中 \(\ell\) 为 Gaussian chord 数。各 pair-overlap ratio 的分母零阶项在充分大 \(p\) 时一致远离零，ratio 的变化为 \(O_T(p^{-1})\)。先完成精确 singleton 相消后，每个所需 nilpotent 系数由
\(K_T(1+K_{\rm ov})^r e^{cK_{\rm ov}}\) 控制，\(r,c\) 只依赖固定 slot 数，\(K_{\rm ov}\) 是这些固定 \(p\)-sets 的总交叠数。前驱的 hypergeometric overlap 指数矩界提供统一支配；三重交叠仍以 \(O(p^3/N^2)\) 消失。全局的 \(\cosh(t/(2p))^N\) 在紧 \(t\) 窗一致有界。

同一估计适用于加入任意固定左右能量幂以及固定线性组合后的正范数。特别在 \(t=\pm T\) 给出电荷指数加权二矩的一致有界，从而控制每个固定 \(X_N^r\) 的尾。也可先在复 \(t\) 紧集用上述局部解析表达及 \(\operatorname{Re}t\) 支配得到局部一致解析收敛，再用 Cauchy 公式取导数。两种方式均证明全部固定 charge-power 混合矩收敛；没有从未加权矩收敛直接交换微分。

令一张两场图中 field \(a\) 的 internal chord 数为 \(i_a\)，open chord 集为 \(O_a\)，\(d_a=|O_a|\)。设
\[
 B=O_1\cap O_2,\quad b=|B|,\quad
 k_a=-4\operatorname{cr}(O_a)-2i_a,
\]
\[
 K_D(B)=2b(b-1)-4\operatorname{cr}(B).
\]
whole-field 源的 joint tilt 两标记权重为
\[
 C^\ell q^{\operatorname{cr}\pi}e^{t^2/(8\lambda)}
 \{\lambda[4i_1i_2+K_D(B)]
 +(\lambda k_1+i_1t)(\lambda k_2+i_2t)\}.\tag{2}
\]
若 \(X\sim\gamma_\lambda=N(0,1/(4\lambda))\)，则 (2) 等价于在 \(X\) 上取平均的条件规则：
\[
 \text{single}_a=\lambda[k_a+4i_aX],\qquad
 \text{connected}_{12}=\lambda K_D(B).\tag{3}
\]
确切地，\(\mathbb E[Xe^{tX}]=(t/(4\lambda))e^{t^2/(8\lambda)}\)，而
\(16\lambda^2 i_1i_2\operatorname{Var}X=4\lambda i_1i_2\)。所以旧的无桥项 \(4i_1i_2\) 被识别为一个共同守恒电荷的协方差，并未被丢掉。

未标记图给出联合标量律
\[
 (X_N,H_N)\ \Longrightarrow\ \gamma_\lambda\otimes\rho_q
\]
及全部固定混合矩。这里说的是 commuting charge/energy 在该态下的联合谱律，不是任意动态场的概率独立性。

## 3. 切图证明联合正交回归

固定 \(H^n\) field 中的内部 partial matching \(\gamma\)，令 \(i=|\gamma|\)、\(d=n-2i\)。去掉内部弦后留下 \(d\) 个开放槽。该组的固定因子为
\[
 C^i q^{c_\gamma},
\]
其中 \(c_\gamma\) 包含 internal--internal 及 internal--open 的全部交叉。给定外部伙伴位置后，开放端点置换的相对权重为 \(q^{\operatorname{inv}\sigma}\)。规则 (3) 在这一组中变成
\[
 \lambda[-4\operatorname{inv}\sigma-2i+4iX].
\]
其 Mahonian 平均是
\[
 \lambda[-4m_d-2i+4iX].\tag{4}
\]
对 \(\gamma\) 求和，\(\sum_\gamma C^iq^{c_\gamma}=a_{n,d}\)，得到 (1)。

该计算也直接给出所有 \(X^rH_k(H)\) 测试下的正交关系。外部 spectral \(q\)-Hermite 的标准 Wick cancellation 只留下 \(k=d\) 的开放匹配，内积恰为 (4) 乘 \(C^dP_d\)。在这一线性内积内使用已证 \(Z_d-H_d(H_N)\to0\) 是合法的；没有对这条误差再施加 \(D\)。Gaussian charge 多项式和紧支撑 energy 多项式在 \(L^2(\gamma_\lambda\otimes\rho_q)\) 稠密，故
\[
 g_n=\operatorname{Proj}_{L^2(X,H)}G_n.\tag{5}
\]
特别注意，(4) 的 mark 不作用于 \(q^{c_\gamma}\)。对 \(a_{n,d}\) 再取 \(q\)-导数会加入本模型没有的项。

低阶为
\[
 g_0=g_1=0,\qquad
 g_2=-\frac{4\lambda q}{1+q}H_2-2\lambda C+4\lambda CX.
\]
这同时恢复 energy 回归和旧 \(G_2\) 的零 chaos 电荷背景。

## 4. 共同原始耦合中的强 Gram 归约

令 \(Z_{d,N}=\mathbb W_d(H_N^d)\) 为最高 coupling chaos，并定义
\[
 W^-_{d,N}=(dp-D)Z_{d,N},\quad
 \kappa_d^- =\lambda[d(d-1)-4m_d],
\]
\[
 J_{d,N}=W^-_{d,N}-\kappa_d^-H_d(H_N),\qquad
 R_{n,N}=G_{n,N}-g_n(X_N,H_N).\tag{6}
\]
\(J_0=J_1=0\)。所有这些对象都在同一原始 \(J\) 样本上定义。

对一对固定的内部 matchings \(\gamma_1,\gamma_2\)，展开 \(R_n\)--\(R_m\) 的四项。用 (3) 完成 charge 平均，再用 (4) 减去各自回归。其剩余单源系数为
\[
 -4\lambda[\operatorname{cr}(O_a)-m_{d_a}],
\]
connected 系数为 \(\lambda K_D(B)\)。所有 \(i_a\) 和 \(X\) 依赖消失。恰好得到 pure Wick fields \(J_{d_1}\)、\(J_{d_2}\) 的中心化源规则。内部弦只留下原来的 \(C^{i_a}q^{c_{\gamma_a}}\)。

同样的计算还适用于 **mixed common-copy** 的 \(R_n\)--\(J_e\) 内积：第二场取 \(i_2=0\)，并施加它的纯 Wick 禁内部配对规则即可。故任意固定外部能量次数下，三种 Gram 项
\[
 \langle R_n,R_m\rangle,\qquad
 \left\langle R_n,\sum_ea_{m,e}J_e\right\rangle,\qquad
 \left\langle\sum_da_{n,d}J_d,\sum_ea_{m,e}J_e\right\rangle
\]
都有同一个图极限。展开同一概率空间中的差范数，得到
\[
 \boxed{\left\|R_{n,N}-\sum_da_{n,d}J_{d,N}\right\|_2\longrightarrow0.}\tag{7}
\]
这一步依赖 mixed common-copy 项，不能仅凭两个独立场族的边缘 Gram 律相同推出。

差场的 coupling Gaussian degree 有固定上界；\(X_N\) 对 coupling 变量是确定性算符。因此既有固定-degree 左右乘法界允许在 (7) 两侧加固定 \(H_N^a,H_N^b\)，仍强收敛。这里没有使用“任意 \(L^2\) 收敛可乘无界算符”的错误规则。

## 5. Pure 残量的范数与无桥相消

记 \(v_d=(q\partial_q)^2\log P_d=\operatorname{Var}_q(\operatorname{inv})\)。等时两场的弦全部是 direct bridges，故直接由中心化源规则得到
\[
 \boxed{\|J_d\|_2^2=C^dP_d
 \{\lambda[2d(d-1)-4m_d]+16\lambda^2v_d\}.}\tag{8}
\]
这里是有限 \(N\) 范数的极限。因为
\(0\le\operatorname{inv}\le d(d-1)/2\)，
\[
 \|J_d\|_2\le K(\sqrt\lambda\,d+\lambda d^2)
                 \sqrt{C^dP_d}.\tag{9}
\]
在 \(d=2\)，(8) 给出
\[
 \|J_2\|_2^2=C^2\left[4\lambda+
                  \frac{16\lambda^2q}{1+q}\right],
\]
恢复旧 \(Y_-\) 的动态范数；没有留下额外的 \(4\lambda C^2\) charge 静态方差。

现在在 \(J_d^\dagger H^aJ_eH^b\) 图中固定外部伙伴位置。若两 field 没有 direct Gaussian bridge，则 \(B=\varnothing\)、\(K_D(B)=0\)。两个 field 内的开放置换独立；在一个连续 field 区间中置换其端点，不改变两端均在区间外的弦与其相交的数目。因此各自中心化 Mahonian 和为零，整组无桥图消失。单张图不必消失。

等价地，对原 \(R_n\) 先在联合 \(X,H\) 上回归，再以固定内部 matching 为组做同样的抵消。(7) 提供更简洁的有限线性组合归约。

## 6. 两能量谱的 Hilbert--Schmidt 界

令 \(p_r=H_r/\sqrt{C^rP_r}\)，对固定 \(J_d,J_e\) 定义
\[
 M_{mn}^{de}=\lim_N\mathbb E\operatorname{tr}
 [J_{d,N}^\dagger p_m(H_N)J_{e,N}p_n(H_N)].
\]
外部两个 \(q\)-Hermite 块取消各自内部配对。被删除的 energy-internal chord 不属于任何 marked incident set；删除它不改变剩余弦之间的交叉或 mark 权重，故标准 multiblock Wick cancellation 在带标记图中仍成立。

取 \(D=d+e\)。无桥组相消后，每张剩余图至少有一条 field bridge。若 \(r\) 条弦连接两个外部 energy 块，端点计数给
\[
 |m-n|\le D,\qquad r\ge\min(m,n)-D.
\]
在循环块顺序 field--energy--field--energy 下，每条 field bridge 与每条 energy bridge 相交。故至少产生 \(q^r\)。固定场端点的选择至多贡献 \(K_D(m+n+1)^D\)，剩余 energy bridge 置换的权重和为 \(P_r(q)\)。又有
\[
 \frac{P_r}{\sqrt{P_mP_n}}\le(q;q)_\infty^{-1}
 \quad(r\le\min(m,n)).
\]
场 mark 和回归系数由固定 \(d,e\) 控制，因而
\[
 \boxed{|M_{mn}^{de}|\le K_D(m+n+1)^D
 q^{\max(0,\min(m,n)-D)}\,
 \mathbf1_{|m-n|\le D}.}\tag{10}
\]
此界对外部 \(m,n\) 可平方求和。常数可以依赖固定场、\(C,\lambda,q\)，没有声称它对场次数或 \(\lambda\to0\) 一致。

有限 \(N\) 谱测度约定为
\[
 \nu_{de,N}(A,B)=\mathbb E\operatorname{tr}
 [J_{d,N}^\dagger P_{H_N}(B)J_{e,N}P_{H_N}(A)],\tag{11}
\]
其中 \(E\) 是 input/column，\(E'\) 是 output/row。对有限场列表，这是正矩阵测度。全部固定两能量矩由上述图规则给出；固定 mark 的绝对值界、高矩一致可积性和正性给出弱极限，且支撑于
\[
 [-R_H,R_H]^2,\qquad R_H=2\sqrt{C/(1-q)}.
\]
可显式补出支撑步骤：固定场 degree 后，所有 mark/回归系数的绝对值有仅依赖该 degree 的上界；去掉 Wick 禁配对只增加普通非负 q-弦图和。因此两个坐标分别满足
\[
 \int E^{2m}\,d\nu_{dd}\le K_dR_H^{2m},\qquad
 \int(E')^{2m}\,d\nu_{dd}\le K_dR_H^{2m}.
\]
固定场附加次数的 \(R_H\) 幂已吸收入 \(K_d\)，它不依赖外部 \(m\)。任意 \(\epsilon>0\) 的谱尾由 \(K_d[R_H/(R_H+\epsilon)]^{2m}\) 控制，令 \(m\to\infty\) 即为零；正矩阵测度把对角支撑扩到交叉项。矩的唯一性来自紧支撑。由 (10)，其密度为
\[
 f_{de}(E,E')=\sum_{m,n\ge0}M_{mn}^{de}p_n(E)p_m(E'),
 \qquad d\nu_{de}=f_{de}\,d\rho_q(E)d\rho_q(E'),\tag{12}
\]
级数在 \(L^2(\rho_q\otimes\rho_q)\) 收敛。它与原测度有相同多项式矩，故是同一个测度。\(\rho_q\) 无原子，所以能量对角线 \(\Delta=\{E=E'\}\) 的质量为零。

任意固定 \(R_n\) 通过 (7) 继承此结论，固定左右能量多项式也继承。对 neutral 阶段的残量，下文改记为 \(N_j\) 以免与 \(R_n\) 冲突；其对角谱已证绝对连续。共同源图给出 combined \((J_d,N_j)\) 的全部 mixed moments，有限 \(N\) 正性随极限保留。对任意 Borel 集 \(S\)，谱 Cauchy--Schwarz 给
\[
 |\nu_{d,j}(S)|^2\le\nu_{d,d}(S)\nu_{j,j}(S),
\]
更一般地给出相应 variation 的控制；因此 mixed 谱也绝对连续。其矩阵 density 正性进一步给出
\[
 |f_{ab}(E,E')|^2\le f_{aa}(E,E')f_{bb}(E,E'),\qquad
 \|f_{ab}\|_2^2\le\|f_{aa}\|_2\|f_{bb}\|_2.
\]
各固定 charged/neutral 场的对角 density 已在 \(L^2\)，所以 combined 固定场 matrix density 的每个条目都在 \(L^2\)，也都有 (12) 形式的 Hermite 展开。存在性来自共同源矩；绝对连续性和这个 \(L^2\) 推论不要求另猜独立的 probe 模型。

mixed 单图中应使用 whole charged field 的 \(O_D\) 与较小 neutral marked Wick subblock 的 incident 集 \(I_L\)，
\[
 B_{DL}=O_D\cap I_L,\qquad
 K_{DL}=|B_{DL}|(|B_{DL}|-1)-4\operatorname{cr}(B_{DL}).
\]
\(I_L\) 可能包含 neutral 较大 composite field 的内部弦；不能用该大 field 的 open 集替代。对 pure \(J_d\)，charged internal 数为零，故 mixed source 不需要额外的 charged internal 背景项。

## 7. 精确 pinching 的传递

定义完整退化能量块的精确 pinching
\[
 \mathcal P_HY=\sum_eP_eYP_e.
\]
由 (11) 与无对角质量，闭集 \(\Delta\) 上的 Portmanteau 不等式给
\[
 \limsup_N\|\mathcal P_{H_N}R_{n,N}\|_2^2
 =\limsup_N\nu_{R_nR_n,N}(\Delta)
 \le\nu_{R_nR_n}(\Delta)=0.\tag{13}
\]
没有将无限时间平均与 double scaling 极限交换。因为 \(g_n(X_N,H_N)\) 精确守恒，
\[
 (1-\mathcal P_H)G_{n,N}
   =R_{n,N}-\mathcal P_HR_{n,N}.
\]
故联合电荷回归后的动态残量正是实际 pinching 后的强极限。

实际 charged force 使用 ordinary 系数
\[
 K_{n,N}=G_{n+2,N}-\{H_N,G_{n+1,N}\}+H_NG_{n,N}H_N.
\]
设
\[
 k_n(X,H)=g_{n+2}-2Hg_{n+1}+H^2g_n,
\]
\[
 L_{n,N}=K_{n,N}-k_n(X_N,H_N)
 =R_{n+2,N}-\{H_N,R_{n+1,N}\}+H_NR_{n,N}H_N.\tag{14}
\]
每个 \(k_n(X_N,H_N)\) 也精确守恒。极限的联合正交投影是 \(H\)-双模 map：在内积中把左/右 \(H\) 移到 commuting 测试函数即可证明。因此 \(k_n\) 是 \(K_n\) 的联合回归。\(L_n\) 的谱由有限线性组合和能量乘法获得，仍无对角质量，故 (13) 同样适用。

## 8. 先 ordinary Taylor 截断，再去除截断

已封存阶段给出
\[
 \|K_{n,N}\|_2\le K'(n+3)^4L^{n+2}\sqrt{n!},
\]
以及实际 charged 场的紧时间强替换
\[
 A_N(s)\sim(1-\mathcal P_H)e^{-iH_Ns/2}
 \sum_{n\ge0}\frac{(is)^n}{n!}K_{n,N}e^{-iH_Ns/2}.\tag{15}
\]
有限 \(N\) 层面使用 (15) 的统一尾及 pinching 的收缩性；不先假定 \(g_n(X_N,H_N)\) 是有限 \(N\) 正交投影。

先固定截断 \(n\le M\)，精确减去 \(k_n(X_N,H_N)\)，再用 (13)--(14) 去掉残量的 pinching。最后才令 \(M\to\infty\)。极限共同 Hilbert 空间中回归是正交投影，故
\[
 \|L_n\|_2\le\|K_n\|_2
\]
继承同一阶乘包络。这样完整 charged 时间场由绝对收敛的 Hilbert 空间级数构造，不只是一列 Taylor 系数。

若记 \(\nu^{L}_{nm}\) 为 (14) 的两能量谱，则
\[
 \lim_N\mathbb E\operatorname{tr}
 [A_N(s)^\dagger e^{iTH_N}A_N(u)e^{-iTH_N}]
\]
\[
 =\sum_{n,m\ge0}
 \overline{\frac{(is)^n}{n!}}\frac{(iu)^m}{m!}
 \int e^{\frac i2(s-u)(E+E')+iT(E'-E)}\,
             d\nu^{L}_{nm}(E,E').\tag{16}
\]
charged 的 (16) 没有 \(4C^2\) 前因子；该因子来自 neutral 分量的两次 \(2C\)。级数尾在紧 \(s,u\) 窗对所有实 \(T\) 一致，因为相位模为一、测度 variation 受场范数乘积控制。有限 \(N\) 的谱极限仍只对固定 \(T\) 成立，并可由高能矩紧性提升为紧 \(T\) 窗局部一致；没有全实轴 \(T\) 一致收敛结论。

截断场的谱在去截断时按 total variation 收敛，因此完整场谱仍绝对连续，且无静态对角质量。这里不宣称无限级数密度仍在 \(L^2\)；其 \(L^1\) 与绝对连续性已足够。

## 9. 只在极限 Hilbert 空间中作 Hermite 重求和

(7) 允许进一步简化，但这一步发生在 finite-\(N\) ordinary-tail/pinching 证明之后。由于 \(a_{n,d}\ge0\)，正交多项式展开给
\[
 \sum_d a_{n,d}\sqrt{C^dP_d}
 \le\sqrt{n+1}\left(\sum_d a_{n,d}^2C^dP_d\right)^{1/2}
 =\sqrt{n+1}\|H^n\|_{\rho_q}
 \le\sqrt{n+1}R_H^n.
\]
结合 (9)，
\[
 \sum_{n,d}\frac{S^n}{n!}a_{n,d}\|J_d\|_2
 \le K\sum_n\frac{(SR_H)^n}{n!}(n+1)^{5/2}<\infty.\tag{17}
\]
固定时间导数只增加 \(n\) 的固定幂，仍可和。令
\[
 c_d(s)=\frac{\int e^{isE}H_d(E)\,d\rho_q(E)}{C^dP_d}
       =\sum_n\frac{(is)^n}{n!}a_{n,d}.
\]
由 Tonelli/Fubini 的绝对 Hilbert 可和性，
\[
 \sum_n\frac{(is)^n}{n!}R_n=\sum_dc_d(s)J_d.\tag{18}
\]
两侧能量乘法在此极限空间中有界，范数至多 \(R_H\)。因此 charged 场的谱 amplitude 可写成
\[
 D_d(s;E,E')=-c_d''(s)+i(E+E')c_d'(s)+EE'c_d(s).\tag{19}
\]
共同双侧相位仍是 \(e^{-is(E+E')/2}\)。式 (19) 是 (16) 的合法重求和，不是把未证的 finite-\(N\) body/Jacobi bridge 当成输入。

低阶直接给
\[
 A(0)=J_2,
 \qquad A'(0)/i=J_3-\tfrac32\{H,J_2\}.
\]
若 \(\mathsf J_HY=\{H,Y\}/2\)，旧的
\(A_1=R+a\mathsf J_HY_-\)、\(a=-(1-q)^2/(1+q^2)\)
应等价地核对为
\[
 J_3=R+(3+a)\mathsf J_HY_-,\qquad
 3+a=\frac{2(1+q+q^2)}{1+q^2}.
\]
这一行给出与旧命名场比较的明确目标；其旧 \(R\) 的识别及相应 mixed 低阶校准需用共同图规则或旧精确 field 恒等式检查，不能仅凭上述符号改名当作新的独立核验。

## 10. 与 neutral 分量合并的全核记号

将既有 neutral 动态残量写为 \(N_j=V_j-v_j(H)\)，保持原 source、OU 和 descendants 定义。其场系数是 **正的** \(+2C c_j(s)\)，并带相同的双侧相位。取组合场族
\[
 \mathcal X_a\in\{J_d,N_j\},\qquad
 f_{J_d}(s;E,E')=D_d(s;E,E'),\quad
 f_{N_j}(s;E,E')=2C c_j(s).
\]
用 (11) 约定的共同正矩阵谱 \(\nu_{ab}\)，则完整 \(A+B\) 两点核为
\[
 \boxed{\sum_{a,b}\int
 e^{\frac i2(s-u)(E+E')+iT(E'-E)}
 \overline{f_a(s;E,E')}f_b(u;E,E')\,d\nu_{ab}(E,E').}\tag{20}
\]
这同时保留 AA、AB、BA、BB。charged 端由 (17)--(19) 绝对可和，neutral 端由已封存的 \((j+1)^2X_S^j/\sqrt{j!}\) 包络绝对可和；谱 Cauchy--Schwarz 给 mixed 双和的绝对收敛。能量依赖 amplitude 中的 \(E,E'\) 在极限谱支撑有界，不带来新尾问题。

(20) 的 \(\nu_{ab}\) 由逐图算法与 (12) 的 Hermite 展开确定，而非假设 AB/BA 为零。实际低阶比较必须恢复旧的非零
\(\nu_{A_0,B_2}=d_0(w+\eta)\mu_{q^2}\)，其中旧 \(d_0=4\lambda Cq^2\) 与这里的 Hermite 次数 \(d\) 是不同记号。本文证明其 general common-source 构造、静态剥离及时间级数存在；低阶命名场映射的独立枚举应另存新的科学输出，由该输出记录参数、版本及哈希。

本文新增闭合点是：联合 charge/energy 正交回归、共同原始耦合的强残量归约、全部固定 charged 残量的无对角谱、精确 pinching 传递，以及完整时间 charged/combined 核的收敛构造。它们是该 DSSYK 模型内结果；de Sitter 体质量/位置字典与有限 \(G\) 四点动力学仍是后续问题。
