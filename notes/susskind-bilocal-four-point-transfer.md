# 四点量的实际算符接口、pinching 与全时间边界

2026-11-07。接续 [完整 bilocal 强重建](susskind-bilocal-reconstruction-proof.md) 和 [joint-regression](susskind-memory-joint-regression-proof.md)。本稿核查新四点图究竟对应哪些实际算符，并给出可执行的非交换诊断。旧文件及旧科学输出不改写。

模型、态和极限保持：standard CAR、同一原始 Hermitian Gaussian 耦合、overlapping complex \(p\)-body Hamiltonian、\(\beta=\mu=0\) normalized annealed full trace；偶 \(N\)、\(p^2/N\to\lambda>0\)，固定 \(J_s>0\)、\(\hbar=1\)、string 时间。记 \(C=J_s^2/\lambda\)、\(q=e^{-4\lambda}\)，所有 \(L^r\) 均指含 disorder 与全 charge 权重的 Schatten 范数。

## 1. J2 与旧多 singlet 场的精确共同样本接口

令
\[
 \mathcal A=\sum_i[H,c_i]^\dagger[H,c_i],\quad
 DX=QX-\sum_i c_i^\dagger Xc_i,\quad C_Q=\mathbb E_JH^2,
\]
\[
 Z_2=H^2-C_Q,\qquad
 O=\mathcal A-\mathbb E_J\mathcal A,\qquad
 c=\kappa_2=2\lambda\frac{1-q}{1+q}.
\]
逐 sector 的 \(C_Q=C_mI\) 是确定性 charge 函数。由 \(DH=pH\) 展开两个 commutators，得到有限 \(N\) 精确式
\[
 \mathcal A=(2p-D)H^2.\tag{1}
\]
\(D\) 与 coupling Gaussian chaos 投影可交换，故
\[
 O=(2p-D)Z_2.\tag{2}
\]
新场和旧 [singlet-multimark](susskind-singlet-multimark.md) 的动态场分别为
\[
 J_{2,N}=O-cH_2^{(q,C)}(H)=O-c(H^2-CI),
 \qquad Y_N^{\rm old}=O-cZ_2.
\]
因此
\[
 \boxed{J_{2,N}-Y_N^{\rm old}=c(C-C_Q)I.}\tag{3}
\]
这不是两个边缘极限律的比较，而是同一原始耦合上的精确差式。\(C_m\) 一致有界，且由中央 charge 渐近与二项集中，在全迹下 \(C_Q\to C\)；于是 (3) 在每个固定 \(L^r\) 都趋零。

所以旧 \(Y^{\rm old}/\sqrt w\) 与 \(H/\sqrt C\) 的所有固定混合矩、固定高阶范数和实时多点结果可直接回接，其中
\[
 w=C^2\left(4\lambda+\frac{16\lambda^2q}{1+q}\right).
\]
尤其 \(\tau J_2^4\to(2+q^4)w^2\) 是旧结果的接口检查，不是本阶段的新四点发现。charge 回归的不同规范已由 (3) 明确处理。

## 2. Actual A0 的 pinching 不需要另一个 rank 桥

取原来的 charge-removal dilation，\(G=R^\dagger R\)、\(h=H_m\)、\(\alpha=1-p/m\)。直接由压缩块平方得
\[
 mG=\sum_i c_i^\dagger H_{m-1}^2c_i-m\alpha^2h^2
    =\mathcal A|_m-\frac{p^2}{m}h^2.\tag{4}
\]
零时 actual memory/charged coefficient 为
\[
 A_{0,N}=\mathcal F_N(0)=(1-\mathcal P_h)mG.
\]
\(\mathcal A-O=\mathbb E_J\mathcal A\)、\(C_Q\)、\(h^2\) 均与 \(h\) 对易。故在包括全部简并块的原 pinching 下
\[
 \boxed{A_{0,N}=(1-\mathcal P_h)O
                    =(1-\mathcal P_h)J_{2,N}.}\tag{5}
\]
在 \(m<p\) 上实际场为零；(5) 仍成立，因为剩余 \(J_2\) 只有 charge 标量。\(m=0\) 单独置零，不在 (4) 中使用 \(1/m\)。

前驱已证 \(\|\mathcal P_hJ_{2,N}\|_2\to0\)。旧多 singlet 定理与 (3) 已给每个固定 \(L^{2k}\) 一致界。pinching 对每个 Schatten 阶数收缩，所以对例如 \(r=4\)
\[
 \|A_{0,N}-J_{2,N}\|_4
 =\|\mathcal P_hJ_{2,N}\|_4
 \le\|\mathcal P_hJ_{2,N}\|_2^{1/3}
       \|J_{2,N}\|_8^{2/3}\longrightarrow0.\tag{6}
\]
选择任意更高偶矩，同样得每个固定有限 \(L^r\) 的强替换。这里的高矩不是从二点谱推断，而是使用已有多标记矩定理。

## 3. W1、rank 版本与 actual neutral 一次系数

局域 half-Clifford size 记为 \(\mathcal L\)，定义
\[
 W_{1,N}=(p-\mathcal L)H_N-\lambda H_N,
 \qquad N_{1,N}^{\rm loc}=4\lambda q W_{1,N}.\tag{7}
\]
前驱 size bridge 只直接给 actual rank 版本与局域版本的 \(L^2\) 强桥；不能把它未经论证提升为高点接口。但在这个 **一次 Gaussian** 系数上可单独闭合高矩。

设 \(X_m\) 是 Hermitian、\(U(N)\)-equivariant 的线性 Gaussian 场。\(\wedge^m\mathbb C^N\) 的不可约性与 Schur 引理给
\[
 \mathbb E X_m^2=v_mI,\qquad
 v_m=\mathbb E\operatorname{tr}_mX_m^2.
\]
写 \(X=\sum_\alpha g_\alpha A_\alpha\)，\(\sum A_\alpha^2=v_mI\)。Gaussian 分部积分展开 \(\mathbb E\operatorname{tr}X^{2k}\)，在 \(X\) 的本征基中对每个
\(|x_i|^l|x_j|^{2k-2-l}\) 使用加权算术几何平均，再对 \(\alpha,i,j\) 求和，得到
\[
 \mathbb E\operatorname{tr}_mX_m^{2k}
 \le(2k-1)v_m\mathbb E\operatorname{tr}_mX_m^{2k-2}
 \le(2k-1)!!v_m^k.\tag{8}
\]
这依赖 scalar covariance；不是一个对任意矩阵 Gaussian 多项式成立的 \(L^2\to S_{2k}\) hypercontractivity。

actual rank 场 \((p-\mathcal R)H-\lambda H\) 满足 (8)。旧 broad-central rank 二矩界使 \(v_m\) 一致有界；远区 \(\|(p-\mathcal R)H-\lambda H\|_{2,m}\le(p+\lambda)\sqrt{C_m}\)，其任意固定幂都被二项外尾压住。因此它有全迹的所有固定高矩界。局域 \(W_1\) 的所有固定高矩由本阶段[任意固定源数证明](susskind-bilocal-multisource-proof.md)第 3 节给出；该估计适用于任意固定 mark 数，而已保存有限检查只覆盖其声明的次数。于是旧 rank/size \(L^2\) 桥通过插值升为每个固定 \(L^r\) 强桥。

actual neutral 一次时间系数也可直接处理。保留旧定义
\[
 B_{1,N}=m\epsilon\alpha(1-\mathcal P_h)B_Z,
 \quad B_Z=\Phi(h)-\bar\kappa_m h,\quad \epsilon=1-\alpha^2,
\]
其中 \(\Phi\)、\(\bar\kappa_m\) 均为确定性 sector 参数/线性映射。未 pinching 的
\(T_{1,m}=m\epsilon\alpha B_Z\)
仍是上述 invariant Hermitian linear Gaussian 场。旧 two-component 稿已给 broad-central 的
\[
 \|B_Z\|_{2,m}\le C_m^{3/2}\frac{12p}{N}\sqrt{K_1},
 \qquad m\epsilon\le2p,
\]
所以 \(\|T_{1,m}\|_{2,m}\le24\lambda_N C_m^{3/2}\sqrt{K_1}\)；外区的粗界为 \(O(p)\)。用 (8) 和二项尾得到所有固定高矩的一致界，再用 \(1-\mathcal P_h\) 的 Schatten 范数 \(\le2\) 控制 actual \(B_1\)。

既有 full-charge 强接口是
\[
 \|B_{1,N}-8\lambda Cq W_{1,N}\|_2\to0.
\]
两边的高矩现已分别控制，故
\[
 \boxed{\|B_{1,N}-8\lambda Cq W_{1,N}\|_r\to0
        \quad\text{对每个固定有限 }r.}\tag{9}
\]
这也说明 \(N_1=4\lambda qW_1\)、\(B_1=2CN_1\) 的归一化。这里只对一次 Gaussian 场使用 (8)；更高 rank 后代一般是非线性 Gaussian 多项式，不能照搬该论证。

## 4. 怎样把四点和其他固定词传到实际场

对四个因子写 telescoping difference，并使用 tracial Hölder：
\[
 |\tau(X_1X_2X_3X_4-Y_1Y_2Y_3Y_4)|
 \le\sum_{j=1}^4\|X_j-Y_j\|_4
       \prod_{i\ne j}\max(\|X_i\|_4,\|Y_i\|_4).\tag{10}
\]
(6)、(9) 与一致高矩界使右侧趋零。有限个 \(e^{itH}\) 插入不会改变 Schatten 范数；固定 \(H\) 多项式插入则选用更高阶 Hölder 和既有能量高矩。因而本阶段 \(J_2,W_1\) 的四点确实是 actual \(A_0\) 与归一 actual neutral slope
\(B_1/(8\lambda Cq)\)
的四点。也可把该论证用于任意固定长度的词，只要每个所需局域矩已按任意固定 mark 数的定理构造。

实际低阶生成元的固定实时多点同样成立。对偶整数 \(r\)，实轴 Taylor 余项满足
\[
 \sup_{|t|\le S}\left\|e^{itH_N}-
       \sum_{n=0}^{K-1}\frac{(itH_N)^n}{n!}\right\|_r
 \le\frac{S^K}{K!}[(rK-1)!!M^{rK/2}]^{1/r}\longrightarrow0,
\]
且对 \(N\) 一致。固定插入数时，选定一个足够高但固定的 Hölder 阶数，与 markers 的固定高矩组合：先对有限 Taylor 词取 DS 极限，再解除 \(K\) 截断。因此这里的“实时四点”是 \(A_{0,N}(t)\)、\(B_{1,N}(t)\) 等固定生成元由 \(H_N\) 演化后的关联；它不等于下文尚未闭合的完整 relative-time bilocal 四点。

仅证明两个场各自的四矩存在不足以做 same-copy 替换；这里所用的是精确差式 (5) 及已建立的共同耦合强桥 (9)。rank/body-loss/charge 反项没有被当成可互换的有限 \(N\) map。

## 5. 三种场的逐 site 配对解释

\(W_1\) 插入只含一个 Hamiltonian slot。若 pair collision 只涉及其 incident chord 的一个 \(p\)-set，该 field 上的局域词为 odd \(c\) 或 \(c^\dagger\)。size 源只给标量缩放；它与 source-dependent singleton 分母逐项精确抵消。因此含 \(W\) 源的项只能来自该 Gaussian chord 自己的 creator/annihilator 两集合相碰。

own collision 最多接触该 chord 的两个端点 field，因而不可能产生带第三个 field source 的连通项。其单 \(W\) 系数为 \(\lambda\)，两端均为 \(W\) 时 mixed 系数为 \(\lambda\)；若另一端是 charged two-slot body field，mixed 系数为零。减去 \(\lambda H\) 后单源消失，所以 \(W\) 必须与另一个 \(W\) 沿同一 Gaussian chord 成对。

剩余 two-slot \(J_2\) 的相邻端点交换抵消正是旧 singlet-multimark 证明；\(W\) 弦作为 spectator 不改变交换轨道的比值。于是三个归一化 species
\[
 H/\sqrt C,\qquad W_1/\sqrt{C\lambda},\qquad J_2/\sqrt w
\]
按同种字母配对，交叉权由有效弦长度 \((1,1,2)\) 给出：
\[
 Q_{HH}=Q_{HW}=Q_{WW}=q,\quad
 Q_{HJ}=Q_{WJ}=q^2,\quad Q_{JJ}=q^4.\tag{11}
\]
这段解释限定于上述 one-slot \(W\) 与 two-slot \(J_2\) 家族；不把一般更高后代的多源式预设成相同的三 species 模型。

## 6. 可执行的非交换诊断

记 \(v=C\lambda\)。由 (11)，
\[
 \tau J_2^2W_1^2=wv,\qquad
 \tau J_2W_1J_2W_1=q^2wv,\qquad
 \tau W_1^4=(2+q)v^2.\tag{12}
\]
两种 ordered products \(J_2W_1\)、\(W_1J_2\) 的归一 Gram 矩阵为
\[
 \Gamma=\begin{pmatrix}1&q^2\\q^2&1\end{pmatrix},
 \qquad\det\Gamma=1-q^4>0.\tag{13}
\]
对任何有限矩阵自伴 \(A,B\)，循环迹给精确恒等式
\[
 \tau([A,B]^\dagger[A,B])
 =2\{\tau(A^2B^2)-\tau(ABAB)\}.\tag{14}
\]
因此这里
\[
 \boxed{\|[J_2,W_1]\|_2^2=2wC\lambda(1-q^2)>0.}\tag{15}
\]
同时 \(\tau[J_2,W_1]=0\)，且二点内积 \(\tau J_2W_1=0\)。所以这个非零 commutator 不可能是标量：若为 \(cI\)，其零迹强制 \(c=0\)，与 (15) 矛盾。这是算符值的非交换诊断；它不构成一个真空 CCR 常数，也不从二点正交推断可交换或独立。

由 (10)，同一结论适用于 actual \(A_0\) 与 \(B_1/(8\lambda Cq)\)。未归一的 \([A_0,B_1]\) 范数平方再乘 \((8\lambda Cq)^2\)。

新 [exact driver](../code/bilocal_four_point_transfer_check.py) 已首次运行并保存 [输出](../results/susskind-bilocal-four-point/transfer-checks.json)：三个完整 Fock fixtures \((N,p)=(2,1),(3,2),(4,2)\)，全部 Fraction 算术，33 个矩阵恒等式、3030 个逐 entry 比较。它独立核 (1)--(5) 的不经数值本征向量的代数接口、Hermiticity/charge preservation，以及 (14) 和 commutator 零迹。最后一个 fixture 的 commutator 范数平方为 \(73567/144\)，前两个为零 controls。

这些 fixture 用确定的实 Hermitian 整数耦合；\(w_N=1\)、formal \(q=1/2\) 和 \(\lambda=p^2/N\) 仅用于有限代数反项，不强制 \(q=e^{-4\lambda}\)。它们不是 Gaussian 抽样、DS 拟合或物理参数下的四点极限数值。输出保存完整参数、耦合、Python 3.14.5 和 driver SHA256
`33bbe8b8e5a0d5d29758b7d4946f3d98cf9c011a8a9179abcf1e9c21a51fdbdf`。
本阶段 root 的独立 source 图输出负责检验 (12)--(13) 的模型权重；不要把这份有限代数检查当成该渐近定理的数值证明。

## 7. 强返回本身可以升到任意固定 Schatten 阶数

旧 bilocal 强返回只在正文显式写了 \(L^2\) 速率，但其输入力矩定理对任意固定阶都成立。设 \(r\ge2\) 固定，在 good 区 \(m\ge N/4\) 保持
\[
 \|\delta M\|_2=O(p^{-2}),\qquad
 \|\delta M\|_{4r}=O(p^{-1}),
 \quad\|F-U_h\|_{2r}+\|F_0-U_h\|_{2r}=O(p^{-1}).
\]
插值指数由
\(1/(2r)=\theta/2+(1-\theta)/(4r)\)
给出 \(\theta=1/(2r-1)\)，故
\[
 \|\delta M\|_{2r}=O(p^{-1-1/(2r-1)}).
\]
在 exact resolvent 差中，Hölder 加 \(m\le N\asymp p^2\) 给
\[
 \boxed{\sup_{|t|\le S}
 \|\mathcal B_N(t)-\widetilde{\mathcal B}_N(t)\|_r
 =O_S(p^{-1/(2r-1)})+O_S(p e^{-N/(8r)}).}\tag{16}
\]
bad active 区直接用 \(\|m\delta M\|_{r,m}\le K_rp\) 和二项尾；inactive 场仍严格零。这里 \(\widetilde{\mathcal B}\) 是前驱对 **实际 memory 场** 的双重线性滤波。

特别 \(r=4\) 使用旧力的第 16 阶矩，得到 \(L^4\) 返回误差 \(O(p^{-1/7})\)。这是新的直接 Hölder/插值推论；不必对原 \(L^2\) 返回误差使用一个不足以消除增长的粗高范数插值。

## 8. 为什么仍未完成 full-time actual bilocal 四点

(16) 本身并不给 full-time bilocal 四点收敛。四点 telescoping (10) 还需要四个实际/滤波场的 \(L^4\) 一致界。当前旧力矩只给
\[
 \sup_{|t|\le S}\|\mathcal B_N(t)\|_r+
 \sup_{|t|\le S}\|\widetilde{\mathcal B}_N(t)\|_r=O_S(p)
\]
的粗界，并另有统一 \(L^2\)。它们不能推出统一四矩。标量反例 \(X_p=p\mathbf1_{A_p}\)、\(\mathbb P(A_p)=p^{-2}\) 已满足 \(\|X_p\|_2=1\)、每个固定 \(\|X_p\|_r\le p\)，但 \(\mathbb E X_p^4=p^2\) 发散。因此这不是可以省掉的技术句子。

一条足够且可继续执行的下一条件是，对有限 field/ordinary-degree 截断 \(\widetilde{\mathcal B}_N^{(K)}\) 证明
\[
 \lim_{K\to\infty}\limsup_N\sup_{|t|\le S}
 \|\widetilde{\mathcal B}_N(t)-
          \widetilde{\mathcal B}_N^{(K)}(t)\|_4=0.\tag{17}
\]
各固定截断的高矩由有限源图控制，(17) 才允许把共同四点从有限截断传到全时间级数。新“任意固定 mark/degree”矩定理不自动提供对 field degree 的可和常数；旧 \(L^2\) 阶乘尾也不能直接改写成 \(L^4\) 尾。若先建立滤波场统一 \(L^4\) 而尚无 (17)，(16) 可先传递 actual 与 filtered 四点的差，但仍需独立识别 filtered 四点的极限。

当前可无条件写出的 actual bilocal 小时间四点是**先小时间、后 DS** 的曲率结果。例如有限 \(N\) 有
\(\mathcal B_N(t)=-t^2A_{0,N}/2+O_N(t^3)\)，所以
\[
 \lim_{N\to\infty}\lim_{t\to0}
 \frac{16}{t^8}\mathbb E\operatorname{tr}
 [\mathcal B_N(t)^\dagger\mathcal B_N(t)]^2
 =(2+q^4)w^2.\tag{18}
\]
此处不能交换两个极限，也没有把只有二点定义的完整时间场自行当成已构造的多点代数。当前新四点诊断 (12)--(15) 完整适用于实际曲率与实际 neutral slope；full-time actual bilocal 四点的可和尾仍是明确的下一接口。
