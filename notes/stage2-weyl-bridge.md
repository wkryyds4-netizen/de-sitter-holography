# 压缩自由场的自伴性、附属关系及 Weyl 混合导数

日期：2026-10-29。对象是第二阶段的自由场—全线钟表示；不加入有限 G 相互作用。本稿新增明确证明，保留 [原物理审查](stage2-physics-audit.md) 中此前尚未证明的历史状态。空间局域化需另接 [静态管近似](stage2-localization-bridge.md)；这里不能仅凭正世界线谱替代实际局域场论。

结论：在下列可核查的充分假设下，共同钟压缩场在包含任意 L² 钟波函数的有限粒子域上本质自伴，闭包附属于指定角代数。至多四次插入的有界压缩 Weyl 中心差商，在任意固定粒子数图范数中联合收敛到压缩场乘积。单粒子范数收敛的空间/时间逼近可与这些差商交换极限。这是对既有自由场构造的严格性补充，不声称新的全息对偶或数学定理首创。

## 1 对象、归一化及所需假设

设

\[
\mathcal H=\Gamma_s(\mathfrak k)\otimes\mathcal K_c,
\quad \mathcal K_c=L^2(\mathbb R),
\quad N=d\Gamma(1)\otimes1,
\quad E_n=1_{[0,n]}(N).
\tag{1}
\]

\(\mathfrak k\) 可以是热准自由 GNS 的加倍单粒子空间。此时 Fock 向量 \(\Omega\) 实现 BD/KMS 参考态，并非静态物理零温真空；\(N\) 是辅助 Fock 粒子数，不是已识别的局域可测粒子数。定义实 Segal 场

\[
A(v)=a(v)+a^\dagger(v),\qquad
\|A(v)\Omega\|^2=\|v\|^2.
\tag{2}
\]

“实”指自伴场及实线性标记；\(v\) 自身不必在某个坐标基中处处为实。不同的 \(1/\sqrt2\) 约定须先换算到 (2)。

所需条件为：

1. \(U\) 为 \(\mathcal H\) 上的酉算符，\(\Pi\) 为正交投影，二者均与 \(N\) **强对易**，即对每个 \(n\) 与 \(E_n\) 对易。
2. \(\widehat{\mathcal M}\) 是全线钟表示中的指定 von Neumann 代数，\(\Pi\in\widehat{\mathcal M}\)。在角空间 \(\mathcal H_+=\Pi\mathcal H\) 上令 \(\mathcal M=\Pi\widehat{\mathcal M}\Pi\)，其单位为 \(I_+=\Pi|_{\mathcal H_+}\)。
3. 对每个使用的 \(v\) 和**所有实数** \(\lambda\)，有

\[
W_v(\lambda)=e^{i\lambda A(v)}\in\mathcal A_0,
\qquad U^\dagger W_v(\lambda)U\in\widehat{\mathcal M}.
\tag{3}
\]

本项目与静态管证明拼接时，在实际 BD 自由场 Fock 表示上取 \(\mathcal A_0=\mathcal A(D)\)，\(H\) 是该场的真实静态等距生成元，\(N\) 是同一表示中的辅助粒子数。静态管笔记的式 (28) 显式给出共形实现中的单粒子酉群，其第二量子化保持 \(N\)。静态补丁的标准 BD/KMS 模流识别作为物理 crossed-product 输入；这里不以只有世界线谱的循环 Gaussian 子模型替代整个 \(\mathcal A_0\)。也不任意限制到尚未证明受 \(N\) 保持的局域循环子空间，再沿用下述有限粒子域。

\(U\) 本身不必属于 \(\widehat{\mathcal M}\)。在物理实现中 \(H=d\Gamma(h)\)、\(H\Omega=0\)，\(h\) 可有正负谱；先在全线钟上定义 \(U=e^{i\tau H}\)、\(\Pi=1_{[0,\infty)}(\pi)\)。由于 \(H\) 与 \(N\) 强对易且 \(\Pi\) 只作用于钟，条件 1 成立。\(\tau,\pi\) 在全线是正则共轭对；这里没有把 \(\tau\) 当作正半线能量空间上的自伴正则时间算符。

条件 2、3 的 crossed-product 入口来自 [Han–Hu, 2610.14175v1，§II 式 (8)、(12)—(14)](https://arxiv.org/html/2610.14175v1)。其适用对象首先是有界 Weyl 算符。局域测试函数逼近若满足 \(v_\delta\to v\)，且每条 \(W_{v_\delta}(\lambda)\) 属于同一个 \(\mathcal A_0\)，则 Weyl 强连续性和 von Neumann 强闭性推出 (3)。这一步的实际 dS 输入由上述静态管笔记负责。

仅有某一个 \(e^{iA}\in\mathcal A_0\) 不够：取二维 \(A=\operatorname{diag}(0,2\pi)\)，则 \(e^{iA}=I\) 属于标量代数，但 \(A\) 不附属于标量代数。因此必须写出整条一参数族，不能省掉“对所有 \(\lambda\)”的条件。

## 2 任意平方可积钟的共同域

在角空间定义

\[
\mathcal D_{\rm fin}^+
=\bigcup_{n\ge0}\Pi E_n\mathcal H,
\quad T(v)=U^\dagger A(v)U,
\quad B(v)=\left.\Pi T(v)\Pi\right|_{\mathcal D_{\rm fin}^+}.
\tag{4}
\]

这里用的是完整的有限粒子扇区，不只有限个张量积的代数线性包。每个扇区允许任意 L² 钟波函数及场—钟纠缠；\(E_n\Pi=\Pi E_n\) 保证该域在 \(\mathcal H_+\) 中稠密。没有要求钟的能量矩或时间导数有限。

由创造、湮灭算符标准估计，\(\psi\in E_n\mathcal H\) 时

\[
\|A(v)\psi\|\le2\|v\|\sqrt{n+1}\,\|\psi\|,
\qquad A(v)E_n\mathcal H\subset E_{n+1}\mathcal H.
\tag{5}
\]

因为 \(U\) 保持粒子数、\(\Pi\) 是保持粒子数的压缩，(5) 原样适用于 \(T(v)\) 及角空间上的 \(B(v)\)。因此 \(B(v)\) 是稠密定义的对称算符，并将 \(\mathcal D_{\rm fin}^+\) 映回自身。对于 \(\psi\in\Pi E_n\mathcal H\)，逐次使用 (5) 得

\[
\|B(v)^k\psi\|\le
(2\|v\|)^k\sqrt{\frac{(n+k)!}{n!}}\,\|\psi\|.
\tag{6}
\]

对任意有限 \(t\ge0\)，级数 \(\sum_k t^k\|B(v)^k\psi\|/k!\) 收敛：其控制项之比为 \(2t\|v\|\sqrt{n+k+1}/(k+1)\to0\)。每个有限粒子向量都是整解析向量。

Nelson 解析向量定理因而给出

\[
\boxed{B(v)\text{ 在 }\mathcal D_{\rm fin}^+\text{ 上本质自伴。}}
\tag{7}
\]

下文以 \(\overline B(v)\) 表示唯一自伴闭包。所用定理的精确可查版本为 [Neeb, Ann. Inst. Fourier 61 (2011), p.1857, Theorem 4.8](https://www.numdam.org/item/10.5802/aif.2660.pdf)：对称算符若有稠密的解析向量集合，即本质自伴。原始来源是 [Nelson, “Analytic Vectors”, Ann. Math. 70 (1959), 572–615](https://doi.org/10.2307/1970331)。此处只应用该标准定理，关键模型检查是 (5)、(6) 和域的完整性。

若 \(U\) 或 \(\Pi\) 不保持粒子数，则本证明失去依据；一般自伴算符的任意压缩并不因“压缩”二字自动本质自伴。

## 3 从有界代数元素证明附属关系

在 \(\mathcal H_+\) 上定义

\[
C_v(\lambda)=\Pi U^\dagger W_v(\lambda)U\Pi\in\mathcal M,
\quad C_v(0)=I_+,
\quad C_v(\lambda)^*=C_v(-\lambda).
\tag{8}
\]

这些是压缩后的缩并算符，一般不是酉算符。对 \(\epsilon\ne0\) 定义有界、自伴的中心差商

\[
S_\epsilon(v)=\frac{C_v(\epsilon)-C_v(-\epsilon)}{2i\epsilon}
=\Pi U^\dagger\frac{\sin(\epsilon A(v))}{\epsilon}U\Pi
\in\mathcal M.
\tag{9}
\]

谱定理与 \(|\sin(\epsilon x)/\epsilon|\le|x|\) 保证，在 \(\mathcal D_{\rm fin}^+\) 上

\[
S_\epsilon(v)\psi\longrightarrow B(v)\psi.
\tag{10}
\]

下面直接证明所需强 resolvent 极限，不借未经核实的定理编号。记 \(\overline B=\overline B(v)\)。对 \(\psi\in\mathcal D_{\rm fin}^+\)、\(\eta=(\overline B-i)\psi\)，有

\[
(S_\epsilon-iI_+)^{-1}\eta-\psi
=(S_\epsilon-iI_+)^{-1}(\overline B-S_\epsilon)\psi\to0,
\tag{11}
\]

因为自伴性给出 resolvent 范数 \(\le1\)。\(\mathcal D_{\rm fin}^+\) 是 \(\overline B\) 的核，所以 \((\overline B-i)\mathcal D_{\rm fin}^+\) 在 \(\mathcal H_+\) 中稠密。一致范数界将 (11) 扩展到所有 \(\eta\)，得到

\[
(S_\epsilon-iI_+)^{-1}\xrightarrow{\rm strong}
(\overline B-iI_+)^{-1}.
\tag{12}
\]

每个左端由有界 Borel 泛函演算属于 \(\mathcal M\)，故其强极限也属于 \(\mathcal M\)。若 \(Q\) 是 \(\mathcal M'\) 中任意酉算符，它与该 resolvent 对易，因而保持其值域 \(D(\overline B)\)，并在该域上与 \(\overline B\) 对易。这正是附属关系：

\[
\boxed{\overline B(v)\ \eta\ \mathcal M.}
\tag{13}
\]

因此 \(\overline B(v)\) 的谱投影及 \(e^{it\overline B(v)}\) 属于 \(\mathcal M\)。附属算符与有界 Borel 泛函演算的对应也可参见 [Budde–Landsman, arXiv:1508.06772v1, Theorems 3–4](https://arxiv.org/html/1508.06772)。本节的 resolvent 证明已经给出当前构造所需的域控制。

## 4 粒子数图范数中的多重差商

定义图范数尺度

\[
\|\xi\|_s=\|(N+1)^{s/2}\xi\|,\qquad s\ge0.
\tag{14}
\]

固定有限插入数 \(r\le4\)，\(\|v_j\|\le M\)，为方便取 \(M\ge1\)。以下乘积保持显示的顺序，右端算符先作用。对 \(\psi\in\Pi E_n\mathcal H\)，\(k_j\ge0\)、\(K=\sum_j k_j\)，记

\[
X_{\mathbf k}\psi
=\Pi T(v_1)^{k_1}\Pi T(v_2)^{k_2}\Pi\cdots
\Pi T(v_r)^{k_r}\Pi\psi.
\]

在任意两组幂之间插入 \(\Pi\) 不会增加范数或粒子数，故

\[
\|X_{\mathbf k}\psi\|_s
\le (n+K+1)^{s/2}(2M)^K
\sqrt{\frac{(n+K)!}{n!}}\,\|\psi\|.
\tag{15}
\]

这统一控制所有有序幂串，未假定场互相对易。

对每个非零 \(\epsilon_j\)，\(S_{\epsilon_j}(v_j)\) 的正弦 Taylor 级数在有限粒子向量上绝对收敛。其有序乘积可写成

\[
\prod_{j=1}^r S_{\epsilon_j}(v_j)\psi
=\sum_{m_1,\ldots,m_r\ge0}
\left[\prod_{j=1}^r
\frac{(-1)^{m_j}\epsilon_j^{2m_j}}{(2m_j+1)!}\right]
X_{(2m_1+1,\ldots,2m_r+1)}\psi.
\tag{16}
\]

合法性不仅来自形式微分。对 \(|\epsilon_j|\le1\)，利用

\[
\sum_{k_1+\cdots+k_r=K}\frac1{k_1!\cdots k_r!}
=\frac{r^K}{K!}
\]

及 (15)，多级数在 \(\|\cdot\|_s\) 中被以下收敛级数控制：

\[
D_{n,r,s,M}=\sum_{K\ge0}
(n+K+1)^{s/2}(2Mr)^K
\frac{\sqrt{(n+K)!/n!}}{K!}<\infty.
\tag{17}
\]

任意固定有界 \(\epsilon_j\) 区间也只需在 (17) 中换掉常数。范数收敛时外层 \(S_\epsilon\) 是有界算符，可以先作用于内层的有限部分和；每个部分和仍是有限粒子向量。再用绝对多级数收敛及闭算符 \((N+1)^{s/2}\) 完成逐层展开与重排，便得到 (16) 在每个固定图范数中的等式。无需假定 Weyl 算符保持有限粒子域；它通常不保持。

(16) 的 \(m_j=0\) 项正是 \(B(v_1)\cdots B(v_r)\psi\)。所有其余项总共至少含两个 \(\epsilon\) 的幂，因此

\[
\boxed{\left\|
\prod_{j=1}^r S_{\epsilon_j}(v_j)\psi
-\prod_{j=1}^r B(v_j)\psi
\right\|_s
\le D_{n,r,s,M}\,\max_j|\epsilon_j|^2\,\|\psi\|,
\quad |\epsilon_j|\le1.}
\tag{18}
\]

这是各 \(\epsilon_j\) 独立趋零的联合极限。它比“各算符有强 resolvent 极限”更强；后一句单独不足以推导无界乘积的极限。(17) 的常数故意保守，不用于数值最优化。论证事实上适用于任意固定有限 \(r\)，本阶段只使用二点和四点。

## 5 空间逼近、时间涂抹逼近与交换极限

设真实局域测试函数给出的 \(v_{j,\delta}\) 满足

\[
v_{j,\delta}\to v_j\text{ in }\mathfrak k,
\qquad \sup_{j,\delta}\|v_{j,\delta}\|\le M.
\tag{19}
\]

这里 \(U,\Pi,N\) 固定，不随 \(\delta\) 改变。由场的实线性和 (5)，

\[
\|(B(v')-B(v))\psi\|
\le2\|v'-v\|\sqrt{n+1}\|\psi\|.
\tag{20}
\]

将有限个场的乘积作逐项 telescoping，在 \(\psi\in\Pi E_n\mathcal H\) 上得到

\[
\left\|\prod_{j=1}^rB(v_{j,\delta})\psi
-\prod_{j=1}^rB(v_j)\psi\right\|_s
\le L_{n,r,s,M}
\sum_j\|v_{j,\delta}-v_j\|\,\|\psi\|,
\tag{21}
\]

其中可取

\[
L_{n,r,s,M}=2(2M)^{r-1}
(n+r+1)^{s/2}\sqrt{(n+r)!/n!}.
\]

将 (18) 对有界族 \(v_{j,\delta}\) 一致应用，再接 (21)，有

\[
\boxed{\begin{aligned}
&\left\|\prod_{j=1}^r S_{\epsilon_j}(v_{j,\delta})\psi
-\prod_{j=1}^rB(v_j)\psi\right\|_s\\
&\quad\le\left[D_{n,r,s,M}\max_j|\epsilon_j|^2
+L_{n,r,s,M}\sum_j\|v_{j,\delta}-v_j\|\right]\|\psi\|.
\end{aligned}}
\tag{22}
\]

所以这些逼近可以联合取极限，也可以按任意次序取迭代极限。为核实先固定 \(\epsilon\) 再取 \(\delta\to0\) 的图范数极限，对 (16) 每个有限幂串作与 (21) 相同的逐因子差分；各项连续，(15)—(17) 给出与 \(\delta\) 无关的可和控制，故级数的支配收敛给出每个固定图范数中的极限。有界 Weyl 的强连续性也可单独核对其 Hilbert 范数极限。对于实时间平移 \(v_j(t_j)=e^{it_jh}v_j\)，同一个酉群使 \(\|e^{it_jh}(v_{j,\delta}-v_j)\|=\|v_{j,\delta}-v_j\|\)；因此当前估计的常数还可对全部实 \(t_j\) 一致。改变动力学或进行复时间延拓不能直接沿用这一结论。

这也给出真正的场算符图范数控制。由全局创造—湮灭界以及 \(U,\Pi\) 与 \(N\) 强对易，

\[
\|B(u)\xi\|\le2\|u\|\|(N+1)^{1/2}\xi\|
\quad (\xi\in\mathcal D_{\rm fin}^+).
\]

对 \(\xi\in D((N+1)^{1/2})\cap\mathcal H_+\) 作粒子数截断，再利用 \(\overline B(u)\) 的闭性，可扩展为

\[
D((N+1)^{1/2})\cap\mathcal H_+\subset D(\overline B(u)),
\qquad
\|\overline B(u)\xi\|\le2\|u\|\|\xi\|_1.
\tag{23}
\]

于是 (22) 的 \(s=1\) 收敛也蕴含任一固定 \(\overline B(u)\) 的通常图范数收敛。我们不声称 \(D(\overline B(u))=D((N+1)^{1/2})\)。

(19) 是必须证明的条件；未涂抹的时空点场一般不满足它。对 Gaussian 的时间截断、cos² 的光滑化、空间管半径趋零，各步骤都要通过实际单粒子范数。仅有等时方差有界而没有交叉范数收敛，不足以应用本节。

## 6 真混合导数与有界角代数的矩识别

令 \(\psi=\Omega\otimes\zeta\)，\(\zeta\) 为正能量半线钟空间中的任意归一化 L² 波函数，零延拓到全线；则 \(\Pi\psi=\psi\)、\(N\psi=0\)。定义普通有界代数向量态

\[
\omega_\zeta(X)=\langle\psi,X\psi\rangle,
\quad
F(\lambda_1,\ldots,\lambda_r)
=\omega_\zeta\big(C_{v_1}(\lambda_1)\cdots C_{v_r}(\lambda_r)\big).
\tag{24}
\]

每组实参数下，括号内的有界算符确实属于 \(\mathcal M\)。为了说明不仅中心差分收敛，而且真正混合导数存在，对每个指数使用完整 Taylor 级数。相应向量级数为

\[
\prod_j C_{v_j}(\lambda_j)\psi
=\sum_{k_1,\ldots,k_r\ge0}
\frac{i^{K}\lambda_1^{k_1}\cdots\lambda_r^{k_r}}
{k_1!\cdots k_r!}\,X_{\mathbf k}\psi.
\tag{25}
\]

(15)、(17) 证明此级数在每个有限参数多圆盘、每个固定粒子数图范数中绝对且一致收敛；各阶导数的级数只增加 \(K\) 的有限幂，依旧收敛。因此它给出参数的整向量函数，实轴上等于原有界算符乘积作用在 \(\psi\) 上。这里的整延拓是向量函数，不是宣称所有复参数下算符有界。

取每个参数的一阶导数，唯一保留的系数为 \(k_1=\cdots=k_r=1\)，从而

\[
\boxed{
i^{-r}\partial_{\lambda_1}\cdots\partial_{\lambda_r}
F(0,\ldots,0)
=\langle\psi,B(v_1)\cdots B(v_r)\psi\rangle.}
\tag{26}
\]

定义第 \(j\) 个变量的中心差分算符

\[
(\Delta_{j,\epsilon}F)(\boldsymbol\lambda)
=\frac{F(\boldsymbol\lambda+\epsilon e_j)
-F(\boldsymbol\lambda-\epsilon e_j)}{2\epsilon}.
\]

于是严格恒等式为

\[
i^{-r}(\Delta_{1,\epsilon_1}\cdots
\Delta_{r,\epsilon_r}F)(\mathbf0)
=\omega_\zeta\left(\prod_jS_{\epsilon_j}(v_j)\right).
\tag{27}
\]
(18)、(22) 给出其收敛及逼近交换的明确界。对单位钟态该界与 \(\zeta\) 无关，所以不需要在钟积分中假设尚未证明的高阶能量矩。

对于本项目的物理 \(U\)，用 \(U\Pi U^\dagger=1_{[0,\infty)}(\pi-H)\)、\(U\psi=\psi\)，上述共同域进一步保证

\[
\langle\psi,B_1\cdots B_r\psi\rangle
=\int_0^\infty |\zeta(p)|^2
\langle\Omega,A_1P_pA_2P_p\cdots P_pA_r\Omega\rangle\,dp,
\quad P_p=1_{(-\infty,p]}(H).
\tag{28}
\]

右端有与 \(p\) 无关的有限粒子界 \((2M)^r\sqrt{r!}\)，足以对归一化 \(|\zeta|^2\) 积分。它是 [四点谱测度笔记](stage2-wick-measures.md) 的共同钟产品入口。因而在实际局域化 (3)、(19) 已建立后，该笔记的两点/四点矩既是压缩自伴场在共同域上的矩，也是原有界角代数中 Weyl 产品的混合导数。完整中间投影链保留，没有把穿衣后四点重新 Wick 分解。

## 7 必须保留的区别

\(C_v(\lambda)\) 与 \(e^{i\lambda\overline B(v)}\) 是一般不同的有界角代数元素。最简单的结构反例在真正 Fock 空间已经存在：取 \(U=I\)、\(\Pi=|\Omega\rangle\langle\Omega|\)（它与 \(N\) 对易），则角空间上一阶压缩场 \(B(v)=0\)，但

\[
\Pi e^{i\lambda A(v)}\Pi
=e^{-\lambda^2\|v\|^2/2}\Pi
\ne\Pi=e^{i\lambda\overline B(v)}
\quad(\lambda\|v\|\ne0).
\tag{29}
\]

此 \(\Pi\) 是说明运算结构的真空投影，不冒充物理正钟投影。二阶导数中 \(\Pi T(v)^2\Pi\) 与 \(\Pi T(v)\Pi T(v)\Pi\) 已经不同。式 (26) 每个插入使用不同变量并各微分一次，所以在插入之间保留全部 \(\Pi\)，得到的是后者类型的共同钟乘积。

本稿证明了当前自由准自由表示、任意归一化正钟 L² 向量、固定有限阶产品的算符接口。它没有证明任意正常态的所有无界矩存在，没有将压缩 Weyl 变成酉群，没有证明未经核查的无限时长窗口属于某个固定有限时空管，也没有加入观察者有限质量反冲或有限 G 引力动力学。代数附属、态的 KMS 性、读时分辨率与引力约束仍是不同问题。

## 8 可复用结论与依据

- 本质自伴性：完整有限粒子域、(5)—(7)，允许任意 L² 钟。
- 附属关系：整条局域 Weyl 族是输入；(9)—(13) 通过有界代数内的 resolvent 强极限得到。
- 二点/四点字典：(15)—(27) 控制有序多重差商、真混合导数和单粒子逼近交换，均保留所有中间投影。
- 与实际 dS 场拼接：由静态管笔记核验 (3)、(19)，由原始场谱及 Wick 测度笔记确定实际关联函数；不能由抽象谱直接推断整个静态补丁代数。

来源定位：Han–Hu 2610.14175v1 为既有 crossed-product 与共享钟构造；Nelson 1959／Neeb 2011 为标准解析向量判据；Budde–Landsman 1508.06772v1 为附属与泛函演算的标准对应。本稿的详细范数估计和交换极限是对本项目当前模型的直接推导。没有为这些标准工具或推论主张首创性。有限振子数值脚本只能核验实现和公式结构，不替代上述无限维证明及 dS 局域化。