# 记忆二次通道的 rank 滤波与输出尾界

2026-11-06，本地。[六页PDF](../output/pdf/de-sitter-susskind-memory-rank-blocks.pdf)、[LaTeX](../manuscript/susskind-memory-rank-blocks/main.tex)、[来源](../literature/susskind-memory-rank-blocks-sources.md)。不 commit/push/ZIP。

沿上一轮 W2Ecal 精确算符继续：先读现有不可约张量乘法，再分解已有 N5p2m3 的残量。这次没有增加一串新的小系统来猜 DS；新大参数结果来自一般 N 的协方差复合与精确 Casimir 矩。

## 有限残量的具体成因

标准 CAR、含指标重叠 complex Gaussian p-body H，固定 charge m，annealed normalized trace，ℏ1。有限 witness 使用 w1。对 output 中性算符，C2(X)=Σij[Eij,[Eji,X]]，Eij=c_i†c_j，rank r 特征值 c_r=2r(N+1−r)，维数d_r=binom(N,r)²−binom(N,r−1)²。

用精确多项式投影 P_r=∏_(s≠r)(C2−c_s)/(c_r−c_s)。在 p-particle coupling 空间的两个 Gaussian 输入指标上也分别投影。二次 tensor 的 (a,b;r)、a≤b 分块彼此正交，但保留块内重数，不能默认一个 rank 三元组只有一个 intertwiner。

已有 N5p2m3,w1 的四个方向 B=H²−C、O=A−a0、J_HF、W2Ecal（C18，F=W1T3+86H），18 个 block 中7个零、11个非零，全部非零4×4 Gram 的秩恰为1。各块求和严格恢复上一轮完整Gram；去r0后严格恢复逐样本traceless Gram。

令 H_a 为原同一耦合的 rank-a 部分；Z_aa;r=P_r(H_a²)，Z_ab;r=P_r{H_a,H_b}（a<b）。存在表中精确标量使

\[
(H^2,A,J_HF,Q_{Ecal})=\sum_{a\le b,r}(1,o_{ab;r},f_{ab;r},g_{ab;r})Z_{ab;r}.
\]

完整11行表见稿第2页及 rank-reconstruction.json。比如 output r1 的 g 在 input(0,1),(1,1),(1,2),(2,2) 分别334/3、88、595/4、−20，不能压成单一 output-rank 函数。g 也不是所有 r>0 块上统一的 1/o/f 线性组合，所以旧三方向全局投影后仍有 norm²=248945355/14903。

原 finite witness 的动态残量仍真实，但本轮表明在这个系统里其来源是通道权重不同；并未发现块内额外乘法方向。不能把这解释成一般 U(N) tensor product 无重数，也不能外推任意 N,p,m 均 collinear。完整Gaussian coefficient Gram的零Schur范数证明 tensor恒等式；另外两个整数耦合样本逐矩阵恢复8式作独立核对。

## 正的乘积输出 rank 分布

考虑同一 charge m 空间上独立、各向同性、Hermitian rank-a/rank-b 算符 U,V，EU²=EV²=I。定义

π_ab(r)=E||P_r(UV)||²_(tr_m)，d=binom(N,m)，R=min(m,N−m)。

π≥0、Σπ=1、r>a+b 时π0。纯rank协方差 Phi_a(X)=EUXU 的 rank-s 特征值为前驱已有

\[
\chi_a(s)={}_4F_3\left(\begin{matrix}-a,a-N-1,-s,s-N-1\\-m,m-N,-N-1\end{matrix};1\right).
\]

由 E(UV)†X(UV)=Phi_b Phi_a(X)，有 χ_a(s)χ_b(s)=Σ_r π_ab(r)χ_r(s)。协方差 Choi 矩阵的 rank projector 以及 reshuffling 的 HS 等距给

Σ_s d_s χ_a(s)χ_b(s)=d² δ_ab/d_a，故

\[
\pi_{ab}(r)=\frac{d_r}{d^2}\sum_{s=0}^R d_s\chi_a(s)\chi_b(s)\chi_r(s).
\]

这是已有 Racah/covariance 结构的应用，正性来自真实乘积，而不是假定求和每项非负。不需要选择重数幅的相位，但也不包含对称乘积的完整交换相位。

N5m3 的 full CAR 独立检查：π12=(0,5/27,22/27)，π22=(1/75,176/675,98/135)。各使用全部10000个有序 Gaussian 坐标对及所有输出rank，另核 sym strength≤4π。

## 精确矩递推与真正的双标度结论

x_r=r(N+1−r)，C_f=m(N−m)(N+1)/N。χ 的一阶式直接给

E x_r=x_a+x_b−x_a x_b/C_f。

令 L_j(x)=∏_(ell=0)^(j−1)[x−ell(N+1−ell)]，
A_nj=(-n)_j(n−N−1)_j/[(-m)_j(m−N)_j(-N−1)_j j!]。
则 χ_r(n)=Σ_j A_nj L_j(x_r)，M_j=E L_j(x_r) 满足

M0=1，M_n=[χ_a(n)χ_b(n)−Σ_(j<n)A_nj M_j]/A_nn，条件 n≤R。

对 N=p²/λ，m=ζN，a=αp，b=βp，λ固定正、αβ∈[0,1]、ζ处于(0,1)固定紧集。令 t=a+b，S=t−r，Dx=x_t−x_r，

θ=αβλ[1/(ζ(1−ζ))−2]。

保留全部符号计算 n1/2/3，得到

\[
(E D_x/N,E D_x^2/N^2,E D_x^3/N^3)
\longrightarrow(\theta,\theta^2+\theta,\theta^3+3\theta^2+\theta).
\]

有理分母只含N,m，缩放后在上述紧集上远离零；分子对αβ为固定阶多项式，因此该收敛对两个输入rank比例与填充率一致。这不只是若干数值序列拟合。

由 0≤S≤t≤2p，Dx=S(N+1−2t+S)，Dx/(NS)=1+O(p/N) 在S>0上一致。于是实际 S 的前三矩也有上述极限。特别固定 input deficits a=p−u,b=p−v、半填充时：

E S→2λ，Var S→2λ，E S³→8λ³+12λ²+2λ。

这些矩与Poisson(2λ)吻合，但本轮没有证明完整Poisson分布。

## 与 same-J 二次场的接口

对于K≥1，

limsup P(S≥K)≤(θ²+θ)/K²；
limsup E[(1+S)²1_(S≥K)]≤4(θ³+3θ²+θ)/K。

后一式用尾上(1+S)²≤4S³/K，可控制至多线性增长滤波的平方尾。

H=Σ_a H_a 的各 rank 是独立 Gaussian，E H_a²=C_a I。对 a<b，有

E||P_r{H_a,H_b}||²≤4C_aC_b π_ab(r)。

对 a=b，引入独立同rank副本 H_a'，Gaussian chaos decoupling 给

E||P_r W2(H_a²)||²=(1/2)E||P_r{H_a,H_a'}||²≤2C_a²π_aa(r)。

因此没有删去 same-J 的 exchange interference；独立乘积只提供正上界。在项目半填充／中央电荷窗内，旧 ω_a=C_a/C 的 input-deficit factorial tail 与本轮 conditional output tail 一起给未滤波 W2(H²) 的总 deficit 2p−r 紧性。

对真实 memory filter，还须证明其 general-N 重数幅及在固定input窗口内的增长。如果幅至多线性依赖S，则本轮 O(1/K) 平方尾足以支持相应极限交换；当前没有从有限11行表推出这个增长条件。也没有把独立乘积π当完整W2Ecal谱。

## 实际证据与保存

5份最终JSON均passed，不覆盖：

- rank-blocks.json：原N5系统18 block、所有projector精确恒等式、完整及traceless Gram恢复；最大整数中间量710400。
- rank-reconstruction.json：11块精确共线及2个整数样本8个算符重建。
- product-rank-moments.json：4个全rank Racah transforms，合法范围的矩递推；fixed-deficit三个符号极限。
- deficit-limit-checks.json：general α/β/ζ 的三个符号极限；λ1/4，p16/32/64/128、input deficits(0,0)/(1,2)共8组80/120位计算。p16/p128、u=v0的 ES² 上界约.7377264557/.7480984651，极限.75；有限p偏差保留，不是区间舍入认证。
- product-covariance-checks.json：两个有序乘积，full CAR全部20000坐标对独立检查π与sym bound。

保留3份初版代码和失败说明：list高级索引错误；SymPy未展开时结构相等失败（精确展开后全0）；在R2系统误用三阶probe递推（合法范围修正，大系统仍核三阶）。未放宽容差或改原数据。第三方PDF截图工具cache miss，原文文本(34)/(35)实际可读，未声称看过截图。

六页LaTeX首次双遍编译通过、95dpi逐页核版完成。本轮计算全终态：rank session89424、reconstruction60272、moment7602均exit0；初版moment29696 exit1已经记录。完整goal active。

下一步求Q_Ecal在general N的reduced filter amplitudes，明确是否仍collinear、保留重数与交换相位；用三deficit(a,b,output)近顶窗口和本轮尾控制计算实际残量DS权重。然后恢复W4后代、raw Z2联合谱、有限memory两参数协方差与体空间/质量字典。
