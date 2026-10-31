# 创建部分的直接位置配对与联合绝对收敛

日期：2026-10-29。正式证明及全部记号见 [LaTeX position.tex](../manuscript/position.tex)，编入主稿第6节。本文是中文推导及边界记录；它更新此前“未调节位置积分全部开放”的状态，但只解决以下明确的创建配对。

模型固定 dS₂、L=1、G_N→0、ν>0、σ>0、s*>0 与 proper SO⁺(2,1)。同一原模式、精确 m(s)=√(s²+1/4)、实际observer逐次投影和负号时间演化均保留。内积对第一变量共轭线性。文献核查见[来源记录](../literature/stage3-position-action-sources.md)。

## 1. 输入与测试域

已证明的输入是：

- [真实种子共同域](stage3-positive-frequency-orbit-core.md)：紧质量包的双端Schwartz径向尾及轨道积分。
- [分数阶轨道界](stage3-ambient-spectral-test-core.md)：α=1/2时，产品轨道范数受两个 b_α 控制。
- [最大弱谱空间](stage3-maximal-weak-frequency-space.md)：物理态 k 具体嵌入环境 FP 截面；Rk∈L¹loc 且对轨道测试绝对可配对。
- [弱场构造](stage3-regulated-local-descent.md)：各粒子层有界创建块 C_N，范数至多√((N+1)𝔐)，原模式识别和共同cutoff下降。

取 P_N 为有限原observer产品包与有限物质占据基向量的有限和。左磁指标0，右observer与物质总磁指标0，两个质量包均在(s*,∞)内光滑紧支撑。这就是此前生成定理的小种子类；B P_N在完整物理N粒子空间中稠密，不要求P_N在完整运动学空间稠密。

E_R(u)=exp(iuH_R)保持该种子类；物理 U_N(u)在已识别的实际右质量a上乘 exp(iu m(a))，N=0时a=d。它酉且满足 U_N(u)Bψ=B E_R(u)ψ。这一步只用物理谱识别，不依赖局部场的存在。

## 2. 转置后的测试是真正轨道测试

物质湮灭部分记为

\[
\mathfrak a_\nu(x)v=\sum_n\psi^\nu_n(x)a_nv .
\]

对有限占据模式的固定 v，此和有限。固定 n 的原模式双端尾是 O((cosh t)^(-1/2))，紧t区间光滑，故全局有界。于是

\[
\|\mathfrak a_\nu(x)v\|\le
\sum_{n\in F_v}\sup_x|\psi^\nu_n(x)|\,\|a_nv\|=:c_v<\infty .
\]

这里没有随最大磁指标一致的界。定义

\[
\zeta_{\phi,u}(y)=\mathfrak a_\nu(x_R)(E_R(u)\phi)(y),
\qquad \phi\in P_{N+1}.
\]

其范数逐点受原observer产品包模长乘 c_v 控制。untwist 中物质表示酉，故这个支配直接传给轨道绝对范数。

对紧质量包 f，质量相位满足

\[
\sum_{j=0}^2\|\partial_s^j(e^{iu m(s)}f(s))\|_1
\le C_f(1+|u|)^2 .
\]

双端模式展开两次分部积分给加权径向一致界

\[
\sup_t(1+|t|)^2\sqrt{2\pi\cosh t}\,
|u_{n,e^{ium}f}(t,0)|\le C_f(1+|u|)^2 .
\]

α=1/2时尾部权为(1+|t|)^(-3/2)，可积。径向L²范数受谱相位酉性控制。因此

\[
\|a_{\zeta_{\phi,u}}\|_2\le C_\phi(1+|u|)^2,\quad
\|a_{E_R(u)\psi}\|_2\le C_\psi(1+|u|)^2 .
\tag{D1}
\]

ζ不必正频，但属于环境D_orb。这正是使用完整环境位置测试而非只在物理子空间内转置的地方。

## 3. 先独立定义直接积分

以实正规化 Gaussian g_(σ,τ)(u)取中心τ、宽σ。对任意物理输入k定义

\[
\mathcal I_N^+(\phi,k)=
\int_{\mathbb R}du\,g_{\sigma,\tau}(u)
\int_Ydy\,\langle\zeta_{\phi,u}(y),R U_N(u)k(y)\rangle .
\tag{D2}
\]

输入和测试均已限制实际observer谱，D2是所选投影场的弱位置配对，不是声称原位置乘法自动保持正频。两边时间相位给 exp(-iu[m(a)-m(b)])，与原(5.6)一致。

R的绝对配对界与D1首先给

\[
\int du\,g_{\sigma,\tau}(u)\int_Ydy\,
\|\zeta_{\phi,u}(y)\|\|R U_N(u)k(y)\|
\le C_\phi\|k\|_{\rm FP}\int du\,g_{\sigma,\tau}(u)(1+|u|)^2<\infty .
\tag{D3}
\]

这是在证明算符身份之前得到的绝对收敛。U_N(u)强连续、R连续进L¹loc，可选(u,y)联合可测代表；也可由种子逼近在每个紧u区间和紧位置集的L¹极限构造。不同几乎处处代表不改变D2。

## 4. 真正允许交换次序的界

先取k=Bψ，ψ∈P_N，ψ_u=E_R(u)ψ。coarea和Haar换元在非负范数级别给

\[
\int_Gdg\int_Ydy\,\|\zeta_{\phi,u}(y)\|\|U_N(g)\psi_u(y)\|
=\int_Zd\mu\,a_{\zeta_{\phi,u}}a_{\psi_u}
\le C_\phi C_\psi(1+|u|)^4 .
\tag{D4}
\]

Gaussian有全部多项式矩，故u、g、y三重积分联合绝对可积。可先取群平均成R U_N(u)Bψ，也可先做u和位置积分。这比此前“先恢复有界算符，再取群cutoff极限”更强，但只针对已写出的测试配对。

固定g时，输出测试中的物质湮灭只留有限个物质基向量；它们和变换后的输入物质态取有限个内积。左角积分选一个输入磁指标，右角积分按 n_b=n_a+n 选有限个输入observer模式。相应质量系数是原紧包乘固定g的光滑主系列矩阵系数。

因此固定g只需已经证明的有限三模式积分；Gaussian积分产生原h_(σ,τ)(a,b)。这里不把群变换后的整个无限模式和直接移入时空积分。对既有有界C_N，从原Hilbert基逼近输入，有限选择规则给同一个矩阵元。得到

\[
\mathcal I_N^+(\phi,B\psi)
=\int_Gdg\,\langle\phi,(1_L\otimes C_N)U_N(g)\psi\rangle_{\rm kin}
=\langle B\phi,\mathsf C_N B\psi\rangle_{\rm FP}.
\tag{D5}
\]

D3给固定φ时输入物理范数连续，B P_N稠密，创建块已有界。故对所有物理k，

\[
\boxed{\mathcal I_N^+(\phi,k)
=\langle B\phi,\mathsf C_N k\rangle_{\rm FP}.}
\tag{D6}
\]

这进一步给真正的物理测试范数界√((N+1)𝔐)‖Bφ‖‖k‖，并证明Bφ=0时D2为零。不能先把依赖具体测试代表的Cφ当成这个物理范数常数。

## 5. 位置cutoff误差与结论边界

0≤ξ_j≤1为紧位置cutoff、逐点趋1。插入ξ_j后的误差至多

\[
\|k\|_{\rm FP}\int du\,g_{\sigma,\tau}(u)
\|a_{(1-\xi_j)\zeta_{\phi,u}}\|_2\longrightarrow0 .
\tag{D7}
\]

固定u时用轨道L²支配收敛，再以D1作u积分的支配。对固定测试φ，这在全部物理输入单位球上一致。时间截断也由D3的Gaussian尾控制；没有给出未经推导的具体位置收敛速率。

得到的是创建块的直接位置弱配对，覆盖任意物理输入，测试类的物理像稠密，且无需辅助位置/群cutoff。它没有移除Gaussian或observer硬下限，不是逐点Hilbert向量乘法。

物理湮灭块作为C_N*已经存在；但若要其独立的直接位置积分，将场转到测试侧会出现 a†(x)v 的无限模式和。上面的有限和估计不适用。后续需研究积分后的one-particle核、适当物质Sobolev测试或其他分布域，不能只引用伴随就宣布同样的位置积分收敛。两个物理不变态在未规约全Y上相乘积分，也不由D6授权。

## 6. 本轮自查

- 相位：测试侧E_R(u)处于内积第一项，贡献exp(-iu m(a))；输入贡献exp(iu m(b))。
- 物理时间：使用实际observer a，不替换为总表示d；N=0另取a=d。
- 群体积：只群平均输入；测试是原种子，未积两个不变态。
- 可积性：D3是任意物理输入；D4只在真实种子输入使用。未把两类估计混同。
- Bose归一化：使用标准占据基的a_n和已有C_N；不新增√N或跨species因子。
- 结论级别：解析域与积分识别；本轮无新数值扫描，不以编译通过证明数学命题，不作未经充分查新的原创性声明。
