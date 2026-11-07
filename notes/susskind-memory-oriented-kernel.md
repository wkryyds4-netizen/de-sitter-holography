# 两种费米子力与完整投影二次记忆谱

2026-11-06，本地。[七页PDF](../output/pdf/de-sitter-susskind-memory-oriented-kernel.pdf)、[LaTeX](../manuscript/susskind-memory-oriented-kernel/main.tex)、[来源](../literature/susskind-memory-oriented-kernel-sources.md)。不commit/push/ZIP。

本轮把上一阶段的等时双chaos残量推进到完整双能量谱。关键是先将输出秩缺额化成湮灭/产生两种费米子力的平均，再证明与Q4的能量分辨cross相消。

## 固定态与算符

标准CAR、完整overlapping complex Gaussian p-body，同J相邻charge，h=H_m、m=N/2，annealed normalized HS，ℏ=1。先λ_N=p²/N→λ>0，q=e^(−4λ)，C_N=Etr_mh²→C=Js²/λ；后小λ。W_j为coupling Gaussian chaos，不是谱q-Wick。

d_i=[H,c_i]，A_−=Σd†d，A_+=Σdd†，O_±=W2A_±，B=W2h²。D X=QX−Σc†Xc、Dbar X=(N−Q)X−ΣcXc†，全部composition先在全Fock作，再限制charge。

## 精确恒等式与强桥接

半填充Casimir Cx=(m+1)D+mDbar−D Dbar，x_a=a(N+1−a)，DH=pH。h_a为neutral rank分量，u=p−a，则

\[
(p-\bar D)h_a=\left(2u+\frac{u(u-1)}{m-p+1}\right)h_a.
\]

D_in h=Σu h_a，E_in=Σu(u−1)h_a/(m−p+1)，U_in,2={h,D_in²h}+2(D_inh)²。barA=(2p−Dbar)H²=A_++{H,(p−Dbar)H}，E_dbl=(2p−D)barA。

旧S_N=Σ(a+b−r)W2Z⁺_ab;r，proxy S_tilde=Σ[x_(a+b)−x_r]W2Z⁺/N。精确有限N恒等式：

\[
N[\widetilde S-(O_-+O_+)/2]
=(1-2p)O_--2pO_+-W_2\{h,D_{\rm in}h\}
-W_2U_{{\rm in},2}+W_2E_{\rm dbl}
+(m-2p)W_2\{h,E_{\rm in}\}.
\]

double mark的界有实际局部理由：对H²原始单项，a_i/b_i数该site原始creation/annihilation数，总和均2p。局部表在(I,c,c†,n,1−n)上D_i=(0,0,c†,n,−n)，Dbar_i=(0,c,0,n−1,1−n)。故2p−D=Σ(a_i−D_i)，2p−Dbar=Σ(b_i−Dbar_i)，两者在singleton site都为零，只在≤K个重叠site作用；两个操作的绝对系数最多const·K²。第一次操作不改变singleton，第二次不产生隐藏p因子。

W2删内部Gaussian配对，固定Wick图中不同p元集有全部固定指数交叠矩Eexp(aK)一致界。展开occupation的绝对系数为poly(K)exp(cK)，所以O±、W2E_dbl及所需高矩有界。旧半填充迹桥接R≤exp(const p²/N)、平衡singleton与fixed-K窗将此界转到固定charge。input factorial tail给D_in^j h各固定矩有界、E_in norm O(N^-1)；Gaussian四阶界及Holder处理乘积。除N后桥接误差O(p/N)=O(1/p)，旧proxy误差也O(1/p)：

\[
S_N=(O_-+O_+)/2+O_{L^2}(p^{-1}).
\]

90个有限矩阵恒等式是独立验证，不替代上述尾控制。

## 两方向的全阶谱矩阵

单mark alpha=2−4·incident-cross，两方向相同。同向双mark delta=+4仅反向双桥，异向delta=−4仅平行双桥。局部determinant normal entry−2、reverse entry+2给此规则；全图q^cr(λδ+λ²αβ)。

两相邻端点交换论证对任意两弧长度成立。异向双桥仅crossed成员留下，所以Σq^crδ=−4qP_ab；同向为4P_ab。alpha beta仍为4[(1−q)/(1+q)]²S_n+16qP_ab/(1+q)，single alpha为2[(1−q)/(1+q)]S_n，S_n=m_(n+4)−2m_(n+2)+m_n。P_ab是已有r=q² probe矩。

c=2λ(1−q)/(1+q)，Y±=O±−cB。减同一守恒背景后

\[
d\nu_{(Y_-,Y_+)}=\begin{pmatrix}w&x\\x&w\end{pmatrix}d\mu_{q^2},
\quad w=4C^2\lambda+\frac{16C^2\lambda^2q}{1+q},
\quad x=-4C^2\lambda q+\frac{16C^2\lambda^2q}{1+q}.
\]

w±x严格正。v=2λ(1−q)+16λ²q/(1+q)，d=2λ(1+q)，rho=v/(v+d)。真正finite-N Pi投影span(I,B,O−,JhF)；旧Gram连续性与JhF系数→0给

\[
T_N=(1-\Pi_N)S_N=(Y_-+Y_+)/2-\rho Y_-+o_{L^2}(1),
\quad d\nu_T=V_{\rm new}\mu_{q^2},
\]

Vnew=C²vd/(v+d)=Js4(1+q)(1−q²+8λq)/[λ(1+q+4λq)]。ν_T,Y−=0是能量分辨结论。

## Q4交叉与完整结果

复用旧Q4整数q多项式权重，四槽对两槽Y+，k=0/1/2桥、未配端点A/B色和两弧split共128模板。每orbit 4!2!/k!；6072字全部λ/λ²多项式相消。此次唯一改变是已独立验证的reverse delta=−4·平行双桥。任意外部spectator在同模板orbit中公共，故证书覆盖任意固定双能量矩。旧Q4–Y−证书加新Q4–Y+证书给ν_Q4,T=0。

有限N正matrix spectral measure的对角矩给紧性及单项式一致可积，off-diagonal由spectral CS控制；候选compact support确定全部矩。强替换对谱测度TV误差≤||X−X'||₂(||X||₂+||X'||₂)，故无需为L²差再强求全部能量矩。

hatPi是旧八方向span(I,B,O−,JhF,E0,M0,S0,Z0)实际finite-N投影。X_N=(1−hatPi)mZ2，沿旧raw reduction得到

\[
\boxed{d\nu_X=V_4\mu_{q^4}+V_2\mu_{q^2},\qquad C_X(t)=V_4G_{q^4}(t)+V_2G_{q^2}(t).}
\]

V4见旧quartic及稿中显式式，V2=(8λCq²)²Vnew。ν_X,Y−=0，包括全部固定左右能量后代。正μ为已有q-Hermite/Mehler核；指定投影后的raw系数没有DS静态atom。

后小λ：V2→1024Js8，V4~768Js8/(5λ)，V2/(V4+V2)~20λ/3，归一X核→sech⁸(2Js t)。20/3只是低chaos权重份额首项，不能当完整time核一阶修正；G_q4本身也有形变。

## 已执行验证与下一步

4JSON：oriented-marks（1688minor、640CAR字/2560比较，115非零double mark/3负）；force-bridge（5系统15非零integer Hermitian样本，90恒等式）；reverse-cross（128模板6072字、5符号式及2小λ系数）；spectral-checks（9双能量矩、3λ两核混合）。Jacobi截断解析界<1e−13、+12 refinement；λ.25/.1另128/192-node Mehler，实际最大差3.38e−14。λ.025不是asymptotic fraction的精确近似，保留有限值。依赖版本、种子、输入hash全部在各JSON；旧科学文件保持原样。

七页LaTeX实际双遍编译、95dpi逐页核版。新计算首次全部通过，无覆盖旧JSON。来源核对继续区分ordinary coupling chaos与谱q-Wick。仅本地保存，无commit/push/ZIP。

下一步：X与旧三阶memory场L的混合谱。已完成X–Y−及所有能量后代，无需重做；仍需Y+–R、Y+–F、Q4–F等接口，不能从等时或不同chaos推出0。随后中央charge窗的uniform扩展、有限两参数memory余项，以及Susskind体singlet/介子空间质量字典。
