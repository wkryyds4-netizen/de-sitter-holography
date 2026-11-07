# Slater轨道密度的一致矩界与单激发中的精确pinching

2026-11-08。前驱neutral-slope已证明实际归一t_N与u_N在未pinch的全部selected words/time中相接，并给纯协变词内的exact pinching小量。本轮补上单个线性selected制备的有限N估计；不是用极限密度ρ1倒推有限N性质。

## 1. 有限N的一般线性制备命题

H_m=wedge^m C^N，d_m=binom(N,m)，固定Slater |O〉。原p-body Hermitian GUE耦合写成独立实标准Gaussian g_a；物理U(N)在这个实noise空间正交作用R(U)。P_O(J)为齐次线性算符族，要求对O的稳定子K=U(m)×U(N−m)满足P_O(k·J)=k_m P_O(J)k_m†。M†、h、u及固定线性组合均满足；a_N=E||P_O(J)|O〉||²。

定义q_U(J)=U_m†P_O(U·J)|O〉，Z(J):L²(U(N))→H_m的核为sqrt(d_m) q_U(J)。只需有限维matrix-coefficient支撑。R_P=ZZ†=d_m∫q_Uq_U†dU≥0。

耦合旋转不变与Y协变给精确式E〈O|P†Y P|O〉=E tr_m(R_PY)。全U平均加Schur给E ZZ†=a_N I_m。真正需要补的是另一侧E Z†Z≤a_N I；只有E R=aI并不够。

## 2. 权重轨道给第二侧variance

写q_U=Σ_a g_a f_a(U)，Ψ_U=Σ_a e_a⊗f_a(U)=ρ(U^−1)Ψ，其中ρ=R⊗wedge^m、||Ψ||²=a_N。在ordered complex J标签基中noise representation是conjugate Ad(wedge^p U)，与实Hermitian基酉等价，不能误写成同向noise weight。对最大torus D，稳定子协变使ρ(D)Ψ=χ_O(D)Ψ，χ_O=∏_{i∈O}D_ii。

任何含Ψ分量的irrep都含weight (1^m,0^(N−m))。置换矩阵生成其binom(N,m)个不同权重；这些weight空间互相正交，所以每个这样的irrep维数d_α≥d_m。把ρ分解为⊕_α V_α⊗M_α，Schur twirl给

d_m∫|Ψ_U〉〈Ψ_U|dU=⊕_α (d_m/d_α) I_{Vα}⊗Tr_{Vα}|Ψα〉〈Ψα| ≤ a_N I。

multiplicity部分保留，trace上界不需要multiplicity-free。这个frame与E Z†Z是同一个analysis map的TT†和T†T，非零谱相同。于是E Z†Z≤a_N I、Tr(E Z†Z)=d_m a_N。

## 3. Gaussian矩界及选态小量转移

Hermitian dilation S=[[0,Z],[Z†,0]]是线性Gaussian。V=E S²=diag(a_N I,E Z†Z)≤a_N I，TrV=2d_m a_N。Gaussian分部积分及在S本征基中的加权AM–GM给ETrS^(2r)≤(2r−1) ETr[V S^(2r−2)]≤(2r−1)a_N ETrS^(2r−2)。起点TrV，而非ambient dimension，故

E tr_m R_P^r ≤(2r−1)!! a_N^r，r≥1。

这是常见matrixGaussian矩证明在具体orbit variance上的应用，不宣称新Khintchine一般定理。若P=M†，a_N=1精确，特别E tr_m R_M²≤3，对N,p,m都成立（合法p≤min(m,N−m)）。对任意covariant Y，|E〈O|P†Y P|O〉|≤sqrt(3)a_N ||Y||_(2,m)。对covariant X有E||X P|O〉||²≤sqrt(3)a_N ||X||_(4,m)²。

后置协变制备B_N也可处理：R_{B P}=B R_P B†逐样本精确，故||R_BP||r≤||B||_(4r)² ||R_P||_(2r)。所有固定covariant多项式/实际unitary phase已有相应高矩或norm界。这里的顺序是B P|O〉，不能无推导变成P B|O〉。

## 4. 实际pinched斜率在单selected态闭合

前驱central-sector证明X_N=P_H t_N在每个固定S_r趋0，包括全部简并energy blocks。故E||X_N M_N†|O〉||²≤sqrt3 ||X_N||4²→0。实际归一exact-pinched b_N=(1−P_H)t_N因此在此制备态中可替成u。更一般，令W为任意固定covariant word，含h_N、u_N、t_N、b_N及固定实际unitary phases。trace Hölder用高矩证明W(b_N)−W(t_N)在S2趋0，再用R_M的S2界，得到

lim_N φ_N[M_N W(h_N,b_N) M_N†]=τ[ρ1 W(h,u)]，
ρ1=(γh−√λγu)²+[1−sqrt(q)(1+λ)]I。

例如φ_N(M_N b_N² M_N†)→1+λsqrt(q)(1+q)。这是finiteN先exact pinching后的实际极限，不再是先N后finite-window的替代。允许任意有限个b_N插入和固定外部h-time，均位于该M...M†制备夹心内部。

本轮没有覆盖任意多次交错selected插入或M†B|O〉。一般P多项式时Z(J)不再linearGaussian；degree本身不提供这里的矩递推。已完成单激发制备及covariant后置dressing，下一步需多个selected制备或直接bilocal/J2混合结构，不能把该范围写成全joint exact pinching。

## 5. 精确有限检查

seed fixtures (N,m,p)=(4,2,1),(4,2,2),(6,3,1),(6,3,2),(8,4,2)：59个噪声+物理权重等式、122个permutation frame计数；N4两种p各核两种stabilizer变换，共624个矩阵元素，保留noise conjugation和fermion signs。

另用m=p=1的projective Haar积分独立核完整R。对任意Hermitian单粒子J：

R(J)=[N J²−2Tr(J)J+((TrJ)²+TrJ²)/(N+1) I]/[(N−1)(N+2)]。

N2,3,4各两份exact complex fixtures，Haar三阶monomial公式独立展开核58矩阵元素和全部principal minors。N2只有此例同时半填充；其余是有限N一般lemma检查，不是DS渐近。

实Gaussian坐标独立Wick核N2..5的mean I和second moment，共978坐标pair；精确second moments在1与3之间。N2可化R=(χ_3²/3)I，r1..8的exact moments∏_{j<r}(3+2j)/3^r均≤(2r−1)!!。

首轮checker在比较未canonicalize复数表达式时assert失败，改成逐entry simplify及principal minors后首次成功写science JSON。保留初失败诊断；成功后不重跑。研究范围和PDF输入输出hash另由read-only manifest保存，仅本地无Git/issue。
