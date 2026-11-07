# 同一 J 的协变大小涨落：从选态代数接回平稳动力学

2026-11-08。模型同前，standard CAR、半填充固定Slater annealed态、overlapping Hermitian Gaussian p-body、p²/N→λ>0、h=H/√(wD)，q=e^(−4λ)、γ=q^(1/4)、s=√q。新加入实际微观u_N=(p−L−λ_N)h_N/√λ_N，L是full-Fock half-Clifford-size，先作用后限制中央sector。它来自旧bilocal的W1，同一J，非新增独立probe。h时间无量纲，物理string能量乘J_s/√λ。

## 从局部源得到全部混合词

每个C_A†C_B单项有r个重叠site；p−L在这些site将n替成I/2后求和。形式源exp[z(p−L)]在重叠site为n+(exp(z)−1)I/2。每endpoint一个nilpotent z，所以仅n+zI/2。singletons无源。保留旧Slater局部CAR与Jordan符号后，一个G弦自身双membership在空/占据求和的系数为1、z_i/2、z_j/2、z_i z_j/2；selected弦不允许own overlap。不同弦双membership无源。

旧归一化defect证明在有限nilpotent系数范数下照样成立：singleton不变，所有局部系数被固定B控制，activity m(p/m)^|S|；size≥3消失、backward需≥p个defect亦消失。2λ乘own-double系数，减λ每W endpoint，得到G边上的exp(λ z_i z_j)，selected边W端的exp(−λz_i)。除√λ即完整四字母规则：HH=UU=1、HU=0、HM混合γ、UM混合−√λγ、MM†=1，只保留selected forward，乘q^(crossings+inner-selected nestings)。奇词零。这是任意固定word的解析推导，不由六阶检查外推。

u_N是Hermitian线性Gaussian且U(N)-equivariant；E u_N²=v_N I，二点规则给v_N→1。因此实Gaussian系数row/column norm一致有界，接前轮逐Wick图noise-slot收缩与阶乘尾，实际有限N的所有固定插入在紧复时间窗收敛。全joint quenched集中未另证明。

## 新方向和具体混合矩

记c=−√λγ。φ(hu)=0、φ(u²)=1、φ(uM†)=c。新H/U/P传递的各grade正交；旧hM cyclic空间由H/P字串穷尽，故uΩ在旧空间的投影只需第一层：

\[
P_{hM}u\Omega=\frac{\sqrt\lambda s}{1-s}h\Omega-
\frac{\sqrt\lambda\gamma}{1-s}M^\dagger\Omega,
\quad \|u\Omega-P_{hM}u\Omega\|^2=1-\frac{\lambda}{e^{2\lambda}-1}>0.
\]

余量被h-annihilator杀掉，是新的Hamiltonian primitive；λ→0时范数平方趋1/2，并非可丢的小方向。四点：φ(Mu²M†)=1+λs(1+q)，φ(uM†Mu)=λs，φ(MuM†u)=q(1+λs)，φ(MM†u²)=1。

## 全传递与平稳态

字串按oldest→newest，append H/U/P分别是a_H†/a_U†/b†。令c_H=γ、c_U=−√λγ。a_i删同色第j端权q^(n−j)，删P权c_i q^(n−1)，异色零。b删H/U/P权分别c_H/c_U/1乘q^(n−1)。h=a_H+a_H†、u=a_U+a_U†、M=b。Gram来自切图/微观正性，字串不正交。a_i a_j†−q a_j†a_i=δij；a_i b†−q b†a_i=b a_i†−q a_i†b=c_i q^N；bb†−q b†b=q^N。

finite grades≤3^n，creator norm界[n+1]_q与(n+1)q^n使h,u bounded而M仍compact。h沿primitive q-Jacobi towers有pureAC谱；qcomm幂抽出PΩ，joint weak closure仍B(H)。这些事实并不要求证明每grade恰3^n。

令D=C*(1,h,u)、A+=C*(1,h,u,M)=D+K(H)。φ在D上是标准两色qGaussian trace。于是lim_|t|∞φ(alpha_t(B+K))=φ(B)=τ∞(B+K)，B∈D、K compact。该limit保证定义无歧义、正且归一；τ∞为A+上的tracial stationary state，杀M但不杀u。它的GNS是Γ_q(R²)，按Ricard已知定理是II1 factor。未声称整个A+/K与某个特定C* algebra同构；未分类所有stationary states。相对原Type I表示，该态非normal，不能当B(H)普通density。

τ∞([h,u]†[h,u])=2(1−q)>0。一般有限β可取τβ(A)=τ∞(e^(−βh)A)/Zβ，是真正α-KMS态；β数值尚未赋予de Sitter温度。标准qMehler核K_q(E,E′)=Σq^n p_n(E)p_n(E′)，C_u(t)=∫e^{it(E′−E)}K_q μ_q²=Σq^n|〈n|e^{ith}|0〉|²，二阶为1−(1−q)t²。finite-G解释仍未建立。

## 实际rank/bilocal范围

旧size-bridge在精确半填充给||(R−L)h||_(2,m)=O(N^-1)。该差也是同J线性Gaussian equivariant场，row/column covariance是O(N^-2)I。因此用相同逐图收缩，可把u替换成(p−R−λ_N)h/√λ_N，所有固定Slater mixed words和紧时间关联差趋零。这个结论补上旧tracial rank/size桥到选态混合词的接口。

尚未把pinched实际neutral slope或完整bilocal放入这里的全joint定理：全迹S_r强差不能无证明在M激发态中使用，pinching也破坏linear-Gaussian条件。下一任务应从实际线性未pinch slope与rank谱比较起，再控制pinched项在选态mixed words中的质量，随后接J2和完整bilocal。已有完整full-trace bilocal四点不是待重做任务。

## 检查与封存

新site-checks：136704有限CAR矩阵比较、7880局部source系数；word-checks：8192六字、1458旧限制、720代数、62投影正交、1638Gram元素。精确有理γ,d独立测试点只核代数，不冒充物理q=e^(−4d²)曲线。另三物理q、18时间行比较80/120 Jacobi、解析cutoff界与曲率；浮点非interval。成功science不重跑，三个JSON/代码/参数/版本/hash与最终PDF独立绑定。当前仅本地，无Git/issue操作。
