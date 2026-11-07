# 二次记忆 core 的强滤波极限

2026-11-06，本地。[六页PDF](../output/pdf/de-sitter-susskind-memory-filter-limit.pdf)、[LaTeX](../manuscript/susskind-memory-filter-limit/main.tex)、[来源](../literature/susskind-memory-filter-limit-sources.md)。本轮不commit/push/ZIP。

完成固定半填充下整个 W2 Ecal 的算符 L² 渐近式。关键改变是用正乘积权重的六阶矩控制增长窗口外的 A，不必先求下charge交换二阶Casimir，也不必先证明A在固定输出窗口内紧。前一轮的相干范数抵消仍用于窗口内。

## 模型、态与有限 N 定义

标准CAR、完整overlap complex Gaussian p-body系综，ℏ=1，固定m=N/2的maximally mixed归一迹，annealed HS范数。h=H_m、k=H_(m−1)同耦合，C_N=Etr h²；W2为普通Gaussian二次chaos，B=W2h²，A=Rem(k²)+(2p−m)h²，O=W2A。Ecal=B22^f+[h,[h,A]]/4，为先前quartic memory core。

λ_N=p²/N→λ>0，q_N=exp(−4λ_N)，C_N→Js²/λ。以下C=C_N。领先p项必须保留C_N、q_N，不能在没有收敛速率的条件下把它们直接换成极限常数。

## 正乘积权重的六阶一致界

独立单位方差rank-a/b的Hermitian矩阵U,V给正权π_ab^(M)(r)=Etr|P_r(UV)|²。和为1，支撑r≤a+b。已有covariance composition给Eπ χ_r(n)=χ_a(n)χ_b(n)，其中

\[
\chi_a(n)={}_4F_3[-a,a-N-1,-n,n-N-1;-M,M-N,-N-1;1].
\]

x_j=j(N+1−j)，P_j(x)=∏_(l<j)(x−x_l)，c_nj=(-n)_j(n−N−1)_j/[(-M)_j(M−N)_j(-N−1)_j j!]。则F0=1，F_n=[χ_a(n)χ_b(n)−Σ_(j<n)c_njF_j]/c_nn是精确EπP_n递推。

取ε=N^(−1/2)、a=A0/ε、b=B0/ε、M=1/(2ε²)+δ，Y=ε³x_r。在QQ[A0,B0,δ]上展开至ε^12。scaled c_nn抽去ε^n后的首项为(-4)^n，其余有理分母在紧参数集远离0；不依赖非零A0/B0。因此形式计算给真实一致解析展开。

X=x_(a+b)−x_r。1≤j≤6的Eπ(X/N)^j所有负幂逐一消去，共21项；常数项为Touchard T_j(2A0B0)，与δ无关。T6(t)=t6+15t5+65t4+90t3+31t2+t。故M=m,m−1、全部a,b≤p的Eπ(X/N)^6一致有界。

S=a+b−r≥0，X=S[N+1−2(a+b)+S]≥S(N+1−4p)，最终括号为正。所以两charge的EπS6≤Kλ。这只使用六个矩，不宣称完整Poisson律；δ有界也不等同中央charge窗已完成。

## 增长窗口外可以使用粗界

Z⁺=P_r(h_a⊙h_b)，Z⁻=P_rRem(k_a⊙k_b)，a<b用反对易子，a=b用平方。A_ab;r=Z⁻+(2p−m)Z⁺。n=m+1，Ca−/Ca=(m−p)²n/[m(mn−xa)]≤Kλ，Rem的normalized HS范数增益≤m。

独立对称乘积强度≤4π；a=b的中心平方经Gaussian decoupling保留绝对常数界。上下charge各用正π六阶矩，再用范数三角不等式平方，得到

\[
\sum_{a\le b,r:S>R}\|W_2 A_{ab;r}\|^2
\le K_\lambda C^2m^2/R^6\le K_\lambda C^2p^4/R^6.
\]

两个O(m)项此处分别粗界没有被隐藏：损失p4明确保留，随后让R增长压过它。它不能给固定R的A尾紧性。不同unordered input pairs在二次Gaussian chaos正交；不同output ranks在HS正交。

## g_A 的一致界与强极限

精确g_A/C=U0+Ua κa/C+Ub κb/C+Ur κr/C。四个U完整写于稿件/JSON。a=αp、b=βp、r=ρp，αβ∈[0,1]、ρ∈[0,2]时一致为

\[
(U_0,U_a,U_b,U_r)=(3/2,-1,-1,1/2)
+\lambda_N(8,-5/2,-5/2,3)/p+O(p^{-2}).
\]

κ的binomial product给相邻rank差≤KλC/p，κp/C=q_N+O1/p、κ2p/C=q_N²+O1/p。设u=p−a、v=p−b、γ_N=(3−4q_N+q_N²)/2，则全rank上|gA/C|≤Kλ、|gA/C−γ_N|≤Kλ(1+u+v+S)/p。

前阶段已证明Σ(1+u+v)²||W2A_ab;r||²≤KλC²。误差E_A=Σ(g_A−Cγ_N)W2A_ab;r：窗内用此相干范数界及S≤R；窗外用g_A一致界和上述六阶尾。因此

\[
\|E_A\|^2\le K_\lambda[p^{-2}+R^2/p^2+p^4/R^6].
\]

R=p^(3/4)给||E_A||=Oλ(p^(−1/4))。不需要另截input log窗口；已有加权范数直接处理u/v。

## 完整 W2 core 与投影残量

S_N=Σ(a+b−r)W2Z⁺，D_N=Σ(p−a)h_a，bar s_N=Σ(p−a)ωa→λ。旧H部分加上本轮gA，再用W2{h,D_N}=Linput+2bar s_NB和JhF−4λ_Nq_N C Linput→0，得到

\[
W_2\mathcal E=C[2p(1-q_N)^2+\bar K_N]B+C\gamma_NO
+2(q_N-1)J_hF+8\lambda_NCq_N^2S_N+o_{L^2}(1),
\]

\[
\bar K_N=\lambda_N(1-10q_N+9q_N^2)-16\lambda_N^2q_N^2.
\]

这里只给gA部分O(p^−1/4)速率；旧H余项仍是oL²。所有估计对正λ紧集一致，保留N下标兼容一般λ_N→λ序列。

令Π_N为真正有限N的span(B,O,JhF)正交投影。前阶段Gram非退化，S_N的投影残量范数平方→Vnew。投影精确消除p项，收缩性保留oL²，所以

\[
\lim\|(1-\Pi_N)W_2\mathcal E\|^2
=\frac{64J_s^8q^4(1+q)(1-q^2+8\lambda q)}{\lambda(1+q+4\lambda q)}>0.
\]

先DS再λ→0得到1024Js8。这是指定三个方向之外的等时二次chaos残量，尚不是完整raw Z2、新粒子或动力学谱。

## 计算证据与后续

3JSON首次运行均通过且未覆盖：sixth-moment-certificate（精确多项式21负幂相消/六个常数项）；sixth-moment-checks（4完整π表28矩及48组80/120位exact Fraction后求值）；filter-certificate（130有限式及组装/小λ核对）。

λ.5,u=v0,上charge的p16/32/64/128六阶X/N proxy约75.6892/123.1213/157.8280/178.9193，极限203；对应真实S六阶的有理上界小数约166.4086/180.7793/190.8031/196.6133。最后一组下charge/u1v2的proxy偏差约−37.34285，原始记录保留。非实际大系统S矩，非区间舍入证书。

六页PDF双遍首次通过、95dpi逐页核版，公式/引用无溢出。构建脚本code/build_memory_filter_limit.ps1；输入输出哈希见results/susskind-memory-filter-limit/local-evidence.json。

下一步恢复第四chaos后代的low-chaos组分，组装raw Z2并扩中央charge窗，再求联合双能量谱/finite-memory/体空间质量字典。不再把固定输出紧性或下charge交换二阶Casimir当作本阶段的必经关卡。
