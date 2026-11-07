# 完整双局域的单选态均值

2026-11-08。本轮把上一轮有限词投影候选推进为任意固定次数的证明，并借助既有实际场尾与本轮投影尾，完成完整双局域的一阶矩表示。[七页 LaTeX](../manuscript/susskind-selected-bilocal-mean/main.tex)、[PDF](../output/pdf/de-sitter-susskind-selected-bilocal-mean.pdf)、[来源接口](../literature/susskind-selected-bilocal-mean-sources.md)。当前仅本地，不新增 Git/issue。

## 对象和结论

模型沿用同 J、overlapping complex Gaussian p-body CAR。先有限 N exact pinching，再半填充 fixed-lambda double scaling；本轮不取小 lambda。C=Js²/lambda，q=exp(-4lambda)，h=H/√C，W=√lambda U；h,U是正交单位方差qGaussian。本文 t=√C t_phys、T=√C T_phys，前稿 string变量 u=2√lambda T；hbar=1。rho=1−√q(1+lambda)+√q(h−W)² 是已接实际制备的单选态密度。

令 P 为共同迹 Hilbert 空间到 h,U 多项式闭包的正交投影。C=1下，S_v=(vp−L)Z_v，T_(v,w)=S_v diamond Z_w，n=v+w。一般式为

P T_(v,w)=kappa_v H_n(h)+√lambda Σ_(i=0)^(v−1) Wick_q(h^i U h^(n−1−i))，

kappa_v=lambda[v²−4q∂q log([v]_q!)]。这不是整个高阶场等于该多项式。每个中心化纯 charged J_d 的投影为零。

证明利用任意有限 disjoint D/L 源的已证局域性：单槽W的singleton中心化为零，W-W只在同弦两端相连，L-W只有W弦接入L块时系数为1，D-W为0，三阶以上源项消失。唯一L标记要么独立贡献已知Mahonian能量回归，要么连接一个W，恰好变成含一个W的whole-Wick块。其余W两两配对，故结论对任意测试词成立。

## 高阶递推与两能量表示

有限core上的左右h乘法可转移到测试词；不同coupling Gaussian chaos正交，whole-Wick图的外出弦计数给投影与chaos/OU相容性。写 P eta_j=g_j+√lambda E_j，保留全部低chaos，得到

E_(j+1)=hE_j−[j]_q E_(j−1)−4lambda Σ_(v<j)q^v Σ_(i<v) W_(i,j−2−i)，E0=E1=E2=0；

P N_j=√lambda Xi_j，Xi_j=4lambda q^j Σ_(i<j) W_(i,j−1−i)+(q^j−q^Ncal)E_j。

其中 W_(a,b)=Wick_q(h^a U h^b)，其degree为a+b+1。左乘h的三项分别W_(a+1,b)、[a]W_(a−1,b)、q^(a+1)[b]W_(a,b−1)。Xi3含4lambda q²(1−q²)U，不能丢掉。

以 dnu=K_q(x,y)dmu_q(x)dmu_q(y) 表示one-U空间，h^a U h^b对应x^a y^b。W_(a,b)的符号P_ab正是Casper–Kolb–Yakimov双变量连续q-Hermite的r=q规范；这里只用其递推/显式式，不照搬另一正交测度。

c_j(s)=∫exp(isx)H_j(x)dmu/[j]_q!，Psi_s=Σc_j(s)Xi_j(x,y)。完整实际记忆中charged部分的单选态均值为零；neutral部分给

sigma(alpha_T Fhat(s))=−2√q lambda ∫(x+y) exp[iT(x−y)−is(x+y)/2] Psi_s dnu。

实际Volterra滤波B(T,t)=−∫_0^t ds∫_−(t−s)/2^((t−s)/2) dv alpha_(T+v)Fhat(s)，t≥0。因此完整均值是

sigma B(T,t)=2√q lambda ∫(x+y)exp[iT(x−y)] Σ_j Xi_j L_j dnu，

L_j=∫_0^t(t−s)sinc((t−s)(x−y)/2)exp[−is(x+y)/2]c_j(s)ds。

物理Fhat=F_phys(s/√C)/C，两个时间微分抵消C，B不额外乘C。

## 尾界与短时结果

利用H_j实根上的插值和Hermite–Genocchi公式，直接得到所有实s的 |c_j(s)|≤|s|^j/j!。这是标准插值工具的应用。r=(1−q)^−1/2，Fq=4lambda q/(1−q)，one-U向量范数≤r^degree给 ||Xi_j||₂≤4lambda j(qr)^j+Fq(3r)^j。稿内式(18)给可计算的级数尾R_J(S)；均值误差≤t²||rho||₂ R_J(S)，||rho||₂²=1+q(1+q)(1+lambda)²，对实中心时间T一致。实际finiteN到局域场的识别仍用旧统一尾及sector/preparation桥，不能只凭极限场的新界倒推。

记 Hq(T)=∫(x+y)²exp[iT(x−y)]dnu，则

sigma B(T,t)=[√q lambda²q(1−q)/3]Hq(T)t⁴+O(t⁶)，

sigma B(0,t)=(2/3)√q lambda²q(1−q²)t⁴+O(t⁶)。

这是正的首项。能量单独密度对完整该一阶矩给零。高阶源的自伴性与parity保证这里的even时间结构。

## 执行范围

projection-checks首次成功：196 prefix、40 charged、44 chaos、32原始descendant、8完整selected边界，共320项精确检查。形式有理点(q,lambda)=(1/16,1/3),(16/81,1/5)，不在物理指数曲线。独立时间驱动核72个CKY显式多项式、10个原始场force/commutator Taylor系数，并接bilocal t4/t6。旧114项探索输出只保留未重跑。

物理q=.25,.5取lambda=−log(q)/4，共12个(t,T)点；(J,L,ns)从(24,48,24)到(32,64,32)，最大复差2.78e−17以下，refined解析degree尾最大2.79e−17。相消因子最大约37.7，数值为float64/complex128诊断；不把加密吻合或解析degree尾说成包含roundoff、能谱与时间求积的总认证。证据清单另外从保存值导出相对差与相对虚部，原成功输出不改写。q=.25,t1,T0均值约.0074594258627641，T=.5约.0066385092055227。

下一步是两插入投影：L-L、D-L连接项及one-U正交补会贡献方差，不能将P(XY)替成P(X)P(Y)。低阶必须恢复旧neutral-selected-response。完整选态双局域高点、空间/质量生成元、有限G仍开放。本轮长期goal为progress。
