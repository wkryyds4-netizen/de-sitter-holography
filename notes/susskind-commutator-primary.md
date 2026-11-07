# 混合增长、排序反对称复合态与三点顶点

2026-11-06，本地阶段。[六页PDF](../output/pdf/de-sitter-susskind-commutator-primary.pdf)、[LaTeX](../manuscript/susskind-commutator-primary/main.tex)。不启代理、不commit/push/ZIP。

## 规范与已有理论接入

标准CAR，含重叠集合的complex Gaussian p-body H_s，退火归一trace、β=μchem=0、ℏ=1。先DS固定λ=p²/N>0，再固定string时间取小λ；h=4λ，q=e^(−h)，Hsrc=Hs/(2Js)，τHsrc²=1/h，u=2Js t。沿用归一X和Rhat，d=2、e=3。

原文Lin–Stanford(205)保留ΔVΔW，直接用于不等species。β→0即c→1后，crossed connected/h为

−de g_d(x)g_e(y) cosh A/(cosh x cosh y)，x=u1−u3，y=u2−u4，A=u1+u3−u2−u4，g_d=sech^(2d)。扣除的是完整有限h两点乘积。G.2(388)给uncrossed energy exchange/h=−de g_d(x)g_e(y)tanh x tanh y；组合commutator时注意反转pair方向。

所以固定时间 **D23(t)=48λcosh(4Js t)+O(λ²)**。这是已有半经典展开的代入，不再仅靠上一轮t6猜函数；O(λ²)未认证为统一数值余项，也不用于随log(1/λ)增长的时间。

更一般B_h(u,v)=−i[X(u),Rhat(v)]/√h的二点极限为

12 sech⁵(u−u′) sech⁷(v−v′) [cosh(u+u′−v−v′)+sinh(u−u′)sinh(v−v′)]。

对角给12cosh2(u−v)，共同时间线u=v、u′=v′给12sech¹⁰(u−u′)。四时间射线、20系数和60 Laurent相消独立核对该式。

## 精确复合态与乘法

K=−i[X,Rhat]/√[2(1−q⁶)]，S={X,Rhat}/√[2(1+q⁶)]。两个排序态的Gram是[[1,q⁶],[q⁶,1]]，因此K/S正交、各范数1，而且

XRhat=√[(1+q⁶)/2] S+i√[(1−q⁶)/2] K。

排序奇分量强度(1−q⁶)/2~12λ。原文§4.8及F.1已经明确讨论这类commutator primary；此处是d2/d3的具体应用，不宣称新发现该家族。k0的归一不需要一般k的猜测公式。

任意外部H只放在两个复合block之间，不在block内。交换block内部XR次序只改变互交叉q⁶，其余H交叉不变。因此所有固定a,b满足τ(KH^aKH^b)=τ(SH^aSH^b)=σ^(a+b)P_ab(q,q⁵)，τ(KH^aSH^b)=0。紧矩确定性给完整Mehler两能量测度r=q⁵及零cross。K时间核在小λ趋sech¹⁰(2Js t)。25两弧矩至8个外H核验。

同样交换三时间word的首个相邻block，证明精确

τ[K(t0)X(t1)Rhat(t2)] = i√[(1−q⁶)/2] τ[Rhat(t0)X(t0)X(t1)Rhat(t2)]。

右边为已有uncrossed四点。固定时间领先为i√(12λ)sech⁴(u0−u1)sech⁶(u0−u2)+O(λ^(3/2))。50三弧矩至6个外H验证；等时顶点严格i√[(1−q⁶)/2]，反转X/R变号。

## 高点与体构建

精确τK⁴=2+q²⁵+(1−q⁶)²(1+q⁴+q⁹)/2。最后项>0，领先54h²=864λ²，不能把K作为任意高点独立q⁵探针。τK²/K⁴/K⁶的小h极限1/3/15只是已算矩，不是完整joint Gaussian定理。

K是label5的额外排序副本，适用此前已证明的规范化fixed-descendant收缩，得到D⁺_5和正Q_h；极限Q种子谱Gamma(shape10,rate2)。它与普通两粒子乘积具有相同传播label而不同乘法。这为体映射提供额外约束，尚不选出物理质量。

下一步优先把两种排序的整套相对态内积写成实际算符。LS F.1的领先norm因子提示h[k(k+9)+12]，即h(C_total−C_2−C_3)。必须区分原文一般有限q归一公式的证据层次，优先从原Gram/cross递推证明领先二次型。然后检查已有Q和平移候选如何作用在这额外扇区，及乘法顶点与质量算符是否相容，不再反复单核拟合。

## 计算证据

exact.json：25两弧、三个偶矩、20covariance Taylor、60 Laurent相消。vertex.json：50三弧、Gram、简化第四矩。numerics.json：λ=.02/.01/.005、cutoff244/497/1032、四时间点，各cutoff+16；最大绝对变化2.78e−15。四个50位原非对称BINT和式与正Jacobi形式比较，最大绝对1.97e−17、相对5.94e−15。非对称和条件数≤1.000004；D=2(1−F)相消条件数≤16.69。解析Jacobi尾界在每组记录中；浮点控制非区间认证。

λ=.005,t=.35的D/λ=87.789905，较领先103.243126低14.97%；两点Richardson仍低3.52%。保留真实偏差，不把表示精度当渐近误差。初次LaTeX因长行溢出被构建检查拦截，保留initial log，拆开数学列表后重编通过并六页核版。
