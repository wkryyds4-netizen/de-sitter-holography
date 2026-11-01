# 正内积局部表示与尚未解决的五边形复现

2026-10-31 至 2026-11-01。Alan Day / Chengheng。[LaTeX](../manuscript/unitary-flip/main.tex)、[PDF](../output/pdf/de-sitter-unitary-flip.pdf)、[来源](../literature/stage4-unitary-flip-sources.md)。

本轮从CLP真实缝合缺口转向其§6引用的Kim–Scarinci已有构造。它提供明确Lebesgue正内积和局部酉算符，但不是已完成的有质量观察者字典。标准公式复现不能包装成新全息结果。

对hbar>0，z=x+i hbar y，Phi^{i hbar}(z)=psi_{exp(-pi hbar)}(exp z)/psi_{exp(-pi/hbar)}(exp[z/(i hbar)])；F=Phi/conj(Phi)几乎处处模长1。孤立零极点可任意指定模长1的值，不影响连续乘法算符。每条边x、y强对易时可用联合谱定理定义F(x,y)。CLP double-sine通过其(3.48)与Phi相差明确Gaussian前因子，不能据此把其固定拓扑模平方认作同一个酉算符。

四穿孔球取边12,13,14,23,24,34，U行=(1,0),(0,1),(-1,-1),(-1,-1),(0,1),(1,0)，J=[[0,1],[-1,0]]；定向交换矩阵E=UJU^T，rank E=2。V为端点关联矩阵，V^T U=0。L²(R²)中p=-pi i derivative，x=Up、y=UJt给[x_i,y_j]=pi i Eij及零穿孔约束。R=V(V^TV)^-1后，x→x+R xi、y→y+R eta给任意实中心值，仍为自伴Weyl表示。此为标准Heisenberg表示的显式实例与代数扩展；尚未证明xi/eta等于锥角或观察者质量，也未完成改变三角剖分的全套交织子。

测试采用等价规范x1=x,y1=y,x2=p_y,y2=p_x。令A=F(x,y)、B=F(p_y,p_x)、M=exp(i xy/pi)，则C=F(x+p_y,y+p_x)=M^-1 B M。需复现原文AB=BCA，而不只是A、B各自酉。态取归一exp[-(x²+y²)/(2pi)]。两个次序ABψ、BAψ的Gram对角为1、非对角z，距离D=2−2Re z。Gram正性是这个明确Hilbert模型内事实，不是观察者四点的正性结论。

21组统一协议计算包含FFT的(256,16),(512,24),(1024,32)和Hermite-DVR的64,128,256,512；各取hbar=.5、1和独立Lambda=0对照。高阶FFT的五边形残差为.02955658、.12168478，Hermite为.02946719、.12156956。平直对照分别3.01e-12、4.81e-15。两种方法的离散范数均保持，但正宇宙常数五边形复现未通过。不得把实算残差写成截断已认证的连续反例。

36项标量/代数核查通过，包括55位独立q-Pochhammer、反演式、3个倾斜Barnes积分、精确穿孔矩阵及chirp共轭式。乘积尾界只控制精确算术的函数截断；不控制浮点舍入、FFT环绕、Hermite离散化或算符极限。原始不同态/表示的prototype及失败全部保留，不能把它们的数值与统一协议混表。audit只核数据和算术一致性，其passed不表示五边形passed。

下一步先控制连续空间：用真正的Hermite Galerkin矩阵元替代DVR节点函数演算，以已知乘法模长1从投影范数差估计每步漏出；隔离Phi零极点邻域并给误差界，再结合有界算符的逐步误差传播控制一个投影残差。需要独立实现与原始域/解析延拓复核后，才讨论定理层面的结论。不得将复现困难投射成de Sitter全息不可能。实际质量映射、边界态、全翻转及物理四点缝合仍开放。
