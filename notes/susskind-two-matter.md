# 两 singlet 通道与首个相互作用信号

2026-11-06，本地阶段。见[六页 PDF](../output/pdf/de-sitter-susskind-two-matter.pdf)和[LaTeX](../manuscript/susskind-two-matter/main.tex)。Alan Day、Chengheng 署名；不提交、推送 Git，不生成 ZIP。

## 这轮实际算出了什么

把已有的短时间对易子平方系数，精确拆成相互正交的弦通道。设 h=4λ，q=e^(−h)，r=q²，χ=q⁴=r²，Hs=2J_s Hsrc，τHsrc²=1/h。完整 annealed trace，β=μ_chem=0；先 complex DSSYK DS，再固定阶/时间的小 h。

F=Hsrc X²−2XHsrc X+X²Hsrc，满足 [dot X,X]=2iJ_s F，因此 C(t)=4J_s²t²||F||²+O(t⁴)。精确有限 h 结果：

||F||²={4(1−r)²+2(1−r)⁴+4(1−r²)²}/h。

三项分别来自零物质 H 投影、双物质 k0 外后代、双物质 k1 相对 primary。对 C(t) 而言，它们的领先小 λ 系数为256、0、1024，单位 J_s²λt²。合计1280，与旧完整四点精确短时间计算相同。这里多出的是通道内容，不是另一个 Lyapunov 指数。

## 直接推导

一根 H 相对两个 open X 的位置按 L,M,R 排列。direct Gram G=[[1,r,r²],[r,1,r],[r²,r,1]]，cross Gram Gx=χ[[1,r,r²],[r,r²,r],[r²,r,1]]。特别是交叉 middle-middle 元是χr²。完整同种二开放弦内积为 G+Gx；内部 X 收缩另入零物质 sector。

direct 归一 primary u1=(-r,1+r²,−r)/√(1−r⁴)，两个外 annihilator 都湮灭它，d_eff=5。s0=(1,0,1)，v=(1,−2,1)。严格投影 v=c s0+a u1，其中 c=(1−r)²/(1+r²)，a=−2(1−r²)/√(1−r⁴)。G 与 Gx 中两个方向均正交，u1 的 physical norm²=1−r⁴。因此两通道范数恰为2(1−r)⁴和4(1−r²)²。

Hsrc X²、X²Hsrc 中两 X 收缩各留 Hsrc；XHsrc X 留 rHsrc。故零物质部分=2(1−r)Hsrc，贡献4(1−r)²/h。它与二物质 Wick sector 正交。

k1 的 direct norm 为1，physical norm²=1−e^(−8h)~8h；投影 a²~8h，乘上 F 自带的1/h后留下 O(h)。不能先把趋 null 的相对方向丢掉，再求首阶相互作用。这解释了领先系数中4/5来自何处。该分解对应明确的 Wick/primary 通道，不能自动读作空间散射概率。

## 两物质 Q_h 与相对坐标

LS已有 ordered 两物质表示 R2(d1,d2)=⊕_k R1(Dk)，Dk=d1+d2+k。每个归一 primary 的存在可由 grade 维数与外 raising image 的正交补构造，不依赖一般闭式归一化猜想。对每个固定k套用上轮范数证明：Q_h=E+(K−−K+)/(2i)正、自伴，固定k,n反对称后代收缩为D+_(Dk)。这只控制所选固定后代扇区，不能抹掉其余奇异重标度方向。

标准正离散系列张量积给同一极限直和。用P=p1+p2、x=p1/P，酉坐标变换 Uf=√P f(Px,P(1−x))，测度dP dx。相对基χk=N_k^(−1/2) x^(d1−1/2)(1−x)^(d2−1/2) Jacobi_k^(2d2−1,2d1−1)(2x−1)。channel basis映为φ_n^(Dk)(P)χk(x)，Q=P，B=i(P∂P+1/2)。正文给N_k及相对 Sturm–Liouville算符；三个一般参数微分恒等式独立核验了完整变换。这是极限表示的酉实现，不是finite-h coproduct忽然变成primitive sum。

若外加自由两粒子 Poincare 动力学，M²=m²/[x(1−x)]。d1=d2=2时，相对空间的〈χ0|M²/m²|χ0〉=14/3，〈χ0|M²/m²|χ2〉=√11/3≠0。这里不是含不同径向波函数的完整Ωk矩阵元。即使自由质量也混合k；离散sl2标签不等于介子质量。一般正的相对算符 Mx²配P−=Mx²/(μP)都可满足运动学，需独立微观动力学来选。

## 文献公式的使用范围

LS附录F的ψk00归一化和Rk跨序重叠通式由原文猜测并核到k25，不能在我们的稿里升格为已证全k定理。我们的k0,k1结果直接解析证明；新代码把k0..4、三组参数同该通式核对。Rk的模小于1，描述内积而不是unitary S matrix；fixed h的顺序信息不能直接换成普通Bose交换投影。

Xu §3.2已实际读完相关factorization/measure公式：普通张量空间须更改内积才给isometry。这支持保留Gram权，不给独立粒子动量算符。Koelink–Van der Jeugt及Van der Jeugt–Jagannathan给标准Jacobi/Laguerre耦合表示，我们明确复用。

## 验证、修正与下一问题

validation.json保存38项精确恒等式，包括完整六字母Wick配对，与Gram通道计算独立；tensor-check.json另3项一般参数微分检查。三组(h,d1,d2)=(.05,2,2),(.01,2,2),(.02,1.3,2.1)，k0..4，用85/115位。最大相对精度变化2.42e−78，115位cross overlap残差≤1.71e−108；这些是有限点经验精度检查，不是全k证明。

初版mpmath无pivot QR在k3稀疏排序上遇到精确零pivot；改小矩阵高精度正规方程，保留完整annihilation residual与85/115位比较。未生成失败版validation，不覆盖旧阶段输出。编译先修复转义产生的form-feed，再修一处文件名溢出，最终六页完成核版。

下一实质方向：回到单singlet的独立导数双线性族，先分开Hermitian分量、守恒投影、Liouvillian后代，再算新的独立Gram/谱矩阵。理由是本轮两个open X描述两个既有singlet的乘积，主要约束相互作用；MSS单介子谱需要不同的独立双线性方向。本轮x不能重命名成单介子内部夸克动量，k也不能重命名成介子编号。继续优先文献及已有chord/embedded-ensemble恒等式，保留同一无序J。
