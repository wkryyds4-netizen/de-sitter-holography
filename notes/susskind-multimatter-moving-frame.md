# 共同背景下的两物质移动框架

2026-11-08。本地阶段：[六页PDF](../output/pdf/de-sitter-susskind-multimatter-moving-frame.pdf)、[LaTeX证明](../manuscript/susskind-multimatter-moving-frame/main.tex)。先fixed-lambda半填充DS，再q=exp(-epsilon)、epsilon=4lambda→0；K=H/(2Js)=h/sqrtepsilon，hbar=1。对象是权d1,d2>0的**ordered direct**两物质模块，物理charged例子d1=d2=2。

旧two-matter阶段已经做过primary分解和Jacobi/Laguerre耦合，本轮不重复称新结果。新读HVX2511.08743v2补齐此前遗漏的量子群factorization前驱；该文q_HVX²=q。本轮工作是把这类明确因子与实际LS大共同弦数生成元连接。

## 精确三角框架

在普通ell²_n⊗ell²_k1⊗ell²_k2上，a|n>=sqrt([n])|n-1>，R_i=q^(d_i+N_i)，b_i|k>=sqrt([k](1-q^(2d_i+k-1)))|k-1>。定义

\[
A_0=a,\quad A_1=aR_1+b_1,\quad
A_2=A_1R_2+b_2=aR_1R_2+b_1R_2+b_2.
\]

各A_i为q振子，creators互相对易。raw chord |l,m,r>对应(A0†)^l(A1†)^m(A2†)^r|000>，外aL=A0、aR=A2，总nbar=N0+N1+N2+d1+d2。该映射逐grade三角且非零对角，满射；外lowering和puremiddle norm [m]!确定全部direct Gram，与实际chord一致，所以完成后为幺正F。不是仅检验代数就假设内积正确。

在共同因子位移V=F(D(-sqrtm)⊗1⊗1)，取固定有限j,k1,k2核心。m为位移参数，区别于raw中间弦数。令E_i=d_i+N_i、L_i-|k>=sqrt(k(2d_i+k-1))|k-1>、B_i=(L_i++L_i-)/2、P_i=i(L_i+-L_i-)/2。在m→∞、epsilon m→0下，对每个固定非交换多项式，

\[
(E,\widehat B,\widehat P,B_{tr})
\longrightarrow(E_1+E_2,B_1+B_2,P_1+P_2,-B_1-B_2)
\]

的含义是p(actual)Vv−Vp(limit)v范数→0。保留共同oscillator multiplicity。真实heat端点与coherent的关系仍用旧epsilon m²→0；不将本轮全部coherent后代改称真实热后代。

## 无大数相减的作用公式及证明

置w=d1+d2+N1+N2，R=R1R2，C=b1R2+b2，c=q^((N0+w)/2)，Z=a†RC+C†Ra+C†C，则

\[
E=\frac{[w]-(1-R)^2[N_0]/2-Z/2}{c},\quad
\widehat B=\frac{(1-R^2)[N_0]-Z}{2c\sqrt{1-c^2}},
\]
\[
\widehat P=\frac{i(a^\dagger C-C^\dagger a)}{2\sqrt{1-c^2}},\quad
B_{tr}=\frac{(a+a^\dagger)(1-R)-(C+C^\dagger)}{2\sqrt\epsilon}.
\]

bounded internal levels时C/sqrtepsilon→L1-+L2-。共同n∈[m/2,2m]上E跨项和B对角、Btr共同项均O(sqrt(epsilon m))，B的C†C项O(sqrt(epsilon/m))。共同移位在D(-sqrtm)|j>上趋-I，确定符号。固定词只改变有限步坐标；与旧证明同型指数系数包络乘固定多项式Poisson尾可积并消失，直接控制词范数。该论证也容许穿插固定internal generators，从而得到下面个体协变。

Q_i=E_i-P_i为正自伴且不同因子强对易；Qsum=Q1+Q2，Tsum=-B1-B2满足[Tsum,Qi]=iQi。实际LS E-Phat在V核心上逼近Qsum。还可定义有限q的正强对易候选 mathfrakQi=F Qi F†，其与Btr的协变缺陷在V核心上趋零，和逼近实际LS Q。F中的普通内部Qi并非有限q LS Q；后者不保证正。候选不是唯一性或单侧局域算符声明。

## 同种交换与质量

相同X的physical two-open Wick内积还加cross Gram。旧已证d2前两个direct unit primaries的physical norm²分别1+q⁴与1-q⁸~8epsilon；相同因子随整个外后代模块传播。共同位移不会取消它。ordered factor算符要接同种物质，必须保留crossed metric与奇异归一化，不能直接把裸ordered空间当Bose空间。

HVX(82)–(83)已有一般channel swapping权重的推导，更新了旧LS有限级猜测的文献状态；尚未在本项目完全核对所有能量态归一化。下一步可直接从这一文献结果及原Gram递推入手，不再重新猜全k公式。

总Q/boost即使已确定，也不选择质量：在Qsum=P、Tsum=i(P∂P+1/2)的谱表示，任何正且强对易的反向协变Pminus都为A/(mu P)，A为relative/multiplicity空间上的正自伴算符。由谱分解后的算符场与logP平移协变直接推出。自由两质量取A=m1²/x+m2²/(1-x)只是一个动力学选择。另一个LS组合S=Esum+Psum满足[Qsum,S]=-2iTsum，不是对易的Poincare动量；这不要求全局de Sitter平移对易。

## 验证与继续位置

checks.json（SHA 4f3621ca842f4fac8db04b276567991204fb3ac65bf825053ee0500fddf084db）保存1,484项rational Gram和280项outer lowering，共1,764精确比较；q=1/4,4/9，(d1,d2)=(2,2),(1,2)，总级0..4。12个移动组为m8/32/128/512、epsilon=m^-3、d1=d2=2和(j,k1,k2)=(0,0,0),(1,1,0),(2,1,2)，每组8词及BP/TQ对易子。

最大截断词差1.527e-24；输入norm² Poisson尾3.178e-48，不当全词证书。最大相消BP358.56、TQ528.98。m128,j2,k1=1,k2=2用60位稀疏外ladder直接复算E/B/P/T/Q/BE以及两个对易子残差本身，maxrelative5.923e-13；不是区间证书。成功driver/helper已封存，不重跑/改写。

下一首任务：把HVX已导出的crossing channel权重接入本移动框架，检验epsilon首阶的Casimir依赖，保留physical奇通道的epsilon^-1/2归一化，再算真实同种singlet四点中的影响。尚未声称完成这一首阶计算或确定质量。selected原排序、完整四bilocal、finiteG继续开放。本轮无Git/issue。
