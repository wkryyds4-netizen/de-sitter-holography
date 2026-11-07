# 全一物质弦模块的正交框架与移动极限

2026-11-08。本地阶段稿：[七页 LaTeX](../manuscript/susskind-chord-module-frame/main.tex)、[PDF](../output/pdf/de-sitter-susskind-chord-module-frame.pdf)。沿用 Lin–Stanford 表示；正交化与标准离散系列不作为新的表示论发现。当前推进是把实际弦生成元接到明确的移动框架，建立固定有限词的范数极限。

## 设定与精确映射

先 fixed-lambda 半填充双标度，再 q=exp(-epsilon)→1，epsilon=4lambda；K=H/(2Js)=h/sqrt(epsilon)，hbar=1。tracial one-insertion module 权 d>0、r=q^d，实际 charged 曲率通道 d=2。

记 [n]=(1-q^n)/(1-q)，(r²;q)_k=prod_(j<k)(1-r²q^j)。左湮灭算符的真空线可直接正交化：

\[
\xi_k=\frac{\prod_{j<k}(a_R^\dagger-rq^j a_L^\dagger)\Omega_d}
 {\sqrt{[k]!(r^2;q)_k}},\qquad
e_{n,k}=\frac{(a_L^\dagger)^n}{\sqrt{[n]!}}\xi_k.
\]

证明见稿件 Proposition 1：aL xi_k=0，同一总级不同 n 有不同 aL†aL 特征值 [n]；每级维数给完备性；新真空范数递推 s_k²=1+q s_(k-1)²-r²q^(2k)=[k+1](1-r²q^k)。因此 F:|n,k>→e_(n,k) 为真正乘积内积的幺正映射。

\[
a_L=a_q\otimes1,\quad a_R=a_q\otimes q^{d+K_0}+1\otimes b_q,
\quad b_q|k\rangle=\sqrt{[k](1-r^2q^{k-1})}|k-1\rangle.
\]

旧端点的 f_n=e_(n-1,1)。固定左能量后，右 Hamiltonian 是参数 r exp(±i theta_L) 的 Al-Salam–Chihara Jacobi 算符，与既有 Mehler matter measure 相容。Xu/Okuyama 原线性坐标的内积含 matter 权；不能将其未经说明视为普通乘积内积。

## 固定有限词，而非只有种子态

取普通振子位移 D(-sqrt(m))，Ucal=F(D⊗1)，Dcore=finite span{|j,k>}。定义第二因子上的

\[
L_-|k\rangle=\sqrt{k(2d+k-1)}|k-1\rangle,\quad
E_0=d+K_0,\quad B_0=(L_++L_-)/2,\quad P_0=i(L_+-L_-)/2.
\]

在 m→∞、epsilon m→0 下，对每个固定非交换多项式 p 和 v∈Dcore，

\[
\|p(E,\widehat B,\widehat P,B_{tr})\mathcal Uv
-\mathcal U p(E_0,B_0,P_0,-B_0)v\|\to0.
\]

这里 Btr=(K_L-K_R)/2，负号不能丢。真热端点态与该相干态的联系仍用上一阶段更强的 epsilon m²→0；本定理本身是移动相干后代核心，不擅自把所有后代换成真热态。

证明路径：位移后 N/m 的均方偏差为 (2j+1)/m+j²/m²；系数是固定次数 Charlier 乘 Poisson 振幅。单边移位在这些向量上趋 -I。n∈[m/2,2m] 的算符系数一致收敛；固定词只作有限步移位。精确系数包络为 C(1+epsilon(n+L))^L exp(L epsilon(n+k+d+L))，平方包络的尾由指数倾斜 Poisson 控制。因此直接控制有限路径展开，不从单算符强收敛非法推乘积收敛。

各有限 q 的 E/Bhat/Phat/Q 是有限总级 Hermitian 块的自伴直和，Btr 有界；极限 Jacobi 满足 Carleman 条件。核心收敛给各共轭自伴算符的 strong-resolvent 收敛，以及固定时间幺正群的强收敛。未声称无界算符的全空间算符范数收敛。

## 正的零方向生成元候选

令 Q0=E0-P0、T0=-B0。则 [T0,Q0]=iQ0，且

\[
\langle v,Q_0v\rangle=\tfrac12\sum_{k\ge0}
|\sqrt{2d+k}v_k+i\sqrt{k+1}v_{k+1}|^2\ge0.
\]

零模递推不在 ell² 中。相位 v_k=i^k w_k 把 Q0 化成标准 Laguerre Jacobi；种子谱为 Gamma(shape=2d,rate=2)，均值 d、方差 d/2。该谱及表示是标准结果；微观弦生成元的移动逼近才是本阶段的连接工作。

有限 q 的 Q=E-Phat 不保证正。极限满足 exp(i eta T0)Q0 exp(-i eta T0)=exp(-eta)Q0。E0+P0 与 Q0 不对易，不能直接当另一个 Poincare 光锥动量；M²/(mu Q0) 仍把质量 M 手动放入。d=2 不决定体质量。

## 完整选态的指定插入协议

保留已知全代数 rho=a+s(h-W)²，s=sqrt(q)，a=1-s(1+lambda)>0，W=sqrt(lambda)U，||W||≤sqrt(epsilon/(1-exp(-epsilon)))。将 one-matter 模块嵌入全 tracial 空间，对单位 zeta=Ucal v 定义 zeta_rho=L_sqrt(rho) zeta/||L_sqrt(rho) zeta||。这是左向量插入协议，不自动等同原 tau(rho X† alpha_t(X)) 的排序。

\[
\eta=\frac{\|(h+2\sqrt m)\zeta\|+\|W\|}{2\sqrt m}
+\frac{\sqrt a}{2\sqrt{sm}},\quad
\|\zeta_\rho-\zeta\|\le2\eta\to0.
\]

谱演算给 ||(|h-W|-2sqrt(m))zeta||≤||(h-W+2sqrt(m))zeta||；sqrt(a+sA²) 与 sqrt(s)|A| 的算符距离≤sqrt(a)。普通位移及 q 振子误差 O_j(epsilon m^(3/2)) 给 h 浓集。对有界 A 的期望差≤4eta||A||，插入范数平方/(4sm)→1。这保留了 U，却只说明领先阶制备信息在该协议里消失；原 selected 关联的排序与次阶 U 仍需计算。

## 证据与数值范围

- frame-checks.json：1,224 Gram、864 ladder、21 Laguerre 等式，共 2,109；q=1/4,4/9，d=1,2,3，总级0..7。
- 48个移动参数组：m=8,32,128,512，epsilon=m^-3，d=2，j=0..2、k=0..3；每组9个词、3个对易子、Casimir及boost匹配。两截断相差48，最大词向量差1.075e-23；输入范数平方尾≤3.178e-48，不把输入尾当全部无界词误差证书。
- m512,j2,k3 的8个词以独立50位稀疏系数计算，最大相对向量差1.395e-13。BP相消max733.45；EP和BE分别可到4.027e8、1.518e9，后两残差未单独高精度认证。解析证明承担收敛主张，末位数不作区间证书。
- null-checks.json：3个d的精确Q正性形式与8个TQ-iQ移动残差；null相消max376.62。m512的残差为0.022581（j0k0）和0.172697（j2k3）。
- m512的归一选态向量距离**上界**为0.088390（j0）及0.143017（j2），不是实际计算出的向量差。
- 输入、输出、版本及哈希见两份JSON与local-evidence.json。成功科学driver/helper已封存，不重跑或改写。

下一步先精读多物质弦 coproduct，确定正Q在实际双物质模块的组合律及另一零方向的微观来源；同时保留实际selected排序。质量、四bilocal和finite-G尚未完成。本轮本地保存，既有Git/issue请求已完成。
