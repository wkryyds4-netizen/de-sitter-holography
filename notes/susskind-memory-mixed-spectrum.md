# 三个记忆场的完整混合两点谱

2026-11-06，本地。[六页PDF](../output/pdf/de-sitter-susskind-memory-mixed-spectrum.pdf)、[LaTeX](../manuscript/susskind-memory-mixed-spectrum/main.tex)、[来源](../literature/susskind-memory-mixed-spectrum-sources.md)。不commit/push/ZIP。

本轮完成上一阶段留下的Y+–R、Y+–F、Q4–F三项。由此，投影后二次memory系数X与一次memory场L的完整双能量cross为零。Y/L/X正谱矩阵闭合；有限memory核仍需恢复投影掉的动态后代和余项。

## 约定和新结论

标准CAR、完整overlapping complex Gaussian p-body，固定半填充m=N/2，annealed归一HS，h=H_m，ℏ1。先λ_N=p²/N→λ>0，q=e^−4λ，C_N→C=Js²/λ，之后才小λ。coupling chaos Wj与spectral q-Wick不混用。

Y=Y−；一次场L_N=mZ1−C1,m=R+aJ_hY+F+oL²，a=−(1−q)²/(1+q²)。二次X=(1−hatPi)mZ2，hatPi为既定八span(I,B,O−,JhF,E0,M0,S0,Z0)，X=Q4+8λ_N C_N q_N²T_N+oL²。上轮已证X–Y所有能量后代cross0。

本轮证明ν_R,Y+=ν_F,Y+=ν_F,Q4=0，结合旧结果给ν_T,R=ν_T,F=0及ν_X,L=0。五场(Y,T,F,R,Q4)两点谱对角，其中Y与T共享归一q²核；不由此推出混合高点独立。

## reverse cubic端点证书

三槽R按lex marks(01,02,12)的numerator是[−q(1+q+q²),D3,−q(1+q+q²)]，D3=(1+q²)(1+q+q²)，unmarked为−2λ(1−q)(1+q³)。Y+ numerator为mark(1+q)、unmarked−2λ(1−q)。

三槽对两槽k=0,1,2桥，Σbinom(5−2k,3−k)(6−2k)=74模板，3!2!/k!每orbit，共876字。reverse delta=−4·平行双桥，normal delta=+4·反向双桥，alpha beta不变。两方向均完整计算，148模板1752字λ/λ²多项式全0；每个minor另以determinant核对，共5256次。外部spectator在同orbit公共，因此任意固定双能量阶被覆盖，非低矩外推。谱正性、紧性和候选compact support给完整ν_R,Y+=0。

## 线性rank filter的上电荷提升

沿用真实算符L²桥接f=F/√VF、fhat=√λ_N Σ[R1(p−a)−1]h_a/√C_N，R1(u)=u(N−2p+u+1)(N−m+p)/[(m−p+1)p²]。f−fhat→0由input factorial tail证明，不只是协方差总变差。

设d=m−p+1、g(n)=(n−p+1)/d，U为参考m归一的p−1-body ensemble。g(Q)U是p-body lift，其rank charge因子identity对任意允许n成立；这轮明确检验m→m+1的93076个covariance tensor entries。原来的m→m−1还不足以单独覆盖reverse vertex。

直接定义M−(U,V)=−Σ[U,c†][V,c]、M+(U,V)=−Σ[U,c][V,c†]。任意neutral U,V，限制m后有精确式

\[
M_+(gU,V)=(1+d^{-1})M_+(U,V)-d^{-1}U\bar DV,
\]
\[
M_+(V,gU)=(1+d^{-1})M_+(V,U)-d^{-1}(\bar DV)U.
\]

normal同式为(1−d^-1)M−加d^-1 UDV或d^-1(DV)U。对未shift p-body V，DV=pV，DbarV=pV−DeltaV，旧half input deficit给DeltaV各固定矩有界。所以reverse附加项为O_L²(p/N)=O1/p，插固定字长同样用Holder控制。

不能把D机械替成Dbar的自然body-count式：

\[
\bar D[f(Q)B]_m=f(m+1)\bar DB+(N-m)[f(m)-f(m+1)]B
\]

会留下O1项；实际交换子identity保留了正确相消。

Q4的normal mark用一般k式

\[
[(kp-D)(gB)]_m=(1-d^{-1})[(kp-1-D)B]_m+(k-1)pB_m/d,
\]

B中一因子p−1，其余p，k4适用。非相邻pair mark由连续子产品total marks有限包含排除得到。四个不同非交换因子、每个lift位置都已实际验证。

每个含一个fhat的允许Wick图只有一条signed edge，等于lifted p−1 covariance减p covariance。BNR generic-length overlap规则、固定mark局部determinant及absolute exponential tail给两图相同极限。W2/W4先删内部pair，lift误差O1/p，故逐图相消。fhat替换误差由||f−fhat||₂乘已有B自谱高能量矩控制。正matrix spectral CS再从所有矩得到ν_F,Y+=ν_F,Q4=0。

## 正3×3谱与指定投影多项式

e=(E+E′)/2，W=4Js4[1/λ+4q/(1+q)]，VF=64Js6q²；VR沿旧memory-fluctuations，V4/V2沿上一阶段。完整谱matrix为

\[
d\nu_{(Y,L,X)}=
\begin{pmatrix}
W\mu_{q^2}&aeW\mu_{q^2}&0\\
aeW\mu_{q^2}&V_F\mu_q+V_R\mu_{q^3}+a^2e^2W\mu_{q^2}&0\\
0&0&V_4\mu_{q^4}+V_2\mu_{q^2}
\end{pmatrix}.
\]

相对共同product marginal密度Kj，factorization Mdiag(WK2,VFK1,VRK3,V4K4,V2K2)Mᵀ，M三行(1,0,0,0,0)、(ae,1,1,0,0)、(0,0,0,1,1)，行列式WK2(VFK1+VRK3)(V4K4+V2K2)>0。

明确过滤的Phi(s)=Y+isL−s²X/2有

\[
d\nu_{\Phi(s),\Phi(u)}=W[1+iae(u-s)+a^2e^2su]\mu_{q^2}
+su(V_F\mu_q+V_R\mu_{q^3})+s^2u^2\nu_X/4.
\]

ordinary中心时间反射(E,E′)→(−E′,−E)消去奇e项，得WG_q2(T)+su C_L(T)+s²u² C_X(T)/4。X额外投影了非守恒后代，Phi不是raw memory Taylor polynomial，也没有由此得到有限s,u余项。

## 实际证据与修正

3JSON：cubic-force-cross(2×74模板、2×876字、5256minor)；charged-lift(93076上charge协方差entries、81直接CAR identities、80非零，另保留旧15normal controls)；matrix-assembly(正factorization、行列式、9个多项式系数、共轭与反射)。依赖版本、种子、输入输出hash保存。

保留失败：首checker要求每个witness非零，N6p2m3一traceless one-body U使辅助Dbar(gU)=0；其清分母式−3tr(J)I本来应零。现明确记录index4为零，其余reverse/normal/quartic仍要求非零，不放宽精确等式。LaTeX首编译转义使frac变form-feed，修正后又发现0.57117pt行溢出，改display；旧源码和日志保留。最终六页双遍编译/逐页95dpi核版。

下一步先恢复raw二次记忆系数被投影掉的动态后代，求其与(Y,L,X)的谱接口；已有raw五字W2表无需重算。随后扩中央charge窗并控制有限memory两参数余项，再接体介子空间字典。完整研究goal仍active，当前只本地保存。
