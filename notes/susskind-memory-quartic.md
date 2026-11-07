# 二次记忆系数与第四阶动态singlet

2026-11-06，本地。[六页PDF](../output/pdf/de-sitter-susskind-memory-quartic.pdf)、[LaTeX](../manuscript/susskind-memory-quartic/main.tex)、[来源](../literature/susskind-memory-quartic-sources.md)。无commit/push/ZIP。

## 模型与本轮推进

标准CAR、完整overlap complex Gaussian p-body H=Hs，w=4^p(Js²N/p²)(p!/N^p)²，β=μ=0 annealed归一迹、ℏ1。先N,p→∞且p²/N→λ>0，q=e^(−4λ)、σ²=Js²/λ，后小λ，string时间固定。

本轮得到二次居中memory Z2的有限N精确式，分离首返/漂移的守恒贡献与消失余项；在其第四coupling chaos中确定一个独立动态Q4。807模板全阶相消证明完整两点为已有r=q⁴ Mehler核，另128/354模板证明与Y/R的完整双能量cross为0。原始Z2仍含第二chaos，尚未全部求完，不能把raw Z2等同Q4。

## 精确二次系数

固定m，h=H_m，Vψ=Σ|i〉c_iψ/√m，K=1⊗H_(m−1)，r=p/m，Jh={h,·}/2，Lh=[h,·]。B_n=V†K^nV、b_n=B_n−Σ_(j=1)^(n−1)b_j B_(n−j)。

b4=B4−B1B3−B2²−B3B1+B1²B2+B1B2B1+B2B1²−B1⁴，Z2=b4−2Jh b3+Jh²G，Mhat(s)=G+isZ1−s²Z2/2+O(s³)。此为标准ordered first-return/interval-partition结构，不假设原H或charge blocks具有Boolean独立性。

B^f_ab=Σ(L^a c_i)†L^b c_i，L B^f_ab=−B^f_(a+1,b)+B^f_(a,b+1)，A=B11、M_n=B0n、T_n=(M_n+M_n†)/2。于是

\[
M_4=B_{31}^{f}-3B_{22}^{f}+3B_{13}^{f},\quad
T_4=B_{22}^{f}+2L^2A,\quad M_4-T_4=2LT_3.
\]

中心bilocal第四系数U4=T4−(3/2)L²A=B22^f+(1/2)L²A。代入renewal精确得到

\[
mZ_2=B_{22,m}^{f}+\frac{(1+r)^2}{4}L_h^2A_m
+2rJ_hT_{3,m}+3r^2J_h^2A_m-\frac{A_m^2}{m}-mr^4h^4.
\]

定义旧守恒KA=a_m+c_m(h²−C_m)，Y=A−KA，U=T3−P_m。新守恒多项式与剩余core：

\[
D_{2,m}=2rJ_hP_m+3r^2J_h^2K_A-K_A^2/m-mr^4h^4,
\qquad \mathcal E=B_{22}^{f}+L^2A/4.
\]

精确余项

\[
mZ_2-D_{2,m}-\mathcal E_m
=(r/2+r^2/4)L_h^2Y+2rJ_hU+3r^2J_h^2Y
-[\{K_A,Y\}+Y^2]/m.
\]

中央窗r=O(1/p)、a_m=O(p)，旧固定阶矩界给右侧O_L²(1/p)；尾由binomial指数界压过多项式。特别是不能先删A²/m：其O(1)守恒部分已经放入D2。该分离还没有扣掉E自身的守恒和第二chaos部分。W4D2→0，且Gaussian chaos投影为L²收缩，所以第四chaos的有效core确实为W4E。

## 第四chaos的五维空间

记D_op X=QfX−Σc†Xc、Tcal_k=(kp−D_op)H^k。五个raw Hermitian表达式：

- E0=(AH²+H²A)/2；M0=HAH；
- S0=(H Tcal3+Tcal3 H)/2−E0−M0；
- C0=Tcal4−2E0−M0−2S0；Z0=H⁴。

对全部取W4后记E,M,S,C,Z。在四槽中它们对应marks：(01+23)/2、12、(02+13)/2、03、无mark。此处E/M等与旧memory的记号不同，均限本节basis。

W4 E_core=C−2S+(E+M)/2；W4(J_HT3)=S−(E+M)/2，W4(J_H²A)=(E+M)/2，W4(L²A)=2(E−M)。故E,M,S,Z正是旧场的第四chaos后代方向，不全是守恒量。

Gram矩阵由两个四槽W4块间的24配对求得。定义d3=q²−q+1、d5=1+q+q²+q³+q⁴，扣除四个后代得到

\[
Q_4=C-\alpha E-\beta M-\gamma S-\delta Z,
\]

\[
\alpha=-2q^3/(d_3d_5),\quad
\beta=-q^2(1+q+q^2)/(d_3d_5),\quad\gamma=2q/d_3,
\]

\[
\delta=\frac{2\lambda(1-q)^3(1+q)(1+q+q^2)}{(1+q^2)d_5}.
\]

Gram正交是DS结论；定义在有限N用同样λ依赖系数，不声称有限N精确正交。W4E_core中Q4系数为1，其余是(α+1/2)E+(β+1/2)M+(γ−2)S+δZ。

## 正范数、全阶核与交叉

P(q)=1+q+2q²+q³+q⁴，T(q)=2+q+2q²+q³+2q⁴，

\[
V_4=\sigma^8\frac{4\lambda(1-q)^2(1+q)^2[P(q)+4\lambda q^3T(q)]}{(1+q^2)d_3d_5}>0.
\]

独立加两个H的216+216配对给μ2=2σ²(1−q⁴)，但全核不靠此低矩外推。采用旧single-mark determinant局部规则：每图q^cr[λδ+λ²αβ]，α=2−4×incident-cross，δ=4仅两标记反向双桥连接。

两个四槽块间k=0..4条桥，未配完的8−2k端点以A/B色及两弧split分类。模板总数Σ binom(8−2k,4−k)(9−2k)=807；每orbit24²/k!，共452760字。消分母后的整数q多项式证书显示k<4的λ与λ²系数全部恒0，k4幸存项与V4精确一致。外部spectator与连续块内部置换的cross不变，故可提出；这穷尽任意固定两弧矩，不是只检有限阶。

唯一幸存四桥穿一条H弦权q⁴，正双能量测度/高矩紧性/紧候选矩确定性给

\[
\langle Q_4^\dagger Q_4(T)\rangle_{DS}=V_4G_{q^4}(T).
\]

四槽对两槽Y有128模板、对三槽R有354模板，全部λ和λ²系数为0。同一spectator论证与spectral Cauchy–Schwarz给ν_Q4,Y=ν_Q4,R=0，连带所有固定左右能量后代cross为0。未声称混合高点独立，也未完成Q4–F或Q4–Y²。

核是已有Mehler公式。0<q⁴<1使候选绝对连续、没有diagonal atom；Portmanteau给有限N归一exact static权重趋0。后小λ：V4∼(768/5)Js⁸/λ，归一关联→sech⁸(2JsT)，μ2→32Js²。chord Δ4尚不是体介子质量。

## 实际验证与保留的问题

- quartic.json：finite初版15 identities、5×5 Gram及projection，24/216/216配对，324局部minor比较。原finite用弱Gaussian J取整，两个p2系统成零H；这些原值保留，不当充分验证。
- independent-car.json补上真正非零integer Hermitian J：3系统12个精确式，m=p保留为G=0控制；另64 Wick固定字、1600对照，239非零、8负trace，核五个primitive mark与局部determinant的独立有限CAR一致性。
- orbit-certificate.json：807模板、452760配对字、360独立minor核对；cross-kernel.json另482 cross模板、survivor与Gram精确一致、Laurent768/5及4λ曲线。Jacobi/+12解析截断<1e−13，λ.25/.1用128/192-node独立Mehler；舍入非interval认证。
- 独立checker首轮未展开complex表达式，SymPy无法判断正性；保留初版与说明，改为先expand。第二次表达式膨胀运行，确认live后主动Ctrl-C，保存源码并在每个矩阵moment/renewal步立即展开，session14677终态通过。没有以无输出误判结束或启动重复进程。
- cross脚本首次解析保存的lambda变量触发Python关键字SyntaxError；原源码/说明保留，parser token改lam后通过。PDF初版缺mathscr宏，换mathcal；后有3.16pt行溢出，改display公式，最终六页双遍和逐页核版通过，旧log保留。

下一步优先W2 E_core、其H²守恒投影与动态rank-filtered残量，并把W4E/M/S/Z后代恢复成完整算符，补回低chaos收缩。随后求实际Z2完整谱及Q4与F/Y²混合，再控制有限s,u的二插入余项，接体空间字典。完整goal active，不把纯Q4全核当raw Z2或完整全息已完成。
