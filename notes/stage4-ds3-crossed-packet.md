# dS₃真实制备的交叉四点与质量下界见证

2026-10-31。Alan Day / Chengheng，内部研究推导。接续[交叉字典](stage4-ds3-crossing-dictionary.md)及[两物质](stage4-ds3-two-matter.md)。[来源定位](../literature/stage4-ds3-crossed-packet-sources.md)。

## 1. 明确态、算符与所计算的量

半径ℓ=1、Newton G=0、proper SO⁺(3,1)，观察者标量BD表示标签≥s=.2；两个不同自由KG物质species p=.6、q=.9；局部耦合1。ρ(c)=2c²/π，h(c)=√(1+c²)。初态是实非负系数f(c)，

\[
\rho(c)|f(c)|^2=\pi_0(c)=\frac{630}{c_1-c_0}x^4(1-x)^4\mathbf1_{[c_0,c_1]}(c),\quad
x=\frac{c-c_0}{c_1-c_0},\quad(c_0,c_1)=(.4,.8).
\]

这是真实归一初态，f在端点按二次方消失，零延拓为C¹；不称C∞。其质量矩有限。观察者谱在中间仍为[s,∞)，没有额外内部滤波。

使用前轮已证明有界的半周期创建块tilde A。两种次序为v_s=tilde A_q(t)tilde A_p(0)f和v_t=tilde A_p(0)tilde A_q(t)f，时间生成元是真实h而非标签。因物质species不同，两次作用于物质真空都落在两粒子扇区。

令X(t)=〈v_s,v_t〉，N_s=||v_s||²、N_t=||v_t||²，D(t)=N_s+N_t−2ReX(t)=||v_s−v_t||²。它是调节后两次创建的次序差；没有把tilde A说成自伴场，也没有把D称为有限G混沌诊断。

沿用hatR_p(x,y)、g(x)及U_00。完整实谱表达为

\[
X(t)=\int_{c_0}^{c_1}\!dc\,\pi_0(c)\int_s^\infty\!da\,db\,dd\,
\frac{g(a)}{g(c)}
\sqrt{\widehat R_q(a,b)\widehat R_p(b,c)\widehat R_p(a,d)\widehat R_q(d,c)}
U_{00}(b,d;a,q,c,p)e^{it[h(a)+h(c)-h(b)-h(d)]}.
\]

按前轮Hilbert截断极限解释，不预设整个四重积分绝对收敛。正范数N_s=∫π₀(c)g(a)/g(c)·hatR_q(a,b)hatR_p(b,c)，N_t交换p,q。范数与t无关，|X|≤√(N_sN_t)。

## 2. 保留完整谱的解析截断界

数值取a,b,d≤L，c仍覆盖全部初态。k_p(u)=[8π(coshπu+coshπp)]⁻¹，0<hatR_p≤k_p≤e^(−π|u|)/(4π)。记δ(x)=1/(√(1+x²)+x)，

\[
g_{\min}=\frac{c_0e^{\pi\delta(c_1)}}{1+e^{-2\pi c_0}},\quad
C_s=e^{\pi\delta(s)},\quad C_0=\frac{C_s-1}{2\delta(s)}.
\]

有g(c)≥g_min。由exp在[0,δ(s)]上的弦界及xδ(x)≤1/2，g(x)≤x+C₀；若x≥L，另有g(x)≤x exp[πδ(L)]。

每个次序态被删去的平方范数ε_L=||v−v_L||²以如下union界控制（L>c₁）：

\[
\epsilon_L\le T_a+T_b,
\quad T_a=\frac{e^{\pi\delta(L)}e^{-\pi(L-c_1)}}{16\pi^3g_{\min}}
\left[L^2+L(3/\pi-c_0)+3/\pi^2-c_0/\pi\right],
\]
\[
T_b=\frac{e^{-\pi(L-c_1)}}{8\pi^4g_{\min}}(L+C_0+2/\pi).
\]

证明：a>L部分用全实b卷积∫e^(−π|a−b|−π|b−c|)db=e^(−π(a−c))(a−c+1/π)，再以c₀、c₁分别控制多项式和指数。b>L部分先积a，用∫_{a≥s}(a+C₀)e^(−π|a−b|)da≤2(b+C₀)/π+2/π²，再积b。它们允许无限谱，不将L当新物理参数。

令A_p=∫k_p=p/(4πsinhπp)，B=1/(2π³)≥∫|u|k_p(u)du，则

\[
N_s,N_t\le\mathcal N_{\max}
=\frac{(c_1+C_0)A_pA_q+(A_p+A_q)B}{g_{\min}}.
\]

因此有不依赖数值范数的严格质量截断界

\[
|X-X_L|\le2\sqrt{(T_a+T_b)\mathcal N_{\max}}.
\]

当前参数Nmax≈8.42931e−4；L=12的平方范数尾≤1.16353e−16，交叉绝对尾≤6.26346e−10。D的对应谱截断误差≤2(Ta+Tb)+4√((Ta+Tb)Nmax)。这些界不包含盒内求积、有限精度或hypergeometric渐近尾误差。

## 3. 实际球分量的精确odd-convolution字典

局部顶角的000振幅可以直接化为

\[
E_{abp}=C_p\sqrt{\frac{\cosh\pi b}{\cosh\pi a}}H_p(a,b),\quad
C_p=\frac1{4\sqrt\pi}\sqrt{\frac{\sinh\pi p}{p}},\quad
H_p(a,b)=r_p(a-b)-r_p(a+b),\quad
r_p(x)=\frac1{\cosh\pi x+\cosh\pi p}.
\]

它由已固定的BD/KG顶角精确约分，不是重新定义插入。r_p是偶L¹函数；对odd延拓，半轴核H_p是普通卷积。因此

\[
\int_0^\infty H_q(a,b)H_p(b,c)db
=\int_0^\infty H_p(a,b)H_q(b,c)db.
\]

取两种创建态的四个外部球模全0。因为J_b(0000)=1，其零时刻差分量（不含共同f(c)及半周期因子）为

\[
\Delta_{0000}(a,c)
=-C_pC_q\sqrt{\frac{\cosh\pi c}{\cosh\pi a}}
\int_0^s[H_q(a,b)H_p(b,c)-H_p(a,b)H_q(b,c)]db.
\]

保留半周期后乘exp{π[h(a)−h(c)]/2}，仍是正实乘数。由此，这个特定分量的非交换性恰为删掉[0,s)谱段的边界项。**不能从该分量在s→0消失推出完整态差消失**；其他K分量尚未如此约化。

## 4. 非零符号与有限范数见证

令P=coshπp，Q=coshπq，u=coshπa，v=coshπb。D_P(u,v)=u²+v²−1+2Puv+P²，有H_q/H_p=D_P/D_Q。若q>p且PQ≥sinh²πs，则0<b<s时

\[
\frac{\partial}{\partial u}\frac{D_P}{D_Q}
=\frac{2(Q-P)[vu^2+(P+Q)u+v(1+PQ-v^2)]}{D_Q^2}>0.
\]

于是a>c时方括号积分严格正，Δ0000<0。当前条件左减右约28.1211>0。

取运动学有限K测试向量χ，四球模全0，质量系数χ(a,c)=𝟙{1.1≤a≤1.4}f(c)。其group average的物理平方范数为

\[
N_\chi=\int_{1.1}^{1.4}da\int_{.4}^{.8}dc\,|f(c)|^2K_4(a,q,c,p),
\quad K_4=\int_0^\infty\sinh^2r\prod_{n=a,q,c,p}\frac{\sin(nr)}{n\sinh r}dr.
\]

此前四体群平均保证该紧质量/有限K测试范数有限。其与实际调节后次序差的配对为

\[
Z_\chi=\int_{1.1}^{1.4}da\int_{.4}^{.8}dc\,|f(c)|^2
e^{\pi[h(a)-h(c)]/2}\Delta_{0000}(a,c)<0.
\]

物理Cauchy不等式给D(0)≥|Zχ|²/Nχ>0。这是模型内非零性的解析见证，不依赖完整6j积分的近似误差。数值给Zχ≈−1.1370592534e−5、Nχ≈.69678827298，对应下界约1.85552e−10；数值常数仍未区间认证。Nχ同时用独立径向群积分核对。

固定a,c,p,q，H_p(a,b)=2πb sinhπa/(coshπa+P)²+O(b³)，所以该球面差分量在s→0为O(s³)，系数由两核乘积给出。这里只说明该见证的缩放，没有声称全D为s⁶或完整s→0极限为零。

## 5. 完整6j的加速与精度边界

原LPRS每个₄F₃的参数excess严格为1，因此其第n项可写

\[
t_n=Cn^{-2}\exp\left(\sum_{k\ge1}\ell_kn^{-k}\right),\quad
\ell_k=\frac{(-1)^{k+1}}{k(k+1)}
\left[\sum_i B_{k+1}(a_i)-\sum_j B_{k+1}(b_j)-B_{k+1}(1)\right],
\]

其中C=∏Γ(b_j)/∏Γ(a_i)，B_k为Bernoulli多项式。令c₀=1，c_k=k⁻¹Σ_{j=1}^k jℓ_j c_{k−j}，用CΣ_{k=0}^K c_k ζ(k+2,N)近似n≥N尾，并显式累加前N项。使用complex128批处理，保留完整两Omega项、两shadow项、相位及相消指标。常规N=128,K=12，精度变化N=256,K=16。

这是已知Gamma比渐近与Hurwitz zeta求和的应用，没有证明复参数一致的余项界。45位参考检查含6个高波包权重敏感点、1个最坏相消点和6个固定随机点。最坏相消点相对误差约1.97%，其余高权重点最高约1.96e−11。不能把少数好点泛化到全谱，也不把condition·机器epsilon当认证误差。

## 6. 数值保存与当前解释

完整复kernel网格以NPZ保存；JSON记录初态、模型、各次求积节点、软件版本、源码和网格哈希、真实复X、绝对贡献和、相消指标及敏感节点。五组完整网格共3,442,832个复kernel点，[汇总](../results/stage4/ds3-crossed-packet/summary.json)保存独立E/W/ρ组装重放、完整复差及Cauchy比；最终56³×12=2,107,392点。最高盒L=12、每段8节点，初态12节点。正范数用L=20、每段32节点、初态16节点另算。

N_s≈6.053283797898185e−5，N_t≈6.0662442257299775e−5。

| t | Re X | Im X | D(t) | Cauchy比 |
|---|---:|---:|---:|---:|
|0|6.02722781982e−5|−6.87e−19|6.50723839841e−7|.99463135|
|.5|5.88726862884e−5|−5.61586481373e−6|3.44990765951e−6|.97594497|
|1|5.46710731852e−5|−1.11360697956e−5|1.18531338658e−5|.92072470|
|2|3.85245427342e−5|−1.98191785175e−5|4.41461947678e−5|.71494027|

这些小数位是保存值，不代表精度认证。L12的6→8节点比较同时改变初态10→12节点及series设置，最大完整复差1.01718e−11；L10→12的6节点比较同时改变series设置，最大差8.53e−15。同一pilot网格只改series参数，四个X最大差2.03309e−20。不能把混合加密当成各误差源独立认证。

[顶角、谱尾和见证检查](../results/stage4/ds3-crossed-packet/bounds-validation.json)24项通过，见证范数独立径向积分T40附解析尾。所有五组完整raw网格用原始E、W、ρ重新组装，重放差低于2e−17验收阈值。13点45位kernel对照保留最坏相消点约1.97%的相对误差；未以整体series微小变化来抹掉这一精度限制。

当前结论分级：谱截断界与非零球面见证是模型内解析推导；完整X(t)数值是实际无限谱模型的受控截断加经验盒内求积，尚非总误差认证；有限G、全Fock、同species Bose以及完整de Sitter全息仍开放。
