# 横向质量、有效 Bessel 参数与中心局域化的限制

2026-10-30；Alan Day、Chengheng。作者自查证明草稿，未作独立认证或优先权声明。[LaTeX](../manuscript/recoil-transverse/main.tex)、[PDF](../output/pdf/de-sitter-recoil-transverse.pdf)、[来源](../literature/stage4-recoil-transverse-sources.md)、[计算](../code/recoil_transverse.py)、[诊断](../results/stage4/recoil-transverse/diagnostics.json)、[后果计算](../results/stage4/recoil-transverse/consequences.json)。

## 1. 对象、精确核与本轮结论

对象仍为固定质量先 κ→∞、半径1的内部手征 Barnes 图；随后做 ν、M 大质量极限。沿用上一轮的内部轮廓解析继续，ω=iα，α≥0，u、v实。没有新给出物理态、代数、外部插入或断开归一化；不是有限G四点函数。

质量为 m₁,₃=M+ω/2±u、m₂,₄=M−ω/2±v，边次序为 (m₁,ν,m₂,m₃,ν,m₄)。CLP(3.46)四个三元集合给 a_j=M+ν/2+(±u±v)/2；三个四元集合给 M+ν±ω/2、2M。注意 CLP(3.47)本身是 V 的定义，不是另一个质量集合公式。

令 K=ν²/(2M)，ρ=K³/ν²，固定 0≤θ<1、Y<∞，且 |u|+|v|≤θν、α≤YK。定义
\[
w_{\epsilon\eta}=\frac{\nu-\epsilon u-\eta v}{2},\quad
k=\frac{2\sqrt{\prod w_j}}M
=\frac{\sqrt{[\nu^2-(u+v)^2][\nu^2-(u-v)^2]}}{2M},
\]
\[
\lambda=\frac18\sum_jw_j^{-1}
=\frac14\sum_{\epsilon,\eta=\pm1}(\nu-\epsilon u-\eta v)^{-1}.
\]
于是 (1−θ²)K≤k≤K，w_j≥(1−θ)ν/2，λ≤1/[(1−θ)ν]。对称点回到 k=K、λ=1/ν。

b=1/2−i(M−ν)、c=1/2−i(M+ν)，归一因子和核为
\[
C(u,v)=e^{-\pi(M+\nu)}\prod_j\Gamma(-iw_j)\frac{\Gamma(b)}{\Gamma(c)},
\]
\[
Q(t)=e^{-i\pi t}\prod_j\frac{\Gamma(t-iw_j)}{\Gamma(-iw_j)}
\frac{\Gamma(b-t)\Gamma(c)}{\Gamma(b)\Gamma(c+t)}
\Gamma\left(\frac{1+\alpha}{2}-t\right)
\Gamma\left(\frac{1-\alpha}{2}-t\right).
\]
I=(2πi)⁻¹∫Q=𝔅/C。轮廓从−N−i∞向上到−N+i，再向右。κ极限中的顶角已经在原分解中处理，这里没有再往完整6j上乘一次顶角。

**联合估计。** 对固定 θ、Y，K充分大、ρ充分小时，存在只依赖这两个固定参数的正数 C₀,c₀，使
\[
|I-e^{i\lambda S_k}J(k,\alpha)|
\le C_0(\rho+e^{-c_0K})E(k,\alpha),\qquad S_k=k^2+\alpha^2-1,
\]
\[
J(k,\alpha)=-e^{i\pi\alpha}kK_\alpha(k)-i\pi kI_\alpha(k),
\qquad E=kK_\alpha(k)+\pi kI_\alpha(k).
\]
这是上一轮对称内部图估计的横向扩展，**常数不对 θ→1 一致**。包含首项相消处的包络绝对误差；只在 |J|≥εE 的区域转成相对估计。未认证具体最小质量门槛。

## 2. 保留横向的乘子与精确递推

沿用 G、V 的定义，将其中 K 换为 k；令
\[
R(t)=e^{i\pi t}(k/2)^{-2t}
\prod_j\frac{\Gamma(t-iw_j)}{\Gamma(-iw_j)}
\frac{\Gamma(b-t)\Gamma(c)}{\Gamma(b)\Gamma(c+t)}.
\]
有 Q=R exp(−2πit)G 和完全相同的精确分支身份
\[
e^{-2\pi it}G=-e^{i\pi\alpha}G+2\pi e^{i\pi\alpha/2}V.
\]
在 |t|≤CK 上，从移位0连续选对数，Gamma比展开给
\[
\log R=4i\lambda(t^2-t)+O_{C,\theta}(\rho).
\]
具体地，四个大Gamma的首幂是 (∏w_j)^t exp(−2πit)，二次项是 i(t²−t)Σ1/w_j /2。余项由 Σ(|t|+|t|³)/w_j²≤CθK³/ν² 控制。b、c首幂相对 M⁻²ᵗ exp(iπt) 的对数偏差至多 C|t|(ν²/M²+1/M)，其移位余项至多 C(|t|+|t|²)/M，全部为 O(ρ)。因此不需要将 λK² 当成小量；矩形内只用 |R|≤exp[Oθ(ρ+K²/ν)]=exp[o(K)]。

Gamma递推还给精确身份
\[
\frac{Q(t+1)}{Q(t)}=
\frac{\prod_j(t-iw_j)}{
[(t+\frac12-i\nu)^2+M^2][(t+\frac12)^2-\alpha^2/4]}.
\]
这两项从原始 Barnes 核逐点核对，归一化中的相位没有通过重新取模来替代。

## 3. 尾部常数必须重新选择

**正尾。** 取 η=1、d₀=max(α²/4−1,0)、A=max(2,8K,2d₀)、L₀=M−ν+1、a=max_j(w_j−1)。充分大ν时所有 w_j>1，a≤ν。若 ν≥16K，则 L₀≥ν²/(4K)，且
\[
(1+1/(2A))^2(2a^2+a^4/A^2)
\le \frac{75}{128}\frac{\nu^4}{16K^2}\le L_0^2.
\]
证明将左式除以 ν⁴/(16K²)，用 A≥8K、ν≥16K 得 (25/16)(32K²/ν²+16K²/A²)≤75/128。旧的 A≥2K 不足以直接照搬；这里已增加到8K。

令 q=x+1/2。两分母下界仍是 q²+L₀²、q²−d₀，四个分子模之积≤(x²+a²)²。与旧代数链相同，对 x≥A 有 |Q(x+1+i)/Q(x+i)|≤[x/(x+1)]^(3/2)。

令 D=max(8,Y²/2)、ν²≥32(D+1)K³、R₀≥(32+2Y)K。在 R₀≤x≤A+1 的过渡区，
\[
|Q(t+1)/Q(t)|\le2[(A+1)/L_0+a^2/(R_0L_0)]^2\le1/8.
\]
两项分别≤1/8；如果该区间为空无需使用。格点递推和远尾相接给 ∫_{R₀}∞|Q(x+i)|dx≤10(R₀+1) sup_{[R₀,R₀+1]}|Q|。故正尾仍从O(K)基段起控。

**负尾。** 取 L=√(k²+α²)/2、N=ceil L、d=1/2−L。K充分大保证 N>α/2+1/2。对每个 w>0，精确整数Gamma模长式给
\[
\left|\frac{\Gamma(-N-i(w+s))}{\Gamma(-N-iw)}\right|
\le e^{-\pi s/2},\qquad s\ge0.
\]
其余因子使用与上一轮相同的半整数乘积及小Gamma正实部模长界，得到
\[
|Q(-N-is)|\le|Q(-N)|e^{-\Lambda s},\quad
\Lambda=2\pi-N[(M-\nu)^{-1}+(M+\nu)^{-1}].
\]
这是四个不同w分别成立的界，不依赖极点合并。−N实数上的领先相位模长1，故基点由 k 的G鞍点控制，|Q(−N)|≤Cθ,Y E/√K。取s=HK后尾部相对于E指数小。

## 4. 有限变形、两支矩与联合余项

选r在 [(32+2Y)K,(32+2Y)K+1] 内，距两组实极点格点至少1/8。取 T=Hθ,Y K；充分大K、小ρ保证 T<min w_j 且 T<M−ν。

G支改走有限竖线 d−iT 到 d+iT；V支在下面顺时针闭合，只收 t_n=(1+α)/2+n<r，权重
\[
\omega_n=\frac{(k/2)^{1+\alpha+2n}}{n!\Gamma(1+\alpha+n)}>0.
\]
不得将真实R沿竖线无限延伸穿过min w_j，也不得闭入 b−t 的远端极点。上下边的Stirling指数为 −πHK+Oθ,Y(K[1+log(1+H)])，H充分大时压过 E≥cθ√K exp(−K)。右边经反射、极点间隙和r≫α,k给指数小界；V在闭合底边及右边Im t≤1处满足 |V|≤C|G|。乘子exp[o(K)]不改变指数小性质。正尾基段用同样反射估计，再接第3节；负尾用第3节的精确界。

故
\[
I=-\frac{e^{i\pi\alpha}}{2\pi}\int_{-T}^{T}R(d+is)G(d+is)ds
-2\pi i\sum_{t_n<r}\omega_nR(t_n)+O(e^{-cK}E).
\]
这是有限竖线、有限留数和与全部连接边的分解。

令 A(t)=4(t²−t)，概率 p_n=2ω_n/[kIα(k)]。固定α对k求导的Bessel方程给 E_p(A−S_k)=0、E_p(A−S_k)²=8k² E_p t_n=Oθ,Y(K³)。实相位Taylor余项乘 λ² 后为 Oθ,Y(ρ)。

负支令 q(s)=A(d+is)−S_k=−4s²−8iLs。Gamma模长积给有限线上的 exp(−c s²/K) 控制、全辅助线的 ∫|G|≤C kKα(k)，及复一阶矩 ∫qG=0。Taylor余项含 λ²|q|²exp(8λL|s|)；完成平方的代价是 exp[O(λ²L²K)]=exp[Oθ,Y(ρ)]。于是两支各自的主误差为 O(ρ) 倍各自Bessel模长；辅助多项式尾和截断权重尾指数小。相加证明第1节结论。旧相消条件数问题没有通过除以J规避。

## 5. 自然横向尺度：两支的曲率相反

设 u=√M ξ、v=√M ζ，r²=ξ²+ζ²，α=Ky。固定ξ、ζ、y时，
\[
k^2=K^2-Kr^2+(\xi^2-\zeta^2)^2/4,\qquad
k=K-r^2/2-\xi^2\zeta^2/(2K)+O(K^{-2}).
\]
固定α改变Bessel自变量，记 Φ(x,α)=√(x²+α²)−α asinh(α/x)，则 ∂ₓΦ=√(1+α²/x²)。标准正实轴Bessel鞍点展开因此给
\[
\frac{kI_\alpha(k)}{KI_\alpha(K)}\longrightarrow
e^{-\sqrt{1+y^2}r^2/2},\qquad
\frac{kK_\alpha(k)}{KK_\alpha(K)}\longrightarrow
e^{+\sqrt{1+y^2}r^2/2}.
\]
对固定有界y及 r≤B、B≥1、B=o(K^(1/4))，两式对数余项为 O_Y(B⁴/K)。此处对α求导会得到错误的曲率，不能沿α=yx同时求导。y=0可直接用大自变量展开；有界正y用均匀鞍点，连接到y=0不出现奇点。

这只直接描述Bessel支模长。有限解析支另有 Oθ,Y(ρ) 误差；完整和在相消点不保证相对比。对于固定有界r，λS_k−(K²+α²−1)/ν=O_Y(Kr²/ν+r⁴/ν)→0，但C的相位仍须精确保留。

Gamma纯虚模长给精确式
\[
\left|\frac{C(u,v)}{C(0,0)}\right|
=\frac Kk\frac{(1-e^{-\pi\nu})^2}{\sqrt{\prod_j(1-e^{-2\pi w_j})}}.
\]
自然尺度时它趋1，固定θ盒内等于(K/k)[1+Oθ(exp(−π(1−θ)ν))]，不会供应抵消增长支的横向Gaussian。但其复相位通常不小，不能由模长趋1写C比→1。复质量的两手征因子也不能擅自替换为新模方。

## 6. 一个可证的非局域化后果

令 F(y)=√(1+y²)−y asinh y，正根y_c≈1.50887956。固定y>y_c以及p=u/ν、h=v/ν，0<|p|+|h|≤θ<1，定义 q=k/K<1、φ_y(q)=√(q²+y²)−y asinh(y/q)。因为φ'_y(q)=√(q²+y²)/q>0，Δ=φ_y(1)−φ_y(q)>0，且φ_y(q)<F(y)<0。两处都是K支指数主导，包络相对误差可用于真实内部图。由第1节及C精确模长：
\[
\left|\frac{\mathfrak B(u,v)}{\mathfrak B(0,0)}\right|
=\left(\frac{1+y^2}{q^2+y^2}\right)^{1/4}
e^{K\Delta}\{1+O_{\theta,y}(\rho+K^{-1}+e^{-cK})\}.
\]
此为κ先取极限的内部图模长，不是物理OTOC增长。它说明在这些参数上，整个实横向盒不能由中心的衰减Gaussian统一控制；**并不证明物理外积分发散**，也不排除物理轮廓、两手征组合或外腿改变结论。

y=2、p=.3、h=.2时，q≈.861684397、Δ≈.328336461。复用旧对称积分与本轮横向有限分支计算、乘精确C模长，K=8/16/32得到比值约13.8842、192.471、36871.8；渐近式为14.0119、193.756、37048.6。它们支持该指数，但未认证有限K的全连接边常数。

## 7. 数值证据与下一步

主程序67项诊断通过：3个有效参数乘积、18个原核归一/递推/分支身份、3个C模长身份、16个尾界样本、4个有限轮廓/细化、12个联合加阶、8个精度/高度扩展、3个自然尺度速率比较。后果程序另16个乘子残差增精度检查；K=16至1024的残差/ρ记录为样本，不包装成一致界认证。两程序共83项。

联合路径ν=8K^(7/4)，K=8/16/32、y=0/2，横向取(√M,√M/2)或(.3ν,.2ν)。使用20/32阶及55/75位，最大K另95位与高度8K→10K。最大加阶包络变化1.7543e−18，最大增精度3.5523e−70，最大高度变化2.3381e−45；它们是经验误差，不是总误差证书。负/正支条件数最大约1.0804/1.0263。有限轮廓带全部连接边的最大身份误差1.011e−44，网格.2，未放宽阈值。原始复数、参数、版本和哈希均保存。

下一步应先重读完整物理外部循环、两手征因子及N归一，判断当前虚ω延拓上的K支如何进入物理量；再决定是否需要θ→1的极点逼近层或更远横向域。不能先把增长支扔掉，再用单鞍点Gaussian宣称已完成外积分。此前单支局部Hessian结论保留，适用范围不升级。有限κ、有限G四点、实际态和完整de Sitter对偶仍开放。
