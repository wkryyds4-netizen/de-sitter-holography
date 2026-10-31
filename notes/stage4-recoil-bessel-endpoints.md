# 联合端点控制与中央段的相消障碍

2026-10-30。Alan Day、Chengheng。作者自查推导；数值为有限点诊断，不是区间认证。英文完整推导见 [LaTeX](../manuscript/recoil-bessel-endpoints/main.tex)，[PDF](../output/pdf/de-sitter-recoil-bessel-endpoints.pdf)。

## 1. 对象与结论边界

沿用固定质量先 κ→∞ 的内部手征 Barnes 图积分。半径 1，u=v=0，ω=iα，α=Ky，0≤y≤Y<∞，M=ν²/(2K)。态、物理外插入、算符代数和四点归一化尚未识别，不能称为有限 G OTOC。

精确 Qν、辅助 U 及 C 的定义沿用[完整固定 K 笔记](stage4-recoil-bessel-limit.md)。轮廓从 −N−i∞ 向上到 −N+i，再向右；积分带 1/(2πi)。辅助全积分为
\[
J=-K[e^{i\pi\alpha}K_\alpha(K)+i\pi I_\alpha(K)],\qquad
\mathcal E=K[K_\alpha(K)+\pi I_\alpha(K)]>0.
\]
负 α 利用完整和的偶性改用 |α|。以下用 E 归一，不除以存在相消零点的 J。

本轮补齐上一轮缺少的正尾基段 H_R/E 及 N=O(K) 负竖直轮廓。**中央段仍未完成**：在固定高度轮廓上，即使除以 E，其绝对积分也可指数增长。不能用两端已控推断全积分余项已控。

## 2. 保留相位的局部乘子

取 K→∞、ρ=K³/ν²≤ρ_Y，ρ_Y 为充分小的固定正数。对 |t|≤C_YK，在取消相同的小 Gamma 因子后，以连续对数定义 Qν/U：
\[
\log(Q_\nu/U)=4i(t^2-t)/\nu+O_Y(\rho).
\]
不是直接套用固定移位的 Gamma 比公式：对大参数的安全扇区使用 Stirling，再对 w=O(K) 作 Taylor，
\[
\log\Gamma(z+w)-\log\Gamma(z)
=w\log z+\frac{w(w-1)}{2z}
+O\big((|w|+|w|^2+|w|^3)/|z|^2\big).
\]
整个线段 z+θw 仍处在统一安全扇区，因为 K/ν→0。四个 ν/2 Gamma 给首相位，大 b,c 的剩余量由
\[
C_Y(K^2|t|+K|t|^2+|t|^3+K)/\nu^2
\]
控制。K²/ν 未必小，不能展开掉首相位。

t=−N 为实数时，首项纯虚，故 |Qν/U|≤exp(C_Yρ)。t=x+i 时首项实部为 −4(2x−1)/ν；正基段因此也有相同上界。中央固定高度条带有双向模长比较 exp[±C_Y(ρ+K/ν)]。

## 3. 优化左基点

令 F(y)=√(1+y²)−y asinh y，L=K√(1+y²)/2，N=ceil L。K 充分大时 N>α/2+1/2，两个小 Gamma 实部至少 c_YK。移动左角点不跨极点：小 Gamma 的极点族在其右侧，分子极点高度为 ν/2>1，b−t 的极点实部为正；下方闭合边由固定参数竖直衰减消失。实 Stirling 给
\[
|U(-N)|=2\pi\exp[-KF(y)+O_Y(K^{-1})].
\]
把主指数视为 n 的函数，
\[
g(n)=(n+\alpha/2)\log[2(n+\alpha/2)/K]
+(n-\alpha/2)\log[2(n-\alpha/2)/K]-2n.
\]
g(L)=−KF、g'(L)=0、g''=O_Y(1/K)，所以取整只付出 O_Y(1/K)。

为覆盖 y→0，E 的统一下界不用大阶数渐近，而用标准正积分
\[
K_\alpha(K)=\int_0^\infty e^{-K\cosh u}\cosh(\alpha u)\,du.
\]
在 u*=asinh y 到 u*+K^(-1/2) 上积分，以 cosh(αu)≥e^(αu)/2 和极大点 Taylor 控制，得到
\[
\mathcal E\ge c_Y\sqrt K e^{-KF(y)},\qquad
|Q_\nu(-N)|/\mathcal E\le C_YK^{-1/2}e^{C_Y\rho}.
\]

## 4. 精确负竖直衰减

写 B=M−ν>0、C₀=M+ν、v=ν/2。Gamma 递推与标准模长身份给
\[
|\Gamma(-N-iw)|^2=
\frac{\pi}{w\sinh(\pi w)\prod_{j=1}^N(j^2+w^2)},\quad w>0,
\]
\[
|\Gamma(N+\tfrac12+iw)|^2=
\frac{\pi\prod_{j=0}^{N-1}[(j+\tfrac12)^2+w^2]}{\cosh(\pi w)}.
\]
负半整数的倒 Gamma 同样由递推获得。相对 s=0 的模长比中：

- 显式指数提供 e^(−πs)；
- 四个负整数实部 Gamma 提供至多 e^(−2πs)，所有有理乘积因子≤1；
- 两个大半整数 Gamma 的双曲部分至多 e^(πs)，因为 log cosh 的导数模≤1；
- 多项式部分至多 (1+s/B)^N(1+s/C₀)^N；
- 两个小 Gamma 的比≤1，由 a>0 时 |Γ(a+is)|≤Γ(a)。

所以对所有 s≥0 有精确不等式
\[
|Q_\nu(-N-is)|\le |Q_\nu(-N)|e^{-\lambda s},\quad
\lambda=2\pi-N(B^{-1}+C_0^{-1}).
\]
在上述层级内 λ≥π 最终统一成立，因而
\[
\mathcal E^{-1}\int_{s_0}^{\infty}|Q_\nu(-N-is)|ds
\le C_YK^{-1/2}e^{C_Y\rho-\pi s_0}.
\]
−N 到 −N+i 的短连接段，用 |Γ(a−iτ)|≤Γ(a)、e^(2πτ)≤e^(2π) 和乘子界，同样得到 O_Y(K^(-1/2))，额外因子 exp[C_Y(ρ+K/ν)]。

**这些估计仅含竖直段，不含中央负水平段。**

## 5. 正尾归一估计

R=c_YK，c_Y=8+2Y。在 x∈[R,R+1]、t=x+i 上反射两个小 Gamma：
\[
U(t)=\frac{\pi^2e^{-2\pi it}(K/2)^{2t}}
{\cos\pi(t-\alpha/2)\cos\pi(t+\alpha/2)
\Gamma(\tfrac12+t-\alpha/2)\Gamma(\tfrac12+t+\alpha/2)}.
\]
cos 因子模有正下界。r=x/K，Stirling 给 |U(x+i)|≤C_Y exp[−Kh_y(r)]，其中
\[
h_y(r)=(r+y/2)\log[2(r+y/2)]
+(r-y/2)\log[2(r-y/2)]-2r.
\]
它对 y 偶且凸，h_y≥h₀=2r(log(2r)−1)。因此
\[
H_R\le C_Ye^{C_Y\rho-Kh_0(c_Y)}.
\]
结合上一轮的联合精确递推引理（ν²≥32[max(2,Y²/2)+1]K³），以及 F≤1，
\[
\mathcal E^{-1}\int_R^\infty|Q_\nu(x+i)|dx
\le C_Y\sqrt K\,e^{-[h_0(c_Y)-1]K+C_Y\rho}.
\]
H_R 是真正 supremum；此解析上界不以数值采样最大值代替。

## 6. 中央段的绝对条件数

固定 y>1，对 t=Kx+i、−y/2<x<y/2，仅反射负实部小 Gamma。任何内紧区间上
\[
|U(Kx+i)|\asymp_y e^{KA_y(x)},\quad
A_y=u(\log(2u)-1)-v(\log(2v)-1),
\quad u=y/2-x,\ v=y/2+x.
\]
高度 1 上的周期因子上下有正界，但在复积分中不能抛弃其相位。A'_y=−log(y²−4x²)，极大点
\[
x_*=-\tfrac12\sqrt{y^2-1},\quad
A_y(x_*)=B(y)=y\operatorname{arcosh}y-\sqrt{y^2-1},\quad
A_y''(x_*)=-4\sqrt{y^2-1}.
\]
在 t 实部宽 O(√K) 的窗口 W_K 上，绝对积分至少 c_y√K exp[KB(y)]。

对固定 y≥y_c=1.50887956…，标准大阶数 Bessel 渐近给 E≤C_y√K exp[−KF(y)]，所以
\[
\mathcal E^{-1}\int_{W_K}|U(t)|\,|dt|
\ge c_y e^{K\Delta(y)},\qquad \Delta(y)=B(y)+F(y)>0.
\]
证明正性：Δ'=acosh y−asinh y<0，Δ(y)→0（首项 1/(2y)）。在 y=2，Δ≈0.2506620134229251。

双向乘子模长界将该下界传递给 Qν，仅付出 exp[−C_Y(ρ+K/ν)]。保留的二次相位也只以 exp[O_Y(K/ν)] 改变此条带的模长。

这不是完整联合渐近的反例。它排除了“点态多项式余项乘以当前轮廓绝对积分，就得到 E 尺度多项式误差”的直接路线。下一步必须利用相位、解析分支或可证明的轮廓变形。

## 7. 数值与继续位置

[代码](../code/recoil_bessel_endpoints.py)、[原始记录](../results/stage4/recoil-bessel-endpoints/diagnostics.json)保存版本、参数、精度与哈希。ν=32K^(3/2)，ρ=1/1024 固定，是联合界测试，不是 ρ→0 收敛序列。

K=8、32、128、512，y=0、1、2；竖直采样 s=0、.1、1、10、K、ν、M−ν、M、M+ν、2M。独立实乘积公式与复 log-Gamma 对比，避免巨大指数溢出；保留完整对数比。

300 项诊断通过：120 个竖直身份、120 个衰减不等式、12 个 λ≥π、36 个反射身份、10 个中央绝对求积加阶和 2 个增精度。另保存基点和局部乘子样本，不给未知常数编造数值界。

中央窗口为 |Re t+√3K/2|≤√K/4、高度 1、y=2。20/32 阶 Gauss–Legendre、单位或更短网格、55/75 位，最后一行再用 95 位。最大相对加阶变化<7.4e−27；不是严格误差证书。

|K|∫窗口 \|U\| / E|∫窗口 \|Qν\| / E|log(∫窗口 \|U\| / E)/K|
|---|---|---|---|
|8|9.33072e2|1.00797e3|.854810|
|16|5.65331e3|5.94872e3|.540000|
|32|2.77261e5|2.86796e5|.391647|
|64|7.95113e8|8.12952e8|.320219|
|128|7.15645e15|7.25825e15|.285209|

这是绝对范数诊断，不计算全复振幅，也不将这些大数解释为物理增长。常数 prefactor 在有限 K 仍明显影响对数速率。

下一步：在不跨错极点的前提下分解/变形中央解析分支，建立保留相位的 E 归一误差。随后才恢复 u/v、横向尾部、外部能量循环与物理插入归一化。有限 G 与完整全息仍开放。
