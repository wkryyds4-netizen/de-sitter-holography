# 保留相位的两分支联合估计

2026-10-30。Alan Day、Chengheng。作者自查证明草稿，未经过独立审查，不作优先权声明。[LaTeX](../manuscript/recoil-bessel-joint/main.tex)、[PDF](../output/pdf/de-sitter-recoil-bessel-joint.pdf)、[代码](../code/recoil_bessel_joint.py)、[数据](../results/stage4/recoil-bessel-joint/diagnostics.json)。

## 1. 本轮主张与对象

固定质量先 κ→∞ 的内部手征 Barnes 图，半径1，K,ν>0，M=ν²/(2K)，u=v=0，ω=iα。没有新增物理态、算符代数、外插入或断开归一化，不能称有限G四点。精确 Qν 及 Iν=(2πi)^(-1)∫Qν 沿用[固定K笔记](stage4-recoil-bessel-limit.md)；旧图归一 C 独立于α。

令
\[
A_K=KK_\alpha(K),\quad A_I=KI_\alpha(K),\quad
J=-e^{i\pi\alpha}A_K-i\pi A_I,\quad\mathcal E=A_K+\pi A_I.
\]
α≥0；负α利用完整表达的偶性。

**联合领先估计。** 对每个有限Y≥0，存在K_Y、ρ_Y、C_Y、c_Y>0，对K≥K_Y、0≤α≤YK、ρ=K³/ν²≤ρ_Y，有
\[
\left|\mathcal I_\nu-e^{iS/\nu}J\right|
\le C_Y(\rho+e^{-c_YK})\mathcal E,\qquad S=K^2+\alpha^2-1.
\]
因此K→∞、ρ→0足够，**不要求K²/ν→0**。常数没有认证数值门槛；这是领先联合估计，不是把旧固定参数高阶展开的全部系数升级为一致展开。J=0处只声明包络绝对误差。

## 2. 精确分支身份

记a±=(1±α)/2，b=1/2−i(M−ν)，c=1/2−i(M+ν)，并定义
\[
G(t)=(K/2)^{2t}\Gamma(a_+-t)\Gamma(a_--t),
\]
\[
V(t)=e^{-i\pi t}(K/2)^{2t}
\frac{\Gamma(a_+-t)}{\Gamma(\frac12+t+\alpha/2)},
\]
\[
R_\nu(t)=e^{i\pi t}(K/2)^{-2t}
\frac{\Gamma(t-i\nu/2)^4}{\Gamma(-i\nu/2)^4}
\frac{\Gamma(b-t)\Gamma(c)}{\Gamma(b)\Gamma(c+t)}.
\]
Qν=Rν exp(−2πit)G。利用
\[
e^{-2\pi it}+e^{i\pi\alpha}
=2e^{-i\pi(t-\alpha/2)}\cos\pi(t+\alpha/2)
\]
和Gamma反射，得到精确身份
\[
e^{-2\pi it}G=-e^{i\pi\alpha}G+2\pi e^{i\pi\alpha/2}V.
\]
V只剩a₊一族实极点，α为整数时也无发散连接系数。这一步把上一轮固定高度上的指数相消转化为两个分别可控的解析对象，而非删除快速相位。

旧端点推导的增长区域Gamma估计给，在任何固定C的|t|≤CK上，
\[
R_\nu=e^{iA(t)/\nu}(1+\epsilon_\nu),\quad
A(t)=4(t^2-t),\quad |\epsilon_\nu|\le C_C\rho e^{C_C\rho}.
\]
取消共同小Gamma后取连续对数。K/ν→0保证大参数的安全扇区；矩形内|Rν|≤exp[C_Cρ+C_CK²/ν]=exp[o(K)]。实t上的领先因子模长为1。不能把真实Rν竖线积分无限延伸穿过ν/2处极点。

## 3. 有限轮廓与所有边界

L=√(K²+α²)/2，d=1/2−L，N=ceil L，r₀=8+2Y。在[r₀K,r₀K+1]选择r，与两组a±+整数格点均相距至少1/8；两个坏集合总长不超过半个单位区间，因此可选。

固定充分大H=H_Y，T=HK。增大K_Y后保证T<ν/2、T<M−ν、N>α/2+1/2。ρ_Y同时包括旧正尾引理所需ν²≥32[max(2,Y²/2)+1]K³。

原完整轮廓可截为
\[
A=-N-iT\to B=-N+i\to C_r=r+i
\]
而误差≤C_Y exp(−cK)E：正尾和负下尾使用上一轮已证界。

G支将路径改为A→d−iT→d+iT→r+iT→C_r，不跨极点。a±极点在实轴右侧；分子极点高度ν/2>T，b−t极点在−T下方。V支用C_r→r−iT→A闭合为顺时针矩形，只包围t_n=a₊+n<r，留数为
\[
\operatorname{Res}_{t_n}(R_\nu V)
=-e^{-i\pi a_+}w_nR_\nu(t_n),\quad
w_n=\frac{(K/2)^{1+\alpha+2n}}{n!\Gamma(1+\alpha+n)}>0.
\]
因此
\[
\mathcal I_\nu=
-\frac{e^{i\pi\alpha}}{2\pi}\int_{-T}^{T}R_\nu(d+is)G(d+is)\,ds
-2\pi i\sum_{t_n<r}w_nR_\nu(t_n)
+O_Y(e^{-cK}\mathcal E).
\]

边界不能省略，具体控制如下：

- 顶/底横线x±iHK，−N≤x≤r：Stirling中K的幂抵消，
  log|G|≤K{−πH+C_Y[1+log(1+H)]}+O_Y(logK)。H大于有界缩放实部的两倍后安全扇区及C_Y不依赖H，取H充分大可使主指数<−3K。
- 右竖线：反射和极点间隙给|G(r)|≤C_Y exp[−Kh_y(r/K)]，h_y沿用上一轮，至少h₀=2r(log(2r)−1)。对任意实s，|G(r+is)|≤|G(r)|，由Gamma移位到正实部及有限递推因子推出。
- V底线及右竖线−T≤Im t≤1上，精确身份给|V|≤C|G|。G的短底连接段也在上述横线界内。

乘子只有exp[o(K)]，路径长O(K)，E≥c_Y√K exp(−K)，故所有连接边均为包络尺度的指数小量。没有穿过b−t远端极点，也没有对原高度1的中央绝对范数作三角估计。

## 4. 正留数支：用正权重消去首误差

标准Bessel级数给2Σw_n=A_I。令D=(K/2)∂K，**固定α**，则Dw_n=t_nw_n；Bessel方程给4(D²−D)A_I=SA_I。

以p_n=2w_n/A_I定义概率权重，q_n=A(t_n)−S，有精确矩：
\[
\mathbb E q_n=0,\qquad \mathbb E q_n^2=8K^2\mathbb E t_n.
\]
二阶身份来自对4(D²−D)与S作交换，DS=D²S=K²。另由级数移位，
\[
\mathbb E[n(n+\alpha)]=K^2/4,\quad
\mathbb E n\le K/2,\quad
\mathbb E t_n\le(1+\alpha+K)/2.
\]
故二阶中心矩O_Y(K³)。实数q上|exp(iq/ν)−1−iq/ν|≤q²/(2ν²)，一阶矩精确抵消，得到
\[
\left|2\sum_{n\ge0}w_ne^{iA(t_n)/\nu}-e^{iS/\nu}A_I\right|
\le C_Y\rho A_I.
\]
真实乘子只在t_n<r有限和中替换。省略的**辅助**权重之比为K²/[4(n+1)(n+α+1)]，到r₀K后统一小于1，首项Stirling给指数小尾。保留有限Rν，不使用其全无穷级数。此支总误差为O_Y(ρA_I+exp(−cK)E)。

## 5. 负竖直支：复中心矩与Gaussian模长

t=d+is处两Gamma实部A±=L±α/2∈[c_YK,C_YK]；实Stirling给G(d)=2πexp[−KF(y)+O_Y(1/K)]。

标准Gamma模长乘积
\[
\log\frac{|\Gamma(a+is)|}{\Gamma(a)}
=-\frac12\sum_{j=0}^\infty\log\left(1+\frac{s^2}{(a+j)^2}\right)
\]
在|s|≤HK给|G(d+is)|≤G(d)exp(−c_{Y,H}s²/K)。例如log(1+z)≥z/(1+z)，保留首O(K)项即可。大于HK的尾部，用同一递减正和的积分下界，其衰减率有正下界；任意固定多项式矩仍为exp(−cK)小。

结合Kα正积分的下界，∫|G|/(2π)≤C_YA_K，没有上一轮中央段的指数条件数。标准Mellin–Barnes表示给∫G/(2π)=A_K；先在局部固定竖线、固定α下求导，再移回d，Bessel方程给
\[
\frac1{2\pi}\int_{\mathbb R}[A(d+is)-S]G(d+is)\,ds=0.
\]
这是复积分身份，不是正权重断言。

选择d=1/2−L使
\[
q(s)=A(d+is)-S=-4s^2-8iLs.
\]
复Taylor余项由|q|² exp(8L|s|/ν)/(2ν²)控制。Gaussian完成平方后，
\[
\frac1{2\pi}\int_{-HK}^{HK}|G|\,|q|^2e^{8L|s|/\nu}\,ds
\le C_YK^3A_K e^{C_Y\rho}.
\]
平移Gaussian中心代价L²K/ν²=O_Y(ρ)，所有多项式矩有同样一致控制。利用全辅助积分的一阶中心矩为零，截断的多项式尾指数小，故该支是
\[
e^{iS/\nu}A_K+O_Y[(\rho+e^{-cK})A_K].
\]
Rν到首相位的点态误差可在这个良性竖线上积分，贡献同阶O_Y(ρA_K)。两支加上连接边即证第1节估计。

## 6. 对研究窗口的意义

K=ν²/(2M)→∞和ρ=ν⁴/(8M³)→0等价于
\[
\nu\to\infty,\qquad \nu^{4/3}\ll M\ll\nu^2.
\]
此前只是系数尺度支持的候选窗口，现在有指定对称内部图的完整领先联合估计。并未处理有限κ与质量联合极限，也未处理真实物理外部归一化。

α=0处两支正交，|J(K,0)|≥E(K,0)/√2。于是当误差充分小时，
\[
\frac{\mathcal I_\nu(K,\alpha)}{\mathcal I_\nu(K,0)}
=e^{i\alpha^2/\nu}\frac{J(K,\alpha)}{J(K,0)}
+O_Y[(\rho+e^{-cK})\mathcal E(K,\alpha)/\mathcal E(K,0)].
\]
旧C独立α因而抵消。这仍只是内部手征图比值。靠近辅助零点不除以J，不声称求出了完整零点位置或小抬升的相对误差。

## 7. 验证、失败保留与下一步

最终64项检查：12个精确分支/乘子身份，12个正权重矩，9个竖线矩，6个有限轮廓核对，6个冻结全积分比较，12个联合求积细化，6个增精度/扩范围，1个辅助零点核对。

三个有限轮廓故意取小高度T=1.2，边界修正/E分别约.0018048、.0390218、.0111411，不可忽略。初版两组身份误差约3.49e−25和4.54e−29，未达1e−35阈值；保存[初版输出](../results/stage4/recoil-bessel-joint/initial-contour-grid.json)和[对应代码](../results/stage4/recoil-bessel-joint/initial-contour-grid.py)，将网格由最多.5缩为.25后重跑，最终最大约2.26e−39，没有降低阈值。

新的两支表示与旧6组全积分比较，最大包络差约5.88e−22。联合序列ν=8K^(7/4)、ρ=1/(64√K)，K²/ν=K^(1/4)/8会增长；不同于上一轮ρ固定的界测试。

|K|包络误差/ρ，y=0|y=1|y=2|
|---|---|---|---|
|8|.392828|1.002373|3.615498|
|16|.364228|.974595|3.665900|
|32|.349064|.959404|3.691384|
|64|.341276|.951441|3.705275|

误差使用完整复差/E。20/32阶、55/75位，最大K另95位及T=8K→10K。最大加阶包络差约5.84e−20，增精度约2.36e−70；扩范围在该工作精度下未分辨出变化，不报为严格零误差。数值T=min(8K,ν/4)、r近12K，未认证为定理常数给出的最小截断门槛。

两支内部条件数最大分别约1.063和1.026。α=15.5的辅助零点K≈10.2731额外检查，实际包络误差约.009615，合并相消条件数约104；保留复和与绝对误差，不用首项零点做分母。这些都是经验诊断，不是区间认证。

下一主线恢复u/v及横向尾部，之后接外部能量循环、插入/断开归一化和有限G/微观字典。完整de Sitter全息尚未完成。
