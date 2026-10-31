# 全局反冲积分：第二分支与 Barnes–Bessel 辅助模型

2026-10-30；Alan Day、Chengheng。内部研究记录，作者自查草稿。配套[LaTeX](../manuscript/recoil-global/main.tex)、[PDF](../output/pdf/de-sitter-recoil-global.pdf)、[诊断](../results/stage4/recoil-global/diagnostics.json)、[来源](../literature/stage4-recoil-global-sources.md)。

## 本轮改变了什么

上一轮[真实局部积分](stage4-recoil-local.md)只控制右鞍点附近。本轮保留完整积分的负方向后，在 y=2、K=4、8、16 的有限参数计算中发现远大于局部项的贡献。因此下一步不能预设“对全部有界 y，右局部项主导”。局部证明不被这些结果否定；有限数值也还不是联合渐近反例的证明。

对象仍是固定质量先取 κ→∞ 后的内部手征图积分。de Sitter 半径为1，u=v=0，K=ν²/(2M)，ω=iKy，M>ν，y实。没有指定新的物理态、算符代数、外部插入和能量循环；未得到有限G的物理OTOC。复杂质量时保留手征解析继续，不能改成新的模方。

## 精确负方向反射

沿用[内部积分](stage4-recoil-internal.md)及[鞍点定义](stage4-recoil-saddle.md)：
\[
B(z)=e^{-i\pi z}
\frac{\prod_a\Gamma(z+ia)\prod_r\Gamma(\frac12-z-ir)}
{\Gamma(z+\frac12)},\qquad
a_j=M+\nu/2,\quad r=(M+\nu+\omega/2,2M,M+\nu-\omega/2).
\]
令w=−z，逐个使用Gamma反射公式得
\[
B(-w)=8\pi^3e^{4\pi iw+\pi\sum r}
\frac{\Gamma(w+\frac12)\prod_r\Gamma(w+\frac12-ir)}
{\prod_a\Gamma(1+w-ia)}
\frac{1+e^{2\pi iw}}{\prod_a(1-e^{2\pi i(w-ia)})}.
\]
这一步精确，未删去极点因子。实a和共轭不变r下，兼容对数支给
\[
\Psi(w)=2\pi iw-\overline{\Phi(\bar w)},\qquad
D=\Re\Phi(z_+)-\Re\Psi(\bar z_+)=2\pi c_++2\Re\Phi(z_+).
\]
其中z₊=x₊−ic₊。D是右减负方向的实指数差，不是已经确定的鞍点系数。ε=ν/M、δ=M+ν−c₊时，对称小ε展开为
\[
D/K=2F(y)-2\pi\delta/K+O(\epsilon^2),\quad
F(y)=\sqrt{1+y^2}-y\,\operatorname{arsinh}y.
\]

## 用于数值的完整轮廓

取L=x₊>K|y|/2，c=M+ν−d，满足max a<c<min Re r。左尾沿x=−L向下变形；扫过的左下区域无极点：第一组极点在高度−a>−c，另一组在−L右边。固定参数反射/Stirling估计去掉远弧。向上左射线的定向给
\[
\int_{\mathbb R-ic}B\,dz
=\int_{-L}^{\infty}B(x-ic)\,dx
+i\int_0^\infty B(-L-i(c+t))\,dt.
\]
代码除以g₊expΦ₊后求积，最终除以sqrt(2π/(−Φ″₊))得到完整项/局部首项。固定参数变形不意味着质量一致的尾界。

## 紧集上的核极限与精确辅助积分

设t=z+i(M+ν)、α=Ky，并取
\[
C_{M,\nu}=e^{-\pi(M+\nu)}
\frac{\Gamma(-i\nu/2)^4\Gamma(\frac12-i(M-\nu))}
{\Gamma(\frac12-i(M+\nu))}.
\]
固定K、α，ν→∞且M=ν²/(2K)，远离极点的紧集上
\[
B(t-i(M+\nu))/C_{M,\nu}\to
U_{K,\alpha}(t)=e^{-2\pi it}(K/2)^{2t}
\Gamma(\tfrac12-t+\alpha/2)\Gamma(\tfrac12-t-\alpha/2).
\]
这只是核的紧集收敛，还没有交换完整积分与极限。

辅助轮廓从d₀−i∞向上到d₀+iη，再向右到+∞+iη；d₀<1/2−|α|/2、η>0。下尾具有e^(−3π|Im t|)衰减，右尾反射后有倒Gamma的阶乘衰减。顺时针向右下闭合，非整数α的两族极点分别为1/2±α/2+n。两族留数之和为
\[
\begin{split}
J={1\over2\pi i}\int U\,dt
={}&-e^{-i\pi\alpha}\Gamma(-\alpha)(K/2)^{1+\alpha}
{}_0F_1(1+\alpha;K^2/4)\\
&-e^{i\pi\alpha}\Gamma(\alpha)(K/2)^{1-\alpha}
{}_0F_1(1-\alpha;K^2/4)\\
={}&-K[e^{i\pi\alpha}K_\alpha(K)+i\pi I_\alpha(K)]
=-KK_\alpha(Ke^{-i\pi}).
\end{split}
\]
最后一个K函数是负轴下方的解析继续。整数α必须继续整个和，不能分别继续发散的留数系数。这是标准特殊函数关系的应用，不作原创性宣称。

固定y>0，K→∞，标准大阶渐近给
\[
I_{Ky}(K)\sim\frac{e^{KF(y)}}{\sqrt{2\pi K\sqrt{1+y^2}}},
\quad K_{Ky}(K)\sim\sqrt{\frac{\pi}{2K\sqrt{1+y^2}}}e^{-KF(y)}.
\]
两项模长比K_(Ky)/(πI_(Ky))约e^(−2KF(y))。F'(y)=−asinh(y)，唯一正零点y_c=1.5088795615383199289。它只标记辅助模型的指数模长交叉，不能称为物理临界点或完整图已证实的Stokes跃迁。

## 首阶紧集修正

Gamma比值展开给
\[
B/C=U[1+(4i/\nu)(t^2-t)+O_{K,\alpha,t}(\nu^{-2})].
\]
四个−iν/2大参数贡献4i(t²−t)/ν；另外两个参数为O(M)，在固定K下从ν⁻²阶进入。令D_K=(K/2)∂_K，必须固定α，则D_K U=tU，辅助积分的尾衰减允许求导。Bessel方程精确给
\[
(D_K^2-D_K)J=(K^2+\alpha^2-1)J/4.
\]
因此形式上，归一比B(K,α)/B(K,0)的首阶修正为iα²/ν。数值中使用exp(iα²/ν)保留相位；不是全阶重求和定理，也尚未证明Barnes余项可一致积分。

## 数值记录和真正支持的主张

取ν=4K³。所有点u=v=0：

| K | y | 完整项/右局部首项（约值） |
|---|---|---|
|4|0|1.0377303−0.0016775i|
|4|2|0.7574552−182.4665112i|
|8|0|1.0169041−0.0001656i|
|8|1|1.0014509−0.0003609i|
|8|2|−4.4614748−33541.8312194i|
|16|0|1.0081070−0.0000205i|
|16|2|−22852.0717−1.1232773362×10⁹i|

以R为完整归一手征比，R_J为J(K,Ky)/J(K,0)：

|K|y|未修正相对差|乘exp(i(Ky)²/ν)后的相对差|
|---|---|---|---|
|4|2|0.24889358|0.00389292|
|8|1|0.03124995|0.00007430|
|8|2|0.12488755|0.00049349|
|16|2|0.06248789|0.00006189|

没有拟合相位。诊断共34项：12项负反射/指数身份、6项辅助Bessel/微分身份、7项完整轮廓求积细化、9项高度/范围/精度比较。求积20/36阶，50/65位，d=1，X=max(32,8K)，T=18；y=2增加d=.75、X=max(32,10K)、T=24和90位比较。各段相消比最大约60.5，比较实际复和的相对变化。完整积分实际采用有限范围近似；扩展范围变化不构成无限尾的区间认证。

## 下一步的可证伪目标

先证明固定K、α的全Barnes→Bessel极限：建立质量一致的两端尾界，尤其正水平尾，考虑x≪ν、ν≲x≲M、x≳M三个区域及交叠。随后才允许K增长、建立积分首阶修正和两分支统一误差；两项可相消处先控制绝对误差再取相对误差。再恢复u、v及横向盒外积分。物理插入、外部能量循环、断开基准、有限G和微观字典继续开放。
