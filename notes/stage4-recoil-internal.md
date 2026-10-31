# 固定质量纯反冲图的 Barnes 积分

2026-10-30。作者自查证明草稿，无独立认证或优先权声明。
[四页PDF](../output/pdf/de-sitter-recoil-internal.pdf)、[LaTeX](../manuscript/recoil-internal/main.tex)、[来源核查](../literature/stage4-recoil-internal-sources.md)。接续[合并顶角](stage4-recoil-vertex.md)，本轮补的是完整内部积分，不是外部三能量积分或物理OTOC。

## 精确变量与轮廓

三维、de Sitter半径1、CLP固定拓扑图。六边质量μⱼ为固定正实数；b²=i/κ、L=b⁻²=−iκ、iPⱼ=−1/(2b)−ibμⱼ/2。顶点集合V={123,156,246,345}，四边集合R={1245,1346,2356}，定义aᵥ=(Σᵥμ)/2、rρ=(Σρμ)/2，满足Σa=Σr。要求固定严格间隙max a<min r，c=(max a+min r)/2。

在CLP (3.46)的内部变量中设iP=−3/(2b)+bz，dP=−ib dz。被积函数精确为

\[
K_\kappa(z)=\frac{\prod_vS_b(b(z+ia_v))}
{\prod_\rho S_b(b(z+\frac12+ir_\rho+L))\,S_b(b(z+\frac12-L))}.
\]

可能的极点族为z=−n+iκm−iaᵥ、z=1/2+n−iκm−irρ及z=1/2+n−iκ(m+2)，n,m≥0；有些会被零点消去。CLP (3.58)的s轮廓经s=−8πi−4π(z−1/2)/κ变成**z∈ℝ−ic**，不是一开始就竖直的Barnes线。大κ时前族在上，其余在下。方向取左到右，连同dP Jacobian跟踪；相反约定只改整体符号。

## 精确 q-Gamma 表达及积分极限

q=e⁻²π/κ、τ=e⁻²πκ。CLP A.5精确给

\[
S_b(bw)=C_\kappa(1-q)^{w-1}
e^{i\pi w/2+\pi(w^2-w)/(2\kappa)}
\Gamma_q(w)(\tau e^{2\pi iw};\tau)_\infty .
\]

Cκ与w无关。做三次正逆b移位和一次负逆b移位，归一
\(\mathcal B_\kappa=-8\pi^3e^{-\pi\kappa}(1-q)^2K_\kappa\)。
其精确q-Gamma、三角和模对偶乘积见英文式(9)。四个小分子参数A=z+ia和四个分母参数D=z+1/2+ir、z+1/2的和差为−2，平方减一次项之差Δ₂中z系数为−4。紧集上的极限是

\[
\mathcal B(z)=e^{-i\pi z}
\frac{\prod_v\Gamma(z+ia_v)\prod_\rho\Gamma(\frac12-z-ir_\rho)}
{\Gamma(z+\frac12)}.
\]

本轮不止取点态极限。q=e⁻ʰ，f_q=logΓ_q，在实x>0有

\[
\log\frac{1-q^x}{1-q}-\frac{h}{e^{hx}-1}\le f'_q(x)
\le\log\frac{1-q^x}{1-q},
\]

而沿有界复平移的f″由1/(x−A)+1/(x−A)²控制。对参数和差为−2的四比四Gamma_q商作Taylor展开，得到C[(1−q)/(1−e⁻ʰˣ)]²。精确式的剩余指数为e⁻ʰˣ，乘积≤C/x²。负尾先反射所有八个原S_b、令w=−z，再移位；其小参数和差仍为−2，独立得到同一界。不能凭正尾对称猜负尾。

固定质量及固定间隙下，因而
\(|\mathcal B_\kappa(x-ic)|\le C/(1+x^2)\)，C对大κ一致但依赖质量。
支配收敛给

\[
-\frac{4\pi^2(1-q)^2e^{-\pi\kappa}}bJ_b
\longrightarrow\mathfrak B(\boldsymbol\mu)
=\frac1{2\pi i}\int_{\mathbb R-ic}\mathcal B(z)\,dz
=\frac1{2\pi i}\int_{d-i\infty}^{d+i\infty}\mathcal B(z)\,dz,\quad0<d<\tfrac12.
\]

最后一步是对**极限函数**的轮廓旋转：右端向上、左端向下，无极点穿越，连接弧O(R⁻¹)趋零。竖线两端分别指数衰减e⁻²πy和e⁻⁴π|y|乘多项式。未宣称有限κ积分也可以用同一竖线。

## 可计算表达与完整图比

不同aᵥ时，向左闭合得四项₄F₃(1)，系数及参数完整保存在英文式(17)—(18)。每项的参数过剩为1，收敛；多个a相等时单项系数发散，只能对完整和取极限或直接用竖线积分。Barnes/Racah表达有既有前驱，不把该函数命名为新的6j，也未建立它与某个群的酉主系列、测度和归一化字典。

结合上一轮T=|F|的合并顶角，对同样满足间隙且B不为零的参考质量μ*，

\[
\lim_{\kappa\to\infty}\frac{Z_{\rm tet}(\mu;\kappa)}{Z_{\rm tet}(\mu^*;\kappa)}
=\prod_v\frac{T(\mu_v)}{T(\mu_v^*)}
\left|\frac{\mathfrak B(\mu)}{\mathfrak B(\mu^*)}\right|^2.
\]

这给整个固定拓扑图在固定实质量下的归一极限。模方只在实Euclidean数据上使用；复杂外部能量须分别继续两手征，不可重新取模。没有指定物理四点外腿或断开基准。

取边映射(μ₁,…,μ₆)=(m₁,ν,m₂,m₃,ν,m₄)，其中m₁,₃=M+ω/2±u、m₂,₄=M−ω/2±v，则四个a=M+ν/2+(±u±v)/2，独立于ω；三个r=M+ν+ω/2、2M、M+ν−ω/2。u=v=0会使四族左极点全部合并。此图重标号可能与CLP的ω号相反，此前Gram结论为偶函数；未来非偶量须明确用当前映射。

## 证据与局限

[34项诊断](../results/stage4/recoil-internal/diagnostics.json)由[代码](../code/recoil_internal.py)保存：4个符号平衡、9个精确有限κ移位及各自精度复核、3组极限积分复核、1个四项超几何对照、8个双侧尾复核。45位工作、指定量65位复核；未直接求完整有限κ积分，也未作区间认证。点态趋近误差和算术一致性分开保存。

首个积分样本无意中有三族相同，因此最初32项没有执行不同族的超几何比较。检查原始参数后补入真正四族不同的第三组，最终34项通过，原先两组保留。第三组四项绝对值和/最终绝对值约268134.37；积分与完整和仍相符，不能只看单项或只报最终数值精度。所有输入、依赖和哈希保存。

下一步在同一Barnes表示做M、ν大参数的内部鞍点与涨落，保留u、v，再接横向边缘积分。C/(1+x²)的C尚未对增长质量或复杂能量统一，不能把固定质量κ极限自动升级为联合有限G误差。物理插入、外部能量轮廓及微观态字典仍开放。
