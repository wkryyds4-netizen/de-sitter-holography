# 固定电荷的记忆返回：守恒比较与全有限范数

2026-11-08。本阶段稿为 manuscript/susskind-sector-memory-return/main.tex，阅读版为 output/pdf/de-sitter-susskind-sector-memory-return.pdf。沿用同一 static Gaussian overlapping p-body CAR Hamiltonian、归一 annealed sector Schatten 范数、固定 Js、ℏ=1 和 string 时间。p²/N→固定正 lambda，|m−N/2|≤L0p、|t|≤S，L0、S 和有限范数阶 r≥1 固定。

## 本轮结果

实际记忆 M(t)=R†exp(itDperp)R，普通传播 Lambda(t)=R†exp(itK)R。保留逐样本完整简并块能量投影，定义 Fscr=(1−P_H)m exp(−iHt/2)M exp(−iHt/2)，X0 用 Lambda 替换 M。本轮得到该中央电荷窗口内所有固定有限 r 的

\[
\sup_{m,|t|\le S}\|\mathcal F(t)-X_0(t)\|_{r,m}=O(p^{-1/2}).
\]

它补齐 total memory 到普通传播的逐扇区返回接口，尚不是 X0 的两分量极限或完整 selected bilocal 定理。证明不将 full-Fock 平均界直接限制到固定电荷。

## Force defect 的精确逆式与时间尾

K_n=(p Ncal−D)H^n，其中 D 先作用完整电荷族。F_n=Gamma_(H,m)(H_(m−1)^n)−2p Phi_m(H_m^n)。Gamma 与 Phi 的中间输入来自不同扇区，不能省略这一差别。旧 operator-draft 式(21) 给 K_(n+2)=Σ_(a+b+c=n)H^a F_b H^c。

先减 HK_(n+1) 消掉 a≥1，再减 K_(n+1)H−HK_nH 消掉 c≥1，得到精确式

\[
F_n=K_{n+2}-HK_{n+1}-K_{n+1}H+HK_nH.
\]

K0=K1=0，因而 n=0、1 边界同样成立。这只是已有 product identity 的有序代数反演。sector-gradient-tail 的 (O) 已给 K_n 的 polynomial(n)L_r^n sqrt(n!) 包络；也可从其 G_n 加回 Gaussian number counterterm，用旧 OU/Chebyshev multiplier 界得到。逆式只用固定 2r/4r Hölder 与两次外部 H，额外 sqrt((n+2)!/n!) 吸收到多项式。故指数 force defect 的所有固定 r 范数在紧时间窗一致 O(1)，不止逐阶极限。其尾由 Σ_(n>L)poly(n)(L_rS)^n/sqrt(n!) 支配。

## 有限参数守恒比较

旧 Jacobi 恒等式和 T_j=p(Phi−C_m theta_j)Q_j(H) 给

\[
\Phi_m(e^{itH})=C_m f_m(t,H)+p^{-1}Z_m(t),\qquad
f_m=\sum_jc_j(t)\theta_jQ_j(H),\quad Z_m=\sum_jc_j(t)T_j.
\]

theta、Q、c 都保留 finite-N 系数。f 在每个样本上是 H 的函数；CP 收缩及已证 T 尾给 ||f||_r≤1+O(1/p)，无需一致 samplewise operator-norm 界。

定义 b_m(t,H)=(2p²/m)C_m f_m(t,H)。由 D_F=KV−VH、R=D_F+(p/m)VH，以及精确 D_F†exp(itK)D_F=Gamma(exp(itH_(m−1)))/m，得到

\[
\frac p m\Gamma(e^{itH_{m-1}})-b_m
=\frac p m F_m(t)+\frac{2p}{m}Z_m(t)=O_r(p^{-1}).
\]

旧逐sector Gaussian force Gram 高矩给 ||D_F||_(2r)=O(p^-1/2)、||(p/m)VH||_(2r)=O(p^-1)。两交叉项乘 p 后 O(p^-1/2)，平方项 O(p^-1)。于是 p Lambda−b_m=O_r(p^-1/2)，所有 fixed finite r 直接成立，不从带 p 增长的粗界插值。

## 返回展开

有限 Volterra 恒等式为 M=Lambda+J1+E1，其中 J1=Lambda*U*Lambda、E1=M*U*Lambda*U*Lambda，U(t)=exp(i alpha Ht)。构造

\[
\mathcal K_m(s,H)=\frac m{p^2}e^{-iHs}
\int_{v,w\ge0,v+w\le s}b_m(v,H)e^{i\alpha H(s-v-w)}b_m(w,H)\,dv\,dw.
\]

逐样本 [Kcal,H]=0。两次 2r Hölder 给中心化 mJ1−Kcal=O_r(p^-1/2)。E1 有六个 R/R† 因子，以统一源维数归一并取 6r Hölder，得 ||E1||_r≤s⁴/24 ||G||_(3r)³，G=R†R=O_r(p^-1)，所以 mE1=O_r(p^-1)。负时间用伴随与反向单纯形。最后 (1−P_H)Kcal=0 精确，且 ||1−P_H||_(r→r)≤2，得到主结果。

在半填充与已控制固定多项式制备 P 中，轨道密度满足 ||R_P||2统一有界，故 E||[Fscr−X0]P O||²≤||R_P||2 ||Fscr−X0||4²=O(p^-1)。这不自动控制任意 actual total-memory 产品中的其余因子，也不解决任意 nonpolynomial selected 交错。

Kcal 的有限参数初始两项是 2p²C_m² s²/m，以及 i[2p²/(3m)]C_m²(2theta1+alpha−3)Hs³；半填充极限再现旧 4lambda C²s²+i(8lambda C²/3)(q−1)Hs³。没有把旧 calibration 当作新谱公式。

## 检查与保存

新 pure helper code/sector_memory_return.py 和新 driver code/check_sector_memory_return.py 检查了21个精确整数 free-word 逆式；三 CAR/Hermite 系统共9行，全部有非零低chaos贡献；15个有正有负的实时间点及90条 r=1,2,4 Schatten 不等式。block matrix exponential 直接实现有序卷积，没有时间求积或有限时间Taylor截断。float64 CAR最大相对误差9.50e-16、绝对误差2.71e-15，Volterra最大相对误差3.70e-16。浮点检查不是区间证书或双标度证明。

科学输出首次成功即封存，SHA256 05599623446e6af7dc9a9cc00f1ab175703db3ae53ce8bf6705c45c1b6adf3b6。参数、版本和代码哈希在该JSON中。旧成功driver与manifest未重跑。正式稿双遍编译与全页核版见本轮local-evidence。

下一步：X0 的 charged alpha 接口与 central four-u 条件回归，再做实际 bilocal Volterra 返回；空间/质量字典和finite-G仍开放。本轮仅本地，不重复上一轮686b199/issuecomment-6712633148的发布。
