# 同种 singlet 的固定时间首阶四点

2026-11-06。主稿见 manuscript/susskind-first-order/main.tex，五页 PDF 已编译至 output/pdf/de-sitter-susskind-first-order.pdf。本阶段沿 Susskind complex DSSYK∞ / Rindler meson 主线。

采用此前守恒分量扣除后的 X=(O−cB)/√W_d，完整 annealed 归一迹 β=μ=0；先固定 λ=p²/N>0 双标度，再固定 string 时间取小 λ。设 h=4λ、u_j=2J_s t_j、d=2，g_ij=sech^{2d}(u_i−u_j)、T_ij=tanh(u_i−u_j)。连通函数扣除完整 finite-h 两点的三个乘积，而非只扣领先 sech。

一般四时间结果为 κ₄=hK₄+O(h²)，

\[
K_4=-d^2\left[g_{12}g_{34}T_{12}T_{34}+g_{23}g_{41}T_{23}T_{41}
+g_{13}g_{24}\frac{\cosh(u_1+u_3-u_2-u_4)}{\cosh(u_1-u_3)\cosh(u_2-u_4)}\right].
\]

这是将 Lin–Stanford / Streicher 现成半经典公式应用到当前已归约 singlet 的结果，不是新的普适四点定理。O(h²)表示固定时间渐近，没有证明统一余项界。

因此，C(t)=τ([X(t),X]†[X(t),X]) 为

\[
C(t)=32\lambda[\cosh(4J_st)+(4\tanh^2(2J_st)-1)\sech^8(2J_st)]+O(\lambda^2),
\]

而相邻时间连通函数为

\[
\kappa_{22}(t)=16\lambda(2\tanh^2(2J_st)-1)\sech^8(2J_st)+O(\lambda^2).
\]

两个排序的 disconnected 部分同为1+2G_h(t)²，所以两点自身的首阶修正在 C 中严格消去。三个收缩通道都必须保留。前者短时给1280J_s²λt²，增长重叠给16λexp(4J_s t)；后者给κ22(0)=−16λ、κ22″(0)=768J_s²λ。领先κ22零点J_s t=.4406868，不直接等于有限λ零点。

推导的两个关键位置：crossed由Lin–Stanford(205)在c→1、τ+=i(u3−u1)、τ+'=i(u4−u2)、τ−=i(u1+u3−u2−u4)给出；uncrossed由附录G.2能量交换及τ(Hsrc XiXj)=−i∂sG/2给出。第二个uncrossed用有方向的T41，不能改成T14。互交叉q^{d²}已在crossed连通首阶中，不能再乘一次。

代码 GeneralKernel 将旧交替时间实现推广到任意四个实时间，crossed保留复双线性积，不能普遍取绝对平方。独立50位非对称BINT指标和在L14、两个λ和两组非交替时间下核对12个通道，误差均<1e−12。精确弦多项式另核对22个Taylor系数及72个非正Laurent次幂相消；四条时间射线覆盖至u¹⁰或u⁸。

四组λ=.04,.02,.01,.005采用L169,332,672,1383，包含16个交替时间点和12组一般四时间；其中λ=.02,.01的交替点直接复用旧输出。选点L+16最大通道差1.23e−15。一般crossed内部最大相消条件2.20，连通减法373。这些是经验稳定性与表示核对；原cutoff输出中的解析界针对旧交替核，一般四时间总界尚未完成证明，不将它挪作认证误差。

有限λ偏差明确存在：J_s t=.5时，λ=.005的C/h为24.850065，首阶系数30.426106，低18.33%；以.005和.01作一次Richardson外推得28.771611，仍低5.44%。相邻κ/h对应外推.01986431，对照.01991637，差约.26%。一般四时间第三组在整个λ扫描中并非单调，但从.01到.005三组系数误差均缩小。没有为美化收敛而删点。

物理意义：已得到具体微观四点靶标，尚未把它变成此前空间KG波包的乘法保持字典。任意正质量仍可由径向伸缩配合，下一步必须加入独立空间/平移生成元及交叉关联，检验Rindler平移代数和质量观测量。不要仅由d=2报meson质量，也不要把已有时间核的扫描继续当作主线。

构建：powershell -NoProfile -ExecutionPolicy Bypass -File code/build_singlet_first_order.ps1；pdftoppm -scale-to 1200 -png。最终五页全部检查，引用无未定义、无公式溢出，输入输出哈希见本阶段local-evidence.json。无commit/push/ZIP。
