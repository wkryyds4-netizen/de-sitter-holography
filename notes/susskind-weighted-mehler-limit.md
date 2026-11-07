# 带权 q-Mehler 极限：把选态预测补成余项受控的结论

2026-11-08。正式证明见 [七页 LaTeX](../manuscript/susskind-weighted-mehler-limit/main.tex) 和 [PDF](../output/pdf/de-sitter-susskind-weighted-mehler-limit.pdf)。本阶段承接 selected-rindler-response；旧稿的条件性措辞保留，不追改历史证据。

## 对象、单位和结论

取 q=exp(−epsilon)，epsilon=4lambda，e=H/√C，K=H/(2Js)=e/√epsilon，u=2Js t，hbar=1。先做有限 N 的实际精确投影，再取半填充、同 J、固定 lambda 的双标度极限，最后 epsilon→0。本文处理这个极限核，没有给 N 一致速率。

对任意固定正整数 d，归一两能量 q-Mehler 测度的参数为 q^d。C_d 为迹态关联，Ctilde_d 在输入能量乘正权 f(e)=1−√q+√q e²。d=2 时，这正是此前接通的实际曲率单选态；一般 d 此处只定义能量倾斜，不增添微观选态解释。

令 g_d(u)=sech^(2d)(u)，R_d(u)=(d/2)g_d(u)[u tanh(u)−(2d+1)tanh²(u)]。固定 d、S 和有限导数阶 k 时，

- C_d=g_d−epsilon R_d/2+O_{C^k([-S,S])}(epsilon²)；
- Ctilde_d=g_d+epsilon R_d/2+O_{C^k([-S,S])}(epsilon²)。

因此 d=2 的选态差为 epsilon sech⁴(u)[u tanh(u)−5tanh²(u)]+O(epsilon²)，上一阶段尚有条件的预测获得直接积分证明。真实正态混合 (tau+sigma1)/2 在此曲率通道消掉一阶偏差，平均关联为 sech⁴(u)+O(epsilon²)。这不是体真空或质量的识别。

## 证明中的关键步骤

已有 q-Mehler 乘积先作整数有限移位，然后用 theta 模变换。坐标 a=theta+phi−pi=√epsilon y、delta=theta−phi=epsilon z 将主密度分成标准 Gaussian phi(y) 与

f_d(z)=2^(2d−1) z product_{j=1}^{d−1}(j²+z²)/[(2d−1)! sinh(pi z)]。

反射像必须成对处理：T(x)=Gaussian·[1−exp(−2pi b(x)/epsilon)]·Qim(x)，其中 b=min(x,2pi−x)，|Qim−1|≤4exp(−2pi²/epsilon)。delta 的边缘因子已吸收到 sinh(pi z)，不能重复。剩余的 sin((theta+phi)/2)/F(theta+phi) 一致有界，得到全谱 p≤C_d phi f_d，包括靠近端点的薄层。

首阶相对密度 A_d=(4d+1)/8−(2d+1)y²/8−z²/2。对每个固定多项式权 w，有 integral w|p−phi f_d(1+epsilon A_d)|=O(epsilon²)。证明在 |y|≤4√log(1/epsilon)、|z|≤8log(1/epsilon)/pi 展开；中心余项被 epsilon² 乘固定多项式支配，外部 Gaussian/指数尾为 epsilon^8 乘对数幂，因此没有保留对数损失。

对称化能量插入 W0=y²−1、W1=z²+y²/2−y⁴/12；相位为 −2z[1+epsilon(1/4−y²/8)]。Gaussian 积分的首项消去，首阶留下 z²−d/2+iuz/2。傅里叶变换给 −g_d''/4−dg_d/2−ug_d'/4=R_d。没有插入时为 −R_d/2。全谱多项式控制也覆盖固定阶时间导数。

## 已执行检查与保存

modular-checks.json：27 个原始乘积/模表示比较、6 个边缘/对角层、18 个局部展开、一般 d 的符号 Gaussian 相消、12 个独立 Fourier 积分。乘积比较最大相对差 1.36e−68，边缘比较 6.83e−67。Fourier 使用 50/70 位及 |z|≤32，最大差 2.76e−34、显式尾界最大 2.83e−33。数值积分不声称区间认证。

trace-mixture-checks.json：到 u^14 的 16 个精确迹态系数（主项和一阶项），只读复用旧 28 个时间点；混合误差除 epsilon² 的最大值约 0.446075。未重跑成功科学 driver。

双遍 pdflatex、90dpi 七页渲染全部检查。首次编译因集合左括号缺转义失败，原日志留 tmp；修正后无未定义引用及公式溢出。当前科学输入、PDF、编译与渲染哈希由 local-evidence.json 绑定。

## 下一步

全双局域态响应需保留 chi=e−W 的混合源。先从已完成的共同代数和 bilocal 重建提取中性项，再研究其小 epsilon 极限；曲率通道的混合相消不能直接推广。非整数 d、增长时间、联合 N/epsilon 极限、任意 selected 交错、微观空间/质量生成元及有限 G 尚未完成。
