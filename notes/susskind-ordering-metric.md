# 排序奇内积、相对 Casimir 与平移相容性

2026-11-06，本地阶段。正式推导见 [LaTeX](../manuscript/susskind-ordering-metric/main.tex) 与 [PDF](../output/pdf/de-sitter-susskind-ordering-metric.pdf)。延续已有署名，不作对外发表。本阶段不 commit/push，不创建 ZIP。

## 从已有四点核确定内积

沿用 complex Gaussian p-body DSSYK、标准 CAR、退火归一迹，β=μchem=0，ℏ=1。先取固定 λ=p²/N 的双标度极限，再令 h=4λ→0；q=e^(−h)，u=2Js t 固定。X 与 Rhat 是此前已归一的 d=2、e=3 导数 singlet。

Lin–Stanford 已给出领先 commutator 核及 block 权重。令 B_h(u,v)=−i[X(u),Rhat(v)]/√h，将此前得到的四时间协方差写成 D⁺_d⊗D⁺_e 中的矩阵元，直接得到

\[
\langle f,g\rangle_{\rm odd}=\langle f,\mathcal A g\rangle_0,
\quad \mathcal A=C_{\rm tot}-C_1-C_2,
\quad \mathcal A_k=k(k+2d+2e-1)+2de.
\]

这是已有领先核的算符实现，不把排序 primary 或标准张量积分解当作原创。对于 d=2、e=3，A_k=k(k+9)+12。证明只用相干态矩阵元与标准表示分解，不依赖原文由低阶推测的有限 q 全 k overlap 式 (355)。有限时间轨道的差商在 form norm 中产生后代；有限后代张成 form core，因此确定闭合的领先二次型。该论证不声称有限 h 算符的一致收敛。

令 P=p1+p2，x=p1/P，抽取 beta 权重后，H0=L²(dP w(x)dx)，w=504x³(1−x)⁵，且

\[
\mathcal A=-x(1-x)\partial_x^2-(4-10x)\partial_x+12,
\qquad
\|f\|_{\rm odd}^2=\int dP\,w(x)dx\,[x(1-x)|\partial_xf|^2+12|f|^2].
\]

Hodd=D(A^(1/2)) 配以上述范数，U=A^(1/2):Hodd→H0 是满射酉映射。总 Q=P 与 B=i(P∂P+1/2) 保持不变。归一奇扇区中的乘法顶点为 i√h A^(1/2)/2；k=0 给 i√(12λ)，与此前三点相符。它要求输入处于 form domain；hk²≈1 时不能沿用固定通道小 h 展开。

## 保持范数的平移及其更强约束

采用标准自由两粒子 light-front 质量平方 M(x)=m1²/x+m2²/(1−x)。P、x 无量纲，μ、m_j 为能量单位，P⁺=μP、P⁻_0=M/(μP)。普通 H0 中的 M 自伴正定，谱为 [(m1+m2)²,∞)。

在 Hodd 中直接乘 M 不对称，因为 [A,M] 的一阶导数项是 −2x(1−x)M′∂x。修正为 M_odd=A^(−1/2) M A^(1/2)，定义域通过 U 拉回，即得到自伴、正且同谱的算符；P⁻_odd=M_odd/(μP) 与 P⁺ 强对易，并满足 boost 伸缩关系。

精确诊断：m1=m2=1，固定 P fiber 中取 p0=1、p1=(10x−4)/√(24/11)，f=(p0+ip1)/√2。则

\[
\|e^{-isM}f\|_{\rm odd}^2=17+\sqrt{66}s+181s^2.
\]

这证明未修正的自由平移破坏奇扇区范数；不是有限矩阵截断效应。

更强地，若要求乘法顶点 T=A^(1/2) 是从普通两粒子空间到闭合奇扇区的平移交织映射，即 e^(−isN)T=T e^(−isM)，其中 N 在归一目标空间自伴，则取范数强制右边保持 A form。上式给反例。因此任意选 N 都无法满足这一特定要求。此结论以闭合目标扇区、自由源质量和交织要求为前提；微观乘法插入尚未被证明是 on-shell 过程。

F(A)>0 可以同时保持范数和交织，但质量函数仍任意。下一步必须从体模型确定空间插入或相互作用生成元，而非将 F 的任意选择称作介子谱。

## 证据与误差

exact.json 保存 50 条符号恒等式、6×6 质量矩阵和精确范数多项式。checks.json 保存五组独立相对坐标积分（55/85 位）、三组 h=.04/.02/.01 的直接/交叉 Gram（k≤4，85/115 位）和四组相位范数积分（70 位）。原始结果及输入哈希保留，未重跑覆盖。

85 位核积分残差最大 5.15e−85；Gram 提精度相对变化≤4.37e−79。但 h=.01,k=4 的有限 h 值为 54.77019259，对领先 64 仍偏低 14.42%。这与数值舍入不同，不省略。有限 q 检查不是全 k 证明。

编译首次遇到已有命令 AA 重名，续编发现两处轻微行溢出；原日志分别保留为 build-initial.log 和 build-overflow.log，修正命令名与措辞后双遍编译。正式阅读版及渲染核查记录见 local-evidence.json。
