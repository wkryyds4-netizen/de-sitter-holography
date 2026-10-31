# 任意紧转移、物理谱尾与两次创建的 Bose 极限

2026-10-31，Alan Day / Chengheng。自查研究稿，未作原创性声明。[七页 LaTeX](../manuscript/ds3-compact-global/main.tex)、[PDF](../output/pdf/de-sitter-ds3-compact-global.pdf)、[来源](../literature/stage4-ds3-compact-global-sources.md)。承接 [physical-local](stage4-ds3-physical-local.md)；旧科学文件保持原字节。

## 模型与本轮结论

曲率半径 ℓ=1、Newton G=0、proper SO⁺(3,1)，两个可区分 BD 标量观察者，标签下界固定为 s>0，KG 物质 p,q>0，单位局部耦合及半周期调节创建块。静能量仍为 h(x)=√(1+x²)。初态条件 c=M 由归一窄质量包极限定义。b=M+u、a=M+u+v、d=M+v+r；条件测度 da db dd=du dv dr，能量变量另有 ∏h(x)/x，不能重复插入谱密度。

本轮将旧小窗口核定理扩展到任意固定紧转移集，再对实际两条创建链分别截取。得到先 M→∞、后截取半径 K→∞ 的物理配对极限。在这个指定模型的两次创建扇区，同种物质两种次序的 Bose 范数及相互重叠均趋于 Lp(0)²+|Lp(t)|²，次序差平方范数 DB→0，紧实时间一致。它不包含四个自伴场的中间真空收缩，不是全 Fock 或有限 G 定理。

## 为什么能扩大紧窗口

令 Aⱼ=1/2+iαⱼ，α=(v+q,p−u,−u−p,v−q)/2，z=−ir/2，|r|≤R。对任意实 α，cosh 无限乘积给

\[
\frac{|(1/2+i\alpha)_n|}{n!}
\le\sqrt{\frac{\cosh(\pi\alpha)}{n+1}}.
\]

其中中央二项式因子的递推由 (n+1)³−(n+1/2)²(n+2)=(3n+2)/4>0 保证。原小参数不等式只需增加常数 Cser=cosh²(πA*)，A*=max|αⱼ|+R/2；不再要求 A*≤1/2。

取 B=max(|u|,|v|)+R，m=M−B>0。两个小参数级数项满足

\[
|f_n|,|h_n|\le\frac{C_{\rm ser}}{(n+1)^2}
\prod_{j=1}^n\frac{j^2}{j^2+m^2},
\quad |f'_n|\le6H_n|f_n|,\quad |h'_n|\le H_n|h_n|.
\]

导数估计只用实部 1/2 与 1，虚部扩大不改变可和性。令 𝔞=zyA、𝔟=−zyB，两者在零点为1，F(0)=H(0)，先合并 Y=(𝔞F−𝔟H)/z。新的 C0、Q、RF、RY 和 TY(N) 在稿式(3)–(5)及代码 constants/merged_y 中明确给出。系数界变大，但阶仍是 O_K(logM/M²)。

完整伴随 G 的逐项 sinh 乘积界及十次首一多项式尾证书不依赖旧 A*≤1/2 条件；每个多项式系数对 M 至多八次，整个伴随贡献仍为 O_K(M⁸e^(−πM/2))。保留完整 G，不以首项代替级数。

小 Gamma 参数的实部为1/2或1，大参数留在远离负实轴的固定扇形。移位 Stirling/digamma 的 C¹ 余项因而在任意固定紧移位集合一致。去掉 M^±z 后的两支误差在 r=0 共享 ε0=O(M⁻²)；分离 ε0(e^(−irlogM)−1) 后没有虚假的 1/r 发散。因此

\[
U_{00}=\frac{\sin(r\log M-\delta(r))}{\pi r}+O_{\mathcal K}(M^{-1}),
\quad
\delta=\Im\left[2\log\Gamma(1-ir/2)+\sum_j\log\frac{\Gamma(A_j+ir/2)}{\Gamma(A_j)}\right].
\]

相位连续从零取值；r=0 用连续值。Dirichlet 论证给任意固定紧支撑、C² 范数一致有界测试族的 dr 分布极限，系数为1。紧集常数随 K 增长，不能任意令 K 与 M 同时增长。

## 实际链截取与全物理配对

原始 E/W/ρ 约分后 s 链的平方范数密度为 g(a)/g(M)·hatR_q(a,b)hatR_p(b,M)。旧真实权重定义给 uniform M≥s 的上界

\[
\frac{B_s}{16\pi^2}\left(1+\frac{|u|+|v|}{s}\right)e^{-\pi(|u|+|v|)},
\qquad B_s=2e^{\pi/(\sqrt{1+s^2}+s)}.
\]

选平滑0≤ηK≤1，|x|≤K时1，|x|≥2K时0。两条链分别乘 ηK(a−M)ηK(b−M) 与 ηK(a−M)ηK(d−M)，交叉积分中的乘子必须是 χK=ηK(u+v)²ηK(u)ηK(v+r)。它支撑于 |u|≤2K、|v|≤4K、|r|≤6K；这是确切的链乘子，不能任意指定 χ 再假称为投影。

被删平方范数≤εK，全部平方范数≤Nbar：

\[
\epsilon_K=\frac{B_se^{-\pi K}}{4\pi^2}
\left[\frac K\pi+\frac1{\pi^2}+\frac1s\left(\frac{K^2}\pi+\frac{2K}{\pi^2}+\frac2{\pi^3}\right)\right],
\quad \overline N_s=\frac{B_s}{4\pi^4}+\frac{B_s}{2s\pi^5}.
\]

理由是删除区域蕴含 |u|+|v|>K；四象限积分化为径向 ∫K∞x(1+x/s)e^(−πx)dx。完整重耦合为 Hilbert 空间酉变换，标量压缩为 contraction，由 Cauchy 得 |XM−XM^K|≤2√(εKNbar)，对所有实时间有效。不声称逐点 |U|≤1，也不要求 |AMU| 的全谱绝对可积；全配对由 Hilbert 极限定义。

固定 K 后，用真实权重的紧 C² 展开，峰面极限为 kp(u)kq(v)，其中 kp(w)=[8π(coshπw+coshπp)]⁻¹。于是 XM^K→IK，IK→Lp(0)Lq(0)。用“尾 + 固定K误差 + 极限尾”三角不等式，先取 M 的 limsup 再 K→∞，得到不同 species 的完整交叉基准。Lp(t)=sin(pt)/(4πsinhπp sinh t)，Lp(0)=p/(4πsinhπp)。

同种物质记 N、Y、E1、E2、Z 如旧 Bose 稿。自交换相位 exp[±it(hd−hb)] 在峰面趋 exp[±it(v−u)]，crossed 相位趋1。所以 N,Z→Lp(0)²，Y,E1,E2→|Lp(t)|²。√2P+ 给 Bose 范数 N+Ei、重叠 Y+Z、DB=2N+E1+E2−2Re(Y+Z)→0。没有额外2!。

此外，对 [s,∞) 上归一且质量概率逃向无穷的任意制备均成立：每个固定C，πM([s,C])→0。以统一范数界分开 c≤C 和 c>C 即可，不要求质量矩。固定不移动初态包、移动阈值和增长时间均不在本结论内。

## 数值证据、失败和修复

[初始数据](../results/stage4/ds3-compact-global/study.json)：两组 (u,v,p,q)=(1.4,−.8,1.2,1.7)、(−.9,1.1,2.2,.6)，M16/64/256、r−.6/0/.6，90位、merged N160。18个完整核点均满足有限 M 的 Gamma-exact 界。没有新的全物理积分数值求积；有限点不能认证全积分极限。

6个120位未合并参考中 index16 失败：M256、第二组、r0，稳定值约1.88534345914，参考约1.87011200511，绝对差.0152314540，相对差.00807887。原尾约4.79e−134，未含舍入。[诊断](../results/stage4/ds3-compact-global/collision-diagnosis.json)复现了相同失败：90位解析再120位序列化各标签，使精确有理重合残差 a+c−b−d=10⁻¹¹⁷，绕过碰撞分支；两项相消condition约3.32e110。同字符串210位、精确有理标签120/170位均对原稳定值差约4.93e−88。原数据、初始failed状态及容差1e−65全部保持。诊断记录完整复值、精度、条件数和依赖哈希，非区间认证。

[独立审计](../results/stage4/ds3-compact-global/audit.json)初版442项中441通过，失败来自 SymPy 因式式与展开式的结构 == 比较。[修复审计](../results/stage4/ds3-compact-global/audit-refined.json)以精确 expand(lhs−rhs)==0 核相同恒等式，444项通过，数据和阈值未变。检查包括 Decimal 常数、完整复数组装、多项式导数正性、失败参考重放、谱尾原函数、20个有限Pochhammer核对及输入哈希。它不是定理证明器或总舍入认证。

## 后续

下一步先核完整 Fock 共同域及四个自伴场中的真空收缩，把两次创建结论接到明确算符有序积；再研究时间增长及有限 G。当前给迭代极限，没有优化全局 M 速率。双方署名 PDF 已实际双遍编译、桌面编译成功及七页逐页渲染核查。
