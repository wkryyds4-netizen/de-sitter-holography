# 第一次记忆返回的守恒背景

2026-11-06，本地阶段。[六页PDF](../output/pdf/de-sitter-susskind-memory-return-background.pdf)、[LaTeX](../manuscript/susskind-memory-return-background/main.tex)、[来源](../literature/susskind-memory-return-background-sources.md)。沿用 Alan Day 与 Chengheng 的研究记录署名。未 commit/push/ZIP。

## 结果与规范

沿用含指标重叠的 complex Gaussian p-body H、标准 CAR、共同 J、β=μ=0 全归一 annealed trace、ℏ=1、string 时间。偶 N，λ_N=p²/N→λ>0 且在正紧集，Js>0 固定。q=e^(−4λ)，σ²=Js²/λ；固定任意有限时间窗，不交换小 λ 或长时间极限。

令 Λ(s)=R†exp(iKs)R，M(s)=R†exp(iQKQs)R，h=H_m，a=(1−p/m)h。P_h 是逐样本完整简并能量块 pinching。实际场 Fscr=(1−P_h)m exp(−ihs/2)M exp(−ihs/2)，普通项 X0 用 Λ 替换 M；m<p 块设零。

本轮证明固定 S 上：

\[
\sup_{|s|\le S}\|p\Lambda(s)-b_s(H_N)\|_{2,\rm ann}\to0,
\qquad
\sup_{|s|\le S}\|\mathsf F(s)-X_0(s)\|_{2,\rm ann}\to0.
\]

其中 b_s=4Js² P_q(exp(is·)) 是文献已有 Mehler 条件算符作用后的能量函数。有限 N 能谱可能出界，所以只把 Mehler 函数的外部参数截到极限谱区间 [−E*,E*]；后面的真实时间相位仍用实际能量 x。

## 从两个带标签的收缩到强极限

设 Fi=[H,ci]，Q_N(f)=(2p/N)Σ Fi†f(H)Fi。计算 Eτ Q_N(f)†Q_N(g) 的固定次数 Gaussian 配对。每条弦独立求和两个 p-集合；一个 force pair 的两端若未内部配对，它的共同 site 必须落在两条独立弦的集合交集 A∩B 中。两 tag 总计数至多 p|A∩B|，即加权 O(p)，而两个内部 pair 给 p²。共享 tag、共享弦、重合交集以及外部 H 收缩均不能解除约束。

固定字长时，以 Z 表所有集合两两交叠之和，CAR 相对迹权由 C exp(CZ) 控制。超几何阶乘矩给 E exp(tZ)、E|A∩B|exp(tZ) 统一有界。双内部 pair 中剔除 tag 碰到其他集合也只损 O(p)；各孤立 tag 去掉一个 trace=1/2 占据因子，故各乘2，剩下两条指定且不交叉的 bubble。

已有 BKS 条件期望在 q-Hermite 基上满足 E_H[B W_n(H) B]=q^n W_n(H)。依次条件到(H,C)再到H，得到两个 Mehler 像的乘积。辅助 B、C 只编码指定收缩，不能把实际 Fi 改成独立 probe。两条给定 bubble 在常数函数下贡献1；同一独立 probe 的四次矩则为2+q，后者多出的配对必须先由 tag 论证压掉。

令 c=4Js²，则双插入矩趋 c²〈P_q f,P_q g〉；已有单插入混合矩趋 c〈r,P_q f〉。再用 H 的矩收敛，展开平方即得固定多项式 Q_N(f)−c(P_q f)(H_N)→0 的 L² 强极限。标准 Mehler 接口不是本项目新公式；本轮新增的是实际共同 J 对象的双 tag 桥接及其应用，不作优先权判断。

## 有限时间、全电荷与乘积

旧 return-expansion 给 Fi 两侧 covariance≤v_N，v_N=4(Js²/λ_N)(p/N)(1+p/N)^(2p−1)，并给每个 sector Etr A^k≤(2k−1)!!(Nv_N)^k。旧 charge-sum 的 C_m≤(Js²/λ_N)e^(2p/N) 使 H 的 Gaussian 矩界统一。

全 Fock Schatten Hölder 给 ||Q_N(f)||₂≤(2p/N)Σ||Fi||₈²||f(H)||₄≤C||f(H)||₄。实轴 Taylor 余项：

\[
\|e^{isH}-T_{k-1}(s,H)\|_4
\le \frac{S^k}{k!}[(4k-1)!!M^{2k}]^{1/4}\to0.
\]

先固定 k 取 N 极限，再取 k→∞；不需要所有次数一致的 diagram 常数。有界连续目标 extension 与多项式目标之间由谱弱收敛、高矩一致可积性比较；有限 s 网格与一致余项给紧时间一致。

D_F=KV−Vh，R=D_F+(p/m)Vh，精确 pD_F†exp(iKs)D_F=(N/(2m))Q_N(exp(is·))|m。中央 m≥N/4 上 D_F 的四范数 O(p^−1/2)，(p/m)Vh 的四范数 O(p^−1)，故 R 替换误差乘 p 后 O(p^−1/2)。Q_N 的 sector 二范数统一有界，binomial 集中给 prefactor 误差 O(N^−1/2)。活跃下尾 m≥p 只有多项式增长，被 exp(−N/8) 压住。

由 G≤A/m 得 ||pΛ||₈统一有界，插值把 L² 收敛升为 L⁴，指数为1/3。这一步才允许在第一次返回的乘积里替换两个 pΛ。矩形 Schatten 统一按源维数 d_m 归一。

## 第一次返回的具体函数

对 s≥0，Δ_s={(v,w):v,w≥0,v+w≤s}、t=s−v−w。J1=∫Λ(v)exp(iat)Λ(w)。先构造保留 m、a 的精确守恒代理：

\[
\mathcal K_{N,m}(s)=\frac m{p^2}e^{-ihs}
\int_{\Delta_s}b_v(h)e^{iat}b_w(h)\,dv\,dw.
\]

m/p²≤1/λ_min，L⁴ 替换给 m exp(−ihs/2)J1 exp(−ihs/2)−K_N→0 的 L² 收敛。令 g_v(x)=exp(−ivx)b_v(x)/(4Js²)，合并相位后 K_N 只多 exp(−i(p/m)ht)。去掉该相位时 (m/p²)(p/m)=1/p 精确抵消，误差≤c²S³||H||₂/(6p)，无需提前作中央近似。

再以电荷集中取 m/p²→1/(2λ)，得到

\[
\kappa_s(x)=\frac{8J_s^4}{\lambda}
\int_{\Delta_s}g_v(x)g_w(x)\,dv\,dw,
\quad |\kappa_s(x)|\le\frac{4J_s^4}{\lambda}s^2.
\]

负时间按 κ_(−s)=conj κ_s；J1 的两次定向积分反序也给伴随关系。谱内小 s 系数为 4Js⁴s²/λ+i[8Js⁴/(3λ)](q−1)x s³+O(s⁴)。旧一次返回余项集体范数 O(1/p) 接上后，完整 m exp(−ihs/2)(M−Λ)exp(−ihs/2)→κ_s(H_N)。精确 pinching 消去它，得到本轮场等价。

## 此处尚不能省去的下一步

两个场的统一 annealed L² 有界性等价；本轮没有证明其中任何一个有界。已有粗范数只有 O(p)。因此 o(1) 场差不能直接推出 o(1) 绝对协方差差。只有补出统一 L² 界后，Cauchy–Schwarz 才给对任意外部中心时间 T 一致的协方差等价。此前两次返回六产品的 O(1/p) 绝对协方差归约仍无条件成立。

下一任务集中于 X0 的次领先共同 J 双 bubble 涨落、两个 memory 参数和精确能量去对角。固定次数的 Y/L/U2 谱可作低阶检查，不能代替有限时间余项控制。实际体介子内部核与有限 G 字典继续在后续。

## 保存的计算

`bubble-and-return-checks.json` 已运行一次，不覆盖：三个有理 q 各791个五 H 块字、总 H 次数≤8，共2373个精确恒等式。它验证指定 bubble 的标准条件期望，不验证有限 N tag 压低。

κ 的27个 float64 值：Js=1，λ=.25/.5/1，能量为谱边缘的0/.4/.85，s=.5/1/2；Jacobi谱阶64/96、单纯形Gauss阶20/32，最大阶数变化1.3255300427371534e−13。非区间认证误差。Python3.14.5、numpy2.4.6、scipy1.18.0；哈希及双遍编译方式见 local-evidence。最终六页95dpi逐页核版，首轮三处 overfull 日志保留。
