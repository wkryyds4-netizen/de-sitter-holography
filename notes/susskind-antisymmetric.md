# 左右弦差的固定阶极限与正生成元

2026-11-06，本地阶段。主证明见 [六页 LaTeX 稿](../manuscript/susskind-antisymmetric/main.tex) 与 [PDF](../output/pdf/de-sitter-susskind-antisymmetric.pdf)。署名延续 Alan Day、Chengheng。本轮无 Git 提交、推送或 ZIP。

## 物理推进

前一轮的 Q=p 由 boost 核构造出来。这一轮在一物质 chord Hilbert 空间里给出实际算符 Q_h，并证明其正性、固定后代空间的代数极限及初态谱测度收敛。因此候选不再只有标量两点拟合，而有与真实 Liouvillian 相连的弦表示来源。质量与能量标尺仍要靠独立物理信息确定。

态：完整 annealed normalized trace，β=μ_chem=0。先 complex DSSYK double scaling，固定 λ=p²/N，再 h=4λ→0；q=exp(−h)，物质标签 d=2。证明可用于任意固定 d>0。这里是 DS 后的一物质弦表示，不是有限 N 样本算符证明。

E_h=nbar，grade n 上 E_h=d+n；C=aL†−aR†，D=C†；K+=C/√h、K−=D/√h；B_h=(K++K−)/2，Π_h=(K−−K+)/(2i)，Q_h=E_h+Π_h。B_h=(Hs,L−Hs,R)/(4J_s)，η=4J_s t，全部新生成元无量纲。

## 一般固定阶证明

从 Lin–Stanford 的确切弦代数出发，δ=1−q，A=aL†aL+aR†aR≥0：

1. 每个 annihilation grade block 的范数平方≤[n]_q≤n，由 a a†=1+q a†a 归纳。
2. [D,C]=2(1−q^E)−δA。令 t_n=||C:H_n→H_(n+1)||²，则 t_n≤t_(n−1)+2(1−q^(d+n))，故 t_n≤h(n+1)(2d+n)。这是整个 grade 空间上的界，包含 Gram 退化方向。
3. w_n=K+^n Ω 的范数平方≤n!(2d)_n。又 [aL,K+]=(1−q^E−δaL†aL)/√h，在 grade k 上范数≤(d+2k)√h。逐次交换得 ||aL/R w_n||=O(√h)，进而 ||A w_n||=O(√h)。所有常数依赖固定 n,d。
4. 因而 ([K−,K+]−2(d+n))w_n=O(√h)。再展开 [K−,K+^n]，得 K−w_n=n(2d+n−1)w_(n−1)+O(√h)。取内积递推即得 ||w_n||²→n!(2d)_n。
5. 归一后代 e_n(h) 正交，升/降系数分别趋 √((n+1)(2d+n)) 与 √(n(2d+n−1))；降算符在此族以外的分量也受范数控制。所有固定多项式 word 在有限后代 span 上收敛。

这解决了最初真正的难点：仅把文献的 c_k(n) 展为 h n(2d+n−1)+O(h²)，不足以控制 h^(−n/2) 归一化后几乎 null 的向量。现在使用整个 grade 的算符范数界，不需要假定坐标误差等于 Hilbert 范数误差。

## 正性、谱与域

标准 Laguerre 正交函数 g_n=√(n!/(2d)_n) ψ_d L_n^(2d−1)(2p)，φ_n=i^n g_n，ψ_d=2^d p^(d−1/2)e^(−p)/√Γ(2d)。在此基中 B=i(p∂p+1/2)，Q=p，E 的本征值为 d+n。B、Q 的 Jacobi 系数线性增长，Carleman 给有限基 span 上本质自伴；Q 的闭包为正乘法算符。

对有限 grade 向量 f=Σf_n，r_n=||f_n||，c_n=√((n+1)(2d+n))：

〈f,Q_h f〉≥Σ(d+n)r_n²−Σc_n r_n r_(n+1)=∫p|Σr_n g_n(p)|²dp≥0。

有限 h 下 aL/R 有界，Π_h 有界，故 Q_h 自伴于 D(E_h)，有限 grade 是 core。正性不只是渐近性质。

协变缺陷精确为 [B_h,Q_h]−iQ_h=i{(1−q^E)/h−E−δA/(2h)}，在固定后代上→0。它是双侧 chord 算符，且有限 h 不满足精确伸缩协变，因此不与前轮有限迹单侧 obstruction 冲突。

初态 Q_h 的所有固定矩趋 (2d)_k/2^k。正性和矩界给 tightness/一致可积性，Gamma law 的矩确定性给弱谱收敛：ν(dp)=2^(2d)p^(2d−1)e^(−2p)dp/Γ(2d)，特征函数 (1+is/2)^(−2d)。有限 h 精确均值 d，方差 (1−q^d)/(2h)。没有由此推出 Q_h 的逆收敛；p=0 奇异端点须另查。

## 文献与验证

二项式后代不是我们发现的：Ambrosini 等 2412.15318 §3.1 (81) 已给出，(91) 给精确降系数。他们亦指出有限 q 从第四级起 binomial 不再是精确 Krylov 基。其 fixed tilde-q heavy limit、fixed hn 大级数极限不等于本轮 fixed d,n。一般证明是本轮从已有代数整理的推导，未宣称文献优先权。

已有 probe.json 的 d2、n1..6 精确 Gram 多项式不重跑、不覆盖。新 validation.json 使用独立数值递推：d=.7,2；h=.1,.02,.005；80/110 位；六阶后代、五个整级 creation block 范数、八阶 Q_h 矩。最大相对精度变化 2.04e−69 以下；这是经验算术核对，不是 interval certification。另有 n0..8 的18个精确 Laguerre 微分/递推恒等式。

d2、h=.005、n6 归一范数比 .93482790596，第八谱矩相对极限偏差 .01268639139。没有把精度稳定写成小 h 近似已经精确。精确 probe 见到各受测 n 的归一 lowering residual 为 O(h)，正文只采用已对所有固定 n 证明的 O(√h) 界。

## 下一步

精读实际两物质表示与 coproduct，在两物质内积上构造相应 Q_h，检查相对后代、交换顺序与相互作用。不能直接采用独立单粒子动量 primitive 相加。与已得三 channel 首阶四点对接，再研究 P−、能量标尺和介子质量。暂不把任意 m 的 P−=m²/(μQ) 当微观动力学。
