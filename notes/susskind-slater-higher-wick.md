# 高阶 Wick 后代的单选态与逐阶精确去静态

2026-11-08。接 selected-state 阶段与集中发布 b6b2a45。正式六页稿见 manuscript/susskind-slater-higher-wick/main.tex，阅读版 output/pdf/de-sitter-susskind-slater-higher-wick.pdf。本轮只本地，不新增 Git/issue。

## 本轮闭合的接口

单选态的旧 rho1 不只覆盖 h,u,x：对任意固定长度、互不重叠的 whole-Wick D 标记与 L 标记子块，仍有

\[
\sigma_1(B)=[1-s(1+\lambda)]\tau(B)+s\tau((h-W)B(h-W)),
\quad s=\sqrt q=e^{-2\lambda}.
\]

这里是多项式矩恒等式，不宣称所有高阶生成元有界。允许普通外部 H、较大的 whole-Wick 禁内部配对区，以及整体 chaos 投影；不允许带内部 Gaussian pair 的标记区，也不涵盖相互重叠的多次 loss。常数中心化与 Hermite 回归由线性展开处理。

模型为 standard CAR、overlapping Hermitian complex Gaussian p-body、同原 J、annealed 半填充固定 Slater，lambda_N=p²/N→固定正lambda。D/L 先作用 fullFock，再限 m=N/2；归一用目标 sector 确定常数 C_m，不能用 C_Q 替换。h=H_m/√C_m，W=(p−L−lambda_N)H|m/√C_m，h 时间无量纲、hbar1。

## 为何能从有限表到任意次数

两标签碰撞至多看四个端点，其他槽位在该 site 都是 I。删去 spectator 后，disjoint marked interval 仍为连续区间；空mark满足 D_i I=L_i I=0；singleton分母及JW sign不变。触极端的F仍触极端。因此两端点own-edge与四端点cross-edge的有限穷举覆盖任意原始长度。与旧固定Slater normalized-defect domination结合：固定layout只增加有限系数界B，activity仍为B m(p/m)^|S|，先控制尾再删triple；禁配对和chaos条件只删图。

own-edge全部G的source dependent项为 ell1+ell2+ell1 ell2，ell只在L端点取对应z；F的项为0。boundary selected F至少有一个端点在全部内部mark之外，故只减内部L线性源。不同edge的cross所有源系数不变，且没有inner-F nesting。因此 K_F=K_G−Σ_(L marks) n_a z_a，n_a是该mark接边界端点数。

另一侧加入allG两个边界L singleton a,b。除边界自身线性a+b后，分开接内端点时只剩 a ell0+b ellL；直接外配只剩ab。用h−W对应(1−∂a)(1−∂b)，得到(1−lambda ell0)(1−lambda ellL)，正匹配F源；外配给1+lambda，真实selected外配给1，两条mixed边归一带gamma²=s，合并即上述density。

248个局部layout中，241个cross-source恒等式、489个own-source恒等式、241个边界W cross零式、489个边界W own式全部exact Fraction通过，保留所有nilpotent masks。有限表是计算辅助的局部代数证书，一般次数结论另依赖上述压缩证明，不来自对低阶矩外推。

## charged tower 与实际 pinching

无L内部mark时ell0=ellL=0，所以整个固定次数charged塔有

\[
\sigma_1(B)=\tau([1-s+s h^2]B).
\]

取 monic q-Hermite H_d，m_d=q∂q log[d]_q!，

\[
J_{d,N}=C_m^{-d/2}(dp-D)\mathcal W_d(H^d)|_m
-\lambda_N[d(d-1)-4m_d(q_N)]H_d^{(q_N)}(h_N).
\]

d2恢复旧j_N。全部G的source逐图等于旧fulltrace源：local ratio=2Σ_oF_o，两者background都−2k²。Schur covariance再直接给central trace全部固定矩；没有从全charge平均反推central。

于是旧 memory-joint-regression-proof §§5–7 对J_d的固定degree两能量谱结果可逐项移植。外部normalized Hermite矩M_ab满足有限带宽 |a−b|≤d+e、polynomial(a+b)×q^min(a,b)界，谱有L²(mu_q²)密度、紧能量支撑、对角线零质量。central高矩、Portmanteau和pinching contraction给 ||P_h J_d,N||_(r,m)→0，任意固定有限r。随后旧fixed-polynomial orbitdensity把此误差传到所有固定多项式制备向量；single-selected协变块含固定多个actual a_d=(1−P_h)J_d和固定h-time也可替换。严格先finiteN exactpinch，再DS。

以E为input/column、E'为output/row，

\[
\lim_N\varphi_N(M a_d(t)a_e M^\dagger)
=\int [1-s+sE^2]e^{it(E-E')}d\nu_{de}(E,E').
\]

这是原有charged谱的选态重加权。D-only限制态平稳；不声称每个J_d边缘谱等于trace（此前只对曲率证明了该更强结论）。没有新增dS温度或bulk质量解释。

## 全时间接口的剩余位置

已可直接使用逐sector的pure-Wick tagged-deficit、coupling OU projection、normalized CP covariance path界；其固定r包络为poly(n)L_r^n√n!。普通G_n=[pNcal−D−lambda_N Ncal(Ncal−1)]H^n有mark内部pair，不在本源定理内。旧ordinary全次数尾的Clifford Riesz只在fullFock迹；fixed-degree preparation常数随degree增长，不能直接求和。

下一步优先补central ordinary charged强桥及uniform actual尾，然后selected完整bilocal重求和，再空间/质量生成元和finiteG。可以先独立重建含内部pair的central source，找共同守恒charge背景的固定sector替代；不要将fulltrace charge-tilt结果简单设X=0。已完成fulltrace完整四点继续保留。

## 计算与封存

pure helper code/slater_higher_wick_sources.py 提供validate/allowed、nilpotent exponential、collision_part、moment。moment的mode=trace/selected/chi；marks与wick_regions/chaos_regions明确分开，shifts是实际标量中心，不自动按旧W/J猜。

新driver首次完整通过，科学结果 results/susskind-slater-higher-wick/exact-checks.json SHA256为27b3557c6ce6c68605eb1a89ecbb7e23ebe7f11d00fc40e87f84f8425f845143。除局部表，还有两有理gamma/lambda点的30个高阶layout density比较（最高8内部槽）及18个旧H/W/J回归比较。有理参数不在物理曲线上，只验证代数。成功driver不得重跑覆盖。

六页PDF已双遍编译，90dpi逐页目检；首轮括号命令拼写错误保存diagnostics，修正后无undefined/overfull。编译命令 code/build_slater_higher_wick.ps1。依赖、父科学证据、最终PDF/LaTeX及渲染哈希见本阶段local-evidence.json。
