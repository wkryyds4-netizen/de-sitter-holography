# 固定半填充的普通 charged 源与逐阶强桥

2026-11-08。接 higher-Wick 阶段。正式稿 manuscript/susskind-central-ordinary/main.tex，六页阅读版 output/pdf/de-sitter-susskind-central-ordinary.pdf。本轮只本地，已有 Git b6b2a45/issue6711589107 不重复。

## 新闭合点

对于任意固定 ordinary D 标记区（允许内部 Gaussian pair）和无内部pair的L子区，精确投影到m=N/2之后，源指数为

\[
K_{\pi,\mathrm{cen}}(z)=K_{\pi,\mathrm{full}}(z)-2I(z)^2,
\qquad I(z)=\sum_a i_a z_a.
\]

i_a为第a个D区的内部弦数；L区取0。因此DD二源系数扣去4i_ai_b；single源不改，DL/LL不改。该结论从有限N Fourier charge投影及整个积分区间的一致可积界推出，不从fulltrace联合矩直接令X=0。

standard CAR、overlapping Hermitian complex Gaussian p-body、各charge同一原J、annealed normalized half-filled sector；lambda_N=p²/N→固定lambda>0，h_N=H_m/√C_m。D/L先作用fullFock再限sector，C_m是目标sector确定常数。h时间无量纲、hbar1。所有次数/mark数先固定，不含增长degree、增长time或联合小lambda极限。

## Fourier证明的要点

内部弦的两个singleton在rho=e^(iθ)/(1+e^(iθ))下给(1−z)(1+e^(iθ)z)，p个站点后为exp_nil[p(e^(iθ)−1)z]。open chord的odd singleton与pd源精确抵消。必须保留wrong-occupancy source项：D(1−n)=−n，所以原来零的占据贡献可变成纯z。

有k弦时，标量singleton和全charge迹合成cos(θ/2)^(N−2kp)。目标sector归一化的单弦标量A_N=4^(−p)binom(N,p)²/[binom(m,p)binom(m+p,p)]→e^(−lambda)，有界。固定t=pθ时，旧collision展开推广任意固定mark数，给

\[
F_{\pi,N}(t/p,z)\to e^{-t^2/(8\lambda)}
\exp_{nil}[\lambda K_\pi(z)+itI(z)].
\]

局部一致不够，额外控制如下。固定小δ<π/2，对|θ|≤δ，source singleton分母的零阶均非零。R为有限个p-set的总pair overlap，missing singleton memberships≤2R；消去extensive p源之后，碰撞/缺额系数≤K(1+R)^r A^R。hypergeometric指数矩与Holder给R的所有固定指数矩有界。于是每个source coefficient≤K(1+p|θ|)^r cos(θ/2)^(N−2kp)，在t变量由K(1+|t|)^r exp(−t²/(16lambda+))支配。

对δ≤|θ|≤π不用complex singleton分母。每配置support U≤2kp，support外D_iI=L_iI=0，精确保留cos(θ/2)^(N−U)；直接展开loss/number最多p^r，归一标签计数最多e^(Kp)。故远区界Kp^r e^(Kp)cos(θ/2)^(N−2kp)，乘p后积分仍指数小。两段共同给真正Fourier dominated convergence。

将分子积分除以未标记binomial Fourier积分，Gaussian变量t的方差为4lambda，所以E_t exp(itI)=exp(−2lambda I²)。完整保留nilpotent多源系数，得到central source。Wick/chaos filters只删有限图，不影响该证明。这里常数仍依赖固定图，不是degree-uniform估计。

## ordinary强桥与实际投影

定义目标sector归一G_n=C_m^(−n/2)[pNcal−D−lambda_N Ncal(Ncal−1)]H^n|m。monic q-Hermite展开h^n=Σa_ndH_d，m_d=q∂qlog[d]_q!，中央守恒反项

\[
g_n^{cen}(h)=\sum_d a_{nd}[-4\lambda m_d-\lambda(n-d)]H_d(h).
\]

有限N使用q_N、lambda_N。内部matching切图后single系数−4lambda inv−2lambda i、connected只保留桥弦项；减Mahonian平均正得旧pureJ_d源。RR、RJ、JJ同原J的Gram同时一致，因此

\[
\|G_{n,N}-g_n^{cen}(h_N)-\sum_d a_{nd}(q_N)J_{d,N}\|_{r,m}\to0,
\]

任意固定n和有限r。先二范数差，任意固定更多source给高矩有界，再插值升r；Holder允许固定左右H幂。不是两个边缘律相同就宣布same-copy强收敛。

g与h精确对易，上一轮P_h J_d在fixedS_r趋0，所以actual (1−P_h)G_n→Σa_ndJ_d。完整简并能量块保留，先finiteN exactpinch后DS。fixed-polynomial orbitdensity把误差传到已有固定多项式制备向量、selected单协变块，固定h-time phases依旧先用central Gaussian Taylor尾控制。所有固定ordinary truncations因此接上selected态，限制态仍由charged塔f(h)=1−√q+√q h²给出。

二阶精确校准：m=N/2，C_(m−1)/C_m=(m−p)(m+p+1)/(m(m+1))，归一零chaos b_Np=−p(p+1)/(m+1)。整个R2−J2精确等于(b_Np+2lambda_N)I。N=4p²时误差=(1−2p)/(4p²+2)→0。两个分别internal G2的fulltrace图为4lambda²+4lambda，central为4lambda²；新规则确实删除电荷方差而保留固定charge背景。

## 独立检查及覆盖缺口记录

三个新输出均保存后未重跑。

- exact-checks：2186 Gaussian Fourier nilpotent恒等式、20energy regressions、32Gram layouts/96equalities、2个full/central G2非零差、5finiteN scalar bridges；另存p16/32/64/128/256、N4p²的30个exact singleton条件系数和真实有限偏差。SHA54286ec6f07967d3c60770c44ffd93799dd4eac26cb346dbdcab49165280f45c。
- 首driver的local循环因递归tail条件写错只生成单标记，打印source comparisons0；该输出不作为local证据，原样保留。另起check_central_ordinary_local.py，显式生成区间并强制2/3/4source计数非零，核3496/459/1265案例，共5220；1380含internalpairs、4010负nonconstant系数。SHAe7321618b31f4bffe489a4f2eb61f379ee56e443dcfc84dbb912c372195b7add。
- Fourier数值独立比较exact binomial系数：10组finite-support case，60/90digits，最大精度差4.58e−59，analytic |t|>16尾最大1.01e−41。保存absolute integral/cancellation ratio；积分离散/舍入没有interval认证。不是完整有限N模型仿真。SHAeceb91cd3d59939dad487a968f5a8370ade15aaff8d7690bbe317662a27961f3。

六页PDF双遍成功，90dpi全部目检。补A_N显式归一及support措辞后再双遍，只有page2/3改变并重检，其余4页render hash不变。封存manifest绑定当前、父科学文件与编译渲染。

## 下一步

central ordinary固定degree强桥已解决，剩实际全次数tail。旧fullFock Riesz不能直接限制sector，本轮Fourier domination常数随图数/degree变化，不能将其用作统一尾；fixed-degree preparation常数也不够。优先从旧charged-proof的具体cross-product估计找sector替代，或证明真正可行的双截断论证。不能以逐阶收敛直接交换无限和/N极限。随后完整selected bilocal、bulk空间/质量、finiteG。旧fulltrace完整四点不重列缺口。
