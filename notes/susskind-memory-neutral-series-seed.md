# 下一计算入口：neutral分量的有限参数全次数分解

2026-11-07，two-component阶段后的继续种子。下式由旧Jacobi桥的代数重排得到；下一轮先独立核有限多项式版本，再用于联合涨落核。此处尚没有新的极限谱或数值输出，不能将界当作收敛定理。

沿 uniform-bound 记号，`Z_j=W_j(H^j)`、`E_j=Z_j−Q_j(H)`、`e^(isH)=Σ c_j(s)Q_j(H)`；`θ_j=μ(jp)`。Φ/C=μ(Rcal)，θ(Ncal)为coupling chaos乘子。每个Z_j纯degree-j，因此精确重排给

\[
\frac1C(1-P_h)\Phi(e^{isH})
=(1-P_h)\left[
(\mu(\mathcal R)-\theta(\mathcal N))e^{isH}
+\sum_{j\ge0}c_j(s)(\theta_j-\theta(\mathcal N))E_j
\right].
\]

理由：`θ(N)e = Σ c_j θ_j Q_j(H) + Σ c_j(θ_j−θ(N))E_j`，第一项逐样本是能量函数，pinching精确消去。原Jacobi估计和1/j!系数尾保证该重排在annealed L²中绝对收敛；放大p后 `p(p/N)=λ_N` 仍有可和上界。

这把neutral场的真正涨落分成两部分：rank–chaos multiplier差，以及Jacobi递推的低chaos误差。后者不能删掉，即使Z_j在最高chaos层已受控。设η_j=pE_j，则

\[
\eta_{j+1}=H\eta_j-Ca_j\eta_{j-1}-p r_j,
\quad
r_j=\W_{j-1}\sum_{v<j}[\Phi(H^v)-C\mu(vp)H^v]H^{j-1-v}.
\]

下一步可先求固定j的这些marked对象与charged transfer defect的joint moment：首个非零η_3来自r_2，必须恢复本轮B1/B2校准。对于第一项，固定chaos/固定deficit的斜率倾向 `4λq^k D`；但actual full exponential所需的D²统一可积性不能只由一个二矩上界推出。要给定量rank尾或高矩，再交换rank-cutoff与N极限。

charged侧仍须保留reverse source；尚未证明A/B全时间cross核或收敛。先核本接口与现有文献中的marked-chord/多probe方法有无直接可用公式，再安排新计算。不要重新证明已完成four-u归约，不重跑旧科学JSON。
