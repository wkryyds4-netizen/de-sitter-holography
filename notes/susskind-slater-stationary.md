# 从选择激发回到观察者态

2026-11-08。本轮接[全开弦传递](susskind-slater-transfer.md)，先核[文献](../literature/susskind-slater-stationary-sources.md)，再写[LaTeX证明](../manuscript/susskind-slater-stationary/main.tex)。固定 λ>0，q=e^(−4λ)∈(0,1)，半填充 annealed Slater 态，同原 Gaussian J，标准 CAR/overlapping p-body。h 按方差归一，时间与 β 均为其无量纲共轭量；物理 string H=(J_s/√λ)h。先双标度极限，再讨论下述长时间，不声称有限 N 长时间一致。

**主结果。** 在既有 h,M 的 cyclic GNS 空间中，h=J_q⊗I 有纯绝对连续谱，M compact。令 A=C*(1,h,M)，则

\[
A=C^*(h)+\mathcal K(\mathcal H),\qquad
A/\mathcal K(\mathcal H)\simeq C([-E_*,E_*]),\quad E_*=2/\sqrt{1-q}.
\]

证明先由Q=MM†−qM†M=q^N的幂范数逼近真空投影，再由cyclicity生成所有rank-one。h谱的每个非空开区间均有无限维投影，故没有非零连续f(h)是compact，商范数就是sup|f|；和空间闭合。这比仅知W*(h,M)=B(H)多给了真正的norm algebra。

B(H)上没有normal stationary state：其密度ρ若与h对易，则ρ的任一正特征值给非零有限维h不变子空间，产生h本征向量，矛盾。e^(−βh)普通trace对任意有限β发散。更强地，B(H)上连singular有限β KMS也不存在：有界inner h的KMS态经e^(βh)倾斜成为tracial state，而两个正交isometry排除B(H)的归一trace。不要将这个结论误报为“没有任何singular invariant state”。

任意A上的stationary state限制到compact ideal是trace-class positive functional，其density同样被迫为零。因此全部stationary state恰为ω_ν(f(h)+K)=∫f dν，ν为谱区间上任意概率测度；每个也是任意有限β的KMS态，因商动力学为恒等。所有含至少一个M或M†的有限有界word均被杀掉；GNS是乘法代数L∞(ν)，不是非交换静态物质理论。β不在这里选出唯一热分布。该分类只限当前A和当前flow。

**可构造的极限和有限时间状态。** 纯AC谱使e^(−ith)ψ弱趋零，故任何normal初态ρ在A上逐observable趋ω_(νρ)，νρ是其守恒能量分布。这是pointwise weak*收敛，不是normal态空间或B(H)上处处收敛。对选择/欧氏制备，νβ=(1−√q+√q E²)e^(−βE)μ_q/Z_M；能量倾斜保留，M关联消失。

以Ω为初态，ρ_T=(2T)^−1∫_[−T,T] e^(−ith)PΩe^(ith)dt是具体normal有限时间态，任意u的stationarity defect ≤min(2,|u|/T)。设R(t)=∫e^(itE)dμ_q，κ=π∫ρ_q(E)²dE，则

\[
\operatorname{Tr}(\rho_TP_\Omega)=T^{-1}\int_0^T|R(t)|^2dt\sim\kappa/T,
\quad \operatorname{Tr}\rho_T^2=T^{-1}\int_0^{2T}(1-t/(2T))|R(t)|^2dt\sim\kappa/T.
\]

有限资源的近似平稳化与有限rank信号流失可同时量化；这不是固定finite-N的热化定理。有限Jacobi exact pinching会留下Σw_j²平台，不能把有限cutoff恒定项当连续谱normal equilibrium。

**建设性接口。** 对任何不依赖固定O的联合U(N)-covariant算符族B_a(J)，每个可积word F的Gaussian平均在wedge^m不可约表示上为scalar，因此E〈O|F|O〉=Eτ_m(F)。这是旧H-only Schur事实的直接推广，包含适当的协变bilocal、时间移位和joint-annealed Gibbs夹心；不包含固定O的M。有限N协变子代数中的Slater限制态已经tracial/stationary，可先建立该子代数的GNS，再求与M的mixed moments。这里只是固定charge等式，不能跳过central/full-charge桥。

扩大joint observable family时，旧h,M的cyclic空间只是更大空间的reducing子空间。M在旧子空间上compact不保证在新增扇区上compact，故本轮结论不阻止构建更大的非平凡平稳可观测代数。这是下一主线。若候选generic算符始终留在当前norm A且其normal参考态限制stationary，long-time compact消失迫使其物理表示只剩commuting能量函数；要承载非平凡bilocal动力学就必须越过这个范围。

**检查。** 新独立driver三q、800/1100节点、36时间记录，最大差1.23e−15；θ密度两网格/两product截断、三直接时间积分，以及五小系统occupation projector精确计数。浮点诊断不证明无限维定理、不算区间证书。输出含相对渐近比和有符号求和条件数，所有参数/版本/输入hash保存。已成功的科学driver不重跑。下一步优先明确旧bilocal中最小非平凡协变生成元的同态/同J mixed接口，不再尝试用ρ∝e^(−βh)在当前B(H)内修补平稳态。无新Git或issue动作。
