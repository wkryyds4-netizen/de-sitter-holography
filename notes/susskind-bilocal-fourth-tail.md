# 完整有限时间双局域四点：统一高范数尾与实际算符传递

2026-11-07，本地阶段稿。[十一页 PDF](../output/pdf/de-sitter-susskind-bilocal-fourth-tail.pdf)、[LaTeX](../manuscript/susskind-bilocal-fourth-tail/main.tex)。本阶段完成指定 complex DSSYK 模型内完整有限相对时间的共同四点构造。下一步转向态、极化与微观体生成元。

## 结论与适用范围

标准 CAR、允许 I/J 重叠的复 Hermitian Gaussian p-body Hamiltonian，各电荷扇区共用原始 J；态为 beta=mu=0 的归一 annealed full-Fock trace。偶 N、p→∞、lambda_N=p²/N→固定 lambda>0，lambda_N 位于正紧区间，Js 固定，hbar=1，string 时间。所有范数阶和实时间窗固定。

本次证明实际两分量表示的带电普通系数 G_n、K_n 及中性系数 T_j 满足有限 N 一致包络

\[
 \|X_{n,N}\|_r\le K_r(1+n)^{b_r}L_r^n\sqrt{n!}.
\]

相应时间级数尾被 sum_(n>K) polynomial(n)L_(r,S)^n/sqrt(n!) 控制，先在有限 N 建界，再交换截断与 DS 极限。完整 compressed-bath memory 用已有并已升阶的强近似接入；没有另称其自身 Taylor 系数已逐项证明同一界。

## 证明链与可复查入口

1. [标记图证明](susskind-bilocal-fourth-tail-diagram-draft.md)：受限 p-set covariance 的精确占据公式，配对端点的 row/column ancilla contraction，保留交叉和 CAR 符号。pure Wick 内每个标记站点约束两个独立标签，N(p/N)²=lambda_N。重复 D/L 缺额通过 graded body product rule 控制；给任意固定标记数的全次数阶乘包络。
2. [实际中性证明](susskind-bilocal-fourth-tail-operator-draft.md)：OU/Chebyshev 同阶投影估计；Phi_m/C_m 精确为 p 次 S 后 p 次 T 的 normalized CP 通道。逐步 telescoping 给实际 A_v=p(Phi−Ctheta_v)Z_v 的包络；展开受迫 Jacobi 递推给 eta_j，最终 T_j=A_j+(Ctheta_j−Phi)eta_j。这避免将单系数点值界当作矩阵和的范数界。
3. [普通带电证明](susskind-bilocal-fourth-tail-charged-proof.md)及[独立审计](susskind-bilocal-charged-tail-audit.md)：精确 D/L 与 normalized CP 的交换式，保留 charge 系数差分。Clifford Riesz 只用于 full-Fock trace，控制乘积 cross/p；因此实际 eta_j、Deta_j/p、Leta_j/p 同时受控。逆 Jacobi 展开保留五项，获得 ordinary G_n=[p Ncal−D−lambda_N Ncal(Ncal−1)]H^n 的有限 N 包络。
4. [完整四点传递](susskind-bilocal-full-four-point-transfer.md)：普通带电源不直接满足 pure 多源定理，使用旧同样本 R_n→sum_d a_nd J_d 的 L² 强桥与新高矩插值。中性后代的整个 composite Wick 区禁止内部配对，标记子区符合原多源规则。先固定次数传递 exact pinching，再解除截断，再接旧 S4 返回误差 O(p^-1/7)。

关键文献输入是 Ben Efraim–Lust-Piquard 的 Theorem 6.8(1),(3) 与 Lemma 5.6，[原文与规范核对](../literature/susskind-bilocal-clifford-riesz-sources.md)明确 N_Clifford=2L。其它矩阵 Gaussian/Clifford 工具及不能直接使用的推断见[文献备忘](../literature/susskind-bilocal-fourth-tail-sources.md)。已有定理是输入，CP 路径及本模型的重组是本地推导，不声称完成全领域原创性排查。

## 四点如何计算

固定次数先把全部后代展开为真正的共同源配对图；保留真实 Wick 禁配对区、chaos 筛选、外部普通 H 幂与所有反项。每张 2ell 槽图的有限源系数来自

\[
 (Ce^{-\lambda})^{\ell}
 \exp\{\lambda\sum_{u<v}(R_{uv}(z)-1)\}.
\]

R_uv 保留源依赖 singleton 归一和费米符号。一般后代不能套用仅 H/W1/J2 的三 species mixed-q 规则；两能量谱本身也不提供四点数据。

令 R_n^infty=sum_d a_nd J_d、L_n=R_(n+2)^infty−{H,R_(n+1)^infty}+H R_n^infty H，N_j 为中性回归残量，c_j 为 monic q-Hermite 指数系数。共同极限为

\[
 \mathcal F(s)=e^{-iHs/2}\left[
 \sum_{n\ge0}\frac{(is)^n}{n!}L_n+2C\sum_{j\ge1}c_j(s)N_j
 \right]e^{-iHs/2},
\]
\[
 \mathcal B(t)=-\int_0^t ds\int_{-(t-s)/2}^{(t-s)/2}
 e^{iHu}\mathcal F(s)e^{-iHu}du,
 \qquad \mathcal B(-t)=\mathcal B(t)^\dagger.
\]

上述级数在每个固定 S_r 和时间紧窗绝对可和；滤波范数至多 S²/2。四点按给定算符顺序，取四份截断场与中心时间共轭的共同迹，再去截断。四因子 Hölder 及实际 S4 返回界完成 actual-operator 传递。N 极限只在固定/紧中心时间成立，不声称全实轴或增长时间一致。

## 新科学检查与保存

- tagged-covariance-checks.json：42,336 covariance 比较，564 networks（452 crossed），94,743 label assignments，565 ancilla permutations；1,776 exact rational sector bounds，113 非零 crossed CAR-sandwich cases。首跑 1.83 秒。
- coherence-channel-checks.json：340 coherence identities、265 CP-path、75 inactive、2,609 kernel entries、1,305 site bounds；132 非常数输入及14非平凡 CAR gauges。首跑 0.091 秒。
- charged-reassembly-checks.json：348 neutral matrix units，3,312 matrix identities、9,949 entry comparisons、84 polynomial identities、16非零 charge-difference fixtures。首跑 0.528 秒；单实 Gaussian 坐标仅核代数，不充当完整 ensemble 的矩检查。

三份新输出均独占创建，未重跑旧科学 main。参数、版本、输入/输出哈希、三份父证据、双遍编译和最终十一页渲染由[本地证据](../results/susskind-bilocal-fourth-tail/local-evidence.json)绑定。封存器只复核保存的标量/计数和文件来源，不把 provenance 检查称为解析证明。

构建：`powershell -NoProfile -ExecutionPolicy Bypass -File code/build_bilocal_fourth_tail.ps1`；95 dpi 渲染，最终十一页逐页目检。初次缺宏和两处行溢出的编译日志保留在本阶段 diagnostics，最终没有 undefined/overfull。当前仅本地，无 Git/外部回复。

## 接下来的物理问题

完整模型四点已成为可使用的输入，接下来先重读 Susskind/Itakura/FLZ 的态与 canonical bilocal 构造，明确各自 Hilbert 空间、真空和正规序。检验当前 tracial 集体场可采用何种极化，以及微观 H 能产生哪些实际空间/内部生成元；尤其不能用中心能差或 memory 频率直接命名为体质量。有限 G 引力散射和完整 de Sitter 字典仍待建立。Issue #1 属 KS 支线，其旧9890端点失败与21143单元未独立复积分的状态保持。
