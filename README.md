# de Sitter 全息研究

## 2026-11-08 集中快照：选态密度与曲率动力学

按本轮要求集中发布 polynomial-density、composite、composite-transfer、selected-state 四个阶段，附四份已编译 PDF、LaTeX 与封存证据。[阅读入口与发布范围](notes/research-release-2026-11-08-selected-state.md)。下文“本地”“未 Git”保留各阶段记录时点的含义。后续研究接高阶双局域后代及固定电荷扇区的一致时间尾。

## 当前进展：选态密度与曲率时间关联（2026-11-08，本地）

[七页 PDF](output/pdf/de-sitter-susskind-slater-selected-state.pdf)、[LaTeX](manuscript/susskind-slater-selected-state/main.tex)、[中文证明](notes/susskind-slater-selected-state.md)和[来源](literature/susskind-slater-selected-state-sources.md)已保存。边界源比较将单selected激发密度推广到完整 h,u,x 协变代数；其 h,x restriction 为平稳密度 1−√q+√q h²，曲率等时谱仍与迹态相同，时间关联却不同。实际 exact-pinched 曲率的固定时间极限已接入；q=1/4、t=2 的关联负号由独立有理谱矩及余项界严格认证，迹态同期严格为正。

2186 密度词、2953 边界源、1022 reduced 词及三物理参数的双方法动力学核对通过。七页双遍编译/逐页目检；[证据清单](results/susskind-slater-selected-state/local-evidence.json)绑定三科学输出、15 当前/98 父科学文件。下一步高双局域后代的 central-sector 与 selected 全时间接口；本轮只本地，无新 Git/issue。

## 当前进展：曲率的显式开放双弦递推（2026-11-08，本地）

[六页 PDF](output/pdf/de-sitter-susskind-slater-composite-transfer.pdf)、[LaTeX](manuscript/susskind-slater-composite-transfer/main.tex)、[中文推导](notes/susskind-slater-composite-transfer.md)和[来源](literature/susskind-slater-composite-transfer-sources.md)已完成。归一化碰撞源的局域性给完整几何分类；五种开放端点及两类链接将全部固定 h,W,j,M 混合词化为显式递推，闭弦历史可合并。可达状态空间的正 Gram 商恢复共同 GNS，保留同 J 关联。

2774 源指数、994 矩比较、4096 旧六字词及五个长词已精确核对；八曲率插入例最多保留 21 个中间状态。双遍编译、六页目检和 11 当前/86 父科学文件核验见[证据清单](results/susskind-slater-composite-transfer/local-evidence.json)。下一步单选态激发在完整协变三场代数中的密度，再接高双局域后代与体生成元。当前仅本地；上一批 Git/issue 已完成，不重复发布。

## 当前进展：实际曲率与 selected 态的同 J 共同代数（2026-11-08，本地）

[八页 PDF](output/pdf/de-sitter-susskind-slater-composite.pdf)、[LaTeX](manuscript/susskind-slater-composite/main.tex)、[中文推导](notes/susskind-slater-composite.md)及[来源](literature/susskind-slater-composite-sources.md)已完成。分开空/占据局部矩阵元的多源公式构造 h,u,J₂,M 的全部固定混合矩；全部G弦时严格恢复旧tracial源指数，直接补齐central-sector高矩与exact pinching。新density界把实际曲率和中性斜率传到扩大后的多项式制备；协变三场在共同GNS上为有界正常表示。曲率与双selected激发的非零重叠说明同J收缩不能用独立probe替代。

592704项有限CAR、326源指数、340旧selected词及混合Gram已核；八页双遍编译/逐页目检，[证据清单](results/susskind-slater-composite/local-evidence.json)核11当前/72父科学文件。下一步高效composite transfer及完整bilocal后代；本轮仅本地，无新Git/issue。

## 当前进展：固定多项式制备中的精确投影（2026-11-08，本地）

[八页 PDF](output/pdf/de-sitter-susskind-slater-polynomial-density.pdf)、[LaTeX 完整证明](manuscript/susskind-slater-polynomial-density/main.tex)、[中文笔记](notes/susskind-slater-polynomial-density.md)及[来源](literature/susskind-slater-polynomial-density-sources.md)已保存。半填充 DS 中，near-Slater 权重维数下界控制所有部分 tensor flattenings，再经固定 Schatten 高斯迭代与条件 Wick 去耦，得到任意固定 Gaussian degree 的一致轨道密度矩界。实际 exact pinching 因而可传到所有固定 selected 多项式制备 core，以及多项式准备之间的协变观测块。二次制备 M†h、(M†)² 的归一斜率响应已显式计算。

[证据清单](results/susskind-slater-polynomial-density/local-evidence.json)绑定首次成功精确检查、9 当前/11 父文件和双遍编译；八页已逐页目检。多块非多项式 pinching 与 selected 因子交错的所有原始矩仍需额外控制。下一步实际 J₂ 与 M 的共同表示，再接完整 bilocal 和体生成元。本轮仅本地保存，上一批集中发布 a5be605 与 issue 回复 6710690876 不重复。

## 2026-11-08 集中研究快照：Slater 态与精确投影

本次集中整理 stationary、covariant、neutral-slope、orbit-density 四个阶段，附四份已编译 PDF、LaTeX、代码与封存证据。[发布范围与阅读入口](notes/research-release-2026-11-08-slater-density.md)。下文各阶段的“本地”“未 Git”保留其记录时点含义；当前已完成范围止于线性制备及协变后置 dressing，非线性制备仍为下一任务。

## 当前进展：单激发中的精确能量投影（2026-11-08，本地）

[六页 PDF](output/pdf/de-sitter-susskind-slater-orbit-density.pdf)、[完整 LaTeX](manuscript/susskind-slater-orbit-density/main.tex)、[中文证明](notes/susskind-slater-orbit-density.md)和[来源](literature/susskind-slater-orbit-density-sources.md)已保存。有限N Slater轨道的权重结构同时控制Gaussian synthesis map的两侧variance，给线性制备密度的一致矩界 E tr R_P^r≤(2r−1)!! a_N^r。对M†制备，a_N=1，二阶矩≤3。由此把真实exact pinching误差传入M†|O〉态，完成该制备内任意固定协变斜率关联；归一二点回到1+λ√q(1+q)。下一步处理二次及多次交错selected制备，再接J₂与完整bilocal。

[证据清单](results/susskind-slater-orbit-density/local-evidence.json)绑定精确有限检查、11当前/14父文件、最终双遍编译与六页目检。非线性selected制备尚未包含在本轮定理中；本轮仅本地保存。

## 当前进展：实际中性斜率与选态中的有限分辨率极限（2026-11-08，本地）

[八页 PDF](output/pdf/de-sitter-susskind-slater-neutral-slope.pdf)、[LaTeX](manuscript/susskind-slater-neutral-slope/main.tex)、[中文推导](notes/susskind-slater-neutral-slope.md)和[来源](literature/susskind-slater-neutral-slope-sources.md)已保存。实际未投影中性斜率通过精确rank系数接入u，覆盖全部固定selected混合词。新的共同primitive分解证明h,u在整个Slater空间上是标准qFock表示的放大，并识别A+/K≅C*_q(R²)。单个M†激发对应显式正密度ρ₁=(γh−√λγu)²+[1−√q(1+λ)]I；由此控制时间平均，迹方差为κ_u/T渐近，得到先双标度、后长平均窗的实际选态极限。精确有限N pinching已覆盖纯协变词，含M的接口仍为下一任务。

三份科学输出、14当前/15父文件哈希、双遍编译及八页目检见[证据清单](results/susskind-slater-neutral-slope/local-evidence.json)。本轮无新Git或issue回复。

## 当前进展：同 J 协变通道与非平凡平稳态（2026-11-08，本地）

[九页 PDF](output/pdf/de-sitter-susskind-slater-covariant.pdf)、[LaTeX完整证明](manuscript/susskind-slater-covariant/main.tex)、[中文笔记](notes/susskind-slater-covariant.md)及[文献来源](literature/susskind-slater-covariant-sources.md)已完成。加入实际微观大小涨落 u=(p−L−λ)h/√λ 后，归一化局部源展开给全部同 J 混合词及紧时间窗收敛。uΩ 在旧 h,M cyclic 空间外的范数平方为 1−λ/(e^(2λ)−1)>0。扩大的代数存在保留 h,u 非平凡动力学的平稳迹态，其 GNS 为标准两色 q-Gaussian II₁ 因子；这是指定态的表示，未识别整个 C* 商代数。三个科学输出、双遍编译、九页目检及父阶段哈希由[证据清单](results/susskind-slater-covariant/local-evidence.json)绑定。

下一步处理实际未 pinching neutral slope 的 rank 系数和 selected 态中的 pinching 控制，再接 J₂、完整 bilocal 与空间/质量生成元。此前 Git 推送及 issue 回复已完成；stationary、covariant 两个后续阶段先本地保存，沿用低频发布安排。

## 当前进展：平稳态分类与协变可观测量接口（2026-11-08，本地）

[八页 PDF](output/pdf/de-sitter-susskind-slater-stationary.pdf)、[LaTeX证明](manuscript/susskind-slater-stationary/main.tex)、[中文笔记](notes/susskind-slater-stationary.md)和[文献](literature/susskind-slater-stationary-sources.md)已保存。当前 h,M 的范数代数是 C*(h)+K(H)；所有平稳态只保留能量商，含 M 的有界关联为零，全 B(H) 上没有正常平稳密度。有限时间平均给近似平稳态，其真空信号与纯度同为 κ/T 渐近。与此同时，U(N)协变算符族的Slater联合平均精确回到固定扇区迹，为回接bilocal提供路径。下一步计算实际协变生成元与M的同J混合Gram，建立更大的共同表示；完整体字典仍开放。

本轮科学输出及[证据清单](results/susskind-slater-stationary/local-evidence.json)本地封存，八页已双遍编译、核版。上一轮集中发布为[560d683](https://github.com/wkryyds4-netizen/de-sitter-holography/commit/560d683503beb9b6d85832e0bbe6b65f45f54d31)，[issue回复](https://github.com/wkryyds4-netizen/de-sitter-holography/issues/1#issuecomment-6709573165)已完成；本轮无新增Git或评论。

## 当前进展：同 J 开弦传递与选择激发谱（2026-11-08）

[九页 PDF](output/pdf/de-sitter-susskind-slater-transfer.pdf)、[LaTeX](manuscript/susskind-slater-transfer/main.tex)、[中文推导](notes/susskind-slater-transfer.md)及[文献接口](literature/susskind-slater-transfer-sources.md)已完成。全 same-J word 规则现已化为带正 Gram 度量的显式 H/P 开弦传递；Hamiltonian 分解为 q-Jacobi 阶梯，选择激发的谱权重为 1−√q+√q E²。三份科学输出、最终双遍编译与九页核版由[证据清单](results/susskind-slater-transfer/local-evidence.json)绑定。下一步研究该表示中的静态观察者态，再接微观空间与质量生成元。

用户本轮要求集中 Git 并回复 issue，本次整理 joint-fourth、joint-words、transfer 三个阶段。[发布说明](notes/research-release-2026-11-08.md)记录范围；下方“仅本地”等文字保留历史时点，不代表本次发布状态。Issue #1 的完整重积分仍未完成。

## 当前进展：全部 same-J 联合 word 与有限时间构造（2026-11-07，本地）

[九页PDF](output/pdf/de-sitter-susskind-slater-joint-words.pdf)、[LaTeX完整证明](manuscript/susskind-slater-joint-words/main.tex)、[中文推导](notes/susskind-slater-joint-words.md)与[来源](literature/susskind-slater-joint-word-sources.md)已保存。精确site多项式保留CAR负号，归一化缺陷级数控制高重交集和含H的backward配对，证明任意固定joint word的γ^b q^ν规则；有限N行列covariance收缩给全次数尾，接正GNS、自伴H和紧时间窗实际关联。1,461,120六点标签、1,215系数、1,814局部表、729六字与time界核对完成，九页双遍编译/逐页目检。下一步全层transfer/谱核及state preparation，再物理态和体字典；全joint自平均与finite-G仍开放。本轮仅本地，不Git/issue，见[NEXT](NEXT.md)。

## 当前进展：同一 J 的联合四点与 Hamiltonian 前两层（2026-11-07，本地）

[七页 PDF](output/pdf/de-sitter-susskind-slater-joint-fourth.pdf)、[LaTeX](manuscript/susskind-slater-joint-fourth/main.tex)、[推导](notes/susskind-slater-joint-fourth.md)与[原文接口](literature/susskind-slater-joint-sources.md)已保存。同 J 交叉收缩接既有 EGUE covariance，给全部二/四次矩及七维 degree≤2 正空间；第一层 H 作用含范数平方 1−q 的新交换子方向。实际 disorder complement 的归一四矩为 2+√q+2q，不能补成标准 q 模式。四个有限 CAR 系统、62 个符号恒等式及12个文献标量比较已核，PDF双遍编译/逐页目检完成。下一步三弦/六阶，再全联合表示与真实演化；本轮尚非全阶闭包或体字典。

本次用户授权的集中 Git 已推送 [6e0f40a](https://github.com/wkryyds4-netizen/de-sitter-holography/commit/6e0f40af31f2ebf0ff0cf8f842244aaa04d82462)，[issue #1 回复](https://github.com/wkryyds4-netizen/de-sitter-holography/issues/1#issuecomment-6708558378)已发送并读回核对，issue保持开启。下方历史阶段的“未发布/未回复”保留为当时记录。新增 joint-fourth 阶段仅本地保存，继续位置见 [NEXT](NEXT.md)。

## 当前进展：同一 J 的 Slater 极化与微观 canonical 阶梯（2026-11-07，本地）

[九页PDF](output/pdf/de-sitter-susskind-slater-ladder.pdf)、[LaTeX](manuscript/susskind-slater-ladder/main.tex)、[主推导](notes/susskind-slater-ladder.md)、[完整证明](notes/susskind-slater-ladder-proof.md)及[文献接口](notes/susskind-slater-ladder-literature-interface.md)已保存。先核纯态DSSYK前驱，再从原H取半填充Slater适配creation块；forward精确计数与backward余项给全部固定word的dressed boson律。实际占据数去dressing得到指定态固定激发core的渐近CCR。完整H首步仍有Poisson(2lambda)缺额，单ladder只捕获sqrt(q)范数平方；下一步同J联合(H,M)动力学。1,537,310个精确标签项、一份科学JSON、双遍编译及九页目检完成。完整体字典仍开放；仅本地、不Git，见[NEXT](NEXT.md)。

## 当前进展：完整有限时间双局域四点（2026-11-07，本地）

[十一页PDF](output/pdf/de-sitter-susskind-bilocal-fourth-tail.pdf)、[LaTeX](manuscript/susskind-bilocal-fourth-tail/main.tex)、[主推导](notes/susskind-bilocal-fourth-tail.md)及[实际四点传递](notes/susskind-bilocal-full-four-point-transfer.md)已保存。受限covariance/ancilla界、精确CP通道链及Clifford Riesz接普通带电重组，补齐有限N一致的全次数高范数尾；共同源图与旧强返回因此构造完整finite-time bilocal四点。三组科学检查、双遍编译和十一页目检完成。下一步转到态/极化与微观体生成元，有限G与完整体字典仍开放。本轮仅本地；GitHub仍为已核issue #1，暂无新评论。见[NEXT](NEXT.md)。

## 当前进展：低双局域通道的共同四点代数（2026-11-07，本地）

[九页PDF](output/pdf/de-sitter-susskind-bilocal-four-point.pdf)、[LaTeX](manuscript/susskind-bilocal-four-point/main.tex)、[推导](notes/susskind-bilocal-four-point.md)、[多源证明](notes/susskind-bilocal-multisource-proof.md)、[实际算符接口](notes/susskind-bilocal-four-point-transfer.md)已保存。核清旧J2定理的共同样本接口，将中性W1接入长度(1,1,2)的三字母mixed-q代数；实际曲率/中性斜率的高范数强桥给固定实时间四点及严格正交换子。四科学输出、双遍编译及九页目检完成。返回误差升至L4的O(p^-1/7)，下一步完整时间级数L4可和尾，再态/极化和微观体生成元；完整bilocal四点仍未宣称完成。仅本地、不Git，见[NEXT](NEXT.md)。

## 当前进展：完整双局域涨落的强重建（2026-11-07，本地）

[七页PDF](output/pdf/de-sitter-susskind-bilocal-reconstruction.pdf)、[LaTeX](manuscript/susskind-bilocal-reconstruction/main.tex)、[推导](notes/susskind-bilocal-reconstruction.md)、[完整证明](notes/susskind-bilocal-reconstruction-proof.md)、[来源](literature/susskind-bilocal-reconstruction-sources.md)已保存。静态memory的收缩dilation加旧正力八阶矩/集体二阶界，给实际bilocal到完整memory滤波的O(p^−1/3)全电荷强误差；charged幅度精确化简，neutral写成重节点连续的指数差商。两科学输出、双遍编译及七页目检完成。下一步共同J四点/交换子、适用态和canonical结构，再构造体空间/质量生成元；仅本地、不Git，见[NEXT](NEXT.md)。

## 当前进展：完整有限时间联合记忆核（2026-11-07，本地）

[九页PDF](output/pdf/de-sitter-susskind-memory-joint-regression.pdf)、[LaTeX](manuscript/susskind-memory-joint-regression/main.tex)、[推导](notes/susskind-memory-joint-regression.md)、[完整证明](notes/susskind-memory-joint-regression-proof.md)、[来源](literature/susskind-memory-joint-regression-sources.md)已保存。联合电荷/能量回归与same-copy强归约补齐charged谱；无对角质量传递exact pinching，ordinary尾后再Hermite重求和，完成指定模型内AA/AB/BA/BB的收敛构造。三份科学输出、双遍编译和九页目检完成。下一步从记忆核重建完整bilocal涨落，接体空间/质量与canonical dressing；体字典和有限G动力学仍开放。仅本地、不Git，见[NEXT](NEXT.md)。

## 当前进展：完整带电局部源与电荷静态背景（2026-11-07，本地）

[八页PDF](output/pdf/de-sitter-susskind-memory-oriented-transfer.pdf)、[LaTeX](manuscript/susskind-memory-oriented-transfer/main.tex)、[推导](notes/susskind-memory-oriented-transfer.md)、[来源](literature/susskind-memory-oriented-transfer-sources.md)已保存。精确full-Fock loss接口接上全时间加权替换与可和尾；带源singleton归一给所有内部Gaussian收缩的闭式图规则。识别无直接桥仍存在的charge静态背景，扣除后恢复旧Y₋方差。五份科学输出、双遍编译及八页目检已完成。下一步联合电荷/能量回归、AA/AB/BA谱，再接体字典。当前仅本地、不Git，见[NEXT](NEXT.md)。

## 当前进展：完整中性记忆核与精确去静态（2026-11-07，本地）

[PDF](output/pdf/de-sitter-susskind-memory-wick-descendants.pdf)、[LaTeX](manuscript/susskind-memory-wick-descendants/main.tex)、[推导](notes/susskind-memory-wick-descendants.md)、[来源](literature/susskind-memory-wick-descendants-sources.md)已保存。有序Wick后代及chaos投影接入共同J图；切图回归使无桥组抵消，剩余桥给所有外部Hermite次数的指数界，因此残量谱绝对连续，能传递exact pinching。旧时间级数尾给完整neutral B协方差的收敛构造。42有限OU/Hermite、6非空抵消、59矩多项式已核。下一步charged A的forward/reverse及AA/AB/BA，完整联合核和体字典仍开放。只本地、不Git，见[NEXT](NEXT.md)。

## 当前进展：秩标记的局域大小桥与离散分布（2026-11-07，本地）

[十页PDF](output/pdf/de-sitter-susskind-memory-size-bridge.pdf)、[LaTeX](manuscript/susskind-memory-size-bridge/main.tex)、[推导](notes/susskind-memory-size-bridge.md)、[来源](literature/susskind-memory-size-bridge-sources.md)已保存。接上已有非Abelian大小代数与Clifford semigroup，证明固定次数rank/half-length桥，得到两标记共同J弦图、全部等时矩和完整离散缺额分布；实际rank law在TV下收敛。5716合法图的25884精确比较及有限CAR/分布检查通过，PDF含分布图并逐页核版。下一步有序Wick后代、全时间中性核及零频质量；完整联合核和体字典仍开放。只本地、不Git，见[NEXT](NEXT.md)。

## 当前进展：中性记忆场的全时间标记级数（2026-11-07，本地）

[九页PDF](output/pdf/de-sitter-susskind-memory-neutral-marked-series.pdf)、[LaTeX](manuscript/susskind-memory-neutral-marked-series/main.tex)、[推导](notes/susskind-memory-neutral-marked-series.md)已保存。固定次数高矩与实际指数的rank尾接上有限参数Jacobi递推，将整个中性分量表示为最高chaos标记场及其全部收缩后代；有显式可和尾和全电荷L²误差控制。87完整CAR系数比较及450精确斜率支配通过。下一步求与charged defect的共同极限矩，并控制精确pinching的零频质量；完整联合核仍待计算。仅本地，无Git发布，见[NEXT](NEXT.md)。

## 当前进展：实际记忆场的两分量归约（2026-11-07，本地）

[八页PDF](output/pdf/de-sitter-susskind-memory-two-component.pdf)、[LaTeX](manuscript/susskind-memory-two-component/main.tex)、[推导](notes/susskind-memory-two-component.md)已保存。三色条件期望与全实时间可和尾证明four-u修正趋守恒背景，pinching后消失；实际field归约为charged defect+neutral residual。neutral前两系数接旧谱，二阶mixed cross严格非零，作为全时间核的校准。9828有理q恒等式、24CAR积分、12有限差式及360rank斜率已核。下一步共同J全时间联合核，见[NEXT](NEXT.md)；仅本地，无Git发布。

## 当前进展：实际有限时间记忆场的统一二阶界（2026-11-07，本地）

[九页PDF](output/pdf/de-sitter-susskind-memory-uniform-bound.pdf)、[LaTeX](manuscript/susskind-memory-uniform-bound/main.tex)、[推导](notes/susskind-memory-uniform-bound.md)已保存。加权charged传递、同副本Dyson展开与有限参数Jacobi桥证明实际X0的全charge统一L²界；结合首返守恒结果，实际记忆场与X0的两点函数可渐近替换。12完整Hermite案例、9实时间余项、105有理式及6高精度谱外展开已核。下一步求三个有界分量的共同J联合涨落核及cross；完整核和体字典仍开放。仅本地保存，见[NEXT](NEXT.md)。

## 当前进展：真实相邻电荷传播的OU近似（2026-11-07，本地）

[八页PDF](output/pdf/de-sitter-susskind-memory-charged-transfer.pdf)、[LaTeX](manuscript/susskind-memory-charged-transfer/main.tex)、[推导](notes/susskind-memory-charged-transfer.md)已保存。正反向缺陷与polar传递接全部Gaussian chaos，证明实际F_m=T_m exp(isH_{m−1})到同一原J的OU平滑传播子的全charge L²误差为O_S(1/N)。六完整CAR系统、230有理源方差及三JSON已核。下一步加时间导数/能量权重，求m放大后的力核缺陷；OU参考尚非samplewise conserved。仅本地保存，见[NEXT](NEXT.md)。

## 外部复核：issue #1（2026-11-06，本地）

已读[issue #1](https://github.com/wkryyds4-netizen/de-sitter-holography/issues/1)及8页附件，核对KS原文，独立有理算术重组支持所提供区间的普通/整体相位优化正下界。9890项旧端点重现失败保留，未在本机重积分21143单元。见[三页核查PDF](output/pdf/de-sitter-issue-1-review.pdf)、[记录](notes/issue-1-review.md)。无外部回复/提交/推送；主线继续下面的有限记忆计算。

## 当前进展：完整普通传播子的全次数rank界（2026-11-07，本地）

[八页PDF](output/pdf/de-sitter-susskind-memory-wick-rank.pdf)、[LaTeX](manuscript/susskind-memory-wick-rank/main.tex)、[推导](notes/susskind-memory-wick-rank.md)已保存。正rank链的全k二矩接Gaussian creation/annihilation和逐样本Duhamel导数，证明实际exp(isH)的rank–chaos缺额在固定实时间窗统一有界，覆盖所有低chaos。独立线性charged probe四点及精确能量回归给X0的正充分上界。下一步接实际相邻charge force；μ(pNcal)尚不是逐样本守恒函数，未宣布X0已紧。三JSON、八页核版与哈希本地保存，见[NEXT](NEXT.md)。

## 前一进展：第一次返回的守恒背景（2026-11-06，本地）

[六页PDF](output/pdf/de-sitter-susskind-memory-return-background.pdf)、[LaTeX](manuscript/susskind-memory-return-background/main.tex)、[推导](notes/susskind-memory-return-background.md)已保存。双tag计数接已有q-Mehler条件期望，把普通核的领先能量回归升为固定有限时间L²强近似；第一次返回趋守恒函数，精确去能量块后实际场与普通项X0的L²差趋零。下一步集中于X0次领先涨落的统一二阶界与极限；协方差替换仍需该界，旧六产品归约保留。2373精确恒等式、27数值和六页核版本地保存，见[NEXT](NEXT.md)。

## 前一进展：有限记忆的两次返回归约（2026-11-06，本地）

[六页PDF](output/pdf/de-sitter-susskind-memory-return-expansion.pdf)、[LaTeX](manuscript/susskind-memory-return-expansion/main.tex)、[推导](notes/susskind-memory-return-expansion.md)已保存。Gaussian力矩与有限Volterra迭代给集体场O(p^-2)、协方差O(p^-1)误差，保留完整时间指数；下一步计算六个共同J、pinched普通传播产品的极限。旧KS支线issue初核见上方独立记录。仅本地保存，见[NEXT](NEXT.md)。

## 当前进展：全电荷记忆谱与精确去静态接口（2026-11-06，本地）

[六页PDF](output/pdf/de-sitter-susskind-memory-charge-sum.pdf)、[LaTeX](manuscript/susskind-memory-charge-sum/main.tex)、[推导](notes/susskind-memory-charge-sum.md)已保存。全sector正回归分母、增长窗余项及binomial尾将(Y,L,U2)谱接到β=μ=0全迹；charge Gaussian与该双能量matrix因子化。精确samplewise pinching保留前三个Taylor极限系数。6766精确方差案例、15CAR sectors、12全rank组及4尾求和已核。下一步实际finite-memory covariance，仅本地、不commit/push/ZIP，见[NEXT](NEXT.md)。

## 当前进展：中央原始记忆系数与线性滤波场（2026-11-06，本地）

[六页PDF](output/pdf/de-sitter-susskind-memory-central-raw.pdf)、[LaTeX](manuscript/susskind-memory-central-raw/main.tex)、[推导](notes/susskind-memory-central-raw.md)已保存。五个raw后代新增电荷项均为守恒B修正；中央F的斜率、范数及双向lift接上完整(Y,L,U2)谱。123个有限系数参考、36组全rank高精度及旧证书hash已核。下一步窗外charge求和，再finite-memory与体字典。仅本地，不commit/push/ZIP，见[NEXT](NEXT.md)。

## 当前进展：中央双向力谱与投影二次core（2026-11-06，本地）

[六页PDF](output/pdf/de-sitter-susskind-memory-central-spectrum.pdf)、[LaTeX](manuscript/susskind-memory-central-spectrum/main.tex)、[推导](notes/susskind-memory-central-spectrum.md)已保存。中央迹相对fullFock每个固定marked图有O(1/p)误差；一般charge force桥接补N−2m项，得到中央双向力谱和实际投影W2core的V2 μ_q²两能量测度。324高精度迹点、108精确CAR恒等式和旧证书hash已核。下一步中央raw U2后代/守恒项与F接口，只本地、无commit/push/ZIP，见[NEXT](NEXT.md)。

## 当前进展：中央四方向Gram与正记忆残量（2026-11-06，本地）

[六页PDF](output/pdf/de-sitter-susskind-memory-central-gram.pdf)、[LaTeX](manuscript/susskind-memory-central-gram/main.tex)、[推导](notes/susskind-memory-central-gram.md)已保存。相邻charge极化加neutral矩O(1/(Np))漂移独立证明S–O交叉，完成偶N中央紧窗四Gram及正W2core投影残量；不依赖此前待证的有限Ncharged候选式。758精确比、36全秩组及3217完整CAR坐标对已核。下一步中央energy insertions/raw谱，仅本地、无commit/push/ZIP，见[NEXT](NEXT.md)。

## 当前进展：中央产品传递与输出缺额范数（2026-11-06，本地）

[六页PDF](output/pdf/de-sitter-susskind-memory-central-cross.pdf)、[LaTeX](manuscript/susskind-memory-central-cross/main.tex)、[推导](notes/susskind-memory-central-cross.md)已保存。相邻charge等距传递与两端六阶尾给输出缺额O(p^−2/3)传递，完成偶N中央紧窗(B,L_in,S)Gram。S–O候选通过60精确CAR恒等式及36全秩组；general charged gradient解析步骤仍待补，完整四Gram未标为完成。继续仅本地、不commit/push/ZIP，见[NEXT](NEXT.md)。

## 当前进展：中央窗口的相干范数与core强展开（2026-11-06，本地）

[六页PDF](output/pdf/de-sitter-susskind-memory-central-norm.pdf)、[LaTeX](manuscript/susskind-memory-central-norm/main.tex)、[推导](notes/susskind-memory-central-norm.md)已保存。一般charge相干A范数与输入尾接上中央输出六阶矩，完成gA的L²强替换、(B,L_in,O)三方向Gram和含charge守恒修正的完整W2core展开。四JSON、31466完整CAR坐标对及36全秩高精度组本地保存。下一步S_N相关Gram及动态插入，见[NEXT](NEXT.md)。

## 当前进展：中央电荷窗口的秩尾与守恒修正（2026-11-06，本地）

[六页PDF](output/pdf/de-sitter-susskind-memory-central-window.pdf)、[LaTeX](manuscript/susskind-memory-central-window/main.tex)、[推导](notes/susskind-memory-central-window.md)已保存。正rank乘积六阶矩扩至m=N/2+z√N的紧窗口；共同方差留下p(Cm/C0−1)→4λ²x(1−x)，raw gH另增4λx(1−q)²守恒项。三JSON、60精确filter恒等式和六页核版仅本地保存。下一步实际memory的charged范数及Gram，见[NEXT](NEXT.md)。

## 当前进展：恢复只扣守恒项的二次记忆谱（2026-11-06，本地）

[六页PDF](output/pdf/de-sitter-susskind-memory-raw-spectrum.pdf)、[LaTeX](manuscript/susskind-memory-raw-spectrum/main.tex)、[推导](notes/susskind-memory-raw-spectrum.md)已保存。恢复旧投影去掉的动态后代，得到U2完整正谱及Y/L交叉谱；有序小λ归一形状为(9sech⁸−2sech⁶)/7，恢复的Y后代占领先强度11/35。三个科学JSON、六页核版与哈希仅本地保存。下一步中央charge窗及finite-memory余项，见[NEXT](NEXT.md)。继续减少Git频率，本轮无commit/push。

## 当前进展：三个记忆场的完整混合谱（2026-11-06，本地）

[六页PDF](output/pdf/de-sitter-susskind-memory-mixed-spectrum.pdf)、[LaTeX](manuscript/susskind-memory-mixed-spectrum/main.tex)、[推导](notes/susskind-memory-mixed-spectrum.md)已保存。74类reverse端点结构和上电荷lift补齐Y+–R/Y+–F/Q4–F，得到X–L完整双能量cross0与(Y,L,X)正3×3谱。下一步恢复raw二次系数被投影掉的动态后代，再中央charge窗/finite-memory。3JSON、失败证据与六页核版本地保存，无commit/push/ZIP，见[NEXT](NEXT.md)。

## 当前进展：双向费米子力与完整投影记忆谱（2026-11-06，本地）

[七页PDF](output/pdf/de-sitter-susskind-memory-oriented-kernel.pdf)、[LaTeX](manuscript/susskind-memory-oriented-kernel/main.tex)、[推导](notes/susskind-memory-oriented-kernel.md)已保存。输出缺额强桥接到产生/湮灭force平均，128模板证书补齐Q4–Y+全能量cross；raw mZ2在指定八方向span外的完整两点谱为V4 μ_q4+V2 μ_q2，并与Y−全部能量后代cross0。下一步求与三阶memory场L的mixed谱，再中央charge窗、finite-memory与体字典。4JSON/七页核版本地保存，无commit/push/ZIP，见[NEXT](NEXT.md)。

## 当前进展：raw 记忆系数的双 chaos 残量（2026-11-06，本地）

[七页PDF](output/pdf/de-sitter-susskind-memory-descendants.pdf)、[LaTeX](manuscript/susskind-memory-descendants/main.tex)、[推导](notes/susskind-memory-descendants.md)已保存。恢复四个raw后代低chaos，证明旧W4 Gram的半填充迹桥接；完整mZ2在指定八方向span外保留Q4与二次输出缺额两部分，等时方差V4+V2。下一步求两者联合能量谱，不能由chaos等时正交推出时间cross0。3600收缩、3JSON和七页核版本地保存，无commit/push/ZIP，见[NEXT](NEXT.md)。

## 当前进展：完整二次记忆 core 的强极限（2026-11-06，本地）

[六页PDF](output/pdf/de-sitter-susskind-memory-filter-limit.pdf)、[LaTeX](manuscript/susskind-memory-filter-limit/main.tex)、[推导](notes/susskind-memory-filter-limit.md)已保存。正rank乘积六阶矩加增长窗口给g_A的L²替换，完成固定半填充W2core四方向分解；扣(B,O,JhF)后残量严格正，后小λ趋1024Js8。下一步恢复W4后代低chaos、组装rawZ2并扩中央charge窗。3JSON及六页核版本地保存，继续降低Git频率，不commit/push/ZIP，见[NEXT](NEXT.md)。

## 当前进展：半填充范数与正记忆残量（2026-11-06，本地）

[六页PDF](output/pdf/de-sitter-susskind-memory-residual-norm.pdf)、[LaTeX](manuscript/susskind-memory-residual-norm/main.tex)、[推导](notes/susskind-memory-residual-norm.md)已保存。独立求固定半填充A范数，补齐四方向等时Gram；扣H²/A/JhF后输出缺额仍有正残量，小λ趋16Js4。输入尾已控制，下一步A的输出rank尾及g_A极限。68544全CAR坐标对、符号相消和12组全rank高精度验证本地保存，不commit/push/ZIP，见[NEXT](NEXT.md)。

## 当前进展：相邻电荷的相干交叉项（2026-11-06，本地）

[六页PDF](output/pdf/de-sitter-susskind-memory-coherent-cross.pdf)、[LaTeX](manuscript/susskind-memory-coherent-cross/main.tex)、[推导](notes/susskind-memory-coherent-cross.md)已保存。求出真实S–A等时cross；扣H²后cross恰等于S残量范数，下一步求固定半填充的A−S范数及A加权尾。10完整charged maps、26227全CAR坐标对和12组全rank高精度核对通过。只本地保存，不commit/push/ZIP，见[NEXT](NEXT.md)。

## 当前进展：输出缺额的奇偶矩与 Gram（2026-11-06，本地）

[五页PDF](output/pdf/de-sitter-susskind-memory-parity-gram.pdf)、[LaTeX](manuscript/susskind-memory-parity-gram/main.tex)、[推导](notes/susskind-memory-parity-gram.md)已保存。由既有Racah差分式求signed前二矩，得到(B、输入缺额、输出缺额)完整等时DS Gram。输出缺额扣H²后仍有正范数，小λ趋16Js4，并与JhF正交；A cross尚未求出。417矩、1468全CAR坐标对及五页核版本地保存，无commit/push/ZIP，见[NEXT](NEXT.md)。

## 当前进展：记忆 core 的一般双电荷公式（2026-11-06，本地）

[六页PDF](output/pdf/de-sitter-susskind-memory-two-charge.pdf)、[LaTeX](manuscript/susskind-memory-two-charge/main.tex)、[推导](notes/susskind-memory-two-charge.md)已保存。普适charged transfer给完整Q_Ecal两个显式有限N振幅；合并O(m)项后，H²部分的三个deficit线性极限具有L²控制。下一步求output-deficit奇偶加权Gram、A加权尾及cross，再组装raw Z2谱。解析证明、87个完整linear maps、522个收缩核对及六页核版本地保存，无commit/push/ZIP，见[NEXT](NEXT.md)。

## 当前进展：相邻电荷打开乘积禁戒通道（2026-11-06，本地）

[四页PDF](output/pdf/de-sitter-susskind-memory-charge-lift.pdf)、[LaTeX](manuscript/susskind-memory-charge-lift/main.tex)、[推导](notes/susskind-memory-charge-lift.md)已保存。半填充选择规则和一般解析族证明：普通 H² 的 rank 滤波不足以表示记忆算符；N6p2m3 精确见证 P1H²=0、P1Q=99D，独立225方向收缩核对。改用两个电荷的相干乘法幅，另得 charge-removal 精确范数传递。下一步求包含 Ψ 的一般幅及 DS cross。四页核版/计算本地保存，无 commit/push/ZIP，见[NEXT](NEXT.md)。

## 当前进展：rank 滤波与双标度尾控制（2026-11-06，本地）

[六页PDF](output/pdf/de-sitter-susskind-memory-rank-blocks.pdf)、[LaTeX](manuscript/susskind-memory-rank-blocks/main.tex)、[推导](notes/susskind-memory-rank-blocks.md)和[来源](literature/susskind-memory-rank-blocks-sources.md)已保存。既有有限残量被分成11个共线rank块，差别来自通道权重；协方差复合另给正输出rank律和一般填充的前三矩极限，半填充偏移均值/方差均2λ。尾界接same-J对称乘积上界，下一步求真实memory的general reduced filter及DS残量权重。5份JSON及六页核版本地保存，不commit/push/ZIP，见[NEXT](NEXT.md)。

## 当前进展：记忆 core 的二次收缩（2026-11-06，本地）

[五页PDF](output/pdf/de-sitter-susskind-memory-quadratic-contraction.pdf)、[LaTeX](manuscript/susskind-memory-quadratic-contraction/main.tex)、[推导](notes/susskind-memory-quadratic-contraction.md)及[来源](literature/susskind-memory-quadratic-contraction-sources.md)已保存。相邻charge协方差给W2Ecal有限N精确式；完整Gram和逐样本scalar分离显示N5p2m3仍有旧方向之外的动态残量，而两个小系统闭合。下一步求该残量的DS强度与rank分解，不由有限N断言新粒子。四个计算JSON、五页核版与哈希本地保存；不commit/push/ZIP，见[NEXT](NEXT.md)。

## 当前进展：二次记忆系数与第四阶singlet（2026-11-06，本地）

[六页PDF](output/pdf/de-sitter-susskind-memory-quartic.pdf)、[LaTeX](manuscript/susskind-memory-quartic/main.tex)、[中文推导](notes/susskind-memory-quartic.md)及[来源](literature/susskind-memory-quartic-sources.md)已保存。精确求Z2及首返项分离；第四chaos扣旧后代后得到正动态Q4，807端点模板全阶相消接已有q⁴ Mehler核，另482模板给Q4–Y/R双能量cross0。小λ归一核sech⁸。下一步补core第二chaos与低阶收缩，完成raw Z2谱；未把纯Q4当完整memory。六页核版、计算及失败证据本地保存，无commit/push/ZIP，见[NEXT](NEXT.md)。

## 当前进展：能量背景与首对记忆涨落（2026-11-06，本地）

[六页PDF](output/pdf/de-sitter-susskind-memory-fluctuations.pdf)、[LaTeX](manuscript/susskind-memory-fluctuations/main.tex)、[中文推导](notes/susskind-memory-fluctuations.md)和[来源](literature/susskind-memory-fluctuations-sources.md)已保存。已有Mehler条件算符给能量回归背景及Vbg=G_(q²)−G_q²，证明集体尺度只扣scalar mean不足。明确零/一次守恒counterterms，得到首对memory系数的完整正谱矩阵；居中去掉raw三阶时间后代，smallλ tangent核为sech⁶。下一步二次memory、第四singlet及有限memory余项。六页核版、计算和失败证据仅本地保存，无commit/push/ZIP，见[NEXT](NEXT.md)。

## 当前进展：平均记忆核与费米子首阶修正（2026-11-06，本地）

[七页PDF](output/pdf/de-sitter-susskind-memory-force.pdf)、[LaTeX](manuscript/susskind-memory-force/main.tex)、[中文推导](notes/susskind-memory-force.md)及[来源](literature/susskind-memory-force-sources.md)已保存。同J带电force的标记收缩与两次Duhamel界给能量居中平均记忆p Etr Mhat(s)→4Js²G_q(s)，固定紧时间一致；积分给原始fermion首个1/p项，小λ恢复既有log cosh。精确renewal与一阶记忆系数接旧导数singlet。下一步求两个memory参数的算符涨落协方差，再接体介子内部核。本轮计算、七页核版与证据仅本地保存，不commit/push/ZIP，见[NEXT](NEXT.md)。

## 当前进展：完整bilocal的压缩演化与记忆核（2026-11-06，本地）

[六页PDF](output/pdf/de-sitter-susskind-bilocal-channel.pdf)、[LaTeX](manuscript/susskind-bilocal-channel/main.tex)、[中文推导](notes/susskind-bilocal-channel.md)和[来源](literature/susskind-bilocal-channel-sources.md)已保存。charge-removal等距构造给完整双时singlet的正核、组合偏差和精确记忆方程；记忆初值减守恒部分接旧Y通道。Duhamel界证明普通归一bilocal在固定紧时间DS趋单位，重建须保留涨落。下一步求完整正交记忆核及减背景涨落；本轮核版与证据本地保存，无commit/push/ZIP，见[NEXT](NEXT.md)。

## 当前进展：双局域到介子的重建接口（2026-11-06，本地）

[六页PDF](output/pdf/de-sitter-susskind-bilocal-interface.pdf)、[LaTeX](manuscript/susskind-bilocal-interface/main.tex)、[中文推导](notes/susskind-bilocal-interface.md)和[来源](literature/susskind-bilocal-interface-sources.md)已保存。接入Itakura已有投影，计算内部源的介子质量混合与高矩尾；推导complex SYK双时CAR导数的精确次数权重，集体零时曲率在DS保留。下一步求完整singlet投影及极化/dressing，不将未投影范数当作代数不闭合证明。计算、六页核版与哈希证据本地保存，继续降低Git频率，无commit/push/ZIP，见[NEXT](NEXT.md)。

## 当前进展：接入已知体介子谱与空间响应靶标（2026-11-06，本地）

[五页PDF](output/pdf/de-sitter-susskind-bulk-target.pdf)、[LaTeX](manuscript/susskind-bulk-target/main.tex)、[中文记录](notes/susskind-bulk-target.md)与[来源](literature/susskind-bulk-target-sources.md)已保存。MSS重整化质量字典明确接到FLZ的a=0有隙介子方程；两种独立表示复现前六级，差约6.93e−13。正二次型比较给平面谱相对位移≤π/p²，构造反向光锥响应的质量方差判据。下一步求微观bilocal/空间插入映射；既有体谱不当作DSSYK质量结果。五页核版、原始失败和证据保存，无commit/push/ZIP，见[NEXT](NEXT.md)。

## 当前进展：排序奇内积与空间平移检验（2026-11-06，本地）

[六页PDF](output/pdf/de-sitter-susskind-ordering-metric.pdf)、[LaTeX](manuscript/susskind-ordering-metric/main.tex)、[中文推导](notes/susskind-ordering-metric.md)与[来源](literature/susskind-ordering-metric-sources.md)已保存。从已有commutator核直接确定领先Casimir内积 A=Ctot−C1−C2，d2/d3相对谱为k(k+9)+12。构造保持范数的正自伴平移，并以精确范数变化证明普通自由两粒子质量与既有乘法顶点的条件性交织不相容。下一步回到Susskind体规范理论的实际空间字典。50条符号检查、独立积分及有限h Gram保存，六页核版完成；无commit/push/ZIP，见[NEXT](NEXT.md)。

## 当前进展：混合增长与排序反对称复合态（2026-11-06，本地）

[六页PDF](output/pdf/de-sitter-susskind-commutator-primary.pdf)、[LaTeX](manuscript/susskind-commutator-primary/main.tex)、[推导](notes/susskind-commutator-primary.md)与[来源](literature/susskind-commutator-primary-sources.md)已保存。应用LS现有公式得固定时间D23=48λcosh(4Js t)+O(λ²)；归一commutator K具有精确q⁵两点、非零三点顶点和可计算第四矩额外项。将文献既有排序奇primary落实到两个导数singlet。精确弦矩、三λ实时与50位独立和式核对完成；下一步把排序扇区内积接相对Casimir和体平移。无commit/push/ZIP，见[NEXT](NEXT.md)。

## 当前进展：二阶/三阶联合弦规则与混合对易子（2026-11-06，本地）

[六页PDF](output/pdf/de-sitter-susskind-mixed-cubic.pdf)、[LaTeX](manuscript/susskind-mixed-cubic/main.tex)、[推导](notes/susskind-mixed-cubic.md)与[原文接口](literature/susskind-mixed-cubic-sources.md)已保存。14类局部组件相消给(H,X,Rhat)任意固定阶的已有多probe规则，互交叉q⁶；接BINT不等长核得到mixed OTOC和D23。有限CAR6656项、局部精确多项式、三λ实时核对及六页核版完成。下一步求固定时间首阶并接两species体字典；无commit/push/ZIP，见[NEXT](NEXT.md)。

## 当前进展：扣除守恒部分后的三阶两点谱闭合（2026-11-06，本地）

[七页PDF](output/pdf/de-sitter-susskind-filtered-spectrum.pdf)、[LaTeX](manuscript/susskind-filtered-spectrum/main.tex)、[推导](notes/susskind-filtered-spectrum.md)和[来源](literature/susskind-filtered-spectrum-sources.md)已保存。精确charged lift把rank-filter桥接接入导数mark，证明ν_FY=ν_FR=0；结合旧ν_RY=0，得到U=M₃−P_N完整正双能量测度。小λ归一time核为(47sech⁶−36sech⁴)/11，二阶时间后代占领先等时强度9/11。下一步转mixed三/四点与体字典。CAR、full-rank、Jacobi/Mehler/Fourier及七页核版完成；无commit/push/ZIP，见[NEXT](NEXT.md)。

## 当前进展：原始三阶归约与交叉谱（2026-11-06，本地）

[七页PDF](output/pdf/de-sitter-susskind-cubic-cross.pdf)、[LaTeX](manuscript/susskind-cubic-cross/main.tex)、[推导](notes/susskind-cubic-cross.md)和[来源](literature/susskind-cubic-cross-sources.md)已保存。精确求出S/Z线性反项及1/p极限，把raw M₃扣守恒多项式后归为F、R和二阶singlet的能量/时间后代；177端点模板相消与谱绝对连续性给R–Y完整交叉谱为零。F–R及能量分辨F–Y仍待算，下一步直接接rank-filter与marked chord。有限交叉、全rank高精度及独立矩阵恒等式保存；按降低Git频率的要求，无提交/推送/ZIP，见[NEXT](NEXT.md)。

## 当前进展：三阶一次残量的完整q核（2026-11-06，本地）

[七页PDF](output/pdf/de-sitter-susskind-linear-kernel.pdf)、[LaTeX](manuscript/susskind-linear-kernel/main.tex)、[推导](notes/susskind-linear-kernel.md)和[来源](literature/susskind-linear-kernel-sources.md)已保存。以p/p−1/p−2系综协方差总变差桥接，控制任意固定Wick图，证明(H,F)正交q-Gaussian联合矩；F完整两点接已有r=q Mehler核，小λ为sech²(2Js t)，DS无静态atom。下一步求F与R/S/Z实时cross，完成raw三阶谱再做混合相互作用。精确比、six-slot收缩、联合四次矩与双核计算保存；无Git提交/推送/ZIP，见[NEXT](NEXT.md)。

## 当前进展：三阶一次残量的动态谱矩（2026-11-06，本地）

[六页PDF](output/pdf/de-sitter-susskind-cubic-dynamics.pdf)、[LaTeX](manuscript/susskind-cubic-dynamics/main.tex)、[中文推导](notes/susskind-cubic-dynamics.md)和[来源](literature/susskind-cubic-dynamics-sources.md)已保存。用已有Racah求和工具得到纯rank交叉协方差，完整保留同J收缩，证明F的归一动态曲率为2Js²(1−q)/λ，小λ趋8Js²。第四矩统一界进一步证明动态频率权重在DS中保留。下一步求F高阶双能量矩及全核，再接raw三阶交叉与混合相互作用。精确CAR、组合数和全rank高精度核对完成；无Git提交/推送/ZIP，见[NEXT](NEXT.md)。

## 当前进展：三阶算符的一次 Gaussian 分量（2026-11-06，本地）

[六页PDF](output/pdf/de-sitter-susskind-cubic-linear.pdf)、[LaTeX](manuscript/susskind-cubic-linear/main.tex)、[中文推导](notes/susskind-cubic-linear.md)与[来源](literature/susskind-cubic-linear-sources.md)已保存。求出W₁T₃的有限N不可约通道系数，精确反例证明它不总守恒；扣除逐电荷H投影后的F在DS中范数为64Js⁶q²，相对R小λ弱2λ。下一步求F的动态谱与同J交叉项，不能由正范数直接认定新粒子。42个通道、31个charged map、199个权重核对和六页核版完成；无Git提交/推送/ZIP，见[NEXT](NEXT.md)。

## 当前进展：三阶 singlet 的完整两点核（2026-11-06，本地）

[六页PDF](output/pdf/de-sitter-susskind-cubic-kernel.pdf)、[LaTeX](manuscript/susskind-cubic-kernel/main.tex)、[中文推导](notes/susskind-cubic-kernel.md)及[来源](literature/susskind-cubic-kernel-sources.md)已保存。任意图归约为177端点模板，独立整数枚举证实所有不足三桥通道相消，由此证明全阶q³探针等价并接已有Mehler核；固定时间小λ得到sech⁶(2Js t)。五组有限λ曲线与精度证据已保存，λ=.005,t1仍偏极限7.09%。下一步raw T₃一次chaos与d2/d3混合高点；无Git提交/推送/ZIP，见[NEXT](NEXT.md)。

## 当前进展：三阶导数 singlet 的独立方向（2026-11-06，本地）

[五页PDF](output/pdf/de-sitter-susskind-cubic-singlet.pdf)、[LaTeX](manuscript/susskind-cubic-singlet/main.tex)、[中文推导](notes/susskind-cubic-singlet.md)与[来源](literature/susskind-cubic-singlet-sources.md)已保存。精确分离 M₃ 的既有时间后代，在第三 coupling chaos 中扣 S、Z 后得到正范数动态 R；μ₂=2Js²(1−q³)/λ，μ₄及九个双边矩至四个外部H均匹配q³探针。下一步证明任意阶归约，接已有完整核；当前未把低阶矩当质量谱。独立精确CAR、3072个标记对照和五页核版完成；继续降低Git频率，无提交、推送或ZIP，见[NEXT](NEXT.md)。

## 当前进展：两 singlet 通道与对易子增长（2026-11-06，本地）

[六页PDF](output/pdf/de-sitter-susskind-two-matter.pdf)、[LaTeX](manuscript/susskind-two-matter/main.tex)、[中文推导](notes/susskind-two-matter.md)与[文献对应](literature/susskind-two-matter-sources.md)已保存。把领先对易子曲率1280分成256的H投影与1024的第一相对通道，独立完整Wick验证通过；构造两物质正Q及标准相对坐标，实现中自由质量也混合k通道。下一步回到独立导数fermion bilinear谱矩阵，区分单介子种类与两个既有X的乘积。41项精确核对、三组高精度Gram及六页核版完成；无Git提交或推送，见[NEXT](NEXT.md)。

## 当前进展：反对称弦后代与正生成元（2026-11-06，本地）

[六页PDF](output/pdf/de-sitter-susskind-antisymmetric.pdf)、[LaTeX](manuscript/susskind-antisymmetric/main.tex)、[中文推导](notes/susskind-antisymmetric.md)和[文献对应](literature/susskind-antisymmetric-sources.md)已保存。从既有 chord 代数证明所有固定阶后代的规范化极限，并构造有限 h 下正自伴 Q_h，其初态谱测度趋于上一轮的 Gamma 分布；boost 使用真实 Liouvillian。下一步研究两物质表示、实际 coproduct 和相互作用。精确 Gram、独立高精度核对与六页核版完成；只本地保存，继续降低 Git 频率，见[NEXT](NEXT.md)。

## 当前进展：候选空间平移与质量族（2026-11-06，本地）

[六页PDF](output/pdf/de-sitter-susskind-translations.pdf)、[LaTeX](manuscript/susskind-translations/main.tex)、[中文推导](notes/susskind-translations.md)及[文献对应](literature/susskind-translations-sources.md)已保存。构造共享singlet boost核的任意正质量Poincare表示，得到可区分质量的Bessel位移核及自由两粒子连续谱；有限q种子测试表明LS长度保持boost不能直接认同当前Liouvillian。下一步从左右弦差的双侧生成元与Gram极限构造微观null平移。9个精确算符恒等式、8个高精度核和六页核版完成；无commit/push/ZIP，见[NEXT](NEXT.md)。

## 当前进展：固定时间首阶四点（2026-11-06，本地）

[五页PDF](output/pdf/de-sitter-susskind-first-order.pdf)、[LaTeX](manuscript/susskind-first-order/main.tex)、[中文记录](notes/susskind-first-order.md)与[来源对应](literature/susskind-first-order-sources.md)已保存。应用Lin–Stanford/Streicher已有公式，保留同种singlet三个通道，取得任意固定四时间连通系数及C(t)、κ22(t)的首阶λ闭式。22个精确Taylor系数、72个Laurent相消、12个50位通道比较通过，有限λ实际偏差完整保留。下一步构造独立空间/平移观测量以推进介子质量字典；无commit/push/ZIP，见[NEXT](NEXT.md)。

## 当前进展：Rindler 波包构造与四点轮廓（2026-11-06，本地）

[五页PDF](output/pdf/de-sitter-susskind-rindler-packet.pdf)、[LaTeX](manuscript/susskind-rindler-packet/main.tex)、[中文记录](notes/susskind-rindler-packet.md)和[来源对应](literature/susskind-rindler-packet-sources.md)完成。用已有KL逆变换构造任意正质量的调节体波包，midpoint两点一致误差≤ε，径向范数有闭式；质量与空间尺度的自由度明确。普通等间隔热四点的零实时间极限为33，微观第四矩为3，下一步沿Lin–Stanford已有fake-circle多点映射确定插入方式。91项主核对、19项独立核对及相消精度已保存。仅本地记录，无commit/push/ZIP，见[NEXT](NEXT.md)。

## 当前进展：singlet 的 OTOC 与 scrambling（2026-11-06，本地）

[六页PDF](output/pdf/de-sitter-susskind-singlet-otoc.pdf)、[LaTeX](manuscript/susskind-singlet-otoc/main.tex)、[中文记录](notes/susskind-singlet-otoc.md)及[原文对应](literature/susskind-singlet-otoc-sources.md)已完成。直接调用已有BINT/BNR四点核，保留两个uncrossed及q⁴加权crossed通道，算出完整同种singlet OTOC与对易子平方；小λ长时间接到已有Tricomi scrambling函数，在string时间中增长率4J_s。780项有理通道比较、独立高精度和解析截断界已保存。下一步将这个具体微观四点接到Rindler体波包观测量。只作本地记录，无commit/push，见[NEXT](NEXT.md)。

## 当前进展：多 singlet 联合矩与实时连通四点（2026-11-06，本地）

[七页PDF](output/pdf/de-sitter-susskind-singlet-multimark.pdf)、[LaTeX](manuscript/susskind-singlet-multimark/main.tex)、[中文推导](notes/susskind-singlet-multimark.md)和[文献接入](literature/susskind-singlet-multimark-sources.md)已完成。扣除守恒分量后，任意固定阶联合矩化到HH:q、HX:q²、XX:q⁴的已有多弦规则；得到κ₄(t)=(1+q⁴)G_{q⁴}(t)−2G_{q²}(t)²，并计算五组λ的实时曲线。95组精确联合矩、3584项独立CAR等式及Jacobi/Mehler核对已保存。下一步使用现成crossed四点核及实际Rindler字典。遵照降低Git频率的要求，本轮无commit/push，见[NEXT](NEXT.md)。

## 当前进展：首个 singlet 的完整双标度两点核（2026-11-06，本地）

[八页PDF](output/pdf/de-sitter-susskind-singlet-kernel.pdf)、[LaTeX](manuscript/susskind-singlet-kernel/main.tex)、[中文推导](notes/susskind-singlet-kernel.md)及[来源](literature/susskind-singlet-kernel-sources.md)已完成。任意阶端点交换证明给出“守恒二次通道＋已有r=q²探针核”的完整两点谱，并确定零频权重；小λ接到sech⁴(2J_s t)。4320项精确CAR对照及独立Jacobi/Mehler核校验保存。下一步进入实际Rindler体传播字典与高点相互作用；遵照用户要求，本轮不提交或推送Git，见[NEXT](NEXT.md)。

## 当前进展：首个 singlet 的双标度谱矩（2026-11-05，本地）

[七页PDF](output/pdf/de-sitter-susskind-singlet-double-scale.pdf)、[LaTeX](manuscript/susskind-singlet-double-scale/main.tex)、[中文记录](notes/susskind-singlet-double-scale.md)和[来源](literature/susskind-singlet-double-scale-sources.md)已完成。通过局部Wick行列式与交叠极限，求得中心化singlet的等时强度和第二频率矩；归一初始曲率在随后小λ极限趋16J_s²。一般N,p系数式、独立CAR与全协方差检查已保存。下一步求零频投影、第四频率矩和完整谱；当前仅本地保存，未提交/推送，见[NEXT](NEXT.md)。

## 当前进展：singlet 动态谱矩的对称性化简（2026-11-05，本地）

[五页PDF](output/pdf/de-sitter-susskind-singlet-moments.pdf)、[LaTeX](manuscript/susskind-singlet-moments/main.tex)、[中文推导](notes/susskind-singlet-moments.md)及[来源](literature/susskind-singlet-moments-sources.md)已保存。用标准 Gaussian Wick 方法与已有 U(N) 无重数分解，把首个动态谱矩化为纯三次张量正范数；完整有限系统收缩、独立15配对及打破对称性的求积控制已核对。下一步将张量和压缩为一般 N,p 的组合数表达，求双标度行为。当前为本地阶段稿，尚未提交或发布；继续位置见[NEXT](NEXT.md)。

## 当前进展：规范已统一，singlet谱化为带电四点问题（2026-11-05）

[六页推导PDF](output/pdf/de-sitter-susskind-dictionary.pdf)、[LaTeX](manuscript/susskind-dictionary/main.tex)、[中文记录](notes/susskind-dictionary.md)与[来源](literature/susskind-dictionary-sources.md)完成。直接复用BNR的complex chord谱与embedded-ensemble二阶矩，统一CAR、方差和时间单位。首个候选A₂已写成相邻电荷扇区的四点谱密度；18个样本及独立相干重建完成校验。下一步研究该特定四点对象的大N矩与重求和，见[NEXT](NEXT.md)。

## 当前主线：从 Susskind 静态补丁全息到体粒子（2026-11-05）

研究重心转向 **静态补丁 → DSSYK∞ → 尺度分离 → singlet 动力学 → 体粒子与相互作用**。首先精读已有论证，再求 complex DSSYK 的 singlet 关联矩阵及其体谱映射。

[七页研究稿 PDF](output/pdf/de-sitter-susskind-mainline.pdf)、[LaTeX](manuscript/susskind-mainline/main.tex)、[中文规划](notes/susskind-mainline.md)与[文献地图](literature/susskind-mainline.md)已保存。首步给出导数复合算符的守恒低项、候选 A₂=Σ[H,cᵢ]†[H,cᵢ]，以及它与双时四点核的关系。下一步统一模型规范并推小正λ下的实时核，见[NEXT](NEXT.md)。

以下为阶段历史；表示接口及五边形计算作为独立技术补充保存。

## 最新进展：参数连续性与正规表示接口（2026-11-05）

[四页PDF](output/pdf/de-sitter-representation-interface.pdf)、[LaTeX](manuscript/representation-interface/main.tex)、[中文推论](notes/stage4-representation-interface.md)及[来源核查](literature/stage4-representation-interface-sources.md)已完成。应用已有支配收敛和 Stone–von Neumann 定理，说明上一阶段的非零残差延续到一个未定量的正参数邻域，且在同一 Weyl 关系、相位函数与普通谱演算下不能靠换正规表示消除。它是旧计算辅助证书的解析后果，未独立重算证书，也不涉及有限G物理结论。下一步返回观察者反射与实际内积/缝合字典，见[NEXT](NEXT.md)。

## 最新进展：明确连续表示中的非零五边形残差（2026-11-01）

[十页PDF](output/pdf/de-sitter-residual-certificate.pdf)、[LaTeX](manuscript/residual-certificate/main.tex)、[中文记录](notes/stage4-residual-certificate.md)与[来源核查](literature/stage4-residual-certificate-sources.md)已完成。直接采用已有Gaussian矩、Mehler和严格求积定理，认证两步chirp及五个全平面F重叠；对明确hbar=1连续L²表示和归一Gaussian，得到||(AB−BM⁻¹BMA)f|| > .04489。

21,143叶子、63,551项覆盖/算术/区间核查及独立局部求积对照已保存。结果与原文五边形的字面谱实现有实质张力，下一步独立核其定义及延拓框架，再回有限G物理缝合；不是de Sitter全息否定结论。继续位置见[NEXT](NEXT.md)。

## 最新进展：严格热截断界与尚未认证的矩阵元（2026-11-01）

[六页PDF](output/pdf/de-sitter-heat-matrix.pdf)、[LaTeX](manuscript/heat-matrix/main.tex)、[中文记录](notes/stage4-heat-matrix.md)及[来源边界](literature/stage4-heat-matrix-sources.md)已完成。固定hbar=1和归一Gaussian测试，用挖孔Sobolev近似控制完整相位的热截断误差；纯整数/有理区间证书给ε=1e−10时严格上界0.000982。标准方法不作原创性声明。

37条FFT/Hermite矩阵元、12条独立左侧积分和188项证书/算术核查已完成。矩阵元模约0.0057–0.0059，但求积误差未认证，不能据此宣布连续反例。平直对照也表明热截断本身会破坏恒等式。下一步先核已有严格核公式能否覆盖当前表示，再补真正缺失的右侧积分误差；继续位置见[NEXT](NEXT.md)。以下保留阶段历史。

## 最新进展：有界配对的 Fourier 字典与热核误差（2026-11-01）

[六页PDF](output/pdf/de-sitter-paired-kernel.pdf)、[LaTeX](manuscript/paired-kernel/main.tex)、[中文记录](notes/stage4-paired-kernel.md)和[文献对照](literature/stage4-paired-kernel-sources.md)已完成。核用GK已有Gaussian-group框架，明确当前测度、符号及三个相近构造的适用差别；给完整有界相位的绝对可积截断核和Gaussian矩阵元误差式，避免拆分无界因子。

240条嵌套Sobol记录、30项字典/控制检查及1290项保存算术/哈希重放完成。热误差的数值仍为经验估计，尚未计算或认证截断五边形矩阵元。下一步直接推进该矩阵元及其严格误差；完整恒等式、物理缝合和全息对偶仍开放，以[NEXT](NEXT.md)为准。

## 最新进展：字面弱积分与移动极点（2026-11-01）

[五页PDF](output/pdf/de-sitter-weak-continuation.pdf)、[LaTeX](manuscript/weak-continuation/main.tex)、[中文记录](notes/stage4-weak-continuation.md)及[来源核查](literature/stage4-weak-continuation-sources.md)已完成。应用已有 Cauchy 核分布恒等式，明确一个固定紧测试的字面弱积分具有非零参数源项，不能自动等同于全纯延拓。例子在第一象限内部；未测试最终有界配对五边形，也不是 de Sitter 全息否定结论。

13条主记录、3个归一正态差分、34项诊断及62项算术/哈希重放通过；数值误差为经验比较。下一步先查直接适用于有界配对函数的积分/Fourier 核与严格恒等式，逐项核对文献假设；继续位置见[NEXT](NEXT.md)。下列条目保留阶段历史。

## 最新进展：定义域、相位奇性与文献证明条件（2026-11-01）

[六页PDF](output/pdf/de-sitter-flip-domain.pdf)、[LaTeX](manuscript/flip-domain/main.tex)、[中文推导](notes/stage4-flip-domain.md)与[证明依赖图](literature/flip-proof-dependency-map.md)整理了已有理论的适用条件。在明确L²(R²)表示中，形式对易不推出强对易；单因子极点阻止Gaussian入域，配对虽有界却使该态失去H¹正则性。移动极点还排除了字面单因子族在参数开集上的共同稠密域。这些诊断不是有界五边形反例，也不作数学原创性声明。

11组过采样Galerkin计算及每步向量已保存，322项算术/哈希重放通过；原22项诊断的2项失败与另5项refinement完整保留，新增31项局部正则性控制通过。正margin仍未认证求积误差。下一步核弱/轮廓延拓与配对谱算符的确切关系，继而回到实际观察者插入、质量和缝合字典。完整全息对偶仍开放，当前以[NEXT](NEXT.md)为准。

## 最新进展：正内积局部表示与未解决的五边形复现（2026-11-01）

[五页PDF](output/pdf/de-sitter-unitary-flip.pdf)、[LaTeX](manuscript/unitary-flip/main.tex)、[中文记录](notes/stage4-unitary-flip.md)保存了对Kim–Scarinci已有有限形变量子化的具体复现。明确四穿孔球的正Hilbert表示、相位配对与非零中心值的代数扩展；尚未建立观察者质量字典或物理缝合。

21组统一态FFT/Hermite计算中，平直五边形对照达到约1e-12/1e-15，但正宇宙常数复现仍有非零残差。36项标量/代数检查和161项证据算术核查通过，不代表连续五边形通过或已被反驳。下一步按稿中投影泄漏界控制连续误差；完整研究仍开放。当前以[NEXT](NEXT.md)为准。

## 最新进展：四段质量滤波与相关时间（2026-10-31）

[五页PDF](output/pdf/de-sitter-mass-filter.pdf)、[LaTeX](manuscript/mass-filter/main.tex)、[中文推导](notes/stage4-mass-filter.md)完成。复读CSTY确认原文已有四段质量滤波；精确变换给相关Gaussian时间，中心时间方差是旧独立涂抹的一半。相应小宽度首阶修正可以反号，这不意味着引力增长项反号。

64组条件比较、16组指定质量平均、23项诊断与444项原始算术重放通过。质量平均保持G固定，包含断开归一权重；有限κ映射仍是比较模型的额外假设。下一步核实际插入、Gram缝合及反射正性，完整物理迹与全息对偶尚未建立。当前继续位置见[NEXT](NEXT.md)。

## 最新进展：实 principal 探针与带相位四点（2026-10-31）

[六页 PDF](output/pdf/de-sitter-principal-eikonal.pdf)、[LaTeX](manuscript/principal-eikonal/main.tex)及[中文推导](notes/stage4-principal-eikonal.md)完成。保留真实principal Bessel波函数，通过原动量的相反旋转确定两侧解析延拓；Gaussian平均的各阶矩有正谱表示，但完整四点仍有相位与相消。数值扫描显示涂抹后虚部可以变号，最晚点符号未分辨，不能沿用共形全时间负号定理。

144次完整scrambled积分、60项公式/独立求积检查、14项原复轮廓/改变重要性检查及727项原始算术重放完成。误差为经验估计，未认证完整引力误差。下一步回到实际制备、有限κ插入归一和Gram缝合；完整全息对偶仍开放。旧条目保存历史，当前从[NEXT](NEXT.md)继续。

## 最新进展：四时间涂抹与有限引力 eikonal 接口（2026-10-31）

[六页 PDF](output/pdf/de-sitter-gaussian-eikonal.pdf)、[LaTeX](manuscript/gaussian-eikonal/main.tex)和[中文推导](notes/stage4-gaussian-eikonal.md)已完成。在已有双共形探针近似中，四个独立Gaussian涂抹化为三个时间变量，必须同时改变断开两点基准。证明重观察者条件下两条候选轮廓的非零虚部在涂抹后仍存在；有限时间分辨率本身没有解决轮廓实性问题。不是完整有限G或全息否定定理。

8组主求积及加阶，原四坐标独立计算，52项诊断和11项进一步检查通过；最大32/48阶差约1.38e−9。下一步精读观察者边缘模/引力零模来源，并推进真实principal探针的Bessel核；本轮共形complementary特例不替代上一轮全Fock模型或最终目标。


## 最新进展：受调节 Fock 场与完整四点（2026-10-31）

[九页阶段稿 PDF](output/pdf/de-sitter-ds3-fock-vacuum.pdf)、[LaTeX](manuscript/ds3-fock-vacuum/main.tex)及[中文推导](notes/stage4-ds3-fock-vacuum.md)已完成。固定 s,p,σ>0、G=0，完整 CG 重数空间上的受调节创建块给全有限粒子数共同域与本质自伴场；补上真空返回后得到完整四点基准。真实 Gaussian 时间涂抹保留重质量极限的热详细平衡与非零自由场对易子；半周期创建加伴随是另一种 observable，不能混同。

20组半周期与16组 Gaussian 双截取求积、36与24项诊断、154与115项独立复核通过。原序列化失败与符号结构比较失败保留；没有新的完整四场全谱数值求积。下一步回读有限 G 插入、态与缝合测度，再比较微观全息字典；完整对偶仍未建立。以下旧条目保存历史，当前以 [NEXT](NEXT.md) 为准。


## 最新进展：任意紧转移与两次创建的 Bose 极限（2026-10-31）

[七页阶段稿 PDF](output/pdf/de-sitter-ds3-compact-global.pdf)、[LaTeX](manuscript/ds3-compact-global/main.tex)及[中文推导](notes/stage4-ds3-compact-global.md)已完成。用 cosh 乘积扩大核的紧转移窗口，再由实际两链截取的 Hilbert 范数尾控制去掉截断；在 G=0、固定观察者下界、半周期调节、紧时间的两次创建模型中得到完整 Bose 自由观察者基准。18个新完整核点、原失败参考及精确标签/加精度修复全部保留，444项独立核查通过；未做新的全物理积分数值认证。

下一步是完整 Fock 共同域与四场中的真空收缩，然后继续有限 G 和微观全息比较。当前是模型内的自查推导，完整 de Sitter 全息对偶仍开放。以下旧阶段进度保留当时状态，当前以 [NEXT](NEXT.md) 首段为准。

这里是 Alan Day 与澄衡的独立研究记录。我们关心：宇宙内部的有限资源观察者，如何获得自洽的量子可观测量，以及这些量如何与视界熵、反冲和引力动力学联系起来。完整 de Sitter 全息对偶仍是长期开放目标。

公开仓库：[wkryyds4-netizen/de-sitter-holography](https://github.com/wkryyds4-netizen/de-sitter-holography)。这里保存工作稿及可复算证据；私人约定与过程日志留在本机，上传范围见 [REPOSITORY](REPOSITORY.md)。

## 当前进展（2026-10-31）

最新[真实谱权重与局部配对 PDF](output/pdf/de-sitter-ds3-physical-local.pdf)共4页，附[LaTeX](manuscript/ds3-physical-local/main.tex)和[中文推导](notes/stage4-ds3-physical-local.md)。实际创建幅满足固定紧窗的统一C²界，接上局部物理积分；97项诊断及165项独立Decimal重放通过。下一步扩大核定理的紧转移适用范围，再控制物理谱尾。完整Bose自由极限、有限G和完整全息对偶仍未完成。

以下为各阶段原有入口及历史进展，旧“运行中/未冻结/下一步”以[NEXT](NEXT.md)首段为准；旧科学文件保持原字节。

## 从这里开始

- [跨零统一极限与局部测试积分 PDF](output/pdf/de-sitter-ds3-uniform-local.pdf)：5页，双方署名；[LaTeX](manuscript/ds3-uniform-local/main.tex)。给穿过r=0的一致O(1/M)误差与局部dr分布极限。21组跨零计算、11个加精度比较、144个完整核积分节点及446项核查通过。真实物理谱积分仍开放。[Git保存范围](REPOSITORY.md)说明本机私人记录与第三方缓存的排除规则。

- [完整核横向相位 PDF](output/pdf/de-sitter-ds3-transverse-phase.pdf)：4页，双方署名；[LaTeX](manuscript/ds3-transverse-phase/main.tex)。得到去零紧窗口的带相移正弦式及完整核逐点sinc(λ)极限；16组计算、203项独立核查、4个加精度检查通过。旧100位验证已完成，以下旧阶段“未冻结/仍运行”文字保留当时状态，当前以NEXT首段为准。

- [完整交叉核峰面极限 PDF](output/pdf/de-sitter-ds3-peak-limit.pdf)：5页，Alan Day 与 Chengheng 署名；[LaTeX](manuscript/ds3-peak-limit/main.tex)。保留完整物理相位，给出 logM/π、常数、1/M 修正和紧转移一致余项；12组计算、修复后191项审计通过。仅精确峰面，未证明波包积分极限；尚未冻结。

- [Wilson对偶与稳定交叉核 PDF](output/pdf/de-sitter-ds3-shadow-duality.pdf)：4页，双方署名；[LaTeX](manuscript/ds3-shadow-duality/main.tex)。用已知Wilson对偶建立Omega反射式，解析消去原shadow系数的指数误差放大；18组完整核计算、365项独立审计，旧失败与高精度修复保留。完整波包重质量极限仍开放，尚未冻结。

- [完整尾界与shadow指数放大 PDF](output/pdf/de-sitter-ds3-full-tail.pdf)：4页，双方署名；[LaTeX](manuscript/ds3-full-tail/main.tex)。给出完整核截断界及exp(2πM)误差放大机制，407项独立重放检查通过；失败近似全部保留。尚未证明完整重质量极限，尚未冻结。

- [可去奇点与重质量因子 PDF](output/pdf/de-sitter-ds3-collision.pdf)：5页，Alan Day与Chengheng署名；[LaTeX](manuscript/ds3-collision/main.tex)。当前未冻结阶段稿；额外100位独立核对仍运行，最终审计未完成，继续位置见NEXT。

- [同种粒子与Bose归一 PDF](output/pdf/de-sitter-ds3-bose.pdf)：5页，双方署名；[LaTeX](manuscript/ds3-bose/main.tex)。

- [真实交叉波包与谱下界见证 PDF](output/pdf/de-sitter-ds3-crossed-packet.pdf)：5页，双方署名；[LaTeX](manuscript/ds3-crossed-packet/main.tex)。

- [交叉重耦合物理字典 PDF](output/pdf/de-sitter-ds3-crossing-dictionary.pdf)：5页，Alan Day、Chengheng署名；[LaTeX](manuscript/ds3-crossing-dictionary/main.tex)。

- [dS₃两物质与非交叉极限 PDF](output/pdf/de-sitter-ds3-two-matter.pdf)：6页，Alan Day、Chengheng署名；[LaTeX项目](manuscript/ds3-two-matter/main.tex)。

- [dS₃群配对与物理两点 PDF](output/pdf/de-sitter-ds3-group-pairing.pdf)：7页，Alan Day、Chengheng署名及AI助手脚注；[LaTeX项目](manuscript/ds3-group-pairing/main.tex)。

- [dS₃质量谱与BD投影 PDF](output/pdf/de-sitter-ds3-mass-spectrum.pdf)：5页，Alan Day、Chengheng署名；[LaTeX项目](manuscript/ds3-mass-spectrum/main.tex)。

- [dS₃局部插入与CG系数 PDF](output/pdf/de-sitter-ds3-local-vertex.pdf)：6页，Alan Day、Chengheng署名；[LaTeX项目](manuscript/ds3-local-vertex/main.tex)。

- [全局重耦合与归一不变检验 PDF](output/pdf/de-sitter-global-racah.pdf)：5页，Alan Day、Chengheng署名；[LaTeX项目](manuscript/global-racah/main.tex)。

- [缝合与恒等极限 PDF](output/pdf/de-sitter-recoil-identity.pdf)：5页，双方署名；[LaTeX项目](manuscript/recoil-identity/main.tex)。

- [两点谱校准与归一化自由度 PDF](output/pdf/de-sitter-recoil-spectral.pdf)：5页，双方署名；[LaTeX项目](manuscript/recoil-spectral/main.tex)。

- [夹逼极点消去与两点校准 PDF](output/pdf/de-sitter-recoil-poles.pdf)：5页，双方署名；[LaTeX项目](manuscript/recoil-poles/main.tex)。

- [两手征配对与外部接口 PDF](output/pdf/de-sitter-recoil-pairing.pdf)：6页，Alan Day、Chengheng署名；[LaTeX项目](manuscript/recoil-pairing/main.tex)。

- [横向质量与中心局域化 PDF](output/pdf/de-sitter-recoil-transverse.pdf)：6页，双方署名；[LaTeX项目](manuscript/recoil-transverse/main.tex)。恢复u/v的联合估计，并说明两支横向曲率相反。
- [研究稿 PDF](output/pdf/de-sitter-observer-foundations.pdf)：55 页，署名 Alan Day 与 Chengheng，含目录、证明、四点函数与两个附录；[LaTeX 主文件](manuscript/main.tex)和[编译检查记录](results/manuscript/build-review.md)可用于继续修改。
- [两分支联合估计 PDF](output/pdf/de-sitter-recoil-bessel-joint.pdf)：6 页，双方署名；[LaTeX 项目](manuscript/recoil-bessel-joint/main.tex)。保留整体相位，建立对称内部图的完整领先联合误差。
- [联合端点与中央相消 PDF](output/pdf/de-sitter-recoil-bessel-endpoints.pdf)：5 页，双方署名；[LaTeX 项目](manuscript/recoil-bessel-endpoints/main.tex)。两端已得到包络控制，中央段仍须保留相位。
- [相消零点与联合正尾 PDF](output/pdf/de-sitter-recoil-bessel-cancellation.pdf)：6页，双方署名；[LaTeX项目](manuscript/recoil-bessel-cancellation/main.tex)。
- [完整固定K极限与首阶修正 PDF](output/pdf/de-sitter-recoil-bessel-limit.pdf)：5页，双方署名；[LaTeX项目](manuscript/recoil-bessel-limit/main.tex)。
- [全局轮廓与两分支Bessel PDF](output/pdf/de-sitter-recoil-global.pdf)：5页，双方署名；[LaTeX项目](manuscript/recoil-global/main.tex)。
- [作用量恒等式与局部积分 PDF](output/pdf/de-sitter-recoil-local.pdf)：5页，双方署名；[LaTeX项目](manuscript/recoil-local/main.tex)。
- [反射鞍点与极点振荡 PDF](output/pdf/de-sitter-recoil-saddle.pdf)：4页，双方署名；[LaTeX项目](manuscript/recoil-saddle/main.tex)。
- [固定质量内部积分 PDF](output/pdf/de-sitter-recoil-internal.pdf)：4页，双方署名；[LaTeX项目](manuscript/recoil-internal/main.tex)。
- [纯反冲合并顶角 PDF](output/pdf/de-sitter-recoil-vertex.pdf)：4页，双方署名；[LaTeX项目](manuscript/recoil-vertex/main.tex)。
- [Gram与横向涨落 PDF](output/pdf/de-sitter-gram-prefactor.pdf)：5页，署名 Alan Day 与 Chengheng；[LaTeX项目](manuscript/gram-prefactor/main.tex)。
- [三维WKB比较 PDF](output/pdf/de-sitter-wkb-overlap.pdf)：4页，给重探针相对窗口、临界缩放及振幅检验；[LaTeX项目](manuscript/finite-g-bridge/main.tex)。
- [研究规划](PLAN.md)：方向依据、阶段成果与推进条件。
- [下一步](NEXT.md)：当前问题、阅读顺序和验收标准。
- [研究日志](journal/2026-10-31.md)：决策、尝试、失败与修正；历史段落保留当时状态，当前进度以本页和 NEXT 为准。
- [观察者与代数](literature/observers.md)、[候选微观模型](literature/models.md)：领域发展与选题依据。

## 当前研究到哪里

最新[跨零局部极限](notes/stage4-ds3-uniform-local.md)通过先合并Gamma极点再取导数，消去旧1/|r|界的假发散；利用两支在零点的共同O(M⁻²)误差，得到完整核对相移正弦式的一致O(1/M)近似。C²紧测试的dr积分趋于f(0)，速率O_f(1/logM)+O(||f||₁/M)。这尚未纳入实际创建幅、质量Jacobian和非紧尾。下一步核对真实物理权重能否满足局部定理的测试族条件。

最新[横向相位推导](notes/stage4-ds3-transverse-phase.md)说明物理归一相位使单因子的sinc(λ/2)变为完整核的sinc(λ)。在固定非零r上得到U=sin(rlogM−δ(r))/(πr)+O(1/M)，并给Gamma系数精确保留时的有限M界；穿过r=0的联合余项及物理积分仍未证明。90/120位复核保留慢而不单调的收敛。旧47841正常退出，12组比较及collision独立审计通过，当前没有旧活跃计算句柄。

最新[峰面完整核推导](notes/stage4-ds3-peak-limit.md)在指定 G=0、ℓ=1、BD/KG 与半周期创建模型中，解析证明两系数共轭关系和完整归一相位，控制大参数伴随级数整体和，得到精确峰面的重质量展开。初始审计因 π 精度不足失败，原始文件保留；独立 Machin 重放191项通过，五页 PDF 编译及逐页核版完成。下一步横向邻域、振荡尾与真实创建幅上的积分极限；旧100位验证仍运行，四个新增阶段尚待整体登记和冻结。下列条目保留历史阶段状态，当前继续位置以 NEXT 首段为准。

最新[同species Bose研究](notes/stage4-ds3-bose.md)从Fock创建和物理群平均推出√2P+及直接项加交换项，证明同时间次序差严格为零；对指定紧初态，球面有限K见证证明D_B(t)=Q_Bt²+O(t⁴)、Q_B>0。三组完整网格共1,216,480点，独立E/W/ρ重放通过，Q_B数值约3.17558e−5。保留一个低质量近重合点160%相对精度失败及55/75位修复；网格加密差大于谱尾，总误差未认证。下一步稳定完整交叉的重质量极限及Bose自由基准，继而全Fock、有限G；完整对偶仍开放。

最新[真实交叉波包](notes/stage4-ds3-crossed-packet.md)已为明确Beta(5,5)初态完成四重质量积分，保留完整复相位及无限谱尾控制。五组raw网格共3,442,832个kernel点，零时刻调节后次序差平方范数约6.50724e−7，两个细网格的交叉复差约1e−11；L12质量谱截断绝对界6.26346e−10。另由真实顶角的odd-convolution表达得到缺失[0,s)谱段的确定符号，并用有限K群平均测试证明此模型D(0)>0。24项解析/见证诊断及全部网格独立E/W重放通过；五页PDF完成双编译与视觉核查。远端严重相消点约1.97%相对误差保留，整体数值未作总误差认证。下一步同species Bose，继而一般制备/重质量及有限G；完整对偶仍开放。

最新[交叉字典](notes/stage4-ds3-crossing-dictionary.md)已将实际局部创建链接到LPRS标准6j，固定000球周期、shadow反射、32π³ quotient比例与物理U_00的相位。独立Euclidean复平面积分在两组40阶样本的复相对差约8.7e−11和4.6e−8，另给三个puncture区域的保守解析尾界。33+10项诊断通过；五页双方署名PDF双编译及逐页核查完成。实际四重波包积分尚未计算，不能把scalar子块当完整酉变换；下一步沿已明确的创建振幅和正范数尾控制计算真实交叉态重叠。有限G和完整对偶仍开放。

最新[两物质推导](notes/stage4-ds3-two-matter.md)在明确的G=0、两不同物质species模型中，以完整CG理论确定局部观察者所达标量通道的四体逆范数，建立1→2闭块及二阶质量矩充分域。真实静质量半周期延拓后创建块有界；固定下界、条件外质量M→∞的非交叉四点恢复独立反极探针乘积，并给显式绝对界与i/M相位修正。56+10项诊断通过，含非球面群积分；六页双方署名PDF编译并逐页核查。交叉全自旋重耦合、同species Bose、全Fock及有限引力仍开放。

最新[群配对与两点推导](notes/stage4-ds3-group-pairing.md)在明确G=0局部模型中，用SO⁺(3,1)群平均确定ρ_G=2ν²/π与三体逆范数，线性左右配对选出sech(πν)态，得到完整单物质角不变扇区的正核Q=μ(a)μ(b)/(8D)。建立0↔1闭块及真实静质量相位两点，并在统一大质量规范后得到旧recoil核的n_loc=[(1+e⁻²πa)(1+e⁻²πb)]⁻¹。群密度与sech态均有明确前驱，未作原创声明。81+16项诊断通过，七页署名PDF经双编译及逐页核查。完整Fock/复合gluing、有限κ和有限G仍开放，下一步是两物质四点扇区。

本轮[dS₃质量谱与投影](notes/stage4-ds3-mass-spectrum.md)以Darboux极分解给完整连续双通道、束缚态及零阈值共振，推得时空Gram和μ(ν)=νtanh(πν)/π。BD与其反极像的归一重叠为−sech(πν)，真正正交补另有其线性组合。得到保持proper群的观察者投影、两腿重标及有界压缩插入的两点恒等式。150项诊断通过，初始短时盒两项失败完整保留；五页双方署名PDF编译及视觉核查完成。物理规约内积、态、无限物质模域及有限G仍开放。

本轮[dS₃局部插入](notes/stage4-ds3-local-vertex.md)从BD/KG归一模式得到收敛三点时空积分、有限Beta求和、future/shadow字典及完整复CG系数。球面周期以Funk–Hecke和Dougall独立核定规范；与旧C₀配对的差别是显式腿因子。249项诊断及24项比较通过，六页双方署名PDF完成编译和逐页检查。这是G=0指定局部模型结果；观察者质量谱、物理两点n、有限G缝合及完整对偶仍开放。

本轮[全局SL(2,C)字典检验](notes/stage4-global-racah.md)核清主系列权重、shadow反射和含自旋求和的测度。对明确的拓扑保持标量映射，以四点矩形比值消去全部局部行/列归一后，H和配对图Z与已知全局6j仍分别相差约5.06%和4.59%；加精度与独立Euler积分支持这一数值反例，尚非区间认证。另给非负连续核不能独立构成完整标量酉变换的正性证明。44项诊断通过。下一步从G=0 dS₃观察者局部场矩阵元推导完整三标签系数，不再仅凭图形相似猜字典；有限G与完整对偶仍开放。

最新[缝合与恒等极限研究](notes/stage4-recoil-identity.md)核清实b理论的块测度及恒等腿分布；CLP共同角路径的两权重和恒为1，ν=±i不是双真空。约化两点正规化后有严格δ极限、显式跨极点补项及有限M误差界，但全部nλ仍能通过该端点检验。因此不以identity校准宣称归一唯一。106项诊断通过（包含重复真空标签核对）。下一步核G=0全局SL(2,C) partial-wave的反射、内积和测度，再与反冲H建立或否定完整字典；物理有限κ缝合、有限G和完整对偶仍开放。

最新[两点谱校准](notes/stage4-recoil-spectral.md)把直接C₀和观察者权重组装为显式有限质量因子R，区分固定平均质量、固定初态及阈值以上迹积分三种协议。构造nλ族证明领先QFT匹配仍不唯一确定有限质量归一；在明确n=1模型中计算了完整谱与阈值t⁻²尾。71项诊断通过。下一步核实际缝合、谱测度和identity-leg分布极限，并检独立自由场的断开四点基准；物理内积、有限G和完整对偶仍开放。

主线是**有限制备怎样影响观察者可测量的关联函数**。在 dS₂、G→0 的指定模型中，正规化有限质量包的交叉四点计算已完成；指定弱实现的算符基础已形成记录，交叉四点显式渐近窗口已形成自查证明草稿；当前转向三维有限G比较所需的共同参数域与完整振幅。

前期通过完整环境谱测试及分数阶径向矩估计，把**明确弱正频测试定义下的最大有限FP空间**与真实种子生成闭包、完整谱空间识别。LaTeX 稿汇集了各有限粒子数的识别、共同群cutoff下的Gaussian弱场与本质自伴性，并把有限包四点接到指定压缩闭场的有序乘积上。

这个结果限定于dS₂、G→0、ν>0、σ>0、实际observer下界s*>0及proper群。这里给出的是明确测试拓扑的弱实现，不代表任意形式求和规则下的最大空间。创建块已有明确种子测试、任意物理输入上的绝对位置配对。新增第7节证明：湮灭的有限模式位置积分在固定测试下以物理对偶范数收敛到伴随；在稠密种子输入上另有Fejér点极限及绝对未截断积分。任意输入的点值推广仍开放，Gaussian和observer下限保留。下一主线转向交叉大质量及增长时间极限；有限G和完整全息对偶尚未建立。

此前[第9节推导](notes/stage3-crossed-limit-reduction.md)区分固定质量下界与随制备移动的下界，证明正范数谱尾的质量一致控制。当时精确交叉核的紧测试分布极限（9.8）尚未证明，命题9.2以它为假设；本轮第11节补齐了该步骤。独立比较模型中，移动边界产生严格正的静态对易子贡献，六组制备的19项数值诊断通过；这不是精确6j收敛验证，也不是Lyapunov增长。原文与障碍见[文献记录](literature/stage3-crossed-large-mass-sources.md)。

此前[第10节](notes/stage3-hypergeometric-endpoint.md)从精确Euler积分证明单个₄F₃因子的两支展开，紧参数一致余项O(Λ⁻¹⁺ε)，包括δ=0的端点合并。第二支在纯虚参数下不衰减，须在四项组装时处理。当时15项有限点诊断通过，完整四支抵消及(9.8)留待下一步。

此前[第11节](notes/stage3-assembled-fixed-time.md)给出新增端点支的解析抵消、无相反标签极点的U±，并以小复管及最大模原理证明合并余项O(Λ⁻²⁺ε)。恢复外腿与精确能量Jacobian后，得到(9.8)，接上两种阈值协议的固定时间有限包极限。移动阈值保留正的静态边界贡献，固定阈值的重复时间对易子极限为零。21项代数诊断通过，另存8个精确核点和4个领先表达测试积分；这些不替代证明。该成果为作者自查的证明草稿，随后接续对数增长时间。

此前[第12节](notes/stage3-logarithmic-time.md)在t_M=log(2M)+T下给delta加PV局部极限，用对所有实时间有效的谱尾界恢复F/Q两种协议及任意固定L¹概率外包。先M→∞再T→±∞，早期回到固定时间平台，晚期交叉项趋零；移动阈值的早期静态D_edge仍保留。精确Gamma系数核对定位了KL主文/附录的一单位时间零点差异，未称作者确认勘误。25项系数及领先表达求积诊断通过，但不是物理包或精确6j的积分认证。当时保留真实制备的曲线与早期速率任务。

此前[第13节](notes/stage3-preparation-asymptotics.md)给完整外腿后的Gaussian横向分离。固定阈值F的早期对易子范数为2C exp(2T)+O(exp(hT))，2<h<4，C>0；该指数已有KL前驱，本轮保留了指定涂抹。移动硬阈值Q对π₀(r)=β exp(−βr)则出现角点：减去静态D_edge后，早期领先为2B/|T|³，B>0。结论均先M→∞再T→−∞，当前为作者自查证明草稿。37项诊断通过，F真实极限曲线已计算，Q只核局部角点，完整曲线及其他密度待继续；有限M窗口、有限G和完整对偶仍开放。

此前[第14节](notes/stage3-moving-threshold-curve.md)直接计算指数外包Q的完整复值极限曲线，按真实max边界分区，并独立加密横向与Fourier求积。静态D_edge约3.88408×10⁻⁴；T=−128的动态差分约1.09044×10⁻¹⁰，与逆三次领先项的比值约.96205。122项经验诊断通过，小信号最大相对加阶变化约1.9×10⁻⁷；这不是总误差认证。数据、曲线和实际频谱保存，下一步转向密度端点消失/支撑间隙及有限M窗口。

此前[第15节](notes/stage3-density-endpoints.md)处理整数阶Gamma密度及固定支撑间隙。复关联的首项在奇数消失阶数时纯虚，实对易子要看下一阶；β=1的p=1、2给负的早期五次幂修正，因此D_edge是早期平台而非普适下界。移位指数密度给振荡四次幂尾。三组完整积分、复频谱及384项经验检查保存；最远端小实部精度损失明确列出。当前为作者自查证明草稿；下一主线转向有限M窗口，有限G和完整对偶继续开放。

此前[第16节](notes/stage3-two-branch-window.md)保留两支及真实硬投影，得到领先包的镜像恒等式。由真实能量紧盒比较与正范数尾界，证明全时间轴的一致绝对逼近；其增长盒常数尚未定量。Gamma幂律尾的快速/慢速比趋于(s/k)^n，s为log时间原点的提前量，k=2log(2M)−s；故固定比例提前量下不能相对忽略快速支。F指数尾的对应比为exp(−4t)。300项代数与实际领先包求积诊断通过，均非精确有限M积分认证。下一步定量增长盒常数，建立可计算的相对窗口；有限G及完整全息仍开放。

本轮[第17节](notes/stage3-growing-box-windows.md)追踪增长盒常数，得到指数乘多项式包络。在外包具有正指数矩的条件下，精确有限M与两支表达的全时间误差比任意固定逆log幂更小，因而s=sqrt(log(2M))等显式渐近窗口成立；固定比例提前量仍须保留两支。290项有限点诊断通过。结果为作者自查证明草稿，未认证具体百分比精度所需的最小质量。下一主线恢复具有共同适用条件的有限G/微观比较，数值常数及精确6j包认证保留为验证任务。

[WKB重叠笔记](notes/stage4-wkb-overlap.md)已给CSTY约化纯反冲₃F₂到₁F₂的显式相对界：Q=o(ν^(2/3))时比值趋一，Q=cν^(2/3)时趋于exp(c³/3)。约化作用量逆变换与反冲函数在t=0分别归一后仍相差趋向cosh²t的prefactor。42项数值检查通过；这是自查的约化模型推导，不是有限G对应。下一步从完整6j、顶角、外腿和谱权重推导领先振幅，跟踪外部轮廓及退化极限。

[Gram与横向涨落笔记](notes/stage4-gram-prefactor.md)给出精确Gram因式分解和等探针经典横向Hessian。CSTY反极点时间关联函数须对两个额外能量方向积分；仅取对称截面遗漏其高斯因子。在明确一圈延拓假设下，Gram与横向因子的乘积趋向此前所需的1/(1+y²)。83项检查通过。这是局部涨落的条件性解释，顶角、物理外腿、外部轮廓及退化一致性仍须处理，尚未建立完整有限G匹配。

[纯反冲顶角笔记](notes/stage4-recoil-vertex.md)从完整C和V精确约分，用Barnes格点三阶导数及平衡矩给固定质量纯反冲极限。合并顶角化为普通Gamma/双曲函数，四顶角在指定横向盒内归一后趋1，附显式指数修正界。42项检查通过。该因子属于内部积分表示，不能直接再乘到完整6j的Gram近似上；下一步转向同一表示的内部double-sine积分、轮廓及涨落。物理四点归一化和有限G匹配仍未完成。

[内部积分笔记](notes/stage4-recoil-internal.md)把CLP完整内部double-sine积分在固定正实质量、严格极点间隙下退化为Barnes积分。精确q-Gamma表达与两端统一可积界允许交换κ极限和内部积分；结合顶角得到整个固定拓扑图的归一纯反冲极限。34项检查通过，含极点合并样本及严重相消的四项₄F₃比较。下一步保留u、v做大质量鞍点与涨落；物理外部能量循环、插入归一化和有限G匹配仍开放。

[反射鞍点笔记](notes/stage4-recoil-saddle.md)保留u、v，从同一Barnes表示得到二次方程、Δ=−det(CM)/32及内部局部行列式，避免重复顶角。经典作用量与横向Hessian通过有限点核对，一般角度恒等式尚未证明。精确反射因子在部分重质量尺度不逐点趋一；W^(4,1) Fourier引理给出与极点间隙无关的积分界，因此不能据逐点误差否定整个窗口。57项诊断通过，完整轮廓、真实局部导数界及物理外部接口仍开放。

[作用量与局部积分笔记](notes/stage4-recoil-local.md)将一般六质量的作用量匹配提升为解析恒等式：对数参数比给外二面角，齐次性恢复CSTY作用量。指定虚能量和保留横向盒内，真实Gamma商的四阶导数界与非抵消高斯下界已建立，极点因子的局部积分修正为O(K⁻²)，即使极点间隙趋零也一致；局部整体相对误差O(K⁻¹)。50项诊断通过。下一步必须处理完整积分与局部贡献之差，不能把局部定理当全局或物理OTOC结果。

[全局反冲与Bessel笔记](notes/stage4-recoil-global.md)保留负方向精确反射，给出紧集核退化与精确双分支Bessel辅助积分。完整轮廓的有限参数数值在y=2显著偏离右局部首项；由紧集展开推得的首阶相位修正改善辅助模型比较。34项诊断通过。下一步转为两分支与质量一致尾界，不再预设全部有界y的右鞍点主导；全Barnes→Bessel极限尚未证明。

[完整固定K极限笔记](notes/stage4-recoil-bessel-limit.md)用精确Gamma递推建立正水平尾的质量一致x⁻³/²界，并用固定实部Gamma估计控制负竖直尾。完整Barnes→双分支Bessel极限及积分首阶O(ν⁻²)绝对余项在固定K、α下得到证明草稿。69项诊断通过；联合K增长、横向变量与有限G物理字典仍开放。

[相消与联合正尾笔记](notes/stage4-recoil-bessel-cancellation.md)分类两支Bessel和的实参数零点，用固定参数的二阶导数项证明这些辅助零点在完整积分中以ν⁻²阶抬升。三阶相位和O(ν⁻⁴)绝对余项已形成证明草稿；联合正尾在ν²≥C_Y K³下可从O(K)基段控制。65项诊断通过。完整联合误差须用绝对或两支模长和，不能在相消处除以首项J。

[联合端点与中央相消笔记](notes/stage4-recoil-bessel-endpoints.md)取 N=ceil[K√(1+y²)/2]，以精确 Gamma 模长控制整个负竖直尾及短连接段；反射与 Stirling 补足正尾基段相对于 E 的界。另一方面，对固定 y≥y_c 的中央窗口，绝对积分/E 至少按 exp[KΔ(y)] 增长，Δ(2)≈.250662。300 项诊断通过，y=2、K=128 的辅助中央绝对范数/E 约 7.16e15。这不是物理振幅增长，说明中央余项必须保留相位；完整联合渐近仍未证明。

[保留相位的两分支联合估计](notes/stage4-recoil-bessel-joint.md)用精确反射拆出负竖线和正权重留数，两支均以Bessel方程消去一阶中心矩，完整有向积分满足 |Iν−exp(iS/ν)J|≤C_Y(K³/ν²+exp(−c_YK))E，覆盖有界α/K及相消点。此前候选ν^(4/3)≪M≪ν²窗口现有对称内部图的领先证明草稿。64项检查通过，初版网格失败保留；下一步恢复u/v及横向尾部，物理外部归一化、有限G与完整对偶仍开放。

| 当前材料 | 作用 |
|---|---|
| [湮灭积分模式极限](notes/stage3-integrated-annihilation.md) | 固定左包传递强收敛；有限截断不假定完整群等变；种子输入有Fejér正核支配 |
| [直接位置创建配对](notes/stage3-direct-position-creation.md) | 有限模式测试侧转置、时间/群/位置联合绝对界及所有物理输入上的算符识别 |
| [完整环境谱测试](notes/stage3-ambient-spectral-test-core.md) | 双连续通道、全部径向束缚态，以及低谱截止在强轨道范数中的收敛 |
| [最大弱正频空间](notes/stage3-maximal-weak-frequency-space.md) | 独立位置测试定义、四个频率sector及实际observer阈值；不对内部CG谱加下界 |
| [真实正频共同域](notes/stage3-positive-frequency-orbit-core.md) | 原模式双端衰减、六轨道整体可积、FP＝GA；另证FP有限截面定义局部可积位置分布 |
| [径向散射](notes/stage3-radial-scattering-core.md) | 非共振、两通道ds谱和自相关L²界，固定角向归一化 |
| [产品GA幅度与稠密性](notes/stage3-product-seed-ga-density.md) | 左质量先积分后的绝对Fubini、固定左包连续界、完整通道稠密性 |
| [有限磁指标延拓](notes/stage3-finite-k-ga-extension.md) | 正确对偶相位、相干交叉项与大小产品域闭包一致 |
| [全有限粒子数空间](notes/stage3-all-n-generated-space.md) | PP/PD归纳、Bose投影及真空独立处理，明确全部J_N |
| [局部弱场及自伴性](notes/stage3-regulated-local-descent.md) | 原模式核经共同群cutoff下降，输出、伴随、√N块界和唯一自伴闭包 |

前置文献与规范见[群平均查新](literature/stage3-group-averaging-update.md)、[完整CG树](literature/stage3-two-particle-cg-map.md)及[运动学创建审查](notes/stage3-creation-descent-audit.md)。此前51项高精度有限点公式诊断通过；这些不是区间认证，也不替代解析证明。

## 已有成果入口

| 阶段 | 研究内容 | 记录 |
|---|---|---|
| 一 | 指定有限谱模型的资源比较、优化与有限读时方差 | [阶段一](notes/stage1.md) |
| 二 | dS₄ 指定自由场谱、共同钟四点、局域化与压缩算符接口 | [阶段二](notes/stage2.md)、[局域化](notes/stage2-localization-bridge.md)、[算符接口](notes/stage2-weyl-bridge.md) |
| 三 | dS₃ 反冲基线；dS₂ 有限包二点、非交叉四点与大质量首修正 | [反冲](notes/stage3-baseline.md)、[有限包](notes/stage3-packets.md)、[非交叉](notes/stage3-fourpoint-study.md)、[速率](notes/stage3-clock-rate.md) |
| 三 | 正规化有限外包的交叉四点与远端尾界；顶角归一化核查 | [交叉报告](notes/stage3-packet-crossing-study.md)、[尾界](notes/stage3-packet-tail-envelope.md)、[顶角](notes/stage3-vertex-analytic-audit.md) |
| 四 | 有限 G 双共形基线、微观时间修正与 TQFT 半经典字典 | [有限 G](notes/stage4-conformal-setup.md)、[时间项](literature/stage4-time-corrections.md)、[TQFT 字典](literature/tqft-recoil-dictionary.md) |

有限外包算例取 ν=.6、σ=1、内部谱[.2,3]、外包 d∈[.7,.9]。四组内外加阶及最终93项诊断通过，15,560个原始点核已保存。内部最大加阶变化约1.17×10⁻⁶；外层差低于内部经验变化，不能作为认证外积分误差。cap20以外的交叉配对尾严格小于2×10⁻¹¹，但盒内及有限区间求积尚无严格总误差条。完整复值、参数和失败历史见报告。

## 文件怎样安排

| 位置 | 内容 |
|---|---|
| `literature/` | 原文版本、公式位置、适用范围与查新 |
| `notes/` | 推导、阶段结果、条件命题与候选路线 |
| `code/` | 复算、验证、来源下载与归档脚本 |
| `results/` | 参数、原始输出、环境及代码和输入哈希 |
| `sources/` | 主清单中的26份来源文件及 SHA-256；不是26篇论文；另有四项在线核读但未成功缓存，以及仅在线核读的JT前驱、DLMF页面、Meng–Yang定理1.3及MTV附录B.1 |
| `manuscript/`、`output/pdf/` | 正式 LaTeX 研究稿、BibTeX 文献与编译后的阅读版 |
| `journal/` | 研究过程、决策与修正 |
| `snapshots/` | 不覆盖的阶段 ZIP 与校验值 |

原始论文版权归各自权利人，缓存用于本地研究复核。项目尚未对外发表，新研究只写本独立目录。

## 保存与恢复

每轮更新报告、PLAN、NEXT 和日志；确认计算与编辑停止后，在本项目根目录执行：

```powershell
python code/project_audit.py
python code/project_audit.py --snapshot --label new-checkpoint-name
```

[审计记录](results/project-audit.json)检查必需文件、本地链接、来源与已绑定代码/输入的哈希；[文件清单](results/project-files.json)记录当前字节数和 SHA-256。它检查记录一致性，不替代数学证明或误差分析，也不独立认证未绑定的 CSV 与图片。

本轮检查点名称为 `2026-10-30-recoil-transverse`。ZIP 写入后逐文件读回核对，校验值另存；既有快照不覆盖，后续必须换新名称。归档包含研究材料、代码、原文缓存、原始结果和失败记录，不递归打包旧快照与临时编译文件。恢复时先用同名 `.sha256` 核对选定 ZIP，再解压到一个空目录，从 README 和 NEXT 接续。
