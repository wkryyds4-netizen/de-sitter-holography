# 双时算符到介子的重建接口

2026-11-06，本地阶段稿；正式证明与公式见 manuscript/susskind-bilocal-interface/main.tex，六页阅读版见 output/pdf/de-sitter-susskind-bilocal-interface.pdf。

## 文献给出的入口

Itakura hep-th/9604032v2 §5 已给出光前双局域 boson 到介子的投影。我们保留其 √P Jacobian 和 L²(dx) 内积，接上一轮 FLZ 的 a=0 质量算符。原文使用 M²=2P⁺P⁻，本项目用 M²=P⁺P⁻，没有直接移植 Hamiltonian 前因子。原文的 canonical B 与原始夸克双线性之间还有非线性 boson expansion；这是构建微观映射的重要接口。

## 一个直接可用的体计算

令 2x−1=cosθ，内部源 f_k=sin(kθ)/√(Bkk/2)。这里的 k 是内部基底序号，没有认同为 SYK 时间导数阶数。使用已存 256 维 FLZ 模式，投影至前六介子；400 点角积分独立核对。

完整矩可直接计算：〈M²〉/τ=πk/Bkk，〈M⁴〉/τ²=4k²Ckk/Bkk，其中 Bkk=1+1/(4k²−1)，Ckk=∫_{−1}¹U_{k−1}(z)²dz。k=1 的最低介子权重 .9947003267，但质量矩比 32/(3π²)=1.08075929；前六级之外仅 .0002089268 的范数，贡献 .2127687761 的四次质量矩。k=2 对应比值 512/(45π²)=1.15280991，尾部四次矩4.09931829/40。不能凭很大的单级重叠把源视为单质量。

## 微观侧的精确接口

采用标准 CAR、β=μchem=0 归一迹、H=2^pΣJ_IJ C†_I C_J，E J_IJ J_KL=v δ_IL δ_JK，v=Js²N/p²(p!/N^p)²。

定义 D_i=[H,c_i]，T_ij={D_i,c†_j}，U_ij={D_i,c_j}。有限N的 unequal-time CAR 首阶为 δ_ij+itT_ij 与 itU_ij。对四个奇算符有精确式

`[ab,cd] = a{b,c}d − c{d,a}b + ca{d,b} − {a,c}bd`。

应用于双时 singlet 时必须保留全部四项。算符值的反对易子不能预先替换成标量。

按 Majorana/Clifford 次数 2k 正交分解 H，w_k=||H_k||²。这个 k 与旧 U(N) 协方差不可约 rank、coupling Gaussian chaos 均不同。中性 2k 次 exterior 向量有 k 个创建与 k 个湮灭形式；一次 contraction 系数 √2，再分别收缩 k 或 k−1 个形式，得到

- Σ_i||D_i||²=2Σ_k k w_k；
- Σ_ij||T_ij||²=4Σ_k k²w_k；
- Σ_ij||U_ij||²=4Σ_k k(k−1)w_k；
- Σ_ij|tr₀T_ij|²=4w₁。

对重叠 r=|I∩J| 展开共同占据 n_i=(1+z_i)/2。选择 s 个 z 的次数为 k=p−r+s；经 2^p 规范化后每项 HS 范数平方为1。Wick 收缩与 Vandermonde 给精确式

`E w_k = v C(N,p) C(p,k) C(N−p+k,k)`。

设缺额 ℓ=p−k，其相对最高次数权重为 (p)_ℓ²/[ℓ!(N)_ℓ]。相邻比≤p²/[(N−p)(ℓ+1)]，控制归一化与每个固定矩。在 p²/N→λ>0 的 DS 极限，归一化缺额趋 Poisson(λ)，E||H||²→Js²/λ。因而 mixed、扣标量后的 mixed、anomalous 三种张量范数除N均趋4Js²。

这是先 t→0 取曲率，再 DS 的结论。它不提供固定有限时间的一致误差；也没有求完 singlet 双局域对易子四项之间的相消。下一步必须处理该投影，不能把未投影张量的正范数当成投影后不闭合的证明。

## 已验证和下一步

micro.json：5组精确 Hamiltonian/CAR、77组整数恒等式、54个单项协方差及15组80位有限参数记录。projection.json：四个内部源，解析完整矩与独立积分。参数/依赖/输入SHA256都已保存；本轮两脚本无失败，不重跑覆盖。

下一任务是 S(t₁,t₂)=Σc†_i(t₁)c_i(t₂) 的完整投影代数与协方差，寻找可供 boson expansion 的极化/态与非线性 dressing。可对照 Cavicchi–Di Vecchia–Pesando 的体 collective action，逐项匹配二次内部核及相互作用。Jevicki–Suzuki–Yoon 的红外 SYK 双局域工具可参考，但不直接移植其 conformal/AdS saddle 到当前无限温 complex DS 问题。内部夸克分数、单介子 boost 与两介子相对分数保持区分。
