# 固定多项式制备的 Slater 轨道密度与精确 pinching

2026-11-08。本轮从已封存的线性 orbit-density 结果继续，先核矩阵 Gaussian chaos 前驱，再推固定次数制备。正式逐步证明在 manuscript/susskind-slater-polynomial-density/main.tex。不是将线性 Gaussian 矩递推直接用于多项式。

## 模型和命题

沿用 standard CAR、overlapping complex Hermitian Gaussian p-body、同一原 J、半填充 m=N/2、annealed Slater O。实际观测量 h,u,M,M† 的归一及 string time 不变。一般 P_O(J) 的 Gaussian degree≤固定 K，且对 O 的稳定子 U(m)×U(N−m) 协变。按 Wick chaos 分解 P_O O=Σ_k P_k O，a_k=E||P_k O||²、a=Σa_k。系数允许随 N 变化。

设 λ_N=p²/N≤Λ，p≤N/(12K)，C_K=exp(36K²Λ+1)，B_0r=1，B_kr=[k(2r−1)]^k C_K^(1−1/r)/k!。存在正的协变 R_P，使 covariant Y∈joint L² 满足

E〈O|P†YP|O〉=E tr_m(R_PY)，
||R_P||_r≤(Σ_k sqrt(B_kr a_k))²≤(K+1)max_k B_kr a。

这里 r≥1 为整数，tr_m 为归一扇区迹。常数很粗，但每个固定 K,r,Λ 下不依赖 N。DS 固定 λ>0 时条件最终成立；不是随次数增长的一致估计。旧线性界更强且适用于非半填充，保留不变。

## 系数张量与部分权重

q_U=U_m†P_O(U·J)O，Zf=√d_m∫q_U f(U)dU，R_P=ZZ†。只需有限维 matrix coefficient 的 Haar 支撑。k-chaos 用有序对称张量 A_i 写 Σ A_i :g_i1...g_ik:，Ψ∈W^⊗k⊗H_m、a_k=k!||Ψ||²。轨道系数 Ψ_U=ρ(U^-1)Ψ。稳定子 torus 给整体权重 w=(1^m,0...)。

复化实 Hermitian noise 并换正交 weight basis；其表示与 conjugate Ad(wedge^p) 酉等价，每条 noise weight 是两个 p-subset 指示向量之差，总和0、l1≤2p。概率估计仍用实独立 Gaussian；酉基变换只用于不改变奇异值的系数 flattenings。

trace 去掉 ℓ 条 noise legs 后，T_I=Σβ |ψβ〉〈ψβ|，每个 ψβ 的权重 z=w−νβ 满足 sum z=m、||z−w||1≤2ℓp。不能说部分种子仍是精确 Slater 权重。

## near-Slater 维数下界

令 b=#(z_i不属{0,1})、S 为这些坏值之和。z−w 正负总质量各≤ℓp，所以 b≤2ℓp，−ℓp≤S≤b+ℓp，|S|、|S−b|≤3ℓp=D。零/一的计数 n0=m−b+S、n1=m−S；δ0=S−b、δ1=−S，和为−b，绝对值≤D≤m/2。

正 δ 用 log(1+x)≤x，负 δ 用 −log(1−x)≤2x，得到 log[(m+δ)!/m!]≤δ logm+(δ²+|δ|)/m。坏值各 multiplicity 的阶乘积≤b!，故 Weyl orbit / d_m≥m^b/b! exp[−(2D²+2D)/m]≥exp[−36ℓ²Λ−1]≥C_K^-1。每个支持 irrep 含整个置换权重轨道，所以 dim Vα≥d_m/C_K。

## 所有 flattenings 与概率界

F_I 将物理输出及 noise I 放行、Haar输入及 noise J 放列。列侧 noise 的酉矩阵求和相消，得到 F_I F_I†=d_m∫ρ_I(U^-1)T_Iρ_I(U^-1)†dU，HS²=d_m||Ψ||²。将每个 ψβ 投影到 isotypic 块仍保留其权重，故所有支持块适用上节下界。Schur twirl 保留 multiplicity matrix，其 norm≤trace≤||Ψ||²，给 ||F_I||op²≤C_K||Ψ||²；进而 ||F_I||_(2r)^(2r)≤d_m C_K^(r−1)||Ψ||^(2r)。

线性矩阵 Gaussian X=Σg_i A_i 经 Hermitian dilation、IBP 和 Schatten Hölder 给 E||X||_(2r)^(2r)≤(2r−1)^r/2 [Tr(ΣAA†)^r+Tr(ΣA†A)^r]。对 decoupled k-chaos 逐个条件 Gaussian 应用，每步对应 row concatenation / column stacking，最终恰好 2^k 个 orientations，系数 (2r−1)^(kr)2^-k。结合 F_I 界得归一 E||Z_dec||_(2r)^(2r)/d_m≤(2r−1)^(kr)C_K^(r−1)(a_k/k!)^r。

令 g=Σ_j g^(j)/√k。条件均值 g/√k、不同 slot covariance=−I/k。生成函数逐slot微分给 E[Π_j g^(j)_ij|g]=k^(-k/2):Π_j g_ij:，重复坐标也成立。故 Z_k=k^(k/2)E[Z_dec|g]，Schatten Jensen 给所述 B_kr。k=0 时 P_0 O=cO、R_0=|c|²I。最后对 Z=ΣZ_k 用 Minkowski，不能逐样本删去跨 chaos 干涉项。

## 实际精确能量投影与二次响应

对 covariant X，E||XP O||²≤||R_P||2 ||X||4²。固定 selected 多项式 P(h,u,M,M†) 的 a_N 由既有同J全word定理有界。旧 t_N−u_N 与 P_H t_N 在所有固定扇区 S_r 趋零，因此实际 b_N=(1−P_H)t_N 满足 E||(b_N−u_N)P_NO||²→0，覆盖整个固定多项式制备 core。

对固定 covariant block W(h,u,b) 及内部固定 h-time phases，trace Hölder/telescoping先给 W(b)−W(u) 的 L² 小量，再用 R_P 和极化得到 φ_N(P_N† W_N(b_N) Q_N)→φ(P†W(h,u)Q)。P,Q 可以各有任意固定数量 selected 因子，但本定理不直接覆盖多个非多项式 pinched blocks 之间插 M 的所有 raw moments。不能由 core 收敛偷换这一结论。

P=M†h 的内部收缩配反向 hop，右侧 annihilates E、杀死 O；P=(M†)² 的两个 selected 标签不能互为逆。因此两个向量都是 pure chaos2。transfer 给 a1=q(1+√q)、a2=2q。归一实际斜率二点分别为

P=M†h：1+λ√q q(1+q)/(1+√q)；
P=(M†)²：1+2λ√q q(1+q)。

这两个式子先在有限N exact pinching 后再 N→∞。物理参数 q=e^-4λ，固定 λ>0，不称热平衡态或体质量；小有限N 的双激发可能为零，该处不归一。

## 证据与下一步

五组实际半填充 CAR fixtures：1610 full seed weights、6440 partial weights；20199 rational factorial steps、342 near-weight计数；126 degree≤6条件Wick恒等式；8个小型矩形tensor的全orientation/decoupling矩比较。独立 raw pairing 和 symmetrized tensor 核 k!；四个极限 transfer 多项式另核。参数、版本、代码与父manifest hash均在 exact-checks.json。首次运行成功后封存，未重跑。有限检查不代替一般证明，不拟合DS。

下一步优先实际二次协变 J2 与 M 的同J混合结构，使用本轮密度控制处理其 tracial approximation / pinching 接口；先复用已完成 fulltrace J2 结果，不能把 fulltrace 四点再当空缺。更高 nonpolynomial interlacing 的矩尾仍单列。随后 fullbilocal、空间/质量生成元与 finite-G。仅本地保存，无新 Git/issue。
