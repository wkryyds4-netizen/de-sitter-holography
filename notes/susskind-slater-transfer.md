# 同 J 开弦传递、Hamiltonian阶梯与选择激发谱

2026-11-07。[完整LaTeX证明](../manuscript/susskind-slater-transfer/main.tex)、[PDF](../output/pdf/de-sitter-susskind-slater-transfer.pdf)、[原文核读](../literature/susskind-slater-transfer-sources.md)。沿用standard CAR、半填充annealed Slater态、same original Gaussian J与overlapping p-body H，q=e^(−4λ)、γ=q^(1/4)、s=√q；固定λ>0及固定时间，h=H/√(wD)无量纲、ℏ=1/string H单位。此页为joint-words之后的新stage，旧科学数据不覆盖。

## 1. 全transfer与正metric

按打开时间由老到新记录w=t₁…tₙ，t∈{H,P}。a† append H，b† append P，b=M。删除第j个endpoint得到w\j：

\[
a|w\rangle=\sum_{t_j=H}q^{n-j}|w\setminus j\rangle
+\gamma q^{n-1}\sum_{t_j=P}|w\setminus j\rangle,
\]
\[
b|w\rangle=q^{n-1}\left(\gamma\sum_{t_j=H}|w\setminus j\rangle
+\sum_{t_j=P}|w\setminus j\rangle\right),\quad h=a+a^\dagger.
\]

首次关闭edge与newer open形成cross，与older open形成inner nesting；仅selected closing对后一项加q，因此得到exact weights，不需在打开时猜最终edge是G还是S。逐word空串系数等于旧allword定理。

String不能当正交单位基。不同length正交，等length的Gram用所有跨cut matching，或按bra的反序letters逐a/b关闭ket。图反射给对称性，移动一个endpoint给creation/annihilation伴随。递归χ_wP=M†χ_w、χ_wH=hχ_w−Σ a_vwχ_v，把formal strings送入实际GNS，Gram=φ(χ_u†χ_v)≥0。商null后各level维度≤2^n；未证明每level正好2^n。固定r个M†插入的h传播每level至多Σ_(j≤r)C(n,j)个状态，可稀疏计算。

## 2. 全graded代数及Type I

\[
aa^\dagger-qa^\dagger a=I,\quad
ab^\dagger-qb^\dagger a=ba^\dagger-qa^\dagger b=\gamma q^{\mathcal N},
\]
\[
bb^\dagger-qb^\dagger b=q^{\mathcal N}.
\]

每个creator升总grade一。逐level正性递推给||a†_n||²=[n+1]_q、||b†_n||²=(n+1)q^n，纯H/P链达到界。因此a,b有界，b=M在full joint空间compact（finite level + level norm→0）。Q=MM†−qM†M=q^N正，Q^r范数收敛到rank-one vacuum projection。该projection与joint cyclicity生成所有rank-one，W*(h,M)=B(H)，指定向量态纯。这是RX state-adapted Type I机制的本模型实现。

β=q^(−N/2)b、β†=b†q^(−N/2)在finite-level core满足[β,β†]=I，aβ†−√qβ†a=γq^(N/2)。逆dressing无界；尚未把fulljoint N无推导认同微观Q_E/p。

## 3. h的正交primitive towers

K_d=ker(a|H_d)，H_n=a†H_(n−1) ⊥ K_n。每个primitive ξ生出(a†)^r ξ/√[r]_q!，不同primitive towers正交，穷尽空间。因此H≅ℓ²(N₀)⊗K，h≅J_q⊗I_K；K=⊕K_d的具体维数无需假定。J_q offdiag√[r+1]_q，||h||=2/√(1−q)，标准q-Gaussian μ_q为每tower底态谱。

Π|H_n=Σ_(j=0)^n (−1)^j q^(j(j−1)/2)/( [j]_q! ) (a†)^j a^j投影到ker a；v=Σ_r(a†)^r Πa^r v/[r]_q!。h谱虽相同，M仍耦合不同primitive towers；不能删掉multiplicity后计算多插入。

## 4. 具体谱核与态制备

η=(b†−γa†)Ω为primitive，norm²=1−s，故b†Ω=γa†Ω+η。对任意bounded Borel f：

\[
\varphi(Mf(h)M^\dagger)=\int(1-s+sE^2)f(E)d\mu_q(E).
\]

设R(t)=∫e^(itE)dμ，C_M(t)=(1−s)R(t)−sR″(t)，φ(e^(ith)M†)=−iγR′(t)。Euclidean向量ψβ∝e^(−βh/2)M†Ω的Z_M=(1−s)Z+sZ″，meanE=−∂βlogZ_M；这只是明确向量制备，不是全代数KMS态。

vacuum h-cyclic tower的单插入压缩P₀M†P₀为weighted shift，d_n=γ(n+1)q^n/√[n+1]_q。其两能量HS核B(E,E′)=Σ d_n p_(n+1)(E)p_n(E′)。尾Σ_(n≥L)d_n²≤s r^L[(L+1)²/(1−r)+2(L+1)r/(1−r)²+r(1+r)/(1−r)³]，r=q²。压缩只捕获||P₀M†Ω||²=s，完整为1；不能直接迭代该核当full四点。

## 5. 计算和失败记录

- exact-checks.json：1458保存sixth-word比较、232独立8字、248q关系、124projection、62tower、680adjoint、170symmetry、54 residual/compressed/spectral恒等式。
- canonical-checks.json：62grading、62CCR、62same-J cross、18sharp norm。
- spectral-checks.json：q=.0625/.25/.64，Jacobi160/240；15time差≤1.34e−15、12Euclidean preparation、9HS tail行，另8193theta/256product谱密度归一与二矩检查。t=0的log尾字段保存为Python JSON的−Infinity扩展，表示解析尾为零；其余数据正常有限。解析quadrature截断界不包括eigensolver roundoff，不称区间证书。

初checker在1/qfac(0)处误引入float，HHP分解出现8.7e−19伪余项，未产scienceJSON；失败source/说明留diagnostics。改用Fraction(1)后精确比较全部通过，成功outputs不重跑。数值图selected-spectrum.png由scientific driver生成，哈希同存。

下一物理任务应优先检查stationary observer态在该B(H)、连续谱表示中的存在条件，再明确grading/微观occupation及原bilocal的共同state接口。谱和全transfer已经足够开展这些测试，不再重复固定矩推导。finite-G、空间/质量与有限视界熵仍待建立。本轮仅本地，未Git或issue回复。
