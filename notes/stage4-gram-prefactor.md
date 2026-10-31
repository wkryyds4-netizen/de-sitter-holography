# Gram 因子与横向能量涨落

日期：2026-10-30。作者自查推导草稿；未独立认证或确认首创性。
英文推导和证明见[LaTeX](../manuscript/gram-prefactor/main.tex)，阅读版见[五页 PDF](../output/pdf/de-sitter-gram-prefactor.pdf)。前一步为[WKB重叠](stage4-wkb-overlap.md)，原文位置见[来源记录](../literature/stage4-gram-prefactor-sources.md)。

## 对象与结论边界

这里比较三维 CLP 固定拓扑四面体图与 CSTY 固定平均观察者质量的混合表示反极点关联函数。取 de Sitter 半径为1、观察者 proper 时间 t。该对象不是二维有限真空包，也没有在本轮建立新的物理算符代数或态字典。CSTY 使用指定经典纯反冲鞍点；使用 CLP 的 Gram 一圈渐近则额外依赖解析延拓及退化一致性假设。

本轮精确结果是 Gram 行列式的因式分解、外部能量变换的 Jacobian 与相位、以及等探针质量时经典作用量的横向 Hessian。将它们解释为完整量子振幅仍有条件。两项局部涨落组合给出此前需要的能量依赖，但未证明它们穷尽物理振幅。

## 1. Gram 行列式及极限次序

按 CLP 顶点三元组123、156、246、345排列面法向 Gram 矩阵，记 D=det G。令

\[
\theta_{1,4}=H+d,\quad\theta_{3,6}=H-d,\quad
\theta_2=\pi-\beta',\quad\theta_5=\pi-\beta,
\quad d=\frac{\pi\omega}{2\kappa},\quad
(\beta,\beta')=\frac\pi\kappa(\nu,\nu').
\]

交换矩阵的1、3及2、4指标，分解到两个二维子空间，精确得到

\[
D(\omega)=16[\cos^2(\beta/2)\cos^2(\beta'/2)-\cos^2H\cos^2d]
[\sin^2(\beta/2)\sin^2(\beta'/2)-\sin^2H\sin^2d].
\]

固定 H 远离0、π，ν/κ、ν′/κ→0，A=πνν′/(2κ sin H)，有 D(Az)/D(0)→1−z²，对有界复 z 一致；相对使用须远离±1。这个结论只关于精确行列式。

另一条路线先令 H=π(1−M/κ)，固定 M、ν、ν′、ω，κ→∞：

\[
D(\omega)=\frac{\pi^6}{4\kappa^6}
(4M^2-\nu^2-\nu'^2+\omega^2)(\nu^2\nu'^2-4M^2\omega^2)
+O(\kappa^{-8}).
\]

记 B=4M²−ν²−ν′²>0、A₀=νν′/(2M)，归一行列式趋于 (1+ω²/B)(1−ω²/A₀²)。还需 ν/M、ν′/M→0、ω/A₀ 有界才能舍弃第一个因子。不能将固定亏角公式直接代到 H=π。

条件性双手征 Gram 因子为 R_G=√(D(0)/D(ω))，先定义在实 Euclidean 分支，再沿解析路径延拓；不能在复 ω 上重新取模平方。纯反冲再重观察者极限下，ω=iA₀y 给 R_G→(1+y²)^(−1/2)。单独它只补足所需因子的一半。

## 2. 时间关联函数需要边缘积分

从 CSTY (2.3) 出发，设

\[
m_1=\bar M+\omega/2+u,\quad m_3=\bar M+\omega/2-u,
\qquad m_2=\bar M-\omega/2+v,\quad m_4=\bar M-\omega/2-v.
\]

绝对 Jacobian 为4。反极点四段时间间隔为 π/2−it、π/2+it、π/2−it、π/2+it，相位的线性组合为2πM̄−2itω。平均质量 delta 积分后，原定义中的热因子抵消固定部分，留下

\[
\mathcal F(t|M)=\frac4{Z_{\rm sphere}(M)}
\int d\omega\,e^{-2it\omega}\int_{\mathcal D_\omega}du\,dv\,
\boldsymbol{\mathcal F}(m_1,m_2,m_3,m_4).
\]

积分域和轮廓继承原定义。u=v=0的截面不等于这个边缘积分。横向高斯给2π/√det I″⊥，它与内部6j积分的涨落不是同一积分，不能混为一次计算或重复计数。CSTY脚注19已有正交能量积分的说明，因此这里不是发现原文遗漏，而是恢复先前比较中未计算的步骤。

## 3. 等探针质量的经典 Hessian

令 ν′=ν>0，B=4M²−2ν²>0，A₀=ν²/(2M)，先在 |ω|<A₀ 的实非退化分支。采用 CSTY (2.54)、(B.5)—(B.7) 的六边作用量 I=Σmᵢⱼdᵢⱼ，顶点顺序与英文稿一致。利用在壳微分 dI=Σdᵢⱼdmᵢⱼ，直接对原余弦公式求导可得

\[
I''_\perp=h(\omega)\mathbf1_2,\qquad
h(\omega)=\frac{4\sqrt{(\nu^4-4M^2\omega^2)(4M^2-2\nu^2+\omega^2)}}
{(4M^2-\nu^2)(\nu^2-\omega^2)}.
\]

平方根在零点取正。延拓只到根、极点或分支障碍出现以前；局部正性不证明 Lorentzian 轮廓必经该鞍点。归一横向因子为

\[
R_\perp=\frac{h(0)}{h(\omega)}
=\frac{1-\omega^2/\nu^2}{\sqrt{1-\omega^2/A_0^2}\sqrt{1+\omega^2/B}}.
\]

与纯反冲 Gram 比相乘：

\[
R_GR_\perp=\frac{1-\omega^2/\nu^2}
{(1-\omega^2/A_0^2)(1+\omega^2/B)}
\longrightarrow\frac1{1+y^2},\quad\omega=iA_0y,\quad\nu/M\to0.
\]

对有界实 y 一致；在此前外鞍点 y=sinh t 给 sech²t，正好满足上一轮的必要振幅条件。这是从几何和边缘积分得到的条件性解释，不是把目标函数拟合回测度。

h(0)~2/M，横向宽度约√M。要让该宽度小于探针尺度，需要√M≪ν≪M。结合前一轮有界时间下 A₀ cosh t≪ν^(2/3)，存在兼容层级 ν^(4/3)≪M≪ν²。这只是尺度兼容性检查，尚未给完整积分的一致高斯余项。

## 4. 仍缺什么

CLP (3.57)只处理领先指数，不能推出顶角 C 的一圈项抵消。精确顶角函数已给出，其联合大参数渐近尚待计算。CLP (3.37)的物理两点插入归一化只给固定探针数据下的领先项；不能单凭它确定增长探针的联合极限。英文稿给出的乘性例子仅说明这种渐近信息不足，不主张任意补全符合全部物理一致性。

此外须控制退化球面上的6j一圈余项、两个横向方向、外部三能量轮廓、其他鞍点及相位。Meng–Yang固定双曲角定理不能直接提供上述全部控制。未得到完整有限 G 的 OTOC 或微观字典。

## 5. 数值与保存

[脚本](../code/gram_prefactor.py)的[83项检查](../results/stage4/gram-prefactor/diagnostics.json)全部通过：3项符号检查、5个一般判别式、9组固定角、6组纯反冲、6组仅带Gram因子的逆积分，以及9组横向Hessian。Hessian通过直接对完整六角作用量求导与解析式比较，含实与纯虚能量样本。工作65位；指定比值与积分用95/90位复核，未声称每个Hessian都做了加精度复核。

诊断保留相对误差、参数、依赖版本与4个输入哈希；它们不是区间认证。仅带Gram的逆积分仍保留趋向cosh t的剩余差异，未隐去未匹配结果。新[证据记录](../results/stage4/gram-prefactor/evidence.json)绑定源码、来源、中文记录、PDF与最终日志。二维55页稿及前一轮4页稿保持冻结。
