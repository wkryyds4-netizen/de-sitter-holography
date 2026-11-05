# 从complex DSSYK的现有解到首个singlet谱

2026-11-05。沿Susskind主线继续，上一轮规划已成为本轮明确的算符与谱密度问题。[六页LaTeX](../manuscript/susskind-dictionary/main.tex)与[PDF](../output/pdf/de-sitter-susskind-dictionary.pdf)保存完整公式；原文逐项位置见[来源表](../literature/susskind-dictionary-sources.md)。

## 已有解可以直接用在哪里

BNR（2006.13983v2）已经给complex双标度SYK的Hamiltonian chord矩和通用随机探针关联。我们将其binom(N,p)⁻²方差与MSS的(p!/Nᵖ)²方差精确匹配，得到J_B=(J_s/√λ)∏(1−j/N)。乘积在固定λ下趋于e^(−λ/2)，与BNR矩的尺度因子相消。Hamiltonian二阶矩在双标度极限为J_s²/λ，弦交叉权重q=e^(−4λ)。这是已有解的规范转写。

代码采用{c,c†}=1，原文ψ=√2c，因此2p次Hamiltonian系数要乘2ᵖ，方差乘4ᵖ。时间单独按t_s=p t_c、H_c=p H_s换算；导数复合算符的幅度也要一起换。固定电荷二阶矩直接复用Tall–Tomsovic（2604.14522v2）及其所引embedded-ensemble恒等式。

## 真正要构建的对象

令T(X)=Σcᵢ†Xcᵢ，标准CAR给

**A₂=T(H²)−(Q−2p)H²。**

在m粒子扇区，T(H²)连接到同一样本的H_(m−1)²；后一项守恒。因此A₂的非零频谱只需T(H²)的非对角矩阵元。展开中间态可把它写成四个c/c†矩阵元及两个中间能量平方权重的相干组合。稿中Ξ_m正是这一带电四点密度，给下一轮大N计算明确的外部插入和能量权重。

特别地，m=p扇区有A₂=pH²；p=1时整个A₂都守恒。这两个控制比只检查正性更直接。与之相对，N=4–6的交互样本有非零动态谱权重，数据保存为原矩阵及所有离散跃迁，没有将其外推为双标度介子质量。

## 本轮计算与接续

5组完整协方差矩阵求和与精确组合公式一致；另一套粒子数扇区公式独立复核一致。18组固定seed样本完成178项代数及谱矩检查，最大scaled residual约1.30e−15。另由相邻扇区能量本征态的相干和核查四点重建，输出单列，不改写原校验记录。

下一步先读Tall–Tomsovic §3、§4.2与附录A，再对照BNR §3及附录C，算Ξ_m的低阶能量加权矩并尝试重求和。现有闭式探针解采用独立随机系数；A₂使用同一J，需保留这类额外收缩。暂不扩大同类小N样本，也不回旧五边形支线。

复算入口：`python code/run_complex_syk_dictionary.py`只允许新输出，`python code/audit_complex_syk_dictionary.py`与`python code/check_singlet_charge_bridge.py`首次写审计。已有输出后可导入`audit(root)`、`verify(root)`只读重放。环境Python3.14.5、NumPy2.4.6，详细参数及哈希在results/susskind-dictionary/calibration.json。
