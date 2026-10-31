# 双共形微观四点：时间修正与一维积分核

2026-10-29。接续 [条件性系数比较](stage4-microscopic-match.md)。本轮得到两项可复算结果：保留 C.39 的线性时间和 digamma 后，各阶出现不同的修正因子，原来的 m=2s₀ 对应只在 leading-log 层成立；同一反极点配置下，C.6 的三个频率积分可降为一个无实轴奇点的积分，因而能够直接检查 C.39 额外丢弃的有限时间项。

输入是 [Narovlansky 2506.02109v2](https://arxiv.org/html/2506.02109v2) 的 C.6、C.9、C.13、C.18、C.21、C.39。下文是对这些既有公式的特例整理与验证，未作首创性判断。

## 1 保持原来的模型、归一化和极限

采用原文固定谱态 |s₀⟩₁|s₀⟩₂、等能投影、双拷贝复合算符。两个探针均取 Δ₁=Δ₂=1/4，另一拷贝参数为 3/4，对应 dS₃ 共形标量。

原文先取 λ→0，再取高于基态的无量纲能量大；C.6 已经使用这组 holographic-semiclassical 近似。以下“一维完整核”仅指**C.6 这条已近似公式的全部四个 sector**，不是有限 λ 的完整 DSSYK 模型，不是有限 G 引力结果，也不是有限观察者质量字典。

令 uᵢ=2πT tᵢ、t=2πT t_phys，保持

\[
u_1=i\pi/2,\quad u_3=-i\pi/2,\quad
u_2=t,\quad u_4=t-i\pi.
\tag{T1}
\]

这里公式中的 u₁、u₂ 是时间变量，非物质谱标签 ν。两个同种插入的差均为 −iπ。可用

\[
u_1=i(\pi-\epsilon)/2,\quad u_3=-i(\pi-\epsilon)/2,
\quad u_2=t,\quad u_4=t-i(\pi-\epsilon),\quad \epsilon\downarrow0
\tag{T2}
\]

从原负虚时间条带接近。所有结果保持同一个投影尺度 𝔍，并仍定义

\[
R=\frac{\langle O_2O_1O_2O_1\rangle}
{\langle O_2O_2\rangle\langle O_1O_1\rangle},\qquad
C_0=\frac{\mathfrak J\log s_0}{4\pi^2T},\qquad \mathcal R=R/C_0.
\tag{T3}
\]

没有重新选择 𝔍 把某个时间点归一成 1；下面 n=0 项不等于 1，是保留原 C₀ 后的实际结果。

## 2 C.39 中 digamma 积分可以精确完成

设 a=1/4、b=3/4。反极点使 C.39 的两个独立频率权重指数归零。所需 Barnes 积分为

\[
I(z)=\int_{\mathbb R}dw\,
\Gamma(a+z+iw)\Gamma(b+z+iw)\Gamma(a-iw)\Gamma(b-iw)
=\frac{2\pi\Gamma(z+1/2)\Gamma(z+1)^2\Gamma(z+3/2)}{\Gamma(2z+2)}
=\pi^{3/2}4^{-z}\Gamma(z+1/2)\Gamma(z+1).
\tag{T4}
\]

Barnes 第一引理的实部条件在 z=n≥0 和 z=±iω、ω 实时都满足。I(0)=π²，I(n)/I(0)=n!(1/2)ₙ/4ⁿ。另一个频率因子是 w→−w 后的同一积分。

C.39 每个频率上的两项 digamma 正好是对 z 求导，因此积分后的平均为

\[
D_n=\partial_n\log I(n)
=\psi(n+1/2)+\psi(n+1)-2\log2.
\tag{T5}
\]

可在固定 n 的邻域交换求导与积分：四个 Gamma 给出指数衰减，digamma 只添加对数增长。这里不需要把 n 当作离散差分，也没有忽略任一 digamma 项。

实际时间差 u₄₃=t−iπ/2，使原括号的 −π 与欧氏时间部分恰好抵消；外部指数的 (−1)ⁿ 也被相位抵消。由此，原 C.39 保留的早时渐近块为

\[
\boxed{\quad
\mathcal R_{C.39}\sim\sum_{n\ge0}(1/2)_n^2
\left(\frac{e^{2t}}{64s_0^2}\right)^n
\frac{B_n(t,s_0)}{\log s_0},\qquad
B_n=\log(8s_0)-t-\psi(n+1/2).
\quad}
\tag{T6}
\]

令 B=B₀=log(32s₀)+γ−t，则

\[
B_n=B-2\sum_{k=0}^{n-1}\frac1{2k+1},\qquad
B_1=B-2,\quad B_2=B-\frac83.
\tag{T7}
\]

具体到前两阶：

\[
\mathcal R_{C.39}\sim\frac{B}{\log s_0}
+\frac{e^{2t}}{256s_0^2}\frac{B-2}{\log s_0}
+\frac{9e^{4t}}{65536s_0^4}\frac{B-8/3}{\log s_0}+\cdots.
\tag{T8}
\]

这些是指定渐近块的系数，不是对 C.6 全积分剩余误差的完整展开。对旧 leading-log 结论，安全层级仍是固定 n、1≪t≪log s₀；此时 Bₙ/log s₀→1，原系数身份保持正确。

## 3 一个尺度重标度不能吸收这项改动

保持原 C₀，零阶已经出现 (γ+5log2−t)/log s₀。n=1、2 又有不同常数，不能把整个变化叫作同一个归一化因子。

即使仅作为诊断，另外除以本块的 n=0 值，并用首阶重新选择一个有效尺度，次阶仍不匹配。令 x=e²ᵗ/(64s₀²)，则诊断性的级数为

\[
1+\frac14x\frac{B-2}{B}+\frac9{16}x^2\frac{B-8/3}{B}+\cdots.
\]

单一 x 重标度不能改变的比值，本来是 9，现在成为

\[
\frac{[x^2]\mathcal R_{C.39}/\mathcal R_{C.39}^{(0)}}
{([x]\mathcal R_{C.39}/\mathcal R_{C.39}^{(0)})^2}
=9\frac{B(B-8/3)}{(B-2)^2}
=9\left[1+\frac4{3B}+O(B^{-2})\right].
\tag{T9}
\]

这里 [xⁿ] 指在固定 t、s₀ 的形式组织中提取幂项；没有把 B(t) 当成与时间无关的物理常数。B=3 的孤立代数例外不构成整个早时区间上的重标度。

因此，m=2s₀ 的全阶身份是 leading-log 身份；保留 C.39 时间/digamma 项后，不再是同一单 sector 纯反冲函数的全阶身份。这个判断没有给 s₀ 赋予观察者质量。

## 4 C.6 的四项合并后降为一维

只读 C.39 还不足以称为“完整次 leading-log”：原文对另外三个 sector 也做了大时间/大 s₀ 的局域化。因此进一步从 C.6 本身保留全部四项。

取

\[
\omega=\omega_1+\omega_3-\omega_2,\quad
\Omega_1=\omega_1,\quad\Omega_2=\omega_2-\omega_1,
\]

Jacobian 为 1。两个非混合项的内部频率积分给 I(∓iω)²；混合项则给

\[
M(\omega)^2,\qquad
M(\omega)=2\pi\Gamma(1/2-i\omega)\Gamma(3/2+i\omega)
=\frac{2\pi^2(1/2+i\omega)}{\cosh(\pi\omega)},
\tag{T10}
\]

以及 M(−ω)²。把 L₂=log(2s₀)，组合核为

\[
\begin{aligned}
K(\omega)={}&e^{2iL_2\omega}\Gamma(i\omega)^2I(-i\omega)^2
+e^{-2iL_2\omega}\Gamma(-i\omega)^2I(i\omega)^2\\
&+\Gamma(i\omega)\Gamma(-i\omega)
\left[M(\omega)^2+M(-\omega)^2\right].
\end{aligned}
\tag{T11}
\]

这一步先对非零 ω 做 Barnes 积分，保持四项的共同调节，再相加取消调节。不能分别把各个 1/ω² 项当成普通实积分。两项非混合核各有 −π⁴/ω²，混合项合计 +2π⁴/ω²；合并后双极点消失，奇次项也消失。共同的小轮廓位移可在不跨越其他 Gamma 极点的范围内撤去。更直接地，在 C.6 原 integrand 的每一个 Δ 括号中，ω=0 时两项的其余分子相同，而 Γ(iω) 与 Γ(−iω) 的简单极点留数相反；因此在内部频率积分之前已点态相消。这说明共同调节的撤除没有另外添加接触项。

令

\[
\ell=\log(8s_0),\quad
\Phi(\omega)=2\omega\ell-2\operatorname{Im}\log\Gamma(1/2+i\omega),\quad
A(\omega)=\frac{(1-4\omega^2)\tanh(\pi\omega)}{\pi\omega}.
\tag{T12}
\]

由 Gamma 反射公式，K=2π⁶(A−cosΦ)/(sinh²πω coshπω)。总前因子可以直接从原 C.6 与二点式 A.3 核对，避免仅靠 leading-log 首项反推。令 D=Γ(1/2)Γ(3/2)=π/2；每个独立谱积分恢复一个相同的 𝔍，则原 C.6 前因子及反极点二点分别为

\[
P_4=\frac{\mathfrak J^3(2\lambda s_0)^4}{(16\pi^3T)^3D^2},\qquad
Q_2=\frac{\mathfrak J(2\lambda s_0)^2}{16\pi^3T D}I(0),\qquad I(0)=\pi^2.
\tag{T12a}
\]

这里使用与 A.3/C.6 相同的算符约定；若同时使用正文复合算符的额外常数，其四个算符因子在四点/两个二点之比中抵消。两个同种探针的二点都为 Q₂，所以

\[
R=\frac{P_4}{Q_2^2}\int_{\mathbb R}d\omega\,e^{-2it\omega}K(\omega)
=\frac{\mathfrak J}{16\pi^3T I(0)^2}
\int_{\mathbb R}d\omega\,e^{-2it\omega}K(\omega),
\]
\[
\mathcal R=\frac{R}{C_0}
=\frac{1}{4\pi\log s_0\,I(0)^2}
\int_{\mathbb R}d\omega\,e^{-2it\omega}K(\omega).
\tag{T12b}
\]

K 是实偶函数。把全实轴积分改为两倍正半轴余弦积分，再代入 K 的 2π⁶ 和 I(0)²=π⁴，净系数为 2·2π⁶/(4π·π⁴log s₀)=π/log s₀。这在固定同一个 𝔍 下明确给出

\[
\boxed{\quad
\mathcal R_{C.6}(s_0,t)=\frac{\pi}{\log s_0}
\int_0^\infty d\omega\,\cos(2t\omega)
\frac{A(\omega)-\cos\Phi(\omega)}
{\sinh^2(\pi\omega)\cosh(\pi\omega)}.
\quad}
\tag{T13}
\]

其原点的可去极限是

\[
\lim_{\omega\to0}\frac{A-\cos\Phi}{\sinh^2(\pi\omega)\cosh(\pi\omega)}
=\frac{2[\log(32s_0)+\gamma]^2-4-\pi^2/3}{\pi^2}.
\tag{T14}
\]

无穷端至多为多项式乘 exp(−3πω)，所以这是一条普通绝对收敛的积分。实性和 t→−t 对称来自原四项的代数合并；没有对结果取实部，没有平均两个候选引力轮廓。该特定欧氏插入配置也不能据此确定一般微观/引力的轮廓字典。

式 (T13) 可以计算 C.6 模型在任意有限 t 的函数，不要求先截断 C.39。它不恢复 C.1→C.6 已丢失的有限 λ 或大 s₀ 非增强修正；对这些修正没有误差认证，也没有宣称与引力晚时曲线匹配。

## 5 实际检查到的有限时间余项

下表所有值保持原 C₀。最后一列为全 C.6 核减 C.39 的 n=0、1、2 项，不是有限 G 效应。

| log s₀ | t | C.6 一维积分 | C.39 的 n=0 项 | 差值（扣除 n≤2） |
|---:|---:|---:|---:|---:|
| 20 | 5 | 0.94784702623355 | 0.95214757838506 | −0.00430055215152 |
| 20 | 10 | 0.70208977288891 | 0.70214757838506 | −0.00005780550100 |
| 50 | 6 | 0.96010097348066 | 0.96085903135403 | −0.00075805787336 |

第一行 C.39 的 n=1 项只有 3.1149×10⁻¹⁶；第三行只有 2.1779×10⁻⁴¹。有限时间余项远大于它们。因此不能把“已保留 C.39 的 digamma”误写成“完整微观函数已受控到该指数项的精度”。这些差值包含 C.39 对其他 sector 的额外渐近化及本次有限阶截断；没有逐项宣称其全部来自某一种物理机制。t=10、log s₀=20 只是积分核例子，并不满足旧 leading-log 的 t≪log s₀ 层级。

这同时提供了实际后续工具：现在可以用 (T13) 分离有限时间效应、归一化问题与渐近系数匹配，不必用一条 leading-log 形式级数代替完整 C.6 积分。

## 6 独立验证与复算

[代码](../code/microscopic_time_corrections.py) 与 [结果 JSON](../results/stage4/time-corrections.json) 为本轮新增文件，既有第四阶段代码和输出未修改。

50 位工作精度下 **41 项检查通过，0 失败**：

- n=0、1、2、4 的直接实轴 Gamma 积分和 digamma 加权积分，对照 Barnes 闭式与其连续参数导数；
- n=0…8 的 C.39 原括号、欧氏相位和奇调和数整理；
- 两个不同复移位的混合 Barnes 积分；
- 四 sector Gamma 形式与无奇点三角函数形式的逐点对照；
- ε=0.02、0.002、0.0002 从同一负虚时间条带接近的一阶系数；
- 三个一维全核算例的 40→50 位精度及分段加密检查，最大差约 6.37×10⁻⁴²。

最后的差值是算法一致性，不是原物理近似拥有 42 位精度的声明。JSON 记录源缓存 SHA-256、代码 SHA-256、参数与环境。复算命令：

```powershell
python code/microscopic_time_corrections.py
```

本轮具体解决了“时间/digamma 是否只是整体常数”的问题，答案是否定的；还得到一个可直接计算所有 C.6 sector 的特例核。仍待做的是该核的独立解析交叉审查、与相同微观态/可观测量的引力轮廓比较，以及新推导的查新。有限观察者质量、有限 G 和有限钟制备均未因此接通。
