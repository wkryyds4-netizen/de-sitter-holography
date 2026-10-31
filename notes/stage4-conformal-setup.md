# 第四阶段：双共形探针的有限 G eikonal 基线

2026-10-29。状态：公式与轮廓已从原文推导到可算的一维角积分；数值结果在本文末追加。这里不引入有限内钟态，也不声称得到完整有限 G 量子引力四点函数。

来源：[2607.14042v1](https://arxiv.org/html/2607.14042v1) 式 (2.8)、(2.12)、(2.25)—(2.26)、(2.33)—(2.35)、(2.39)—(2.42)、附录 G 的 (G.1)—(G.2)，并参照 §3 的字典检验。两探针均取 ν=ν′=i/2，即 M_field²L²=3/4。它与第三阶段重—轻公式 R=atanh(q)/q 是不同的参数极限。

## 对象、参数和归一化

仍使用固定平均观察者质量、带热权重的混合表示。L=1，t=τ_proper/L，m=LM_obs，G=G₃/L。设

\[
 \gamma=4Gm,\quad k=1-\gamma,\quad\alpha=2\pi k,\qquad
 0\leq\gamma<1,
\]

这里 γ 仅是本文件的无量纲质量亏角参数，不是第三阶段的轨迹折角。代码使用名字 `deficit` 避免混淆。

取附录 G 所用插入配置

\[
 (\theta_1,\theta_2,\theta_3,\theta_4)
 =(-\pi/2+it+\tau,0,\pi/2+it+\tau,\pi).
\tag{S1}
\]

τ 是欧氏偏移，不能与物理固有时间混用。主体结果在 τ→0 两种指定轮廓下计算。两探针对各相隔 π，因此断开乘积

\[
 \mathcal F_0={1\over16\alpha^2},\qquad R=\mathcal F/\mathcal F_0.
\tag{S2}
\]

这是在同一 eikonal 近似内的 TOC/断开基准；没有计算非增强的 TOC 引力修正。在 α=2π 时，F₀=1/(64π²)，与共形标量的两个反极点传播子乘积一致。

## 从 G.2 到单变量核

定义

\[
 A=e^{-i\theta_4}-e^{-i\theta_2}=-2,\quad
 B=e^{i\theta_1}-e^{i\theta_3}=-2ie^{-t+i\tau}.
\]

G.2 的前因子为 −ie^{-t+iτ}/(4α²)。把 p₋=−v、v>0，先按原文 F₁₂ 轮廓积分 p₊，得到

\[
 \mathcal F_{12}^{(2.32)}
 ={e^{-t+i\tau}\over4\alpha^2}
 \left\langle\int_0^\infty {e^{-2e^{-t+i\tau}v}\,dv\over2+h(\phi)v+i0}\right\rangle_\phi.
\tag{S3}
\]

令 x=2e^{-t+iτ}v 并保留连续轮廓，得到归一化核

\[
 J(c)=\int_0^\infty{e^{-x}\over1+cx}\,dx
 ={e^{1/c}E_1(1/c)\over c},\qquad J(0)=1.
\tag{S4}
\]

离开负实 c 轴时取与积分一致的主支 E₁。对负实 c，设 a=1/|c|：

\[
 J_\sigma(c)=a e^{-a}\operatorname{Ei}(a)
 -\sigma i\pi a e^{-a},\qquad c<0,\quad\sigma=\pm1.
\tag{S5}
\]

σ=+1 对应 c+i0。对 c>0，两支都为 (1/c)e^{1/c}E₁(1/c)。这条显式边值式避免在数值库的切割上默认选边。

从式 (2.33) 的交换和 h→−h，另一个相对时间 sector 的 c 为 he^{-t}/4，i0 相反。因此，原文 (2.25)—(2.26) 所指定的组合为

\[
 R_{12}(t)=\left\langle
 J_+\left({h e^t\over4}\right)+
 J_-\left({h e^{-t}\over4}\right)-1\right\rangle_\phi,
\tag{S6}
\]

\[
 R_{14}(t)=\left\langle
 J_-\left({h e^t\over4}\right)+
 J_+\left({h e^{-t}\over4}\right)-1\right\rangle_\phi
 =R_{12}(-t)=R_{12}(t)^*.
\tag{S7}
\]

这里 t 位于 V 的插入，故 sector 标签与“W 相对 V 的正负时间”相反。用原文公式编号标记 sector，避免只凭 superscript 正负误读。

有限正 ε 的独立边值诊断采用

\[
 R_{12,\epsilon}(t)=\left\langle
 J\left({h e^t e^{-i\epsilon}\over4}\right)+
 J\left({h e^{-t}e^{i\epsilon}\over4}\right)-1\right\rangle,
\tag{S8}
\]

R₁₄,ε 用 ε→−ε。对 h<0，第一项从 c 上半平面逼近，第二项从下半平面逼近；h>0 无切割。取 ε>0、ε≪1，并检查 ε 收敛。本式只被用作 τ=0 两候选轮廓的边值诊断，不声称简单主支表达式覆盖 F₁₂ 在整个欧氏位移条带中的继续路径。

两支都保留复数。循环置换 (2.35) 交换 R₁₂/R₁₄，不是证明某一支满足迹循环性。取它们的实部或平均不会自动恢复 §3 所要求的正性。

## 稳定冲击波、角积分及独立矩检查

利用角反射对称，〈f〉=(1/π)∫₀^πf(φ)dφ。原文 h 可稳定改写为

\[
 h(\phi)=-{\cos[k(\pi-\phi)]\over2m\,\mathrm{sinc}(\pi\gamma)},\qquad
 \mathrm{sinc}(x)={\sin x\over x},\quad 0\leq\phi\leq\pi.
\tag{S9}
\]

此写法在 G→0 时没有 0/0，直接给出 h₀(φ)=cosφ/(2m)。若 k>1/2，h 的零点为 φ₀=π−π/(2k)，角积分在此拆段；若 k≤1/2，h≤0 覆盖全部角区间。

独立几何积分应复现

\[
 \langle h\rangle=-{4\pi G\over\alpha},\qquad
 \langle h^2\rangle={2\pi^2G^2(\alpha+\sin\alpha)\over\alpha\sin^2(\alpha/2)}.
\tag{S10}
\]

第二式的稳定形式是 [1+sinα/α]/[8m²sinc²(πγ)]。G=0 时分别为 0 和 1/(8m²)。有限 G 的冲击波包含纯反冲，不应先删去 ℓ=1 模式。

(S5) 还给出不依赖 E₁ 分支库的虚部验证：

\[
 \operatorname{Im}R_{12}(t)=4\pi\left\langle
 {\mathbf 1_{h<0}\over|h|}
 \left[e^t e^{-4e^t/|h|}-e^{-t}e^{-4e^{-t}/|h|}\right]
 \right\rangle.
\tag{S11}
\]

## 扰动交叉检查与时间区域

将 (S4) 作渐近展开，而非假定收敛 Taylor 级数：

\[
 J(c)\sim\sum_{n\ge0}(-c)^n n!.
\]

由此得到

\[
 R\sim1-{\langle h\rangle\over4}(e^t+e^{-t})
 +{\langle h^2\rangle\over8}(e^{2t}+e^{-2t})+\cdots.
\tag{S12}
\]

它正好是原文 (2.42) 在 ν=ν′=i/2 的连续极限，那里看似 0×∞ 的第一个质量因子应取极限 1/2。纯反冲 G=0 时第一项消失，最先保留 (e^{2t}+e^{-2t})/(64m²)。

拟采用 m=1000，亏角参数 γ={0,0.05,0.25}。γ=0.25 对应 α=3π/2，与原文图 12 的几何相同；对应 G=6.25×10⁻⁵。令 t=log m+u，u 在 −4 到 2，覆盖小增强项至重求和必要的区域，另在验证中检查 t=0；t=0 不属于主 CSV 的 21 组扫描。原文的受控极限是 m 大、G 小而保留 Ge^{|t|} 与 e^{|t|}/m；上述数值是该有效公式的参数示例，不是对丢弃项的严格误差界。

只在 max|h|e^{|t|}/4≪1 的点比较 (S12)。增强项达到 O(1) 后使用完整核，不用截断级数下结论。固定 γ、m→∞、t 固定得到 R→1；固定 t−log m 则仍有非平凡单侧散射。只把 G→0 而固定 m 不会去掉反冲。

## 数值策略与成果边界

主实现采用任意精度的 E₁/Ei 边值核和拆段角积分。独立检查包括：普通有理核积分 (S4)、减去极点后的实主值积分 (S5)、角矩 (S10)、虚部实积分 (S11)、有限 ε 收敛、原文二阶系数和纯反冲连续极限。保存双轮廓值及 TOC/断开基准，不构造未经说明的实平均。

这一步推进到有限 G 的 eikonal 动力学复现，但没有引入第三阶段尚缺的有限内钟开放端点核。即使数值显示两轮廓不同或早时增长，也只能视为原文模型现象的复现，不能据此排除所有 de Sitter 全息方案。

## 已完成的数值复现

本轮已经实现 [conformal_shock.py](../code/conformal_shock.py) 与 [run_conformal_shock.py](../code/run_conformal_shock.py)，结果保存于 `results/stage4/`。主体采用 35 位精度，选点增至 45 位核查。共计算三组亏角、七个 t−log m 的 21 组双轮廓关联函数，并保存 12 组有限 ε 诊断。

**84 项检查通过，0 项失败。** 检查包括有理积分、去极点实主值、切割符号、断开归一化、五组冲击波角矩、独立留数虚部、两轮廓共轭、循环置换的时间翻转、ε 收敛、二阶渐近系数、纯反冲连续极限与固定时间大质量极限。

| 检查 | 最大绝对差 |
|---|---:|
| E₁ 核与独立有理积分 | 1.69×10⁻³⁶ |
| 冲击波解析角矩与积分 | 7.35×10⁻⁴⁰ |
| 关联函数虚部与留数实积分 | 1.89×10⁻³⁷ |
| 35 位主体结果与 45 位复算选点 | 4.07×10⁻³¹ |

最后一行包括 CSV/JSON 保存到 30 位有效数字的截断。这些误差只描述数值复现，不约束丢弃的非 eikonal 物理贡献。

m=1000、γ=0.25、G=6.25×10⁻⁵、α=3π/2 时的部分结果：

| t−log m | Re R₁₂ | Im R₁₂ |
|---:|---:|---:|
| −4 | 1.000770545856659 | −5.45×10⁻¹⁷⁰ |
| −2 | 1.005928245696152 | −9.37×10⁻²³ |
| 0 | 1.066039637322169 | −0.003285830982471 |
| 2 | 0.946095155927024 | −0.567374160149312 |

R₁₄ 的实部相同，虚部反号。早时增长和随后复数结构均是原文有效模型的复现。在后两个点不能用二阶渐近截断替代完整核；Re R 降低也不能单独当成恢复正性或完成混沌检验的证据。

![双共形有限引力eikonal基线](C:/Users/weiji/Documents/ChatGPT/de-sitter-holography/results/stage4/conformal-shock.png)

图中直线连接离散计算点，非额外拟合。虚线保留第二轮廓，未将两条结果取平均。图已打开核查。

固定 m=1000、t=log m，将 γ 从 10⁻²、10⁻³、10⁻⁴ 趋于零，R 与独立 γ=0 纯反冲值的绝对差依次为 1.22617×10⁻³、1.20856×10⁻⁴、1.20682×10⁻⁵，符合连续极限。这个纯反冲值仍不是 1。

复算命令：

```powershell
python code/run_conformal_shock.py
```

参数、版本和两份脚本 SHA-256 见 [study.json](../results/stage4/study.json)，逐项容差与状态见 [validation.json](../results/stage4/validation.json)，原始值见 [correlators.csv](../results/stage4/correlators.csv)、[regulators.csv](../results/stage4/regulators.csv)、[shock_moments.csv](../results/stage4/shock_moments.csv)。依赖为 Python 3.14.5、mpmath 1.3.0、matplotlib 3.11.0。

## 下一项字典比较的条件性线索

[Narovlansky 2506.02109v2](https://arxiv.org/html/2506.02109v2) 的 (6.5)—(6.7) 可在同一场质量处作局部比较。取 Δ=1/4、μ=1/2，并把两个对点差送到无量纲时间 x→−iπ，显式展开得到

\[
 g_\Delta(-i\pi)=g'_\Delta(-i\pi)=\tfrac18.
\]

这里 g′ 是原文另行定义的函数，不是导数。数值从 x=−iπ+ε 逼近也确认该极限。令微观无量纲时间 uᵢ=2πTtᵢ，取 u₁=iπ/2、u₃=−iπ/2、u₂=t、u₄=t−iπ，则 e^{4πTt₄₃}=−e^{2t}。因此原文交叉四点在除去尚未解释的整体因子 C₀=log(s₀)/(4π²T) 后，单侧首项为

\[
 1+{e^{2t}\over256s_0^2}+\cdots.
\]

与本文件 G=0 的单侧首项 1+e^{2t}/(64m²)+… 比较，形式上提示 m=2s₀。这只是同一插入配置下一个系数的条件性匹配；它还不确定 C₀、态密度、观察者能量、完整轮廓或有限 G 的 Ge^t 项。不能把它称为已建立的质量字典，尤其不能据此将内钟 p 替换成 s₀。下一步需查明微观极限中哪一个独立参数承载有限 G 贡献，并在未人为重定归一化之前比较完整观测量。

后续核验已扩至固定全阶形式系数和纯反冲全函数分支，见[微观比较专题](../literature/stage4-microscopic-match.md)。其中明确区分了早时 leading-log 展开与扰乱时间附近的重求和，并保留整体归一化和物理质量字典的缺口。
