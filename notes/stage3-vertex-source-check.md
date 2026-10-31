# A.5 三模积分：独立数值核验与误差记录

日期：2026-10-29。本记录续接[源公式审查](stage3-cg-domain-investigation.md)，但不修改该记录及其证据哈希。原文为 [Kolchmeyer–Liu 2411.08090v1](https://arxiv.org/html/2411.08090v1)。本轮直接计算模式函数积分，不由 A.7 反推积分常数。

## 计算的对象与实现检查

取 (s₁,s₂,s₃)=(0.4,0.7,1.2)，n₁=n₂=n₃=0。A.6 给 e^{iφ₀,s}=i，固定半相位 e^{iπ/4} 后，A.5 的 Gamma 比、正则化超几何 c=1 的因子均消去：

\[
\psi_0^s(t)=\frac{{}_2F_1(1/2-is,1/2+is;1;(1-i\sinh t)/2)}{2\sqrt{\cosh\pi s}},
\quad \Psi_0^s=\sqrt{\frac{s\tanh\pi s}{\pi}}\psi_0^s.
\]

角积分给 2π；由 ψ₀(−t)=overline(ψ₀(t))，三模积分是

\[
I=4\pi\operatorname{Re}\int_0^\infty\cosh t\,
\overline{\Psi_0^{s_1}(t)}\Psi_0^{s_2}(t)\Psi_0^{s_3}(t)dt.
\]

在每个质量及 t=0,0.7,2 上，分别检查原 (3.4) 的 KG 范数与 KG 微分方程；40 位工作精度下最大误差为 4.60×10⁻⁴¹ 与 8.12×10⁻⁴²。超几何导数由相邻参数公式计算。这里不把工作精度等同于所有特殊函数的真实精度。

CG 先用两个途径计算：原 E.22 后的 C_S 显式两项正则化 4F3；以及我们从 E.15/E.17 导出的零 Fourier 卷积。设 α=Δ₁+Δ₃−Δ₂、β=Δ₂+Δ₃−Δ₁、γ=Δ₁+Δ₂−Δ₃，

\[
f_\alpha(0)=\frac{2^\alpha\Gamma(1-\alpha)}{\Gamma(1-\alpha/2)^2},\quad
C_{S,000}=\prod_{x=\alpha,\beta,\gamma}f_x(0)
\left[2\,{}_4F_3\left(1,\alpha/2,\beta/2,\gamma/2;1-\alpha/2,1-\beta/2,1-\gamma/2;1\right)-1\right].
\]

该卷积的项为 O(k⁻³ᐟ²)，绝对收敛，但短直接截断不能当高精度参考。按 E.19 用 Δ=(1/2−is₃,1/2−is₂,1/2+is₁)，再按 E.24 乘 e^{-iπ/4}。两种 4F3 表达式在40位计算中彼此相对差3.09×10⁻⁴¹。

**独立参考揭示了共同误差。** 后来由 [Bros 等 0901.4223v1 (124)](https://arxiv.org/html/0901.4223v1)取得纯 Gamma CG 参考，发现上述两种 4F3 的共同相对偏差约1.28×10⁻²⁷。故“彼此一致到40位”不能解释为真实40位准确。这一偏差远小于下述时间截断误差，但已原样记录，没有覆盖早期结果。Gamma 参考与解析映射见[解析审查](stage3-vertex-analytic-audit.md)。

## 截断、求积及比较

对正时间按宽度2分面板，用任意精度 Gauss–Legendre。第一档工作精度30位、每面板24点，记录 T=16,32,48,64,80；第二档40位、每面板36点，独立重算 T=80。两档 T=80 结果绝对差2.40×10⁻³¹。该阶数/精度比较是经验数值检查，不是严格求积区间包络。

| 量 | 结果 |
|---|---|
| 直接 A.5 积分，T=80 | 0.0111785603498889046822781311122602 |
| literal A.7 的 D.8 预测 | −0.000506744653141198668085 − 0.001537009468258153025657 i |
| 仅把末 Γ 移到分母、用纯 Gamma CG 的预测 | 0.0111785603498889045844547505661441 |
| literal 预测相对差 | 1.05433573768 |
| Gamma 参考与直接积分的绝对差 | 9.78234×10⁻²⁰ |
| T=80 的解析双尾界数值 | 8.49373×10⁻¹⁹ |

后两者相容。反事实版本在本样本的数值支持不单独证明唯一修正或全标签正确性；全标签 n=0 的解析证据另记。literal 预测连应有的实相位也不满足，且此前已有不依赖数值的输入交换矛盾。

## 解析尾包络

令 z=(1−i sinh t)/2、a=1/2−is、b=1/2+is。2F1 的无穷远连接式由两个 (-z)^{-a}、(-z)^{-b} 分支及 1/z 的 2F1(a,a;1−2is;1/z) 类级数组成。其系数相邻绝对比满足

\[
\left|\frac{(a+k)^2}{(1-2is+k)(k+1)}\right|
\le\frac{(k+1/2)^2+s^2}{(k+1)^2}\le B_s:=\max(1,1/4+s^2).
\]

设 θ_T=arg[−(1−i sinh T)/2]，D_s=|Γ(2is)/Γ(1/2+is)²|，当 2B_s/cosh T<1 时，所有 t≥T 都有

\[
|\Psi_0^s(t)|\le L_s(T)\cosh^{-1/2}t,\quad
L_s(T)=\sqrt{\frac{s\tanh\pi s}{\pi}}
\frac{\sqrt2D_s\cosh(s\theta_T)}{\sqrt{\cosh\pi s}\,[1-2B_s/\cosh T]}.
\]

这里 |z|=cosh t/2，θ_t随 t下降到π/2。由 cosh t≥e^t/2，双尾严格满足

\[
|I-I_{[-T,T]}|\le8\sqrt2\pi\,\prod_jL_{s_j}(T)e^{-T/2}.
\]

不等式本身是解析结果；文件中的数值以高精度计算该上界，没有声称进行了向外舍入区间算术。

## 可复算文件

- [独立时域与双 CG 脚本](../code/validate_vertex_source.py)：`--quick` 输出早期 pilot；不带参数重做30/40位两档。
- [冻结的完整原始数值](../results/stage3/cg-domain-source-audit/vertex-check.json)与[早期 pilot](../results/stage3/cg-domain-source-audit/vertex-check-pilot.json)。
- [纯 Gamma 参考脚本](../code/validate_vertex_gamma_reference.py)及[参考结果](../results/stage3/cg-domain-source-audit/vertex-analytic-reference.json)，明确记录4F3的共同偏差。
- 运行环境为 mpmath，当前复用本机 blackhole 项目的 Python 运行时；没有向 blackhole 写入研究文件。