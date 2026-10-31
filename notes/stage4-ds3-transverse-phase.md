# 完整核的横向相位与逐点对数尺度

2026-10-31。Alan Day 与 Chengheng；[四页PDF](../output/pdf/de-sitter-ds3-transverse-phase.pdf)、[LaTeX](../manuscript/ds3-transverse-phase/main.tex)、[来源](../literature/stage4-ds3-transverse-phase-sources.md)。本轮为已知核在指定模型中的推导，不作首创性声明。

## 模型、坐标及结论

ℓ=1、G=0、proper SO⁺(3,1)、BD标量观察者、KG物质、固定正观察者标签下界与半周期调节创建块。U₀₀保留物理CG球周期相位；标量矩阵元不是独立酉变换。M为principal-series标签，真实静能量√(1+M²)。没有做质量积分或能量Jacobian替换。

取 c=M、b=M+u、a=M+u+v、d=M+v+r，Aⱼ=1/2+iαⱼ，α=(v+q,p−u,−u−p,v−q)/2，要求p,q>0、|αⱼ|≤1/2、|αⱼ+r/2|≤1/2。定义从δ(0)=0连续取值的相位

\[
\delta(r)=\Im\left[2\log\Gamma(1-ir/2)+
\sum_j\bigl(\log\Gamma(A_j+ir/2)-\log\Gamma(A_j)\bigr)\right].
\]

在固定非零r及去掉r=0的紧转移集合上，

\[
U_{00}=\frac{\sin(r\log M-\delta(r))}{\pi r}+O(M^{-1}).
\]

对r=0窗口严格内部的固定转移，每个固定实λ满足

\[
\frac{\pi U_{00}(r=\lambda/\log M)}{\log M}\longrightarrow
\frac{\sin\lambda}{\lambda},
\]

λ=0连续取1。单个Omega之前的sinc(λ/2)不是完整核的极限；物理归一相位贡献了另半个振荡相位。逐点缩放不意味着一致收敛或分布极限。

## 推导中的关键结构

令z=−ir/2。LPRS两级数在Wilson反射之后有精确共同因子：X=c_AF+c_BG、Y=y_AF+y_BH。F在两处的四个上参数和三个下参数完全相同，符号比较已保存。Y的两个Gamma系数见稿第2页。

D₀=C₊c_A=D∞M^z[1+O(1/M)]。其完整CG相位的大参数部分为−rlogM/2−r(3r+4v)/(8M)+O(M⁻²)，小参数Gamma相位不应丢弃。将它与y_A、y_B合并后，一支为T(r)e^(−irlogM)，另一支为J(r)。Gamma反射和模恒等式给T=i e^(iδ)/(2πr)、J纯虚；因此2Re后得到上述正弦形，而非单因子轮廓。相位、符号及共轭对的完整清单见audit.json。

设S=y_A+y_B；沿用小上参数比值估计得到|F−1|≤R_F、|Y−S|≤R_Y。大参数伴随G用sinh乘积和精确十次首一多项式给整个和的B_G界，保留N*=O(M⁸)辅助阈值，不实际求和到该阈值。完整有限M界为

\[
|U-2\Re(D_0S)|\le2|D_0|[R_Y+R_F(|S|+R_Y)]
+2|E_0|B_G(|S|+R_Y).
\]

这是Gamma系数精确保留的近似误差界，不是最终正弦式的显式有限M误差界。固定r≠0的级数误差O(M⁻²)，Gamma渐近另带O(M⁻¹)。r=λ/logM、固定λ≠0时系数为O(logM)，相应误差仍消失；λ=0用上一峰面证明单独处理。

δ′(0)=−Re[2ψ(1)−Σψ(Aⱼ)]/2，故正弦式的可去端点与峰面logM及常数匹配。但R_Y直接界含1/|r|，尚无穿过r=0的一致误差；不能对物理波包先换δ再积分。

## 数值证据与局限

[16组数据](../results/stage4/ds3-transverse-phase/study.json)：(u,v,p,q)=(.1,−.1,.4,.4)，M=16,64,256,1024，r=±.08与±1/logM，90位、至少128项。全部满足有限M界，另计保存的系列尾。203项独立Decimal/Fraction及符号核查通过。

[四个120位/双截断复核](../results/stage4/ds3-transverse-phase/refinement.json)通过；M16两点相对变化约1.63e−23、9.29e−24，M1024两个缩放点约1.52e−86、1.42e−86。使用1e−75绝对舍入allowance，仅为经验容差，没有区间认证。

λ=1时πU/logM依次约.9225、.9975、.9977、.9853，极限sin1约.84147：收敛慢且这一段不单调。保留全部值，不据中等M图形宣称已达到极限。小M显式界也很保守。

## 当前继续位置

下一步合并两个分支的可去奇点，建立穿过r=0的联合误差，再对真实创建幅和谱测度控制振荡尾、质量积分及极限交换。固定条件质量、初态包、移动下界、增长时间分别处理。Bose自由极限、完整Fock、有限G和全息微观字典仍开放。

旧47841本轮明确exit 0，12组100位比较及独立collision审计通过，当前没有需续查的旧句柄。旧稿中的pending是当时检查点状态，见[完成补记](../results/stage4/ds3-collision/completion.md)。
