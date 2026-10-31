# 完整交叉核的精确峰面极限

2026-10-31。内部阶段稿，Alan Day 与 Chengheng；[LaTeX](../manuscript/ds3-peak-limit/main.tex)、[五页 PDF](../output/pdf/de-sitter-ds3-peak-limit.pdf)、[来源](../literature/stage4-ds3-peak-limit-sources.md)。以下为指定模型内的推导，不作原创性或完整全息对偶声明。

## 定义与结论

沿用曲率半径 ℓ=1、G=0、proper SO⁺(3,1)、BD 标量观察者、KG 物质、固定正观察者标签下界与半周期调节创建块。U₀₀ 是保留物理 CG 球周期相位的标量矩阵元，不是独立酉标量变换。M 是 principal-series 标签，静能量为 √(1+M²)，未做能量 Jacobian 或质量积分。

仅在精确峰面 c=M、b=M+u、a=M+u+v、d=M+v 上，取 p,q>0、|u|+p≤1、|v|+q≤1。令 Aⱼ=1/2+iαⱼ，α=(v+q,p−u,−u−p,v−q)/2，得到

\[
U_{00}=\frac{\log M}{\pi}
+\frac{2\psi(1)-\Re\sum_j\psi(A_j)}{2\pi}
+\frac{u+v}{2\pi M}+O(\log M/M^2).
\]

余项对满足上述条件的紧转移集合一致。它只描述峰面；峰高不能决定横向宽度、积分权重或 δ 分布极限。

## 证明路线及可复查对象

1. 由上一阶段 Wilson 反射，加上实际三球周期的 Gamma 相位，得到更强恒等式 Ĉ₋=conj(C₊)，故 U=2Re(C₊XY)。六个独立实标签的 Gamma argument/weight 列表及共轭配对证书保存在 study.json 的 symbolic 字段。这是解析恒等式，未用数值实部投影消除误差。
2. 峰面 X=c_AF+c_BG。完整 leading 系数 D₀=C₊c_A=e^{iθ}/(4π)，其中 θ=[−p²v+q²u+u²v−uv²+u−v]/(4M²)+O(M⁻³)。小参数相位共轭抵消，移位 Stirling 展开中的 MlogM、M、logM、常数及 1/M 项全部抵消；连续相位选为趋零。代码 ds3_peak_normalization.py 保存符号核查。
3. 大参数伴随级数 G 不能仅保留首项。sinh 无限乘积给逐项一致界 B_b，十次首一项比多项式给 N*=2+ceil(max|pⱼ|)，进而 |G|≤(3N*+3)B_b。N*=O(M⁸) 仅是解析阈值，数值没有求和至 N*。结合 |C₊c_B|=O(M⁻¹/²e⁻πM)，得到 O(M⁸e⁻πM/²) 的保守整体抑制。
4. 令 q₀=ψ(1−ib)+ψ(1+id)+2ψ(1)−Σψ(Aⱼ)，R_F=(ζ(2)−1)/(1+min(b,d)²)，R_Y 为前期碰撞 digamma 级数余项。稿件式 (10) 给出完整有限 M 绝对界，再由 digamma 展开得到上式。解析界以精确算术为前提，浮点求值没有区间舍入认证。

## 数值、失败与修复

[study.json](../results/stage4/ds3-peak-limit/study.json) 保存两组 (u,v,p,q)=(.3,−.2,.6,.6)、(.2,.1,.7,.8)，每组 M=8,16,32,64,128,256；100 位、128 项，完整复数和解析截断尾分别记录。12 组均满足有限 M 界（另计截断尾）。第一组 M=8 的误差约 .002807 而界约 11286；M=256 的误差约 5.77e−6、界约 6.36e−5。保留小质量界过宽的事实。M²θ 分别趋向 .1625、.04525；数值趋势不能替代证明。

初始独立审计导入仅约 70 位的 π，却要求 85–90 位比较，导致失败。[原审计](../code/audit_ds3_peak_limit.py)及[失败输出](../results/stage4/ds3-peak-limit/audit.json)原样保留。新审计用 180 位 Decimal 的独立 Machin 公式，atan(1/5) 116 项、atan(1/239) 34 项，交错级数截断界约 9.93e−165。未改科学数据或放宽容差，[修复审计](../results/stage4/ds3-peak-limit/audit-refined.json)的 191 项检查通过。该审计核保存算术、精确有理阈值及相位证书，不认证全局舍入。

## 继续位置

下一步先建立离开 a+c=b+d 的紧横向转移一致渐近，再控制振荡尾与真实创建幅上的积分及极限交换。固定条件质量、固定初态包、移动制备与增长时间分别处理。Bose 自由观察者极限、全 Fock、有限 G 与微观全息字典仍开放。

旧参数导数验证 session 47841 续查仍运行，controlled-validation.json 未生成，未计为通过、未重启。四个新增阶段尚待项目整体登记和冻结；最新冻结仍是 ds3-bose。
