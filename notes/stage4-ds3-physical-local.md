# 真实谱权重与局部物理配对

2026-10-31。Alan Day / Chengheng。[四页 LaTeX](../manuscript/ds3-physical-local/main.tex)、[PDF](../output/pdf/de-sitter-ds3-physical-local.pdf)。沿用[跨零局部核](stage4-ds3-uniform-local.md)，先核[两物质](stage4-ds3-two-matter.md)、[交叉制备](stage4-ds3-crossed-packet.md)及[Bose](stage4-ds3-bose.md)的实际幅度；[来源](../literature/stage4-ds3-physical-local-sources.md)。

## 模型与精确测度

G=0、ℓ=1、proper SO⁺(3,1)、两个可区分BD标量观察者、KG物质p,q>0、固定标签下界s>0、半周期调节创建块、局部耦合1。初态条件化在c=M，由归一窄包极限定义。M是表示标签，真实静能量h(x)=√(1+x²)。本轮不是将[.4,.8]固定初态包送到重质量。

原始E/W/ρ组装为

\[
A_M=\frac{\overline{E_{abq}E_{bcp}}E_{adp}E_{dcq}e^{\pi(h_a-h_c)}}{\rho(c)\sqrt{W_sW_t}},\quad
W_s=\rho(b)K_{abq}K_{bcp},\ W_t=\rho(d)K_{adp}K_{dcq}.
\]

约分后恰为g(a)/g(c)乘四个hatR的平方根，公式见PDF(1)。条件质量积分使用da db dd；不应再乘ρ(b)ρ(d)。条件化前为π₀(c)dc，其中π₀=ρ|f|²。转为c=M、b=M+u、a=M+u+v、d=M+v+r时两个相关行列式绝对值均为1。若改用三个真实静能量，测度另有∏h(x)/x，x=a,b,d；该因子趋1但有限M不能省略。

## C²权重与局部极限

令F_M=A_M exp[it(h_a+h_c−h_b−h_d)]。任意固定紧转移、正p/q紧集及|t|≤T上，

\[
F_M=(1+(u+v)/M)F_\infty+O_{C_r^2,T}(M^{-2}),\quad
F_\infty=\sqrt{k_q(v)k_p(u)k_p(u-r)k_q(v+r)}e^{-itr}.
\]

证明直接使用hatR/k=1−[coshπ(x−y)+coshπp]/[coshπ(x+y)+coshπp]的函数及两阶导数指数界、δ=h−x的|δ′|≤1/(2x²)、δ″≤x⁻³，以及g比的真实Taylor式。没有形式微分未知余项。固定下界在足够大的M下不碰到该紧窗口。

对固定C²紧截取χ，若其支撑还满足旧核条件max|αⱼ|+R/2≤1/2，则

\[
\int\chi F_MU_{00}\,du\,dv\,dr
=\int\chi(u,v,0)k_p(u)k_q(v)\,du\,dv
+O_{\chi,T}(1/\log M)+O_{\chi,T}(1/M).
\]

这是先前局部定理对真实权重的推论，包含u/v的固定紧积分；没有去掉截取。核定理仍依赖项目内前轮证明，未宣称独立外部同行验证。

同species的另一交换相位exp[it(h_d−h_b)]在r=0给exp[it(v−u)]，并非常数。将来若全局扩展成立，这两种配对分别产生L_p(0)²与|L_p(t)|²。两者不能混淆；目前不能因此宣布完整Bose自由极限。

## 计算证据与限制

[diagnostics](../results/stage4/ds3-physical-local/diagnostics.json)保存85位算术、12个实际E/W/ρ比较、两种相位的72个导数比较及12个加权积分重组，共97项通过。M4/16/64/256，含不等p/q及非零u/v/r；原始表示最高相对差约2.95e−83。有限M导数的M²缩放残差保存为观察量，没有设置经验“证明阈值”。

加权积分复用旧144个完整复U节点，绝不是144个新核计算。固定(u,v,p,q)=(.1,−.1,.4,.4)，R=.25，χ_r=(1−(r/R)²)³区间外0，24/48阶、t=0/1/2。48阶零时相对于极限的比值仅约.470982/.552526（M64/256），t2则约.465087−.002135i/.545862−.002340i；有限质量尚远离极限。舍入、求积与解析series尾分开记录。

另以不导入mpmath、SymPy或权重实现的Decimal脚本重放实际权重和加权积分，见[audit](../results/stage4/ds3-physical-local/audit.json)。计算并保存24/48网格差。这是独立权重重放，不是重新认证输入U。

## 下一步的真实障碍

原小窗口条件不能穷尽整个(u,v)平面，也限制当前可覆盖的p/q。即使窗口外有正范数尾界，也不能将这个固定小窗口的补集直接称为可忽略尾。下一步须扩大完整核在任意固定紧转移上的一致渐近，再用真实Hilbert范数裁剪两条链并控制其配对差；标量U是压缩核，不能从算符收缩性推出逐点|U|≤1。任意χ也不自动分解为两条链的独立乘子。

固定初态包、移动阈值、增长时间另行处理。全Fock、有限G与完整微观全息目标仍开放。
