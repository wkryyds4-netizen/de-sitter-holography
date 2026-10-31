# dS₃ 同种粒子：Bose 归一与小时间非零性

2026-10-31。Alan Day / Chengheng，内部研究稿。[LaTeX](../manuscript/ds3-bose/main.tex)、[五页PDF](../output/pdf/de-sitter-ds3-bose.pdf)、[来源核查](../literature/stage4-ds3-bose-sources.md)。接续[真实交叉波包](stage4-ds3-crossed-packet.md)。

## 模型、物理空间与归一

ℓ=1、G=0、proper SO⁺(3,1)，两个可区分BD标量观察者质量标签≥s，物质为一个KG归一自由标量species p。保留真实质量h(x)=√(1+x²)、局部耦合1、半周期调节创建块。数值s=.2,p=.6；初态ρ|f|²是[.4,.8]上的Beta(5,5)，f实非负、零延拓C¹。计算仍保留无限内部质量谱。

先将两个同质量物质临时标号。完整四体物理空间上，交换物质槽的F与对角群作用交换；在有限K、紧质量测试的收敛群平均上保持正半定form和零空间，因而下降为物理完成空间的酉自伴对合。P+=(1+F)/2是正交投影，Bose空间为其像。这是完整四体空间的陈述，不要求单个标量Racah子块酉。

有限模中，局部创建写作Σ M_j(t) a†_j，标号链C_jk=M_j(t₂)M_k(t₁)f。从occupation基的√(n_j+1)直接推出B=√2 P+C。故物理两创建的配对为

\[
\langle B_k,B_l\rangle=\langle k,l\rangle+\langle k,F l\rangle.
\]

范数通常不是2||k||²；达到上界要求Fk=k。观察者矩阵一般不交换，不能假定标号链已对称。√2是标准Fock规范，不另添2!、Haar或shadow反射因子。已知调节后标号1→2块有界，因此Bose块范数≤√2倍；不据此声称完整Fock共同域已建立。

取旧规范e_b=√W_b J_b、W_b=ρ(b)K_abp K_bcp，J_b(0000)=1。F把s配对映到t配对，已有物理U_00(b,d;a,p,c,p)=〈e_b,F e_d〉。此核为Hermitian contraction的矩阵元，标量压缩通常不平方为1。

## 两时间次序差

共同实振幅T_b=f(c)E_abp E_bcp exp{π[h(a)−h(c)]/2}/√W_b。令α_b=exp{−it[h(a)−h(b)]}，β_b=exp{−it[h(b)−h(c)]}，k=Tα，l=Tβ。两Bose态是√2P+k和√2P+l。

N=||k||²=||l||²，Y=〈k,l〉，E₁=〈k,Fk〉，E₂=〈l,Fl〉，Z=〈k,Fl〉。于是N_B1=N+E₁，N_B2=N+E₂，Bose overlap=Y+Z，D_B=2N+E₁+E₂−2Re(Y+Z)。交换项用完整复U及旧物理权重，三种相位分别为exp{it(h_d−h_b)}、其逆及exp{it(h_a+h_c−h_b−h_d)}。

稳定计算用e_b^en=2h(b)−h(a)−h(c)，δ_b=2i exp{−it[h(a)−h(c)]/2} sin(te_b^en/2)：

\[
D_B=\|T\delta\|^2+\langle T\delta,F T\delta\rangle.
\]

它在t=0严格为零，不靠两个有限大数相减。两创建扇区不含四个自伴场中的中间真空收缩，不能混称完整四场关联。

## 模型内小时间定理

紧初态及hatR≤exp(−π|a−b|)/(4π)、g(a)≤a+C₀保证T的全部多项式质量矩有限。定义j=e^en T，用正弦三阶余项在Hilbert范数中得到

\[
D_B(t)=t^2Q_B+O(t^4),\quad Q_B=\|\sqrt2P_+j\|^2\ge0.
\]

外部共同相位与F交换，故D_B为偶函数。该推导由六阶矩控制余项，不依赖数值Racah求值。

严格正性使用前轮真实顶角E_abp=Cp sqrt(coshπb/coshπa) H_p(a,b)，H_p(a,b)=r_p(a−b)−r_p(a+b)，r_p(x)=1/(coshπx+coshπp)。全部外部球模为0时，j的系数除去正共同因子后为

\[
J(a,c)=\int_s^\infty H_p(a,b)H_p(b,c)[2h(b)-h(a)-h(c)]\,db.
\]

a=c>s时，将b=c−x与c+x在0<x<c−s配对：r严格递减给H(c,c+x)>H(c,c−x)>0，h严格凸增给h(c+x)−h(c)>h(c)−h(c−x)>0。故每一配对贡献为正，余下b>2c−s也为正，因此J(c,c)>0。由连续性，c取(.4,.8)中的任意紧区间I时，存在统一ε>0使|a−c|<ε内仍正。选质量系数f(c)𝟙{c∈I,|a−c|<ε}的全球模测试向量，其群平均范数有限，与√2P+j的配对严格正。物理Cauchy于是给Q_B>0，进而小的非零t有D_B(t)>0。

这是存在性证明，没有给ε、时间邻域或Q的区间认证数值。与前轮不同species在t=0的谱下界见证不是同一定理，也未声称有限G混沌。

## 解析谱尾和自由基准

沿用前轮q=p的ε_L与Nmax。对a,b≤L的标号链截断，直接重叠尾≤ε_L、交换尾≤2√(ε_L Nmax)。D_B的直接尾≤4ε_L、交换尾≤8√(ε_L Nmax)，总界是其和。L12给ε≈1.16352690246e−16、Nmax≈.00136455222524、D尾≤约3.18767e−9。表达式为解析界，小数为浮点计算；不包含盒内、series或舍入误差，也不优化t→0。

固定观察者自由场基准为||a†(u_t)a†(u₀)Ω||²=L_p(0)²+|L_p(t)|²，L_p(t)=sin(pt)/(4πsinhπp sinh t)，D_free=0。零时为2L_p(0)²。它是标准Fock Gram permanent；当前尚未证明整个dS Bose交换项在重观察者极限趋向它，不能仅引用旧非交叉极限。

## 数值、失败与独立核查

[三网格及重放](../results/stage4/ds3-bose/summary.json)共1,216,480个复kernel点，NPZ全部保留。最终L12、每段6节点、初态12节点、series N128 K12。另有同L12/nc12/series但4节点网格，以及L8/nc8/4节点pilot。固定其他参数4→6的D(2)差2.26752e−8，完整Bose重叠最大差1.64309e−8，均大于解析谱尾，不能作总误差认证。

最终保存值N_B(0)≈2.2345037479e−4，Q_B≈3.1755764941e−5；D_B(.5)≈7.70115634e−6，D_B(1)≈2.83756540e−5，D_B(2)≈8.95584228e−5。展示位数不等于认证精度。t=.001,.01的D/t²趋近Q，t0稳定公式给严格0。

所有网格以原始E/W/ρ独立重放，误差<3e−17；完整复值、Cauchy、正性、展开式与稳定式一致，最终加权Hermiticity残差约6.77e−18。Hermiticity小不能排除对称实部误差。最大权重condition proxy的六个节点做55位比较，最高相对误差7.151e−11；只代表这些点。

[初始诊断](../results/stage4/ds3-bose/validation.json)22项中21项通过、1项失败，原样保留。有限Fock矩阵、C3交换下降和45位球面配对积分均符合公式；C3仅为独立有限模型核查。预期普通的(a,b,c,d)=(.7,.5,.6,.8)恰在a+c=b+d近旁，旧complex128结果相对45位参考错约160.15%，Hermiticity残差约.926。另一个远端最大condition点相对错约.128%。

[精度修复](../results/stage4/ds3-bose/precision-repair.json)对选定病态点用55/75位比较，恢复Hermiticity；[新helper](../code/ds3_bose_precision.py)对condition>1e10的单点回退。该helper未静默改写网格，不能直接当高质量全网格认证算法。最终网格所检a+c−b−d离零最近约1.88e−5；其余潜在相消面仍需审计。下一步重质量前必须处理统一稳定求值，不能只把complex128升采样。

本轮标准Fock归一、对称化与Gram permanent不作原创声明；新模型推导也未作穷尽查新。下一主线：稳定交叉重质量极限及同species自由基准，继而全Fock与有限G。
