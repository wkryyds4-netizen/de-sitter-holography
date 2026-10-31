# dS₃ 群配对与有限质量两点

2026-10-31。Alan Day 与 Chengheng 的内部研究记录。对应[七页署名稿](../manuscript/ds3-group-pairing/main.tex)及[PDF](../output/pdf/de-sitter-ds3-group-pairing.pdf)。本轮为指定模型的推导与经验数值验证；前驱见[来源记录](../literature/stage4-ds3-group-pairing-sources.md)。

## 1. 模型、规范与任务

半径 ℓ=1，Newton G=0，proper 群为 SO⁺(3,1)。两个可区分标量 BD 观察者，主系列标签 ν≥s*>0；一类固定主系列物质标量 p>0，BD 真空；局部耦合取1。观察者质量包以 dν 积分，时空归一模式为 vνlm=√μ(ν)uνlm，μ=ν tanh(πν)/π。物质仍用 KG 模式。本轮继续[质量谱](stage4-ds3-mass-spectrum.md)和[局部顶角](stage4-ds3-local-vertex.md)，没有改变它们的冻结文件。

目标是从群配对得到物理范数，再决定态及局部两点。时空 μ、群密度 ρ_G、三点分布的逆范数是三种不同对象。

## 2. Haar、反射与双观察者物理空间

令 K=SO(3) Haar 概率测度，定义总体规范
\[
\int_G F(g)dg=\int_0^\infty dr\,\sinh^2r\int_Kdk_1dk_2F(k_1a_rk_2).
\]
单位球函数
\[
\phi_\nu(r)=\frac12\int_{-1}^1(\cosh r-x\sinh r)^{-1-i\nu}dx
=\frac{\sin\nu r}{\nu\sinh r}.
\]
正弦正交给正标签上的分布
\[
\int_G\phi_a\phi_b\,dg=K_2(a)\delta(a-b),\qquad
K_2=\frac{\pi}{2\nu^2},\qquad \rho_G=K_2^{-1}=\frac{2\nu^2}{\pi}.
\]
ν>0 已选一个反射代表；不再加额外2。此密度已有 YYZ (2.25) 直接前驱。

BD K 基中的不变分布为
\[
I_\nu=\sum_{lm}(-1)^{l+m}|\nu,lm\rangle_L|\nu,l,-m\rangle_R.
\]
Θ|lm⟩=(-1)^(l+m)|l,−m⟩ 并共轭系数，Θ²=1。旋转由 singlet 恒等式；boost 的实相邻 l 矩阵元经 l 奇偶反号，和反幺正共轭合起来给群 intertwiner。

种子域取两个质量标签内紧支撑光滑、远离0和阈值、有限 K 型。先做质量积分。紧质量 Harish–Chandra 渐近的一阶两体乘积为 e^(−2r) 乘 exp[i(±a±b)r]，分别对质量积分分部得到任意逆 r 幂；余项另有指数衰减。乘 Haar 后可积。有限 K 域本身不假称群不变，可取其群平移张成域。固定标签只能用分布配对，不能先要求普通积分收敛。

设 cψ(ν)=Iν†ψ(ν,ν)，群平均给不变分布系数 f=K₂cψ，故
\[
\langle\eta\psi,\eta\chi\rangle
=\int K_2\overline{c_\psi}c_\chi\,d\nu
=\int\rho_G\bar f g\,d\nu.
\]
这里 Ψf=∫fIνdν 的系数空间为 H₀=L²(ρ_Gdν)。平滑对角数据可由双质量种子实现，紧支撑系数在该空间稠密。不存在 δ² 操作，也不交换 rigging 权重与不变分布范数。

## 3. 三粒子范数与完整角扇区

固定 a,b,p>0。每粒子有限 K 矩阵元为 O((1+r)^N e^(−r))，三体乘积与 Haar 相乘可积。泛型标量主系列不变三线性形式唯一，取 I_abp 的000分量为1。于是 η₃=K₃|I⟩⟨I|。记
\[
S=a+b+p,\quad D=\cosh(\pi S/2)\prod_{n=a,b,p}\cosh[\pi(S-2n)/2].
\]
由三正弦积和标准 ∫sin(kr)/sinh r dr=(π/2)tanh(πk/2)，得到
\[
K_3=\int_0^\infty\sinh^2r\,\phi_a\phi_b\phi_p\,dr
=\frac{\pi}{8abp}\left[\sum_n\tanh\frac{\pi(S-2n)}2-\tanh\frac{\pi S}2\right]
=\frac{\pi}{8D}\prod_n\frac{\sinh\pi n}{n}>0.
\]
不变分布系数 h(a,b)I_abp 的物理范数是
\[
\|h\|_1^2=\int_{s_*}^\infty da\,db\,\frac{|h(a,b)|^2}{K_3(a,b,p)}.
\]
这覆盖完整角不变扇区；000 只选基准规范，不是将物质截成 l=0。独立非球面对角外群积分核查见第7节。

## 4. 线性左右配对态

已有局部顶角的球分量（出观察者 a、入 b、出物质 p）为
\[
B_*=\frac{\pi}{4\sqrt{4\pi}D}\sqrt{\prod_n\frac{\sinh\pi n}{n}},\qquad
E_{abp,000}=\sqrt{\mu(a)\mu(b)}B_*\cosh\pi b.
\]
左作用后不变系数为 f(b)E_abp,000，右作用为 f(a)E_bap,000。三点唯一性及 singlet 相位把其余角分量接上。相等要求
\[
f(b)\cosh\pi b=f(a)\cosh\pi a,\qquad f_{\rm lin}(\nu)=\sech\pi\nu
\]
（整体常数除外，几乎处处）。KL dS₂ 已有同一函数；这里重推 dS₃ 局部模型。称“线性配对态/候选 HH 态”；不将线性关系升格为任意复合算符 gluing 或完整 HH 代数。

\[
Z(s)=\int_s^\infty\frac{2\nu^2}{\pi}\sech^2(\pi\nu)d\nu
=\frac{4s^2}{\pi^2}\frac q{1+q}
+\frac{4s}{\pi^3}\log(1+q)-\frac{2}{\pi^4}\mathrm{Li}_2(-q),\quad q=e^{-2\pi s}.
\]
Z(0)=1/(6π²)仅形式端点检查，实际 s*>0。概率为 ρ_G sech²/Z。大标签权重约 (8/π)ν²e^(−2πν)，不是把有限标签ν当作真实质量。

## 5. 创建块、定义域和两点

一物质范数给
\[
Q(a,b;p)=\frac{|f_{\rm lin}(b)E_{abp,000}|^2}{K_3}
=\frac{\mu(a)\mu(b)}{8D}>0.
\]
Q 对 a,b 对称。定义 A:H₀→H₁，(Af)(a,b)=f(b)E_abp,000，则
\[
\|Af\|_1^2=\int db\,\rho_G(b)|f(b)|^2\mathcal S_p(b),\quad
\mathcal S_p(b)=\frac{\cosh^2\pi b}{\rho_G(b)}\int_{s_*}^\infty Q(a,b;p)da.
\]
因 D≥exp[2π max(a,b)]/16，Q≤2ab exp[−2π max(a,b)]/π²。把 a 积分在 b 分开，再用 cosh²πb≤exp(2πb)，得到
\[
\mathcal S_p(b)\le\frac b{2\pi}+\frac1{2\pi^2}+\frac1{4\pi^3s_*}.
\]
有限一阶质量矩足够；最大域则是上式精确 S 加权可积域。该映射是直接积分中的纤维乘法，图收敛子列的几乎处处极限仍为 f(b)E，故闭；紧质量支撑给稠密域。伴随定义反向块。这只完成0↔1，未证明所有N、共同核心或全场自伴性。

本模型明确选择静质量钟 Hamiltonian h(ν)=√(1+ν²)，不是标签ν。归一两点二次型
\[
C_p(t)=Z^{-1}\int_{s_*}^\infty da\,db\,Q(a,b;p)e^{-it[h(a)-h(b)]}.
\]
态在 A 的域内。绝对收敛、实且偶；推前到 Ω=h(a)−h(b) 的正谱测度说明正定性，不说明任意t点值非负。半径恢复后频率 Ω/ℓ，典范标量两点单位ℓ⁻¹。平方盒 a,b≤L 的未归一外尾界为
\[
\frac{2e^{-\beta L}}{\pi^2}
\left[\frac{L^3-s_*^2L}{\beta}+\frac{3L^2-s_*^2}{\beta^2}
+\frac{6L}{\beta^3}+\frac6{\beta^4}\right],\quad \beta=2\pi.
\]
除 Z 得归一界，对所有实t成立；盒内求积误差另算。

## 6. 旧 recoil 两点的有限质量补全

记 M=(a+b)/2，ω=a−b，A₀=cosh2πM，B₀=coshπω，C₀=coshπp。D=(A₀+C₀)(B₀+C₀)/4。旧定义
\[
k_p(\omega)=\frac1{8\pi(B_0+C_0)},\quad
R_{\rm old}=(1-\omega^2/4M^2)\frac{A_0-B_0}{A_0+C_0},\quad
Q_{\rm old}=M^2e^{-2\pi M}k_pR_{\rm old}.
\]
本轮直接代数给
\[
\frac Q{\rho_G(M)\sech^2\pi M}=k_pR_{\rm phys},\qquad
R_{\rm phys}=R_{\rm old}\frac{A_0+1}{A_0+B_0},\quad0<R_{\rm phys}\le1.
\]
用大质量匹配统一整体规范，
\[
\frac{\pi}{8}Q=Q_{\rm old}n_{\rm loc},\quad
n_{\rm loc}(a,b)=\frac1{(1+e^{-2\pi a})(1+e^{-2\pi b})}.
\]
这用新增局部插入、态和物理内积输入确定了此模型，未解决有限κ的真实插入/缝合。旧 nλ 自由族在仅有领先匹配时仍成立；旧 n(M,M)=1 是诊断选择而非物理定律。

固定ω,p且 M→∞，R_phys→1。真实能差另满足
\[
|h(a)-h(b)-(a-b)|\le|\omega|/[2\min(a,b)^2].
\]
证明是积分 1−h′(x)≤1/(2x²)。紧时间范围相位误差≤|t|乘此界；没有据此声称任意增长时间或全实轴相对收敛。

## 7. 数值、独立性及局限

[主诊断](../results/stage4/ds3-group-pairing/diagnostics.json)81项：球函数9，双观察者包2，三体群积分9，高角gluing15，Q/对称/旧核比9，Z3，创建率界3，平均谱/能差15，时间求积16。[非球面诊断](../results/stage4/ds3-group-pairing/offdiagonal.json)16项。全部通过，不替代解析证明。

主三标签样本 (.7,1.1,.6)、(1.4,.9,1.2)、(3.2,.4,.5)。角三元组 (0,0,0)、(1,1,0)、(2,1,1)、(2,2,2)、(3,2,1)。主要50位。三体直接积分截到65，尾界2 atanh(e⁻⁶⁵)/(abp)。双体测试 g=(ν−.4)⁴(1.6−ν)⁴、质量系数νg，径向R=40/80；Schur目标(π/2)(1.2)¹⁷B(9,9)，四阶导数Parseval给R⁻⁸尾。

另以旋转singlet S_l=(2l+1)^(-1/2)Σ_m(-1)^(l−m)|lm⟩|l,−m⟩|00⟩，I(S_l)=√(2l+1)T_ll0/T_000。直接群积分
\[
\frac{(-1)^l}{\sqrt{2l+1}}\int_0^{70}\sinh^2r\,
D^a_{l0,00}(r)D^b_{l0,00}(r)\phi_p(r)dr
\]
独立比较 K₃I(S_l)。boost 从真实紧边界角积分及有限多项式原函数计算；前两组质量、l=1,2，45位，外尾≤2√(2l+1)e⁻ᵀ(T²+2T+2)/[p(1−e⁻²ᵀ)]。它检验非000分量而非重复归一定义。

Z(.2)=.01552699344879428，Z(.8)=.002500278866978929，Z(1.5)=.00009085555337748436。

| t | C，p=.6，s*=.2 | C，p=1.2，s*=.8 |
|---|---:|---:|
| 0 | .005995878861534646 | .001281989884583574 |
| .5 | .005942866399361680 | .001261299720256264 |
| 2 | .005231147771197946 | .001002967854031784 |
| 5 | .002966801878925909 | .000418357991895478 |

时间平方盒 L=s*+8；NumPy双精度80/120 Gauss阶。两组归一外尾分别≤5.11e−20、8.94e−21，盒内加阶变化更大，未区间认证。参数、依赖版本、源哈希及实际复值均在JSON保存。

## 8. 继续位置

下一步先核两物质/四点的完整CG自旋通道、真实规约范数和重耦合，并建立0→1→2作用域。独立自由物质反极断开基准为 pp′/[16π²sinhπp sinhπp′]，须在其大质量条件下检验；同场需额外Wick项。复合gluing、完整Fock代数、有限κ缝合、有限G及完整全息均开放。本轮没有新代理、外发或正式发表。
