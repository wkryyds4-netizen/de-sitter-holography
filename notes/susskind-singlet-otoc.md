# 导数 singlet 的乱序关联与 scrambling

2026-11-06。沿同一 Susskind complex DSSYK 主线，已把上轮多点归约接到现成四点核，实际计算完整同种算符 OTOC 和对易子平方。[六页 LaTeX](../manuscript/susskind-singlet-otoc/main.tex)、[阅读 PDF](../output/pdf/de-sitter-susskind-singlet-otoc.pdf)。本轮仍仅本地保存，不 Git 提交、推送或打包。

## 现在可比较的动力学

模型、态、算符保持不变：标准 CAR，完整 Hermitian Gaussian G，指标集可重叠，β=μ=0全归一迹，H=H_s生成string time，ℏ=1，w=4ᵖ(J_s²/λ)(p!/Nᵖ)²，先N,p→∞固定λ=p²/N>0。q=e⁻⁴λ，σ²=J_s²/λ，X=(O−cB)/√W_d，c=2λ(1−q)/(1+q)，W_d=4J_s⁴[1/λ+4q/(1+q)]。O=A−EA(Q)，A=Σ[H,cᵢ]†[H,cᵢ]；B=H²−EH²(Q)守恒。

上轮已经证明(H/σ,X)的所有固定阶联合矩及固定实时关联对应HH:q、HX:r=q²、XX:χ=q⁴的已有多弦规则。本轮直接应用BINT/BNR既有公式；没有重新解R-matrix。

真正的乱序量与对易子平方为

\[
F_{\rm OTO}(t)=\tau(X(t)XX(t)X),\qquad
\mathcal C(t)=\tau([X(t),X]^\dagger[X(t),X])
=2\{1+(1+\chi)G_{r^2}(t)-F_{\rm OTO}(t)\}.
\]

τ是annealed tracial state。F为实数但不保证正；C≥0且C(0)=0。直接配对求得

\[
\mathcal C(t)=2\sigma^2(1-r)^2(5+2r+3r^2)t^2+O(t^4).
\]

t²系数的小λ首项为1280J_s²λ。这里是β=0导数singlet的微观不对易动力学，不是凭先前相邻插入的κ₄推测OTOC。

## 三个channel的正确组装

令J为对称Jacobi矩阵，J_{n,n+1}=σ√[(1−q^{n+1})/(1−q)]，U±=e^{±itJ}，D_r=diag(rⁿ)，u=U−|0⟩，v=U+D_ru。

两个uncrossed channel在本时间排列下相同，均为 **U(t)=vᵀv**；这里是转置，不是共轭转置，因而该项可负。crossed channel（暂不含互探针χ）可写为

\[
f_k=(q;q)_k,\quad b_i=\frac{(r^2;q)_i}{(q;q)_i}r^i,\qquad
V_{li}=\sqrt{b_i}\,r^l\sqrt{f_{l+i}/f_l}\,u_{l+i},
\]
\[
\mathcal R(t)=\|D_r^{1/2}U_+V\|_{\rm HS}^2\ge0,
\qquad F_{\rm OTO}=2\mathcal U+\chi\mathcal R.
\]

这是BINT(3.9)/(4.10)作对角相似变换后的正项重排，原文§4.2省略的互探针因子在这里恢复χ=q⁴。t=0得到U=R=1、F=2+χ，正确重现⟨X⁴⟩。r=q²时比值(q⁴;q)_i/(q;q)_i=(1−q^{i+1})(1−q^{i+2})(1−q^{i+3})/[(1−q)(1−q²)(1−q³)]，避免了大阶乘相似矩阵。

## 数值结果及真正的误差范围

五组λ=.25,.1,.05,.02,.01，每组31个J_s t∈[0,1.5]时间点，Jacobi尺寸53、103、190、459、921。J_s t=.5处C分别为2.47621、2.08900、1.79281、1.26429、.837141；随着λ降低，增长推迟。λ=.25的采样F有约−.00365小负值，由uncrossed项贡献；它不违反C的正性。

独立核对：

- 全部四弧Hamiltonian次数和0/2/4/6、q=2/5与3/4：260组、780项有理channel比较，直接弦图与原非对称转移和一致。
- 从2τ(H²X⁴)+6τ(HX²HX²)−8τ(HXHX³)独立推出短时t²系数；五组 t=.0001/Js 数值与之核对。
- 4组50位、同有限L=18的原非对称矩阵指数/指标求和，与对称正项式比较，差≤1.6e−15。它检验相似变换和组装，不冒充无限截断验证。
- 8个参数/时间点的独立Mehler三能量积分核对uncrossed，96/192节点，差≤5.0e−14；完整相消条件数与虚部保存。
- 25个L+16主曲线复核点，全F最大差≤1.1e−14，相对差约6.5e−14；过零附近只用绝对误差判断。

解析尾界单独处理非真空传播，不能沿用仅vacuum Taylor尾。设R_H=2σ/√(1−q)、a=|t|R_H、B_d=aᵈ/d!/[1−a/(d+1)]，源截断0≤l<M<L。ε₀=2B_L、ε_M=2B_{L−M}，C_q=1/[(1−q)(1−q²)(1−q³)]，K=C_q/(1−r)。由||V||HS≤√K、源尾≤√K rᴹ，取δ=ε₀+ε_M+2rᴹ，得到

\[
|F-F_L|\le(4+2\chi K)\delta.
\]

全部主曲线区间上实际界3.6–4.9e−11；它同时控制Jacobi无限维和i/l弦索引尾。C的界为2[2(1+χ)ε₀+(4+2χK)δ]。不含浮点舍入，后者由上述独立比较经验核查。没有从1.5后的未知曲线宣称普适晚时plateau。

另有独立70位尾界重放，核全部8个主曲线/移动时间截断；164条组装算术重放。155个主点的uncrossed内部相位相消条件数最大约411，而三完整channel之间的相消条件数最大1.84。数据保存在error-replay.json，避免只看最终F掩盖内部相消。

## 接到已有的scrambling极限

精读Mukhametzhanov2303.03474v3 §5.1–5.2，下载原PDF并查看印刷页12（PDF第14页）核式(5.20)/(5.22)的前因子。原文J=1；字典是λ_src=4λ、J_src=2J_s、u_src=2J_s t、Δ=2。由原saddle(3.11)先取β→0，λL_src→2，于是物理string增长率为4J_s，且

\[
z=\frac\lambda2 e^{4J_s t},\qquad
\mathcal R(t)\longrightarrow\mathcal F(z)=z^{-4}U(4,1,1/z),
\quad\lambda\to0,\ z>0\ \text{fixed}.
\]

采用标准归一Tricomi，即F(z)=∫₀∞y³e⁻ʸ(1+zy)⁻⁴dy/Γ(4)，与文献(5.24)的F(0)=1相符；三点用60/50位hyperu与gamma积分交叉检验。

在此已有scrambling近似的领先阶，χ→1，uncrossed贡献消失，完整F→F(z)，C→2[1−F(z)]。z≪1时F=1−16z+O(z²)，故late string time与small z重叠区有C≈16λe^{4J_s t}、t*=[log(1/λ)]/(4J_s)+O(1/Js)。它不替换t→0的t²起始，也不是低温Schwarzian式直接代β=0。

实际用λ=.02,.01,.005、z=.025,.05,.1测试移动时间。以z=.1为例，R=.411526,.358189,.333041，对应极限.309611。有限λ偏差明显，大于数值尾；不能说已经等于极限。最大L1559并加16对照。完整同种F的finite-λ差有通道间抵消，不用它的偶然更接近来替代crossed本身的偏差。原值与各channel均保存。

## 接下来接哪一段物理

我们已有一个具体singlet的传播核与完整交替时间四点。下一步优先明确能与之比较的Rindler体波包四点：态、径向/空间涂抹、string与boost时间，以及哪种介子叠加对应X。微观增长率尚不能直接认作某个体观察者的proper-time Lyapunov率。另求固定时间的首阶λ修正，包括前轮κ₄。完整对偶、finite curvature和在壳体散射仍开放。
