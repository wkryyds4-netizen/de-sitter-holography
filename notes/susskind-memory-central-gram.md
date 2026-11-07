# 中央四方向 Gram：用极化恒等式闭合交叉项

2026-11-06。LaTeX：`manuscript/susskind-memory-central-gram/main.tex`；六页PDF：`output/pdf/de-sitter-susskind-memory-central-gram.pdf`。同J、完整重叠complex Gaussian p-body系综、归一最大混合sector迹、普通Gaussian W2、hbar=1。偶N，p²/N=λ_N→λ>0，m=N/2+xp，x固定紧集，λ_N正紧集。

## 完成的结果

沿用B=W2H²、L_in=W2{h,D−bar u h}、S=Σ(a+b−r)W2Z_ab;r、O=W2[Rem(k²)+(2p−m)h²]。

q=exp(−4λ)，z=1+q，b=2λ(1−q)，l=2λ(1+q)，s=4λ²(1+q)+2λ(1−q)，a*=4λ+4λ²(1+q)。本轮证明S–O/C_m²→s，故完整Gram/C_m²为

    [z, 0, b, b]
    [0, l, 0, 0]
    [b, 0, s, s]
    [b, 0, s, a*]

这次不再以候选general charged Casimir式为前提。那个更强的有限N公式仍保留为待证命题，但它已不阻挡所需中央极限。

## 为什么改用极化

前一轮已建立相邻charge等距I_m、同一Gaussian输入标号和一步归一产品差O(1/N)。每个输入/输出block设U=W2Z_m/C_m、V=I_m(W2Z_(m−1)/C_(m−1))、E=V−U，t=R_m ell_r，h0=m−2p，则真实算符恒等式O_block/C_m=tV−h0U给

    ΣW<U,O/C_m> = 1/2 ΣWt (||V||²−||U||²)
                  +ΣW(t−h0)||U||² −1/2 ΣWt||E||².

取W=X/N，X=x_(a+b)−x_r=S[N+1−2(a+b)+S]。W是charge无关的标号权重；所有系数Hermitian，内积实。不是把相位取实以掩盖误差。

最后项以W≤Kp、t≤KN及Σ||E||²≤K/N²给O(1/p)。第一项拆t=h0+δ。δ≤K(1+u+v+S)部分由带权Cauchy和输出四阶/输入矩给O(1/N)。只剩h0(F_(m−1)−F_m)/2，其中F_m=ΣW||U||²是已知neutral第一矩。

## 核心估计：归一后的neutral矩只漂移O(1/(Np))

F_m=Σω_aω_b f_ab(m)，其中

    f_ab=[−2ab(1+χ_ab^(m))+x_a x_b/Cf(m)]/N,
    Cf(m)=m(N−m)(N+1)/N.

neutral terminating series具有|c_j|≤K0^j/j!，且

    c_j(m−1)/c_j(m)=m(n−j)/(n(m−j)), n=N−m+1.

其差比为j(n−m)/(n(m−j))，中央|n−m|≤Kp，因此固定输入charge差|f_ab(m−1)−f_ab(m)|≤Kp/N²。Cf的差也满足这一界。

对输入rank，a降1使每项乘(a−j)(N+2−a)/[a(N+2−a−j)]。a≥p/2时该比与1之差≤Kj/p；若a<p/2，u/p≥1/2直接用有界性。故|f_ab−f_pp|≤K(u+v)/p。

归一权重精确比

    ζ_a=ω_a^(m−1)/ω_a^m=(m−p)n²/[(n+p)(mn−x_a)],
    ζ_a−1=[(n+p)x_a−pn(N+1)]/[(n+p)(mn−x_a)].

a=p处分子p²(m−p)，改成p−u带来≤KuN²的变化，故|ζ_a−1|≤K(1+u)/N。两扇区权重总质量都是1，所以权重变化项里可以先减f_pp。然后用factorial input尾，得|F_(m−1)−F_m|≤K/(Np)。乘h0=O(N)后仍为O(1/p)。这一步的归一相消才是闭合关键。

中间项利用t−h0=(u+v+S)−2λ_N+O((1+u+v+S)/p)，得到

    [||S||²+<S,L_in>+(2bar u−2λ_N)<S,B>]/C_m²+O(1/p) → s.

前轮已证(B,L_in,S)Gram及bar u→λ_N，因此交叉极限成立。最后用||S−X/N||/C_m=O(1/p)、O/C_m有界，换回真实S。

## 物理接口：正投影残量

v=s−b²/z=2λ(1−q)+16λ²q/(1+q)，Δ=a*−s=2λ(1+q)。两者严格正，四Gram行列式zl vΔ>0。

S扣span(B,L_in,O)后为T=S−v/(v+Δ) O−bΔ/[z(v+Δ)]B，||T||²/C_m²→vΔ/(v+Δ)>0。实际有限N projector同极限，因为旧三Gram逆在紧参数集一致有界。

接上一轮完整W2core强展开，(1−Π_N)W2Ecal范数²/C_m⁴→64λ²q⁴vΔ/(v+Δ)。原central charge守恒修正4λx(1−q)²仍在raw core的B方向；这里是被指定投影去掉，不能说raw correction不存在。C_m→Js²/λ后，先DS再小λ，方差趋1024Js⁸。只为equal-time W2，不包含raw memory其余chaos/能量后代或有限时间完整谱。

## 实际验证与保存

`neutral-drift.json`：758精确weight ratio、四个有限有理矩、36全秩邻sector差（λ=.1/.25/.5，p16/32/64/128，x−2/0/2），80/120位未序列化前比较最大差3.554e−74。无rank/series截断。λ=.5,x2时放大drift为−.184286、−.0987564、−.0512180、−.0260982，有限差保留。数值扫描不替代uniform不等式证明，也不是区间认证。

`polarization-car.json`：N5p1m3、N6p1m3、N6p1m4，3217有序Gaussian坐标对、每系统6个a,b,r blocks；实际CAR force与极化式逐块核对，最大误差2.78e−17（float64，阈值1e−12）。独立neutral矩亦匹配exact rational。

`gram-assembly.json`：符号验证Gram行列式、投影系数、残量、core尺度及有序小λ极限。

首版SymPy结构比较把−4(−λ+u)与4(λ−u)视为不相等，保留源码/失败记录后改用代数差factor=0；科学恒等式未更改。首版LaTeX遗漏转义产生一个formfeed，原稿/log保留，修正后双遍成功、六页95dpi逐页检查。代码/JSON/稿/PDF/旧阶段/渲染hash由local-evidence记录。

下一步转中央energy insertions与raw U2谱，再全charge尾、finite-memory两插入余项及Susskind体介子空间/质量字典。本轮只本地保存，无commit/push/ZIP、不启代理；整体goal active。
