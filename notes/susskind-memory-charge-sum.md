# 全电荷求和与有限记忆的精确去静态定义

2026-11-06。六页稿 `manuscript/susskind-memory-charge-sum/main.tex`，PDF `output/pdf/de-sitter-susskind-memory-charge-sum.pdf`。标准 CAR、同一 overlapping complex Gaussian p-body ensemble，hbar=1、string units；偶 N，λ_N=p²/N→λ>0，λ_N在正紧集。这里完成的是 β=μ=0 的 annealed full normalized trace，不是有限温或有限 G 的结果。

## 全sector定义与统一界

全迹权重 b_N(m)=2^(−N)binom(N,m)。m<p 时 H_m、H_(m−1)、memory 全为零，三个字段及反项定义为零。对 m≥p，保留此前 central-raw 的有限sector C0/C1/C2 公式，用实际 a_m、c_m、bar t_m、bar z_m、e_Q。窗外把该公式用作明确定义，不声称每个远离中央的sector有相同的渐近式。

共同方差的更强上界是

    C_m/(Js²/λ_N)
      = product_j 4(m−j)(N−m+p−j)/N²
      ≤ product_j [1+(p−2j)/N]² ≤ exp(2p/N)。

每因子的完成平方差是 (2m−N−p)²，连续m最大点(N+p)/2。相邻比 C_(m−1)/C_m=(m−p)(N−m+p+1)/[m(N−m+1)]≤p+1。

回归分母没有active-sector奇点：||B_m||²=C_m²(1+χ_m)，χ_m=Σω_a κ_a/C_m ∈[0,1]。所以包括 m=p、m=N 在内，c_m=〈O,B〉/||B||²都直接定义。端点不使用带可去分母的rank递推式。非活跃sector只因C_m=0而单独置零。

归一Rem/m是unital、normalized-trace-preserving CP map，因此Schatten收缩。Gaussian integration by parts给Etr h^(2k)≤(2k−1)!! C_m^k。由原始有限多项式及显式t_a的binomial比界，得到所有sector的粗界

    ||Y_m||_(2k) ≤K_k N² M_N，
    ||L_m||_(2k) ≤K_k N² M_N^(3/2)，
    ||U2,m||_(2k) ≤K_k N⁴ M_N²，
    M_N=(Js²/λ_N)exp(2p/N)。

这些界仅用于远尾，不用于判断中央减背景涨落强度。K_A²/m仍保留；不把它的守恒贡献删掉。

## 增长窗与全迹

取 R_N=√(A0 log p)，A0固定。不能单凭此前固定x定理推出新窗口；本轮逐项审核了余项：

1. 六阶rank矩的递推长度固定，负幂相消是z的多项式恒等式。相消后rational remainder为N^(−1/2)乘有限次z多项式，z=x√λ=O(√logp)，仍消失。
2. coherent A范数在u,v,j≤ceil(logp)展开。此前p²,p的相消及有限项x无关都是多项式恒等式，余项至多polylog(p)/p。窗外factorial尾压过原p⁴粗界，统一weighted A范数成立。
3. gA取S≤p^(3/4)；gH先精确减p ell与2λx ell，再对输出取较小p^b。有限次数Taylor余项和六阶尾给强替换。F及cubic remainder同用factorial weighted bounds。
4. central marked trace误差至多p^(−1) exp(KR_N)(1+R_N)^D=p^(−1+o(1))。一般force桥接中的x/p也消失，直接接S_N=(O−+O+)/2，避免重复charge transport。

因此 sup_(|x|≤R_N)Σ||V_m||²≤K_A0，V=(Y,L,U2)。取A0=20/λ_min，binomial Chernoff尾给窗口外 ≤K N⁸p^(−40)=O(p^(−24))。固定L与R_N之间的环带≤K exp(−2λ_min L²)。于是

    lim_(L→∞) limsup_N Σ_(|x_m|>L) b_N(m) Σ_i||V_i,m||² =0。

结合固定中央窗谱收敛和spectral Cauchy–Schwarz，得到full (Y,L,U2)双能量谱等于已求中央正matrix。交叉测度的tail也有TV界；相位未被丢弃。

更进一步，x=(Q−N/2)/p的characteristic function为cos(t/(2p))^N，极限Normal(0,1/(4λ))。joint charge/two-energy matrix因子化为

    sqrt(2λ/π) exp(−2λx²) dx × dν_*(E,E')。

这只是指定两能量谱matrix的因子化，不等于完整算符代数的有限N独立性。旧有序小λ核随fulltrace转移，不能交换λ和N极限。

## 有限memory的新明确接口

每个sample、每个sector定义P_h X=Σ_e P_e X P_e（完整能量简并块pinching）。它是投向commutant的HS正交投影，也等于finite-dimensional无限中心时间平均。

    Fscr_N,m(s)=(1−P_h)[m Mhat_m(s)]。

此定义精确减所有静态部分，无需先猜任意s的确定性背景。其前三个Taylor系数是(1−P_h)(Y,L,U2)。full正谱极限在E=E'线上质量为零；闭集上界给||P_hV||²→0。因此前三个系数保留原谱matrix。此处静态atom消失是先finiteN定义pinching、再N极限的结论；没有把无穷时间与N极限随意交换。

下一问题具体化为固定实s,u,T的

    Etr Fscr_N(s)† exp(iTH) Fscr_N(u) exp(−iTH)。

先从R†exp(iK_perp s)R的两个memory插入、共同coupling和pinching计算或建立uniform bound；不能只用前三个导数代替有限s余项。这个samplewise谱投影本身也不等于已建立体局域字典。

## 证据

- `global-sector-checks.json`：6766个精确binomial/product/adjacent-ratio案例；N3p1、N4p2、N5p2的全部15个sector以直接Gaussian整数CAR核B variance、A mean、BO与regression residual，包括inactive/full端点。
- `growing-window-checks.json`：12组x=±ceil√(4logp)、λ=.25/.5、p32/64/128的全部rank、完整terminating sum，80/120位；另4个full positive binomial tail sums。保存65位一致不等于interval。λ.5,p128,x5的F方差相对极限差仍+.17866，未隐藏。
- 取λ_min=.25、A0=80，p128的精确positive tail约1.8073e−86，N⁸tail≈6.1500e−48；解析尾界另存。
- `charge-sum-assembly.json`核6个prior/current证书source hashes、21个负幂相消、A norm两项相消、平方分解及charge characteristic expansion；未重跑旧证明驱动。
- LaTeX初次缺mathscr宏，改用mathsf；第二次有长inline系统列表溢出，拆开后双遍无警告。失败日志保留，最终六页95dpi全部核版。

仅本地保存，无commit/push/ZIP；长期goal仍active。下一位置见NEXT。
