# 中央电荷窗口的双向力谱与二次core谱

2026-11-06。六页LaTeX/PDF分别在 `manuscript/susskind-memory-central-spectrum/main.tex` 和 `output/pdf/de-sitter-susskind-memory-central-spectrum.pdf`。同J、完整重叠complex Gaussian p-body、standard CAR、hbar=1；归一最大混合charge迹。偶N，λ_N=p²/N→λ>0，m=N/2+xp，x固定紧集、λ_N正紧集，C_m→C=Js²/λ。

## 结果与适用对象

A_−=Σ[H,c_i]†[H,c_i]，A_+=Σ[H,c_i][H,c_i]†，O±=W2 A±，B=W2h²。c=2λ(1−q)/(1+q)，Y±=O±−cB，q=e^(−4λ)。中央窗口的完整两能量矩阵为

    ν_(Y−,Y+) = [[w,η],[η,w]] μ_(q²),
    w=4C²λ+16C²λ²q/(1+q),
    η=−4C²λq+16C²λ²q/(1+q).

旧q² Mehler核本身来自已有弦图框架；本轮证明其在此中央charge态下仍适用。不是由等时Gram猜测动态谱。

Π_N为span(B,L_in,O−)实际有限N投影，T_N=(1−Π_N)S_N。v=2λ(1−q)+16λ²q/(1+q)，Δ=2λ(1+q)，ρ=v/(v+Δ)，则

    T_N=(1/2−ρ)Y−+Y+/2+oL²,
    ν_T=C²vΔ/(v+Δ) μ_(q²), ν_(T,Y−)=ν_(T,Q4)=0.

实际投影W2core R2,N=(1−Π_N)W2 Ecal=8λ_N C_m q_N² T_N+oL²，因此

    ν_R2=V2 μ_(q²), V2=64λ²C⁴q⁴vΔ/(v+Δ),
    C_R2(t)=V2 G_(q²)(t).

固定λ的极限谱无静态atom。极限谱的归一二阶频率矩是2C(1−q²)。先DS再λ→0，有V2→1024Js⁸，G_(q²)(t)→sech⁴(2Js t)。这是投影后的W2core，尚未称为完整raw U2；其余chaos及动态后代需要继续恢复。

## 中央迹转换：定量O(1/p)

一个局域occupation项规定L个site，其中t占据，s=L−t空，d=2t−L，则sector/fullFock精确比

    R=2^L (m)_t (N−m)_s/(N)_L = 2^L binom(N−L,m−t)/binom(N,m).

取ε=(m−N/2)/N，L≤kp、k固定。logR的首阶falling-factorial近似A在稿第2页；显式绝对余项是E(m,t)+E(N−m,s)+E(N,L)，其中E(a,j)=j(j−1)(2j−1)/[12a²(1−(j−1)/a)]。另有全局logR≤L log(1+2|ε|)+L(L−1)/[2(N−L+1)]≤K。

展开得

    logR = L/2 log(1−4ε²)+d/2 log[(1+2ε)/(1−2ε)]
           +(L−d²)/(2N)+O(|ε|L²/N+L³/N²).

故|logR|≤K[(1+|d|)/p+d²/N]。利用logR的正部分有界，|R−1|≤e^K |logR|，无需R的下界。

每个Gaussian chord的两个p-set有相反singleton occupation且原大小相等；除去重叠site后失衡|d|≤c_k K_ov。正常/反向mark以及有限次数body-deficit mark只作用重叠site，Wick块禁止internal contraction；局域occupation展开的绝对系数≤poly(K_ov)e^(cK_ov)。hypergeometric factorial矩与Hölder给所有固定指数的K_ov矩有界。于是每个固定marked graph的中央迹与fullFock迹差≤K/p（采用Gaussian coupling和singleton归一）。保留全部符号，没有对相消后总值强行取绝对值再冒充相对误差。

同样应用于H²，得到C_m/C_Fock=1+O(1/p)，不用先做全charge积分。因而任意固定阶外部H插入、有限mark数、bounded scalar coefficients的Wick-block组合可复用旧fullFock结果。不能将O(1/p)误差无推导乘上raw counterterm中的p；也不覆盖增长word length、增长x。

## 一般charge的精确force桥接

D X=QX−Σc†Xc，Db X=(N−Q)X−ΣcXc†。复合必须先fullFock后restrict。精确

    C_x X=(N−m+1)DX+m DbX−D DbX,
    (Db h)_a=[x_a−p(N−m)]h_a/(m−p+1).

设u=p−a，E_in=Σu(u−1+N−2m)h_a/(m−p+1)，则(p−Db)h=2D_in h+E_in。中央E_in仅O(1/p)，不能继续用half的O(1/N)。

barA=(2p−Db)H²=A+ +{H,(p−Db)H}；E_double=(2p−D)barA，U_in2={h,D_in²h}+2(D_in h)²。Sproxy用Casimir X/N代替output deficit，精确

    N[Sproxy−(O−+O+)/2]
      =(N/2−m+1−2p)O−+(m−N/2−2p)O+
        +(2m−N−1)W2{h,D_in h}−W2U_in2+W2E_double
        +(m−2p)W2{h,E_in}.

双mark仍只作用重叠site，系数poly(K_ov)可积。输入factorial尾和Gaussian高矩控制其余乘积。除N后每项O(1/p)，再用已有S−Sproxy=O_L²(1/p)，得到S=(O−+O+)/2+O_L²(1/p)，中央紧窗一致。

## 从矩到谱

对marked fields使用所有固定双能量矩、正matrix measures、tightness与compact候选谱的moment uniqueness，得到中央uniform weak limit。旧Q4 reverse/forward endpoint certificates随迹转换也成立，故Q4–Y±全谱cross0。S/T/R2的实际强替换用谱Cauchy–Schwarz给总变差控制，足以传递bounded spectral tests和固定实时间；对被替换字段的无界能量多项式，还须单独有uniform integrability，不能只引用L²替换。

## 执行证据

- `trace-ratio.json`：24 exact binomial checks，324中央样本（λ .1/.25/.5；p16/32/64/128；x−2/0/2；L/p2/4/8；d0/2/6）。每个点满足显式log余项与上界，90/130位变化最大2.225e−85。λ.5,x2,L4p,d2的R为.7985081/.9027791/.9522632/.9763479，有限偏差保留。
- `force-bridge.json`：六个非half系统、每组三个整数Hermitian输入，108 exact cleared-rational CAR/rank identities。直接核对一般E_in与fullFock复合。
- `spectrum-assembly.json`：符号matrix/projection/V2，核查旧5份科学输出与各自源码hash。复用128 reverse templates/6072words和原有任意固定energy阶端点证明，未重跑旧扫描。

实际双遍编译成功，最终六页95dpi全部查看，无溢出/未定义引用。三JSON、稿、PDF、渲染与前阶段hash记录在local-evidence。只本地，无commit/push/ZIP，不启代理。

下一步：中央raw U2的low-chaos descendant table与精确守恒counterterm，再filtered F的中央lift/cross，组装完整(Y,L,U2)谱。随后全charge窗外可积性、finite-memory两插入余项、体介子空间/质量字典。整体goal active。
