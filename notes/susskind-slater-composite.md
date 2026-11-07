# 同一 J 的曲率复合场与 Slater selected 代数

2026-11-08。上一轮 polynomial-density 给任意固定 Gaussian degree 的制备态一致轨道密度。本轮构造实际 J₂ 与 h,u,M,M† 的全部固定混合矩，并直接补 central-sector 的高矩/精确 pinching 接口。正式证明：manuscript/susskind-slater-composite/main.tex。

## 定义须先于限制扇区

standard CAR、overlapping Hermitian complex Gaussian p-body、同原 J、m=N/2固定Slater O、annealed φ_N，λ_N=p²/N→固定λ>0。D X=QX−Σc†Xc，L X=(NX−Σ(c†Xc+cXc†))/2。C_Q=E H²，C_m为m扇区值，Z₂=H²−C_Q。

κ_N=2λ_N(1−q_N)/(1+q_N)，q_N=e^-4λ_N；
j_N=C_m^-1 (2p−D−κ_N)Z₂|_m，v_N=4λ_N+16λ_N²q_N/(1+q_N)，x_N=j_N/√v_N；
h_N=H_m/√C_m，W_N=(p−L−λ_N)H|_m/√C_m=√λ_N u_N。

所有loss先在fullFock作用，再限制central；C_m是固定目标扇区的确定常数，置于loss之后。不能先把H换成H/√C_Q后再做D，因为D连接相邻charge，结果会改变。M†仍为同J的纯occupied→empty块/√L_N。无量纲string能量、hbar1；真实H时间需乘√C_m。

## selected 多源规则

j插入展开为相邻两个H槽的whole D mark，禁内部Gaussian配对（Z₂的纯chaos）；W为单H槽L mark；其他H为ordinary。Gaussian edge右端C_P†C_Q、左端其adjoint。G不限标签；forward F要求P⊂empty,Q⊂occupied；backward B反之。nilpotent z_a各只出现一次。

每site先乘mark块完整2×2 CAR词，再作用1−zD_i或1−zL_i。D_i[[a,b],[c,d]]=[[0,b/2],[c/2,d−a]]，L_i=[[a−d,b],[c,d−a]]/2。保持字段原次序后取occupied o=0或1的对角元，再乘(−1)^Σ_cross t_e t_f；不合selected标签则置零。得到f_o(S;z)。每个P/Q标签的正确singleton分别在o=0/1，s_i=f_correct({i};z)零阶1。无mark内部配对使每个singleton的loss只有odd缩放1−z/2。因此exp_nil[pΣd_a z_a]Π_i s_i^p=1，完整多源的extensive项精确抵消。

collision归一F_o(S;z)=f_o(S;z)/Π_(i∈S)s_i，分母可逆。令k=#Gaussianedges、k_G=#G、b=#H/selected混合edge，ν=#crossings+#inner-selected nestings，γ=e^-λ=q^1/4。定义

Kπ(z)=2Σ_(i<j,o)F_o({i,j};z)−2k(k−1)−2k_G，Kπ(0)=−4ν。

任意固定word极限=Σ_forward π γ^b [Πz_a] exp{λKπ−λΣ_W z_a−κΣ_j z_a}，κ=2λ(1−q)/(1+q)。使用u,x时再除λ^(#W/2)v^(#j/2)。每张图保留全部source masks，不用旧tracial mixed-q简单补独立物种。

## 一般收敛证明

精确finite coefficient来自f_0^m f_1^m。令multisite计数n_(o,S)、标签defect d_i、占用site U_o；除L_N^k并取消singletons后，每项是(m)_U0(m)_U1/(m)_p^(2k) ×Π_i(p)_di ×ΠF_o(S)^n/n!。沿旧Slater falling-factorial估计，所有F_o的系数norm≤固定B，将activity m(p/m)^|S|乘B，仍有一致指数级数尾。|S|≥3总activity→0；先支配再去高重碰撞和defect尾。backward强迫某正确singleton数为零，即d_i=p，落在同一大defect尾中。source保持odd singleton parity，不会救活它。

固定defect的falling background=e^-2λk(k−1)，每个pair-membership活动度p²/m→2λ。H归一用D_N=C(m,p)C(m+p,p)，D_N/L_N→e^2λ，n_H=2k_G+b，另给e^-λnH；合并即上式。任意固定word及P†P的矩存在且一致有界，有限N正性传到共同GNS。

## central-sector 与 exact pinching

全部G时旧local trace T=(f0+f1)/2，旧singleton t_i=s_i/2，所以旧R_ij=T/(t_i t_j)=2Σ_o F_o。旧background−C(2k,2)−k=−2k²恰等于本式全部G的background。整个多源指数一般性相同，而不只是若干数值矩吻合。covariant word的annealed矩在wedge^m上由Schur精确等于固定Slater矩。由旧端点交换相消，得到central (h,u,x) 的既有tracial mixed-q律，Q=[[q,q,q²],[q,q,q²],[q²,q²,q⁴]]。所有固定central高矩因此成立，未从fullcharge平均反推central。

x的正两能量测度ν_N弱收敛到K_(q²)(E,E′) μ_q(dE)μ_q(dE′)。全mixed矩+高矩给tightness/UI，compact极限确定测度；q²<1、μ_q无原子，因此没有对角质量。Portmanteau给||P_H x_N||2,m→0，pinching contraction及更高矩插值升全部fixed S_r。完整简并block保留，不交换长时极限。

实际A0=(1−P_H)(2p−D)Z₂|m；回归项κ_N Z₂|m与H对易，故归一实际a_N=A0/(C_m√v_N)=(1−P_H)x_N精确。原固定κ的J₂与C_m j_N只差趋零系数乘Z₂和κ(C−C_m)标量，fixedS_r消失。此处a_N是observable，不是旧prep norm。

新任意fixed P(h,u,x,M,M†) 的Gaussian degree仍固定，且a_P=E||PO||²由新source theorem有界。因此可调用上轮density界：E||X P O||²≤C_P||X||4,m²。实际a_N−x_N和b_N−u_N均在这些prep向量上趋零。对单个covariant block W内任意固定a,b及h-timephases，trace Hölder/telescoping、density和极化给φ_N(P†W(h,a,b)Q)→φ(P†W(h,x,u)Q)。时间尾用centralGaussian h高矩和fixed插入高矩，先在trace L²控制再转移到selectedpreps。多块nonpolynomial pinching与M交错仍需新论证。

## 共同GNS上的协变代数

对core ψ=PΩ、任意selfadjoint covariant polynomial A，density界和central矩给||A^nψ||²≤C_P[τA^(4n)]^1/2≤C_P||A||_D^(2n)，C_P不依赖n。positive moment logconvexity给||Aψ||²≤||ψ||^(2−2/n)C_P^(1/n)||A||_D²，n→∞即协变作用受标准tracial C*norm控制。真空restriction=τ使表示faithful；|〈ψ,Aψ〉|≤C_P||A||2给L² density，densecore/polarization导出vN上的normal表示。

h,u范数≤2/√(1−q)，x≤2/√(1−q⁴)。没有据此分类整个selected代数，未证明M在新空间compact，未给j的闭式开弦transfer。旧小cyclic空间的compactness不能直接搬来。

## 两个同J信号与验证

φj²=v，φjh²=0；φj(M†)²=−2κγ²q，φjM†h=−2κγq，φjhM†=0；φjWM†=φjM†W=0；φMjM†=0，φMj²M†=v。

第一非零信号也可直接看：与双selected creation配对时所有存活hop不交叠，2p−D给零，只剩−κh²项；两图给φh²(M†)²=2γ²q。独立中心化probe会把这个重叠错误置零。实际归一curvature由新pinching桥给φ_N[a_N(M†)²]→−2κγ²q/√v、φ_N[M a_N² M†]→1。前者是off-diagonal矩阵元，别误称归一double-excited期望。

精确finiteCAR在N2p1/N4p1/N4p2的fullFock矩阵上直接施wholeD/L，与local source逐项核592704次；含40个负项，max整数16，forward/backward/unrestricted全部保留。326allGsource恢复旧指数，centralwords恢复旧helper；340无j词回到原HU/M transfer。十个symbolic式与两个14维rationalGram LDL核验；rationalγ,λ点不在γ=e^-λ物理曲线上，只作精确代数诊断，不推无限维新sector。科学JSON首次成功封存，不重跑。

下一步用source规则构造高效j/open-transfer与低模块，再接更高charged/neutral后代及完整bilocal。固定高阶wholeWick mark的source部分可按同法处理，但实际rank/charge-removal桥、selected全时间尾、nonpolynomial interlacing不能省略。fulltrace bilocal四点保持已完成。体空间/质量生成元与finiteG仍开放；本轮仅本地，不Git/issue。
