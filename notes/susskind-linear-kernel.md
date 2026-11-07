# 三阶一次残量的完整q核：协方差桥接

2026-11-06，本地。[七页PDF](../output/pdf/de-sitter-susskind-linear-kernel.pdf)、[LaTeX](../manuscript/susskind-linear-kernel/main.tex)、[来源](../literature/susskind-linear-kernel-sources.md)。无Git提交、推送或ZIP。

## 本轮结果

沿用完整overlap complex p-body Gaussian Hs，标准CAR、β=μchem=0、string时间、ℏ1。F=W₁T₃−tbar(Qf)H，W₁为普通coupling chaos。先DS λ=p²/N固定正，再小λ；q=exp(−4λ)，σ²=Js²/λ，V=〈F²〉DS=64Js⁶q²。

本轮证明：h=H/σ与y=F/√V的任意固定阶联合矩，等于两种正交q-Gaussian变量的配对矩；只有同字母pair保留，每次弦交叉权q。不是有限N统计独立。完整两点为已有Mehler核r=q，固定时间后小λ趋sech²(2Js t)，频率谱S1=πω/[4Js²sinh(πω/(4Js))]，S1(0)=1/Js。

原“F是否只有正动态矩而完整谱未知”的任务因此推进至明确正测度及唯一全核。有限q测度绝对连续，没有零频atom，有限N归一exactly static权重也趋0。该结论仍不把d1变成体质量。

## 为什么同J算符可以接到独立probe文献

用纯rank单位方差covariance Q_a，已有unit variance r-body ensemble K_r=Σω_a^(r)Q_a。ω公式沿用Kota和前两轮。设c=t/C_m，δ=c−Σωc，W=Σωδ²，u=δ/√W；标准化h_N=H/√C_m，y_N=F/√(C_m³W)。

精确HH/HF/FF三个covariance rank权重为ω、ωu、ωu²。HF通常非零，不能从〈HF〉0推出有限N独立。

固定j，s=p−a，精确比

R_j(s)=ω_(p−s)^(p−j)/ω_(p−s)^p
=s下降j × (N−2p+s+1)上升j × (N−m+p−j+1)上升j
/[(m−p+1)上升j × (p下降j)²]。

半填充R_j→(s)下降j/λ^j，有统一多项式界。前轮已证明s→Poisson(λ)、u_(p−s)→(s−λ)/√λ，并控制factorial tails。因此，在rank权重总变差范数中：

K_hh=K_p；
K_hy=√λ(K_(p−1)−K_p)+o(1)；
K_yy=λK_(p−2)+(1−2λ)K_(p−1)+λK_p+o(1)。

λ>.5时系数可以负，这只是Wick图中的有符号协方差替换，不构造虚假的正Gaussian surrogate。

## 整图误差与全阶相消

固定n弦图Dπ，给每根弦一个独立、单位方差的纯rank Hermitian Gaussian A_j。每个变量出现两次；ensemble×normalized trace的Hölder及EtrA_j^(2n)≤(2n−1)!!给|Dπ|≤(2n−1)!!。对有符号rank协方差展开，界乘上各边总变差范数。

若最大替换误差ε_N，B=max(1,2√λ,2λ+|1−2λ|)，整张n弦图误差≤n(2n−1)!! B^(n−1) ε_N→0。这一步使单边协方差逼近能用于任意固定阶，而不是把单crossing误当多弦因子化。

替换后每根弦是已有独立embedded ensemble r=p−j，j固定。应用文献的generic-body弦规则：交叉exp[−r_i r_j/(Nρ(1−ρ))]，ρ=m/N；半填充所有固定body shift均趋q。因此每种body-shift选择给相同q^cr，系数逐边求和：HH1，HF√λ(1−1)=0，FFλ+(1−2λ)+λ=1。任何含HF边的图消失，其余图给标准同色q配对。

full trace用|m−N/2|≤N^(2/3)一致中心窗，外面binomial电荷指数尾压过有限阶Gaussian矩及t_a多项式界。阶数固定，未把随N增长的word degree放进结论。

## 核与更高矩

〈FH^uFH^v〉DS=Vσ^(u+v)P_uv(q,q)。总H次数4时P40=2+q，P31=q(2+q)，P22=1+q²(1+q)，故μ4=2Js⁴(1−q)²(3q+5)/λ²→256Js⁴。另外〈y⁴〉=2+q、〈hyhy〉q。

已有ρ_q(θ)=(q,e^{±2iθ};q)∞/(2π)、Eθ=2σcosθ/√(1−q)，K_r=(r²;q)∞/∏_(ε,ε′)(r exp[i(εθ+ε′φ)];q)∞。取r=q，ρρK提供正两能量测度。所有混合矩收敛、偶高矩一致可积及compact候选矩确定性给弱谱收敛，继而固定实时间C_F。没有把有限N的Taylor展开直接交换极限。

数值使用C_F=Σq^n|〈n|exp(−itT_q)|0〉|²及(Tq)_(n−1,n)=σ√[(1−q^n)/(1−q)]。小λ取d1，沿用已证明任意d固定后代收缩和已有probe小λ公式，得到sech²；不是新发明Mehler/Fourier公式。

## 多项式rank滤波的推论

任意固定多项式P_i(s)，Y_i=ΣP_i(p−a)H_a/√C_m。边权含P_iP_j，展开有限个下降阶乘，用更多固定body shifts重复同一论证，得到q-Gaussian联合矩，Gram Gij=E_Poisson P_iP_j。正交Charlier多项式给一族正交变量；H对应P0=1，F归一极限对应P1=(s−λ)/√λ。此推论是微观rank滤波的算符构造，哪些属于局部导数singlet塔、能否产生不同物理species，需要另算。

## 实际核对

- cubic_covariance_bridge.py：2286个精确有理shift ratio；λ=.1,.25,.5,1，p16/32/64/128/256，20组全rank误差，80/120位保存65位比较。λ1包括负混合系数。
- cubic_sixth_wick.py：三个有限系统各五种six-slot word，每种15配对；保存同J混合pair总和，并用四slot重放已知B20/B11。N5p2m3,w1时V180，M4=314280，归一1746，未声称等于DS。整数einsum有绝对求和上界，所有分母用Fraction。
- check_linear_joint.py：九组全rank〈(h+zy)^4〉多项式，保留奇次混合系数；与(2+q)(1+z²)²的有限差保存，80/120位比较。
- run_linear_kernel.py：五组λ，Jacobi n33/52/138/302/574及+16；解析截断<1e−13，Mehler128/192双实现、四个Fourier点50/80位；λ.005在t∈[0,1]相对小λ核仍差约.560%。没有把数值精度当极限误差。

首次数值绘图因mathcal J的轴标签不合Matplotlib语法退出，错误记录保存；改成mathcal{J}后正常完成，未改变科学容差。工作进程均终态。

下一实质任务：求F与第三chaos的R/S/Z实时cross，完成raw M3谱；本轮(H,F)联合定理不自动覆盖R/S/Z。之后推进X/R/F混合高点、相互作用及体平移字典。完整de Sitter对偶目标保持active。
