# 完整双时 singlet 的压缩演化与记忆核

2026-11-06。六页 LaTeX/PDF：manuscript/susskind-bilocal-channel/main.tex；output/pdf/de-sitter-susskind-bilocal-channel.pdf。主线继续微观bilocal到体介子的构建，本轮不是完成全息字典。

## 完整算符的表示

标准CAR、旧complex p-body Gaussian H_s，时间为string时间。固定电荷m>0，h=H_m、H_−=H_(m−1)。定义 Vψ=m^(−1/2)Σ_i|i〉c_iψ，K=1_N⊗H_−，则V†V=1。令W(t)=e^(iKt)V e^(−iht)，完整B_ab=S_m(t_a,t_b)/m=W_a†W_b。

因此B_aa=1，B_ab†=B_ba，||B_ab||≤1，所有多时间块矩阵[B_ab]半正定。共同中间时间的精确组合为

`B_ab B_bd = B_ad − W_a†(1−W_b W_b†)W_d`。

一般四时间乘积为W_a†W_bW_c†W_d，不需要忽略上一轮四项反对易子结构。能量表示是

`S_c(δ)_ab=Σ_iα (c^i_αa)* c^i_αb exp{iδ[E_α−(E_a+E_b)/2]}`。

先相干相加中间态再求平方，才能得到正确的Gram协方差。这个内部频率是相邻电荷能量组合，尚不是介子质量或光前内部动量分数。

## 实际构建的有效方程

P=VV†、Q⊥=1−P、R=Q⊥KV、K⊥=Q⊥KQ⊥，h_eff=(1−p/m)h。对F(t)=V†e^(iKt)V，有

`F′(t)=i h_eff F(t)−∫₀ᵗ M(t−s)F(s)ds`，

`M(s)=R†e^(iK⊥s)R`。

这是已有projection-memory方法在本bilocal上的应用，有限矩阵下直接分块积分证明；也由Schur补给完整resolvent。不是新发明Mori–Zwanzig，也没有假定无限维极限的域。

其初值有精确正形式

`G=M(0)=A_m/m−p²h²/m²=R†R≥0`，A=Σ[H,c_i]†[H,c_i]。

m=p时记忆严格为零。p=1也可能有非零G，故有记忆本身不是介子相互作用证据。

完整中心bilocal展开为S_c(δ)=m−ipδh−δ²A_m/2+O(δ³)。偶δ部分的首个非守恒系数就是A；两个中心时间偶bilocal对易子在δ²ε²处为[A(T),A(T′)]/4。此Taylor语句只在有限N，不与DS换序。

令旧a_m=Etr_m A、C_m=Etr_m h²、c_m为旧H²投影系数，Y_m=A_m−a_m−c_m(h²−C_m)。则G的常数及H²守恒项扣除后严格为Y_m/m。这把旧全迹DS两点Wd G_(q²)(T)接到记忆初值的中心时间涨落；它不等于尚未求出的正交演化M(s)。这两个时间不能混同。

## 有限时间的受控DS极限

直接复用固定电荷协方差C_m=w binom(m,p)binom(N−m+p,p)，w=4^p Js²N/p²(p!/N^p)²。相邻比为

`C_(m−1)/C_m=(m−p)/m · (N−m+p+1)/(N−m+1)`。

结合a_m=(2p−m)C_m+mC_(m−1)，得

`Etr_m G=C_m p(m−p)(N+1)/[m²(N−m+1)]`。

半填充固定λ=p²/N>0时C_m→Js²/λ，p Etr_mG→4Js²。由KV−Vh_eff=R的Duhamel式得到所有有限系统、所有实t的界

`E||B_0t−1||²_(2,m) ≤ 2t²[Etr_mG+(p/m)²C_m]`。

故在每个固定紧时间区间，普通B趋于单位算符（annealed HS范数），误差平方上界O(t²/p)。这是实际有限时间界；不需要上一轮先t→0的顺序，不控制t~√p及重新放大的涨落。

对完整β=μ=0迹，用|m−N/2|≤N^(2/3)中央窗内的一致乘积展开，窗外binomial尾≤2exp(−2N^(1/3))；C_m≤(Js²/λ_N)exp(2λ_N)、Etr_mG≤(N+1)C_m控制尾。得到同样O(t²/p)。再由Q数目涨落，E||2S/N−1||²≤8E||B−1||²+2/N→0。这个结论要求重建保留适当减背景/归一后的涨落，并不证明有限δ的完整涨落场已构造完毕。

## 验证与继续

- checks.json：四个确定Hamiltonian的20个整数等式、五个精确单项系综求和、18组80位参数；核正性、全Fock独立演化、记忆卷积24/48点、Schur补及12个有限时间界。
- projection-checks.json：两交互系统的完整四项不等时对易子、相干谱重建、Gram正性与偶δ差商二次收敛。N4p2m3四项范数平方约10.16、7.69、.995、9.62，总对易子3.264，明确有相消；不外推其非零到DS。
- 两个计算脚本首次运行通过。LaTeX首次缺少tr宏，失败源和日志保存；补宏、核正参考作者名后双遍通过，六页95dpi逐页核版，无公式溢出或未定义引用。

下一步直接研究M(s)=R†e^(iK⊥s)R的投影谱/协方差、其守恒扣除与有限相对时间涨落重求和。先找能接入的相邻charge相关resolvent或已有projection-kernel方法。Xu与Rajgadia–Xu提供不同极化/态的具体已有构造，需匹配当前charge-conserving Hamiltonian，不能直接搬用KM态结论。体的Itakura/FLZ质量算符仍作为目标，不任意指定P⁻使之相符。
