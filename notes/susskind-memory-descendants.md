# 恢复 raw 后代，组装二次记忆的双 chaos 残量

2026-11-06，本地。[七页PDF](../output/pdf/de-sitter-susskind-memory-descendants.pdf)、[LaTeX](../manuscript/susskind-memory-descendants/main.tex)、[来源](../literature/susskind-memory-descendants-sources.md)。无commit/push/ZIP。

本轮恢复旧四个quartic后代的W0/W2，并把旧W4 Gram接到同一半填充迹。结果：扣除指定旧方向后，集体二次记忆系数mZ2同时留下Q4与二次输出缺额残量。它们等时正交；完整时间谱仍待算。

## 有限 N 的 raw 字典

标准CAR、完整重叠complex Gaussian p-body，同J相邻charge h=H_m、k=H_(m−1)，m=N/2，归一maximally mixed迹，annealed HS。ℏ1。Rj=Rem(k^j)，A=R2+(2p−m)h²，Tcal_j=Rj+(jp−m)h^j。

E0={A,h²}/2，M0=hAh，S0={h,Tcal3}/2−E0−M0，C0=Tcal4−2E0−M0−2S0，Z0=h4。按五字(R4,{h,R3},{h²,R2},hR2h,h4)：

\[
E_0=\tfrac12\{h^2,R_2\}+(2p-m)h^4,\quad
M_0=hR_2h+(2p-m)h^4,
\]
\[
S_0=\tfrac12\{h,R_3\}-\tfrac12\{h^2,R_2\}-hR_2h+(m-p)h^4,
\quad C_0=R_4-\{h,R_3\}+hR_2h.
\]

所以Ecal=C0−2S0+(E0+M0)/2精确。对五字线性组合应用已有D4/D13/D22/D121/D40即可恢复Q_P=ΔP/2。W2P=Q_P−2eP、W4P=P−Q_P+eP，eP=EP。

若Q_P=Σ(g⁺Z⁺+g⁻Z⁻)，在共同rank域p≤min(m−1,N−m)中，

\[
e_P=\frac12\sum_a[C_a g^+_P(a,a,0)+mC_a^-g^-_P(a,a,0)]I.
\]

这给精确W0。N4p2m3不在该共同域：上charge最高rank1，下charge额外有rank2，不能套只含共同rank的式子。直接CAR五字公式仍有效。半填充DS最终满足共同域。

## 完整后代的 W2 表

λ_N=p²/N、q_N=e^−4λ_N、C_N=Etr h²。写W2P=C_N[p ell_P+k_P]B+C_No_PO+j_PJhF+C_Ns_PS_N+oL²。下表λ/q保留N下标含义：

| P | ell | k | o | j | s |
|---|---|---|---|---|---|
| E0 | 2 | −2λq² | 1 | 0 | 0 |
| M0 | 2 | −2λ(2q−1) | q² | 0 | 0 |
| S0 | 2q | −2λ(q²−q+1) | q | 1 | 0 |
| C0 | 2q² | −2λ(8λq²−3q²+2q+2) | 1 | 2q | 8λq² |
| Z0 | 0 | q²+2q+3 | 0 | 0 | 0 |

全部有理prefactors、固定缺额与全比例展开见symbolic-filters.json。每个A滤波全rank有界、偏差≤K(1+u+v+S)/p，沿上一轮六阶增长窗口得到O_L²(p^−1/4)。H滤波扣p主项后有线性缺额界；既有input factorial/output三阶矩保证L²。总余项只报oL²。

四个旧后代E0/M0/S0/Z0的s全为0。原始Q4组合的W2仍保留8λCq²S_N；它不是纯Q4。

## 固定半填充迹的桥接证明

旧W4 Gram原先在全Fock迹计算，不能未经证明直接混用。对L个受约束模式、t个occupied的局部项，半填充/全Fock迹比精确为

\[
R_N(L,t)=2^L\frac{(m)_t^\downarrow(m)_{L-t}^\downarrow}{(N)_L^\downarrow}
=2^L\frac{\binom{N-L}{m-t}}{\binom Nm}.
\]

L≤kp=O(sqrtN)时，logR=[L−(2t−L)²]/(2N)+O(L³/N²)，全局R≤exp(K_k p²/N)。无内部Gaussian弦的marked Wick块，每条弦的两个p元指标集在singleton模式有相反初始占据。集合大小相同；删去multiply-used sites后不平衡|2t−L|≤c_kK，K为总两两交叠数。

body-rank mark在这样的块内连接不同Gaussian弦，只作用交叠模式，不改变singleton平衡。展开交叠模式的占据projectors后，绝对系数≤exp(cK)乘K的固定多项式。先K≤K0时R→1一致；窗外用已有固定数p元集的交叠指数矩界，再K0→∞。这证明同一marked Wick字的两种trace极限相等，覆盖W4 Gram及固定energy moments，而非靠有限数值猜测。

相同方法给Y/U等中心Wick块的所需高矩界。线性rank过滤Gaussian场X若EX²=vI，Wick＋CS给EtrX4≤3v²；旧input-tail控制v。因此不能只靠Y/U的L²来控制乘积，但这里所需Etrh8/Y4/U4确实有界。

## 完整旧 span 的真实投影

Pi_N投影span(I,B,O,JhF)，T_N=(1−Pi_N)S_N。P_i=E0/M0/S0/Z0经1−Pi后=W4P_i+oL²。除σ8后的旧四方向W4 Gram determinant为

\[
16\lambda^3(1-q)^2(1+q)^4(1+q^2)d_3d_5(1+q+4\lambda q)^2
[1-q^3+4\lambda q^2(2q^2-q+2)]>0.
\]

d3=q²−q+1、d5=1+q+q²+q³+q4。故finite-dimensional投影连续。设hatPi_N投影span(I,B,O,JhF,E0,M0,S0,Z0)，旧Q4=W4(C0−αE0−βM0−γS0−δZ0)，则

\[
(1-\widehat\Pi_N)\mathcal E
=Q_{4,N}+8\lambda_NC_Nq_N^2T_N+o_{L^2}(1).
\]

不是先独立投影W4后冒充raw最优投影：先精确去低阶span，再证明剩余四后代趋W4，并核Gram非退化。Q4与T_N等时严格正交，极限范数平方V4+V2，V4沿用旧显式公式，V2=64Js8 q4(1+q)(1−q²+8λq)/[λ(1+q+4λq)]。后小λ，V4∼768Js8/(5λ)，V2→1024Js8。

## 接回 raw mZ2

旧精确分离mZ2−D2−Ecal=(r/2+r²/4)Lh²Y+2rJhU+3r²Jh²Y−({KA,Y}+Y²)/m，r=p/m。KA=a_m+c_m(h²−C)，P_m为h/h³多项式，D2=2rJhP_m+3r²Jh²KA−KA²/m−mr4h4。

D2精确在span(I,h²,h4)中，属于hatPi的span。r=O1/p、a_m=O(p)、c_m=O1，上述高矩及noncommutative Holder给余项O_L²(1/p)。所以

\[
(1-\widehat\Pi_N)mZ_2
=Q_{4,N}+8\lambda_NC_Nq_N^2T_N+o_{L^2}(1).
\]

这里是集体mZ2、固定半填充、明确八方向span的等时结果。没有把A²/m直接删去，其有限守恒部分已在D2。

## 验证与继续位置

3JSON：symbolic-filters；finite-descendants（720方向3600收缩＋720raw core恒等式，前两系统完整二次系数张量，较大半填充3方向；15个scalar均值）；assembly-certificate（123trace恒等式、10共同rank mean、投影/行列式/两小λ极限）。所有科学输入输出哈希记录于local-evidence.json。

保留两次失败：读取旧JSON中的lambda触发Python保留字解析错误；N4p2m3错误套共同rank均值式，检查暴露lower-only rank，现明确域并记录skip。PDF校对另发现乘号输入成tab+imes，修正后重新编译核版。七页最终双遍通过，逐页95dpi核版。

下一步求T_N的双能量自谱，以及Q4–T_N双能量cross。ordinary chaos等时正交不能推出插入H或exp(itH)后仍正交。之后扩W2中央charge窗、finite-memory两参数及体介子字典。不将当前等时组装称为完整Z2动力学谱。
