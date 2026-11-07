# 四次记忆 core 的二次 Gaussian 收缩

2026-11-06，本地。[五页 PDF](../output/pdf/de-sitter-susskind-memory-quadratic-contraction.pdf)、[LaTeX](../manuscript/susskind-memory-quadratic-contraction/main.tex)、[来源核读](../literature/susskind-memory-quadratic-contraction-sources.md)。本阶段不 commit/push/ZIP。

本轮把上一轮尚未求出的 W2 Ecal 写成有限 N 的精确协方差映射，并给出 H² 投影及旧二次方向的完整有限系统 Gram。没有据有限系统宣布新的 DS 通道：去掉逐样本单位算符部分之后，两个小系统闭合，N5p2m3 则仍有非零动态残量。下一步已具体化为这一个残量的双标度权重，而不是再算纯 Q4 的低矩。

## 状态、规范与 Gaussian 分解

标准 CAR，含指标重叠的 complex Gaussian p-body H，E J_IJ J_KL=w δ_IL δ_JK。固定 charge m≥p，annealed normalized HS 内积，ℏ=1。所有新公式均有限 N 精确，不取渐近。计算使用 w=1；项目物理规范 w=4^p(Js²N/p²)(p!/N^p)²。

Ecal=Σ(L²c_i)†L²c_i+L²A/4，A=Σ(Lc_i)†Lc_i。前轮已经将 mZ2 的显式 first-return/drift 守恒部分分离，余项在 DS 消失；本轮求 core 的低 coupling chaos，不替换完整 finite-memory observable。

对四次齐次矩阵多项式 P，实 Gaussian 坐标方差 v_a，Δ=Σv_a ∂a²，Q_P=ΔP/2，e_P=EP，则

\[
W_0P=e_P,\qquad W_2P=Q_P-2e_P,\qquad W_4P=P-Q_P+e_P.
\]

这里 Δ²P=8e_P，E Q_P=2e_P。特别注意不能把 W2 写成 Q_P−e_P。来源中的谱 q-Hermite 正规序与此处坐标 chaos 定义不同；可能的 DS 等价须另证。

## 有限 N 协方差公式

h=H_m，k=H_(m−1)，c_i:H_m→H_(m−1)。令 Rem(X)=Σc†Xc，R_j=Rem(k^j)。R0=mI，R1=(m−p)h，直接 CAR 展开给

\[
A=R_2+(2p-m)h^2,
\quad Ecal=R_4-2\{h,R_3\}+\frac54\{h^2,R_2\}
+\frac72hR_2h+(4p-3m)h^4.
\]

独立同分布副本 (u,v) 仍在两个相邻扇区共享同一组新耦合。定义

Phi+(X)=E uXu，Phi−(X)=E vXv，Psi(X)=E uXv，
Gamma(X,Y)=ΣPsi(c†X)Yc。

Psi 是 d_m×d_(m−1) 矩阵上的映射，并非独立取两个 charge 的平均。其精确系综表达为 wΣ_IJ T_IJ^(m) X T_JI^(m−1)。中性 Phi 的 rank-ν 特征值是已有 κν(m)=w binom(m−ν,p)binom(N−m+p−ν,p)。

C=w binom(m,p)binom(N−m+p,p)，C−=w binom(m−1,p)binom(N−m+p+1,p)，a=(2p−m)C+mC−。

P+=Phi+(h)，P−=Phi−(k)，
U=Gamma(k,I)+Gamma(I,k)，
V=Gamma(I,k²)+Gamma(k,k)+Gamma(k²,I)，
U_h=Σ[Psi(hc†k)c+Psi(hc†)kc]。

定义 Ctr P=E[t²]P(h+tu,k+tv)。五个精确单收缩：

\[
\begin{aligned}
D_4&=Rem(3C_-k^2+\{k,P_-\}+Phi_-(k^2)),\\
D_{13}&=\{h,Rem(2C_-k+P_-)\}+V+V^\dagger,\\
D_{22}&=2CR_2+2mC_-h^2+hU+U_h+U^\dagger h+U_h^\dagger,\\
D_{121}&=Phi_+(R_2)+Uh+hU^\dagger+mC_-h^2,\\
D_{40}&=3Ch^2+\{h,P_+\}+Phi_+(h^2).
\end{aligned}
\]

依次对应 R4、{h,R3}、{h²,R2}、hR2h、h4。

\[
Q_{Ecal}=D_4-2D_{13}+\tfrac54D_{22}+\tfrac72D_{121}+(4p-3m)D_{40}.
\]

证明的关键排序：B=Rem(kv+vk)，E uB=U，E uhB=U_h，E Bu=U†；R3 的 v 线性项左乘 u 收缩给 V。所有次序保留。五个原始字也张成上一轮 E0/M0/S0/C0/Z0，因此同一套公式能补回其低 chaos，避免把 W4 后代直接当 raw 后代。

## 旧方向与有限系统精确投影

B=h²−CI，O=A−aI，F=W1T3−tbar h，tbar=<h,W1T3>/C。令 U_A=U+(2p−m)(P++Ch)，则

W1T3=Rem(2C−k+P−)+(3p−m)(2Ch+P+)−3[ah+(U_A+U_A†)/2]。

在代码中这提供同一协方差规范下的 F；其与原 rank 公式一致性的前驱结果见 cubic-linear。

对 P=Σx_a x_b P_ab（对 a,b 对称），
<W2P,W2S>=2Σv_av_b tr P_ab†S_ab。实坐标对角方差1，非对角实/虚方差1/2。矩阵系数为 dyadic Gaussian integers，检查后用 Fraction 累加 Gram，不进行 Monte Carlo 拟合。

未去逐样本 scalar 时，W2Ecal 对 B 的投影系数 cE 在三个系统分别为 5、1089/17、4596/41。均值 e_m 分别38、684、2268。

只扣 B/O/J_hF 后的 norm² 分别80/3、322410/49、3189421/65。但这个数字混入随机 scalar I，不能据此判断动态。进一步用 pi X=X−tr_m(X)I 对每个样本去 scalar，再投影到 piB/piO/piJ_hF：

| N,p,m | 独立旧方向 | 精确残量 norm² |
|---|---|---|
| 3,1,2 | piB,piO | 0 |
| 4,2,3 | piB,piO,piJ_hF | 0 |
| 5,2,3 | piB,piO,piJ_hF | 248945355/14903 |

最后一个系统的三个系数为 (2106991/44709,1499753/29806,−77111/89418)。保存的整数耦合样本给其 residual D 的精确正见证

tr_m[h,D]†[h,D]=12991771132439678939/2220994090。

这证明有限 N 下不能用上述有限旧方向和 scalar I 普遍闭合。pi 并非完整 commutant 投影；D 可以同时含其它守恒与动态部分，非零对易子只证明动态部分存在。尚无 DS 存活结论。

物理方差从 w1 恢复后 W2Ecal 缩放 w²、norm² 缩放 w4、对易子 norm² 缩放 w5。不能直接把上述 w1 数字用于大 N 推断。

## 一般 p=1 控制

H=dGamma(J)，A=dGamma(J²)，[H,A]=0，Ecal=dGamma(J4)。基本 GUE 的 E KXK=tr(X)I（此处未归一 trace）给

Q_Ecal=3N A+2tr(J)H+mtr(J²)I，e_m=m(2N²+1)，w1。

每个样本都与 H 对易。这是一般 N 的解析结论；N2..5 的全部14个非空 charge 扇区另有56个恒等式核查。说明不能拿一个 p1 的正 residual norm 当新动态通道。

## 证据与继续位置

4份最终 JSON 均 passed，不应覆盖：

- finite-contractions.json：3系统45/666/5050个对称 quadratic coefficient blocks，完整4×4 Gram、H²及旧方向投影，直接 raw 恒等式/独立 polarization。
- scalar-removed-checks.json：去逐样本 scalar 的完整 Gram；额外45+666+96=807个 independent polarization directions，小两系统穷尽，大系统96方向；N5残量的精确对易子见证。
- onebody-controls.json：14系统56精确式，涵盖 raw core、Q、mean 和 commuting。
- linear-input-audit.json：直接从CAR的三次导数T3做Gaussian极化，3系统全部实坐标加整数样本共148方向，独立核对本轮使用的W1T3/F映射。

初稿 PDF 首次双遍编译有2.48915pt opening paragraph 溢出；保留 initial-main.tex/initial-build.log，拆句后重新双遍通过，最终五页逐页核版。没有计算失败，也未改已有 quartic 阶段文件。

下一步以 general-N 的 bilinear coefficient map 为起点，按两个输入 coupling irreps 和输出 rank 分解 Psi/Phi/Rem 的组合；先合并 O(p)/O(m) 项，再取 rank deficit 及 charge 极限。需要实际计算 residual norm/cross 的 DS 极限，不能直接给 Phi(h²) 填 q² 或从三个有限点外推。随后恢复 W4 后代低收缩，求 raw Z2 联合谱及有限 memory 两插入协方差，继续接体空间与介子质量字典。goal 仍 active。
