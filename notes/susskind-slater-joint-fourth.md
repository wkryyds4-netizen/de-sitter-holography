# 同一 J 的联合四点与 Hamiltonian 前两层

2026-11-07，本地阶段。完整定义、有限配对表及解析证明见 [LaTeX](../manuscript/susskind-slater-joint-fourth/main.tex) 与 [PDF](../output/pdf/de-sitter-susskind-slater-joint-fourth.pdf)；[来源核读](../literature/susskind-slater-joint-sources.md)。

模型延续 slater-ladder：N 偶、m=a=N/2，标准 CAR，原 H=√w ΣJ_AB C_A†C_B，overlapping A/B 允许。φ_N=E_J〈O|·|O〉，固定 Slater O，非 full-Fock 迹、未假设 stationary。p→∞、p²/N→固定 λ>0、J_s>0；ℏ=1、H 为 string 能量，h=H/√(wD) 无量纲。所有字长固定，未交换增长时间/次数极限。

令 L=C(m,p)C(a,p)，D=C(m,p)C(a+p,p)，M† 为原 J 的纯 O→E creation block 除 √L。q=e^(−4λ)，s=√q，γ=q^(1/4)，δ=√(1−s)。精确 φ_N(hM†)=φ_N(Mh)=γ_N=√(L/D)。因此独立 matter probe 的数量选择律不能用在此处。

## 1. 有限 N 两配对表及极限

raw 字母 H,+,− 分别指 H/√w、R†、R。HH edge 是 G；其余 edge 为 S，必须取 pure selected block 并按 M 端点固定方向。forward 表示 creation 先作用（从右读）；++/−− edge 不相容。

记 r_N=C(m−p,p)C(a−p,p)/L，T=C(m−p,p)C(a,p)，W_k=C(m,p)C(p,k)C(a,k)，U=Σ_k W_k C(a−k,p)C(m−k,p)。三配对为 disjoint (01)(23)、crossed (02)(13)、nested (03)(12)。

| 情形 | disjoint | crossed | nested |
|---|---:|---:|---:|
| GG | D² | K | D² |
| SS 均 forward | L² | L²r_N | L²r_N |
| GS，S forward | DL | LT | S outer: DL；G outer: U |
| S outer forward、S inner backward | 0 | 0 | L |
| G outer、S inner backward | 0 | 0 | L |

其他 backward 图为零。每个 word 对三个 pairing 求和，除 D^(n_H/2)L^((4−n_H)/2)。K 是 HHHH crossed 原始计数，用 Kota (28) 可独立预测。全 SS 用旧 even-CAR active-support 证明；GS crossed 用既有 rank-p covariance；GS nested 按 k 数可用 O/E slots；backward 只能全块撤销第一步。这些步骤详写于正式稿。

L/D→s，r_N、T/D、U/(DL)→q。最后一项由 p−k→Poisson(2λ)、固定缺额窗上一致收敛与 [0,1] 支配证明。K/D²→q 用已有 fixed-charge H-only 律及 Schur。存活图的 endpoint weight 为 HH:1、MM†:1、HM†/MH:γ；crossing 添 q，nested 仅 inner 是 S 时添 q。**此处只证明两弦规则**。

## 2. 直接得到的正空间

先由有限 N 的正性取 degree≤2 Gram 极限，再商掉 null vectors，不预设全阶 GNS 表示。MhΩ=γΩ 的精确误差为 E||(M_Nh_N−γ_N)Ω_O||²=1/D_N。

定义

\[
U=(M^\dagger)^2\Omega,\qquad
V=(M^\dagger h-\gamma(M^\dagger)^2)\Omega,
\]
\[
D_v=(hM^\dagger-M^\dagger h-\gamma I)\Omega,
\qquad R_v=(h^2-I)\Omega-sU-2\gamma V-\gamma D_v.
\]

四向量互相正交，亦正交于真空和一次向量；平方范数依次为 2q、q(1−s)、1−q、(1−s)²(1+s+2q)，固定 0<q<1 时严格正。归一为 e₂,e_r2,e_d,e_rr，另设 e₁=M†Ω、e_r=(hΩ−γe₁)/δ，得到七维 degree≤2 空间。第一层作用为

\[
h\Omega=\gamma e_1+\delta e_r,
\]
\[
he_1=\gamma\Omega+\gamma\sqrt{2q}\,e_2
+\delta\sqrt q\,e_{r2}+\sqrt{1-q}\,e_d,
\]
\[
he_r=\delta\Omega+\gamma\sqrt q\,e_{r2}
+\delta\sqrt{1+s+2q}\,e_{rr}.
\]

这是 degree≤1 到 degree≤2 的实际乘法作用。未指定 h 对 degree2 的作用，不是封闭七态 Hamiltonian，也不是全有限时间传播。新交换子方向 D_v 的范数平方为 1−q；不能并入纯阶梯或只调整其标量系数。

## 3. disorder complement 的新校准

g_N=h_N−γ_N(M_N+M_N†) 精确移除原 H 的 selected block 及其伴随。其 Gaussian 坐标与 M 独立，但共同 CAR/state 矩不因此分解。

\[
\varphi(g^2)=1-s,\quad
(g^2-(1-s)I)\Omega=R_v,\quad
\frac{\varphi(g^4)}{\varphi(g^2)^2}=2+s+2q.
\]

q>1/4 时最后一式大于 3；对任何 0<q<1 均不同于标准 q-Gaussian 的 2+q。这否定了把实际 remainder 随手补成另一个标准 q 模式的方案，留下更具体的全联合规律问题。

## 4. 证据与继续位置

- exact-checks.json：四有限系统全部 324 个四字、36 个二字、972 个配对值、2,000,499 个 Gaussian label assignments、116 个 covariance columns；六组 λ=1/2 渐近精确比值。H-only K 在此 driver 中仅枚举。
- gram-checks.json：62 个精确符号恒等式，包括 7×7 Gram、第一层作用、g⁴ 与四个有限 Mh residual。独立代数核对不等于独立全阶证明。
- covariance-checks.json：既有 EGUE 谱/第四矩与保存 CAR 数据作 12 项独立标量比较，无重跑旧 driver。

三个科学输出皆 exclusive-create。参数/依赖版本/driver 哈希在各输出；编译与渲染绑定另见 local-evidence.json。旧成功或失败输出不覆盖。下一任务为三弦/六阶 same-J 有向局部规则，再给 all-degree backward 控制和正表示；不能把两弦因子化猜想升级为已证定理。本阶段本地保存，Git 与 issue 已在本轮开始前按用户授权集中完成，不再次提交或评论。
