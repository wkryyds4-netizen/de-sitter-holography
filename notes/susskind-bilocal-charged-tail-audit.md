# Charged ordinary 尾的独立审计：归一物理导数闭合

2026-11-07，本地证明附录。对象和规范沿用 fourth-tail 两份草稿：共同原 Gaussian 耦合、允许 I/J 重叠的 complex p-body Hamiltonian、标准 CAR、β=μ=0 normalized full-Fock annealed Schatten 范数、string time、固定 J_s，且 λ_N=p²/N 位于正紧区间。本文核查的是有限 N、所有次数同时成立的包络；不是只先取固定次数的 DS 极限。下文 r≥2 固定；非偶整数范数以更高偶整数范数控制。

输入是已另行证明的：(i) pure Wick 重复 D/L 标记的 sector 一致阶乘界；(ii) covariance 的精确 normalized CP 路径；(iii) finite Jacobi 的精确强迫递推及有界递推系数；(iv) Ben Efraim–Lust-Piquard 的 Clifford Riesz 定理及其半群插值引理。文献定理的 full-Fock trace 适用范围必须保留。

## 1. 两类物理导数与 charge 系数

对 neutral sector families 写
\[
 T=\sum_i c_i^\dagger\,\cdot\,c_i,\quad
 S=\sum_i c_i\,\cdot\,c_i^\dagger,\quad
 D=Q-T,\quad L=\tfrac12(N-T-S),\quad E_D=D/p,\ E_L=L/p.
\]
这里 L 是 half-Clifford length，不能与 coupling number \(\mathcal N\) 混同。
若 \(\Delta_-f(Q)=f(Q)-f(Q-1)\)、\(\Delta_+f(Q)=f(Q)-f(Q+1)\)，则
\[
 [D,f(Q)]=(\Delta_-f)T,\qquad
 [L,f(Q)]=\tfrac12[(\Delta_-f)T+(\Delta_+f)S].                 \tag{1}
\]
第二式同时需要前后两个差分。

在远离端点的 charge 区取 \(t=Q^{-1}T\)、\(s=(N-Q)^{-1}S\)。直接 CAR 给
\[
 [D,t]=t-t^2,\quad
 [D,s]=-s+\frac Q{N-Q}ts-\frac{2Q-N}{N-Q},
\]
\[
 [L,t]=\tfrac12\left[-t^2+\frac{N-Q}{Q}st+\frac{2Q-N}{Q}\right],
\]
\[
 [L,s]=\tfrac12\left[-s^2+\frac Q{N-Q}ts-\frac{2Q-N}{N-Q}\right]. \tag{2}
\]
每步 t/s 是相邻 sector 的 normalized-trace Schatten 收缩。在两侧中央区，每个式子的系数有界，所以对长度 \(\ell\le2p\) 的路径 P，\([E_D,P]\)、\([E_L,P]\) 的 sector-sup 算符范数至多 \(K\ell/p\)。估计路径时留出 p+2 个 charge 的余量。

## 2. 普通 H 幂的归一导数

已证 pure 标记要求每个新物理 site 至少约束两个独立 p-set 标签。普通 \(H^w\) 内允许 Gaussian pair；此时一个标记可能只约束一个标签，但至少仍约束一个。对 k 个归一缺额标记、s≤k 个不同物理 sites，受限 covariance 和 ancilla contraction 给
\[
 p^{-k}N^s(p/N)^s=p^{s-k}\le1.                            \tag{3}
\]
不同标记落在同一 site 时先合并 membership 约束。对 D，\((wp-D)H^w\) 是两槽导数和；对 L 还含同一个 Hamiltonian coefficient 内的 I/J 交集标记，该标记要求两个独立标签，故满足同一上界。槽选择只产生 fixed-moment 的 degree 多项式。

由此（中央逐 sector 及聚合全迹）
\[
 \|H^w\|_r+\|E_DH^w\|_r+\|E_LH^w\|_r
 \le K_r(1+w)^{b_r}L_r^w\sqrt{w!}.                       \tag{4}
\]
外 charge 区使用无约束 covariance 粗界和二项尾。这里归一因子 1/p 是关键；(3) 不能用来证明未归一 ordinary charged 缺额本身有界。

## 3. 实际 covariance residual 的归一导数

记
\[
 A_v=p(\Phi-C\theta_v)Z_v,
 \quad\theta_{v,m}=\mu_m(vp),\quad Z_v=:\!H^v\!:.
\]
取 good 区 N/3≤m≤2N/3；先考虑 vp≤N/8。精确 CP telescoping 将 A_v 写成至多 2p 项，每项形式为
\[
 p\,C(Q)\,P_i\,g_i(Q)M_i,
 \quad M_i=(vp-D)Z_v\ \text{或}\ (vp-\bar D)Z_v,
 \quad\bar D=2L-D.                                      \tag{5}
\]
P_i 长度≤2p；g_i 是一个 charge 分母的逆乘至多 2p 个 reference 因子。所有路径所及 charge 仍在 N/4≤k≤3N/4。因 reference 因子位于 [0,1]，
\[
 |g_i|\le K/N,\quad
 |\Delta_\pm g_i|\le K(1+v)/N^2,\quad
 |\Delta_\pm C|\le K(p/N).                              \tag{6}
\]
最后一个不等式可由精确 C 比值直接得出，常数含统一 covariance 上界。

对 (5) 作用 E_D 或 E_L：与 C、g_i 的交换用 (1)，与路径的交换用 (2)，与 M_i 的作用用
\[
 E_\sigma M_i=vM_i-p^{-1}(vp-\sigma)M_i,
 \qquad\sigma=D,L.                                      \tag{7}
\]
右端是固定次数的重复 pure tags。于是每项仍有
\(K(p/N)\operatorname{poly}(v)L_r^v\sqrt{v!}\)
上界；2p 项求和只产生 \(p^2/N=\lambda_N\)。

若 vp>N/8，则
\[
 p\le8\lambda_+v,\qquad N/p<8v,\qquad N<64\lambda_+v^2. \tag{8}
\]
直接对 pC(P−θ_v)Z_v 用 CP 收缩、(1)--(2)、\(|\Delta_\pm\theta_v|\le2\)；所有 p、N/p 或 N 都被 (8) 吸收入 degree 多项式。

外 charge 区 A_v 的粗界为 O(p)，E_DA_v、E_LA_v 的粗界为 O(N)，乘各自 \(L_r^v\sqrt{v!}\) 包络。二项尾及至多一个额外 charge 位移吸收这些 N 多项式。因此得到真正的全 Fock 一致界
\[
 \|A_v\|_r+\|E_DA_v\|_r+\|E_LA_v\|_r
 \le K_r(1+v)^{b_r}L_r^v\sqrt{v!}.                       \tag{9}
\]

## 4. 乘积导数与 η 的闭合

Clifford Riesz 在 full normalized trace 下给
\[
 \|\Gamma_\sigma(X,Y)/p\|_r
 \le K_r\sqrt{\|X\|_{2r}\|E_LX\|_{2r}
                   \|Y\|_{2r}\|E_LY\|_{2r}},            \tag{10}
\]
其中 Γ_σ 是 D 或 L 乘积规则的双梯度项；D 的 complex-CAR 行列梯度只改变固定常数。对 disorder 再用 Hölder 即得 annealed 版本。非自伴输入须保留行、列两种梯度。

因此同时受控的 X,E_DX,E_LX 与 Y,E_DY,E_LY 在相乘后仍受控，使用固定 2r 阶。**不得把 (10) 当作 fixed-charge Riesz 定理。** 如需与中央系数界结合，可先用 χ(Q/N) 作中央 cutoff：χ 在 [1/3,2/3] 为一、支撑于 (1/4,3/4)，固定 Lipschitz 常数，则 [E_σ,χ]=O(1/p)。cutoff 外直接用上述 sector 粗界和二项尾。

记 finite Jacobi 系数 \(\beta_j=C\sum_{v<j}\theta_v\)。有 \(|\beta_j|\le B\)，中央又有
\[
 (N/p)|\Delta_\pm\beta_j|\le K(1+j)^b.                  \tag{11}
\]
低 vp 时逐 reference 因子差分给 \(N|\Delta\theta_v|\le Kv\)；高 vp 时 (8) 和粗差分给 \(N|\Delta\theta_v|\le Kv^2\)。把它们有限求和即足以得 (11)，不需优化几何衰减。

对精确式
\[
 pr_\ell=\sum_{v+w=\ell-1}A_v\diamond Z_w
        =\sum_{v+w=\ell-1}\mathcal W_{v+w}(A_vH^w),
\]
\[
 \eta_{j+1}=H\eta_j-\beta_j\eta_{j-1}-pr_j,
 \qquad\eta_j=p[Z_j-Q_j(H)],                              \tag{12}
\]
先完全展开递推。每项是
\[
 b_{j\ell a}(Q)H^a\mathcal W_{v+w}(A_vH^w),
 \qquad a+v+w\le j-2.                                   \tag{13}
\]
b 是有限 β 路径权重；总绝对系数至多 K^j，(11) 给其归一 charge 差分至多 \(\operatorname{poly}(j)K^j\)。E_D/E_L 与 Gaussian chaos 投影交换；最高 chaos 投影在 degree n 上的同阶 Schatten 代价≤2^n。

对 A_vH^w 用一次 (10)，对 H^a 与其投影结果再用一次 (10)。最多只需固定 4r 阶输入（放宽至 8r 也可），而非随着 j 递归增大 Hölder 阶数。b 的导数单独用 (1)、(11) 处理。阶乘满足 \(\sqrt{a!v!w!}\le\sqrt{j!}\)。外区可直接对 \(p[Z_j-Q_j(H)]\) 用 D,L=O(N) 粗界，吸收入二项尾。于是
\[
 \|\eta_j\|_r+\|E_D\eta_j\|_r+\|E_L\eta_j\|_r
 \le K_r(1+j)^{b_r}L_r^j\sqrt{j!}.                       \tag{14}
\]

## 5. Ordinary charged 系数的精确重组

逆 Jacobi 展开为
\[
 H^n=\sum_{j=0}^n a_{nj}(Q)[Z_j-\eta_j/p].              \tag{15}
\]
a_nj 是长度 n 的有界非负 Jacobi 路径和，故 \(\sum_j|a_{nj}|\le K^n\)。必须保留这些系数的 charge 差分。令 \(X_m=(m-N/2)/p\)，精确比值
\[
 \frac{C_m}{C_{m-1}}=
 \frac m{m-p}\frac{N-m+1}{N-m+p+1}
\]
给中央 \(N|\Delta_-C_m|\le K(1+|X_m|)\)。结合 §4 的 θ 差分，逐路径 telescoping 得
\[
 N|\Delta_-a_{nj,m}|\le K^n(1+n)^b(1+|X_m|).            \tag{16}
\]
这是首阶差分，不能把 charge 因子重复乘 n 次；逐路径乘积的差只有一次差分插入。

由于 coupling number \(\mathcal N\) 与确定性 charge 系数交换，而 D 不交换，由 (1)、(15) 得精确恒等式
\[
 (p\mathcal N-D)H^n=
 \sum_j\left[
 a_{nj}(jp-D)Z_j-(\Delta_-a_{nj})TZ_j
 -a_{nj}\mathcal N\eta_j
 +a_{nj}E_D\eta_j
 +\frac{\Delta_-a_{nj}}pT\eta_j\right].                 \tag{17}
\]
式中第一项是 pure tag；第三项用有限 coupling degree 上的同阶 OU/Chebyshev 乘子界；第四项用 (14)。第二和第五项用 T 的逐 sector 收缩及 (16)，再以 binomial X 的固定矩和已证 Z_j、η_j 的中央 sector 包络聚合。外 charge 区的 N 多项式损失由二项尾吸收。

最后 \(\lambda_N\mathcal N(\mathcal N-1)H^n\) 由同阶 Gaussian degree 乘子和 (4) 控制。故
\[
 \boxed{\sup_N\|[p\mathcal N-D-
 \lambda_N\mathcal N(\mathcal N-1)]H^n\|_r
 \le K_r(1+n)^{b_r}L_r^n\sqrt{n!}.}                      \tag{18}
\]
固定左右 H 能量插入由一次固定高阶 Hölder 和 (4) 处理；固定时间导数只多一个 degree 多项式或有限阶乘位移。因此在每个固定实时间窗，ordinary charged 时间级数的 actual 有限 N 尾满足
\[
 \sup_{N,|s|\le S}\left\|\sum_{n>K}\frac{(is)^n}{n!}G_n\right\|_r
 \le K_{r,S}\sum_{n>K}(1+n)^{b_r}
                  L_{r,S}^n/\sqrt{n!}\longrightarrow0.  \tag{19}
\]

结论范围：这是此前缺失的 charged ordinary finite-N 高范数可和尾。完整实际 bilocal 的四点传递还应在主稿中逐项接回已有返回、two-component、pinching 与固定次数共同多源矩接口；(18)--(19) 本身不声明新的体空间字典或有限 G 引力动力学。
