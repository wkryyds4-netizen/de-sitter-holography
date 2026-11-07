# 有限时间记忆涨落的两分量归约

2026-11-07。本轮正式稿 `manuscript/susskind-memory-two-component/main.tex`，八页阅读版 `output/pdf/de-sitter-susskind-memory-two-component.pdf`。接 uniform-bound 的实际场二阶界；上一轮 goal 为 progress，本轮继续产生新证明、三份科学输出与完整阶段稿。

## 当前结论

模型、态和单位不变：标准 CAR、同 J 的重叠 complex p-body Gaussian 模型，β=μ=0 全 Fock 归一 annealed trace，ℏ=1，string 实时间；偶 N、λ_N=p²/N→λ>0，固定 J_s 和有限实时间窗。只要求 λ_N 留在正紧区间时，场归约由子列论证仍成立；显式谱中的 q=e^(−4λ)、C=J_s²/λ 需选定 λ 极限。

定义精确有限 N 分量（P_h 保留简并能量块）

\[
\mathsf A_m(s)=(1-P_h)m e^{-ihs/2}(\Lambda_m-\Lambda_{{OU},m})e^{-ihs/2},
\]
\[
\mathsf B_m(s)=(1-P_h)m\epsilon_m e^{-ihs/2}\Phi_m(e^{i\alpha_mhs})e^{-ihs/2}.
\]

本轮证明 `sup_|s|≤S ||Fscr−A−B||₂→0`，A、B 各自全迹统一有界。于是实际 memory 两点函数归约为 AA+AB+BA+BB 四项，对紧 memory times、所有外部实时间一致。共同 J 仍保留；没有把 charged force 换成独立 probe，也没有得到这四项的完整大 N 终式。

## 四 u 项为什么可以去掉

先在真正 q-Gaussian 家族中取 h、u、v 正交。BKS 第二量子化的条件期望给 `E_(h,v) W(h,u)=E_h W(h,u)`：u tensor 被投影掉，保留的纯 h tensor 相同。因此两个完整 block 的 trace 满足

\[
\tau_q[W_1(h,u)^\dagger W_2(h,v)]
=\tau_q[(E_hW_1)^\dagger(E_hW_2)].
\]

这不是经典独立性，不适用于任意交错的多色 word。有限 N 下，用两个独立辅助 J_u、J_v 精确重写条件平均的平方范数；颜色只限制同色 Gaussian 配对，允许的 chord 仍有原来的 q crossing。三色 fixed-word 全迹矩收敛与上式，使四 u 多项式条件平均强 L² 趋于能量多项式。

实时间展开按 h 总次数 n 分组：

\[
\mathcal I_s^{(1)}=
\sum_{n\ge0}\frac{s^2(is)^n}{(n+2)!}
\sum_{a+b+c=n}E_u[u h^a u h^b u h^c u].
\]

每个 word 的 L² 范数 ≤ `sqrt((2n+7)!!) M^((n+4)/2)`；三元组数抵消分母后，整层界为 `S^(n+2) sqrt((2n+7)!!) M^((n+4)/2)/(2 n!)`，对任意固定 S 可和。先固定次数取 N 极限，再解除截断，没有使用所有次数统一的 chord 渐近误差。

q 条件平均由 C* 定理属于 C*(h)，在紧谱上连续；clip 到该谱给有限 N 全实轴上的 bounded energy reference。真实相位仍用真实 h。广中央 `|1−α|≤4p/N`，Duhamel+Schatten-10 给差界 `|1−α||s|³ sqrt(9!!) M^(5/2)/2`。`mε²≤16λ`，外区靠 `N² exp(−N/8)`。因此 `(1−P_h)mε² e^(−ihs/2)I_s^(α)e^(−ihs/2)→0`。同副本 Dyson 的三阶余量亦消失，结合首返守恒结果得到两分量归约。

## 显式守恒背景

方差 C、monic q-Hermite W_n；`v=αh+√εu` 且 α²+ε=1。BKS 左右 Wick 递推给完整公式

\[
E_h[uW_n(v)u]=Cq^n\alpha^nW_n(h)
+\epsilon C^2[n]_q[n-1]_q\alpha^{n-2}W_{n-2}(h).
\]

两端 u 分别与内部 v 收缩的第二项不可漏掉。比较沿归一化 α=√(1−ε) 的展开与 frozen-α Dyson，得到

\[
\iota_s=\frac C2P_q[(N_q-x\partial_x)e^{isx}]-C^2D_q^2e^{isx}.
\]

`N_q W_n=nW_n` 是能量 Hermite number，非有限 N coupling chaos number；`D_q W_n=[n]_q W_(n−1)`。`x∂x` 来自归一化漂移，不可漏。背景二、三阶导数为 `(2+q)C²I`、`iC²(1+4q+3q²+q³)h`。这是现有 q-Wick 结构的显式应用，不声称发现新 q-Hermite 定理。

## 两个涨落分量的低阶混合谱

写 `B(s)=is B1−s² B2/2+…`，有精确有限 N 公式

\[
B_1=m\epsilon\alpha(1-P_h)\Phi(h),\qquad
B_2=m\epsilon(1-P_h)[\alpha^2\Phi(h^2)-\alpha\{h,\Phi(h)\}].
\]

令 `Z=h²−C_m I`、`BZ=Φh−barκ h`，则 B2/(mε) 的 pinching 后式可换为 `α²(Φ−κ2p)Z−2α J_hBZ`；未投影的差显式为 `(α²κ2p−2αbarκ)h²+α²C_m(C_m−κ2p)I`。κ2p 用 μ 的正部乘积延拓定义，可能不对应物理 rank。

中央紧窗的有限 rank slope `p[μ(kp−d)−μ(kp)]→4λq^k d`，配 input factorial tail、output sixth moment 与已有 S/Lin/F 桥接，给

\[
B_1\to F,\qquad B_2\to d(Y_-+Y_+)+bJ_HF,
\quad d=4\lambda Cq^2,\ b=2(q-1).
\]

线性 Gaussian 误差乘 h 用 `||J_hX||₂²≤3C_m||X||₂²`，不从一般 L² 收敛裸乘无界 h。pinching 的移除覆盖能量多项式后代，使用旧 Mehler 谱的无对角质量及高能量矩。

设 `e=(E+E′)/2`，w、η 为旧 Y± 矩阵，`VF=64Js⁶q²`。新混合项是 `ν_(B1,B2)=b e VF μ_q`，`ν_B2=2d²(w+η)μ_(q²)+b²e²VF μ_q`。charged 的 `A0→Y−`、`A1→R+aJ_HY−`，`a=−(1−q)²/(1+q²)`，于是

\[
\nu_{A_0,B_2}=d(w+\eta)\mu_{q^2},\qquad
\nu_{A_1,B_2}=ad e(w+\eta)\mu_{q^2}.
\]

`w+η>0`，二阶 cross 非零；一阶与 B1 的 cross 为零。这里只证明系数谱：任何允许逐阶匹配的时间核必须恢复它，尚未交换导数与未证明的有限时间极限。实际 B1/B2 在广中央有 x 无关界，在外区为 O(p)；结合 binomial 尾及旧 Y/L 全迹界，得到上述全迹联合谱，不需要假定所有命名场任意延拓均强收敛。

## 检查及失败记录

- `four-copy-q-checks.json`：6组有理 q/C，2520 conditional trace、7200三色 replica、66 OU Wick、42 Taylor 恒等式，两个独立算法（q-Fock 作用与有色pairings），首次通过。
- `four-copy-finite-time-checks.json`：4个完整 CAR 系统、2个固定 seed、3个正负实时间，共24积分；全部 Gaussian 坐标及三配对，8/12阶 quadrature 差最大1.94e−15。初选3个高对称系统不能见证动态残量，末尾assert失败且无JSON写出；保存源/失败文本后加入 N5p2m3 内部扇区，动态norm最大约.05216。不是DS拟合，不是区间积分证书。
- `neutral-coefficient-checks.json`：12个完整坐标有限恒等式，最大误差1.87e−16；360个有理 slope 支配界，80/120位差已核。有限相对误差完整保留；p128的测试行仍可有约8.67%偏差，不把这些小尺寸当渐近速率拟合。
- 首编译因 LaTeX 内置 AA 宏冲突失败，第二次长公式/行有 overfull，源和日志均保存。改用 MA、分行后最终八页双遍编译干净，95dpi全八页视觉检查通过。PDF marker本阶段成功一次。三个只读协作者分别复核q公式、强回归与低阶谱；没有写文件或Git。

下一步是 A/B 的共同 J 全时间联合核及体介子字典。可从 `susskind-memory-neutral-series-seed.md` 的 exact neutral series 接口出发；先独立核该分解与高阶收缩，再求所有次数的联合谱。当前不是完整理论完成，goal继续active，仅本地保存。
