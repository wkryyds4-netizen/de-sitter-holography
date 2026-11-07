# 完整双局域的制备协方差

2026-11-08。正式证明、公式编号及全部约定见[八页稿](../manuscript/susskind-selected-bilocal-covariance/main.tex)，[PDF](../output/pdf/de-sitter-susskind-selected-bilocal-covariance.pdf)。这一步接上一阶段单插入投影，保留投影之外的共同涨落。

## 1. 对象与归一

同一原始复 Gaussian overlapping p-body CAR 模型，半填充，先有限 N 精确 pinching，再固定 lambda=p²/N 的双标度极限。C=Js²/lambda，q=exp(−4lambda)，zeta=√q，h=H/√C，W=√lambda U；h,U 为正交、单位方差 q-Gaussian 场。时间 t=√C t_phys、T=√C T_phys；string 变量为 u=2√lambda T。

实际单选态为 sigma(Z)=tau(rho Z)，rho=1−zeta(1+lambda)+zeta(h−W)²。能量限制态为 tau_f(Z)=tau(f(h)Z)，f=1−zeta+zeta h²。X,Y 限为各自线性含一个 marked 插入的能量中心化场，允许 h 左右 dressing、chaos 滤波、有限和与已有尾控制的时间极限。P 是联合 L² 中到 h,U 多项式闭包的正交投影。

使用 W 符号 PX=Σ beta_ab h^a W h^b，dnu=K_q(x,y)dmu_q(x)dmu_q(y)，(Ccal beta)(x)=∫K_q(x,y)beta(x,y)dmu_q(y)。因此 ||PX||²=lambda ||beta||²，不可漏掉这一归一。

## 2. 两插入恒等式与残余

稿定理1给

    sigma(X†Y) = tau_f(X†Y) + zeta lambda²(1+q)<Ccal beta_X,Ccal beta_Y>.
    mu_X = −zeta lambda ∫(x+y)beta_X dnu.
    Cov_sigma(X,Y) = tau_f(X†Y) + Gamma(X,Y),
    Gamma = zeta lambda²(1+q)<Ccal beta_X,Ccal beta_Y> − conjugate(mu_X)mu_Y.

证明直接使用既有局域 CAR 源：rho−f=zeta[(W²−lambda)−{h,W}]。相邻两个边界 W 相互配对的项由 −lambda 抵消；其余各接一个不同 L mark，恰等于两个投影的图。D-W 为零，高阶源累积量为零。单边界 W 与中心化两场的项由投影和奇 U 数消失。最后调用已有相邻四 U 的两步 Mehler 恒等式。此论证只用于指定单标记场类，不是 P 的乘法性。

Q_X=X−PX 后，完整协方差为 tau_f(Q_X†Q_Y)+Cov_sigma(PX,PY)，两核都正。特别 PJ₂=0，而 sigma(J₂²)=4lambda+16lambda²q/(1+q)>0；直接乘投影会丢失完整方差的首项。

未扣均值的制备项是正 Gram 核，Gamma 本身可以不正。取 X=S₂−kappa₂H₂，beta=x+y，得到 Gamma(X,X)=√q lambda²(1+q)²(1+q−4√q)。物理曲线上 q>(2−√3)² 时为负。形式有理点 q=16/81、lambda=1/5 独立精确检验给 −1768892/119574225，完整方差仍为正的 806826868/773246655。

## 3. 完整双局域与共同中心时间

沿用上一稿完整 beta_(T,t)=−2 exp[iT(x−y)]Σ_j Xi_j L_j；Xi_j 保留所有低 chaos。基线使用共同原始耦合的正矩阵谱 dnu_ab(x,y)=tau[X_a†P_h(dx)X_bP_h(dy)]，x 为输出、y 为输入。X_a 包含 charged L_n 和 neutral N_j，时间系数分别为 (is)^n/n!、2c_j(s)。实际 Volterra 滤波给 b_a，基线是 Σ_ab∫f(y) conjugate(b_a)b_b dnu_ab。它保留 D-L、L-L 连接及正交残余，仅依赖 T₂−T₁。

完整协方差=上述基线+Gamma。固定内时间下，两中心共同移 R 后，Gamma→0：先在固定 lambda 双标度理论内，条件 Mehler 积分的 Riemann–Lebesgue 衰减给多项式情形，L² 收缩和稠密性推广到完整符号；紧能量支撑控制均值。未给有限 N 递归、衰减速率或联合时间极限。

短时 beta=e^[iT(x−y)]{−(4i/3)lambda q t³−lambda q(1−q)(x+y)t⁴/3+O(t⁵)}。定义 C₂(T)=∫e^[iT(x−y)]K_(q²)dmu²、D_q=(16/9)√q lambda⁴q²(1+q)，则

    Gamma(T,t;T,t)=D_q[C₂(T)t⁶+(1−q)C₂'(T)t⁷/4]+O(t⁸).

中心非零时不能假定内时间方差为偶函数。相同内时间、不同共同中心的完整方差差等于 Gamma 的差，基线精确消去；T≠0 时首项 D_q[C₂(T)−1]t⁶ 为负。完整方差自身从 t⁴ sigma(J₂²)/4 开始。

## 4. 误差与执行证据

稿式(17)–(19)用 r=(1−q)^−1/2、A=qr|t|、B=3r|t| 给显式符号范数和阶乘 degree 尾。尾界不使用浮点计算出的范数，不包含能量/时间求积或舍入误差。完整基线由旧 actual/local L² 尾及有界 f 控制。

- exact-checks.json：386 个完整选态比较、16 奇边界、16 偶边界、2 曲率、2 残余 Gram；两矩阵14个主子式非负。形式有理参数不冒充物理指数曲线。
- gram-checks.json：物理 q=.25,.5，lambda=−log(q)/4；30 对角时间点、两个复 6×6 不等时矩阵、6 个短时点，另只读校准10个旧均值。次数/能量/时间求积从24/48/24增至32/64/32，scalar Taylor cutoff=80。
- 最大对角加密绝对差3.253e−19，相对差8.334e−14；最大解析degree界3.508e−16。均值相消与扣均值条件数保存在原始 JSON。两个 connected excess 矩阵最小本征值约−8.19e−7、−2.21e−7，未裁零。
- sign-witness.json：上述负修正及正完整方差的独立有理数见证。

数值只算制备修正及基线消去后的完整方差差，未积分绝对完整基线。输入哈希、环境版本、最终PDF/双遍编译/渲染及父证据核对见[证据清单](../results/susskind-selected-bilocal-covariance/local-evidence.json)。成功科学输出和脚本不重跑、不回写。

## 5. 继续位置

现在回到 Susskind 静态补丁的态、能量谱、时间单位与空间/质量生成元：先列原文中的算符与规范，再检查同一候选体解释能否同时容纳已算均值及非平稳协方差。Mehler 指数不能直接称为体质量。四个双局域的选态源分区、有限 G 散射和完整对偶仍需另做；两个单标记场的公式不能迭代到多标记乘积。
