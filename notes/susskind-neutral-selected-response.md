# 中性斜率的制备记忆与混合源

2026-11-08。本轮沿用实际 exact-pinched 中性斜率到 U 的已证桥，计算其完整两时间选态关联，以及能量—斜率混合观测。正式推导见[八页 LaTeX](../manuscript/susskind-neutral-selected-response/main.tex)与[PDF](../output/pdf/de-sitter-susskind-neutral-selected-response.pdf)。

## 1. 态与单位

h,U 为同一协变代数中的正交单位方差 q-Gaussian 场，W=√lambda U，epsilon=4lambda，q=e^-epsilon，s=√q。物理 string 时间 u=2Js t，K=h/√epsilon，U(u)=exp(iuK)Uexp(−iuK)，hbar=1。先有限N实际精确投影，再半填充fixed-lambda DS，最后小epsilon；无N一致率。U的迹方差为1，选态初始方差为1+s lambda(1+q)。

单选态密度分解为

rho=f_q(h)+s lambda(U²−1)−s√lambda{h,U}，f_q(h)=1−s+sh²。

后两项正是能量单独倾斜漏掉的部分。已有来源与精确actual接口不重算。

## 2. 两时间精确核

P_r为q-Mehler转移，P_r H_n=r^n H_n。令三能量链为 mu(da)K_q(a,b)mu(db)K_q(b,c)mu(dc)。两条非内部U配对的权分别1、q；对三段h，交叉协方差为q、q、q²，因此

tau[(U²−1)f(h)U g(h)U k(h)]=(1+q) integral f(a)g(b)k(c) dPi。

这是已有q-Wick/Markov结构在当前边界插入上的应用。odd U项对U(u)U(v)为零，所以精确

S_epsilon(u,v)=Ctilde_(d=1,epsilon)(u−v)+s lambda(1+q) B_epsilon(u,v)，

B=integral exp{i[u(a−b)+v(b−c)]/√epsilon} dPi。

B是正Gram核：A_u(b)=exp(−iub/√epsilon)P_q exp(iu·/√epsilon)(b)，B=integral A_u conjugate(A_v) dmu。B(u,0)=C_1(u)，B(u,u)=C_2(u)，故

S(u,u)=1+s lambda(1+q) C_2(u)。

单点均值为零，因而这是真实方差的时间变化。固定q时C_2随|u|→∞消失；这不意味着可以和小epsilon极限交换。

## 3. 紧时间余项

令g(u)=sech²u，R_1=(g/2)[u tanh u−3tanh²u]。取中间角theta1=pi/2+√epsilon x/2、相邻差epsilon z1和epsilon z2，则联合密度精确为

j=p(x+√epsilon z1,z1)p(x−√epsilon z2,z2)/m(x)。

m≤Cphi；两gap边缘有上一阶段的phi f主控。中心矩形|x|≤8√log(1/epsilon)、|zi|≤32log(1/epsilon)/pi内，

j=phi(x)f(z1)f(z2)[1−√epsilon x(z1−z2)]+O_weighted-L1(epsilon)。

外部联合尾由边缘的epsilon^32概率、Holder联合矩及Cauchy–Schwarz控制，避免在端点强行除小m。线性√epsilon项对x为奇，积分消失；相位误差O(epsilon)受多项式权控制，故B=g(u)g(v)+O_Ck(epsilon)。结合上一稿d1定理得

S(u,v)=g(u−v)+(epsilon/2)[R_1(u−v)+g(u)g(v)]+O_Ck(epsilon²)。

所有紧时间及有限阶偏导余项均受控。等权混合tau与sigma后仍有+(epsilon/4)g(u)g(v)，所以曲率通道的一阶相消不能提升为所有两点函数的相消。即使除以初始选态方差，S(u,u)/S(0,0)=1+(epsilon/2)[sech⁴u−1]+O(epsilon²)，记忆仍在。

## 4. 让奇混合源显现

取A(u)={h,U(u)}，其实际对应为归一Hamiltonian与实际中性斜率的反对易子。trace和任何能量单独密度的期望均为零，真正选态给

sigma(A(u))=−s√lambda H_epsilon(u)，

H=integral(a+b)² exp[iu(b−a)/√epsilon] K_q(a,b) dmu²。

中心角中(a+b)²=4y²+O(epsilon P)，并有全谱Cy²主控。上一阶段的带权余项直接给H=4sech²u+O_Ck(epsilon)，因此

sigma(A(u))=−2√epsilon sech²u+O_Ck(epsilon^(3/2))。

初始值精确为−2s√lambda(1+q)。这使混合源有了比“全态不平稳”更直接的可计算检验。这里是既有观测量的制备依赖，不据此识别新的体粒子。

## 5. 验证与后续

330个精确三块式、168个独立selected transfer、90个odd-source恒等式通过。14个独立角度积分/谱比较，最大差3.45e−14；28个两时间小epsilon点、20个混合均值点均核双维数。谱尾单独用实指数Taylor余项、Mehler L2收缩和Gauss多项式精确性控制，最大精细selected exact-arithmetic尾1.44e−35；它不覆盖浮点roundoff，角度quadrature亦非区间认证。记录版本、参数、相位相消与相对截断误差，科学输出成功后封存。

八页双遍编译与全页目检完成。下一步把更高中性残量N_j及charged-neutral共同耦合交叉项送进真实bilocal滤波；不能从当前单个斜率推完整双局域核或有限G。保留上一阶段曲率结果的适用范围，不追改历史稿。本轮仅本地，无Git或issue回复。
