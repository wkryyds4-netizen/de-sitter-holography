# 选态制备如何改变 Rindler 时间核

2026-11-08。本轮先回读 Susskind/MSS 的 singlet 谱目标、项目已复现的体介子谱及平移构造，确认时间 boost 核不能单独选定质量。本轮推进的是制备依赖的可计算接口。[LaTeX](../manuscript/susskind-selected-rindler-response/main.tex)与[PDF](../output/pdf/de-sitter-susskind-selected-rindler-response.pdf)区分精确恒等式与半经典预测。

## 1. 精确的源响应

先 N,p→∞，固定 lambda>0、半填充；同一 disorder、standard CAR，string 时间 t，hbar=1。记 epsilon=4lambda，q=exp(−epsilon)，e=H/√C，C=Js²/lambda，K=H/(2Js)=e/√epsilon，u=2Js t。e 的方差1，K的方差1/epsilon。x 是此前扣除守恒分量并归一的曲率 singlet，Hamiltonian crossing 为 q²。

在 charged 子代数上，已有选态密度 f_q(e)=1−√q+√q e²。定义归一热关联 G_(epsilon,beta)(u)=tau(exp(−beta K)x(u)x)/tau(exp(−beta K))。beta 是 K 的无量纲共轭变量；它对应物理逆温 beta/(2Js)，不是预先指定的体 Rindler 温度。直接微分给精确式

\[
C_1(u)-C_0(u)=\sqrt q\,\epsilon\,
\partial_\beta^2G_{\epsilon,\beta}(u)|_{\beta=0}.
\]

证明只有归一源导数：F_A''(0;Y)=tau(A²Y)−tau(A²)tau(Y)，tau(A)=0。不要求 A 与 Y 对易。这是标准源响应方法在既有选态密度上的应用。

对完整 covariant tower，必须换为源 chi=e−W：sigma1(Y)−tau(Y)=√q F_chi''(0;Y)，tau chi²=1+lambda。它保留 −tau((eW+We)Y)+tau((W²−lambda)Y)。该源不等于实际时间生成元，不可将全 bilocal 的一般非平稳选态称为 Gibbs 态。固定 q 时 chi 有界，实际极限双局域有所有固定矩，因此源在零点的导数合法；不声称全 finite-N quenched partition function 的极限。

## 2. 调用已有热鞍点

Lin–Stanford (137)–(138)、(409)–(410)给 theta=pi v/2、beta=2theta/cos theta，以及未对半分裂的实时间核

\[
g_{\beta,d}(u)=\left[\frac{\cos\theta}
{\cosh(u\cos\theta+i\theta)}\right]^{2d}.
\]

在 beta0 两次微分得 (d/2)sech^(2d)u [u tanh u−(1+2d)tanh²u]。d=2 代入精确接口，得到预测

\[
(C_1-C_0)/\epsilon\longrightarrow
R(u)=\operatorname{sech}^4u\,[u\tanh u-5\tanh^2u].
\]

这里还需要证明热渐近的 beta 二阶导数在零点收敛；beta 点态渐近不能授权微分。足够条件是复 beta 邻域的局部一致解析余项控制，可用 Cauchy 公式。**本轮尚未证明这个条件**，故完整固定时间式是半经典预测；精确源恒等式及有限 Taylor 系数则已验证。epsilon→0 与有限 N→∞不交换。

精确短时曲率差为 √q(1+q)(1−q²)²/epsilon=8epsilon+O(epsilon²)，与 R=−4u²+11u⁴+…相符。若只重新定义 u→(1+a epsilon)u，sech⁴u 给修正 −4a u²+(28a/3)u⁴+…；u² 要求a1，u⁴随即不匹配。归一化固定在u0，所以幅度也不能自由吸收此形状。

## 3. 证据与数值范围

新 helper/checker 从 q-Hermite 递推独立生成有理整数多项式 D_(2k)(q)，对 q=exp(−epsilon)精确展开。时间0..14阶共8系数吻合，36个低阶 Laurent 项精确相消。首几项−4、11、−769/45、6236/315。有限检查不当作所有次数证明。

四个 epsilon=.08,.04,.02,.01，时间u=0,.25,…,1.5；使用 Jacobi sparse exponential action，L=214,364,664,1264，各与L+64比较；epsilon .08另用独立对称三对角 eigensolver。选择直接计算选态减迹态的cross packet，避免相近两数减法，但包内仍有相消，最大条件数约413.48。输入参数、版本、尾界、相对误差均在 response-checks.json。

系数误差 max|(C1−C0)/epsilon−R|依次 .077469、.040022、.020335、.010249。epsilon .01 时 u .25仍相差2.68%，u .5相差1.98%；这是真实有限变形偏差，不能用微小cutoff差遮盖。解析Jacobi尾为exact-arithmetic bound，不包含roundoff，浮点对照不是区间认证。

科学输出 SHA256：814039c5cc5e3366db8d6957bfef66e80efee02adf409c7e1df0ac59e587fee0。首次运行通过后封存，不重跑。既有成功记录未改。

下一首任务是针对实际 q-Mehler 热积分证明带两次 beta 导数的固定时间鞍点余项；随后处理完整 e,W 双源与 bilocal。空间/质量生成元、任意非多项式 selected 交错、finite G仍未完成。本轮仅本地保存，不重复上一轮Git或issue回复。
