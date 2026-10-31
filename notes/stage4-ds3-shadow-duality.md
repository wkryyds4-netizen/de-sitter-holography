# Wilson对偶消去完整交叉核的指数误差放大

2026-10-31，Alan Day / Chengheng。接续[完整尾界](stage4-ds3-full-tail.md)。[四页PDF](../output/pdf/de-sitter-ds3-shadow-duality.pdf)、[LaTeX](../manuscript/ds3-shadow-duality/main.tex)、[来源](../literature/stage4-ds3-shadow-duality-sources.md)。内部自查推导，采用既有Wilson恒等式，不作新特殊函数或优先权声明。

## 范围与关键结果

继续ℓ=1、G=0、proper SO⁺(3,1)、BD标量观察者、KG物质、固定正质量标签下界与半周期创建块。外标签(a,q,c,p)，内标签(b,d)，均为正实数；h_i=(1+iν_i)/2，barh_i=h_i，不是复共轭。U00仍为完整Racah变换中的标量矩阵元，不宣称标量子块完备或酉。

上一阶段确认原负shadow项有exp(2πM)的系数放大。本轮从Groenevelt的Wilson函数对偶性和参数置换对称性导出精确反射关系，把指数小因子与大系数先在解析层面合并。原失败和旧PDF全部保留；改动的是求值表示，没有改物理态、谱测度或核。

## 参数字典与反射证明

对H=h或1−h，Q=h′或1−h′，定义Q0=H+h1+h3−1、z=Q0−Q、

\[
A=(H+h_1+h_2-1,H+h_3+h_4-1,H+h_3-h_4,H+h_1-h_2),
\]
\[
T_j=A_j-z/2,\quad \lambda=z/(2i),\quad x=(Q_0-H-z/2)/i.
\]

以F(T;λ,x)表示Groenevelt Proposition 4.4的对称c函数展开，精确系数见PDF式(1)。S=ΣT，D(T)_j=S/2−T_j；已知对偶及置换对称性给F(T;λ,x)=F(DT;x,λ)。匹配两组LPRS 4F3参数及Gamma系数，而非只匹配上参数，得到

\[
\Omega(H,Q)=K(H,Q)F(T;\lambda,x),\qquad
K=\Gamma(2H)\Gamma(2Q)\Gamma(A_1)\Gamma(A_2)\Gamma(A_3-z)\Gamma(A_4-z).
\]

下标s指同时反射(H,Q)→(1−H,1−Q)，外权重不变。在实标量切片上，λs=x、xs=λ、DT=conj(Ts)，且λ、x为实。由F的对偶性、偶性及复共轭性质，

\[
\Omega(H,Q)=R(H,Q)\overline{\Omega(1-H,1-Q)},\quad R=K/\overline{K_s}.
\]

取消共同Gamma因子后

\[
R=\frac{\Gamma(H+h_1+h_2-1)\Gamma(H+h_3+h_4-1)}
{\Gamma(H+1-h_1-h_2)\Gamma(H+1-h_3-h_4)}
\frac{\Gamma(Q+1-h_1-h_4)\Gamma(Q+1-h_2-h_3)}
{\Gamma(Q+h_1+h_4-1)\Gamma(Q+h_2+h_3-1)}.
\]

所有这些Gamma参数实部为1/2，比例有限且非零；可去碰撞点由连续性包含。本证明使用已知亚纯函数恒等式的解析延拓，不套用Wilson酉变换的正测度域，也不移植其谱密度。

## 完整核与尾界

令P=Ω(h,h′)Ω(1−h,h′)，则

\[
U_{00}=C_+P+\widehat C_-\overline P,\quad
\widehat C_-=C_-R(h,1-h')R(1-h,1-h').
\]

Gamma模恒等式给

\[
|R(h,1-h')R(1-h,1-h')|=
\frac{\cosh(\pi(d-a-p)/2)\cosh(\pi(d-q-c)/2)}
{\cosh(\pi(d+a+p)/2)\cosh(\pi(d+q+c)/2)}.
\]

它正好抵消旧C−/C+的模比，所以对全部正实标签，|Ĉ−|=|C+|。程序将系数与Gamma比的对数合并后指数化，保持复相位，不强取实部或对称化。

若两个正shadow因子的有限和为X̂,Ŷ，精确算术尾分别eX,eY，则完整U绝对尾≤(|C+|+|Ĉ−|)(|X̂|eY+|Ŷ|eX+eX eY)。仍不包含有限和/Gamma舍入，不支持在最终相消点直接除以结果作相对认证，也未证明质量积分极限。

## 验证与保留的失败

[初始33项](../results/stage4/ds3-shadow-duality/validation.json)：24个旧四Omega样本、4个一般标签原hyper比较、4个Gamma比例核对、1个独立Euler配对。符号计算另核六个独立实标签下的三个字典恒等式。Euler在35位的相对差约4.75e−19。

初始33项31通过、2失败：M16,r=.1/.4的旧80位负shadow参考相对残差2.04e−60、7.42e−60，高于预设1e−60。原状态failed、原容差和数据均保留。[新精度复核](../results/stage4/ds3-shadow-duality/refinement.json)对两点各跑100和130位原hyper，全部达到预设1e−75；实际残差约1.24e−80、2.18e−110及3.79e−80、6.95e−110。不是放宽阈值后把旧失败改成通过。碰撞prototype偏差也保留。

[18组稳定完整核计算](../results/stage4/ds3-shadow-duality/study.json)：M2,8,16,32,64,128，N64/128，70/110位。固定p=q=.6,u=.3,v=−.2。N128时M8/16/32对旧完整核参考的复绝对差约1.44e−13、2.80e−23、1.66e−42，对应解析尾9.01e−12、1.11e−21、1.13e−40。M32旧表示误差约4.98e43，改表示的改善可直接比较；旧参考本身未区间认证。

M2直接和仍慢，N128尾约1.52e−4，不能宣称全质量域统一高效。M64/128未提供独立旧公式完整参考；M128的解析尾2.99e−127，70/110位变化却为5.35e−68，明确保留算术精度下限。不能将最小尾数值当总误差。

[365项独立审计](../results/stage4/ds3-shadow-duality/audit.json)通过：Fraction重建多项式与平移系数，Decimal重建遗漏项模、尾及完整复残差，核等模系数、原失败与新修复、哈希。检查的容差按保存精度设置，不把残差最后一位舍入当数学矛盾。

## 后续

优先分析剩下的两个正shadow因子及完整相位的联合渐近，区分固定初态、条件化质量、移动下界及时间窗口，然后对真实创建幅作积分和谱尾控制。单点稳定性未证明完整δ核或Bose自由极限。全Fock、有限G和全息微观字典仍开放。

唯一旧活跃句柄仍为47841（100位参数导数验证），未计入本轮通过项目；本轮63979与3172均已完成。collision、full-tail、shadow-duality三个阶段待整体项目登记与冻结。
