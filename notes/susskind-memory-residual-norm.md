# 半填充的 A 范数与保留下来的记忆方向

2026-11-06，本地。[六页PDF](../output/pdf/de-sitter-susskind-memory-residual-norm.pdf)、[LaTeX](../manuscript/susskind-memory-residual-norm/main.tex)、[来源](../literature/susskind-memory-residual-norm-sources.md)。不commit/push/ZIP。

这次完成上一轮缺失的固定半填充A范数。把H²、A、JhF三个已有方向都加入投影后，输出缺额仍有严格正的残量；smallλ范数平方趋16Js4。证明直接在同一半填充态进行，没有把full-Fock的方差代入。

## 一般电荷的交换 Casimir 恒等式

标准CAR、完整重叠complex Gaussian系综、normalized trace、annealed HS。Cx=Σij[Eij,[Eji,·]]/2，rank r特征值x_r=r(N+1−r)。独立单位方差Hermitian rank-a/b张量U,V，Phi_a(X)=EUXU，Phi_a在rank b的特征值chi_ab。

任意合法charge M，Cf=M(N−M)(N+1)/N：

\[
E\operatorname{tr}(UV)^\dagger C_x(UV)=x_a+x_b-x_ax_b/C_f,
\]
\[
E\operatorname{tr}UV C_x(UV)=(x_a+x_b)\chi_{ab}^{(M)}.
\]

第二式证明很短：Cx(UV)=(xa+xb)UV+Σ[E,U][F,V]。额外项在U平均后为tr[Phi_a(VE)−Phi_a(V)E][F,V]。Phi_a对trace pairing自伴，又Phi_a([F,V])=chi[F,V]，两项逐一抵消。这不是半填充奇偶规则的外推，对一般charge也成立。

## 保留相干的有限 N 范数

回到N=2m，n=m+1，p<m。h=H_m、k=H_(m−1)用同一耦合，A=Rem(k²)+(2p−m)h²，O=W2A，C=Eh²，Ca=Eh_a²，omega_a=Ca/C。Rem rank-r范数增益ell_r²=m(mn−x_r)/n。

Ca−/Ca=R_a=(m−p)²n/[m(mn−xa)]，Cminus_f=(m−1)n(N+1)/N。下charge的removed平方，按每个ordered input pair、除CaCb的norm是

\[
R_aR_b\frac mn[(mn-x_a-x_b)(1+\chi^-_{ab})+x_ax_b/C^-_f].
\]

chi为上charge核，xi上一轮coherent核，D=xi−chi。逐项Pochhammer给chi−=chi+D/n。定义g=(m−p)²/m，h0=2p−m，eta=p²/m，P=g² xaxb/[(mn−xa)(mn−xb)]，精确

\[
\|O\|_m^2/C^2=\sum_{a,b}\omega_a\omega_b f_{ab},
\]
\[
f_{ab}=\eta^2(1+\chi)+P(mn/C^-_f-1-\chi)
+[(g^2-P)/n+2h_0g](\xi-\chi).
\]

每个f都是非负的二次Gaussian系数范数，包括a=b decoupling。单独几项会增长，不能先各自取极限。

## p²/p 大项相消与极限控制

a=p−u，b=p−v，m=p²/(2λ)，chi第j项写c_j，c0_j=(-4λ)^j/j!。

log(c_j/c0_j)=L1/p+L2/p²+O(p^-3)，

L1=−j(u+v+2λ)−j(j−1)，

L2=−j(u²+v²)/2−j(j−1)(u+v)/2−j(j−1)(2j−1)/6+λj(u+v+1)+3λj(j−1)/2−λ²j。

令K=(g²−P)/n+2h0g，c_j乘的系数eta²−P+Kj/(m−j)=p²r2+pr1+r0+O1/p：

- r2=−(j+4λ)/λ；
- r1=4(2j+2λ+u+v)；
- r0=−2j²+(1−20λ)j+4λ(u+v+2)−4uv。

对c0_j求和，用Euler算符z∂z作用e^z，得到三阶系数(0,0,4λ²q)。非chi项eta²+P(mn/Cminus_f−1)→4λ(1+λ)，q=e^-4λ。因此固定deficits的f→4λ+4λ²(1+q)，与u/v无关。

余项不是只在固定j使用：取Kp=ceil(log p)。所有a,b≤p的|c_j|≤Aλ^j/j!，外系数≤Bλ p²(1+j)；j>Kp的绝对尾压过任意p幂。u,v,j≤Kp内logs的Taylor余项≤poly(log p)/p³，乘p²后趋0。input窗外用f≤Bλp4及omega_(p−u)≤Aλ^u/u!，总贡献也趋0。故完整rank和有同一极限。

该证明同时允许任意固定多项式input-deficit权重，解决了A的两个输入缺额的范数尾。**输出缺额S=a+b−r的A范数尾尚未解决。**

## 四方向 Gram 与正残量

theta=2λ，z=1+q，b=theta(1−q)，s=theta²z+theta(1−q)，a=theta²z+2theta。结合前两轮，排序(B,Linput,S_N,O)：

\[
G/C^2\to\begin{pmatrix}z&0&b&b\\0&\theta z&0&0\\b&0&s&s\\b&0&s&a\end{pmatrix}.
\]

v=s−b²/z=theta(1−q)+4theta²q/z，d=a−s=theta z；det=theta²z³v>0。由真实S–A cross得到

\[
\Delta_*:=\lim\|O-S_N\|_m^2=2J_s^4(1+q)/\lambda>0.
\]

令c=b/z，Y=O−cB，Sperp=S_N−cB，rho=v/(v+d)。T_N=Sperp−rhoY在极限中正交于B、O、Linput。旧JhF−4λqC Linput→0的L²桥接现在可以与有界O范数配对，所以也完成相对于(B,O,JhF)的投影。

\[
\lim\|T_N\|^2=
\frac{J_s^4(1+q)(1-q^2+8\lambda q)}{\lambda(1+q+4\lambda q)}>0,
\]

先DS再λ→0为16Js4。这里证明的是明确三方向span之外的独立等时方向，不是完整commutant投影或新体粒子。Y方差恰与旧full-Fock表达式相同，但本轮是独立fixed-sector推导。

## 验证与下一步

4JSON首次运行均通过：exchange-casimir为15rank对、68544全ordered Gaussian坐标对、45恒等式；finite-norm为384contiguous、1200 raw/combined公式及旧完整CAR范数28/81/2010；symbolic-limit为log展开、p²/p相消、Gram/Schur及smallλ16；limit-checks为3旧CAR的27项Gram、2额外exact norms、12组80/120位完整rank和。

数值后两项使用X_N/N代理，不是有限S_N；旧L²代理误差加本轮O有界即可传递极限。λ.5,p128的(A norm、O−proxy norm、proxy投影残量)/C²偏差约(-.0310954,−.0126200,−.0261320)，保存抵消因子和实际finite-p偏差，不当区间证书。

六页LaTeX双遍第一次编译通过，无undefined/Overfull；95dpi六页逐一view_image核对通过。全部计算终态，未覆盖旧科学结果，未提交推送Git。

下一步重点缩为A的输出rank尾，随后证明g_A在实际算符和里可换固定deficit常数；再扩中央charge窗、补W4低chaos/rawZ2及Y/L联合谱，最后finite-memory和体空间质量。不要重算已完成的AA/SA等时Gram。
