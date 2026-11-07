# 下一步seed：charged transfer的精确局部表达

2026-11-07，在wick-descendants阶段收尾时得到；作为下一轮入口保存。以下有限恒等式由定义直接推出并经独立数学复核，尚未运行专门新driver；不将未来的全时间强极限或charged谱写成完成。

同一原J的neutral full-Fock lift X，D X=QX−Σ_i c_i†Xc_i。m>0扇区上

\[
(DX)_m=mX_m-\sum_i c_i^\dagger X_{m-1}c_i.
\]

因此对T_m X=Σc_i†Xc_i/m和F_m(s)=T_m e^{isH_(m−1)}有**精确**关系

\[
F_m(s)=e^{ish_m}-\frac1m[D(e^{isH})]_m.
\]

Gaussian OU和sector restriction可交换。中央m≥p，α_m=1−p/m，Ncal为coupling Gaussian number，令

\[
E_m(s)=m[F_m(s)-\alpha_m^{\mathcal N}e^{ish_m}]
=\left[m(1-\alpha_m^{\mathcal N})-D\right]e^{isH}\big|_m.
\]

这不是把相邻charge Hamiltonian独立化；使用同J full lift正是恒等式成立的条件。m<p时Hamiltonian在相关低charge上为零、F=I，E=0；m=0按零定义。该表达自动包含原forward/reverse源，不能把其中的D改成rank R。

对整数k≥0、u=p/m∈[0,1]，Taylor积分余项给全k界

\[
0\le m[1-(1-u)^k]-kp+\binom{k}{2}\frac{p^2}{m}
\le\binom{k}{3}\frac{p^3}{m^2}.
\]

先保留λ_N=p²/N。在m=N/2+xp、x紧集上，固定k的自然候选为

\[
E_m(s)\simeq
\left[p\mathcal N-D-\lambda_N\mathcal N(\mathcal N-1)\right]
e^{isH}\big|_m.
\]

两个表达之差是纯Ncal乘子，D精确抵消；系数误差≤const_X(k+1)³/p。下一轮可利用已有Gaussian Sobolev/number矩对actual exponential及固定时间导数的界，真正证明全时间与能量权重下的替换，而无需从无导数误差裸求导。这里还没有完成这一全charge估计。

力核接口保持精确。令

\[
\mathscr Q_\alpha f=-f''+i\alpha\{h,f'\}+\alpha^2hfh.
\]

由既有Λ定义，m(Λ−Λ_OU)=Qscr_α E，因此

\[
A_m(s)=(1-P_h)e^{-ihs/2}(\mathscr Q_{\alpha_m}E_m)(s)e^{-ihs/2}.
\]

实际要研究的局部mark是kp−D，和neutral的half-length mark不同。可从D的CAR局部表构建一/二oriented mark的共同J图，再与本轮neutral R_j一起求AA/AB/BA。不能直接把(pNcal−D)e改成Σ_j c_j(jp−D)Z_j；D作用于Jacobi低chaos误差可能留下有限量，必须保留exact E或重新做受控重排。

最先检查：完整finite CAR的T/full-D恒等式与所有chaos数滤波；全k余界；实际Qscr加权替换；随后核A0→Y−、A1及非零A0–B2谱。这一seed不替代任何已有charged transfer/weighted bound，也不需要重新运行旧main。

从D H=pH还可立即核E(0)=E′(0)=0，所以A(0)=−(1−P_h)E″(0)。这是下一轮最直接的符号/归一校准；需要从实际二次词恢复已知Y−，不要只凭方向名称宣布匹配。
