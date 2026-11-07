# 制备字典：模谱与能量窗口详细平衡

2026-11-08。主稿[LaTeX](../manuscript/susskind-preparation-modular-clock/main.tex)、[六页PDF](../output/pdf/de-sitter-susskind-preparation-modular-clock.pdf)。本轮只本地保存。对象仍是半填充、同一原始耦合、先精确pinching后固定lambda双标度的协变代数；不从极限密度倒推有限N模流收敛。

## 已有输入与本轮推进

rho=1−√q(1+lambda)+√q(h−W)²、h/U正交旋转及||[h,rho]||₂²的正值，已在slater-selected-state证明，本轮不再计作新结果。这里用标准模理论与Mehler谱，把该制备的模响应算成显式积分，并提出实际charged能窗探针检验。BKS正交second quantization、Witten模理论及Rindler楔结论均保留归属。

原文回读表明：MSS以原Hamiltonian产生静态补丁boost，微观Boltzmann温度与体GH/Rindler温度分开。我们现有selected state是指定制备，不能仅由其密度对数给原boost换名。已有均值及协方差的非平稳性正是这个制备的数据。

## 精确可计算的模流

令v=(h−√lambda U)/√(1+lambda)、z=(√lambda h+U)/√(1+lambda)，仍是正交单位qGaussian。a=1−√q(1+lambda)>0，b=√q(1+lambda)，rho=p(v)，p(x)=a+bx²，R=2/√(1−q)。在L²(A,tau)的标准表示中，制备向量为sqrt(rho)，Delta=Lrho Rrho^-1；取正温度Heisenberg方向beta_s=Ad rho^(-is)，K=−logrho。它是标准Tomita模自同构的反向参数，sigma对beta_s为逆温1 KMS。

这给完整word制备形式 sigma(A)=a tau(A)+b tau(vAv)=<sqrt(rho),A sqrt(rho)>。标准表示的生成元是LK−RK，不能用单侧LK代替。物理alpha_T=Ad exp(iTh)、T=√C tphys，string u=2√lambda T，既有候选boost eta=2u；s是另一无量纲参数。

one-z谱测度dnu=K_q(x,y)dmu_q²，于是

    Wcal(s)=sigma(beta_s(z)z)=∫p(x)exp{is[k(x)−k(y)]}dnu,
    k(x)=−log p(x), omega=k(y)−k(x).

谱为正，支撑在±Omega，Omega=log[(a+bR²)/a]；交换x,y给dS(−omega)=exp(−omega)dS(omega)，Wcal(s−i)=Wcal(−s)。固定q下K有界，核为整函数。均值频率iWcal'(0)=½∫[p(x)−p(y)]log[p(x)/p(y)]dnu>0，二阶矩为∫p(x)log²[p(x)/p(y)]dnu。相对迹熵为∫p logp dmu。

Hamiltonian在模流下的关联为

    sigma(beta_s(h)h)=[1+b(1+q)]/(1+lambda)+lambda Wcal(s)/(1+lambda).

由于k在正负半轴分别严格单调，pushforward两变量密度仍绝对连续，再取能差得到无原子L¹频谱；Riemann–Lebesgue给Wcal(s)→0，因此上式具有明确平台。此为DS后的模时间结论，不是物理Hamiltonian热化或有限N递归结论。

## 平稳charged子代数上的能窗检验

旧限制态为tau_f，f(h)=1−√q+√q h²，X是normalized charged curvature，crossing q²。它的未滤波X自关联实偶，因此单看该两点的对称频谱不能认证全代数是迹态。

取A_IJ=P_I(h) X P_J(h)，forward= tau_f(A†A)=∫_IxJ f(y)K_(q²)dmu²，reverse=tau_f(AA†)=∫_IxJ f(x)K_(q²)dmu²。这些是正探针范数/跃迁权重，没有引入新散射机制。

对I±=±[.8,1.2]、J=[−.2,.2]，两符号的reverse/forward相等且均≥[1−√q+.64√q]/[1−√q+.04√q]>1。若是原h-flow的常数beta KMS，则reverse=∫e^(−beta(x−y))f(y)dnu_X；正窗口迫使beta<0，负窗口迫使beta>0，beta0又要求比值1。因此不存在任何常数实beta。证明对全部0<q<1成立，不依赖数值。其自身模Hamiltonian仍是−logf(h)，与完整rho的模流须区分。

## 执行与误差

新checker首次成功，22个复时间点用Jacobi与角度乘积两独立表示，各128/192阶；q=.25,.5取物理lambda=−logq/4。最大两表示差2.3315e−14，最大单方法加阶1.8714e−12，相对两表示差2.7542e−14，相消因子至72.914。保留全部复数、sum rules、KMS符号检查。另有正负4个窗口的32/64阶结果：ratio约1.97039949及3.28486964，分别超过解析下界1.57692308、2.32096452。

稿式(17)独立控制q乘积截断的对数尾；192项下两q分别为8.68e−116、1.68e−57量级。它不包含求积、Jacobi截断及roundoff，不能用小product尾替代总误差。KMS交换检查是实现符号检查，不冒充独立理论验证。结果SHAeef7dc91d4280ee0d3955aea75de20a536a243931f2916f3fb941d103da65cc8；环境及输入哈希见原JSON。成功driver不重跑。

## 接续

保留原Hamiltonian为候选boost，把sqrt(rho)制备带入空间生成元矩阵元。下一步具体查Lin–Stanford有限q弦态如何经热Euclidean制备进入大弦数/平滑态窗口。旧beta0零弦种子不能直接认同其sl2生成元，不能再次重复。候选是先有小非零热调节，再研究lambda→0与调节撤除的次序；需由原文鞍点和弦数分布定出窗口，不能先假定它合法。最终检验mixed boost/null-translation Ward identity及反向光锥质量矩，再与已有体介子谱比较。空间生成元、四bilocal及finiteG仍开放。
