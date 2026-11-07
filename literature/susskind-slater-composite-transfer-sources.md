# 开弦递推的文献接口

2026-11-08，先核原文，再从项目已封存的 Slater composite 源公式推导。

- Berkooz, Isachenkov, Narovlansky, Torrents, [Towards a full solution of the large N double-scaled SYK model](https://arxiv.org/html/1811.02584v2)，§2.3，式(2.12)--(2.14)：按切口上的开放弦数累计部分图；闭合一条弦时计入越过其他开放弦的权重。这是本轮传递方法的既有前驱，不作为新方法宣称。本文的同J whole-D复合标记和Slater内嵌F弦需要额外状态。
- Rajgadia, Xu, [Emergent States and Algebras from the Double-Scaling limit of Pure States in SYK](https://arxiv.org/html/2604.14387v3)，§3.2、§4：原文式(69)--(78)、(104)将state-adapted creation/annihilation与dressing接入切图/GNS。其Majorana/KM模型中的独立matter与本文charge-conserving、同原J的曲率不相同；不移用其具体交换关系。
- Josuat-Vergès, Kim, [Touchard--Riordan formulas, T-fractions, and Jacobi's triple product identity](https://arxiv.org/html/1101.5608v1)，引言式(1)：经典Touchard--Riordan有限和用作H^24的独立参照。初始FPSAC与HAL PDF读取超时，随后成功读取arXiv HTML，核对式(1)的二项式差、(-1)^k、q^(k(k+1)/2)和(1-q)^(-n)。该公式并非本项目推导的新成果。

项目内主输入是 code/slater_composite_sources.py 及 manuscript/susskind-slater-composite/main.tex。旧 full-Fock CAR、固定Slater defect支配、central pinching和polynomial-density不重算。本轮新增接口是pair-collision源系数的完整几何分类，以及只记开放端点/两类link的递推。概率/正性只在物理关系 gamma=exp(-lambda)、固定lambda>0上从旧微观极限继承。独立有理点仅检查代数恒等式。
