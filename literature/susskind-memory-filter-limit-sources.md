# 六阶矩与滤波极限的来源

2026-11-06实际打开核对：

- [Kota–Vyas, arXiv:1504.00759](https://arxiv.org/pdf/1504.00759)，§2，PDF第7页，(7)–(10)：有限N covariance的两种分解及Racah换道。使用其框架，不使用稀薄极限。
- [Koekoek–Swarttouw, report 98-17](https://tudelft.roelofkoekoek.nl/documents/as98.pdf)，§1.2，印刷页26／PDF第28页，(1.2.1)：Racah的终止4F3定义。实际定位到该公式，未将前一节Wilson公式当Racah证据。

已有covariance composition与本地rank_product_moments.py给Eπχ=χχ及前三矩。本轮沿同一有限三角递推扩至六阶；Touchard多项式是标准数学对象，不作为原创公式。参数多项式消幂＋解析有理分母给的是一致界，有限采样仅作核对。文献核查没有建立原创优先权。

依赖链：

- [一般双电荷公式](../notes/susskind-memory-two-charge.md)：g_A/g_H有限式、H部分L²极限。
- [rank乘积权重](../notes/susskind-memory-rank-blocks.md)：正π、covariance composition。
- [半填充范数](../notes/susskind-memory-residual-norm.md)：A输入缺额多项式加权范数、四方向Gram、Vnew。
- [奇偶Gram](../notes/susskind-memory-parity-gram.md)：S_N定义、JhF桥接。
- [quartic memory](../notes/susskind-memory-quartic.md)：Ecal、raw Z2及W4后代。

当前贡献是将这些已知或已存推导接成g_A的强滤波极限与完整W2core。尚未把full-Fock中心charge窗、时间谱或体字典补全。
