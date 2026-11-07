# 单选态双局域均值：前驱与规范

2026-11-08在线核对；本轮先使用已知公式，再做本项目的源投影与实际滤波匹配。

- W. Riley Casper、Stefan Kolb、Milen Yakimov：[Bivariate continuous q-Hermite polynomials and deformed quantum Serre relations，v2](https://arxiv.org/html/2002.07895v2)，引言递推及式(1.1)。取该文r=q，坐标缩放√(1−q)x/2、√(1−q)y/2，整体乘(1−q)^−(a+b)/2，得到本稿P_ab。用独立显式有限和检查72个多项式；不将原文另一正交权直接识别成我们的Mehler测度，也不声称发明该多项式。
- Edward G. Effros、Mihai Popa：[Feynman Diagrams and Wick Products Associated with q-Fock Spaces](https://arxiv.org/abs/math/0303045)，前期已核Theorems3.2–3.3的whole-Wick禁配对和有序乘法规则。本轮使用既有局域CAR源作为模型接口，不从抽象qFock公式猜测size mark。
- Carl de Boor：[Divided Differences，arXiv v1](https://arxiv.org/html/math/0502036v1)，Genocchi–Hermite式(52)及紧接其后的插值余项。以monic H_j的实根插值exp(isx)，得到 |c_j(s)|≤|s|^j/j!；此为标准定理的直接应用。本文对应Surveys in Approximation Theory 1 (2005),46–69；作者网页的另一修订排版页码不同，引用固定arXiv版本。

项目内部证明输入：[多源局域性](../notes/susskind-bilocal-multisource-proof.md) §§2–6；[高阶selected密度](../notes/susskind-slater-higher-wick.md)；[有序后代](../notes/susskind-memory-wick-descendants.md)；[neutral actual/local尾](../notes/susskind-memory-neutral-marked-series.md)及其LaTeX §8显式尾；[固定电荷actual bilocal](../notes/susskind-sector-bilocal.md)；[第一neutral选态响应](../notes/susskind-neutral-selected-response.md)。新内容是将单L源的任意测试词结构识别为one-U投影，继而投影整条N_j递推及单选态一阶矩。既有微观桥、数学多项式及通用插值工具均保留归属。

识别边界：投影并非代数同态。下一轮两个高阶插入需要保留共同源连接及投影正交补；本文不能用于直接相乘求高点，也没有识别体质量/空间生成元或有限G。
