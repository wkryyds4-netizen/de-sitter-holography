# 四段质量滤波与条件关联函数：来源记录

2026-10-31。主稿[LaTeX](../manuscript/mass-filter/main.tex)。

1. [CSTY v1](https://arxiv.org/html/2607.14042v1)：本轮重新核§2.1/(2.3)固定四段平均质量变换、§3.1全部正性论证及质量Gaussian滤波、§3.2轮廓问题、附录F。原文已经用四个sδ(mi)=exp[−δ²(mi−M)²]处理从平均质量到迹的比较，并给δ∼M^−1/2及O(1/M)估计。不能声称原文忘记质量制备，不能用平均质量不是态这一点推翻原文结论。本轮只把原文O(δ²ω²)写成精确cycle逆Laplacian，并算其对真实principal核的影响。
2. [CLP v1](https://arxiv.org/html/2610.30247v1)：复核(3.5)、(3.10)、(3.35)–(3.40)。Zobs是固定拓扑提案；物理两点含额外正对称N，(3.37)仅给领先HHL归一。固定平均能量的两侧反极点关联已由原文说明。本轮质量平均数值把半经典α–μ关系额外设为有限κ关系，只定义比较模型，不把它冒充原文已经给定的精确物理字典。
3. 本项目的[正谱矩](../notes/stage4-principal-eikonal.md)及[旧归一不唯一性](../notes/stage4-recoil-spectral.md)保持原字节。Gaussian平均Bn的正谱表示是本轮复用输入；nλ反例不是本轮新发现。连续G=0制备还须区分真实质量sqrt(1+a²)与principal标签a。

有针对性的检索未返回额外相关前驱，结果不用于任何优先权判断。最接近的前驱就是CSTY§3.1本身。精确协方差是标准Gaussian Fourier线性代数，小宽度斜率是标准Bessel方程/Plancherel的直接推论。

输入缓存及前阶段证据的SHA由本轮source-record绑定；不更新旧manifest，不公开第三方全文。
