# 相消与联合尾界：来源说明

2026-10-30。先复读[完整固定K来源](stage4-recoil-bessel-limit-sources.md)及实际推导，再在线核查以下原始公式。

|来源|定位/用途|
|---|---|
|[NIST DLMF 5.11](https://dlmf.nist.gov/5.11)|Gamma比值与Bernoulli多项式展开，得到L₂、L₃；已有[本地缓存](../sources/dlmf-5.11.html)|
|[NIST DLMF 10.28](https://dlmf.nist.gov/10.28)|10.28.2的Wronskian，证明Kα/Iα严格递减及J零点简单|
|[NIST DLMF 10.30](https://dlmf.nist.gov/10.30)|固定阶数、小自变量端点，保证比值从∞开始|
|[NIST DLMF 10.32](https://dlmf.nist.gov/10.32)|正实K、非负α下I和K积分表达的正性|
|[NIST DLMF 10.40](https://dlmf.nist.gov/10.40)|固定阶数、大自变量端点，保证比值趋0|
|[NIST DLMF 10.41](https://dlmf.nist.gov/10.41)|§ii的大阶与导数渐近，用于零点聚集及已知二阶系数的K³尺度|

两支连接、负轴下方继续、完整固定K轮廓沿用上轮已核查内容。DLMF公式不直接提供原Barnes积分的联合误差界。本轮用精确递推补充O(K)正尾过渡，用固定参数的积分余项补充零点抬升；两者不能合并冒称完整联合证明。

本轮没有新增缓存，manifest仍26份。零点分类是针对所选辅助复和的标准特殊函数推论，不代表物理相变或新的全息模型，也没有完成全面优先权检索。
