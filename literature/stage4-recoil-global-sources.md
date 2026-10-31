# 全局反冲与Bessel辅助模型：来源和证据边界

查阅日期：2026-10-30。此表记录本轮直接使用的原始公式；不把标准Bessel连接公式改写成原创结果。

|来源|定位|本轮用途|
|---|---|---|
|[NIST DLMF 5.5](https://dlmf.nist.gov/5.5)|Gamma反射及递推|负实方向精确核分解|
|[NIST DLMF 5.11](https://dlmf.nist.gov/5.11)|Stirling及Gamma比值展开|固定参数尾部、紧集核极限与首阶修正|
|[NIST DLMF 10.25](https://dlmf.nist.gov/10.25)|I函数幂级数与微分方程|两族留数和、辅助微分身份|
|[NIST DLMF 10.27](https://dlmf.nist.gov/10.27)|I与K连接式|辅助积分的两分支表达|
|[NIST DLMF 10.34](https://dlmf.nist.gov/10.34)|10.34.2，m=−1|负实轴下方解析继续的相位与符号|
|[NIST DLMF 10.41](https://dlmf.nist.gov/10.41)|§ii，10.41.3–4|固定y>0的K与阶数同步增长|

5.11已有[本地缓存](../sources/dlmf-5.11.html)，其余本轮在线读取；没有声称新增离线缓存，sources/manifest仍为26份。大阶渐近不是针对原Barnes核的质量一致尾界。

上游物理图、CLP和CSTY原文关系沿用[内部积分来源](stage4-recoil-internal-sources.md)、[反射鞍点来源](stage4-recoil-saddle-sources.md)、[局部作用量来源](stage4-recoil-local-sources.md)。本轮没有完成新的全领域优先权检索，也没有由辅助函数证明任何新的微观全息对应。

证据层次：负反射与辅助Bessel积分为解析恒等式；Gamma展开为紧集渐近；完整Barnes积分为有限参数数值证据；完整极限交换、联合大K、物理四点函数和有限G均未证明。下一步须补质量一致尾界。
