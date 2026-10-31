# 联合端点估计：来源与使用边界

核读日期：2026-10-30。以下为标准特殊函数公式及本项目已冻结结果的应用；没有确认新颖性，不作首创声明。

|来源|本轮使用|范围|
|---|---|---|
|[DLMF 5.4](https://dlmf.nist.gov/5.4)，5.4.3–4|纯虚及半整数 Gamma 的精确模长|结合递推生成整数 N 的乘积；跨 s=M−ν 仍用精确 cosh，未延续错误远端渐近|
|[DLMF 5.5](https://dlmf.nist.gov/5.5)，5.5.1、5.5.3|Gamma 递推与反射|正尾反射两个因子；中央区只反射一个；周期因子保留|
|[DLMF 5.11](https://dlmf.nist.gov/5.11)|带扇区余项的 Stirling|对增长移位自行 Taylor，检查 K/ν→0 的安全扇区；不将固定移位展开直接称为联合界|
|[DLMF 10.32.9](https://dlmf.nist.gov/10.32.E9)|Kα 的正实积分|构造覆盖 y=0 的 E 下界|
|[DLMF 10.41.3–4](https://dlmf.nist.gov/10.41)|Iα、Kα 的大阶数渐近|仅在固定 y≥y_c>0 使用 E 上界；不将其当作 y→0 的一致证明|
|[上一轮联合正尾](../notes/stage4-recoil-bessel-cancellation.md)|∫正尾≤10(R+1)H_R|本轮补足 H_R/E，旧证明及代码冻结|
|[精确核代码](../code/recoil_bessel_limit.py)|Qν 的定义与归一|只读导入 log_q，不重跑旧计算|

上述 DLMF 页本轮在线核读；5.11 另有既存缓存，SHA-256 由新 diagnostics 绑定。本轮未新增原文缓存，主来源清单仍为 26 份文件。

探索中还检索到 Temme, “Uniform Asymptotic Methods for Integrals” ([arXiv:1308.1547](https://arxiv.org/abs/1308.1547))，目前只读搜索摘要，没有引入其任何定理；不把它写成已解决中央轮廓的依据。

本轮的对象仍是固定质量先 κ→∞ 的内部手征 Barnes 图。这里没有新增物理态/代数或微观模型的识别。下一次查新应围绕保持相位的 Barnes 分支分解和移动极点轮廓，再核可用定理条件。
