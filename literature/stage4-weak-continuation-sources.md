# 弱积分与移动极点：来源和使用范围

2026-11-01。先核已有定理，再决定计算；独立推导不等于原创。

| 来源 | 核读位置 | 本轮使用 | 不由它推出 |
|---|---|---|---|
| Kim–Scarinci, *A quantization of moduli spaces of 3-dimensional gravity*, [2112.13329v2](https://arxiv.org/abs/2112.13329v2) | §6.3、Prop.6.6 的产品与极点；§6.4 延拓 | 既有 Phi、极点除子和产品约定；由产品直接提取首留数 | 本文的字面弱积分就是原文意图的完整延拓；最终有界五边形错误 |
| Thierry Ramond, *Distributions and Partial Differential Equations*, version 2.06, 2015-01-19，[大学原站](https://www.math.u-psud.fr/~ramond/docs/cours/distributionsedp.pdf) | Example 4.4.5、Proposition 4.4.6，印刷 p87 / PDF 第88页；该页已渲染查看 | 标准 ∂bar(1/z)=πδ 及 Cauchy 核的分布论基础 | 新的量子双对数定理或物理观察者态 |

Ramond 本机缓存 SHA256：aed544920c4027dde81c3e536718842d295fab1439d4110b8ec3758e112a2f3f。KS 沿用已固定 v2，全部缓存及前阶段证据的哈希见 results/stage4/weak-continuation/source-record.json。第三方全文和核版图片仅留本机，不上传。

本轮新增工作是明确参数映射、Jacobian、单极点邻域及固定紧测试的匹配；Cauchy 源项是标准理论的应用。径向 Cauchy 变换可直接角积分求出，无须当作新延拓理论。40位完整产品计算只是数值交叉核对。

此前 FG/Kim/Pusz–Sołtan 的适用条件仍以[冻结的证明依赖图](flip-proof-dependency-map.md)为准，不回写。Demailly 在线访问未成功，未作为本轮已核依据；没有声称穷尽全部相关文献。

后续文献筛选必须同时记录：具体命题、版本/页码、算符表示及谱、参数域、测试空间/闭包、轮廓/边界值、与最终有界配对 F 的等同性。只找到相近函数名称或形式五边形不足以关闭问题。
