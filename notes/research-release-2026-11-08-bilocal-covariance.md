# 集中研究快照：双局域均值与制备协方差

本次依用户要求，将两个相连阶段一次性提交并推送。日常低频 Git 约定继续；下方材料中的“本地”表示成稿时的研究状态，不是当前仓库可见性。

| 阶段 | 阅读稿 | 证据 |
|---|---|---|
| 单插入投影与完整单选态均值 | [七页PDF](../output/pdf/de-sitter-susskind-selected-bilocal-mean.pdf)、[推导](susskind-selected-bilocal-mean.md) | [清单](../results/susskind-selected-bilocal-mean/local-evidence.json) |
| 完整两点与制备协方差 | [八页PDF](../output/pdf/de-sitter-susskind-selected-bilocal-covariance.pdf)、[推导](susskind-selected-bilocal-covariance.md) | [清单](../results/susskind-selected-bilocal-covariance/local-evidence.json) |

第一个阶段从已有局域源和双变量 q-Hermite 结构得到完整双局域的选态均值。第二个阶段保留同一原始耦合的全部共同涨落，给能量倾斜谱项加制备修正的两点公式。连通制备修正可以为负；共同中心时间衰减发生在固定 lambda 双标度之后。空间/质量字典、四双局域源分区及有限 G 尚未完成。

均值阶段的320项投影、72多项式、10时间系数和12物理点，以及协方差阶段386个选态比较、独立边界/Gram检查和30物理点均已有保存输出。收尾只读核对封存哈希，未重跑成功科学计算。新协方差清单核15个当前及417个父科学文件，八页最终稿双遍编译并逐页渲染检查；七页均值稿沿用封存版本。

数值协方差只涵盖制备项及基线消去后的完整方差差，绝对完整谱基线没有做数值积分。尾界控制degree截断，不是求积与舍入的总误差认证。

issue #1 的[既有核查](issue-1-review.md)保持：所提供区间的独立算术重组支持相应保守下界，但全21,143单元独立重积分及9,890旧端点差异定位未完成。本次没有新增KS证据，issue继续开放。LaTeX、代码、科学输出及公开研究入口属于本次发布范围；私人日志和研究动机仍仅本地保存。
