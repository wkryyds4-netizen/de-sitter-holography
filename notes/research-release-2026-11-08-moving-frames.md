# 集中研究快照：制备、热窗口与物质弦移动框架

本次依用户要求，将四个相连阶段集中提交。此前文稿中的“本地”“未 Git”保留成稿时的含义；日常低频 Git 约定继续。

| 阶段 | 阅读稿 | 封存证据 |
|---|---|---|
| 制备态的模谱与非 KMS 能量窗口 | [六页 PDF](../output/pdf/de-sitter-susskind-preparation-modular-clock.pdf)、[推导](susskind-preparation-modular-clock.md) | [清单](../results/susskind-preparation-modular-clock/local-evidence.json) |
| 热弦数窗口与 boost 单次作用 | [六页 PDF](../output/pdf/de-sitter-susskind-thermal-chord-window.pdf)、[推导](susskind-thermal-chord-window.md) | [清单及精度更正](../results/susskind-thermal-chord-window/local-evidence.json) |
| 全单物质模块的移动框架 | [七页 PDF](../output/pdf/de-sitter-susskind-chord-module-frame.pdf)、[推导](susskind-chord-module-frame.md) | [清单](../results/susskind-chord-module-frame/local-evidence.json) |
| 共同背景中的两物质加法 | [六页 PDF](../output/pdf/de-sitter-susskind-multimatter-moving-frame.pdf)、[推导](susskind-multimatter-moving-frame.md) | [清单](../results/susskind-multimatter-moving-frame/local-evidence.json)、[旧两物质分支补核](../results/susskind-multimatter-moving-frame/prior-two-matter-evidence.json) |

这组工作把制备态、原 Hamiltonian 时间与物质弦生成元接到同一条研究线上。标准模理论应用于已有制备密度后，能量窗口表明其 charged 限制态不是原 Hamiltonian 任意常数温度的 KMS 态。热制备阶段建立充分窗口，再通过共同背景位移研究完整单物质和有序两物质模块。

规范为 q=exp(−epsilon)、epsilon=4lambda、K=H/(2Js)=h/√epsilon、u=2Js t。先取固定 lambda、半填充的双标度极限，再取 epsilon→0、m→∞。固定内部核心上的有限生成元词采用 epsilon m→0；与真实热端点的联系使用更强的 epsilon m²→0。两物质结果在实际 ordered direct 内积上成立，给 LS 生成元的物质因子加法以及正、强对易的个体零方向候选。量子群分解及交换权重明确归于已有文献，见[来源核查](../literature/susskind-multimatter-moving-frame-sources.md)。

四个阶段已有两遍 LaTeX 编译和全页渲染检查记录。本次发布只读核对各清单引用的 **537 个不同文件，零哈希不符**，并读取 PDF 页数为 6、6、7、6；没有重跑封存科学计算或覆盖 PDF。热窗口的精度说明和下溢尾界修正须与原输出连读。单物质阶段极小对易子残差的相消限制也保留在清单中，不把截断差当作总误差认证。

下一步对齐同种物质的 crossed Gram 与已知交换权重，保留物理范数随 epsilon 消失的奇通道，计算其对实际四点函数首阶项的贡献。有序张量加法尚未完成这一交换内积问题，也没有确定质量、另一零方向、完整 selected 四双局域排序或有限 G 动力学。

Issue #1 的[既有核查](issue-1-review.md)独立保留：全 21,143 单元的本机独立重积分和 9,890 项旧端点差异的定位仍未完成。本次没有新增该 KS 支线的计算证据，issue 保持开放。私人日志、研究动机及来源缓存不属于发布材料。
