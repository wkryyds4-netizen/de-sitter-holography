# 双标度singlet谱矩：原文使用记录

核读日期2026-11-05。有限生成式与两项谱矩是本项目对指定算符的推导；其数学工具属于已有Wick、行列式与随机交叠方法，不作原创优先权判断。

| 来源 | 本轮实际阅读 | 在本轮中的位置 |
|---|---|---|
| [Kota–Vyas 1504.00759v1](https://arxiv.org/html/1504.00759v1) | §7定义(101)–(105)、§7.1(106)–(113)、§7.2(114)–(118)及极限说明 | 现成粒子移除矩；(113)仍包含耦合约化矩阵元。§7.2是固定body rank稀薄极限，不直接当本项目双标度。没有逐项计算其Racah系数 |
| [Kota nucl-th/0410069v1](https://arxiv.org/html/nucl-th/0410069v1) | §IV(23)–(31)，承接前轮§II–III | 已有协方差特征值及四阶矩。作为body-rank方法依据，而不把已有谱定理称作新结果 |
| [Ferialdi–Diósi 2110.02920v2](https://arxiv.org/html/2110.02920v2) | §I–III相关定义，§V(43)–(49)及其证明衔接 | 费米排序符号和收缩差；本轮用来组织正常序内部收缩标记 |
| [Klich cond-mat/0209642v1](https://arxiv.org/html/cond-mat/0209642v1) | §2(8)–(12)及有限维说明 | 核查Fock迹到单粒子行列式的既有方法。没有把其二次Hamiltonian指数公式直接应用于相互作用SYK演化 |
| [Berkooz–Narovlansky–Raj 2006.13983v2](https://arxiv.org/html/2006.13983v2) | 本轮重核摘要页作者/版本；具体chord规范沿用已读的§4.1及前阶段dictionary | 无标记极限q=e⁻⁴λ的独立规范核对，不重新计算其通用Hamiltonian谱 |
| [Miyashita–Sekino–Susskind 2506.18054v3](https://arxiv.org/html/2506.18054v3) | 沿用已读§§4–6.1的主线；本轮未新增全文阅读 | singlet谱与体传播的物理动机 |

检索过“complex SYK singlet moments”“double scaled derivative singlet spectrum”“SYK commutator particle removal”；未从本轮结果找到直接给当前A₂完整谱的公式。检索结果混有无关/二手页面，未作为推导来源；这不构成穷尽性的新颖性审查。

实际新增的桥梁是：A=(2p−D)H²作为正常序收缩计数的导数；Gaussian配对后的局部子式生成函数；以指数矩控制的Poisson交叠极限。证明和符号在LaTeX附录，不从文献摘要借用未经核对的系数。所有参数均保持本项目每侧p、总阶2p、标准CAR系数方差w和string时间约定。
