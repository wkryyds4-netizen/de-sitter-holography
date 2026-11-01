# 有界配对的 Fourier 核：文献适用性

2026-11-01。三份新原文已固定版本下载；以下是局部核读，不代表全文独立验证或穷尽检索。完整缓存不公开，哈希由本阶段 source-record.json 保存。

| 来源 | 核读位置 | 可用与不可直接移植的内容 |
|---|---|---|
| Garoufalidis–Kashaev，[1706.08132v2](https://arxiv.org/abs/1706.08132v2)，2018-10-17 | Appendix A eqs.104–132；Appendix B eqs.133–140；PDF30、31页全页核版 | 有界算符五边形与分布 Fourier 核等价；圆周×整数例子必须保留整数谱。 |
| Garoufalidis–Kashaev，[2306.01331v2](https://arxiv.org/abs/2306.01331v2)，2023-06-05 | §2.1 Example2.4 / Remark2.6；§3 Definition3.1 / Lemma3.2 / Theorem3.3及证明开头；PDF4、5、8页全页核版 | R²实方向比值可匹配KS负宇宙常数配对；local-field构造使用乘法群及其对偶。复数域情形并非加法C上的当前模型。Fourier显式式在特殊点仍须分布解释。 |
| Kashaev，[1403.1298v1](https://arxiv.org/abs/1403.1298v1)，2014-03-06 | §2 Lemma1 / Theorem1及相关证明；PDF4页全页核版 | 已有 beta pentagon 构造；这里明确写有平方可积假设，不能把无限平面上的单位模函数直接代入，也不能给函数加阻尼后假设它仍满足五边形。 |
| Kim–Scarinci，[2112.13329v2](https://arxiv.org/abs/2112.13329v2) | 延续§6.3产品/反演、Prop6.15表示，复读§6.5 | 当前要复现的对象和交换关系；并非因找到一般核框架就完成此特定F的证明。 |

本轮核出完整的符号/测度字典：gamma=exp(−ixy/π)，chi=exp[−i(xy′+yx′)/π]，自对偶测度dxdy/(2π²)。文献的p对应−P，q对应Q；负号Fourier变换的逆共轭把该顺序变为当前Q、P。R²实例的变量为dot=x/(2πsqrt(h))、ddot=−sqrt(h)y/π，绝对Jacobian正好是该测度。

圆周函数 C_q(theta,m) 的连续m标量延伸可以因子化当前F，但e^(im theta)在圆周上单值要求m为整数。个别因子的谱定理不能由标量代入继承。两个因子合并后仍有正确反演，故这不是配对五边形反例。

新推导只应用标准Fourier、Plancherel、热半群和Cauchy–Schwarz，给当前表示下的具体误差式；没有提出新的一般数学定理。后续应优先找直接覆盖该配对谱表示的原文；同时可按稿中的绝对可积截断核做一个认证矩阵元。

检索中错误尝试了不存在的2306.01331v3 HTML，后经arXiv摘要确认最新固定版本为v2。Poppler路径下没有pdftotext.exe，改用本机MiKTeX提取成功，其日志权限警告保留在工具输出；PDF原字节与核版不受影响。页面渲染有Type3字体边框警告，实际相关公式均已目视检查。
