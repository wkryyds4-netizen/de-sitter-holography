# 湮灭模式极限的来源与查核范围

2026-10-29。本轮先核原式，再写[湮灭推导](../notes/stage3-integrated-annihilation.md)。没有完成足以支持首创性声明的全领域查重。

1. [Kolchmeyer–Liu 2411.08090v1](https://arxiv.org/html/2411.08090v1)，继续用[本地缓存](../sources/kolchmeyer-liu-2411.08090v1.html)。实际复核(3.14)—(3.15)的共轭方向、(5.4)局部observer矩阵元、(5.6)负号时间演化，以及(5.8)后群体积说明。原式作为模型定义，不把未说明的无限模式交换算作原文已有定理。
2. 固定左包界和扩展GA类使用项目既有[弱下降稿](../notes/stage3-regulated-local-descent.md)L6—L8与[产品种子证明](../notes/stage3-product-seed-ga-density.md)。这一步无需右映射等变，只要求右向量总K=0；完整群等变性仅在最终未截断块下降时使用。
3. [András Vasy, Math 172: Convergence of the Fourier Series](https://math.stanford.edu/~andras/172-3.pdf)，Appendix B，印刷页13—14。在线读到Fejér核的正性、归一化、远离零衰减及连续函数的一致Cesàro收敛。本稿从有限几何和直接重建核，并给Hilbert值推广的短证明，归一化采用dθ/(2π)。此为标准数学输入。浏览工具可读PDF；本地Invoke-WebRequest下载因DNS解析失败，未缓存，不计入26份来源文件。

UPenn讲义链接此次浏览返回403，因此未作为已阅读证据。DLMF 1.8未提供所需Fejér论证，不引用它支持该命题。

本轮真正需要检查的是：滤波只作用新发射腿；有限M保持K但不保持boost；使用固定M绝对Fubini后才作GA识别；用固定左包界传递强收敛；点值Fejér结论只在物理种子域给出。均写入LaTeX和中文推导；没有新增数值扫描或独立审稿。
