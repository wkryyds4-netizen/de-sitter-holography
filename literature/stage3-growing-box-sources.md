# 增长参数估计：核读原文与使用范围

2026-10-30。先核现有Euler端点与组装公式，再建立增长盒估计。本轮为专题续研，未作足以支持首创性声明的系统检索。

- [DLMF 5.11](https://dlmf.nist.gov/5.11)，尤其5.11(i)—(iii)。核Stirling、digamma在避开负实轴的扇区渐近及余项；[本地既有缓存](../sources/dlmf-5.11.html)。本轮不把固定a,b的Gamma比值渐近直接用于增长a,b，而沿平移线积分digamma，单独推出O((1+Q0)²/M)的对数误差。扇区、平移大小与M下界均明确保留。
- [DLMF 15.6.1](https://dlmf.nist.gov/15.6)：本轮在线核Re C>Re B>0及正实路径条件。DLMF写的是正则化Gauss函数；普通₂F₁需恢复Γ(C)。本轮导数上界明确保留该因子。仍仅在线读取，未新增缓存成功声明。
- [DLMF 15.8.4](https://dlmf.nist.gov/15.8)：核连接公式与主值约定；[原缓存](../sources/dlmf-15.8.html)不变。分离δ零点附近的系数极点，再用固定半径圆盘恢复组合的全纯性；不能分别将发散系数当成最终振幅。
- [DLMF 2.3(i)](https://dlmf.nist.gov/2.3)：Fourier端点和分部积分的导数可积条件，本轮在线重读。增长参数下使用固定cutoff和固定导数阶数，显式跟踪参数多项式与Gamma指数；不从一个固定参数O记号推出增长盒一致性。
- KL [arXiv:2411.08090v1](https://arxiv.org/html/2411.08090v1)：[缓存](../sources/kolchmeyer-liu-2411.08090v1.html)及[第11节映射](stage3-assembly-sources.md)沿用，B.21—B.24四项与附录A.14完整外腿必须共同保留。理论仍是严格G_N→0；本轮未扩大到有限G。

本轮推导的结果与限定见[第17节](../manuscript/growingbox.tex)及[中文记录](../notes/stage3-growing-box-windows.md)。增长界的C,c,d,C0存在但没有数值认证；显式渐近时间族不等同于数值可用质量阈值。来源缓存总数未增加。
