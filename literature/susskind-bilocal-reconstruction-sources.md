# Bilocal strong reconstruction：本轮来源与新增接口

2026-11-07。沿用 [上一轮memory-to-bilocal原文笔记](susskind-memory-to-bilocal-sources.md) 定位出的有限Volterra起点和Itakura/CDP接口。本轮先读实际旧证明，再完成此前待证的放大余项；不重复已有bulk谱。

1. [bilocal-channel](../notes/susskind-bilocal-channel.md) 已有 V、K、R、M 的有限维压缩、Schur补/Volterra方程和 S_c/m趋1。投影消去bath和resolvent差式是已有方法，本轮新增的接口是对 M0=P_hM 的逐完整能量块bath直和dilation及其在集体尺度重建中的应用。

2. [return-expansion](../manuscript/susskind-memory-return-expansion/main.tex) 的force moment公式：Gaussian covariance v_N、Etr A_m^k <=(2k-1)!!(Nv_N)^k、G<=A/m。本轮原样继承k=2/4/8，不新做四源任意高矩计算。对annealed Schatten迹使用标准Holder与2/8插值，再结合 [uniform-bound](../notes/susskind-memory-uniform-bound.md) 的实际memory二阶界，才得到O(p^-1/3)强返回率。不能从一般矩阵值Gaussian多项式的二范数裸推四范数。

3. [joint-regression](../notes/susskind-memory-joint-regression.md) 已有完整共同谱nu_ab、ordinary/Hermite可和尾、exact pinching、charged force与neutral幅度。本轮积分保持原mixed矩阵；charged Green消去与neutral差商是作用于这些已证输入的具体推导，不把旧A0B2校准写成新发现。

4. 本轮在线核对 [NIST DLMF §3.3(iii)](https://dlmf.nist.gov/3.3.iii)，(3.3.34)–(3.3.37) 的标准差商、显式插值式和解析contour表示。重节点通过连续延拓；对复指数的t²/2界在本稿用simplex积分直接推导，没有把仅适用于实函数的单一中值点公式套到复函数。差商及三角矩阵函数表示属于既有数值分析工具。

5. 本轮也回看 [BNR 2006.13983v2 §3.3](https://arxiv.org/html/2006.13983v2) 的有限集合交叠/Poisson组织，确认它是后续marked多点扩展的已有入口。最终强返回证明用旧正力高矩已足够，本轮没有建立或使用新的四源Poisson极限定理。

6. 下一物理接口仍是 [Itakura hep-th/9604032v2](https://arxiv.org/abs/hep-th/9604032) §4–5 和 [Cavicchi–Di Vecchia–Pesando hep-th/9306091](https://arxiv.org/abs/hep-th/9306091) 的既有bilocal/collective方法；其具体公式位置与正规序、态、光前规范差异已在上一轮来源笔记核对。Gram metric归一化不等于Holstein–Primakoff canonical dressing。需要先确定本模型多点代数与适用态，再构造纵向变量及实际P^-；不选一个质量算符去拟合已知FLZ谱。

本轮成果是指定complex DSSYK模型内的强重建与谱幅度，不作原创性优先权宣称；既有数学工具和候选体理论的质量结果均保留原始归属。
