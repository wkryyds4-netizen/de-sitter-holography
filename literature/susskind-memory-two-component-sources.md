# 两分量归约：本轮原文接口

2026-11-07，先回读原文再作 strong-regression 计算。与上一轮共用的数学背景不再包装为新定理。

- Bożejko–Kümmerer–Speicher，[funct-an/9604010 PDF](https://arxiv.org/pdf/funct-an/9604010)：Def1.1/Remark1.2（印刷p5）给q-Fock场作用和有界性；Prop2.7 proof/Prop2.9（p10）给Wick左右递推与单方向Hermite；Thm2.11（p13）给正交投影诱导条件期望；Thm2.13（p14）给C*层连续性。root直接打开并读Thm2.11–2.13；两个数学协作者亦核这些接口。三色两block因子化是投影关系的应用，四u背景是左右Wick递推的显式推论。并未将颜色称为经典独立变量。
- Berkooz–Narovlansky–Raj，[2006.13983v2 §3.3](https://arxiv.org/html/2006.13983v2#S3.SS3)，(3.11)–(3.15)以及§4.1：有限CAR trace到固定长度Poisson-overlap/chord的规则。root本轮重新阅读。独立Gaussian颜色只禁止异色配对，允许图的局部trace保持原样；因此可用来扩展三色joint moments。实际charged force没有被这些独立辅助变量替代。
- Tall–Tomsovic，[2604.14522v2 §2.1](https://arxiv.org/html/2604.14522v2#S2.SS1)，(2.12)、(2.17)–(2.21)：中性covariance的rank特征值、有限次数交叉权重。原文明确次数k固定，不能以此直接断言整个实时间指数的uniform误差；本轮另用Gaussian moment给总degree可和尾。也不把其spectral Wick与有限N coupling chaos混同。

本地依赖：uniform-bound（实际场与A/B的有界性、weighted force、Dyson余量），return-background（实际场到X0），central-spectrum（Y±谱与S桥），central-raw（F/BZ和J_hF桥、input/output尾），mixed-spectrum（R/F/Y±交叉），charge-sum（旧Y/L全迹）。本轮新增的是four-u的完整strong regression/pinching消失、实际two-component归约和neutral低系数识别。检索与复核不构成原创优先权声明。
