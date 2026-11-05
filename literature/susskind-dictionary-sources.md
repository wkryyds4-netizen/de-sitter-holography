# complex DSSYK：规范与可复用解

2026-11-05。本轮为Susskind主线第二步，目标是把选定的singlet谱问题接到已有complex SYK方法上。仅在线原文阅读；未新增第三方全文缓存。

| 原文及版本 | 实读位置 | 采用方式 |
|---|---|---|
| [Susskind 2303.00792v1](https://arxiv.org/html/2303.00792v1) | §3.2 (3.9)–(3.15)、§5.2 (5.14)–(5.31)、§6.1 | complex CAR=2与cosmic方差；singlet电荷守恒机制。图形仅读图注，未重算所有图系数 |
| [Miyashita–Sekino–Susskind 2506.18054v3](https://arxiv.org/html/2506.18054v3) | §2 (2.2)–(2.6)、§3 | string方差与归一迹；采用Hermitian Gaussian协方差明确有限N实现 |
| [Sekino–Susskind 2501.09423v2](https://arxiv.org/html/2501.09423v2) | §§3.3–4.2、附录A.1及A.2至(A.24) | Majorana λ=2p²/N、CAR=2、Hc=pHs及闭时间路径；单线self-energy方法参考，尚未完整读A.3 |
| [Rahman–Susskind 2401.08555v3](https://arxiv.org/html/2401.08555v3) | §1时间换算与§6 | cosmic/string时间方向，以及大p与等时导数不可直接交换的例子 |
| [Berkooz–Narovlansky–Raj 2006.13983v2](https://arxiv.org/html/2006.13983v2) | §2.1、§3.2–3.3相关正文、§4.1–4.2相关公式、§7定义及7.1、§8.1至(8.6) | 已有complex chord矩、CAR及组合方差；q=e^(−4λ)。§7–8是独立Gaussian长探针，后续需构造同J的导数singlet插入；未全文复现73页 |
| [Tall–Tomsovic 2604.14522v2](https://arxiv.org/html/2604.14522v2) | 引言、§2.1至(2.11)、§4.1、§6 | 已有固定粒子数二阶矩及strength-density方法；§4独立探针结果的范围需保留。§3证明与§4.2四点细节为下一阅读起点 |

2006.13983最新版本为2020-10-15 v2，2604.14522最新版本为2026-05-11 v2，均由arXiv摘要页核实。一次批量请求连接失败后重新读原站成功，不据失败请求宣称读过内容。

## 直接复用与本轮整理

先匹配BNR的binom(N,p)⁻²方差到MSS的(p!/Nᵖ)²方差；两者比值在双标度极限有有限指数因子，不能用固定p近似抹掉。再将BNR已有chord矩转为本项目单位。固定电荷的二阶矩直接用Tall–Tomsovic (2.8)–(2.9)，以标准CAR复核，不重复宣称新谱密度。

本轮的A₂=Σc†H²c−(Q−2p)H²及相邻电荷扇区公式来自正常序CAR和上一稿恒等式。它给出本项目同J导数singlet的具体插入位置。尚未查明该专门整理是否有逐式相同先例，因此只称标准代数应用。短关键词检索不构成穷尽查新。

## 下一步怎么用现有工作

以固定电荷的四点strength density为工作对象，保留c、c†带来的m↔m−1转移和两个中间能量平方权重。先读2604.14522 §3、§4.2与附录A，再对照BNR §3/附录C的定向弦规则。独立随机中性探针的终式不能直接替代同J、固定指标收缩的A₂，但其正交多项式和收缩方法可复用。Gross–Rosenhaus仍为后续相互作用方法参考；当前优先级让给直接适用于complex模型的原文。

复现代码为`code/complex_syk_dictionary.py`、`code/run_complex_syk_dictionary.py`、`code/audit_complex_syk_dictionary.py`与`code/check_singlet_charge_bridge.py`，分别保存模型、样本生成、独立扇区归一化复核和相干四点重建检查。
