# 固定时间首阶四点：原文与本项目代入

2026-11-06。文献先行；本轮没有重新发明四点R-matrix或宣称新的混沌指数。

1. [Lin–Stanford, A symmetry algebra in double-scaled SYK](https://arxiv.org/html/2307.15725)，当前HTML标明v1。实际读§4.6 (186)–(205)、附录G (375)–(392)，不是只读摘要。此前§4.2 bra/ket fake-circle与§4.3–4.5提供多点插入背景。本文核心引用(205) crossed与G.2 (388) energy-exchange。source λ映射为h=4λ，Hsrc=H/(2J_s)，源时间u=2J_s t，ΔV=ΔW=d=2。

crossed在β→0、c→1后只剩 −h d² cosτ−/(cosτ+ cosτ+')。采用(13)(24)配对的两个中点和差坐标，τ+=i(u3−u1)、τ+'=i(u4−u2)、τ−=i(u1+u3−u2−u4)。等时得到−hd²，与q^{d²}−1一致；互探针权已计入。

G.1 (382)–(387)的A=h(1−c²)/(c²+c√(1−c²)arccos c)、b=c/√(1−c²)+arccos c满足Ab²→h，必须先合并再取极限。G.2 (388)与(392)更直接给leading energy-exchange和高能量多项式的阶数。用全迹两能量测度的交换/符号对称性导出τ(Hsrc XiXj)=−i∂sG/2；第二uncrossed按循环排序使用T41。三个收缩是本同species应用需要的组合。

2. [Streicher, SYK Correlators for All Energies](https://arxiv.org/html/1911.10171)，HTML v1。实际读§3.1 (3.2),(3.3)和§3.2 (3.4),(3.5)，确认有限温度OPE/OTO各自的角色及延拓。该文q指费米相互作用阶数，不等于我们q=e^{-h}。计算系数采用Lin–Stanford已转换到chord规范的公式，避免混淆两套参数。

3. [BINT, Towards a full solution of the large N double-scaled SYK model](https://arxiv.org/html/1811.02584v2)，沿用上一阶段实际核过的(3.9)、(4.9)、(4.10)。本轮把四段Uj=exp(iajJ)保持独立，得到两个uncrossed和一个复双线性crossed；以原非对称T和指标和作独立50位核对。来源定位见susskind-singlet-otoc-sources.md。这里没有宣称重新阅读全文。

本轮结果属于现成半经典公式对已构造singlet的应用与计算验证。它解决固定时间O(λ)微观四点系数，未给独立空间KG插入/质量谱，也没有给固定时间余项的统一界或一般四时间总截断认证。下一步从源文§2–3、§4.1、附录D的生成元和多粒子表示继续，核是否能构造独立平移观测量。
