# 中央产品传递与交叉项：来源接口

2026-11-06 本轮实际在线回读 Tall–Tomsovic，*Double-scaled bosonic and fermionic embedded ensembles, complex SYK, and the dual Hilbert space*，arXiv:2604.14522v2，§2.1，尤其式(2.11)–(2.19)：https://arxiv.org/html/2604.14522v2 。复用其 fixed-charge irreducible covariance 框架。该处 primitive tensor 的 particle-hole symmetry 不能解释成一般同charge乘积的 half parity。

直接数学输入来自本项目已保存推导：

- `manuscript/susskind-memory-coherent-cross/main.tex` 第2节：一般charge T_a(t)、coherent covariance normalization；第3节的 signed completeness 限于 half，本轮不借此证明 general exchanged Casimir。
- `manuscript/susskind-memory-parity-gram/main.tex`：半填充 (B,L_in,S) Gram，作为等距传递的已知参考。
- `manuscript/susskind-memory-central-window/main.tex`：正 output 六阶矩与中央紧窗一致界。
- `manuscript/susskind-memory-central-norm/main.tex`：一般 Rem 范数、实际同 J A 范数、输入加权尾及 (B,L_in,O) Gram。

本轮工作是将已有 Rem norm 归一成等距映射，以实际 A 作为相邻product差的控制量，经 O(p) 步 telescoping 与两个端点六阶尾得到 S 的 O(p^−2/3) 传递。不是把标准 U(N) covariance 或 Racah 恒等式重新命名为原创理论。

候选 exchanged Casimir 已有15纯rank pair的精确CAR核对，但一般解析 gradient contraction 还待完成。624条超几何 contiguous 检查只支持其标量化后的恒等式，不能代替那一算符步骤。后续如找到文献中的现成矩形 recoupling/Casimir 插入定理，先匹配规范再引用。
