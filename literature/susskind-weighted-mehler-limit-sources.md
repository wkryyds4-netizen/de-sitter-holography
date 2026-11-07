# 带权 Mehler 极限的文献接口

2026-11-08。先核已有精确核和特殊函数恒等式，再计算本项目的选态权重。未作穷尽式优先权检索。

| 来源 | 使用位置与范围 |
| --- | --- |
| Berkooz–Isachenkov–Narovlansky–Torrents, [1811.02584v2 §3.2](https://arxiv.org/html/1811.02584v2#S3.SS2)，(3.11)–(3.16) | q-Hermite Poisson/Mehler 核及两能量测度。它是已有公式；本项目先前证明的同 J 曲率归约才允许在这里使用。 |
| DLMF [20.5.1](https://dlmf.nist.gov/20.5.E1) | theta1 乘积；注意本稿 theta nome 为 exp(−epsilon/2)，而 Pochhammer q 为 exp(−epsilon)。 |
| DLMF [20.7.30](https://dlmf.nist.gov/20.7.E30) | tau→−1/tau 模变换，得到反射 Gaussian 格点和。 |
| DLMF [23.18.5](https://dlmf.nist.gov/23.18.E5) | Dedekind eta 变换，给 (q;q)_infinity 的精确对偶乘积。 |
| Mukhametzhanov, [Large p SYK from chord diagrams, 2303.03474v3](https://arxiv.org/html/2303.03474v3)，§§3–4 | 已有态密度和两点函数渐近的前驱。当前不声称该文已提供本稿所需全部带权 C^k 余项，也不声称本文有优先权。 |
| Lin–Stanford, [2307.15725](https://arxiv.org/html/2307.15725)，(137)–(138)、(409)–(410) | 上阶段热鞍点预测的来源；本阶段直接带权积分证明不依赖未经证明的温度求导交换。 |

本项目依赖：selected-state 的协变密度、sector-bilocal 的实际制备接口、selected-rindler-response 的归一源恒等式。旧源响应稿的条件性结论和封存文件保留。当前新增工作集中在同一归一化下的谱边缘支配、带权首阶余项与曲率两态混合；不把既有 Mehler 核、Fourier 恒等式或 leading sech 核写成原创。
