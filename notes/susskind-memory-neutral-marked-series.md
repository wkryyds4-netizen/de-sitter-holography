# 中性记忆场的全时间 marked series

2026-11-07，本地阶段。正式推导见 manuscript/susskind-memory-neutral-marked-series/main.tex，阅读版见 output/pdf/de-sitter-susskind-memory-neutral-marked-series.pdf。沿用共同原始 J 的重叠 complex p-body CAR 模型，β=μ=0 normalized annealed full trace，偶 N、p²/N→λ>0，固定 Js、ℏ=1/string 实时间，固定紧时间区间。这里的结论是模型内强 L² 表示，不是已经求出的联合极限核。

## 本轮闭合的接口

原实际场已有 Fscr=A+B+o(L²)，A 为相邻电荷传递缺陷，B 是 neutral Φ(exp(iαhs)) 的 mε 放大、对称相位包裹和精确能量 pinching。本轮把整个 B 归约为最高 Gaussian chaos 的一个 rank-marked 家族，保留所有低 chaos 收缩。

记 D=pNcal−Rcal，Zj=Wj(h^j)，Uj=DZj，θj=μ(jp)，aj=Σ(v<j)θv。有限参数 Jacobi 多项式 Qj 与系数 cj 满足 exp(ish)=ΣcjQj(h)，Ej=Zj−Qj(h)，ηj=pEj。θ(Ncal)绝不是逐样本能量函数。

精确重排为

    (p/C)(1−Ph)Φ(exp(ish))
    =(1−Ph){p[μ(Rcal)−θ(Ncal)]exp(ish)
             +Σ cj(s)[θj−θ(Ncal)]ηj}.

新增的强 multiplier 线性化保留有限 θ：

    p[μ(Rcal)−θ(Ncal)]exp(ish)
       −4λN θ(Ncal) D exp(ish) → 0.

证明先在 |m−N/2|≤Xp。旧六阶 Racah 证书是在整个 Q[A,B,z] 多项式环中相消，故可把输入 rank 扩到 r≤Kp（K 固定），得到条件 δ⁶ 的统一界。一次 commutator 给实际指数的 D-graph Taylor 尾 TK(S)→0；二次 commutator 由条件 δ⁴ 控制，从而固定 Taylor 多项式有 D² 界。关键新尾式是

    ||D 1(D>R)exp(ish)||₂ ≤ TK(S)+CK(S)/R.

先选 K，再选 R。这控制完整指数的 D² 谱权重统一可积性，不是仅最高 chaos 的二矩，也没有声称整个指数具有统一 D² 范数。固定 k,d 斜率与全 rank 的 16λ+ D 支配随后给乘子极限。

## 低 chaos 的受迫递推

diamond 是带有序矩阵系数的普通 coupling-Gaussian Wick 乘积，不是 q-Hermite 乘积。

    rj=Σ(v<j)[(Φ−Cθv)Zv] diamond Z(j−1−v),
    η(j+1)=hηj−Cajη(j−1)−p rj.

定义 finite proxy source

    rhatj=4λN C Σ(v<j)θv [Uv diamond Z(j−1−v)],
    etahat(j+1)=h etahatj−Caj etahat(j−1)−rhatj,
    etahat0=etahat1=etahat2=0.

固定 j 有 ηj−etahatj→0。Wick 乘积右 creation 界、全次数 rank 二矩和 θv 的几何衰减给 etahatj≤Lhat K*^j C^(j/2)√j!；cj 的 1/j! 系数使整个时间级数绝对可和。固定 j 的 ||DEj||²≤||Ej||||D²Ej||→0，再接 all-j DQj 界，得到 Dexp(ish)−ΣcjUj→0。

最终定义

    Vm(s)=Σj cj(s){4λN θj Uj+[θj−θ(Ncal)]etahatj},
    Btilde_m(s)=(mεm Cm/p)(1−Ph)e^(−ihs/2)Vm(αms)e^(−ihs/2).

广中央之外 Btilde=0。mε/p≤2、C≤M、广中央统一界及 binomial 尾把 compact charge 结论接到全迹，sup|s|≤S Σbm||B−Btilde||²→0。显式可求和尾在稿中 Eq.50；保持 αms 与 mε/p 精确，不提前换成极限常数。旧 Fscr 与 A/B 界意味着两插入协方差可替换成 A+Btilde，对外部实时间 T 一致，四个 cross 项都保留。

## 低阶校准与真正剩下的问题

令 b1=(Φ−Cθ1)h、b2=(Φ−Cθ2)Z2、d1=E tr(hb1)，则

    η3=−p b1,
    η4=−p{h,b1}−p b2+p d1 I.

反对易项已有 2d1 的 scalar contraction，故 η4 净 chaos0 为 −pd1I。η3 的未 pinching 极限仍有 −4λ²Cqh；在下一次乘法或 chaos 操作前不可删掉。最终 pinching 后恢复 −F/2 及 −2λCq²(Y−+Y+)−(1+q)JhF。

现在真正待求的是 h、Uj、有序 Wick 后代和 charged reverse-source defect 在任意固定次数下的共同极限矩及两能量谱。还须控制零频质量：矩收敛本身不允许交换精确能隙零投影与 N 极限。当前定理始终保留 finite-N samplewise pinching，不受这个后续问题影响。

## 执行证据

- coefficient-checks.json：三组 (3,1,2)/(4,1,2)/(4,2,3)，最大次数6/5/4，全部9/16/36个坐标。18 polynomial、12 source，共87比较，最大误差5.6793055687843776e−15。η4 scalar 范数 .25、.11111111111111176、.5185185185185176，明确非零。complex128 对照不是精确算术证明。
- rank-checks.json：4个完整系统，12组 k=3/4/6 六矩、38条件分布、350 Casimir 等式、6 pure-rank mixture；大参数54链及72条件行仅是精确有理上界。450 slope 支配精确通过；80/120位最大变化约1.943e−78。有限 slope 相对极限最大偏差约.9288，完整保留，未把接近极限当测试通过条件。
- 两科学输出均首次运行通过且 exclusive-create。旧 main 未重跑，旧证据不覆盖。九页 PDF 双遍编译并逐页渲染检查；编译方式与输入/输出哈希由 local-evidence.json 绑定。

本轮没有 Git 提交、推送、ZIP 或 issue 回复。issue #1 初核沿用已有记录；9890旧端点重现失败不改，21143单元独立重积分仍未做。
