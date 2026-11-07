# 二阶与三阶 singlet 的联合弦规则和混合对易子

2026-11-06，本地阶段。[LaTeX](../manuscript/susskind-mixed-cubic/main.tex) 与 [六页 PDF](../output/pdf/de-sitter-susskind-mixed-cubic.pdf) 给完整定义和证明。无代理、commit、push 或 ZIP。

## 规范与结果

标准 CAR，含重叠指标的 complex Gaussian p-body Hamiltonian，退火归一 trace、β=μchem=0、ℏ=1。先 N,p→∞、λ=p²/N>0 固定，q=exp(−4λ)、σ²=Js²/λ；小 λ 随后取。

沿用扣守恒二阶 Y 和纯第三 coupling-chaos 残量 R，X=Y/√W、Rhat=R/√VR、h=H/σ。所有固定阶联合矩等于只配对同 species 的已有 multiprobe chord rule；交叉 HH:q、HX:q²、HR:q³、XX:q⁴、XR:q⁶、RR:q⁹。所以 τ(X Rhat²)=0、τ(X² Rhat²)=1、τ(X Rhat X Rhat)=q⁶、τ(Rhat⁴)=2+q⁹。这里是极限矩规则，微观 X/R 仍由同一组 H 耦合构成。

## 任意阶组件证明

二阶/三阶vertex占连续2/3个H槽，每个vertex至多一个导数mark；先删vertex内部Gaussian配对。DS局部指数仅有单mark α=2−4χ和双mark δ；δ=4仅当选中的两弦为嵌套双桥。

至少两桥才在composite vertices间连边。2/3槽使度数至多1，组件只能是单点或两点。固定外部端点、排列组件内槽，spectator交叉总数不变。外mark若两弦都终止于同一组件vertex，它就是双桥邻居、已在组件内；若终止于不同vertex，块次序固定。因此外mark和外δ不随内部排列改变，可以提出局部求和。

14模板、148局部配对穷尽：单Y/单R、双桥YR/RY、仅双桥RR全消，满桥YY/RR给GX/GR。两块余外端点颜色次序及两弧分配全部包括。压缩满桥同种配对，交叉给q^(d_i d_j)，即全固定阶定理。交叉乘积按无序弦对，只算一次。

GX=(1+q)²W/σ⁴，GR=D²VR/σ⁶，D=(1+q²)(1+q+q²)。原始代码X指未归一Y，清去每个X的1+q分母；R清去D，底层每个H槽有σ。不能把代码raw多项式直接当归一矩。

固定实时延拓用Taylor余项的非交换Lp/Hölder界：先固定阶取DS，再用极限H的紧q-Hermite支撑令阶数趋∞；未假设有限N统一算符范数，也未提供增长时间或λ→0的一致界。

## 现成四点核与动力学

读BINT §4.2 (4.10)/脚注8和BNR §8.1；不等长a=q²、b=q³代入已有核，恢复BINT省去的整体q⁶。不同species交替四点仅一个crossed channel，不加同种场的两个uncrossed channel。

F23(t)=τ[X(t)Rhat X(t)Rhat]；τ[X(t)²Rhat²]=1，故D23(t)=τ([X(t),Rhat]†[X(t),Rhat])=2(1−Re F23)。Jacobi正和及解析截断界见稿第4–5页。

F23=q⁶−2q⁶(1−q²)(1−q³)(σt)²+c4(q)(σt)⁴+…。

D23=2(1−q⁶)+4σ²q⁶(1−q²)(1−q³)t²+O(t⁴)。随后小λ：等时48λ+O(λ²)，t²系数384λJs²+O(λ²)。这是两种不同singlet的微观相互作用数据；q-Wick结构不等于普通可交换Gaussian场无相互作用。

## 验证与后续

12词、15392配对、433964局部minor；所有与probe预测的差为零多项式，q/λ作为独立符号。独立旧determinant核对14模板/148配对，共用pairing generator；全阶覆盖由组件证明提供。

独立有限CAR64词/6656选择，N3p1与N4p2，全精确相等，1031非零、2负值，含相干重复指标。核对mark代数，不是有限尺寸外推。

OTOC用λ=.2/.1/.05、Js=1、t=0/.01/.25/.5/1，cutoff52/90/167且各增16。最大绝对cutoff差6.39e−16、species交换差3.34e−16；独立chord Taylor到t6比较差≤1.27e−15。解析Jacobi尾界|t|≤1最高9.87e−12；浮点是经验精度而非区间认证。相对比较及相对尾界存verification-summary。正和主计算无大相位相消，Taylor系数先符号相消再求值。

小λ的t0/t2/t4/t6系数48、384、512、4096/15，吻合48cosh(4Js t)前四项。下一步查已有unequal-probe首阶四点推完整固定时间函数，不由四个Taylor系数宣布成立；再接两species Rindler空间、平移和体顶点字典。F任意混合高点尚未被本轮定理包含。
