# 同一J曲率的开放双弦递推

2026-11-08。正式稿：manuscript/susskind-slater-composite-transfer/main.tex；纯函数：code/slater_composite_transfer.py。沿用同原Gaussian J、半填充固定Slater annealed态、固定lambda=p²/N>0的DS极限。h=H_m/√C_m，W=√lambda u，j=C_m^-1(2p-D-kappa)(H²-C_Q)|m，q=exp(-4lambda)，gamma=exp(-lambda)，kappa=2lambda(1-q)/(1+q)。D/L先fullFock作用、后限制；无量纲能量、hbar=1。j的方差v=4lambda+16lambda²q/(1+q)，实际归一曲率为前一阶段的a_N=(1-P_H)x_N，x=j/√v。

## 1. 完整源几何，而非有限词猜测

一份归一化pair-collision只含两个P/Q标签。未含这两个标签的端点是I；只有一个标签经过的marked块，其odd matrix被1-z/2缩放，恰被singleton分母抵消。因此任何源系数只由两条涉及弦决定，旁观弦不能改变它；JW符号也只涉及这两弦。

同一edge的P,Q碰撞：whole-D作用在该edge一个端点时D(n)=n，与两个singleton因子完全相消；L(n)=n-I/2可留下源。G edge的每个W端给K_W=1，两端都是W时K_WW=1；F edge的P,Q同site被selected约束禁止，所以这些系数都为0。

不同edge的两个标签：L块只有一个incidence，消掉；D源只有该J的两个相邻槽都属于这两弦时才留下。设J槽a,a+1的伙伴分别r,s，epsilon=+1（不交叉）或-1（交叉）。若r,s都在J右方，取右端点较年轻的那条e_*（右端点数值较小）；若都在左方，取扫描时先闭合的那条e_*（左端点数值较大）。K_J=2epsilon 1_{e_*为G}。伙伴分居两侧时K_J=2epsilon，无G指示。两个J源同时留下只可能两弦完整连接两J块；nested给K_JJ=4，crossed为0。WJ、所有三阶及以上源恒零。该论证覆盖任意固定数量的旁观弦。

减去原中心化后，每个W只可能与同G edge上的W成对（权lambda），或与selected端点相连（额外权-lambda）。每个J在源指数的partition展开里，要么单独分配权lambda K_J-kappa，要么与另一个J成对并指定nested双G edge（权4lambda）。这是传递构造的全部额外信息。

## 2. 计算状态与单端点

从右向左扫描。开放列表按右端点从最老到最年轻排列；五种端点颜色：H=右H、预定G；h=右H、预定F；U=右W、预定G；w=右W、预定F；P=右M†、必为F。预定类型只约束未来闭合，最终每张完整图恰计一次。开放H/W分别产生H+h或U+w两个分支，M†只加P，M只闭合。

闭合第i个（0起）开放端点，列表长n：G权q^(n-1-i)；F权q^(n-1)，后者包括older开放弦对inner-F nesting的额外q^i。除这个q权外，允许的闭合为H<-H:1、H<-P:gamma、W<-U:lambda、W<-P:-lambda gamma、M<-h:gamma、M<-w:-lambda gamma、M<-P:1。其他为0。

另有互不共享端点的link，S或R。S是一个已经扫描的singleton-assigned J，其两条弦仍开放；任一端闭合时乘2lambda epsilon 1_G-kappa，较年轻端先闭合epsilon=+1，较老先闭合epsilon=-1，随后删除link，另一端变普通。R是等待partner J的nested双G弦，单端点操作禁止闭合它。

## 3. J的显式作用

先处理singleton分支，将两个H端点从右到左处理；不得在本块内部配对。OO：两端开放，建立S，权留待第一次闭合。OC：右端开、左端关旧edge，乘-2lambda-kappa。CO：右端关、左端开，乘2lambda-kappa。CC：依次关闭原列表i、k，i!=k；epsilon=+1若i>k，否则-1，e_*为max(i,k)的原端点，乘2lambda epsilon 1_G(e_*)-kappa。每次close仍处理它所遇的旧S link。

再加source-pair分支：开两个H端点并建立R；或取一个已有R=(a,b)，a<b，先关b再关a，乘4lambda。两次G闭合的q权为q^[(n-1-b)+(n-2-a)]。R内部不会被其他单端点或singleton J分支拆开。与别的弦的crossing/nesting仍由普通闭弦q权准确计数。

删除端点后重编号，合并完全相同的开放列表和link集合；闭合过的历史全部丢弃。代码中的future pruning只用剩余槽数、剩余M/W/J个数等必要条件，删除不可能完成的状态；它不是数值截断。最坏状态数仍会增长，未证明统一多项式复杂度。

## 4. 与源公式、物理内积的关系

每张合法forward图加一个源partition，唯一决定扫描路径；反之每条回到empty的路径唯一恢复该decorated图。闭弦权还原gamma^b q^nu，S和R还原所有singleton/pair source系数，故任意固定word w严格满足phi(w)=epsilon T_w empty。只在有限图和有限源partition之间重排，不交换N极限或无穷和。

开放列表是auxiliary状态，不赋Euclidean metric。取可达空间V=span{T_P empty}，定义G(T_P empty,T_Q empty)=phi(P†Q)，右边用同一递推计算。若T_P empty=0，则任意左context R有phi(RP)=0，故此定义不依赖P的表示。旧有限N正性使G>=0；商去null后就是共同GNS。T_A的adjoint是T_(A†)，h,W,j symmetric，T_+†=T_-。这只给reachable quotient的物理metric，没有证明任意形式状态都有正metric或已得到最小正交Fock基。

旧covariant bounded-normal theorem仍适用：||h||<=2/√(1-q)，||W||<=2√lambda/√(1-q)，||j||<=2√v/√(1-q^4)。所以极限空间中的h演化可用Taylor多项式及e^z z^(K+1)/(K+1)!余项控制，z=|t| ||h||；fixed preps两侧再乘向量范数。M有界/compact在本扩大空间尚未得出，也未自动补全finiteN interlaced nonlinear pinching的时间尾。

## 5. 已执行核对和继续位置

首次完成科学运行：2774个源指数与局部CAR引擎全mask相等；497词在两组独立有理(gamma,lambda)点给994矩比较；4096个旧六字词回到H/U/P递推。两个14维Gram精确正LDL；5个unpruned/剪枝比较。长词H^24、j^8、jhhjWWjhhj、M j h^4 j M†、M² j h² j (M†)²分别峰值7/21/10/18/22状态，与经典Touchard--Riordan、旧central mixed-q或原selected源引擎相等。状态峰值依赖这些具体词和参数，不是一般复杂度上界。有理参数不在物理gamma=e^-lambda曲线上，只作代数诊断。

初次driver在长H^24的参照计算误选全pairings，主动中断；此前检查无assertion失败、没有科学JSON。改为已知Touchard--Riordan有限和后完整通过，成功输出封存。未重跑旧driver。另一次纯helper探索在gamma=1/2,lambda=1/3处发现phi(M j^n M†)=phi(j^n)对n=1..8成立，未收入本轮科学证据或提升为全阶物理定理；下一步值得先检验单selected制备在covariant三场代数上的完整密度，再用新递推处理高后代/完整bilocal。正式研究仍按原文优先，空间/质量生成元和finiteG保持主目标。
