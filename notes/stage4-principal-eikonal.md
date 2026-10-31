# 实 principal 探针的 eikonal 核与 Gaussian 正谱矩

2026-10-31。详见[LaTeX](../manuscript/principal-eikonal/main.tex)、[PDF](../output/pdf/de-sitter-principal-eikonal.pdf)及[来源核查](../literature/stage4-principal-eikonal-sources.md)。本轮推进原目标所需的真实principal有限G近似接口，不宣称完整对偶、微观Hilbert构造或原创优先权。

## 明确的模型与推导

ℓ=1，m>0，γ=4Gm∈[0,1)，α=2π(1−γ)。两种不同species，ν,ν′>0为实principal，物质质量平方1+ν²、1+ν′²。沿用CSTY固定平均观察者质量混合表示和12/14轮廓，两个sector相加后减去重复断开项。尚未证明它等同上一轮G=0量子观察者Fock真空。

四个独立d_i∼N(0,σ²)，保持原半周期插入位置。X=(d₁−d₃)/2、Y=(d₂−d₄)/2、Z=(d₁+d₃−d₂−d₄)/2相互独立，方差σ²/2、σ²/2、σ²。h=Hγ/m，Hγ=−cos[(1−γ)(π−φ)]/[2sinc(πγ)]，角平均π⁻¹∫₀π。

从原(2.32)做相反±π/4动量旋转，保持p+p−实。两对Bessel参数均有相位z=e^(−iπ/4)，整体前因子为−1，震波项exp(iλxy)。使用标准K的rapidity表示后先积y，再设v=zAx：v积分的射线在正轴下方，负λ时从正极点下方旋到实轴。由此确定η=+，不是数值后取实部。

令A=e^Ycoshr₁+e^(−Y)coshr₂，B=e^Xcoshr₃+e^(−X)coshr₄，W=cosνr₁cosνr₂cosν′r₃cosν′r₄，则单sector

\[
S_\eta(\lambda;X,Y)=(\pi\alpha)^{-2}\int_{\mathbb R_+^4}\frac{W}{A^2B^2}\mathcal K_\eta(\lambda/(AB))d^4r,
\quad\mathcal K=-J'=((1+c)J-1)/c^2.
\]

J(c)=e^(1/c)E₁(1/c)/c，Jη(−1/a)=ae^(−a)Ei(a)−ηiπae^(−a)。ImKη=−ηπa²(a−1)e^(−a)，自身可变号；W也变号。上一轮共形负号定理不可直接推广。

I_nν(X)=∫₀∞x^(n+1)K_iν(xe^X)K_iν(xe^(−X))dx是f_n(r)=e^((n+2)r/2)K_iν(e^r)的自相关。标准Mellin变换给M_nν(ω)=2^(n/2−1−iω)Γ((n+2)/4+i(ν−ω)/2)Γ((n+2)/4−i(ν+ω)/2)，从而

\[
B_{n,\nu}=\E_XI_{n,\nu}(X)=\frac1{2\pi}\int_\mathbb R e^{-\sigma^2\omega^2}|M_{n,\nu}(\omega)|^2d\omega>0.
\]

I₀ν(X)=πsin(2νX)/[2sinhπν sinh2X]可变号，但平均B₀ν严格正。正确断开基准Dσ=B₀νB₀ν′/(πα)²。完整指定模型R=[E〈Sη(he^(t+Z))+S−η(he^(−t−Z))〉−Dσ]/Dσ。不能逐点除以变号的I₀；也不能把Mellin ω无推导地称作模能量或碰撞能量。

固定γ、σ、t、标签时，各固定阶系数为(−1)^n〈h^n〉exp(n²σ²/2)B_nνB_nν′/(n!B₀νB₀ν′)，乘(e^nt+e^(−nt))。γ>0时shock奇矩负、偶矩正，所以这些系数全正；γ=0只有偶阶。这是模型固定阶渐近，不是整个级数收敛或非微扰答案单调。稿中给逐阶余项的存在界，保留n增长和Gaussian矩的依赖。

|K|有有限全实边值上界C_K；由∫drds/A²=X/sinh2X≤1/2，得到完整Gaussian可积性、rapidity和Gaussian连续尾界。未给C_K的认证数值，更不把连续尾冒充离散求积或省略引力项的误差。固定t的大m极限为1，t=logm+u的boosted极限保持非平凡，次序不可混。

## 实际计算和限制

主数据[study](../results/stage4/principal-eikonal/study.json)取ν=.6、ν′=.9、m1000，四组(σ,γ,u)=(.4,0,0),(.4,.25,0),(.8,.25,0),(.8,.25,2)。每组8个scramble、2^18与2^20个点，8维积分完整保留四rapidity、三Gaussian和角。断开项用独立50位谱积分；同时保存raw与断开项control-variate估计。没有用随机分母强制归一。

60项[独立检查](../results/stage4/principal-eikonal/diagnostics.json)通过：小/负参数lateral核、共轭、16个涂抹矩的独立rapidity求积、四个两点闭式Gaussian平均、两个正耦合完整Bessel积分、共形归一和原θ前因子/Jacobian检查。涂抹矩最大相对差1.88e−11，双精度容差非区间认证。

[进一步核查](../results/stage4/principal-eikonal/refinement.json)保留三组相反旋转的原复Bessel积分，正负λ皆核对；另用全新8个seed、重要性率.7、2^21点再算四组，14项明确经验判据通过。主四组R约1.10266−.04401i、1.24688−.08601i、1.19389−.19901i、.30011−.81574i。最后一组分量经验标准误约.000239/.000280，加点并非单调改善，不能当严格置信区间。

[晚时刻扫描](../results/stage4/principal-eikonal/late.json)另存48组seed运行，σ=.4、γ=.25、u4/6/8，各以rate1的2^20和rate.7的2^21计算。u4的ImR约+.02768，SE约.000056，相较u0的负值提供涂抹后变号的数值证据。u6 raw R≈.001038−.000388i，SE≈.000037/.000031，绝对被积函数平均约.1941，已有强相消。u8虚部约+.000034但SE约.000024，符号未分辨，不能宣称稳定振荡。控制变量在late实部未必优于raw，两者原值均保存。

总计144次完整scrambled积分（64主+32加阶+48晚时刻），这些是经验数值，不是完整引力误差认证。负轴固定插入的确定性48/72阶核对另保留，不用它替代角平均与Gaussian平均。所有原代码、参数、版本及输入哈希随结果保存。

## 下一步

本轮得到可与真实principal自由基准比较的有限G候选近似及正谱矩；尚未解决态、有限κ插入归一、Gram逆核、真实缝合、完整四点引力修正或微观对偶。后续优先核实际物理制备和微观插入数据，而不无限增加同一近似中的数值精度。晚时刻相消若用于任何物理主张，必须另建更稳定表示及相对误差控制。
