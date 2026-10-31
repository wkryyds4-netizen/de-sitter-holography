# recoil-global 编译与逐页检查

2026-10-30；5页。Alan Day 与 Chengheng，首页脚注 AI research assistant (OpenAI)。

本机 pdfLaTeX 两遍成功；最终日志无 overfull、未定义引用和重编译请求。初版长路径产生一处 overfull，改用可换行的 path 后重新双编译。桌面 LaTeX 编译器对最终源码返回 success。

最终PDF以 pdftoppm -r 95 转成5张图，逐页视觉检查：第1页标题、双方署名和脚注完整；第2页长相位式与Gamma比未越界；第3页留数和及两支渐近清晰；第4页两张数据表无重叠；第5页路径换行、下一步与参考文献均在页内。未发现裁切、公式重叠或缺字。

检查为作者自查。数值34项诊断通过与排版检查分别记录；没有把编译成功当作数学证明或完整全息对偶认证。
