# 编译与视觉核查

2026-10-31，6页，Alan Day / Chengheng署名及AI研究助手脚注。

实际运行code/build_principal_eikonal.ps1，初版日志检出第2页10.54pt公式溢出，保留build-initial-overflow.log；拆开A/B与W定义后双遍pdfLaTeX exit0。桌面compile_latex_document对当前源返回success。最终main.log无未定义引用、Overfull或Fatal error。

用pdftoppm -scale-to 1100 -png渲染全部6页并逐页显示检查：页1模型、署名；页2相反动量旋转、负轴侧向核；页3正谱表示与归一；页4逐阶符号和余项；页5连续尾与数值表；页6晚时相消、物理边界及文献。未见缺字、裁切、公式重叠或断裂表格。Late u8符号未分辨明确保留。

TeX、最终PDF、编译脚本、最终/初始日志均由evidence.json绑定；桌面编译没有替代本地PDF交付。后完成727项原始算术重放另记录于audit.json，不冒称已写入PDF。
