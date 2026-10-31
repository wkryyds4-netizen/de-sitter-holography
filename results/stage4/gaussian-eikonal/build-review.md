# 编译与逐页核版

2026-10-31，6页，Alan Day / Chengheng，AI research assistant脚注。

实际运行 code/build_gaussian_eikonal.ps1，pdfLaTeX双遍exit0；桌面compile_latex_document返回success。后核原文作者全名，将本轮新稿参考文献中Chen/Tang名字首字母改为Y./H.，再双遍exit0及桌面success。旧科学稿未修改。最终日志main.log无undefined引用、Overfull或Fatal error。MiKTeX有更新提示，编译本身成功。

实际用pdftoppm -scale-to 1100 -png渲染并逐页显示6页：页1范围与不同探针；页2动量积分到三变量；页3尾界与虚部符号；页4固定阶/boosted极限；页5数值及表格；页6来源与下一步。最终作者修正后第6页重新渲染核查；其余页内容不变。无缺字、裁切、公式碰撞或不可读引用。

TeX/PDF/脚本/最终日志由evidence.json绑定哈希。诊断52项、加阶及符号检查11项，检查目的与误差限制见正文；版式检查不验证量子引力解释。
