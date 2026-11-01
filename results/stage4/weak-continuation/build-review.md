# 五页阶段稿编译与检查

2026-11-01。源 manuscript/weak-continuation/main.tex，输出 output/pdf/de-sitter-weak-continuation.pdf；Alan Day / Chengheng 署名及 AI 助手脚注。输入、输出与日志哈希由 evidence.json 绑定。

实际运行 code/build_weak_continuation.ps1：两遍 pdfLaTeX。初版两处溢出日志保存为 build-initial.log；调整长段和展示公式后重新双遍 exit0，最终 main.log 无 Overfull、未定义引用或待重编译提示。桌面 compile_latex_document 返回 success。

使用 pdftoppm -scale-to 1100 -png 渲染五页，已逐页检查署名、公式、数值表和文献，无裁切或重叠。此记录不要求再次运行已成功的编译或 PDF marker。

第2页明确单极点邻域，第3页区别字面积分与亚纯延拓，第4页数值仅作经验验证，第5页保留第一象限内部及未测试有界五边形的边界。Ramond PDF 第88页另行全页核读，截图不公开。
