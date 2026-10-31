# 编译与视觉核查

2026-10-31。源 manuscript/ds3-compact-global/main.tex；交付 output/pdf/de-sitter-ds3-compact-global.pdf，共7页。Alan Day / Chengheng，AI research assistant脚注保持。

本机命令：powershell -NoProfile -ExecutionPolicy Bypass -File code/build_ds3_compact_global.ps1。脚本运行两遍 pdfLaTeX，最终 exit0，main.log无未定义引用、Overfull或重编译请求。首次直接 -File 被本机执行策略阻止；process级命令执行后第一版在第7页发现两处长路径溢出，改用可换行path命令，重新双编译通过。原溢出日志留在 tmp/pdfs/ds3-compact-global-build/initial-overfull.log。

桌面 compile_latex_document 对最终保存源码返回 success。本机 Poppler 用 pdftoppm -scale-to 1100 -png 渲染；7页均实际逐页查看：署名脚注、公式、表格、引用、页码完整可读，无裁切、重叠或缺字。第5页长三角不等式保持在版心内。输入/输出SHA256登记在本阶段evidence.json。

此记录仅说明实际编译和版式检查，不替代科研推导、数值误差认证或同行审查。
