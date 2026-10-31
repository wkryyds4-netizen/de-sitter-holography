# PDF 编译及核版记录

2026-10-31，9页，Alan Day 与 Chengheng，AI research assistant 脚注。源码 manuscript/ds3-fock-vacuum/main.tex；交付 output/pdf/de-sitter-ds3-fock-vacuum.pdf。

实际执行 code/build_ds3_fock_vacuum.ps1：本机 pdfLaTeX 双遍 exit0；随后桌面 compile_latex_document 返回 success。最终日志 main.log，双遍终端输出在本机 tmp/pdfs/ds3-fock-vacuum-build/pass1.txt、pass2.txt。最终无 undefined 引用、Overfull 或 Fatal error。MiKTeX 更新提示不影响此次编译成功。

初稿第2页0.90965pt overfull 已改句；原 initial-overfull.log 留在本机临时目录。最终用 pdftoppm -scale-to 1100 -png 渲染，page-1.png 至 page-9.png 已逐页显示核查：标题、脚注、各页公式、表格、链接和页码可读，无缺字、裁切、碰撞。第8页已含 Gaussian 24项、t=1极限及 KMS 残差；后续 Gaussian 独立115项仅记附属 notes/audit，不声称在 PDF 数值段内。

页1：模型及两种调节；页2：CG synthesis 和单步界；页3：完整重数与绝对群平均；页4：弱下降和解析向量；页5：真空返回及完整四点；页6：重质量谱及 KMS；页7：完整四点极限及边界；页8：数值与失败；页9：来源及有限 G 下一步。全部实际检查完成，无待检查页。

输入 TeX、输出 PDF、编译脚本及最终日志的 SHA256 由本目录 evidence.json 统一绑定。数值审计与版式检查的职责不同，PDF 可读不构成物理证明。
