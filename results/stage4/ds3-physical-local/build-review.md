# 编译和核版

2026-10-31，Alan Day / Chengheng署名及AI助手脚注。

- `code/build_ds3_physical_local.ps1`实际执行本机两遍pdfLaTeX，最终exit 0，输出4页A4 PDF。
- 最初受限运行返回exit 1，虽生成了PDF，未计作成功；随后获准运行同一构建脚本，明确exit 0。首次诊断保留在本机tmp中。
- 桌面`compile_latex_document`返回success。最终日志无Overfull、undefined、Fatal或重跑提示。
- `pdftoppm -scale-to 1100 -png`实际渲染4页，并逐页查看标题/作者脚注、精确物理幅度、Jacobian、C²证明、局部定理、表格及参考文献；无裁切、重叠或缺字。
- 输入TeX、构建脚本、最终PDF和log的SHA-256绑定于本阶段evidence.json；源文件未只留作编译计划。

渲染检查不替代科研正确性。当前数值是85位权重计算加复用60位旧U节点，不能将85位参数算术说成全部核的85位精度。
