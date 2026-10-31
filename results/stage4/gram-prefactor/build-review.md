# 编译与页面核查

2026-10-30。五页英文研究笔记，署名 Alan Day 与 Chengheng；Chengheng 脚注为 AI research assistant (OpenAI)。

初次桌面编译失败：四段质量坐标使用split环境却含两组对齐列，触发Extra alignment tab；另有一处过长行。桌面先尝试的MiKTeX latexmk缺少Perl，回退Tectonic后报告了上述源文件错误。这里保留失败经过，不声称保存了那次完整原始日志，也不是自动审批阻断。

将坐标改为aligned，将Gram因子独立成行之后，桌面compile_latex_document返回success。随后本地build_gram_note.ps1以带引号的绝对输出目录调用pdfLaTeX两遍成功。最终main.log没有Overfull、未解引用或重跑要求，输出5页。阅读版来自本地两遍构建。

最终PDF以pdftoppm、95 dpi渲染全部5页，逐页目视检查：标题、双方姓名、脚注、矩阵、长分式、定理证明、页码和三条参考文献均可读；未见裁切或重叠。渲染图的SHA-256写入evidence.json。不是仅依靠文本提取判断排版。

二维55页稿和WKB四页稿未改。未对外发表。科学内容为自查草稿；页面与编译检查不认证证明或有限G对应。
