# 六页阶段稿编译与检查

2026-11-01，Alan Day / Chengheng，AI助手脚注保留。正式源为 manuscript/flip-domain/main.tex，输出 output/pdf/de-sitter-flip-domain.pdf；哈希由同目录 evidence.json 绑定。

- 实际执行 code/build_flip_domain.ps1，两遍pdfLaTeX均exit0，完整main.log已保存；最终无undefined、rerun、overfull或fatal告警。
- Codex桌面 compile_latex_document 返回kind=success。open_in_codex返回queued，只记录排队，不把它等同于已显示。
- pdfinfo确认6页A4、pdfTeX1.40.29。pdftoppm -scale-to 1100 -png生成六张全页PNG，已逐页查看正文、公式、条件表、数值表、引用与署名，无裁切/重叠。
- 本轮未发生LaTeX源编译失败。首次JS调用因源字符串中反引号而语法失败，未创建项目；移除该字符串冲突后create_latex_project成功，科学文本不受该工具调用错误影响。
- 第2–3页的域论证只针对明确表示；第4页给文献适用边界，第5页明确正margin尚未认证，第6页保留两个初始检查失败。整体状态不是五边形反例或完整全息对偶。

原文核版另检查FG PDF第38页、Kim pseudo-Kähler第21页、Pusz–Sołtan第3页，图像仅存tmp，不公开。
