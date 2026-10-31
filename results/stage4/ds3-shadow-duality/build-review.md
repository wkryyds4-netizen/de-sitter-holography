# Wilson对偶阶段稿：编译与核查

2026-10-31，4页，Alan Day与Chengheng署名，AI研究助手脚注。沿用已读取PDF技能及本会话既有marker，不重跑旧marker。

create_latex_project创建独立目录；code/build_ds3_shadow_duality.ps1完成本机双遍pdfLaTeX，无Overfull、undefined或重跑要求。桌面compile_latex_document明确success；open_in_codex返回queued。

最终四页渲染最长边1100像素并逐页查看。标题、Wilson定义、反射公式、核误差界、失败说明、数值表与参考文献无裁切、碰撞、缺字或重复编号。没有强取实部处理数值；表格明确只显示实部，完整复数保存。

初始33项检查含两项失败，四项新精度检查修复但不覆盖旧结果。365项独立保存值审计通过，18组完整核结果保留。解析尾不含舍入，M128舍入下限明确。旧47841仍活跃，不计入已完成验证；三个新增阶段尚待整体登记和冻结。
