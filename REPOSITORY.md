# Git保存范围

本仓库保存研究代码、数值输出、文献笔记、LaTeX阶段稿和实际编译的PDF。当前成果均为研究草稿，适用模型、近似和未完成事项见各稿件；没有以仓库发布宣称完整de Sitter全息对偶已经建立。

项目约定中的私人交流设置、研究动机和过程日志保留在本机，不进入Git；AGENTS.md、PROJECT.md、journal及完整恢复快照由.gitignore排除。README/PLAN/NEXT中的部分历史链接指向这些本机文件。第三方论文全文缓存也不再分发，版本与来源链接保留在sources/manifest.json及literature。

最新阅读入口为output/pdf中的阶段稿，对应源码在manuscript。当前可直接复核的阶段脚本为code/audit_ds3_uniform_local.py；需安装requirements.txt及该阶段study.json注明的依赖。code/build_ds3_uniform_local.ps1说明实际PDF编译方式。

完整项目审计还核查本机私人日志和第三方原文缓存，因此仅克隆Git内容不等于取得全部本机归档。不可把缺失缓存的审计失败绕过或说成已验证；应按原始版本取得相关文件后核查哈希。数值舍入、求积误差与解析谱尾仍须区分。

未设定远端时，Git提交仅保存在本机。远端发布的目标与可见性以用户明确指定及实际推送结果为准。

2026-10-31：用户明确选择公开，仓库为 https://github.com/wkryyds4-netizen/de-sitter-holography ，默认分支main。初始6d191c3推送成功并核远端哈希。`.gitattributes`关闭自动换行转换，保留证据字节；后续克隆亦须按原文件哈希核对。最新阶段为ds3-physical-local，独立核查脚本为code/audit_ds3_physical_local.py。公开保存工作稿不代表完整对偶已经建立或完成同行评审。

发布复查修正：一份KL论文PDF及两张原文页图位于results/stage3/cg-domain-source-audit，最初未被sources规则覆盖。已扩充.gitignore，将三文件从公开main初始提交移除；本机原文和原始提交保留。科学证据中的原文哈希仍保留，复算者应自行从来源取得文件。
