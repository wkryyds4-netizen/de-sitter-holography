# 编译与核图

2026-10-30，最终4页，署名Alan Day与Chengheng，后者附AI research assistant (OpenAI)脚注。

第一次本地编译因mathscr未定义失败；加入mathrsfs后，桌面编译成功但本地缺mathrsfs.sty。两次本地失败日志分别保存在initial-build-failure.txt、missing-package-build-failure.txt。随后改用已有amssymb提供的mathfrak并移除缺失包，无需安装依赖。最终本地两遍编译及桌面编译成功。

初步五页的第五页只有一条参考文献。压缩数值检查段的重复叙述，并保留34项、精度、抵消率和限制，最终4页。11pt字号、25mm页边距。最终main.log无Overfull、undefined引用、重跑要求或Fatal error。

全部最终页面以95dpi渲染并逐页核查。核对首页双方姓名与脚注、两端界、归一化负号、Barnes积分、复杂分式及全部引文，无裁切重叠。最终渲染哈希写入evidence.json。页面检查不认证数学或物理结论，旧PDF及科学记录保持不变。未对外发表。
