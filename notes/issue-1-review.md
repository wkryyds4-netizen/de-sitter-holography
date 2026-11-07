# GitHub issue #1：本地核查

2026-11-06。用户要求查看新issue；本轮只读与本地核查，未发GitHub回复/邮件，未更改issue状态。

来源：[wc72449 的 issue #1](https://github.com/wkryyds4-netizen/de-sitter-holography/issues/1)，创建于08:58:33 UTC，快照时0评论，针对固定commit714d038bb605cfcffdc351e4dcc62aa95237d25f的KS五边形支线。PDF与ZIP原样缓存于sources/external-reviews/issue-1。

## 已核与判断

这是有实质内容的外部复核，支持已有正残差，并增加整体相位优化后的检验。与当前Susskind记忆核推导无直接冲突。

直接回读[KS原文v2](https://arxiv.org/html/2112.13329v2)普通谱演算定义、Prop6.11证明(74)–(75)、Prop6.15的完整Weyl假设。标准表示(X,Y,Py,Px)、M=exp(iXY/π)和最终词AB−BM⁻¹BMA的对应成立；归一Gaussian为最终有界算符的合法测试态。

强对易异议成立：N=X−iτY、T=Py+iτPx形式对易，但其有界Borel函数exp(iX)、exp(iPx)反对易，commutator范数在u上为2||u||。移动极点的坐标、Jacobian及非零Gaussian局部L²对数发散也支持固定向量全纯性的异议。这与本地flip-domain、weak-continuation和representation-interface记录相符。式(75)应为e^(πih)；补i不解决另外两个问题。未另行验证附件所报原始TeX行号，定位以已读HTML段落为准。

附件8页全部读完并渲染检查。ZIP有44成员，manifest绑定43文件全部hash通过；18项上游code/data与本地完全一致；PDF独立下载与ZIP内副本hash一致。外部代码只读未执行。

新写code/review_issue1_evidence.py，以Fraction精确提升binary64系数、整数isqrt给320bit有理平方根界，从五个重叠及尾界、两个chirp端点独立重组七步误差。7个原保存上界的平方均覆盖独立组合的平方误差上界；最终使用本机重新组合的保守界，确认：

- ||Df|| > 0.04489650872。
- infθ ||ABf−exp(iθ)BM⁻¹BMAf|| > 0.04489650465。

精确有理端点在results/issue-1-review/independent-assembly.json。上述是对所提供区间的独立算术核查，仍以积分包围及误差估计正确为前提；本轮未重新积分21143个单元。旧本地下界约.04489650872660254，外部约.04489650872992831，支持相同保守阈值。相位优化结果排除用常数整体相位消除该例的残差。

## 9890项失败的范围

附件replayed-audit.json明确passed=false，共63551检查、9890失败。逐项看字符串只有bound类：Gaussian255、left1 256、left2 3016、right3 2969、right5 3394，与其portability摘要一致。摘要报cell kind未变。这是对附件记录的检查，不是本机重新执行全部audit。

本地code/audit_tensor_chain.py第74行要求新error区间与旧error端点区间重叠；旧error为向上端点singleton，因此是很强的回归一致性条件。两个有效保守上界未必端点一致。失败不能等同9890个积分错误，也不能删去后当作通过。应分别报告旧端点重现失败与新上界重组有效性；底层积分例程/估计的正确性还需单独核查。具体跨平台差异原因尚未定位，不擅自归因。

## 保存与后续

三页正式稿manuscript/issue-1-review/main.tex与output/pdf/de-sitter-issue-1-review.pdf已双遍编译、95dpi逐页检查。第一次长hash排版溢出日志保留。

完整独立验证还需新输出目录重算所有单元及chirp，审查Gauss/Mehler/产品尾估计；不能运行快速assemble后称完成全量复现。历史audit与证据未修改。主线继续六个共同J有限记忆产品及体介子字典。未commit/push/发布回复。
