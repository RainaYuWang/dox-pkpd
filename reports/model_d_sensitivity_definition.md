# 模型 D：固定 PK 参考值的探索性敏感性分析

在模型 D 基础上，分别把 kEF、kFE、kFB 设为参考值的 0.5、1、2 倍，每次仅改变一项，并重新拟合跨三剂量共享的 kBF、beta。倍数是人为探索性扰动，不是实验不确定性范围。

| 变化项 | 0.5 倍 | 1 倍 | 2 倍 | 单位 |
| --- | ---: | ---: | ---: | --- |
| kEF | 1.75e-6 | 3.5e-6 | 7e-6 | h^-1 |
| kFE | 0.625 | 1.25 | 2.5 | h^-1 |
| kFB | 0.105 | 0.21 | 0.42 | h^-1 |

VI=0.004629 µL，VE=250 µL，kp=0.031094959444983555 h^-1，theta=35143.615284616222 cells 固定。
数据、目标函数、6 h 暴露及理想 washout、ODE 求解器和精度均保持原样。
参考值来源与解释限制见 [D 定义](model_d_definition.md)。

实际执行 9 个情形，每个 8 起点，共 72 个起点。三个 1 倍情形分别重新拟合，并核对一致性；它们是同一条件的重复计算，不是独立实验重复。
所有情形重置同一随机种子 `314159/twister`，使用同一初始猜测 `[0.005,5e-4]` 和同一多起点生成规则，因此起点集合相同。
使用原 fminsearch、MaxIter=650、MaxFunEvals=1800、TolX=1e-5、TolFun=1e-8，探索边界仍为 kBF `[1e-7,10]`、beta `[1e-9,0.1]`。
每个起点保存 MAT 检查点，记录参数、SSE、退出标记。最佳解和各起点的边界诊断沿用 log 距离 0.02（约 2%）规则。

每剂量在 0–12 h 输出 0.05 h 网格，之后至该剂量最末观测时刻输出 501 点，并加入 6 h washout 和全部观测时刻。
输出 CF、CB、N 和有效杀伤速率 `E(t)=beta*CB(t)`；这只是现有方程项的记录，没有增加机制。

使用新 `sensitivity` 模式：`run_joint_pkpd(8,true)`，在既有 A/B/C/D 与 D profile 运行之后执行敏感性分析。quick/full 原入口保持可用。
GitHub Actions 可手动选择 sensitivity，或推送相关代码至 `pkpd-sensitivity` 自动启动；结果沿用已有 Git 回传流程。
主要结果为 `model_d_sensitivity_summary.csv`、`model_d_sensitivity_starts.csv`、`model_d_sensitivity_predictions.csv`、`model_d_sensitivity_residuals.csv`、`model_d_sensitivity_by_dose.csv` 和完整 MAT。

拟合曲线稳定不意味着参数唯一或固定参考值真实。特别是在线性 PK、零初始胞内药物下，kEF 的倍数缩放可以由 beta 的反向缩放精确补偿，使 beta*CB 及细胞数保持相同（只要自由参数边界允许）。实际重拟合将核对这一结构性质，并检查 kFE/kFB 扰动下的早期瞬态与长期细胞数是否有不同敏感性。
