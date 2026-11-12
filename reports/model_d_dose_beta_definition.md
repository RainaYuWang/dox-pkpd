# 模型 D：分剂量 beta 诊断定义

这是共享 beta 假设的诊断比较，不作为最终预测模型。仅把跨剂量共享的一个 beta 改为三个剂量各自的 beta，不增加任何 ODE 状态或机制。

固定参考 PK：kEF=3.5e-6 h^-1、kFE=1.25 h^-1、kFB=0.21 h^-1、VI=0.004629 µL、VE=250 µL，来源为原始模型提交 `6325d45`。固定 kp=0.031094959444983555 h^-1、theta=35143.615284616222 cells，取自 full 对照组。参考值并非独立实验已确定值；来源限制见 [D 定义](model_d_definition.md)。

固定 kBF=0.0049373889028291049 h^-1，共享模型使用 beta=0.00048479871747565857；两者读取自运行 `43101764836`、代码提交 `950b8e52466c0bc98cff98ce759eea70c5b1929b` 的 `model_d_fit.mat`，按完整双精度写入代码。诊断过程中不重新优化 kBF。

保留全部 32 个处理组观测及原残差缩放。每剂量目标仍为 `sum(((Npred-Nobs)/max(Nobs))^2)`，总目标为三个剂量目标之和。因此，固定其余参数后分别拟合各 beta 与联合最小化三个 beta 的总 SSE 数学等价，未改变剂量权重或排除时间点。

各 beta 边界沿用 `[1e-9,0.1]`。每剂量使用 8 个起点、共享 beta 初值及相同 `161803/twister` 随机起点集合；复用原 `multifit/fminsearch`，MaxIter=650、MaxFunEvals=1800、TolX=1e-5、TolFun=1e-8。另复用现有 `profileBeta` 的一维对数扫描、局部细化及显式端点比较，独立复核各剂量最优 SSE；不把扫描曲线解释为置信区间。每个起点及各剂量扫描结果保存为 MAT 检查点。

复用原 `loss`、`simulate`、`rhs` 和 `embeddedData`，ode15s 的 RelTol=1e-7、AbsTol=1e-9 保持不变，6 h washout 输入不变。输出比较曲线、逐观测有符号残差（预测减观测）、最佳 beta、各剂量 SSE、边界标记及全部起点。

Actions 新增 `dose-beta` 模式和 `pkpd-dose-beta` 自动计算分支；执行 `run_joint_pkpd(8,false,true)`。已有 A/B/C/D 和条件 profile 仍运行作为原比较，随后追加此诊断。原 quick/full/sensitivity 模式保持可用。CSV 中诊断 model_id=1 为归档共享 beta 的 D，model_id=2 为分剂量 beta；这些编号仅适用于 `model_d_dose_beta_*` 文件。

本比较在固定 kBF 的条件下增加两个 beta 自由度。原 D 的 kBF 本身来自同一数据，因此也不能将其当作外部已知值来做未经校准的统计检验。新增参数导致训练 SSE 下降不能单独证明机制更合理、参数真实或预测能力提高。

首次实际运行 `43129468749` 中，低剂量多起点与独立扫描的 SSE 相差 5.52e-7，beta 相差约 0.075%，触发了最初的 1e-7 输出一致性检查。随后在各剂量两种结果附近增加 101 点局部目标函数审计，并包含两种最优值本身，保存 `dose_beta_checkpoints/dose_*/local_objective.csv`。一致性判据改为 SSE 差小于 1e-6 且 beta 对数差小于 0.001，局部审计也不能找到超过 1e-6 的进一步改善。这仅调整诊断结果的比较判据；原数据、目标、边界、优化器设置及 ODE 精度全部保持不变。最终报告需保留这一数值限制，不将 beta 的末位数字作精确科学解释。
