# SUM-149PT：联合调整 PK 参数的直接 PK–PD 模型

## 运行

解压后在 MATLAB 打开 `run_joint_pkpd.m`，点击 Run。
也可以在当前文件夹输入 `results = run_joint_pkpd;`。
需要基础 MATLAB（建议 R2019a 或更新版本），不需要 Optimization Toolbox。
默认每个主要模型使用 8 组起点；命令行逐次显示拟合进度。
`run_joint_pkpd(2)` 可用于较快试运行，但较少起点更容易漏掉好的拟合。
数据内嵌在文件最后的 `embeddedData()` 中，不需要另找 CSV。
`verify_equations.py` 只是开发时的独立数值验证，使用 MATLAB 不需要运行它。

云端手动 quick/full 流程、结果下载及日志说明见 [ci/README_zh.md](ci/README_zh.md)。
MATLAB 图窗默认不可见，仍保存所有图片；每完成一个优化起点即保存中间 MAT 文件。

## 方程和固定量

CF' = kEF (VE/VI) CE - (kFE+kFB) CF + kBF CB

CB' = kFB CF - kBF CB

N' = kp N (1-N/theta) - beta CB N

不包含 damage-memory。杀伤项没有 logistic 因子；这一点与论文原人口方程不同。
CE 为已知矩形输入：0–6 h 保持给药浓度，6 h 后为零。
假设理想 washout，没有再摄取胞外残留。换液时胞内状态连续，所有 k 不变。
初始 CF=CB=0，每条曲线的 N0 使用其 t=0 数据点，不使用接种细胞数。
数值求解 CF、CB、log(N)，与上述三个物理方程等价，便于保持 N>0。

固定 SUM-149PT 的 VE=250 µL、VI=0.004629 µL、kEF=3.5e-6 h^-1。
旧参考 PK：kFE=1.25 h^-1、kFB=0.21 h^-1、kBF=0。
这些数值来自 Project_reproduction(2).ipynb；VI 是旧模型校准的有效体积。
它们不是本次由独立 PK 数据测出的量。固定 kEF 是浓度尺度约定，不能称为已验证。
kp、theta 从对照组拟合后固定，共用于三档浓度。

## 数据来源与适用范围

McKenna et al., Scientific Reports 7:5725 (2017), DOI: 10.1038/s41598-017-05902-z。
对照：Supplementary Fig. S3 第一列 SUM-149PT。
处理：正文 Fig. 4 第一列 SUM-149PT，三行分别标注 156、312、625 nM，6 h 暴露。
使用二倍稀释系列的数值 156.25、312.5、625 nM，与图中 AUC 标签对应。
既有代码将 t=0 作为药物施加时刻，负时间的处理前点未用于拟合。
用户提供完整正文及 supplementary 后，已依据正文第 6–8 页核实原文 t=0 为给药开始；当前数字化零刻度与之对应，6 h washout 不重置时钟。详见 [定义核查](reports/paper_definition_audit.md)；这不改变现有初值、观测或残差定义。

沿用上一版的稀疏手工读图数据，不是新的原始数据。
提取的是黑色测量点的近似中心，不是红色拟合线；误差棒、红色叠加和图中文字
遮挡可能造成读图误差。625 nM 被注释框遮挡的部分点未提取，也未插值补点。
两张 `selected_points.png` 用蓝圈标明提取位置；`embeddedData()` 保留像素坐标
和坐标轴换算。正式报告数值前应复核读图。此代码用于初步机制探索。

## 四种比较

A. Old PK / direct PD：旧 PK 三个速率固定，kBF=0，只拟合 beta。
B. Joint irreversible：kBF=0，联合拟合 kFE、kFB、beta。
C. Joint reversible：联合拟合 kFE、kFB、kBF、beta。

D. Fixed reference：固定原参考 PK 和已归档 full 生长参数，仅拟合共享 kBF、beta。
参考值来源与条件 profile 定义见 [模型 D 说明](reports/model_d_definition.md)。

每个模型都用同一批细胞数点和相同误差定义。A/B/C 使用本次对照拟合的 kp/theta；
D 固定使用已归档 full 的 kp/theta。本次 full 中两者完全一致，quick 中可能有微小数值差异。
同一个模型中三档剂量共用参数，不给每档剂量各配一套。
A/B 都不是论文原来的完整 PD 模型；所有比较使用我们新的直接杀伤方程。
C 比 B 多一个名义参数，但可观测的独立组合数也需考虑；不根据误差降低直接
宣称机制正确，不自动给出 AIC、显著性或置信区间。

目标函数：对每条曲线，用该曲线最大观测细胞数除残差，再将所有残差平方相加。
没有使用误差棒权重，没有传播 N0、VI、固定 PK 或生长参数的不确定性。
所有数据都参与拟合，没有独立验证集。
旧 code 的 CF+CB 轨迹只绘图对照，不作为额外的实验约束、不进入目标函数。

## 为什么联合拟合也不能得到唯一的四个参数

令 E(t)=beta CB(t)，在初始胞内药物为零、常参数和已知输入下，可以从两个 PK
方程消元得到：

E'' + a E' + b E = c CE(t)

a = kFE+kFB+kBF
b = kFE*kBF
c = beta*kEF*(VE/VI)*kFB

N 方程只通过 E(t) 感受 PK，所以细胞数轨迹只能通过这些组合区分 PK/PD。
即使连续无噪声地观测 N，固定 kEF 也不能把四个自由参数全部拆开。
这修正了之前“固定 kEF 后其余参数只是可能有相关性”的不完整表述：
这里存在精确的结构性不可辨识，不只是有限数据带来的不确定性。
对照组 kp/theta 已固定；可识别组合在实际稀疏数据中也不保证精确估计。

保持 a,b,c 不变，选一个新的 kBF=x，则：

kFE = b/x
kFB = a-x-b/x
beta = c/[kEF*(VE/VI)*kFB]

只要这些值为正，就产生相同的 E(t) 和 N(t)，但 CF、CB 可能不同。
代码在同一探索边界内生成几套这样的等价参数，直接画出不同 PK、相同 N。
因此 `representative_parameters.csv` 中每一行都是一个拟合代表，而非唯一的
生物学估计。不能将其中某个 kBF 单独解释成真实解离速率。
旧版在固定 kFE/kFB 时得到的 kBF profile 是更强假设下的条件性结果；本版释放
这些参数后，不能沿用旧版的窄谷来主张 kBF 被唯一识别。

## 输出：joint_results 文件夹

- 01_cell_fits.png：对照及三个浓度，三种模型对比。
- 02_residuals.png：各模型随时间的残差，用于发现系统性偏差。
- 03_PK_reference.png：旧模拟 PK 与新拟合代表的 CF+CB、CB；不是实测 PK。
- 04_equivalent_parameters.png：多套参数的胞内 DOX、杀伤率、细胞数对照。
- representative_parameters.csv：三个模型的代表参数及 scaled SSE。
- growth_parameters.csv：kp、theta。
- equivalent_parameter_sets.csv：相同细胞数响应对应的其他参数组合。
- reversible_combinations.csv：a、b、c。
- *_multistart.csv：各起点的参数、误差、退出标记和函数评估次数。
- predictions.csv / residuals.csv / fit_by_dose.csv：分浓度的详细结果。
- search_bounds.csv：人为设置的探索边界，不是生物学置信范围。
- joint_fit.mat：保存 MATLAB 结果和本次使用的数据。

模型编号：1=A，2=B，3=C，4=D。
不保证全局最优；exitflag=1 也不证明唯一性。触及搜索边界要检查。
如果可逆模型明显比重新拟合的不可逆模型更差，先检查优化收敛，而不是直接
解释为生物学结论；代码加入接近不可逆的起点并检查这种情况。

## 验证情况

Codex 容器本身未安装 MATLAB；已在 GitHub-hosted Linux runner 上实际完成 MATLAB quick，
并自动取回、校验和检查结果，见 [云端验证记录](ci/VALIDATION_2026-11-11.md)。
随包的 verification_report.txt / verification_preview.png 来自独立 Python
求解相同三个 ODE，用于检查数值、物质转移恒等式、无药效时的 logistic 极限，
以及等价参数的细胞数轨迹一致性。优化器不同，代表参数可能不同。

本次探索的目标是判断直接动态连接能解释什么行为，以及哪里仍有系统性偏差。
读图误差、共享 theta、线性杀伤、简化 PK 等都可能影响拟合，不能将所有偏差
自动归因于缺少 damage-memory。

模型 D 的实际运行、D/C 对照及条件 profile 见 [D 结果报告](reports/model_d_results.md)。

已完成原文时间形式 A/B 的独立 MATLAB quick/full，保留旧 D、分剂量 beta 和延迟结果而未重跑。A/B 的训练 SSE 为 1.540788 / 0.434762；B 仍比分剂量 beta 的 0.323945 高，625 nM 早期由过低变为过高。原文定义、各剂量参数、残差、局部图和自动取回命令见 [时间形式诊断结果](reports/time_form_results.md)。
