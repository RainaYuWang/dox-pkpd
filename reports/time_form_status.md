# 原文时间形式 A/B 诊断：准备阶段记录

**当前已完成：原文定义核实、实际 MATLAB quick/full、结果自动回传与分析。** A/B full SSE 分别为 1.540788 / 0.434762；完整运行证据、参数和局部残差见 [A/B 结果报告](time_form_results.md)。旧 D、分剂量 beta 和延迟模型未重跑。

**以下全文为收到 PDF 前的准备阶段历史记录；其中“未核实/未运行/待提供”描述当时状态，不代表当前状态。** 当前原文核查见 [定义记录](paper_definition_audit.md)。

本轮核对远端 main=`68f54bc` 和已有归档后，没有重跑 D、分剂量 beta 或延迟模型。以下仅为已保存的旧结果，不是本轮新拟合：

| 比较项 | 已有训练 SSE | 本轮重拟合参数数 |
| --- | ---: | ---: |
| D | 0.99396925495 | 0 |
| 分剂量 beta | 0.32394513275 | 0 |
| 延迟阶段 1 | 0.99395110667 | 0 |
| 延迟阶段 2 | 0.99393720526 | 0 |
| 常数时间形式诊断 A | 未运行 | 计划 3 |
| 暂时时间形式诊断 B | 未运行 | 计划 6 |

本轮没有生成 A/B 的参数、SSE、细胞数拟合曲线或正式科学比较。B 是否优于 A、是否缓解 625 nM 早期过度下降、是否改善 156.25 nM 晚期高估和分剂量 beta 的剩余时间残差，目前均未作结论。

## 原文获取状态与所需材料

当前网络策略仍是已执行的 revision 1，未应用此前允许文献域名的草案。本轮没有重复请求此前返回 CONNECT 403 的入口。新找到的公开仓库 [mckennamt/DoxorubicinPKPD](https://github.com/mckennamt/DoxorubicinPKPD/tree/76952a6b8e6cb412ee0081b89ca1ccba61a00537) 明确引用目标论文，也有相关实现线索，但没有正文或补充材料。详细检索、固定版本链接及代码证据的适用限制见 [原文核查记录](paper_definition_audit.md)。

仍待核实：原文 t=0；D 与浓度/暴露时长/bound-drug 指标的关系；theta(D)；初始细胞数及时间；拟合区间和处理/washout 数据使用方式。不能根据实现线索补全这些原文定义，`paper_time_definitions.json` 仍为 `verified=false`。

需要提供 McKenna et al. (2017) **“A Predictive Mathematical Modeling Approach for the Study of Doxorubicin Treatment in Triple Negative Breast Cancer”**，DOI **10.1038/s41598-017-05902-z** 的完整正文 PDF，以及同一文章页面的 Supplementary Information PDF/附件（含补充 Methods、模型拟合说明与 Fig. S3）。当前没有声称“草案保存即访问恢复”。

## 已完成的独立代码准备

`run_time_form_diagnostics(nStarts,baselineDir)` 是单独入口。它只拟合 A/B，并读取归档 run `43137925939` 的数据、生长参数、D/分剂量 beta/延迟预测与残差。它不会调用 `run_joint_pkpd`、旧 PK 优化器入口或重新积分旧模型来生成比较曲线。MATLAB 优化帮助函数从现有项目逐字复用，旧入口及科学实现未修改；CI 总入口仅添加 time 模式分派。

基线锁定源码 `f05543e978e425785fdf918751185cd1a733d300`、run `43137925939`、attempt 1 和压缩包 SHA-256 `b39c64c22d8ebe1e9fd6bcd07d42c478c2b105e0d7cfdb1df8f8d1c2d5f593d3`。`ci/prepare_time_form_baseline.py` 先检查原文记录，然后通过现有 Git collector 取回并核对该包。未核实时在任何下载或 MATLAB 计算之前停止。

保留 `N'=kp*N*(1-N/theta)-kd*N`，N0 从每条归档数据首点及其原始时刻初始化。生长参数固定为原 full 对照估计 `kp=0.031094959444983555 h^-1`、`theta=35143.615284616222 cells`；所有 32 个观测和按各剂量最大观测值缩放的残差完全保留。效应原点之前仍积分 logistic 并保留观测，在原点延续相同 N，不增加拟合初值或重置状态。

代码区分 `time_h`（当前数据原坐标）、`time_since_dose_h`（减去经核实的给药时刻）、`tau_h`（减去经核实的效应时刻）。这两个时刻都尚未填入。未来需在 `source_evidence` 中为五项定义分别记录 source、location、excerpt；格式完整检查本身不等于文献科学核实。

A 的三个剂量各拟合一个非负 kd_A；B 各拟合 kd_B、r，共六参数。均为 h^-1。A 在效应起点之前为零、之后为常数；B 在 tau≥0 为 `kd_B*r*tau*exp(1-r*tau)`，之前为零。这两个诊断没有动态 CB 输入、不作为跨剂量预测模型、不引入 theta(D) 或把 logistic 因子乘到杀伤项，也不做 AIC 加权。

| 新参数 | 正参数搜索区间 | 依据 |
| --- | --- | --- |
| kd_A、kd_B | `[1e-9,1] h^-1`，另显式比较 kd=0 | 下端接近零效应，上端远大于对照增长率约 0.0311；探索范围，不是论文或实验区间 |
| r | `[1e-4,100] h^-1` | B 的峰时 1/r 覆盖 10000–0.01 h，包含相对数据窗口很慢与很快的效应 |

沿用现有对数参数变换、fminsearch MaxIter=650、MaxFunEvals=1800、TolX=1e-5、TolFun=1e-8 和 ode15s RelTol=1e-7、AbsTol=1e-9。每个形式/剂量分别固定 RNG=`57724+100*form+dose_index`、twister；quick 每组 2 起点（共 12），full 每组 8 起点（共 48）。每个起点保存独立 MAT 和日志，保留失败退出、较差局部解、参数差异。边界判据为与搜索端点的对数距离 <0.02；若零杀伤获选，B 的 r/1/r 输出缺失并标明 selected_zero，不虚构其确定性。

数值检查代码覆盖零杀伤恢复解析 logistic、B 的峰时 1/r 与峰高 kd_B、原点前零杀伤、N 连续及非负性。另用常数杀伤解析解检查效应原点位于初始观测之前、窗口内、窗口之后的三种情况。检查使用显式标注的合成时钟，不给实验数据赋予任意时钟。没有人为细胞数截断。

这些 MATLAB 检查目前**尚未实际执行**，不能把代码存在或静态检查通过称为数值运行成功。

## 工作流和后续执行顺序

新增 `time-quick` / `time-full` 和计算分支 `pkpd-time-quick` / `pkpd-time-full`，沿用官方 MATLAB Actions、超时、失败上传以及 Git 回传。main 更新不触发计算。本轮只推送代码准备到 main，不触发 A/B 或旧模型运行。

原文定义实际核实并填写后，按以下顺序进行；下面是待执行流程，不是已完成记录：

```bash
python3 ci/prepare_time_form_baseline.py --check-only
git push origin HEAD:refs/heads/pkpd-time-quick
python3 ci/collect_results.py --sha HEAD --mode time-quick --wait 5400
# 检查实际 MATLAB quick 成功及数值/输出结果后，再运行 full：
git push origin HEAD:refs/heads/pkpd-time-full
python3 ci/collect_results.py --sha HEAD --mode time-full --wait 5400
python3 ci/analyze_time_forms.py --results <full的joint_results路径> --output reports/time_forms
```

工作流在安装 MATLAB 前执行原文与基线预检查；预检查或拟合失败仍尽量回传已有日志。已有基线通过 Git 获取，不依赖受限的 Actions REST API。只有手动创建全新计算分支时，应先确认分支状态与路径过滤，避免把预检查失败误当作已完成 MATLAB quick。

计划输出包括：A/B 六组参数、总/分剂量 SSE、所有起点和边界标记、原始/给药/效应时间、细胞数与逐时残差、有效杀伤率、1/r、数值检查、图表及 MAT。`ci/analyze_time_forms.py` 只接受状态成功且原文 verified 的新模式结果，再生成六模型比较、156.25 nM 晚期和 625 nM 早期局部图、指标 CSV 与 JSON；没有预先生成虚拟结果。

## 本轮实际验证及限制

- MATLAB 静态解析：新入口、独立 CI 包装器和模式分派均通过 MISS_HIT。
- Python AST、工作流 YAML、所有 push/手动模式映射检查通过；main 不在计算触发分支中。
- 独立入口使用的 multifit/boundedLoss 与旧代码逐字一致。
- 实际执行 `python3 ci/prepare_time_form_baseline.py --check-only` 返回预期退出码 2，列出全部缺失定义与出处，阻止未核实拟合。
- 未执行 A/B MATLAB quick/full，未重跑旧模型，未产生新训练结果或预测验证结论。

新增文件：`run_time_form_diagnostics.m`、`ci/run_time_form_ci.m`、`ci/prepare_time_form_baseline.py`、`ci/analyze_time_forms.py`；修改 workflow/CI 分派以运行独立模式，添加基线缓存忽略项，更新原文来源记录、状态报告和入口文档。待提供原文或访问实际恢复后，剩余工作是文献核实 → quick 数值检查 → full → 自动分析，而不是重新执行已完成的延迟分析。
