# 原文定义核实状态

**当前状态：用户已提供完整正文和 Supplementary Information，五项定义已直接核实。** 原文 t=0 为给药开始；本轮 `tau=当前数据时间`，washout 不重置时钟。完整文件和 SHA-256 保存在 [文献目录](../references/papers/README.md) 及 [定义 JSON](paper_time_definitions.json)。以下未访问/缺失的段落保留为历史记录，已不再构成本轮运行阻塞。

| 项目 | 原文依据与结论 | 本轮诊断处理 |
| --- | --- | --- |
| t=0 | 正文 pp.6–7 Eq.(5A/B) 与相邻段落明确 treatment at t=0，p.8 明确 doxorubicin at t=0 | 给药时刻与效应时刻均映射到当前 Fig.4 数据零刻度；tau=t，不在6 h重置 |
| D | 正文 p.6 定义为综合浓度与暴露时长的处理条件；p.7 以模拟 CB,max 及外部 AUC 描述 | 三个6 h浓度组各自拟合，不把 D 等同于 nM 或动态 CB(t)，不引入模拟 PK 作为数据 |
| theta(D) | 正文 p.7 先拟合对照 kp/theta，再固定 kp、各处理拟合 theta；SI p.4 给出相对对照 theta 的0.7–1.3倍边界 | 仍固定本项目已有 kp、theta，不引入 theta(D)；不照搬原文该边界 |
| 初值 | SI p.3 使用每个重复的首个处理后计数；正文 p.8 的预测也从首个处理后计数开始 | 保留现有三组在当前 t=0 的首点 N0，不拟合/移动/重置；原文未给统一的首个处理后观测小时，不能擅称为washout时刻 |
| 拟合区间 | 正文 p.5 给药前≥3天、给药6/12/24 h、媒体替换后持续每日观测≥30天；p.7区分对照与处理后拟合；SI p.3包含初始处理后计数及随后各点；SI p.8控制数据达平台后截断 | 保留既有32点和各组0–436.36/414.55/533.33 h窗口；当前无负时间观测，不补点、不换拟合区间 |

原文 Eq.(4) 的 logistic 因子作用于净增长率 `(kp-kd)N`，使用处理依赖 theta、逐点相对残差及正则化，并进行 AIC 加权；本轮遵循用户要求使用生长/杀伤分离方程、固定生长参数及既有按剂量最大值缩放的 SSE。故名称仍为“借鉴原文时间形式的诊断对照”，不是完整原文复现。SI pp.3–4 原文 r 的边界为0.001–0.05 h^-1（峰时1000–20 h）；本轮沿用已说明的更宽探索边界0.0001–100 h^-1，明确不把它说成原文边界或实验不确定性范围。

## 以下为原文取得前的历史记录

目标 DOI：10.1038/s41598-017-05902-z。

本次已检查仓库和工作区，没有找到论文正文或 supplementary 文件。仓库有 Fig. 4 / Fig. S3 的选点图片、内嵌数字化数据和 README。README 原有“t=0 按 Methods 中药物施加时刻解释”的断言没有附上可核查的原文摘录，不能用来替代本次原文核实，也不能凭选点截图补全定义。

尝试访问 `https://www.nature.com/articles/s41598-017-05902-z` 和 Europe PMC 的 DOI 查询接口，均在当前云代理的 HTTPS CONNECT 阶段返回 HTTP 403（server: envoy）。因此尚未读到原文或补充材料正文，不把该网络错误当作论文不存在或不公开。

| 需核实的定义 | 当前状态 |
| --- | --- |
| 原文 A/B 的 t=0：给药开始、washout 或其他原点 | 未核实 |
| 原文 D 的定义、单位及与暴露时长的关系 | 未核实 |
| theta(D) 的具体形式及拟合方式 | 未核实 |
| 初始细胞数采用哪个时刻、是否为拟合参数 | 未核实 |
| 拟合区间和处理前观测的处理 | 未核实 |

用户提供的 A/B 形式已实现为可选诊断，但 `paper_time_definitions.json` 的 verified=false 阻止它们进入拟合。只有核实原文并写明原点相对当前数据轴的位置后才启用；当前不默认给药开始，不输出伪造的 A/B SSE，不称为论文复现。

已保存限于相关出版/文献域名的网络配置草案，保留原有预设；平台要求在环境设置中审核、保存并发布，草案本身尚不改变正在运行的网络策略。相关域名为 `www.nature.com`、`nature.com`、`static-content.springer-cdn.com`、`www.ebi.ac.uk`、`europepmc.org`、`pmc.ncbi.nlm.nih.gov`。不需要任何 token。

五状态延迟部分继续沿用用户指定和当前模型既有的输入：从首次给药开始模拟，CE 在 0–6 h 为给药浓度、随后归零，初始 CF/CB/E1/E2 为零，N0 使用当前数据首点。当前数据点和时间坐标不改动；该输入约定的原文依据仍待核对。该部分不依赖未核实的原文 A/B 时间形式定义。

## 后续 A/B 专项核查：新增公开来源，定义仍未核实

本轮首先 fetch 并核对 main=`68f54bc`，已有 D、分剂量 beta、延迟 quick/full 结果均在。没有重复运行任何旧模型。再次检查当前执行环境：observations_current=true、network policy=enforced、revision=1，自定义允许域名仍为空。此前保存草案未改变实际访问，本轮没有重复请求已被拒绝的 Nature/Europe PMC 入口。

当前没有可用的文献检索连接器。通过允许访问的 GitHub 网页进行公开仓库检索，HTTP 200：DOI/DOI+in:readme 无命中，SUM-149PT+in:readme 仅命中本项目，McKenna+doxorubicin+in:readme 找到 [mckennamt/DoxorubicinPKPD](https://github.com/mckennamt/DoxorubicinPKPD)。已只读检出版本 `76952a6b8e6cb412ee0081b89ca1ccba61a00537`，其 [README](https://github.com/mckennamt/DoxorubicinPKPD/blob/76952a6b8e6cb412ee0081b89ca1ccba61a00537/README.md#L5) 明确引用：

> McKenna MT et al. “A Predictive Mathematical Modeling Approach for the Study of Doxorubicin Treatment in Triple Negative Breast Cancer.” Scientific Reports. 7:5725 (2017).

该公开仓库没有正文、supplementary PDF、HTML/XML 正文或 manuscript 文件；它是与目标论文相关的实现来源，而不是已取得的论文正文。以下记录仅陈述这个代码版本的行为，不能替代原文定义核实：

| 实现线索 | 固定版本位置 | 尚不能据此确认的事项 |
| --- | --- | --- |
| 时间轴使用 `times-times(datp)`，给药日程由该时刻与 ExposureTime 设置 | [loadCellLineData.m L207–231](https://github.com/mckennamt/DoxorubicinPKPD/blob/76952a6b8e6cb412ee0081b89ca1ccba61a00537/loadCellLineData.m#L207)；[IJ loader L217](https://github.com/mckennamt/DoxorubicinPKPD/blob/76952a6b8e6cb412ee0081b89ca1ccba61a00537/loadCellLineData_IJ.m#L217) | 原文对 A/B 的 t=0 措辞、当前数字化轴的准确映射 |
| 分别保存 `dox_pss=Yd(end)` 与 `auc=drugAdded*ExposureTime` | [loadCellLineData.m L225–226](https://github.com/mckennamt/DoxorubicinPKPD/blob/76952a6b8e6cb412ee0081b89ca1ccba61a00537/loadCellLineData.m#L225) | 论文符号 D 的正式定义，不能直接把它等同于 nM 浓度、AUC 或 C_B,max |
| 常数形式的参数列表含 kd、theta；暂时形式含 kd、r、theta | [fitEachTreatmentModel_regularized.m L11–37](https://github.com/mckennamt/DoxorubicinPKPD/blob/76952a6b8e6cb412ee0081b89ca1ccba61a00537/fitEachTreatmentModel_regularized.m#L11) | 论文 theta(D) 的精确处理和估计说明 |
| 拟合切片从 `DrugAdded_Tp+1` 开始；每个 well 的 S0 取切片后的首个细胞数 | [fitChemoModel_cellLine_Bootstrap.m L65](https://github.com/mckennamt/DoxorubicinPKPD/blob/76952a6b8e6cb412ee0081b89ca1ccba61a00537/fitChemoModel_cellLine_Bootstrap.m#L65)、[runDamageModel_regularized.m L29](https://github.com/mckennamt/DoxorubicinPKPD/blob/76952a6b8e6cb412ee0081b89ca1ccba61a00537/runDamageModel_regularized.m#L29) / [L106](https://github.com/mckennamt/DoxorubicinPKPD/blob/76952a6b8e6cb412ee0081b89ca1ccba61a00537/runDamageModel_regularized.m#L106) | 正文/补充材料如何说明初值、各条件拟合区间以及处理/washout 数据 |
| 暂时形式归一化峰时为 1/r，另有 logistic 因子作用于总净增长率 | [damageModel_Simplified_TimeMod.m L106](https://github.com/mckennamt/DoxorubicinPKPD/blob/76952a6b8e6cb412ee0081b89ca1ccba61a00537/damageModel_Simplified_TimeMod.m#L106) | 本轮仅借鉴杀伤时间形式，明确不复制完整方程、边界或生长参数 |

未从该仓库导入新的拟合数据、参数或模型代码；当前数据、N0、固定生长条件不变。`verified=false` 保持不变，五项原文定义及两个时间轴映射仍留空。下一步需要下列原始文件，而不是更多代码推断：

1. DOI `10.1038/s41598-017-05902-z` 对应的 **完整正文 PDF**，标题为上述 *A Predictive Mathematical Modeling Approach...*，包含 Methods、模型定义及 Fig. 4 图注。
2. 同一文章页面的 **Supplementary Information 完整 PDF/附件**，包含补充 Methods、模型/拟合说明和 Fig. S3 图注；若这些分在多个附件，请提供对应全部附件。

收到文件或实际恢复出版/文献来源访问后，逐项记录原文页码/节名、短摘录和当前时间轴映射。没有正式原文依据前，不通过设置任意原点来生成 A/B 结果。独立代码准备与本轮执行状态见 [A/B 专项状态](time_form_status.md)。
