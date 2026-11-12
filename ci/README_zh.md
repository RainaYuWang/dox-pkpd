# Linux MATLAB 云端拟合

入口仍为 `results = run_joint_pkpd(nStarts)`，默认 `nStarts=8`。
无绝对路径、外部数据或额外工具箱依赖。图窗不可见，四张 PNG 仍会保存。
ODE、washout 处理、目标函数、数据、参数边界、固定/自由参数、优化迭代上限和容差均保持原样。
随机数生成器显式固定为 `rng(42,'twister')`，串行执行以保持起点顺序。

## 手动运行 GitHub Actions

将 `.github/workflows/matlab-pkpd.yml` 推送至默认分支后：

1. 打开仓库 **Actions → MATLAB PK-PD → Run workflow**。
2. 首次选择 `quick`。它执行 `addpath('ci'); run_pkpd_ci('quick')`，内部调用 `run_joint_pkpd(2)`。
3. 等待 MATLAB 步骤及 artifact 上传结束。只有实际运行成功且输出检查通过，才能认定流程通过。
4. 在该次运行页面底部 **Artifacts** 下载 `pkpd-quick-<run_id>-<attempt>`。
5. 完整拟合仍通过 **Run workflow**，选择 `full`，内部调用 `run_joint_pkpd(8)`。

已登录且有相应权限的 GitHub CLI 也可执行：

```bash
gh workflow run matlab-pkpd.yml --repo RainaYuWang/dox-pkpd --ref main -f mode=quick
gh run list --repo RainaYuWang/dox-pkpd --workflow matlab-pkpd.yml --limit 5
gh run watch RUN_ID --repo RainaYuWang/dox-pkpd --exit-status
gh run download RUN_ID --repo RainaYuWang/dox-pkpd --dir downloaded-results
```

完整拟合将第一条命令中的 `mode=quick` 换为 `mode=full`。
workflow_dispatch 要求工作流首先存在于默认分支；之后可以选择含工作流的其他分支测试。

## 版本、许可与超时

使用 GitHub-hosted `ubuntu-22.04` 和固定发行版 MATLAB `R2024b`（该发行版的安装更新号由官方 Action 选择，实际版本写入日志）。
已核对官方当前文档：

- [setup-matlab@v2](https://github.com/matlab-actions/setup-matlab/tree/v2)：支持 R2021a 及以后版本。
- [run-command@v2](https://github.com/matlab-actions/run-command/tree/v2)：非交互式 MATLAB 执行。
- 公开项目由官方 Actions 自动许可基础 MATLAB；本项目不用 Compiler/Coder，不配置 `MLM_LICENSE_TOKEN`。
- GitHub 官方辅助 Actions 使用 `checkout@v6` 和 `upload-artifact@v7`。

安装上限 30 分钟，quick 拟合及输出检查上限 60 分钟，full 上限 240 分钟，整个 job 上限 300 分钟。
这些是保护性超时，不是已测量的运行时长。实际耗时见日志和 `ci_status.json`。

## 结果与故障定位

`joint_results/` 包含原有 PNG、CSV、`joint_fit.mat`，以及：

- `run.log`：模型运行日志，包含时间戳、每个起点开始/结束、SSE、exitflag、函数评估次数和耗时；长起点每约 60 秒输出一次优化进度。
- `run_inputs.mat`：配置、内嵌数据、起点数量及初始 RNG 状态。
- `checkpoints/*_start_XX.mat`：每完成一个起点立即保存独立文件，包含当前参数、累计完成的起点、最佳结果、边界、优化器输出、耗时和 RNG 状态。
- `ci_status.json`：运行模式、实际 MATLAB 版本、commit/run ID、状态、总耗时及成功时的 SSE。

加入 D 后 quick 共 11 个联合拟合起点（control/A/B/C/D 各 2 个，加 1 个接近不可逆的可逆起点）；full 共 33 个（control/A 各 4 个、B/C/D 各 8 个，加 1 个额外起点）。另逐点计算 D 的条件 kBF profile，保存于独立 profile_checkpoints 目录。
中间 MAT 用于检查和保留已完成工作，当前不自动从中断点续算。

`ci-logs/context.txt` 记录执行命令和提交；`ci-logs/matlab.log` 记录 MATLAB 进程输出，包括错误栈和输出检查。
上传步骤使用 `always()`，失败或拟合步骤超时后尽量保留已有文件，保留期 30 天。
runner 故障、整个 job 超时或强制终止可能使上传无法完成；GitHub 页面上的原始 Actions 日志仍应检查。
优化器未收敛、参数触界及等价轨迹误差会在原有日志中警告；流程通过不代表每个起点都收敛。

本地无图形界面 MATLAB 的等价命令（从仓库根目录执行）：

```bash
matlab -batch "addpath('ci'); run_pkpd_ci('quick')" -logfile matlab-console.log
```

CI 包装器拒绝复用已有 `joint_results`，防止旧文件伪装成本次结果。再次运行前请移动旧结果目录。
GitHub-hosted runner 每次使用全新 checkout。运行输出已加入 `.gitignore`。

## 自动修改 → 运行 → 取回结果

`pkpd-cloud` 是 quick 计算分支，`pkpd-full` 是 full 计算分支；只有 `**/*.m`、`ci/*.py` 或工作流本身改变时才启动对应模式。
`main` 和文档修改不会自动计算；手动 quick/full 仍可用。同一分支只保留一个运行中的任务，后续提交排队（GitHub 只保留最新 pending 任务）。

从仓库根目录，提交代码后执行：

```bash
git push origin HEAD:refs/heads/pkpd-cloud
python3 ci/collect_results.py --sha HEAD --wait 5400
```

全自动 full 使用 `git push origin HEAD:refs/heads/pkpd-full`，随后使用同一个 collector；该分支调用现有 `run_joint_pkpd(8)`，不改变模型或优化设置。

推送仅允许 fast-forward，遇到分支分歧时先检查并合并，不要强推。
collector 每 30 秒查询一次 Git 结果分支；默认不加 `--wait` 时查询一次。
退出码 0 表示本次 MATLAB 输出检查成功，1 表示运行失败并已取回已有结果，2 表示尚未完成或仍在排队。
结果保存在 `cloud-results/<source_sha>/<run_id>-<attempt>/`。包含 `manifest.json`、经过 SHA-256 核对的 `results.zip` 及解压后的输出。
按精确代码提交和运行编号匹配，旧的成功结果不会充当本次结果。

工作流使用 GitHub 自动提供的短期 `GITHUB_TOKEN`，仅为拟合 job 声明 `contents: write`，将状态及已有结果提交到 `pkpd-results` 分支。
无需个人 token。结果分支不触发计算，发布操作也不更改源码分支。标准 Actions artifacts 仍保留。
该分支位于公开仓库，因此结果和模型日志也公开；发布范围只包含 `joint_results`、`ci-logs` 中的结果和日志，不转存进程环境或凭据。
结果分支保留历史运行，适合当前小体积结果；超过 90 MiB 的压缩包会明确报错并保留标准 artifact 路径。

本地传输测试已验证临时 Git 仓库中的运行状态更新、精确 SHA 匹配、压缩包校验和提取、失败退出码，使用的都是合成文件。
实际 MATLAB 是否成功仍以具体云端运行的 `ci_status.json` 和日志为准。

## API 连接限制与最小权限

当前 Codex 云环境的代理在 HTTPS CONNECT 阶段返回 `403 Forbidden`，正文为 `Domain forbidden`；目标域名为 `api.github.com`。
这发生在 GitHub 身份验证之前，因此不能把 `gh auth status` 的失败判断为 token 无效，也不能据此认定组织禁止 Actions。
Git HTTPS 读写已验证可用；公开 GitHub 网页可见运行摘要，但未登录网页要求登录才能查看日志。
Git 回传路径可绕开对外部 Actions API 的依赖，前提是 GitHub runner 能启动且 `GITHUB_TOKEN` 能写入结果分支。

若要恢复原生 API 路径：

1. 在 Codex 当前环境设置 → Internet access / 网络访问 → 自定义允许域名添加 `api.github.com`，保留已有规则。日志/artifact 下载可能重定向至其他域名，应按实际 Location 主机追加，勿把带签名的完整 URL 写入日志。
2. 网络恢复后先用现有 GitHub 绑定测试；Git push 权限不等于 Actions REST 权限。
3. 若 API 明确返回授权不足，GitHub App 仅需该仓库的 Actions 读写（触发、状态、日志、artifacts）；Metadata 只读。已有 Git 代码/工作流写入权限可继续使用。
4. 若平台连接不能授予 Actions 权限，可使用仅选中 `RainaYuWang/dox-pkpd` 的 fine-grained PAT：Actions read/write，Metadata read；不需要额外 Contents write（代码已能通过现有 Git 通道推送）。只需读取结果时 Actions read 即可。仅在云环境安全凭据设置中输入，禁止发到聊天或写入仓库。

Git 回传方案本身不要求上述 API 域名或个人 token。若组织禁止 Actions、禁止官方 Actions 或限制结果分支写入，需要在仓库/组织 Actions 策略或规则集中针对本仓库及 `pkpd-results` 做最小调整，不需要放开所有仓库权限。

## 科学解释与当前证据

quick 用于验证运行流程，不是最终拟合。检查训练数据上的 SSE、优化器退出状态和数值轨迹，不能据此宣称预测能力已获验证或个体 PK 参数唯一可辨识。
cell number 来自论文图像数字化；旧 CF+CB 是模拟参考，不是独立实验数据或验证集。未引入 damage-memory 或阈值机制。

创建工作流或通过静态语法检查，不等于 MATLAB 已运行。本说明不预先宣称 Actions 运行成功；请以具体 run ID、日志、`ci_status.json` 和当次 artifacts 为证据。

实际运行证据及本次 quick 诊断见 [2026-11-11 云端验证记录](VALIDATION_2026-11-11.md)。

固定 PK 敏感性分析使用 `sensitivity` 模式，或推送到 `pkpd-sensitivity`；它调用 `run_joint_pkpd(8,true)`，追加 9 个情形、每个 8 起点的拟合，并输出有效杀伤速率。详见 [敏感性定义](../reports/model_d_sensitivity_definition.md)。

已完成的 MATLAB 敏感性运行、参数/曲线比较及自动取回后重建图表的命令见 [敏感性结果](../reports/model_d_sensitivity_results.md)。分析脚本为 `ci/analyze_model_d_sensitivity.py --results <joint_results路径> --output reports`，依赖 numpy>=2、scipy、matplotlib。

共享 beta 假设的诊断使用 `dose-beta` 模式，或将相关代码推送至 `pkpd-dose-beta`。调用 `run_joint_pkpd(8,false,true)`，保持归档 D 的参考 PK、生长参数和最佳 kBF，追加三个剂量各自的 beta 拟合及一维扫描复核。输出、失败保存及 Git 自动取回沿用相同流程。详见 [分剂量 beta 诊断定义](../reports/model_d_dose_beta_definition.md)。

已完成的实际诊断、各剂量 beta/SSE、残差形状、625 nM 早期下降及数值审计见 [诊断结果](../reports/model_d_dose_beta_results.md)。图表可由 `ci/analyze_model_d_dose_beta.py --results <joint_results路径> --output reports` 重建。
