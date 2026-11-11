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

quick 共 9 个起点（control/A/B/C 各 2 个，加 1 个接近不可逆的可逆起点）；full 共 25 个（control/A 各 4 个、B/C 各 8 个，加 1 个额外起点）。
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

## 首次成功后再启用代码触发

当前工作流仅支持手动触发，不会因文档提交自动计算。
首次实际 MATLAB quick 成功并检查 artifacts 后，在 YAML 的 `on:` 下增加：

```yaml
  push:
    branches: [main]
    paths:
      - '**/*.m'
      - '.github/workflows/matlab-pkpd.yml'
```

自动触发时没有 `inputs.mode`，现有工作流会默认执行 quick。此配置应在首次成功后提交，不要提前启用。

## 科学解释与当前证据

quick 用于验证运行流程，不是最终拟合。检查训练数据上的 SSE、优化器退出状态和数值轨迹，不能据此宣称预测能力已获验证或个体 PK 参数唯一可辨识。
cell number 来自论文图像数字化；旧 CF+CB 是模拟参考，不是独立实验数据或验证集。未引入 damage-memory 或阈值机制。

创建工作流或通过静态语法检查，不等于 MATLAB 已运行。本说明不预先宣称 Actions 运行成功；请以具体 run ID、日志、`ci_status.json` 和当次 artifacts 为证据。
