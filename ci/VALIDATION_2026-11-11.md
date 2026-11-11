# 实际云端验证记录（2026-11-11）

代码提交：`41481687639985e2953803a3f8b96bf40426acf9`。
[GitHub 运行 43099035215](https://github.com/RainaYuWang/dox-pkpd/actions/runs/43099035215)，由 `pkpd-cloud` 分支的代码 push 自动触发。
本次实际完成：修改/推送 → GitHub-hosted Linux MATLAB → 自动发布结果分支 → Codex Git 获取 → 校验、解压及分析。

## 运行证据

- MATLAB `24.2.0.3361896 (R2024b) Update 10`，平台 `GLNXA64`。
- 命令：`addpath('ci'); run_pkpd_ci('quick')`，内部 `run_joint_pkpd(2)`。
- MATLAB 和输出检查耗时 **22.926872 秒**；其中模型及图片输出约 21.6 秒。安装和排队时间不计入。
- `ci_status.json`: `success`，9 个起点全部完成，9 个 checkpoint 的退出标记均为 1。
- 四张 PNG、CSV、最终 MAT、9 个中间 MAT 及日志均已自动取回；已用 Python SciPy 读取全部 checkpoint，并核对图片和 CSV。
- Git 回传包 SHA-256：`ad94ab2ed3d10f184709eb9763d4c74ce7bf7b27dd9600767917d60af75a8b94`。
- [结果分支记录](https://github.com/RainaYuWang/dox-pkpd/tree/pkpd-results/runs/41481687639985e2953803a3f8b96bf40426acf9/43099035215-1)包含 manifest 和 ZIP。
- 标准 Actions artifact 也已生成：`pkpd-quick-43099035215-1`，ID `13890003147`，网页列出约 663 KB。它与 Git 回传包使用不同 ZIP 包装，摘要不应混用。

## 本次 quick 诊断

| 模型 | scaled SSE |
| --- | ---: |
| Old PK / direct PD | 2.5914442953 |
| Joint irreversible | 2.5846219890 |
| Joint reversible | 0.9939678665 |

日志没有 ODE 求解错误。各起点实际耗时约 0.01–3.77 秒，没有重现本机 MATLAB 未响应。
不可逆代表参数 `kFE=9.9994307` 接近人为上界 10，后续正式分析应检查边界敏感性；本轮未改边界。
可逆模型两个普通起点得到接近的 SSE，但不同的个体 PK 参数；额外接近不可逆的起点 SSE=2.5845421，虽然 exitflag=1，仍停在较差解。
等价参数的最大相对细胞数差异为 `3.1550e-6`，低于代码现有 `1e-3` 警告阈值。
图中仍有系统偏差：156.25 nM 后期细胞数被高估，625 nM 的早期下降过深。

quick 仅用于运行流程检查，未运行 full。收敛不代表全局最优或参数唯一；训练 SSE 不验证预测能力。
数据仍是论文图像数字化点，旧 CF+CB 仍仅为模拟参考。未改变科学假设、模型机制、数据、边界或求解精度。

## Actions API 诊断与替代通路

对 `api.github.com:443` 的 HTTPS CONNECT 返回：

```text
HTTP 403 Forbidden
content-type: text/plain
content-length: 16
server: envoy
connection: close

Domain forbidden
```

未看到 GitHub 的 request-id、rate-limit 或 accepted-permissions 响应头。当前云环境限制域名列表不含 `api.github.com`。
POST 工作流触发、GET 运行状态、GET 日志、GET artifact 列表，以及已知 artifact ID 的 ZIP 下载都在网络层受阻；无法据此判断 GitHub API 凭据是否有效或 Actions 权限是否充分。
公开网页可读，但未登录时要求登录查看原始 job 日志；网页摘要还曾滞留在 Queued，因此实际 MATLAB 完成证据使用带准确 commit/run ID 的运行日志及结果 manifest。

已实际证明 Git 读取/推送、push 自动触发、Actions 自带 GITHUB_TOKEN 写入结果分支、Codex 获取完整结果均可用。
完整自动化可继续使用 Git 回传，无需新增 token 或修改用户设置。
若另需原生 REST 下载，需先在云环境网络允许域名中加入 `api.github.com`，再验证现有绑定；当前配置草稿接口返回 INVALID_ARGUMENT，不能由本次工具确认修改。
原生 API 所需的最小仓库权限及安全配置方式见 [运行说明](README_zh.md#api-连接限制与最小权限)。
