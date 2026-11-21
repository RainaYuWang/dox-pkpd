# SUM-149PT 多柔比星 PK–PD

本项目比较了 McKenna 等论文 Fig. 3 中 SUM-149PT 细胞在九种给药条件下的生长和恢复。184 个细胞数观测用于联合拟合九种模型。最终选用**动态可逆结合 PK + 自由 Hill 杀伤**，共 14 个自由参数，按组内最大观测细胞数缩放的总 SSE 为 **0.629771**。

## 下一组从这里开始

做 damage memory 时，可直接使用下表中的 `CB_nM` 作为随时间变化的输入。九组共享四个速率常数和两个体积，浓度曲线的区别来自胞外给药浓度及洗脱时刻。

| 需要什么 | 文件 |
| --- | --- |
| 九组完整的结合药物浓度，从给药起 0–750 h | [下载 CB.csv](reports/fig3_184/cb_input/CB.csv) · [九宫格图](reports/fig3_184/cb_input/CB.png) |
| 四个 k、两个 V、单位及固定或拟合状态 | [PK 参数表](reports/fig3_184/cb_input/pk_parameters.csv) |
| 面板和条件映射、洗脱约定、CSV 字段、MATLAB 任意时刻取值 | [CB 交接说明](reports/fig3_184/cb_input/README.md) |
| 与模型输入配对的细胞数 | [184 点观测](fig3/observations.csv) · [原图取点核对](fig3/point_audit.png) |

## 本轮结果

在同样使用自由 Hill 的情况下，固定平台、保留早期过渡且固定零解离、允许自由解离的训练 SSE 分别为 **3.646524、3.930409、0.629771**。自由解离模型抓住了 g、h 的后期回升，g 的早期下降仍偏浅。结合药物浓度是模型推算值，kBF 由细胞计数间接估计。i 面板有 22 个细胞数观测，最后一个在 509 h。

- [中文汇总 PDF](reports/fig3_184/comparison/fig3_full_comparison.pdf) · [九模型比较、分组 θ、残差与 kBF profile](reports/fig3_184/README_zh.md)
- [同一杀伤函数下比较 CB 输入](reports/fig3_184/comparison/fit_cb.png) · [同一动态 PK 形式下比较杀伤函数](reports/fig3_184/comparison/fit_kill.png)
- [数据、方程和复现说明](fig3/README_zh.md) · [原始 full 输出](reports/fig3_184/full) · [MATLAB 运行记录](https://github.com/RainaYuWang/dox-pkpd/actions/runs/44603901802)
- [论文正文与补充材料](references/papers) · [历史分析索引](README_zh.md)
