# 最终模型的结合药物输入 C_B(t)

[下载浓度 CSV](CB.csv) · [查看九宫格](CB.png) · [下载四个 k 和两个 V](pk_parameters.csv)

`CB_nM` 是供 damage-memory 模型读取的胞内结合药物浓度，单位 nM。来源是 184 点 Fig. 3 full 拟合选定的 `dynamic_hill`（模型 9）。九条曲线用同一套 PK 参数从给药开始计算，各组只改变胞外暴露。细胞数、分组 θ 和 Hill 杀伤参数均不反馈到 PK。[对应的 full 拟合参数](../full/parameters.csv) · [模型比较](../README_zh.md)

## 选择面板

| 给药浓度（nM） | 暴露 6 h | 暴露 12 h | 暴露 24 h |
| ---: | :---: | :---: | :---: |
| 19.53125 | a | b | c |
| 156.25 | d | e | f |
| 1250 | g | h | i |

同一浓度的条件在较早的洗脱时刻前有相同的 CB。初始 CF(0)=CB(0)=0。到给定的 washout 时刻，胞外 CE 从给药浓度切换为零，胞内 CF、CB 连续，时钟继续从给药开始计时。

## 读取 CSV

每行对应一个面板和一个时间。九组均覆盖 0–750 h，每个整小时都有记录，另保留全部实际细胞观测时刻，重合时刻只出现一次。表按面板、时间升序排列，共 6,884 行数据。

| 列 | 含义 |
| --- | --- |
| `panel` | 论文 Fig. 3 面板 a–i |
| `dose_nM`, `wash_h` | 该组给药浓度及从给药到洗脱的时间 |
| `time_h` | 自给药起算的时间 |
| `CE_nM` | 矩形胞外输入。0 ≤ t < wash_h 时等于 dose_nM，从 wash_h 起为 0 |
| `CB_nM` | 该时刻的模型结合药物浓度 |
| `is_cell_observation_time` | 1 表示同面板和时刻有细胞数观测，可与 [184 点数据](../../../fig3/observations.csv) 配对 |
| `after_last_cell_observation` | 1 表示已超过该面板末次细胞观测，该段为模型外推 |

所有行的 CB 均为模拟值。图中实线与虚线分别标记细胞观测窗口内和窗口后的预测，例如 i 在 509 h 后用虚线。九个面板共用刻度，纵轴在 0–1 nM 为线性、1 nM 以上为对数，使低浓度和高浓度曲线可在同一张图上读取。

## 方程与参数

```text
dCF/dt = kEF · (VE/VI) · CE(t) − (kFE + kFB) · CF + kBF · CB
dCB/dt = kFB · CF − kBF · CB
```

| 参数 | 数值 | 单位 | 用途 |
| --- | ---: | --- | --- |
| kEF | 3.5 × 10⁻⁶ | h⁻¹ | 胞外至胞内输入系数，固定 |
| kFE | 1.25 | h⁻¹ | 胞内游离药物外流，固定 |
| kFB | 0.21 | h⁻¹ | 胞内游离至结合，固定 |
| kBF | 0.00715547505854946 | h⁻¹ | 结合至游离，九组共享拟合 |
| VI | 0.004629 | µL | 胞内有效体积，固定 |
| VE | 250 | µL | 胞外体积，固定 |

五个固定量沿用 [参考 PK 校准](../../model_d_definition.md)；kBF 取自所选模型的 full 估计。固定平台对照在推导渐近值时使用 kBF=0。本交接包包含正 kBF 下的建立、解离回流与再结合。

## 在自己的 ODE 中取值

从仓库根目录运行以下 MATLAB 示例。它给出 g 面板在任意非负时间（小时）的 CB。将 `dose` 和 `wash` 换成上表中的其他组合即可使用其他面板。

```matlab
addpath('fig3')
pk = readtable('reports/fig3_184/cb_input/pk_parameters.csv');
kBF = pk.value(strcmp(pk.parameter, 'kBF'));
dose = 1250;
wash = 6;
t = [0; 0.25; 6; 6.1; 24; 267; 750];
[~, CB_nM] = fig3_pk(t, dose, wash, kBF);
```

在 ODE 右端函数中，同一调用也接受求解器传入的标量时间。[拟合时使用的 MATLAB 分段解析求解函数](../../../fig3/fig3_pk.m)

从仓库根目录运行 `python fig3/export_cb.py` 可复现本目录的 CSV 和 PNG，需要 NumPy、SciPy 和 Matplotlib。导出时按 washout 分段传播线性 PK 方程，并与归档的全部 3,784 个模型 9 面板–时间点比较。最大 CB 绝对差为 2.81 × 10⁻¹² nM。另用分段 ODE 独立核对全部 6,884 个导出值，最大绝对差为 1.51 × 10⁻¹¹ nM。[导出脚本](../../../fig3/export_cb.py) · [归档 MATLAB 时程](../full/curves.csv)
