"""Render the Fig. 3 comparison from archived CSVs. Requires matplotlib and a CJK font."""
import csv
import os
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator, FuncFormatter

ROOT = Path(__file__).resolve().parents[1]
FULL = ROOT / 'reports/fig3_184/full'
OUT = ROOT / 'reports/fig3_184/comparison'
LETTERS = 'abcdefghi'
MODEL_NAMES = ['平台·线性', '平台·n=1', '平台·自由 n', '零解离·线性', '零解离·n=1',
               '零解离·自由 n', '自由解离·线性', '自由解离·n=1', '自由解离·自由 n']
COLORS = ['#12658a', '#c36a10', '#98426c']
DARK = '#183448'
PALE = '#f1f4f6'
GAP = '#e4e9ed'
PARAMS = ['kp', 'beta', 'Kmax', 'Khalf', 'n', 'kBF']


def read(file):
    with file.open(newline='') as stream:
        return list(csv.DictReader(stream))


def column(rows, key):
    return np.array([float(r[key]) for r in rows])


def group(rows):
    result = {}
    for row in rows:
        key = int(row['model_id']), int(row['panel_id'])
        result.setdefault(key, []).append(row)
    return result


def load():
    summary = read(FULL / 'summary.csv')
    params = {r['model']: {} for r in summary}
    for r in read(FULL / 'parameters.csv'):
        params[r['model']][r['parameter']] = float(r['estimate'])
    bounds = {r['model']: [] for r in summary}
    for r in read(FULL / 'bounds.csv'):
        if r['near_bound'] == '1':
            value = params[r['model']][r['parameter']]
            side = '上' if abs(np.log(float(r['upper']) / value)) < 0.01 else '下'
            bounds[r['model']].append(r['parameter'].replace('theta_', 'θ') + side)
    curves = group(read(FULL / 'curves.csv'))
    residuals = group(read(FULL / 'residuals.csv'))
    observations = {j: [] for j in range(1, 10)}
    for r in read(ROOT / 'fig3/observations.csv'):
        observations[LETTERS.index(r['panel']) + 1].append(r)
    for r in summary:
        sse = 0
        for j, panel in enumerate(LETTERS, 1):
            rr = residuals[int(r['model_id']), j]
            obs = observations[j]
            np.testing.assert_allclose(column(rr, 'time_h'), column(obs, 'time_h'), atol=1e-9)
            y = column(obs, 'cells')
            np.testing.assert_allclose(column(rr, 'observed_N'), y)
            error = (column(rr, 'predicted_N') - y) / y.max()
            np.testing.assert_allclose(error, column(rr, 'scaled_residual'), atol=1e-12)
            part = error @ error
            assert abs(part - float(r['SSE_' + panel])) < 1e-10
            sse += part
        assert abs(sse - float(r['SSE'])) < 1e-10
    return summary, params, bounds, curves, residuals, observations


def tables(summary, params, bounds):
    main, theta, errors = [], [], []
    for index, r in enumerate(summary):
        p = params[r['model']]
        shared = [f'{p[k]:.5g}' if k in p else '—' for k in PARAMS]
        if index % 3 == 1:
            shared[4] = '1 fixed'
        if 3 <= index < 6:
            shared[5] = '0 fixed'
        main.append([MODEL_NAMES[index], r['free_parameters'], f"{float(r['SSE']):.6f}",
                     *shared, '、'.join(bounds[r['model']]) or '无'])
        theta.append([MODEL_NAMES[index], *[f'{p["theta_"+l]:,.0f}' for l in LETTERS]])
        errors.append([MODEL_NAMES[index], *[f'{float(r["SSE_"+l]):.4g}' for l in LETTERS]])
    return main, theta, errors


def sheet(title, subtitle):
    fig = plt.figure(figsize=(14, 10), facecolor='white')
    fig.text(.06, .95, title, fontsize=24, weight='bold', color=DARK)
    fig.text(.06, .905, subtitle, fontsize=12, color=DARK)
    return fig


def table(fig, rect, headers, values, widths=None):
    ax = fig.add_axes(rect)
    ax.axis('off')
    tab = ax.table(cellText=values, colLabels=headers, colWidths=widths,
                   cellLoc='center', loc='center', bbox=[0, 0, 1, 1])
    tab.auto_set_font_size(False)
    tab.set_fontsize(10)
    for (r, c), cell in tab.get_celld().items():
        cell.set_edgecolor('white')
        cell.set_facecolor(DARK if r == 0 else (PALE if r % 2 else 'white'))
        cell.get_text().set_color('white' if r == 0 else DARK)


def figure_grid(title, subtitle, ids, labels, metric, data, observations):
    fig = sheet(title, subtitle)
    axes = fig.subplots(3, 3)
    fig.subplots_adjust(left=.075, right=.977, bottom=.10, top=.79, wspace=.24, hspace=.54)
    if metric == 'scaled_residual':
        legend = [Line2D([], [], color=c, marker=m, ls='none', label=l)
                  for c, m, l in zip(COLORS, ('o','+','x'), labels)]
    else:
        legend = [Line2D([], [], color=c, lw=2.4, label=l) for c, l in zip(COLORS, labels)]
    if metric == 'N':
        legend.append(Line2D([], [], marker='o', linestyle='none', color='white', markerfacecolor='black',
                             markeredgecolor='white', markersize=7, label='拟合使用的观测'))
    legend += [Line2D([], [], color='#566775', ls=':', label='洗脱'),
               Patch(facecolor=GAP, label='缺少可辨观测 · 模型外推区')]
    fig.legend(handles=legend, loc='upper left', bbox_to_anchor=(.055,.887),
               ncol=len(legend), frameon=False, fontsize=10)
    for j, ax in enumerate(axes.flat, 1):
        obs = observations[j]
        x = column(obs, 'time_h')
        wash = int(obs[0]['wash_h'])
        limit = 720 if j == 9 else x[-1]*1.012
        if j == 9:
            left, right = x[-1], limit
        else:
            left = right = None
        if left is not None:
            ax.axvspan(left, right, color=GAP, zorder=0)
        ax.axvline(wash, color='#526573', lw=1.2, ls=':', zorder=1)
        if metric == 'scaled_residual':
            ax.axhline(0, lw=.8, color='#84939e')
        for model, color, marker in zip(ids, COLORS, ('o','+','x')):
            series = data[model, j]
            xx, yy = column(series, 'time_h'), column(series, metric)
            if metric == 'scaled_residual':
                ax.scatter(xx, yy, color=color, s=14, marker=marker, zorder=3)
            elif left is None:
                ax.plot(xx, yy, color=color, lw=1.8)
            else:
                for mask, style in [(xx <= left,'-'), ((xx >= left) & (xx <= right),'--'),
                                    (xx >= right,'-')]:
                    ax.plot(xx[mask], yy[mask], color=color, ls=style, lw=1.8)
        if metric == 'N':
            ax.scatter(x, column(obs, 'cells'), s=18, facecolor='#181b20', edgecolor='white',
                       linewidth=.55, zorder=5)
            peak = max(np.nanmax(column(data[m,j], 'N')) for m in (3,6,7,8,9))
            ax.set_ylim(0, 1.13*max(peak, max(column(obs,'cells'))))
        elif metric == 'scaled_residual':
            scale = max(max(abs(column(data[m,j], metric))) for m in range(1,10))
            ax.set_ylim(-1.2*scale, 1.2*scale)
        else:
            ax.set_ylim(bottom=0)
        # Label only the actual unrecoverable spans. Ordinary daily samples remain visible.
        if left is not None:
            ax.text((left+right)/2, .94, '模型外推', transform=ax.get_xaxis_transform(),
                    ha='center', va='top', fontsize=9, color='#45596a')
        ax.set_title(f"{LETTERS[j-1]}   {obs[0]['dose_nM']} nM · {wash} h · {len(x)} 点",
                     fontsize=12, loc='left', color=DARK, weight='bold', pad=6)
        ax.set_xlim(-.015*limit, limit)
        ax.xaxis.set_major_locator(MaxNLocator(5))
        ax.yaxis.set_major_locator(MaxNLocator(5))
        if metric == 'N':
            ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v/1000:g}'))
        ax.grid(axis='y', color='#d4dce2', lw=.5, zorder=0)
        ax.spines[['top','right']].set_visible(False)
        ax.spines[['bottom','left']].set_color('#93a2ac')
        ax.tick_params(labelsize=9, length=3)
    ylabels = {'N':'细胞数（千个）','scaled_residual':'归一化残差（预测 − 观测）',
               'CB_nM':'结合药物 CB（nM）','kill_h_inv':'有效杀伤率（每小时）'}
    fig.supylabel(ylabels[metric], x=.012, fontsize=12)
    fig.supxlabel('给药后时间（h）', y=.056, fontsize=12)
    fig.text(.06,.022,'g、h、i 含手动标记的近似读数。i 的观测止于 509 h，灰色为外推区，曲线只画到最后一个输入点。',fontsize=10,color=DARK)
    return fig


def comparison_notes(summary, params, residuals):
    scores = column(summary, 'SSE')
    def improvement(old, new):
        return 100 * (1 - scores[new-1]/scores[old-1])

    changes = {key: float(summary[2]['SSE_'+key]) - float(summary[8]['SSE_'+key])
               for key in LETTERS}
    ordered = sorted(changes, key=changes.get, reverse=True)[:3]
    share = 100 * sum(changes[key] for key in ordered) / (scores[2]-scores[8])
    largest = sorted(LETTERS, key=lambda key: float(summary[8]['SSE_'+key]), reverse=True)[:3]
    g = residuals[9, 7][1]
    turns = next(r for r in read(FULL/'regrowth.csv') if r['model_id']=='9' and r['panel_id']=='7')
    fit = params['dynamic_hill']
    k = fit['kBF']
    t = float(g['time_h'])
    n0 = float(residuals[9, 7][0]['observed_N'])
    net = fit['kp'] - fit['Kmax']
    integral = t if net == 0 else np.expm1(net*t)/net
    lower_g = np.exp(net*t)/(1/n0 + fit['kp']*integral/fit['theta_g'])
    profile = read(FULL/'dynamic_hill_kbf_profile.csv')
    before = max((r for r in profile if float(r['kBF_fixed']) < k*.9999), key=lambda r: float(r['kBF_fixed']))
    after = min((r for r in profile if float(r['kBF_fixed']) > k*1.0001), key=lambda r: float(r['kBF_fixed']))
    fall = [-float(r['last_obs_change_from_peak_pct']) for r in read(FULL/'platform.csv')
            if np.isclose(float(r['kBF']),k,rtol=1e-8,atol=0)]
    return [
        ('分组 θ 和固定平台',f'各组 θ 可调节增长水平和晚期细胞数。固定平台自由 Hill 的 g、h SSE 分别为 {float(summary[2]["SSE_g"]):.3f}、{float(summary[2]["SSE_h"]):.3f}。\n组内杀伤率与 θ 均恒定时，细胞数只能单调变化或保持水平，无法表达 g、h 的先降后升。'),
        ('饱和与解离',f'固定平台下，线性改为 n=1 降低 SSE {improvement(1,2):.1f}%，再放开 n 降低 {improvement(2,3):.2f}%。自由解离下分别降低 {improvement(7,8):.1f}%、{improvement(8,9):.1f}%。\n自由 Hill 的参数数依次为 13、13、14。单独保留早期过渡且固定 kBF=0 时，SSE 为 {scores[5]:.4f}。'),
        ('误差集中在哪些面板',f'固定平台自由 Hill 到动态自由 Hill 的总 SSE 净降幅中，{"、".join(ordered)} 合计贡献 {share:.1f}%。\n动态自由 Hill 的剩余误差以 {"、".join(largest)} 面板最大。逐面板数值见第 2 页。'),
        ('g 的早期高估',f'给药后 {float(g["time_h"]):.1f} h 的观测约 {float(g["observed_N"]):,.0f} 个细胞，动态自由 Hill 预测 {float(g["predicted_N"]):,.0f} 个。\n即使从 0 h 就达到最大杀伤率 Kmax，该时刻也至少有 {lower_g:,.0f} 个。拟合曲线到 {float(turns["predicted_min_t_h"]):.0f} h 才降至最低约 {float(turns["predicted_min_N"]):,.0f} 个。'),
        ('实际窗口与条件 profile',f'最佳 kBF 为 {k:.5g} h⁻¹。九组末次观测时模拟 CB 较各自峰值下降 {min(fall):.0f}%–{max(fall):.0f}%。\n取相邻的 kBF 网格点 {float(before["kBF_fixed"]):.4g} 和 {float(after["kBF_fixed"]):.4g} 并重优化其余参数后，SSE 分别为 {float(before["SSE"]):.4f} 和 {float(after["SSE"]):.4f}，后者 n={float(after["n"]):.2f}。'),
        ('数据能支持的结论','i 的 22 个点约束到 509 h 的变化，后续再增长、消失或长期平台仍缺少数据。手动读数尤其在近零区有误差。\n细胞计数支持这组模型中随时间变化的有效杀伤。真实 CB 浓度和解离速率仍需独立药物测量来验证。')]


def main():
    OUT.mkdir(exist_ok=True, parents=True)
    family = 'Noto Sans CJK JP'
    for path in filter(None, os.environ.get('FIG3_CJK_FONTS', '').split(os.pathsep)):
        font_manager.fontManager.addfont(path)
        family = font_manager.FontProperties(fname=path).get_name()
    plt.rcParams.update({'font.family':[family, 'DejaVu Sans'], 'pdf.fonttype':42,
                         'axes.unicode_minus':False, 'axes.labelcolor':DARK})
    summary, params, bounds, curves, residuals, obs = load()
    main_rows, theta, errors = tables(summary, params, bounds)
    with PdfPages(OUT/'fig3_full_comparison.pdf') as pdf:
        fig = sheet('Fig. 3 九条件联合拟合', 'SUM-149PT  ·  184 个观测  ·  9 个模型  ·  MATLAB full 输入 a04281f')
        fig.text(.06,.84,f"自由 Hill 的总 SSE 在固定平台下为 {float(summary[2]['SSE']):.4f}，允许解离后为 {float(summary[8]['SSE']):.4f}。",fontsize=17,color=DARK)
        table(fig,[.06,.40,.88,.39],['模型','参数数','总 SSE','kp','β','Kmax','Khalf','n','kBF','触界'], main_rows,
              [.14,.07,.10,.10,.11,.10,.09,.08,.10,.11])
        fig.text(.06,.31,'kp、Kmax、kBF 单位为 h⁻¹，β 为 nM⁻¹ h⁻¹，Khalf 为 nM。fixed 表示固定，— 表示不适用。')
        fig.text(.06,.26,'平台值由固定参考 PK 和暴露时长解析给出。固定 kEF=3.5×10⁻⁶、kFE=1.25、kFB=0.21 h⁻¹，VE=250、VI=0.004629 µL。')
        fig.text(.06,.21,'仅列最优解的触界情况。近界定义沿用归档的对数距离 ≤0.01。表内所列 θ 上界均为 10⁸ cells。')
        fig.text(.06,.15,'各组以首观测初始化 N，并以本组最大观测细胞数缩放残差。后续表和图均使用同一批数据。')
        fig.text(.06,.10,'固定平台直接代入常数 CB，所以该模型没有 kBF 参数。零解离模型则从给药起计算 CB(t)，取 kBF=0。')
        starts = column(read(FULL/'dynamic_hill_starts.csv'), 'SSE')
        fig.text(.06,.05,f'动态自由 Hill 的 {len(starts)} 个联合起点中，{sum(starts < min(starts)+1e-5)} 次达到最低 SSE（容差 10⁻⁵），其余起点最高 SSE 为 {max(starts):.4f}。')
        pdf.savefig(fig); plt.close(fig)
        fig = sheet('逐条件参数和误差', '面板 a–i 按 Fig. 3 行优先排列。各行共用上一页的共享参数。')
        fig.text(.06,.85,'各组 θ（cells，按整数展示）',fontsize=15,color=DARK)
        table(fig,[.06,.51,.90,.31],['模型',*LETTERS],theta,[.13]+[.097]*9)
        fig.text(.06,.45,'逐面板归一化 SSE',fontsize=15,color=DARK)
        table(fig,[.06,.10,.90,.32],['模型',*LETTERS],errors,[.13]+[.097]*9)
        fig.text(.06,.04,'i 有 22 个点（含固定初值），其中 20 个是手动近似读数。509 h 以后没有用于拟合的观测。')
        pdf.savefig(fig); plt.close(fig)
        sets = [
            ('fit_cb','三种 CB 输入下的拟合曲线','三条曲线均使用自由 Hill；共享参数及分组 θ 已各自联合估计。',
             [3,6,9],['固定平台','动态 · kBF=0','动态 · 自由 kBF'],'N'),
            ('fit_kill','三种杀伤函数下的拟合曲线','三条曲线均使用动态 CB；每个模型独立估计共享 kBF。',
             [7,8,9],['线性','Hill · n=1','Hill · 自由 n'],'N'),
            ('cb','结合药物的时间过程','自由 Hill 三模型的归档 CB 输入。固定平台同时给出不可逆渐近参考值。',
             [3,6,9],['固定平台','动态 · kBF=0','动态 · 自由 kBF'],'CB_nM'),
            ('kill','有效杀伤率的时间过程','三种 CB 输入下的自由 Hill 杀伤率。纵轴反映输入变化与各模型拟合参数的共同影响。',
             [3,6,9],['固定平台','动态 · kBF=0','动态 · 自由 kBF'],'kill_h_inv')]
        for k, label in enumerate(('固定平台','动态 kBF=0','动态自由 kBF')):
            sets.append((f'residual_{k+1}',f'{label}下的时间残差',
                         '仅在实际取点时刻绘制残差；正值表示模型高估。每个面板在三页残差图使用同一纵轴。',
                         list(range(3*k+1,3*k+4)),['线性','Hill · n=1','Hill · 自由 n'],'scaled_residual'))
        for name,title,subtitle,ids,labels,metric in sets:
            fig=figure_grid(title,subtitle,ids,labels,metric,
                            residuals if metric=='scaled_residual' else curves,obs)
            fig.savefig(OUT/f'{name}.png',dpi=145)
            pdf.savefig(fig); plt.close(fig)

        fig=sheet('比较结果与证据边界','最低训练误差来自动态自由 Hill。结论限于九组可辨点和已比较的模型。')
        notes = comparison_notes(summary, params, residuals)
        for i,(heading,body) in enumerate(notes):
            y=.835-i*.135
            fig.text(.06,y,heading,fontsize=14,color=DARK,weight='bold')
            fig.text(.06,y-.028,body,fontsize=11,va='top',linespacing=1.55,color=DARK)
        pdf.savefig(fig); plt.close(fig)


if __name__ == '__main__':
    main()
