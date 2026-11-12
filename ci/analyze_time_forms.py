#!/usr/bin/env python3
"""Analyze successful, source-verified A/B runs alongside stored old predictions."""
import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--output", default=Path("reports/time_forms"), type=Path)
    parser.add_argument("--quick-results", type=Path, help="Matching successful quick run for comparison")
    args = parser.parse_args()
    source, output = args.results, args.output
    status = json.loads((source / "ci_status.json").read_text())
    paper = json.loads((source / "paper_time_definitions.json").read_text())
    assert status["status"] == "success" and status["mode"] in ("time-quick", "time-full")
    assert paper["verified"] is True and status["old_models_refitted"] is False

    def read(name, directory=source):
        with (directory / (name + ".csv")).open(newline="") as handle:
            return [{key: value if key == "model" else float(value or "nan")
                     for key, value in row.items()} for row in csv.DictReader(handle)]

    def select(rows, model, dose):
        return [row for row in rows if row["model_id"] == model and row["dose_nM"] == dose]

    def write(name, rows):
        with (output / (name + ".csv")).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    def clean(value):
        if isinstance(value, float) and not math.isfinite(value):
            return None
        if isinstance(value, list):
            return [clean(item) for item in value]
        if isinstance(value, dict):
            return {key: clean(item) for key, item in value.items()}
        return value

    summary, parameters, curves, obs, starts = [read("time_form_" + name) for name in
                                               ("summary", "parameters", "predictions", "residuals", "starts")]
    assert len(obs) == 6 * 32
    for row in summary:
        actual = sum(point["scaled_residual"] ** 2 for point in obs if point["model_id"] == row["model_id"])
        assert abs(actual - row["scaled_SSE"]) < 1e-10
    output.mkdir(parents=True, exist_ok=True)
    doses = [156.25, 312.5, 625]
    metrics, start_metrics, effect_metrics, extrema = [], [], [], []
    for dose in doses:
        for model in range(1, 7):
            predicted, points = select(curves, model, dose), select(obs, model, dose)
            late = [p["residual"] for p in points if p["time_h"] >= 300]
            early_residual = [p["residual"] for p in points if 0 < p["time_h"] <= 200]
            early = min((p for p in predicted if p["time_h"] <= 200), key=lambda p: p["N"])
            for index in range(1, len(predicted) - 1):
                left, point, right = predicted[index - 1:index + 2]
                kind = ("minimum" if point["N"] < left["N"] and point["N"] < right["N"] else
                        "maximum" if point["N"] > left["N"] and point["N"] > right["N"] else None)
                if kind and point["time_h"] <= 200:
                    extrema.append(dict(model_id=model, dose_nM=dose, kind=kind,
                                        data_time_h=point["time_h"], N=point["N"]))
            metrics.append(dict(model_id=model, dose_nM=dose,
                                scaled_SSE=sum(p["scaled_residual"] ** 2 for p in points),
                                RMSE_cells=float(np.sqrt(np.mean([p["residual"] ** 2 for p in points]))),
                                early_min_N=early["N"], early_min_data_time_h=early["time_h"],
                                early_count=len(early_residual),
                                early_mean_residual=float(np.mean(early_residual)),
                                early_RMSE=float(np.sqrt(np.mean(np.array(early_residual) ** 2))),
                                late_count=len(late), late_mean_residual=float(np.mean(late)),
                                late_RMSE=float(np.sqrt(np.mean(np.array(late) ** 2))),
                                final_residual=points[-1]["residual"]))
        for model in (3, 4):
            rows = select(starts, model, dose)
            p = select(parameters, model, dose)[0]
            near = [row for row in rows if row["scaled_SSE"] <= p["scaled_SSE"] + 1e-6]
            observed_peak = max(select(curves, model, dose), key=lambda row: row["killing_rate"])
            effect_metrics.append(dict(model_id=model, dose_nM=dose, kd_per_h=p["kd_per_h"],
                                       r_per_h=p["r_per_h"], theoretical_peak_tau_h=p["kernel_peak_tau_h"],
                                       window_peak_rate=observed_peak["killing_rate"],
                                       window_first_peak_data_time_h=observed_peak["time_h"],
                                       theoretical_peak_within_window=(model == 4 and
                                           0 <= paper["origin_in_current_data_hours"] + p["kernel_peak_tau_h"]
                                           <= max(row["time_h"] for row in select(curves, model, dose)))))
            start_metrics.append(dict(model_id=model, dose_nM=dose, starts=len(rows),
                                      successful_exits=sum(row["exitflag"] > 0 for row in rows),
                                      min_SSE=min(row["scaled_SSE"] for row in rows),
                                      max_SSE=max(row["scaled_SSE"] for row in rows),
                                      kd_min=min(row["kd_per_h"] for row in rows),
                                      kd_max=max(row["kd_per_h"] for row in rows),
                                      r_min=min(row["r_per_h"] for row in rows) if model == 4 else None,
                                      r_max=max(row["r_per_h"] for row in rows) if model == 4 else None,
                                      starts_within_1e_minus6_SSE=len(near),
                                      near_best_kd_min=min(row["kd_per_h"] for row in near) if near else None,
                                      near_best_kd_max=max(row["kd_per_h"] for row in near) if near else None,
                                      near_best_r_min=min(row["r_per_h"] for row in near) if near and model == 4 else None,
                                      near_best_r_max=max(row["r_per_h"] for row in near) if near and model == 4 else None,
                                      selected_zero=bool(p["selected_zero"])))
    for row in starts:
        row["kd_near_bound"] = min(math.log(row["kd_per_h"] / 1e-9), math.log(1 / row["kd_per_h"])) < .02
        row["r_near_bound"] = (row["model_id"] == 4 and
                               min(math.log(row["r_per_h"] / 1e-4), math.log(100 / row["r_per_h"])) < .02)
    write("time_form_metrics", metrics)
    if extrema:
        write("time_form_early_local_extrema", extrema)
    write("time_form_effect_metrics", effect_metrics)
    write("time_form_start_metrics", start_metrics)
    write("time_form_start_diagnostics", starts)
    for name in ("summary", "by_dose", "parameters", "residuals", "bounds",
                 "numerical_checks", "peak_checks", "clock_checks"):
        path = source / ("time_form_" + name + ".csv")
        (output / path.name).write_bytes(path.read_bytes())
    for name in ("paper_time_definitions.json", "baseline_manifest.json"):
        (output / name).write_bytes((source / name).read_bytes())
    for name in ("baseline_delay_parameters.csv", "baseline_delay_summary.csv"):
        (output / name).write_bytes((source / name).read_bytes())
    quick_full = []
    if args.quick_results:
        quick = json.loads((args.quick_results / "ci_status.json").read_text())
        assert quick["status"] == "success" and quick["mode"] == "time-quick"
        assert status["mode"] == "time-full" and quick["commit"] == status["commit"]
        for run, ps in ((quick, read("time_form_parameters", args.quick_results)), (status, parameters)):
            for row in ps:
                quick_full.append(dict(mode=run["mode"], run_id=run["run_id"],
                                       model_id=row["model_id"], dose_nM=row["dose_nM"],
                                       kd_per_h=row["kd_per_h"], r_per_h=row["r_per_h"],
                                       scaled_SSE=row["scaled_SSE"], elapsed_seconds=run["elapsed_seconds"]))
        write("time_form_quick_full", quick_full)
    payload = clean(dict(run=status, paper=paper, summary=summary, parameters=parameters,
                         metrics=metrics, effect_metrics=effect_metrics, early_local_extrema=extrema, starts=starts,
                         start_metrics=start_metrics, quick_full=quick_full))
    (output / "time_form_diagnostics.json").write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    styles = {1: ("#15854e", "--", "Archived D"), 2: ("#bd3c2e", ":", "Archived dose beta"),
              3: ("#ba7b00", "-", "Time A diagnostic"), 4: ("#b54a9c", "-", "Time B diagnostic"),
              5: ("#2c8db5", "-.", "Archived delay 1"), 6: ("#41358c", "--", "Archived delay 2")}

    def cells(ax, dose, interval=None):
        for model in range(1, 7):
            rows = select(curves, model, dose)
            if interval:
                rows = [r for r in rows if interval[0] <= r["time_h"] <= interval[1]]
            color, line, label = styles[model]
            ax.plot([r["time_h"] for r in rows], [r["N"] for r in rows], line, color=color, label=label)
        points = select(obs, 1, dose)
        if interval:
            points = [r for r in points if interval[0] <= r["time_h"] <= interval[1]]
        ax.scatter([r["time_h"] for r in points], [r["observed_N"] for r in points], c="black", s=20, zorder=8)
        ax.set(xlabel="Original data time (h)", ylabel="Cell count", title=f"{dose:g} nM")
        ax.grid(alpha=.2)

    fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
    for j, dose in enumerate(doses):
        cells(axes[0, j], dose)
        ax = axes[1, j]
        for model in range(1, 7):
            rows = select(obs, model, dose)
            color, line, _ = styles[model]
            ax.plot([r["time_h"] for r in rows], [r["scaled_residual"] for r in rows], "o" + line, color=color, ms=3)
        ax.axhline(0, color="gray", ls=":")
        ax.set(xlabel="Original data time (h)", ylabel="Scaled residual: prediction - data")
        ax.grid(alpha=.2)
    axes[0, 0].legend(fontsize=8)
    fig.savefig(output / "time_form_comparison.png", dpi=170)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), constrained_layout=True)
    cells(axes[0], 156.25, (300, 500));cells(axes[1], 625, (0, 200))
    axes[0].legend(fontsize=8)
    fig.savefig(output / "time_form_focus.png", dpi=170)
    plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4), constrained_layout=True)
    for ax, dose in zip(axes, doses):
        for model in (3, 4):
            rows = select(curves, model, dose)
            color, line, label = styles[model]
            ax.plot([r["time_h"] for r in rows], [r["killing_rate"] for r in rows], line, color=color, label=label)
        p = select(parameters, 4, dose)[0]
        if math.isfinite(p["kernel_peak_tau_h"]):
            peak = paper["origin_in_current_data_hours"] + p["kernel_peak_tau_h"]
            if min(r["time_h"] for r in rows) <= peak <= max(r["time_h"] for r in rows):
                ax.axvline(peak, color="gray", ls=":", label="B theoretical peak")
            ax.text(.98, .96, f"B: 1/r = {p['kernel_peak_tau_h']:.2f} h", transform=ax.transAxes,
                    ha="right", va="top", fontsize=9)
        ax.set(xlabel="Original data time (h)", ylabel="Killing rate (1/h)", title=f"{dose:g} nM")
        ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.savefig(output / "time_form_killing.png", dpi=170)
    plt.close(fig)
    fig, axes = plt.subplots(2, 3, figsize=(14, 7), constrained_layout=True)
    transient_scores = [math.log10(row["scaled_SSE"]) for row in starts if row["model_id"] == 4]
    for j, dose in enumerate(doses):
        for i, model in enumerate((3, 4)):
            rows = select(starts, model, dose)
            best = select(parameters, model, dose)[0]
            ax = axes[i, j]
            if model == 3:
                ax.scatter([row["start"] for row in rows],
                           [row["kd_per_h"] / best["kd_per_h"] for row in rows], color="#ba7b00")
                ax.axhline(1, color="gray", ls=":")
                ratios = [row["kd_per_h"] / best["kd_per_h"] for row in rows]
                span = max(2e-5, 1.4 * max(abs(value - 1) for value in ratios))
                ax.set(ylim=(1 - span, 1 + span), xlabel="Start number", ylabel="kd_A / best kd_A")
                ax.ticklabel_format(axis="y", style="plain", useOffset=False)
                ax.set_xticks([row["start"] for row in rows])
                score_span = max(row["scaled_SSE"] for row in rows) - min(row["scaled_SSE"] for row in rows)
                ax.text(.03, .96, f"Best kd_A = {best['kd_per_h']:.8f} /h\n"
                        f"Best SSE = {best['scaled_SSE']:.6f}\nSSE range width = {score_span:.2e}",
                        transform=ax.transAxes, va="top", fontsize=9)
            else:
                points = ax.scatter([row["r_per_h"] for row in rows], [row["kd_per_h"] for row in rows],
                                    c=[math.log10(row["scaled_SSE"]) for row in rows], cmap="viridis", s=40,
                                    vmin=min(transient_scores), vmax=max(transient_scores))
                near_ids = []
                for row in rows:
                    if row["scaled_SSE"] <= best["scaled_SSE"] + 1e-6:
                        near_ids.append(str(int(row["start"])))
                    else:
                        ax.annotate(str(int(row["start"])), (row["r_per_h"], row["kd_per_h"]),
                                    fontsize=8, xytext=(4, -12 - 10 * (int(row["start"]) % 3)),
                                    textcoords="offset points")
                ax.set(xscale="log", yscale="log", xlim=(5e-5, 200), ylim=(5e-10, 2),
                       xlabel="r (1/h)", ylabel="kd_B (1/h)")
                ax.text(.03, .04, "Near-best starts: " + ",".join(near_ids),
                        transform=ax.transAxes, fontsize=9)
            ax.set_title(f"Time {'A' if model == 3 else 'B'}: {dose:g} nM")
            ax.grid(alpha=.2)
    fig.colorbar(points, ax=axes[1, :].tolist(), label="log10(SSE)", shrink=.85)
    fig.savefig(output / "time_form_starts.png", dpi=170)
    plt.close(fig)
    print(json.dumps(payload, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
