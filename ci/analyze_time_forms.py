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
    args = parser.parse_args()
    source, output = args.results, args.output
    status = json.loads((source / "ci_status.json").read_text())
    paper = json.loads((source / "paper_time_definitions.json").read_text())
    assert status["status"] == "success" and status["mode"] in ("time-quick", "time-full")
    assert paper["verified"] is True and status["old_models_refitted"] is False

    def read(name):
        with (source / (name + ".csv")).open(newline="") as handle:
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
    metrics, start_metrics = [], []
    for dose in doses:
        for model in range(1, 7):
            predicted, points = select(curves, model, dose), select(obs, model, dose)
            late = [p["residual"] for p in points if p["time_h"] >= 300]
            early = min((p for p in predicted if p["time_h"] <= 200), key=lambda p: p["N"])
            metrics.append(dict(model_id=model, dose_nM=dose,
                                scaled_SSE=sum(p["scaled_residual"] ** 2 for p in points),
                                RMSE_cells=float(np.sqrt(np.mean([p["residual"] ** 2 for p in points]))),
                                early_min_N=early["N"], early_min_data_time_h=early["time_h"],
                                late_count=len(late), late_mean_residual=float(np.mean(late)),
                                late_RMSE=float(np.sqrt(np.mean(np.array(late) ** 2))),
                                final_residual=points[-1]["residual"]))
        for model in (3, 4):
            rows = select(starts, model, dose)
            p = select(parameters, model, dose)[0]
            start_metrics.append(dict(model_id=model, dose_nM=dose, starts=len(rows),
                                      successful_exits=sum(row["exitflag"] > 0 for row in rows),
                                      min_SSE=min(row["scaled_SSE"] for row in rows),
                                      max_SSE=max(row["scaled_SSE"] for row in rows),
                                      kd_min=min(row["kd_per_h"] for row in rows),
                                      kd_max=max(row["kd_per_h"] for row in rows),
                                      r_min=min(row["r_per_h"] for row in rows) if model == 4 else None,
                                      r_max=max(row["r_per_h"] for row in rows) if model == 4 else None,
                                      selected_zero=bool(p["selected_zero"])))
    for row in starts:
        row["kd_near_bound"] = min(math.log(row["kd_per_h"] / 1e-9), math.log(1 / row["kd_per_h"])) < .02
        row["r_near_bound"] = (row["model_id"] == 4 and
                               min(math.log(row["r_per_h"] / 1e-4), math.log(100 / row["r_per_h"])) < .02)
    write("time_form_metrics", metrics)
    write("time_form_start_metrics", start_metrics)
    write("time_form_start_diagnostics", starts)
    for name in ("summary", "by_dose", "parameters", "residuals", "bounds",
                 "numerical_checks", "peak_checks", "clock_checks"):
        path = source / ("time_form_" + name + ".csv")
        (output / path.name).write_bytes(path.read_bytes())
    for name in ("paper_time_definitions.json", "baseline_manifest.json"):
        (output / name).write_bytes((source / name).read_bytes())
    payload = clean(dict(run=status, paper=paper, summary=summary, parameters=parameters,
                         metrics=metrics, starts=starts, start_metrics=start_metrics))
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
            ax.axvline(paper["origin_in_current_data_hours"] + p["kernel_peak_tau_h"], color="gray", ls=":", label="B theoretical peak")
        ax.set(xlabel="Original data time (h)", ylabel="Killing rate (1/h)", title=f"{dose:g} nM")
        ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.savefig(output / "time_form_killing.png", dpi=170)
    plt.close(fig)
    print(json.dumps(payload, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
