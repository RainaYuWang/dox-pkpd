#!/usr/bin/env python3
"""Analyze collected MATLAB delay quick/full outputs without refitting."""
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
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("reports"))
    parser.add_argument("--quick-results", type=Path, help="Optional matching quick joint_results for comparison")
    args = parser.parse_args()
    source, out = args.results, args.output
    out.mkdir(parents=True, exist_ok=True)

    def read(name):
        with (source / (name + ".csv")).open(newline="") as handle:
            return [{key: (value if key in ("model", "status") else float(value or "nan"))
                     for key, value in row.items()} for row in csv.DictReader(handle)]

    def write(name, rows):
        with (out / (name + ".csv")).open("w", newline="") as handle:
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

    status = json.loads((source / "ci_status.json").read_text())
    assert status["status"] == "success" and status["mode"] in ("delay-quick", "delay-full")
    summary, params, by_dose, obs, curves = [read(name) for name in
                                          ("delay_summary", "delay_parameters", "delay_by_dose",
                                           "delay_residuals", "delay_predictions")]
    active = [int(row["model_id"]) for row in summary if row["status"] == "completed"]
    assert len(obs) == 32 * len(active)
    paper = json.loads((source / "paper_time_definitions.json").read_text())
    if not paper["verified"]:
        assert 3 not in active and 4 not in active
    doses = [156.25, 312.5, 625]

    def select(rows, model, dose):
        return [row for row in rows if row["model_id"] == model and row["dose_nM"] == dose]

    metrics, peaks = [], []
    for dose in doses:
        for model in active:
            points, prediction = select(obs, model, dose), select(curves, model, dose)
            rr = np.array([row["residual"] for row in points])
            late = np.array([row["residual"] for row in points if row["time_h"] >= 300])
            early = min((row for row in prediction if row["time_h"] <= 200), key=lambda row: row["N"])
            minimum = min(prediction, key=lambda row: row["N"])
            expected = select(by_dose, model, dose)[0]["scaled_SSE"]
            assert abs(sum(row["scaled_residual"] ** 2 for row in points) - expected) < 1e-10
            metrics.append(dict(model_id=model, dose_nM=dose, scaled_SSE=expected,
                                RMSE_cells=float(np.sqrt(np.mean(rr ** 2))),
                                early_min_N=early["N"], early_min_time_h=early["time_h"],
                                min_N=minimum["N"], min_time_h=minimum["time_h"],
                                late_count=len(late), late_mean_residual=float(np.mean(late)),
                                late_RMSE=float(np.sqrt(np.mean(late ** 2))),
                                final_residual=points[-1]["residual"]))
            if model in (5, 6):
                p = select(by_dose, model, dose)[0]
                cb = max(prediction, key=lambda row: row["CB"])
                e2 = max(prediction, key=lambda row: row["E2"])
                kill = max(prediction, key=lambda row: row["killing_rate"])
                peaks.append(dict(model_id=model, dose_nM=dose, kBF=p["kBF"],
                                  alpha=p["amplitude"], r_per_h=p["r_per_h"],
                                  kernel_peak_delay_h=1 / p["r_per_h"], kernel_mean_delay_h=2 / p["r_per_h"],
                                  CB_peak=cb["CB"], CB_peak_time_h=cb["time_h"],
                                  E2_peak=e2["E2"], E2_peak_time_h=e2["time_h"],
                                  killing_peak_per_h=kill["killing_rate"], killing_peak_time_h=kill["time_h"],
                                  killing_peak_at_window_end=kill["time_h"] == prediction[-1]["time_h"]))
    write("delay_metrics", metrics)
    write("delay_effect_peaks", peaks)
    changes = []
    for dose in doses:
        baseline = select(curves, 1, dose)
        for model in (5, 6):
            prediction = select(curves, model, dose)
            assert [row["time_h"] for row in baseline] == [row["time_h"] for row in prediction]
            changes.append(dict(model_id=model, dose_nM=dose,
                                max_absolute_N_change_from_D=max(abs(a["N"] - b["N"])
                                                                 for a, b in zip(baseline, prediction))))
    write("delay_curve_changes", changes)
    quick_comparison = []
    if args.quick_results:
        quick = json.loads((args.quick_results / "ci_status.json").read_text())
        assert quick["status"] == "success" and quick["mode"] == "delay-quick"
        assert quick["commit"] == status["commit"] and status["mode"] == "delay-full"
        for run in (quick, status):
            for stage in (1, 2):
                p = run[f"delay_stage{stage}_parameters"]
                quick_comparison.append(dict(mode=run["mode"], run_id=run["run_id"], stage=stage,
                                             kBF=p[0], alpha=p[1], r_per_h=p[2],
                                             scaled_SSE=run[f"delay_stage{stage}_SSE"],
                                             elapsed_seconds=run["elapsed_seconds"]))
        write("delay_quick_full", quick_comparison)
    start_diagnostics = []
    start_sets = []
    for stage in (1, 2):
        rows = read(f"delay_stage{stage}_starts")
        start_sets.append(rows)
        bounds = {"alpha": (1e-9, .1), "r_per_h": (1e-4, 100)}
        if stage == 2:
            bounds["kBF"] = (1e-7, 10)
        for index, row in enumerate(rows, 1):
            near = [key for key, (low, high) in bounds.items()
                    if math.log(row[key] / low) < .02 or math.log(high / row[key]) < .02]
            start_diagnostics.append(dict(stage=stage, start=index, **row, near_bounds=near))
    profile = []
    ranges = []
    if (source / "delay_stage2_r_profile.csv").exists():
        profile = read("delay_stage2_r_profile")
        optimum = params[1]["scaled_SSE"]
        for percent in (1, 5, 10):
            accepted = [row for row in profile if row["scaled_SSE"] <= optimum * (1 + percent / 100)]
            if accepted:
                ranges.append(dict(SSE_increase_pct=percent,
                                   r_min=min(row["r_per_h"] for row in accepted),
                                   r_max=max(row["r_per_h"] for row in accepted),
                                   kBF_min=min(row["kBF"] for row in accepted),
                                   kBF_max=max(row["kBF"] for row in accepted),
                                   alpha_min=min(row["alpha"] for row in accepted),
                                   alpha_max=max(row["alpha"] for row in accepted)))
    payload = clean(dict(run=status, summary=summary, parameters=params, by_dose=by_dose,
                         metrics=metrics, peaks=peaks, curve_changes=changes,
                         quick_full=quick_comparison, starts=start_diagnostics,
                         descriptive_profile_ranges=ranges, paper=paper))
    (out / "delay_diagnostics.json").write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    for name in ("delay_summary", "delay_parameters", "delay_by_dose", "delay_residuals",
                 "delay_stage1_starts", "delay_stage2_starts", "delay_stage2_r_profile",
                 "delay_numerical_checks", "delay_matrix_exponential_checks", "delay_bounds", "time_form_starts"):
        path = source / (name + ".csv")
        if path.exists():
            (out / path.name).write_bytes(path.read_bytes())

    styles = {1: ("#15854e", "--", "Archived D"), 2: ("#bd3c2e", ":", "Dose-specific beta"),
              3: ("#ba7b00", "--", "Time-form A diagnostic"), 4: ("#b54a9c", "-.", "Time-form B diagnostic"),
              5: ("#2c8db5", "-.", "Delay stage 1"), 6: ("#41358c", "-", "Delay stage 2")}

    def cell_panel(ax, dose, interval=None):
        for model in active:
            prediction = select(curves, model, dose)
            if interval:
                prediction = [row for row in prediction if interval[0] <= row["time_h"] <= interval[1]]
            color, line, label = styles[model]
            ax.plot([row["time_h"] for row in prediction], [row["N"] for row in prediction],
                    line, color=color, label=label, lw=1.7)
        points = select(obs, 1, dose)
        if interval:
            points = [row for row in points if interval[0] <= row["time_h"] <= interval[1]]
        ax.scatter([row["time_h"] for row in points], [row["observed_N"] for row in points],
                   color="black", s=23, zorder=6, label="Digitized data")
        ax.set(xlabel="Hours", ylabel="Cell count", title=f"{dose:g} nM")
        ax.grid(alpha=.2)

    fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
    for column, dose in enumerate(doses):
        ax, rx = axes[:, column]
        cell_panel(ax, dose)
        for model in active:
            points = select(obs, model, dose)
            color, line, _ = styles[model]
            rx.plot([row["time_h"] for row in points], [row["scaled_residual"] for row in points],
                    "o" + line, color=color, lw=1.3, ms=4)
        rx.axhline(0, color="gray", ls=":")
        rx.set(xlabel="Hours", ylabel="Scaled residual (prediction - data)")
        rx.grid(alpha=.2)
    axes[0, 0].legend(fontsize=8)
    fig.suptitle(f"CB-driven delay study ({status['mode']}); time-form diagnostics "
                 + ("included" if paper["verified"] else "not run: paper origin unverified"), fontsize=12)
    fig.savefig(out / "delay_comparison.png", dpi=170)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), constrained_layout=True)
    cell_panel(axes[0], 156.25, (300, 500))
    axes[0].set_title("156.25 nM: late observations")
    cell_panel(axes[1], 625, (0, 200))
    axes[1].set_title("625 nM: decline and early recovery")
    axes[0].legend(fontsize=8)
    fig.savefig(out / "delay_focus.png", dpi=170)
    plt.close(fig)
    fig, axes = plt.subplots(3, 3, figsize=(14, 10), constrained_layout=True)
    for i, dose in enumerate(doses):
        for j, field in enumerate(("CB", "E2", "killing_rate")):
            ax = axes[i, j]
            for model in ((5, 6) if field == "E2" else (1, 5, 6)):
                prediction = select(curves, model, dose)
                color, line, label = styles[model]
                ax.plot([row["time_h"] for row in prediction], [row[field] for row in prediction],
                        line, color=color, label=label, lw=1.6)
            ax.set(title=f"{dose:g} nM: {field}", xlabel="Hours",
                   ylabel="1/h" if field == "killing_rate" else "nM")
            ax.grid(alpha=.2)
            if i == 0:
                ax.legend(fontsize=8)
    fig.savefig(out / "delay_effect_states.png", dpi=170)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    for model in (5, 6):
        prediction = [row for row in select(curves, model, 625) if row["time_h"] <= 12]
        color, line, label = styles[model]
        axes[0].plot([row["time_h"] for row in prediction], [row["E2"] for row in prediction],
                     line, color=color, label=label + " E2")
        axes[1].plot([row["time_h"] for row in prediction], [row["killing_rate"] for row in prediction],
                     line, color=color, label=label)
    baseline = [row for row in select(curves, 1, 625) if row["time_h"] <= 12]
    for ax, field, label in zip(axes, ("CB", "killing_rate"), ("D CB", "D killing")):
        ax.plot([row["time_h"] for row in baseline], [row[field] for row in baseline],
                "--", color=styles[1][0], label=label)
        ax.axvline(6, color="gray", ls=":", label="Washout")
        ax.set(xlabel="Hours", ylabel="nM" if field == "CB" else "1/h", title="625 nM: first 12 h")
        ax.legend(fontsize=8)
        ax.grid(alpha=.2)
    fig.savefig(out / "delay_effect_early.png", dpi=170)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3), constrained_layout=True)
    for ax, stage, rows, field in zip(axes, (1, 2), start_sets, ("alpha", "kBF")):
        points = ax.scatter([row["r_per_h"] for row in rows], [row[field] for row in rows],
                            c=[np.log10(row["scaled_SSE"]) for row in rows], cmap="viridis", s=45)
        for index, row in enumerate(rows, 1):
            ax.annotate(str(index), (row["r_per_h"], row[field]), fontsize=7, xytext=(3, 3), textcoords="offset points")
        ax.set(xscale="log", yscale="log", xlabel="r (1/h)", ylabel=field,
               title=f"Stage {stage}: all starts (labels = start index)")
        ax.grid(alpha=.2)
        fig.colorbar(points, ax=ax, label="log10(SSE)")
    fig.savefig(out / "delay_starts.png", dpi=170)
    plt.close(fig)
    if profile:
        fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), constrained_layout=True)
        r = [row["r_per_h"] for row in profile]
        axes[0].semilogx(r, [row["scaled_SSE"] for row in profile], "o-", ms=3)
        axes[0].plot(params[1]["r_per_h"], params[1]["scaled_SSE"], "r*", ms=12)
        axes[0].set(ylabel="Profile scaled SSE", title="Stage 2: kBF and alpha re-fitted")
        for ax, field in zip(axes[1:], ("kBF", "alpha")):
            ax.loglog(r, [row[field] for row in profile], "o-", ms=3)
            ax.set(ylabel=field, title=f"Compensation in {field}")
        for ax in axes:
            ax.set_xlabel("Fixed r (1/h)")
            ax.grid(alpha=.2)
        fig.suptitle("Conditional r profile; exploratory bounds, not a confidence interval", fontsize=12)
        fig.savefig(out / "delay_stage2_r_profile.png", dpi=170)
        plt.close(fig)
        fig, ax = plt.subplots(figsize=(8, 4), constrained_layout=True)
        near = [row for row in profile if row["r_per_h"] >= .1]
        best = params[1]["scaled_SSE"]
        ax.semilogx([row["r_per_h"] for row in near],
                    [100 * (row["scaled_SSE"] / best - 1) for row in near], "o-")
        ax.axhline(1, color="gray", ls=":", label="1% descriptive threshold (not a CI)")
        ax.plot(params[1]["r_per_h"], 0, "r*", ms=12, label="Best joint fit")
        ax.set(xlabel="Fixed r (1/h)", ylabel="SSE increase from best (%)",
               title="Stage 2 profile: fast-response plateau", ylim=(-.05, 2.5))
        ax.legend(fontsize=8)
        ax.grid(alpha=.2)
        fig.savefig(out / "delay_r_profile_focus.png", dpi=170)
        plt.close(fig)
    print(json.dumps(payload, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
