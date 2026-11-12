#!/usr/bin/env python3
"""Summarize an actual MATLAB dose-beta diagnostic; no fitting is performed here."""
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
    args = parser.parse_args()
    source, output = args.results, args.output
    output.mkdir(parents=True, exist_ok=True)
    prefix = "model_d_dose_beta_"

    def read(suffix):
        with (source / (prefix + suffix + ".csv")).open(newline="") as handle:
            return [{key: float(value) for key, value in row.items()}
                    for row in csv.DictReader(handle)]

    def write(suffix, rows):
        with (output / (prefix + suffix + ".csv")).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    summary, residuals, predictions, starts = [read(name) for name in
                                              ("summary", "residuals", "predictions", "starts")]
    status = json.loads((source / "ci_status.json").read_text())
    assert status["status"] == "success" and status["mode"] == "dose-beta"
    assert len(summary) == 3 and len(residuals) == 64 and len(starts) == 24
    details, diagnostics = [], []
    for row in summary:
        dose = row["dose_nM"]
        row["beta_ratio_to_shared"] = row["dose_beta"] / row["shared_beta"]
        row["SSE_reduction_pct"] = 100 * (1 - row["dose_scaled_SSE"] / row["shared_scaled_SSE"])
        for model in (1, 2):
            obs = [item for item in residuals if item["dose_nM"] == dose and item["model_id"] == model]
            pred = [item for item in predictions if item["dose_nM"] == dose and item["model_id"] == model]
            errors = np.array([item["residual"] for item in obs])
            scaled_sse = sum(item["scaled_residual"] ** 2 for item in obs)
            expected = row["shared_scaled_SSE"] if model == 1 else row["dose_scaled_SSE"]
            assert abs(scaled_sse - expected) < 1e-10
            assert all(np.isfinite(item["N"]) and item["N"] > 0 for item in pred)
            minimum = min(pred, key=lambda item: item["N"])
            largest = obs[int(np.argmax(abs(errors)))]
            details.append(dict(dose_nM=dose, model_id=model, scaled_SSE=scaled_sse,
                                RMSE_cells=float(np.sqrt(np.mean(errors ** 2))),
                                mean_residual_cells=float(np.mean(errors)),
                                min_residual_cells=float(min(errors)),
                                max_residual_cells=float(max(errors)),
                                max_abs_residual_time_h=largest["time_h"],
                                min_predicted_N=minimum["N"], min_predicted_time_h=minimum["time_h"],
                                initial_N=pred[0]["N"],
                                min_observed_N=min(item["observed_N"] for item in obs)))
        dose_starts = [item for item in starts if item["dose_nM"] == dose]
        diagnostics.append(dict(dose_nM=dose,
                                SSE_min=min(item["scaled_SSE"] for item in dose_starts),
                                SSE_max=max(item["scaled_SSE"] for item in dose_starts),
                                flags=[int(item["exitflag"]) for item in dose_starts],
                                beta_min=min(item["beta"] for item in dose_starts),
                                beta_max=max(item["beta"] for item in dose_starts),
                                near_bound_starts=[int(item["start"]) for item in dose_starts
                                                   if math.log(item["beta"] / 1e-9) < .02
                                                   or math.log(.1 / item["beta"]) < .02]))
    write("summary", summary)
    write("metrics", details)
    for name in ("residuals", "starts", "fixed"):
        (output / (prefix + name + ".csv")).write_bytes((source / (prefix + name + ".csv")).read_bytes())
    for path in sorted((source / "dose_beta_checkpoints").glob("dose_*/local_objective.csv")):
        (output / (prefix + path.parent.name + "_local_objective.csv")).write_bytes(path.read_bytes())
    (output / (prefix + "diagnostics.json")).write_text(
        json.dumps(dict(run=status, starts=diagnostics), indent=2) + "\n")

    fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
    styles = {1: ("#15854e", "--", "Shared beta D"),
              2: ("#bd3c2e", "-", "Dose-specific beta diagnostic")}
    for column, row in enumerate(summary):
        dose = row["dose_nM"]
        ax, rx = axes[:, column]
        points = [item for item in residuals if item["dose_nM"] == dose and item["model_id"] == 1]
        ax.scatter([item["time_h"] for item in points], [item["observed_N"] for item in points],
                   s=23, color="black", zorder=4, label="Digitized data")
        for model, (color, style, label) in styles.items():
            curve = [item for item in predictions if item["dose_nM"] == dose and item["model_id"] == model]
            points = [item for item in residuals if item["dose_nM"] == dose and item["model_id"] == model]
            ax.plot([item["time_h"] for item in curve], [item["N"] for item in curve],
                    style, color=color, label=label, lw=1.8)
            rx.plot([item["time_h"] for item in points], [item["scaled_residual"] for item in points],
                    "o" + style, color=color, label=label, lw=1.5, ms=4)
        ax.set(title=f"{dose:g} nM", xlabel="Hours", ylabel="Cell count")
        rx.set(xlabel="Hours", ylabel="Scaled residual (prediction - data)")
        rx.axhline(0, color="gray", ls=":")
        for panel in (ax, rx):
            panel.grid(alpha=.2)
        if column == 0:
            ax.legend(fontsize=8)
    fig.suptitle("Fixed reference D PK, growth and kBF; diagnostic fit of one beta per dose", fontsize=13)
    fig.savefig(output / (prefix + "comparison.png"), dpi=170)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
    for model, (color, style, label) in styles.items():
        curve = [item for item in predictions if item["dose_nM"] == 625
                 and item["model_id"] == model and item["time_h"] <= 200]
        ax.plot([item["time_h"] for item in curve], [item["N"] for item in curve],
                style, color=color, label=label, lw=2)
    points = [item for item in residuals if item["dose_nM"] == 625
              and item["model_id"] == 1 and item["time_h"] <= 200]
    ax.scatter([item["time_h"] for item in points], [item["observed_N"] for item in points],
               color="black", label="Digitized data", zorder=4)
    ax.set(title="625 nM: decline and early recovery (0-200 h)", xlabel="Hours", ylabel="Cell count")
    ax.grid(alpha=.2)
    ax.legend()
    fig.savefig(output / (prefix + "625_early.png"), dpi=170)
    plt.close(fig)
    print(json.dumps(dict(run_id=status["run_id"], summary=summary, metrics=details,
                          starts=diagnostics), indent=2))


if __name__ == "__main__":
    main()
