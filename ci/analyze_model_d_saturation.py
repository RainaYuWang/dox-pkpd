#!/usr/bin/env python3
"""Reproduce saturation comparison, residual and Ks-profile diagnostics from MATLAB outputs."""
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
    prefix = "model_d_saturation_"

    def read(name):
        with (source / (prefix + name + ".csv")).open(newline="") as handle:
            return [{key: (value if key == "model" else float(value)) for key, value in row.items()}
                    for row in csv.DictReader(handle)]

    def write(name, rows):
        with (output / (prefix + name + ".csv")).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    status = json.loads((source / "ci_status.json").read_text())
    assert status["status"] == "success" and status["mode"] == "saturation"
    summary, parameters, residuals, predictions, by_dose, profile, starts, limits = [
        read(name) for name in ("summary", "parameters", "residuals", "predictions",
                               "by_dose", "profile", "starts", "linear_limit")]
    assert len(summary) == 3 and len(residuals) == 96 and len(starts) == 9
    doses, models = [156.25, 312.5, 625], [1, 2, 3]

    def select(rows, dose, model):
        return [row for row in rows if row["dose_nM"] == dose and row["model_id"] == model]

    metrics = []
    for dose in doses:
        for model in models:
            obs, curve = select(residuals, dose, model), select(predictions, dose, model)
            error = np.array([row["residual"] for row in obs])
            expected = select(by_dose, dose, model)[0]["scaled_SSE"]
            assert abs(sum(row["scaled_residual"] ** 2 for row in obs) - expected) < 1e-10
            assert all(np.isfinite(row["N"]) and row["N"] > 0 for row in curve)
            minimum = min(curve, key=lambda row: row["N"])
            early = min((row for row in curve if row["time_h"] <= 200), key=lambda row: row["N"])
            late = np.array([row["residual"] for row in obs if row["time_h"] >= 300])
            metrics.append(dict(dose_nM=dose, model_id=model, scaled_SSE=expected,
                                RMSE_cells=float(np.sqrt(np.mean(error ** 2))),
                                mean_residual_cells=float(np.mean(error)),
                                min_N=minimum["N"], min_time_h=minimum["time_h"],
                                early_min_N=early["N"], early_min_time_h=early["time_h"],
                                early_drop_pct_from_initial=100 * (1 - early["N"] / curve[0]["N"]),
                                late_observation_count=len(late),
                                late_mean_residual_cells=float(np.mean(late)),
                                late_RMSE_cells=float(np.sqrt(np.mean(late ** 2))),
                                late_min_residual_cells=float(min(late)),
                                late_max_residual_cells=float(max(late)),
                                final_residual_cells=obs[-1]["residual"]))
    write("metrics", metrics)
    for name in ("summary", "parameters", "by_dose", "profile", "starts", "linear_limit",
                 "bounds", "residuals"):
        (output / (prefix + name + ".csv")).write_bytes((source / (prefix + name + ".csv")).read_bytes())
    best = parameters[0]
    finite = [row for row in profile if np.isfinite(row["Ks_nM"])]
    profile_min = min(row["scaled_SSE"] for row in profile)
    ranges = []
    for percent in (1, 5, 10):
        accepted = [row["Ks_nM"] for row in finite if row["scaled_SSE"] <= profile_min * (1 + percent / 100)]
        ranges.append(dict(SSE_increase_pct=percent, sampled_Ks_min=min(accepted),
                           sampled_Ks_max=max(accepted), includes_linear_limit=any(
                               not np.isfinite(row["Ks_nM"]) and row["scaled_SSE"] <= profile_min * (1 + percent / 100)
                               for row in profile)))
    near = []
    for index, row in enumerate(starts, 1):
        for name, low, high in (("beta", 1e-9, .1), ("Ks_nM", 1e-3, 1e8)):
            if math.log(row[name] / low) < .02 or math.log(high / row[name]) < .02:
                near.append(dict(start=index, parameter=name, value=row[name]))
    diagnostics = dict(run=status, best=best, profile_min_SSE=profile_min,
                       descriptive_profile_ranges=ranges,
                       beta_bound_profile_points=[row for row in profile
                                                  if row["beta_near_lower"] or row["beta_near_upper"]],
                       multistart_flags=[int(row["exitflag"]) for row in starts],
                       multistart_near_bounds=near,
                       start_SSE_min=min(row["scaled_SSE"] for row in starts),
                       start_SSE_max=max(row["scaled_SSE"] for row in starts))
    (output / (prefix + "diagnostics.json")).write_text(json.dumps(diagnostics, indent=2) + "\n")

    styles = {1: ("#15854e", "--", "Linear shared beta D"),
              2: ("#bd3c2e", ":", "Dose-specific beta diagnostic"),
              3: ("#235db5", "-", "Shared saturating beta, Ks")}
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
    for column, dose in enumerate(doses):
        ax, rx = axes[:, column]
        points = select(residuals, dose, 1)
        ax.scatter([row["time_h"] for row in points], [row["observed_N"] for row in points],
                   color="black", s=23, zorder=4, label="Digitized data")
        for model, (color, style, label) in styles.items():
            curve, obs = select(predictions, dose, model), select(residuals, dose, model)
            ax.plot([row["time_h"] for row in curve], [row["N"] for row in curve],
                    style, color=color, lw=1.8, label=label)
            rx.plot([row["time_h"] for row in obs], [row["scaled_residual"] for row in obs],
                    "o" + style, color=color, lw=1.4, ms=4)
        ax.set(title=f"{dose:g} nM", xlabel="Hours", ylabel="Cell count")
        rx.set(xlabel="Hours", ylabel="Scaled residual (prediction - data)")
        rx.axhline(0, color="gray", ls=":")
        for panel in (ax, rx):
            panel.grid(alpha=.2)
        if column == 0:
            ax.legend(fontsize=8)
    fig.suptitle("Fixed reference D PK, growth and kBF: PD-form comparison", fontsize=13)
    fig.savefig(output / (prefix + "comparison.png"), dpi=170)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), constrained_layout=True)
    for ax, dose in zip(axes, doses):
        for model, (color, style, label) in styles.items():
            curve = [row for row in select(predictions, dose, model) if row["time_h"] <= 200]
            ax.plot([row["time_h"] for row in curve], [row["N"] for row in curve],
                    style, color=color, lw=1.8, label=label)
        obs = [row for row in select(residuals, dose, 1) if row["time_h"] <= 200]
        ax.scatter([row["time_h"] for row in obs], [row["observed_N"] for row in obs],
                   color="black", zorder=4, label="Digitized data")
        ax.set(title=f"{dose:g} nM: 0-200 h", xlabel="Hours", ylabel="Cell count")
        ax.grid(alpha=.2)
    axes[0].legend(fontsize=8)
    fig.savefig(output / (prefix + "early.png"), dpi=170)
    plt.close(fig)

    ks = np.array([row["Ks_nM"] for row in finite])
    ss = np.array([row["scaled_SSE"] for row in finite])
    bb = np.array([row["beta"] for row in finite])
    linear_sse = next(row["scaled_SSE"] for row in profile if not np.isfinite(row["Ks_nM"]))
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), constrained_layout=True)
    axes[0].semilogx(ks, ss, "o-", ms=3)
    axes[0].axhline(linear_sse, color="gray", ls="--", label="Re-fitted linear limit")
    axes[0].plot(best["Ks_nM"], best["scaled_SSE"], "r*", ms=12, label="Joint optimum")
    axes[0].set(xlabel="Fixed Ks (nM)", ylabel="Profile scaled SSE", title="Entire exploratory range")
    axes[0].legend(fontsize=8)
    axes[1].semilogx(ks, 100 * (ss / profile_min - 1), "o-", ms=3)
    close = ks[ss <= profile_min * 1.1]
    axes[1].set(xlim=(max(min(ks), min(close) / 2), min(max(ks), max(close) * 2)), ylim=(-.5, 10),
                xlabel="Fixed Ks (nM)", ylabel="SSE increase (%)", title="Descriptive valley; not a confidence interval")
    for percent in (1, 5):
        axes[1].axhline(percent, color="gray", ls=":")
    axes[2].loglog(ks, bb, "o-", ms=3)
    axes[2].set(xlabel="Fixed Ks (nM)", ylabel="Re-fitted beta (1/(nM h))", title="Beta compensation across Ks")
    for ax in axes:
        ax.grid(alpha=.2)
    fig.savefig(output / (prefix + "profile.png"), dpi=170)
    plt.close(fig)
    print(json.dumps(dict(diagnostics=diagnostics, summary=summary, by_dose=by_dose,
                          metrics=metrics, linear_limit=limits), indent=2))


if __name__ == "__main__":
    main()
