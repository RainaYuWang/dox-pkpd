#!/usr/bin/env python3
"""Reproduce model-D sensitivity tables and plots from a collected MATLAB run."""
import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True,
                        help="Extracted joint_results directory")
    parser.add_argument("--output", type=Path, default=Path("reports"))
    args = parser.parse_args()
    source, output = args.results, args.output
    output.mkdir(parents=True, exist_ok=True)

    def read(name):
        with (source / name).open(newline="") as handle:
            return list(csv.DictReader(handle))

    def write(name, rows):
        with (output / name).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    summary = read("model_d_sensitivity_summary.csv")
    starts = read("model_d_sensitivity_starts.csv")
    predictions = read("model_d_sensitivity_predictions.csv")
    observations = read("model_d_sensitivity_residuals.csv")
    for row in summary:
        for key in row:
            if key != "varied_parameter":
                row[key] = float(row[key])
    references = [row for row in summary if row["multiplier"] == 1]
    assert len(summary) == 9 and len(starts) == 72 and len(references) == 3
    reference = references[0]
    reference_id = int(reference["case_id"])
    doses = [156.25, 312.5, 625]

    def curve(case, dose):
        rows = [row for row in predictions
                if int(row["case_id"]) == case and float(row["dose_nM"]) == dose]
        array = np.array([[float(row[key]) for key in ("time_h", "N", "beta_CB")]
                          for row in rows])
        assert len(array) > 0 and np.all(np.isfinite(array))
        assert np.all(np.diff(array[:, 0]) > 0) and np.all(array[:, 1] > 0)
        return array

    for row in references:
        assert all(row[key] == reference[key] for key in ("scaled_SSE", "kBF", "beta"))
        for dose in doses:
            assert np.array_equal(curve(int(row["case_id"]), dose),
                                  curve(reference_id, dose))
    original = loadmat(source / "model_d_fit.mat", simplify_cells=True)
    sensitivity = loadmat(source / "model_d_sensitivity.mat", simplify_cells=True)
    assert np.array_equal(original["gD"], sensitivity["g"])
    assert sensitivity["cfg"]["VI"] == .004629 and sensitivity["cfg"]["VE"] == 250

    metrics, diagnostics = [], []
    for row in summary:
        case = int(row["case_id"])
        row["SSE_change_pct"] = 100 * (row["scaled_SSE"] / reference["scaled_SSE"] - 1)
        row["kBF_ratio"] = row["kBF"] / reference["kBF"]
        row["beta_ratio"] = row["beta"] / reference["beta"]
        a = row["fixed_kFE"] + row["fixed_kFB"] + row["kBF"]
        b = row["fixed_kFE"] * row["kBF"]
        row["slow_decay_rate"] = 2 * b / (a + np.sqrt(a * a - 4 * b))
        for dose in doses:
            values, baseline = curve(case, dose), curve(reference_id, dose)
            assert np.array_equal(values[:, 0], baseline[:, 0])
            time = values[:, 0]
            metric = dict(case_id=case, varied_parameter=row["varied_parameter"],
                          multiplier=row["multiplier"], dose_nM=dose)
            for name, column in (("N", 1), ("E", 2)):
                delta = np.abs(values[:, column] - baseline[:, column])
                metric[name + "_max_abs_diff"] = float(max(delta))
                metric[name + "_max_diff_pct_of_ref_peak"] = float(
                    100 * max(delta) / max(baseline[:, column]))
            for name, mask in (("early", time <= 12), ("late", time >= 24)):
                metric["E_" + name + "_diff_pct_of_full_ref_peak"] = float(
                    100 * max(abs(values[mask, 2] - baseline[mask, 2])) / max(baseline[:, 2]))
            metric["E_peak"] = float(max(values[:, 2]))
            metric["E_peak_time"] = float(time[np.argmax(values[:, 2])])
            metric["E_integral"] = float(np.trapezoid(values[:, 2], time))
            metric["E_integral_change_pct"] = float(
                100 * (np.trapezoid(values[:, 2], time) / np.trapezoid(baseline[:, 2], time) - 1))
            metrics.append(metric)
        case_starts = [item for item in starts if int(item["case_id"]) == case]
        assert len(case_starts) == 8
        near = []
        for item in case_starts:
            for name, low, high in (("kBF", 1e-7, 10), ("beta", 1e-9, .1)):
                value = float(item[name])
                side = ("lower" if math.log(value / low) < .02 else
                        "upper" if math.log(high / value) < .02 else "")
                if side:
                    near.append(dict(start=int(item["start"]), parameter=name, side=side))
        errors = [float(item["scaled_SSE"]) for item in case_starts]
        diagnostics.append(dict(case_id=case, min_SSE=min(errors), max_SSE=max(errors),
                                exitflags=[int(item["exitflag"]) for item in case_starts],
                                near_bounds=near))
    write("model_d_sensitivity_summary.csv", summary)
    write("model_d_sensitivity_metrics.csv", metrics)
    for name in ("model_d_sensitivity_starts.csv", "model_d_sensitivity_by_dose.csv"):
        (output / name).write_bytes((source / name).read_bytes())
    status = json.loads((source / "ci_status.json").read_text())
    (output / "model_d_sensitivity_diagnostics.json").write_text(
        json.dumps(dict(run=status, starts=diagnostics), indent=2) + "\n")

    styles = {.5: ("#d87821", "-"), 1.: ("#222222", "--"), 2.: ("#2468b4", "-.")}
    parameters = ["kEF", "kFE", "kFB"]

    def plot_cases(ax, parameter, dose, column, early=False):
        for row in summary:
            if row["varied_parameter"] != parameter:
                continue
            values = curve(int(row["case_id"]), dose)
            if early:
                values = values[values[:, 0] <= 12]
            factor = row["multiplier"]
            color, style = styles[factor]
            ax.plot(values[:, 0], values[:, column], style, color=color,
                    label=f"{factor:g}x reference", lw=1.7)
        ax.set(title=f"Vary {parameter}; {dose:g} nM", xlabel="Hours")
        ax.grid(alpha=.2)

    for column, label, suffix in ((1, "Cell count", "cells"), (2, "beta CB (1/h)", "killing")):
        fig, axes = plt.subplots(3, 3, figsize=(14, 10), constrained_layout=True)
        for j, parameter in enumerate(parameters):
            for i, dose in enumerate(doses):
                ax = axes[i, j]
                plot_cases(ax, parameter, dose, column)
                if column == 1:
                    points = [row for row in observations if int(row["case_id"]) == reference_id
                              and float(row["dose_nM"]) == dose]
                    ax.scatter([float(row["time_h"]) for row in points],
                               [float(row["observed_N"]) for row in points],
                               s=15, color="black", zorder=4, label="Digitized data")
                ax.set_ylabel(label)
                if i == 0:
                    ax.legend(fontsize=8)
        fig.suptitle("Model D: one-at-a-time exploratory PK perturbations, shared kBF and beta re-fitted",
                     fontsize=13)
        fig.savefig(output / f"model_d_sensitivity_{suffix}.png", dpi=160)
        plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), constrained_layout=True)
    for ax, parameter in zip(axes, parameters):
        plot_cases(ax, parameter, 625, 2, early=True)
        ax.axvline(6, color="gray", ls=":", label="Washout")
        ax.set_ylabel("beta CB (1/h)")
        ax.legend(fontsize=8)
    fig.suptitle("Early killing-rate transient (0-12 h); no change to model mechanisms", fontsize=13)
    fig.savefig(output / "model_d_sensitivity_early_killing.png", dpi=170)
    plt.close(fig)
    print(f"Analyzed run {status['run_id']}: 9 cases, 72 starts; tables and figures in {output}")


if __name__ == "__main__":
    main()
