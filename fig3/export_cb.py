"""Export the selected Fig. 3 PK input and verify it against the full MATLAB archive."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]
FULL = ROOT / "reports/fig3_184/full"
OUT = FULL.parent / "cb_input"
MODEL = "dynamic_hill"
PANELS = "abcdefghi"


def read_rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def load_parameters():
    reference = {r["parameter"]: r for r in read_rows(ROOT / "reports/model_d_fixed_parameters.csv")}
    rows = [(name, reference[name]["fixed_value"], reference[name]["unit"], "fixed")
            for name in ("kEF", "kFE", "kFB", "VI", "VE")]
    kbf, = (r["estimate"] for r in read_rows(FULL / "parameters.csv")
            if r["model"] == MODEL and r["parameter"] == "kBF")
    rows.append(("kBF", kbf, "h^-1", "shared_fit"))
    return rows


def bound_concentration(times, dose, wash, p):
    """Bound drug (nM) at nonnegative times (h), with zero initial CF and CB."""
    times = np.asarray(times, dtype=float)
    if np.any(~np.isfinite(times)) or np.any(times < 0):
        raise ValueError("Times must be finite and nonnegative")
    # During dosing the third (constant) state supplies the extracellular input.
    on = np.array([[-p["kFE"] - p["kFB"], p["kBF"], p["kEF"] * p["VE"] / p["VI"] * dose],
                   [p["kFB"], -p["kBF"], 0.0],
                   [0.0, 0.0, 0.0]])
    before = times <= wash
    cb = np.empty_like(times)
    cb[before] = expm(times[before, None, None] * on)[:, 1, 2]
    at_wash = expm(wash * on)[:2, 2]
    after = expm((times[~before, None, None] - wash) * on[:2, :2])
    cb[~before] = after[:, 1, :] @ at_wash
    return cb


def draw(curves):
    fig, axes = plt.subplots(3, 3, figsize=(12, 8), sharex=True, sharey=True)
    fig.subplots_adjust(left=.09, right=.98, top=.84, bottom=.10, wspace=.21, hspace=.48)
    fig.suptitle("Bound intracellular doxorubicin", x=.09, y=.965, ha="left",
                 fontsize=19, weight="bold", color="#20394b")
    blue = "#12658a"
    legend = [Line2D([], [], color=blue, label="Within cell observation window"),
              Line2D([], [], color=blue, ls="--", label="Beyond last cell observation"),
              Line2D([], [], lw=1.2, color="#71573d", ls=":", label="Washout")]
    fig.legend(handles=legend, loc="upper left", bbox_to_anchor=(.082,.925), ncol=3, frameon=False)
    for ax, (panel, dose, wash, end, t, cb) in zip(axes.flat, curves):
        ax.plot(t[t <= end], cb[t <= end], lw=1.9, color=blue)
        ax.plot(t[t >= end], cb[t >= end], lw=1.9, ls="--", color=blue)
        ax.axvline(wash, color="#71573d", lw=1.2, ls=":")
        ax.set_title(f"{panel}    {dose:.8g} nM  ·  {wash:g} h", loc="left", fontsize=11)
        ax.grid(axis="y", alpha=.22)
        ax.spines[["right", "top"]].set_visible(False)
        ax.tick_params(labelbottom=True, labelleft=True, labelsize=9)
    axes[0, 0].set(xlim=(0, 750), xticks=np.arange(0, 751, 150))
    axes[0, 0].set_yscale("symlog", linthresh=1)
    axes[0, 0].set(ylim=(0, 1000), yticks=[0,1,10,100,1000],
                   yticklabels=["0","1","10","100","1000"])
    fig.supxlabel("Time since dosing (h)", y=.035)
    fig.supylabel("Bound concentration (nM; linear to 1, log above)", x=.015)
    fig.savefig(OUT / "CB.png", dpi=175)
    plt.close(fig)


def main():
    params = load_parameters()
    p = {name: float(value) for name, value, _, _ in params}
    observations = read_rows(ROOT / "fig3/observations.csv")
    model_id, = (r["model_id"] for r in read_rows(FULL / "summary.csv") if r["model"] == MODEL)
    archive = [r for r in read_rows(FULL / "curves.csv") if r["model_id"] == model_id]
    curves, output = [], []
    max_error = 0.0
    for panel_id, panel in enumerate(PANELS, 1):
        obs = [r for r in observations if r["panel"] == panel]
        dose, = {float(r["dose_nM"]) for r in obs}
        wash, = {float(r["wash_h"]) for r in obs}
        observed = {float(r["time_h"]) for r in obs}
        end = max(observed)
        times = np.unique([*range(751), wash, *observed])
        cb = bound_concentration(times, dose, wash, p)
        saved = [r for r in archive if int(r["panel_id"]) == panel_id]
        check_times = np.array([float(r["time_h"]) for r in saved])
        reference = np.array([float(r["CB_nM"]) for r in saved])
        assert observed.issubset(check_times) and end == max(check_times)
        replay = bound_concentration(check_times, dose, wash, p)
        max_error = max(max_error, np.max(abs(replay - reference)))
        np.testing.assert_allclose(replay, reference, rtol=1e-10, atol=1e-10)
        assert cb[0] == 0 and np.all(np.isfinite(cb)) and np.all(cb >= 0)
        for t, concentration in zip(times, cb):
            output.append((panel, dose, wash, t, dose if t < wash else 0.,
                           concentration, int(t in observed), int(t > end)))
        curves.append((panel, dose, wash, end, times, cb))

    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "pk_parameters.csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["parameter", "value", "unit", "status"])
        writer.writerows(params)
    with (OUT / "CB.csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["panel", "dose_nM", "wash_h", "time_h", "CE_nM", "CB_nM",
                         "is_cell_observation_time", "after_last_cell_observation"])
        writer.writerows(output)
    draw(curves)
    print(f"{len(output)} rows, {len(archive)} MATLAB comparisons, max |ΔCB| {max_error:.3g} nM")


if __name__ == "__main__":
    main()
