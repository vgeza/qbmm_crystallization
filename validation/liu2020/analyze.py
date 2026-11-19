#!/usr/bin/env python3
"""Compare the Liu et al. (2020) cases with the experiments and the paper's model.
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent

OPENFOAM, EXPERIMENT, PAPER, CSAT = "#2a78d6", "#eb6834", "#1baf7a", "#52514e"


def read_batch(case):
    """Solver output as a dict of columns"""
    f = HERE / case / "postProcessing" / "batch" / "batch.csv"
    if not f.exists():
        return None
    data = np.genfromtxt(f, delimiter=",", names=True)
    out = {name: data[name] for name in data.dtype.names}
    out["t_min"] = out["t"] / 60.0 # time to minutes
    return out


def read_reference(case):
    """Digitized curves: (t [min], C [mg/mL]) of the experiment and the model."""
    rows = np.genfromtxt(HERE / "refdata" / f"{case}.csv", delimiter=",", skip_header=2)
    exp = rows[:, 0:2][~np.isnan(rows[:, 0])]
    model = rows[:, 2:4][~np.isnan(rows[:, 2])]
    return exp[np.argsort(exp[:, 0])], model[np.argsort(model[:, 0])]


def rmse_against(exp, t, c):
    """RMSE of the curve (t, c) at the experiment's times (inside t's range)."""
    keep = (exp[:, 0] >= t[0]) & (exp[:, 0] <= t[-1])
    return float(np.sqrt(np.mean((np.interp(exp[keep, 0], t, c) - exp[keep, 1]) ** 2)))


def main():
    cases = sys.argv[1:] or sorted(p.name for p in HERE.glob("val*") if p.is_dir())
    runs = {c: read_batch(c) for c in cases}
    missing = [c for c, r in runs.items() if r is None]
    if missing:
        print("not run yet (no batch.csv):", ", ".join(missing))
    runs = {c: r for c, r in runs.items() if r is not None}
    if not runs:
        return

    fig, axes = plt.subplots(len(runs), 2, figsize=(10, 3.2 * len(runs)), squeeze=False)
    table = []

    for row, (case, r) in zip(axes, runs.items()):
        exp, model = read_reference(case)

        ax = row[0]
        ax.plot(exp[:, 0], exp[:, 1], "o", ms=4, color=EXPERIMENT, label="LIU 2020 experiment")
        ax.plot(model[:, 0], model[:, 1], "--", lw=2, color=PAPER, label="LIU 2020 model")
        ax.plot(r["t_min"], r["C"], "-", lw=2, color=OPENFOAM, label="OpenFOAM QBMM")
        ax.plot(r["t_min"], r["Csat"], ":", lw=1.5, color=CSAT, label="solubility Csat")
        ax.set_title(f"{case}: solute concentration", loc="left", fontsize=10)
        ax.set_xlabel("time [min]")
        ax.set_ylabel("C [mg/mL]")

        ax = row[1]
        ax.plot(r["t_min"], r["d43"], "-", lw=2, color=OPENFOAM, label="d43 = m4/m3")
        ax.plot(r["t_min"], r["d10"], "--", lw=2, color=OPENFOAM, label="d10 = m1/m0")
        ax.set_title(f"{case}: mean crystal size (OpenFOAM)", loc="left", fontsize=10)
        ax.set_xlabel("time [min]")
        ax.set_ylabel("size [um]")

        for ax in row:
            ax.grid(alpha=0.3, lw=0.5)
            ax.spines[["top", "right"]].set_visible(False)
            ax.legend(fontsize=8, frameon=False)

        table.append(dict(
            case=case,
            rmse_openfoam=rmse_against(exp, r["t_min"], r["C"]),
            rmse_paper=rmse_against(exp, model[:, 0], model[:, 1]),
            C_end=r["C"][-1],
            m0_end=r["m0"][-1],
            d10_end=r["d10"][-1],
            d43_end=r["d43"][-1],
        ))

    out = HERE / "results"
    out.mkdir(exist_ok=True)
    fig.tight_layout()
    fig.savefig(out / "liu2020.png", dpi=130)

    cols = ["case", "rmse_openfoam", "rmse_paper", "C_end", "m0_end", "d10_end", "d43_end"]
    lines = [",".join(cols)] + [
        ",".join(row["case"] if c == "case" else "%.6g" % row[c] for c in cols) for row in table
    ]
    (out / "summary.csv").write_text("\n".join(lines) + "\n")

if __name__ == "__main__":
    main()
