"""Reproduce stage-one tables and figures from the analytic clock models."""
from pathlib import Path
import csv
import json
import math
import platform
import hashlib
import numpy as np
import scipy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from clock_model import optimal_clock, exponential_baseline, SmoothClock, CompactClock

ROOT = Path(__file__).resolve().parent
if ROOT.name == "code":
    ROOT = ROOT.parent
OUT = ROOT / "results"


def main():
    OUT.mkdir(exist_ok=True)
    rows = []
    for Q in (0.5, 1, 2, 5, 10, 20):
        for delta in (0.5, 0.1, 0.01, 0.001, 0.0001):
            k = 1 - delta
            clock = optimal_clock(Q, k)
            baseline = exponential_baseline(Q, k)
            rows.append({"Q": Q, "requested_kernel_error": delta,
                         "optimal_mean": clock.mean, "optimal_KL": clock.kl,
                         "exponential_mean": baseline["mean"],
                         "exponential_KL": baseline["kl_nats"],
                         "mean_ratio": baseline["mean"] / clock.mean})
    with (OUT / "cost-comparison.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    frontier = []
    Q, k = 5.0, 0.99
    unconstrained = optimal_clock(Q, k)
    caps = np.r_[k * Q + np.geomspace(0.001, unconstrained.mean - k*Q, 150), 8, 10]
    for B in caps:
        model = optimal_clock(Q, k, float(B))
        frontier.append({"mean_cap": float(B), **model.record()})
    (OUT / "energy-frontier.json").write_text(json.dumps(frontier, indent=2) + "\n", encoding="utf-8")

    delta = np.geomspace(1e-5, 0.5, 200)
    exact = [optimal_clock(Q, 1-d) for d in delta]
    exp = [exponential_baseline(Q, 1-d) for d in delta]
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.6), constrained_layout=True)
    axes[0].loglog(delta, [e["mean"] for e in exp], label="Unshifted exponential", color="#ba4935")
    axes[0].loglog(delta, [c.mean for c in exact], label="Minimum-KL clock", color="#126d82")
    axes[0].axhline(10, color="#6a52a3", linestyle="--", label="Compact clock, exact window")
    axes[0].set(xlabel="Allowed worst kernel error", ylabel="Mean modular clock energy", title="Cost depends strongly on state family")
    axes[0].legend(fontsize=8)
    axes[1].plot([r["mean"] for r in frontier], [r["kl_nats"] for r in frontier], color="#126d82")
    axes[1].axvline(k*Q, linestyle=":", color="#ba4935", label="Infeasible at / below kQ")
    axes[1].set(xlabel="Mean modular clock energy", ylabel="Relative entropy (nats)", title="Entropy-energy tradeoff, Q=5, k=0.99")
    axes[1].legend(fontsize=8)
    q_grid = np.linspace(0, 18, 600)
    axes[2].plot(q_grid, np.exp(-q_grid/10), label="Exponential, mean=10", color="#ba4935")
    axes[2].plot(q_grid, SmoothClock(5, 1).kernel(q_grid), label="Shifted Gamma, mean=10", color="#126d82")
    axes[2].plot(q_grid, CompactClock(5, 10).kernel(q_grid), label="Compact smooth, mean=10", color="#6a52a3")
    axes[2].axvline(Q, color="#555555", linestyle=":")
    axes[2].set(xlabel="Intermediate modular-energy threshold q", ylabel="K(q)", title="Equal means do not imply equal recovery")
    axes[2].legend(fontsize=8)
    fig.suptitle("Finite-clock spectral recovery in a fixed weak-gravity algebra\nKernel accuracy is not a full observer or finite-G accuracy claim", fontsize=12)
    fig.savefig(OUT / "clock-costs.png", dpi=180)
    fig.savefig(OUT / "clock-costs.svg")
    plt.close(fig)

    smooth_rows = []
    for W in (1, 2, 5, 10, 20, 50):
        smooth_rows.append(CompactClock(5, W).record())
    (OUT / "compact-clocks.json").write_text(json.dumps(smooth_rows, indent=2)+"\n", encoding="utf-8")
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), constrained_layout=True)
    p = np.linspace(0, 26, 1000)
    t = np.linspace(-2, 2, 2000)
    for W, color in ((5, "#ba4935"), (10, "#126d82"), (20, "#6a52a3")):
        model = CompactClock(5, W)
        axes[0].plot(p, model.density(p), label=f"W={W}", color=color)
        axes[1].plot(t, model.time_density(t), label=f"W={W}, RMS={math.sqrt(model.time_variance):.3f}", color=color)
    axes[0].set(xlabel="Modular energy p", ylabel="P(p)", title="Compact energy support, exact K(q<=5)=1")
    axes[1].set(xlabel="Conjugate dimensionless time", ylabel="Fourier POVM density", title="Wider energy support sharpens this time POVM")
    for ax in axes:
        ax.legend(fontsize=8)
    fig.savefig(OUT / "compact-clock-timing.png", dpi=180)
    fig.savefig(OUT / "compact-clock-timing.svg")
    plt.close(fig)

    sample = {"Q": Q, "k": k,
              "exponential": exponential_baseline(Q, k),
              "entropy_optimal": unconstrained.record(),
              "mean_cap_5_2": optimal_clock(Q, k, 5.2).record(),
              "gamma": SmoothClock(5, 1).record(),
              "compact": CompactClock(5, 10).record(),
              "scope": "Fixed vacuum spectral kernels; no finite-G simulation; timing uses a Fourier POVM."}
    code_dir = Path(__file__).resolve().parent
    manifest = {"python": platform.python_version(), "numpy": np.__version__,
                "scipy": scipy.__version__, "matplotlib": matplotlib.__version__,
                "code_sha256": {name: hashlib.sha256((code_dir/name).read_bytes()).hexdigest()
                                for name in ("clock_model.py", "run_study.py")},
                "sample": sample, "comparison_cases": len(rows)}
    (OUT / "study.json").write_text(json.dumps(manifest, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(sample, indent=2))
    print(f"Wrote {len(rows)} comparison cases and {len(frontier)} frontier points to {OUT}")

if __name__ == "__main__":
    main()
