"""Reproduce every result and figure.  ``python experiments/run_all.py [--quick]``."""

from __future__ import annotations

import argparse
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from stabtest import plotting, studies  # noqa: E402
from stabtest.summary import summarise  # noqa: E402
from stabtest.studies import REFERENCE_NOISE, Noise  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="smoke run (< 1 minute), writes to results_quick/")
    quick = ap.parse_args().quick
    res = ROOT / ("results_quick" if quick else "results")
    fig = ROOT / ("figures_quick" if quick else "figures")
    res.mkdir(exist_ok=True)
    fig.mkdir(exist_ok=True)
    t0 = time.time()

    n, n_t = 4, 1
    seeds = [0, 1] if quick else [0, 1, 2, 3, 4]
    size = 6 if quick else 24
    sh_size, sh_trials = (4, 3) if quick else (12, 6)
    budgets = [60, 300, 3000] if quick else [60, 120, 240, 480, 960, 1920, 3840, 7680, 15360]
    ref = REFERENCE_NOISE

    curves = studies.study_curves(n, n_t, ref, seeds, size, budgets, sh_size, sh_trials)
    curves.to_csv(res / "curves.csv", index=False)
    plotting.line_band(
        curves, fig / "fig1_error_vs_copies.png", f"Stabilizer test error versus copies (n = {n}, {n_t} T gate)",
        f"dephasing = depolarising = readout flip = {ref.dephase:g}",
    )
    print(f"[{time.time() - t0:5.1f}s] curves", flush=True)

    deph = [0.0, 0.01] if quick else [0.0, 0.0025, 0.005, 0.01, 0.02]
    ro = [0.0, 0.01] if quick else [0.0, 0.0025, 0.005, 0.01, 0.02]
    heat = studies.study_heatmap(n, n_t, deph, ro, seeds, size)
    heat.to_csv(res / "heatmap.csv", index=False)
    plotting.heatmap_pair(heat, fig / "fig2_noise_heatmap.png", n)
    print(f"[{time.time() - t0:5.1f}s] heatmap", flush=True)

    ns = [3, 4] if quick else [3, 4, 5, 6]
    shadow_ns = [3] if quick else [3, 4, 5]
    grid = [200, 1000, 5000] if quick else [250, 500, 1000, 2000, 4000, 8000, 16000, 32000, 64000, 128000]
    noises = {"noiseless": Noise(), "reference noise": ref}
    sc = studies.study_scaling(ns, shadow_ns, n_t, noises, seeds, size, sh_size, sh_trials, grid)
    sc.to_csv(res / "scaling.csv", index=False)
    pgrid = np.linspace(0.0005, 0.03, 12 if quick else 60)
    pc = studies.study_critical_noise(ns, n_t, seeds, size, pgrid)
    pc.to_csv(res / "critical_noise.csv", index=False)
    plotting.scaling(sc, pc, fig / "fig3_scaling.png")
    print(f"[{time.time() - t0:5.1f}s] scaling", flush=True)

    rounds = 100
    wf = studies.study_waterfall(n, n_t, ref, rounds, seeds, size)
    wf.to_csv(res / "waterfall.csv", index=False)
    plotting.waterfall(wf, fig / "fig4_waterfall.png", rounds, n)

    info = {
        "quick": quick, "seeds": seeds, "ensemble_size_per_class": size, "n": n, "n_t": n_t,
        "reference_noise": vars(ref), "runtime_s": round(time.time() - t0, 1),
        "python": platform.python_version(), "platform": platform.platform(),
    }
    (res / "run_info.json").write_text(json.dumps(info, indent=2))
    summarise(res)
    print(f"[{time.time() - t0:5.1f}s] done -> {res}", flush=True)


if __name__ == "__main__":
    main()
