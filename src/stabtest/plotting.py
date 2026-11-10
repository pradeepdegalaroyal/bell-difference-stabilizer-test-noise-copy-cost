"""Figure helpers (viridis palette, 200 dpi)."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

DPI = 200
VIR = plt.get_cmap("viridis")
STYLES = [("-", "o"), ("--", "s"), (":", "^"), ("-.", "D")]


def _save(fig: plt.Figure, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path, dpi=DPI)
    plt.close(fig)


def line_band(df: pd.DataFrame, path: Path, title: str, noise_label: str) -> None:
    """Balanced error versus copies with mean +/- std bands over seeds."""
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    methods = list(dict.fromkeys(df["method"]))
    for i, m in enumerate(methods):
        g = df[df["method"] == m].groupby("copies")["error"].agg(["mean", "std"]).reset_index()
        c = VIR(0.05 + 0.85 * i / max(len(methods) - 1, 1))
        ls, mk = STYLES[i % len(STYLES)]
        ax.plot(g["copies"], g["mean"], ls=ls, marker=mk, color=c, label=m, ms=5)
        ax.fill_between(g["copies"], g["mean"] - g["std"].fillna(0), g["mean"] + g["std"].fillna(0), color=c, alpha=0.2)
    ax.axhline(0.1, color="grey", lw=0.8, ls=":")
    ax.set_xscale("log")
    ax.set_xlabel("State copies consumed (count)")
    ax.set_ylabel("Balanced error (fraction)")
    ax.set_title(f"{title}\n{noise_label}", fontsize=10)
    ax.set_ylim(-0.02, 0.55)
    ax.legend(fontsize=8)
    _save(fig, path)


def heatmap_pair(df: pd.DataFrame, path: Path, n: int) -> None:
    """Two annotated heatmaps over dephasing and readout error."""
    g = df.groupby(["dephase", "readout"]).agg(
        copies=("copies_needed", "mean"), err=("agnostic_asymptotic_error", "mean")
    ).reset_index()
    xs, ys = sorted(g["readout"].unique()), sorted(g["dephase"].unique())
    cop = g.pivot(index="dephase", columns="readout", values="copies").loc[ys, xs].to_numpy()
    err = g.pivot(index="dephase", columns="readout", values="err").loc[ys, xs].to_numpy()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    for ax, data, ttl, fmt, cb in [
        (axes[0], cop, "Copies to reach 10% error\n(noise-aware threshold)", "{:.0f}", "copies (count)"),
        (axes[1], err, "Asymptotic error, noise-agnostic threshold\n(unlimited copies)", "{:.2f}", "balanced error (fraction)"),
    ]:
        im = ax.imshow(data, origin="lower", cmap="viridis", aspect="auto")
        ax.set_xticks(range(len(xs)), [f"{v:g}" for v in xs])
        ax.set_yticks(range(len(ys)), [f"{v:g}" for v in ys])
        ax.set_xlabel("Readout flip probability per bit")
        ax.set_ylabel("Dephasing probability per qubit")
        ax.set_title(ttl, fontsize=10)
        vmax = np.nanmax(data[np.isfinite(data)])
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                v = data[i, j]
                ax.text(j, i, "inf" if not np.isfinite(v) else fmt.format(v), ha="center", va="center", fontsize=8,
                        color="white" if v < 0.55 * vmax else "black")
        fig.colorbar(im, ax=ax, label=cb)
    fig.suptitle(f"Noise tolerance of the Bell-difference-sampling stabilizer test (n = {n} qubits)", fontsize=11)
    _save(fig, path)


def waterfall(df: pd.DataFrame, path: Path, rounds: int, n: int) -> None:
    """Waterfall of the error contribution of each noise source and threshold repair."""
    order = list(dict.fromkeys(df["step"]))
    piv = df.pivot(index="seed", columns="step", values="error")[order]
    delta = piv.diff(axis=1)
    delta[order[0]] = piv[order[0]]
    mean, std = delta.mean(), delta.std(ddof=1)
    total_m, total_s = piv[order[-1]].mean(), piv[order[-1]].std(ddof=1)
    labels = order + ["final error"]
    fig, ax = plt.subplots(figsize=(8.4, 4.6))
    run = 0.0
    for i, step in enumerate(order):
        d = mean[step]
        bottom = 0.0 if i == 0 else (run if d >= 0 else run + d)
        ax.bar(i, abs(d), bottom=bottom, yerr=std[step], color=VIR(0.1 + 0.8 * i / len(order)), edgecolor="black", capsize=3)
        run += d
    ax.bar(len(order), total_m, yerr=total_s, color=VIR(0.95), edgecolor="black", capsize=3)
    ax.set_xticks(range(len(labels)), labels, rotation=25, ha="right", fontsize=8)
    ax.set_ylabel("Balanced error (fraction)")
    ax.set_title(f"Error contribution of each ablated component (n = {n}, {rounds} rounds = {6 * rounds} copies)", fontsize=10)
    _save(fig, path)


def scaling(sc: pd.DataFrame, pc: pd.DataFrame, path: Path) -> None:
    """Copies needed versus n (left) and critical noise of the agnostic threshold (right)."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    ax = axes[0]
    combos = [(m, z) for m in dict.fromkeys(sc["method"]) for z in dict.fromkeys(sc["noise"])]
    for i, (m, z) in enumerate(combos):
        g = sc[(sc["method"] == m) & (sc["noise"] == z)]
        g = g.replace(np.inf, np.nan).groupby("n")["copies_needed"].agg(["mean", "std"]).dropna(subset=["mean"]).reset_index()
        if g.empty:
            continue
        c = VIR(0.05 + 0.85 * i / max(len(combos) - 1, 1))
        ls, mk = STYLES[i % len(STYLES)]
        ax.plot(g["n"], g["mean"], ls=ls, marker=mk, color=c, label=f"{m}; {z}", ms=5)
        ax.fill_between(g["n"], np.maximum(g["mean"] - g["std"].fillna(0), 1), g["mean"] + g["std"].fillna(0), color=c, alpha=0.2)
    ax.set_yscale("log")
    ax.set_xticks(sorted(sc["n"].unique()))
    ax.set_xlabel("Register size n (qubits)")
    ax.set_ylabel("Copies to reach 10% error (count)")
    ax.set_title("Copy cost versus register size", fontsize=10)
    ax.legend(fontsize=7)
    ax = axes[1]
    g = pc.groupby("n")["p_crit"].agg(["mean", "std"]).reset_index()
    ax.plot(g["n"], g["mean"], "o-", color=VIR(0.3))
    ax.fill_between(g["n"], g["mean"] - g["std"].fillna(0), g["mean"] + g["std"].fillna(0), color=VIR(0.3), alpha=0.25)
    ax.set_xticks(sorted(pc["n"].unique()))
    ax.set_xlabel("Register size n (qubits)")
    ax.set_ylabel("Critical noise strength p (probability)")
    ax.set_title("Noise at which a fixed design threshold fails (>10% error)", fontsize=10)
    _save(fig, path)
