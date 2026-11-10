"""Experiment drivers: each function returns a tidy :class:`pandas.DataFrame`.

All studies share one protocol.  For a seed, ``size`` Haar-random stabilizer
states (null) and ``size`` random Clifford+T states with ``n_t`` T gates and
``M2 >= 0.05`` (alternative) are drawn; the same local noise is applied to both
classes; the figure of merit is the balanced error (mean of false-reject and
false-accept rates over the two ensembles).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .bell_test import (
    COPIES_PER_ROUND,
    acceptance_probability,
    balanced_error_curve,
    calibrated_error,
    copies_needed,
)
from .paulis import FloatArr, apply_local_noise
from .shadow_baseline import shadow_error
from .states import draw_ensembles


@dataclass(frozen=True)
class Noise:
    """Local noise applied to every copy: dephasing, depolarising and readout flip probabilities."""

    dephase: float = 0.0
    depol: float = 0.0
    readout: float = 0.0

    def scaled(self, factor: float) -> "Noise":
        return Noise(self.dephase * factor, self.depol * factor, self.readout * factor)


REFERENCE_NOISE = Noise(0.005, 0.005, 0.005)


def accept_probs(spectra: list[FloatArr], n: int, noise: Noise) -> np.ndarray:
    """Single-round acceptance probability of every state under ``noise``."""
    return np.array(
        [acceptance_probability(apply_local_noise(r, n, noise.dephase, noise.depol), n, noise.readout) for r in spectra]
    )


def noisy(spectra: list[FloatArr], n: int, noise: Noise) -> list[FloatArr]:
    """Apply the state-side noise (dephasing, depolarising) to every spectrum."""
    return [apply_local_noise(r, n, noise.dephase, noise.depol) for r in spectra]


def design_threshold(null: np.ndarray, alt: np.ndarray) -> tuple[float, float]:
    """Noise-agnostic design threshold and gap from the noiseless acceptance rates.

    Returns ``(tau0, gap)`` with ``tau0 = (1 + mean(alt)) / 2`` and
    ``gap = mean(null) - mean(alt)`` (both noiseless).
    """
    return float((1.0 + alt.mean()) / 2.0), float(null.mean() - alt.mean())


def asymptotic_agnostic_error(a_null: np.ndarray, a_alt: np.ndarray, tau0: float) -> float:
    """Balanced error of the noise-agnostic threshold with unlimited rounds."""
    return float(0.5 * ((a_null < tau0).mean() + (a_alt >= tau0).mean()))


def study_curves(
    n: int, n_t: int, noise: Noise, seeds: list[int], size: int, budgets: list[int], shadow_size: int, shadow_trials: int
) -> pd.DataFrame:
    """Balanced error versus copy budget for four testers at a fixed noise level."""
    rows = []
    for seed in seeds:
        rng = np.random.default_rng(1000 + seed)
        null, alt, _ = draw_ensembles(n, n_t, size, rng)
        a0, a1 = accept_probs(null, n, Noise()), accept_probs(alt, n, Noise())
        tau0, gap = design_threshold(a0, a1)
        b0, b1 = accept_probs(null, n, noise), accept_probs(alt, n, noise)
        nn, na = noisy(null[:shadow_size], n, noise), noisy(alt[:shadow_size], n, noise)
        for copies in budgets:
            rounds = max(copies // COPIES_PER_ROUND, 1)
            rows.append((seed, "Bell test, noise-aware threshold", copies, balanced_error_curve(b0, b1, rounds)[0]))
            rows.append((seed, "Bell test, noise-agnostic threshold", copies, balanced_error_curve(b0, b1, rounds, tau0)[0]))
            half = max(copies // (2 * COPIES_PER_ROUND), 1)
            rows.append((seed, "Bell test, reference-calibrated", copies, calibrated_error(b0, b1, half, gap, rng)))
            rows.append((seed, "Pauli shadows, noise-aware threshold", copies, shadow_error(nn, na, n, copies, shadow_trials, rng)))
    return pd.DataFrame(rows, columns=["seed", "method", "copies", "error"])


def study_heatmap(
    n: int, n_t: int, dephase: list[float], readout: list[float], seeds: list[int], size: int
) -> pd.DataFrame:
    """Copies to reach 10 percent error (noise-aware) and asymptotic agnostic error over a noise grid."""
    rows = []
    for seed in seeds:
        rng = np.random.default_rng(2000 + seed)
        null, alt, _ = draw_ensembles(n, n_t, size, rng)
        a0, a1 = accept_probs(null, n, Noise()), accept_probs(alt, n, Noise())
        tau0, _ = design_threshold(a0, a1)
        for pz in dephase:
            for q in readout:
                noise = Noise(dephase=pz, readout=q)
                b0, b1 = accept_probs(null, n, noise), accept_probs(alt, n, noise)
                rows.append(
                    (seed, pz, q, copies_needed(b0, b1), asymptotic_agnostic_error(b0, b1, tau0), float(b0.mean()), float(b1.mean()))
                )
    return pd.DataFrame(rows, columns=["seed", "dephase", "readout", "copies_needed", "agnostic_asymptotic_error", "null_accept", "alt_accept"])


def study_scaling(
    ns: list[int],
    shadow_ns: list[int],
    n_t: int,
    noises: dict[str, Noise],
    seeds: list[int],
    size: int,
    shadow_size: int,
    shadow_trials: int,
    shadow_grid: list[int],
) -> pd.DataFrame:
    """Copies to reach 10 percent error versus register size for the Bell test and the shadow baseline."""
    rows = []
    for seed in seeds:
        for n in ns:
            rng = np.random.default_rng(3000 + 17 * seed + n)
            null, alt, _ = draw_ensembles(n, n_t, size, rng)
            for label, noise in noises.items():
                b0, b1 = accept_probs(null, n, noise), accept_probs(alt, n, noise)
                rows.append((seed, n, label, "Bell test, noise-aware threshold", copies_needed(b0, b1)))
                if n in shadow_ns:
                    nn, na = noisy(null[:shadow_size], n, noise), noisy(alt[:shadow_size], n, noise)
                    found = float("inf")
                    for copies in shadow_grid:
                        if shadow_error(nn, na, n, copies, shadow_trials, rng) <= 0.1:
                            found = float(copies)
                            break
                    rows.append((seed, n, label, "Pauli shadows, noise-aware threshold", found))
    return pd.DataFrame(rows, columns=["seed", "n", "noise", "method", "copies_needed"])


def study_critical_noise(ns: list[int], n_t: int, seeds: list[int], size: int, grid: np.ndarray) -> pd.DataFrame:
    """Noise strength at which the noise-agnostic threshold exceeds 10 percent error, versus ``n``.

    Noise is ``Noise(p, p, p)`` (dephasing, depolarising and readout all equal to ``p``).
    """
    rows = []
    for seed in seeds:
        for n in ns:
            rng = np.random.default_rng(4000 + 17 * seed + n)
            null, alt, _ = draw_ensembles(n, n_t, size, rng)
            a0, a1 = accept_probs(null, n, Noise()), accept_probs(alt, n, Noise())
            tau0, _ = design_threshold(a0, a1)
            errs = []
            for p in grid:
                b0, b1 = accept_probs(null, n, Noise(p, p, p)), accept_probs(alt, n, Noise(p, p, p))
                errs.append(asymptotic_agnostic_error(b0, b1, tau0))
            errs_arr = np.array(errs)
            above = np.nonzero(errs_arr > 0.1)[0]
            if len(above) == 0:
                p_crit = float("nan")
            else:
                i = int(above[0])
                if i == 0:
                    p_crit = float(grid[0])
                else:
                    f = (0.1 - errs_arr[i - 1]) / (errs_arr[i] - errs_arr[i - 1])
                    p_crit = float(grid[i - 1] + f * (grid[i] - grid[i - 1]))
            rows.append((seed, n, p_crit))
    return pd.DataFrame(rows, columns=["seed", "n", "p_crit"])


def study_waterfall(n: int, n_t: int, noise: Noise, rounds: int, seeds: list[int], size: int) -> pd.DataFrame:
    """Error at a fixed budget as noise sources are added and the threshold is then repaired."""
    steps = [
        ("noiseless", Noise(), "agnostic"),
        ("+ dephasing", Noise(dephase=noise.dephase), "agnostic"),
        ("+ depolarising", Noise(noise.dephase, noise.depol), "agnostic"),
        ("+ readout flips", noise, "agnostic"),
        ("noise-aware threshold", noise, "aware"),
        ("reference-calibrated threshold", noise, "calibrated"),
    ]
    rows = []
    for seed in seeds:
        rng = np.random.default_rng(5000 + seed)
        null, alt, _ = draw_ensembles(n, n_t, size, rng)
        a0, a1 = accept_probs(null, n, Noise()), accept_probs(alt, n, Noise())
        tau0, gap = design_threshold(a0, a1)
        for name, nz, mode in steps:
            b0, b1 = accept_probs(null, n, nz), accept_probs(alt, n, nz)
            if mode == "agnostic":
                err = balanced_error_curve(b0, b1, rounds, tau0)[0]
            elif mode == "aware":
                err = balanced_error_curve(b0, b1, rounds)[0]
            else:
                err = calibrated_error(b0, b1, max(rounds // 2, 1), gap, rng)
            rows.append((seed, name, err))
    return pd.DataFrame(rows, columns=["seed", "step", "error"])
