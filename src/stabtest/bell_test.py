"""Exact statistics of the Bell-difference-sampling stabilizer test.

One round of the test consumes six copies of the (noisy) state ``rho``:

1. Bell-measure copies (1,2) and (3,4); the outcomes ``x1, x2`` are ``2n``-bit
   strings (one X-type and one Z-type bit per qubit).  Each bit flips with the
   readout error ``q``.
2. Set ``x = x1 xor x2`` and measure the Pauli ``P_x`` on each of copies 5 and 6.
3. Accept the round when the two outcomes agree (the product is ``+1``).

For ``rho = 2**-n sum_a r_a P_a`` the Bell outcome law is
``B(x) = 4**-n sum_a s_a (-1)**<x,a> r_a**2`` with ``s_a`` the transpose sign,
and a round accepts with probability ``(1 + g(x) r_x**2) / 2`` where
``g(x) = (1 - 2q)**(2 |x|)`` accounts for readout flips of the ``2|x|`` bits
that enter the parity.  Stabilizer states are accepted with probability one in
the noiseless case.
"""

from __future__ import annotations

import numpy as np
from scipy.stats import binom

from .paulis import FloatArr, pauli_index_tables, wht

COPIES_PER_ROUND = 6


def bell_distribution(r: FloatArr, n: int, readout: float = 0.0) -> FloatArr:
    """Law of one Bell measurement on ``rho (x) rho`` indexed ``xx * 2**n + xz``."""
    dim = 2**n
    sign, _, _ = pauli_index_tables(n)
    f = wht(sign * r**2)  # indexed (ux, uz) with ux = xz, uz = xx
    b = (f.reshape(dim, dim).T / 4.0**n).reshape(-1)
    if readout > 0.0:
        t = b.reshape((2,) * (2 * n))
        for axis in range(2 * n):
            t = (1.0 - readout) * t + readout * np.flip(t, axis=axis)
        b = t.reshape(-1)
    return b


def difference_distribution(b: FloatArr, n: int) -> FloatArr:
    """XOR-convolution of two independent Bell outcomes (Bell difference sampling)."""
    d = wht(wht(b) ** 2) / 4.0**n
    return np.clip(d, 0.0, None)


def acceptance_probability(r: FloatArr, n: int, readout: float = 0.0) -> float:
    """Exact probability that one six-copy round accepts."""
    _, _, sw = pauli_index_tables(n)
    d = difference_distribution(bell_distribution(r, n, readout), n)
    g = (1.0 - 2.0 * readout) ** (2 * sw)
    return float(np.sum(d * 0.5 * (1.0 + g * r**2)))


def error_vs_count_threshold(acc_null: np.ndarray, acc_alt: np.ndarray, rounds: int) -> np.ndarray:
    """Balanced error for every count threshold ``k = 0 .. rounds + 1``, exact in the binomial law.

    ``acc_null`` / ``acc_alt`` are per-instance single-round acceptance
    probabilities.  The test declares "stabilizer" when at least ``k`` of the
    ``rounds`` rounds accept; the error averages the false-reject probability of
    the null ensemble and the false-accept probability of the alternative.
    """
    ks = np.arange(0, rounds + 2)
    fr = binom.cdf(ks[None, :] - 1, rounds, acc_null[:, None]).mean(axis=0)
    fa = binom.sf(ks[None, :] - 1, rounds, acc_alt[:, None]).mean(axis=0)
    return 0.5 * (fr + fa)


def balanced_error_curve(
    acc_null: np.ndarray, acc_alt: np.ndarray, rounds: int, threshold: float | None = None
) -> tuple[float, float]:
    """Balanced error of the best (oracle) or a fixed acceptance-rate threshold.

    With ``threshold=None`` the count threshold minimising the ensemble error is
    used (an oracle, state-independent threshold); otherwise
    ``k = ceil(threshold * rounds)``.  Returns ``(error, threshold_rate)``.
    """
    err = error_vs_count_threshold(acc_null, acc_alt, rounds)
    if threshold is None:
        k = int(np.argmin(err))
        return float(err[k]), k / rounds
    k = int(np.clip(np.ceil(threshold * rounds - 1e-12), 0, rounds + 1))
    return float(err[k]), threshold


def copies_needed(
    acc_null: np.ndarray, acc_alt: np.ndarray, target: float = 0.1, max_rounds: int = 2**22
) -> float:
    """Smallest number of copies (6 per round) reaching balanced error ``<= target``.

    Returns ``inf`` when the ensembles overlap so strongly that no budget up to
    ``max_rounds`` rounds reaches the target.
    """
    grid = np.unique(np.round(np.geomspace(1, max_rounds, 60)).astype(int))
    for rounds in grid:
        if balanced_error_curve(acc_null, acc_alt, int(rounds))[0] <= target:
            lo, hi = (int(grid[max(0, np.searchsorted(grid, rounds) - 1)]), int(rounds))
            while hi - lo > max(1, lo // 50):  # refine to ~2 percent
                mid = (lo + hi) // 2
                if balanced_error_curve(acc_null, acc_alt, mid)[0] <= target:
                    hi = mid
                else:
                    lo = mid
            return float(hi * COPIES_PER_ROUND)
    return float("inf")


def calibrated_error(
    acc_null: np.ndarray,
    acc_alt: np.ndarray,
    rounds: int,
    gap: float,
    rng: np.random.Generator,
    n_draws: int = 400,
) -> float:
    """Balanced error when the threshold is calibrated on a noisy reference state.

    A fresh stabilizer reference (drawn from the null ensemble) is run for the
    same number of ``rounds`` as the unknown state; the threshold is placed at
    ``a_ref - gap / 2`` where ``gap`` is the *noiseless* null-minus-alt
    acceptance gap.  The error is averaged over the reference's shot noise.
    """
    ref = rng.choice(acc_null, size=n_draws)
    a_hat = rng.binomial(rounds, ref) / rounds
    err = error_vs_count_threshold(acc_null, acc_alt, rounds)
    ks = np.clip(np.ceil((a_hat - gap / 2) * rounds - 1e-12), 0, rounds + 1).astype(int)
    return float(np.mean(err[ks]))
