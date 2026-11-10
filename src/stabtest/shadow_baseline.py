"""Single-copy classical-shadow baseline for stabilizer testing.

Each copy is measured once in a uniformly random product Pauli basis.  The
Pauli spectrum of the state is estimated from these outcomes (inverse channel
``3**|S|``) and summarised by the statistic ``T = 2**-n sum_a r_a**4``, which
equals one for pure stabilizer states and is smaller for magic states.  Four
independent batches give an unbiased estimate of every fourth power.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np

from .paulis import FloatArr, popcount, wht


@lru_cache(maxsize=None)
def _basis_tables(n: int) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(idx, weight)`` where ``idx[b, S]`` is the Pauli index obtained by
    keeping the factors of basis ``b`` on the subset ``S`` and ``weight`` is ``3**|S|``."""
    dim = 2**n
    nb = 3**n
    digits = np.array([[(b // 3**i) % 3 for i in range(n)] for b in range(nb)])  # 0=X 1=Y 2=Z
    has_x = (digits <= 1).astype(np.int64)
    has_z = (digits >= 1).astype(np.int64)
    bits = ((np.arange(dim)[:, None] >> np.arange(n)[None, :]) & 1).astype(np.int64)  # (S, i)
    ax = (has_x[:, None, :] * bits[None, :, :]) @ (1 << np.arange(n))
    az = (has_z[:, None, :] * bits[None, :, :]) @ (1 << np.arange(n))
    idx = ax * dim + az
    weight = 3.0 ** popcount(np.arange(dim))
    return idx, weight


def _outcome_probabilities(r: FloatArr, n: int) -> np.ndarray:
    idx, _ = _basis_tables(n)
    p = wht(r[idx]) / 2**n
    p = np.clip(p, 0.0, None)
    return p / p.sum(axis=1, keepdims=True)


def _batch_estimate(probs: np.ndarray, n: int, shots: int, rng: np.random.Generator) -> FloatArr:
    idx, weight = _basis_tables(n)
    nb = probs.shape[0]
    per_basis = rng.multinomial(shots, np.full(nb, 1.0 / nb))
    counts = np.stack([rng.multinomial(c, probs[b]) for b, c in enumerate(per_basis)])
    w = wht(counts.astype(float))  # w[b, S] = sum_s counts (-1)^{s.S}
    est = np.bincount(idx.ravel(), weights=(w * weight[None, :]).ravel(), minlength=4**n)
    return est / shots


def shadow_statistic(r: FloatArr, n: int, copies: int, rng: np.random.Generator) -> float:
    """Draw one shadow experiment of ``copies`` single-copy shots and return ``T``."""
    probs = _outcome_probabilities(r, n)
    shots = max(copies // 4, 1)
    prod = np.ones(4**n)
    for _ in range(4):
        prod *= _batch_estimate(probs, n, shots, rng)
    return float(prod.sum() / 2**n)


def shadow_error(
    null: list[FloatArr],
    alt: list[FloatArr],
    n: int,
    copies: int,
    trials: int,
    rng: np.random.Generator,
) -> float:
    """Balanced error of the oracle (best fixed) threshold on ``T``, pooled over instances."""
    t_null = np.array([shadow_statistic(r, n, copies, rng) for r in null for _ in range(trials)])
    t_alt = np.array([shadow_statistic(r, n, copies, rng) for r in alt for _ in range(trials)])
    cand = np.unique(np.concatenate([t_null, t_alt]))
    # stabilizer is declared when T >= c
    fr = (t_null[None, :] < cand[:, None]).mean(axis=1)
    fa = (t_alt[None, :] >= cand[:, None]).mean(axis=1)
    return float(np.min(0.5 * (fr + fa)))
