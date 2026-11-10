"""Pauli-spectrum utilities for small qubit registers.

A Pauli operator on ``n`` qubits is indexed by a pair of ``n``-bit integers
``(ax, az)`` and stored at flat position ``ax * 2**n + az``.  The operator is
``i**popcount(ax & az) * X**ax * Z**az`` (so ``Y = i X Z`` and every Pauli is
Hermitian).  Bit ``i`` of an integer refers to qubit ``i``, which matches the
little-endian convention of Qiskit statevectors.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

FloatArr = NDArray[np.float64]


def popcount(x: NDArray[np.int64] | int) -> NDArray[np.int64]:
    """Number of set bits of each integer in ``x`` (vectorised)."""
    x = np.asarray(x, dtype=np.int64)
    out = np.zeros_like(x)
    while np.any(x):
        out += x & 1
        x = x >> 1
    return out


def wht(v: NDArray) -> NDArray:
    """Fast Walsh-Hadamard transform (unnormalised) along the last axis.

    ``wht(v)[u] = sum_j v[j] * (-1)**popcount(u & j)``.  The last axis length
    must be a power of two.
    """
    v = np.array(v, copy=True)
    n = v.shape[-1]
    h = 1
    while h < n:
        v = v.reshape(v.shape[:-1] + (n // (2 * h), 2, h))
        a, b = v[..., 0, :], v[..., 1, :]
        v = np.stack((a + b, a - b), axis=-2)
        v = v.reshape(v.shape[:-3] + (n,))
        h *= 2
    return v


def pauli_index_tables(n: int) -> tuple[FloatArr, NDArray[np.int64], NDArray[np.int64]]:
    """Return ``(transpose_sign, x_weight, support_weight)`` for all 4**n Paulis.

    ``transpose_sign`` is ``(-1)**(number of Y factors)``, the sign picked up by
    the transpose.  ``x_weight`` counts X or Y factors (dephasing-sensitive) and
    ``support_weight`` counts non-identity factors.
    """
    dim = 2**n
    ax = np.repeat(np.arange(dim), dim)
    az = np.tile(np.arange(dim), dim)
    sign = 1.0 - 2.0 * (popcount(ax & az) % 2)
    return sign, popcount(ax), popcount(ax | az)


def pauli_spectrum(psi: NDArray[np.complex128]) -> FloatArr:
    """All ``4**n`` Pauli expectation values ``<psi|P|psi>`` of a pure state."""
    psi = np.asarray(psi, dtype=np.complex128)
    dim = psi.shape[0]
    n = int(np.log2(dim))
    j = np.arange(dim)
    out = np.empty((dim, dim), dtype=np.complex128)
    for ax in range(dim):
        out[ax] = wht(np.conj(psi[j ^ ax]) * psi)
    ax_grid = np.repeat(np.arange(dim), dim).reshape(dim, dim)
    az_grid = np.tile(np.arange(dim), dim).reshape(dim, dim)
    phase = (1j) ** (popcount(ax_grid & az_grid) % 4)
    full = out * phase
    if n and np.max(np.abs(full.imag)) > 1e-9:
        raise ValueError("Pauli expectation values are not real; input is not a normalised state")
    return full.real.reshape(-1)


def apply_local_noise(r: FloatArr, n: int, p_dephase: float = 0.0, p_depol: float = 0.0) -> FloatArr:
    """Damp a Pauli spectrum by independent single-qubit noise.

    ``p_dephase`` is the probability of a Z flip on each qubit (factor
    ``1 - 2 p`` per X or Y factor).  ``p_depol`` is the probability that each
    qubit is replaced by the maximally mixed state (factor ``1 - p`` per
    non-identity factor).
    """
    _, xw, sw = pauli_index_tables(n)
    return r * (1.0 - 2.0 * p_dephase) ** xw * (1.0 - p_depol) ** sw
