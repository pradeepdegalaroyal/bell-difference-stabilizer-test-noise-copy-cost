"""Instance generators: random stabilizer states and Clifford+T (magic) states."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from qiskit.quantum_info import Operator, random_clifford

from .paulis import FloatArr, pauli_spectrum

_T_PHASE = np.exp(1j * np.pi / 4)


def random_clifford_matrix(n: int, rng: np.random.Generator) -> NDArray[np.complex128]:
    """Dense unitary of a uniformly random ``n``-qubit Clifford."""
    seed = int(rng.integers(0, 2**31 - 1))
    return Operator(random_clifford(n, seed=seed)).data


def random_stabilizer_state(n: int, rng: np.random.Generator) -> NDArray[np.complex128]:
    """Uniformly random stabilizer state, ``C|0...0>``."""
    return random_clifford_matrix(n, rng)[:, 0].copy()


def random_magic_state(n: int, n_t: int, rng: np.random.Generator) -> NDArray[np.complex128]:
    """``C2 (T on n_t random qubits) C1 |0...0>`` with Haar-random Cliffords."""
    psi = random_stabilizer_state(n, rng)
    qubits = rng.choice(n, size=n_t, replace=False)
    idx = np.arange(2**n)
    for q in qubits:
        psi = psi * np.where((idx >> int(q)) & 1, _T_PHASE, 1.0)
    return random_clifford_matrix(n, rng) @ psi


def stabilizer_renyi_m2(r: FloatArr, n: int) -> float:
    """Stabilizer 2-Renyi entropy ``M2 = -log2( sum_a r_a**4 / 2**n )`` of a pure state.

    ``M2 = 0`` exactly for stabilizer states and is strictly positive otherwise.
    """
    return float(-np.log2(np.sum(r**4) / 2**n))


def draw_ensembles(
    n: int, n_t: int, size: int, rng: np.random.Generator, min_m2: float = 0.05
) -> tuple[list[FloatArr], list[FloatArr], list[float]]:
    """Draw ``size`` stabilizer and ``size`` magic Pauli spectra.

    Magic draws whose ``M2`` is below ``min_m2`` (for example a T gate landing
    on a qubit in a computational eigenstate) are redrawn, so every
    alternative is a genuinely non-stabilizer state.
    """
    null = [pauli_spectrum(random_stabilizer_state(n, rng)) for _ in range(size)]
    alt: list[FloatArr] = []
    m2s: list[float] = []
    while len(alt) < size:
        r = pauli_spectrum(random_magic_state(n, n_t, rng))
        m2 = stabilizer_renyi_m2(r, n)
        if m2 >= min_m2:
            alt.append(r)
            m2s.append(m2)
    return null, alt, m2s
