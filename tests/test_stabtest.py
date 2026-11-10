"""Fast correctness tests for the exact Bell-test statistics and the shadow baseline."""

from __future__ import annotations

import numpy as np
import pytest

from stabtest.bell_test import (
    acceptance_probability,
    balanced_error_curve,
    bell_distribution,
    copies_needed,
)
from stabtest.paulis import apply_local_noise, pauli_spectrum, wht
from stabtest.shadow_baseline import shadow_statistic
from stabtest.states import random_magic_state, random_stabilizer_state, stabilizer_renyi_m2

X = np.array([[0, 1], [1, 0]], dtype=complex)
Z = np.diag([1.0, -1.0]).astype(complex)
I2 = np.eye(2, dtype=complex)


def pauli_matrix(ax: int, az: int, n: int) -> np.ndarray:
    m = np.eye(1, dtype=complex)
    for i in reversed(range(n)):
        a, b = (ax >> i) & 1, (az >> i) & 1
        o = I2
        if a:
            o = o @ X
        if b:
            o = o @ Z
        m = np.kron(m, (1j) ** (a & b) * o)
    return m


def test_wht_matches_definition() -> None:
    v = np.random.default_rng(0).normal(size=8)
    ref = np.array([sum(v[j] * (-1) ** bin(u & j).count("1") for j in range(8)) for u in range(8)])
    assert np.allclose(wht(v), ref)


def test_pauli_spectrum_matches_matrices() -> None:
    n = 2
    psi = random_magic_state(n, 1, np.random.default_rng(3))
    r = pauli_spectrum(psi)
    for ax in range(4):
        for az in range(4):
            assert r[ax * 4 + az] == pytest.approx((psi.conj() @ pauli_matrix(ax, az, n) @ psi).real, abs=1e-10)


def test_bell_distribution_matches_direct_bell_measurement() -> None:
    """Closed-form Bell law (including noise) equals <Phi_x| rho (x) rho |Phi_x>."""
    n, dim = 2, 4
    psi = random_magic_state(n, 1, np.random.default_rng(1))
    r = apply_local_noise(pauli_spectrum(psi), n, 0.05, 0.1)
    rho = sum(r[a * dim + b] * pauli_matrix(a, b, n) for a in range(dim) for b in range(dim)) / dim
    big = np.kron(rho, rho)
    direct = np.zeros(dim * dim)
    for ax in range(dim):
        for az in range(dim):
            pm = pauli_matrix(ax, az, n)
            phi = sum(np.kron(np.eye(dim)[j], pm @ np.eye(dim)[j]) for j in range(dim)) / np.sqrt(dim)
            direct[ax * dim + az] = (phi.conj() @ big @ phi).real
    assert np.allclose(direct, bell_distribution(r, n), atol=1e-12)
    assert direct.sum() == pytest.approx(1.0)


@pytest.mark.parametrize("n", [2, 3, 4])
def test_stabilizer_states_always_accepted(n: int) -> None:
    rng = np.random.default_rng(n)
    r = pauli_spectrum(random_stabilizer_state(n, rng))
    assert acceptance_probability(r, n) == pytest.approx(1.0, abs=1e-9)
    assert stabilizer_renyi_m2(r, n) == pytest.approx(0.0, abs=1e-9)


def test_magic_state_rejected_with_known_probability() -> None:
    """One T gate on a random stabilizer state gives M2 = 0.415 and acceptance 13/16 for n = 3."""
    n = 3
    rng = np.random.default_rng(5)
    for _ in range(5):
        r = pauli_spectrum(random_magic_state(n, 1, rng))
        if stabilizer_renyi_m2(r, n) > 0.05:
            assert acceptance_probability(r, n) == pytest.approx(0.8125, abs=1e-9)
            return
    pytest.fail("no magic instance drawn")


def test_noise_and_readout_lower_acceptance() -> None:
    n = 3
    r = pauli_spectrum(random_stabilizer_state(n, np.random.default_rng(0)))
    base = acceptance_probability(r, n)
    assert acceptance_probability(apply_local_noise(r, n, 0.01, 0.0), n) < base
    assert acceptance_probability(r, n, readout=0.01) < base


def test_error_decreases_with_rounds_and_copy_count_finite() -> None:
    null, alt = np.array([0.95, 0.93]), np.array([0.80, 0.82])
    e_small = balanced_error_curve(null, alt, 5)[0]
    e_big = balanced_error_curve(null, alt, 200)[0]
    assert e_big < e_small
    assert np.isfinite(copies_needed(null, alt))
    assert copies_needed(np.array([0.8]), np.array([0.81])) > copies_needed(null, alt)


def test_shadow_statistic_unbiased_for_stabilizer_state() -> None:
    n = 2
    rng = np.random.default_rng(2)
    r = pauli_spectrum(random_stabilizer_state(n, rng))
    vals = [shadow_statistic(r, n, 40000, rng) for _ in range(20)]
    assert np.mean(vals) == pytest.approx(1.0, abs=0.1)
