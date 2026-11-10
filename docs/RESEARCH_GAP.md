# Research gap

## Scope of the provided literature

The literature supplied for this study contains **two** papers (one review, one adjacent paper). The brief therefore
does not allow a ten-paper related-work table drawn from three roles; the table below lists what is available and no
reference outside the brief is cited. Standard techniques used in the study are named in a separate table without a
citation so that the original sources can be added by the author.

## State of the art, as far as the provided papers show

* Wolfowicz et al., 2021 (DOI 10.1038/s41578-021-00306-y), a review of solid-state spin defects, shows that registers
  built from spin defects are limited by dephasing of the electron and nuclear spins and by imperfect optical readout,
  and that these properties are the engineering levers for quantum communication, sensing and computing.
* Umair et al., 2027 (DOI 10.1016/j.bspc.2026.111685), an adjacent study of federated learning (FL) with quantum
  classifiers for electroencephalography (EEG), states that centralised training raises privacy concerns, and that it
  is unclear whether quantum models can match classical deep learning when the data stay on the client.

## Related work table

| Paper | Year | Role | What it does | What it does NOT cover |
|---|---|---|---|---|
| Wolfowicz et al. (10.1038/s41578-021-00306-y) | 2021 | review | Maps spin, optical, charge and material requirements of solid-state spin defects for quantum information | No statistical test of a quantum resource under the dephasing and readout noise it describes; no copy-cost accounting |
| Umair et al. (10.1016/j.bspc.2026.111685) | 2027 | adjacent | Federated training of a quantum support vector machine, variational classifier and quantum convolutional network on five clients | Does not ask how many state copies a client needs to verify a property of its own quantum data; no noise-aware decision threshold |

Techniques used here but **not** part of the provided literature (named without citation): Bell difference sampling and
the six-copy stabilizer test built on it, random Clifford sampling, the stabilizer 2-Rényi entropy, and classical
shadows with random Pauli measurements.

## The gap

A property-testing question sits underneath both papers: a device owner holds *noisy copies* of a state, may not
be able to read it out classically (privacy, or simply no tomography budget), and must decide a yes/no property.
The Bell-difference-sampling test decides "stabilizer state versus not" with a **constant** number of copies in the
noiseless, ideal setting. Under the dephasing and readout noise that dominate spin-defect registers, a noisy stabilizer
state no longer has acceptance probability one, so the textbook decision threshold is miscalibrated. What is missing
is a quantitative, exact account of (i) how much noise breaks a fixed threshold as the register grows, (ii) what a
threshold repaired with a noise model or with an on-device reference costs in copies, and (iii) how this compares
with a single-copy classical-shadow tester that needs no entangling measurement.

**Novelty in one sentence:** this study gives exact acceptance statistics of the noisy Bell-difference-sampling
stabilizer test (dephasing, depolarising, readout), the critical noise at which a fixed threshold fails as a function of
qubit number, and the copy cost of noise-aware and reference-calibrated thresholds against a classical-shadow baseline.

## Hypotheses and success metric

Setting: n = 3 to 6 qubits, exact simulation; null = Haar-random stabilizer states, alternative = random Clifford
circuits containing one T gate (stabilizer 2-Rényi entropy 0.415 bit).

* **H1.** A noise-agnostic threshold placed from noiseless design values yields more than 10% balanced error even with
  unlimited copies once the noise strength p (dephasing = depolarising = readout flip = p) exceeds a critical p\_crit
  that decreases with n. *Falsified if* p\_crit is flat or increasing in n.
* **H2.** A noise-aware threshold restores a finite copy cost, and the Bell test then needs fewer copies than the
  single-copy shadow tester at n = 3 to 5. *Falsified if* the shadow tester reaches 10% error with fewer copies.
* **H3.** A threshold calibrated on a noisy stabilizer reference (no label information) recovers most of the
  noise-aware accuracy at equal total copies. *Falsified if* its error at a fixed budget stays above half of the
  noise-agnostic error.

Success metric: number of state copies needed to reach 10% balanced error (mean of false-reject and false-accept
rates), reported as mean ± standard deviation over five seeds.
