# Method

## Task

Given copies of an unknown n-qubit pure state that has passed through local noise, decide *stabilizer* (null) versus
*not a stabilizer state* (alternative). Instances: the null is a Haar-random stabilizer state `C|0>`; the alternative
is `C2 T C1 |0>` with one T gate on a random qubit and M2 >= 0.05, where M2 is the stabilizer 2-Rényi entropy
(`M2 = -log2(sum_a r_a^4 / 2^n)` over Pauli expectations `r_a`). Each seed draws 24 states per class.

## Noise model (applied identically to every copy)

* dephasing, probability p\_z of a Z flip per qubit (spin-defect-like dominant error),
* depolarising, probability p\_d of replacing a qubit by the maximally mixed state,
* readout flip, probability q per measured bit of the Bell and Pauli measurements.

Noise acts on Pauli expectations multiplicatively, so everything below is computed *exactly* from the 4^n-entry Pauli
spectrum (using a fast Walsh-Hadamard transform, WHT). There is no gate noise inside the Bell-measurement circuit;
readout flips are the only measurement-side error (see Limitations in the README).

## Methods compared (equal total number of state copies)

1. **Bell-difference-sampling test (quantum, standard).** One round uses six copies: two Bell measurements give
   outcomes `x1, x2`; `x = x1 xor x2` (bitwise exclusive-or, XOR); Pauli `P_x` is measured on two more copies; the round accepts if the two
   outcomes agree. Exact acceptance probability per state: `sum_x D(x) (1 + g(x) r_x^2) / 2`, with `D` the XOR-convolution
   of the Bell outcome law `B(x) = 4^-n sum_a s_a (-1)^<x,a> r_a^2` with itself and `g(x) = (1 - 2q)^(2|x|)`. The test
   statistic is the number of accepted rounds out of R = copies / 6.
2. **Pauli-shadow tester (classical, single copy).** Each copy is measured once in a random product Pauli basis; four
   independent batches estimate every `r_a`, and `T = 2^-n sum_a (r_a^(1) r_a^(2) r_a^(3) r_a^(4))` is an unbiased
   estimate of `2^-n sum_a r_a^4` (equal to 1 for pure stabilizer states).

## Decision thresholds (the ablated component)

* **noise-agnostic:** accept if the acceptance rate is at least `(1 + mean noiseless alternative acceptance) / 2`.
* **noise-aware (oracle):** the best state-independent count threshold for the noisy ensemble. It uses ensemble labels
  and is therefore an upper bound on what a deployable rule can do; the shadow tester receives the same oracle treatment.
* **reference-calibrated:** half of the copies run the test on a fresh noisy stabilizer reference; the threshold is the
  reference acceptance rate minus half of the noiseless null-minus-alternative gap. Label-free.

## Metrics

Balanced error = 0.5 (false-reject rate on the null ensemble + false-accept rate on the alternative ensemble), exact in
the binomial law for the Bell test and Monte Carlo (6 trials per state, 12 states per class) for the shadow tester.
Copies needed = smallest budget with balanced error <= 0.10 (Bell: refined to about 2%; shadows: grid of factor 2).
Critical noise p\_crit = smallest p at which the asymptotic (unlimited-copy) error of the agnostic threshold exceeds 10%.

## Experiments

| Study | Sweep | Output |
|---|---|---|
| Error versus copies | 60 to 15360 copies, n = 4, p = 0.005 per source | `results/curves.csv`, Figure 1 |
| Noise heatmap | dephasing x readout flip, 5 x 5 grid, n = 4 | `results/heatmap.csv`, Figure 2 |
| Scaling | n = 3 to 6 (shadows 3 to 5), noiseless and p = 0.005 | `results/scaling.csv`, `results/critical_noise.csv`, Figure 3 |
| Ablation (waterfall) | add noise sources one by one, then repair the threshold, at 600 copies | `results/waterfall.csv`, Figure 4 |

Seeds: 0 to 4 (ensembles are re-drawn per seed). Mean ± standard deviation over seeds throughout.
