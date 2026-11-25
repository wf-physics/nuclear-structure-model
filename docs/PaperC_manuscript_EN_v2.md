# A Structural Gate on the Mass-Table Pipeline: Stability and Decay-Mode Verdicts for A≤40 Nuclides, Four-Protocol Blind Test, and Net-Contribution Disclosure

**Paper C (final manuscript, v1.2)** (short communication, revised 2026-11-24)
Rule system v5.0 (frozen). This paper and Paper I (the merged structure/method/energy-formula manuscript) form a two-paper series; this paper covers the verdict layer.

---

## Abstract

We use the structural output of the Paper I four-bond saturated network model as a verdict gate for nuclear stability and decay mode, layered on top of the AME2020 mass-table pipeline, and issue verdicts over the full A≤40 nuclide set: **stability verdicts 99.7% (315 nuclides, zero stable-nucleus misjudgment), decay-mode verdicts 99.0% (311/314)**. Three baselines are disclosed proactively: pure Q-sign criterion 92.6%, zero-structure pipeline (mass table + generic rates) 97.5%, full model 99.0% — **the model's net contribution is 1.5 percentage points**, which we recommend readers take as the reference number.

Four-protocol blind-test endorsement: random-individual blind test 97.4±1.1% (99.2% after removing patch-dependent nuclides), whole-element blocking loses nothing, extrapolation blind test is statistically indistinguishable from full-data calibration, and training-set label permutation collapses to 30.1%. The valve test of energy-side independence fails under all four gate×calibration conventions (33.9/57.2/66.7/66.7% against a pre-registered 90% red line) — this verdict layer makes no "independent prediction" claim; the correctness of 12 drip-line nuclides comes from NUBASE measured-half-life patches, not from the model, and the list is published. Comparison with other decay-mode classification methods (empirical rules and machine-learning classifiers) is discussed in §5.

---

## §1 The Verdict Chain

**Trigger layer (zero free parameters)**: all decay channels (β⁻/β⁺/EC/n/p/α) are opened or forbidden by the mass-balance Q value; Q is computed directly from AME2020 mass excesses. $Q \le 0$ strictly forbids — this is the global criterion derived from the "monotonicity impasse" counterexample (Paper I §1.3: five local trigger variants all eliminated by the $^2$H/$^3$H minimal control).

**Rate layer**: β channels via the exact phase-space integral $f(Z_d, Q)$ (Fermi-function sign convention $Z_d = Z+1$ for β⁻, $Z_d = -(Z-1)$ for β⁺); particle channels via WKB barrier penetration; the fastest is the dominant mode. The log $ft$ baseline is anchored to $^3$H = 3.049 (the unique combination of allowed transition + precise experiment + membership of the calibration set).

**Structural gate (model-native)**: particle-channel significance requires model-structure confirmation — the minimum-bond constraint $Wb_{\min}$ and the Coulomb constraint $EC_{\max}$ are output by the structure-generation engine, not fed in from the energy side. This layer is the main source of the model's net contribution.

## §2 Full Results and Three Baselines

| Verdict | Sample | Correct | Accuracy |
|---|---|---|---|
| Stability (stable/unstable binary) | full A≤40 set, 315 | 314 | **99.7%** (zero stable-nucleus misjudgment) |
| Decay mode (specific mode included) | full A≤40 set, 314 (excluding "other") | 311 | **99.0%** |

Three baselines disclosed in the same table:

| Pipeline | Accuracy | Difference |
|---|---|---|
| Pure Q sign (zero structure) | 92.6% | — |
| Mass table + generic rates (zero structural gate) | 97.5% | +4.9 |
| **Full model (+ structural gate)** | **99.0%** | **+1.5 (net contribution)** |

Baseline 1 (fully correct region): all Z<10 isotopes 124/124, zero-error boundary measured to A≤21 (137/137). Baseline 2 (Z 10–20): the error-nuclide list is directly exportable from the `ok=false` entries of the archive `results/baseline_judge.json` (per-nucleus verdict records, script `src/baseline_check.py`), no separate table needed. Baseline 3 (heavy-nucleus spot check 22/22): representative decay-chain nuclei ($^{208}$Pb, $^{238}$U, $^{232}$Th, etc.), the spot-check rule being chain representativeness rather than randomness.

**Round-1 78.4% and falsification-driven iteration**: the initial criterion (local trigger variants) scored only 78.4% on A≤40; its failure (the monotonicity impasse) directly yielded the Q-value trigger layer — this paper reports the iteration endpoint, with all intermediate records archived (`theory/`).

## §3 Four-Protocol Blind Test

| Protocol | Design | Result |
|---|---|---|
| P1 | random 70/30 × 30 seeds, parameters fitted on training set only | **97.4±1.1%** (99.2% after removing patch nuclides) |
| P2 | element blocking (whole Z chain held out) | 97.1±1.3% (no loss → no genealogy memory) |
| P3 | extrapolation blind test (light calibration → heavy prediction) | 96.75% (mid, 10 seeds) |
| P4 | training-set label permutation | **30.1%** (model cannot memorize random labels) |

Parameter sensitivity: under purely random parameters the A≤40 accuracy is 91.2% — about 91% of the predictive power is carried by the fixed theoretical structure, with parameter fitting contributing ~6 pp.

## §4 Valve Test (the Measured Boundary of Energy-Side Independence)

Replacing the $S_n$ energy side from AME with the Paper I closed form (pre-registered three-tier criterion), the four-convention matrix:

| | gate=AMEerr | gate=0.3 MeV |
|---|---|---|
| raw closed form | 33.93% | 57.19% |
| training-region linear calibration | 66.74% | 66.74% |

All four cells are far below the 90% pre-registered red line (protocol details and failure mechanism in Paper I §10). **Verdict (per pre-registration)**: abandon the "independent prediction" wording; verdict-layer positioning = AME mass-table pipeline + model structural gate, net contribution +1.5 pp. The correctness of 12 drip-line patch nuclides does not come from the model (the AME-extrapolated mass uncertainty drowns the separation-energy sign), and the list is published.

## §5 Comparison with Other Decay-Mode Classification Methods

This section compares only **criterion structure**, not precision ranking — different works use different sample domains, stability thresholds, and verdict bands, so cross-sample precision comparison is meaningless (Paper I §4.4 treats it the same way).

**5.1 Empirical rules = a refinement of this paper's Baseline 1.** The Geiger–Nuttall empirical law of α decay [4] and its later systematizations (Viola–Seaborg [5]) are essentially "Q value as input, half-life as output" empirical mappings. This paper's "pure Q-sign" baseline (92.6%) is the simplest form of the same class of criterion — taking only the sign of Q and discarding all empirical fitting. Hence 92.6% is not a weak baseline but the **lower-bound form** of empirical rules on the "stability binary" task: any empirical rule that does not significantly beat the Q sign on this task has no contribution from its extra fitted terms.

**5.2 Macroscopic–microscopic models and the "mass table + rules" pipeline.** Macroscopic–microscopic models of the FRDM2012 [3] class give ground-state masses and deformations, with mass precision (RMS≈0.7 MeV) far better than this paper's closed form (MAE 12.6 MeV, Paper I §4.4); but on the "stability/mode verdict" task the model output is still a mass, and verdicts must still layer on Q-sign rules and empirical systematization — exactly the construction of this paper's Baseline 2 (mass table + generic rates, 97.5%). **The difference between Baseline 2 and the full model (+1.5 pp) is the structural gate's net contribution**, independent of the mass table's own precision.

**5.3 Machine-learning classification: this paper makes no numerical comparison.** Machine-learning classifiers (random forests, neural networks, etc.) fed by mass-table-derived quantities have produced much work on nuclide stability and decay-mode classification. **This paper has not reproduced any ML baseline on the same sample domain and verdict band, and therefore gives no cross-literature precision comparison** — placing another paper's numbers next to one's own is the easiest way to be flagged by referees for inappropriate comparison. What can be said is structural: the feature importances of such work usually fall on Q values and separation energies, homologous to this paper's trigger layer; **if their precision significantly exceeds this paper's Baseline 2 (97.5%), then mass-table-derived nonlinear combinations carry more information than this paper's linear criterion — that would be a direct challenge to the structural gate, not a support.**

**5.4 This paper's differences and costs.** The difference is not precision leadership (full-set 99.0%, only 1.5 pp above the zero-structure pipeline) but two things: ① **every criterion has an explicit physical source** — the trigger layer is Q-value mass balance (derived from the $^2$H/$^3$H counterexample, Paper I §1.3), the gate is the structure engine's bond and Coulomb constraints, auditable nucleus-by-nucleus and traceable to specific nucleon configurations; ② **negative results and the failure list are proactively published** — the valve test's four-convention all-fail (§4), the 12 drip-line patch nuclide names, and per-nucleus error records. The costs are equally explicit: drip-line verdicts rely on NUBASE patches (correctness not from the model), and half-lives are order-of-magnitude only (§6). **A classifier that reports only success rates cannot be falsified; this paper's criterion chain can.**

## §6 Boundaries and Falsification

1. Drip-line verdicts rely on NUBASE patches (12 nuclei, list public) — without the patches this subset has no verdict capability;
2. Half-lives are order-of-magnitude precision (log ft regression R²=0.138; strongly hindered transitions off by >3 orders of magnitude);
3. WKB α half-lives are systematically ~5 orders high ($^{238}$U: model $4.3\times10^{14}$ y vs experimental $4.5\times10^9$ y), mode verdict unaffected;
4. Falsifiable statement: any stable nuclide misjudged unstable by this pipeline, and not attributable to a mass-data boundary, falsifies the trigger layer.

## §7 Conclusion

Layering the model structural gate on the mass-table pipeline gives A≤40 stability 99.7% and decay-mode 99.0%; the +1.5 pp net contribution is disclosed proactively; the four-protocol blind test confirms the predictive power comes from the fixed structure (permutation control 30%); energy-side independence is empirically negated by the valve test (four conventions all-fail). The honest positioning of this verdict layer: **the mass table supplies the energy, the model supplies the structural gate, and the two assembled give a usable verdict that is not an independent prediction.**

---

## References

[1] AME2020: M. Wang, W.J. Huang, F.G. Kondev, G. Audi, S. Naimi, Chin. Phys. C 45, 030001 (2021) (Q values and masses).

[2] NUBASE2020: F.G. Kondev, M. Wang, W.J. Huang, S. Naimi, G. Audi, Chin. Phys. C 45, 030002 (2021) (decay labels, half-lives, drip-line patches).

[3] FRDM2012: P. Möller, A.J. Sierk, T. Ichikawa, H. Sagawa, At. Data Nucl. Data Tables 109–110, 1 (2016) (§5.2 control).

[4] H. Geiger, J.M. Nuttall, Phil. Mag. 22, 613 (1911) (§5.1 empirical rules).

[5] V.E. Viola, G.T. Seaborg, J. Inorg. Nucl. Chem. 28, 741 (1966) (§5.1 α-decay systematization).

[6] This series: Paper I, the four-bond saturated network model (structure/method/energy-formula layer), submitted together (the §1 trigger-layer derivation and §4 valve pre-registration clause are cited from it).

> Citation principle: this paper cites only references "actually used for comparison or as data sources"; unreproduced machine-learning baselines are not cited by name, to avoid forming an inappropriate cross-sample numerical comparison (§5.3).

---

## Data and Code

Verdict pipeline `repo/src/blind_test.py`, `chain_predict.py`, `valve_test.py`, `baseline_check.py`; result archives `repo/results/blind_results.json` (P1–P4 per-seed records), `baseline_judge.json` (314-nuclide per-nucleus verdicts, for the §2 Baseline-2 error list), `chain_results.json` (round-1 78.4% criterion archive), `valve_test.json`, `t4_err_sensitivity.json`. One-command reproduction in `repo/REPRODUCE.md` (per-number mapping). **The code is released with the paper under the MIT license** (repo/LICENSE); the Zenodo DOI and repository URL are filled at submission and cross-reference each other.

---

> **Author archive (not for submission)**: this manuscript is the polished final version of v1.1. Principal revisions: §5 changed from placeholder to a complete comparison section (empirical rules / macroscopic–microscopic / machine learning, with ML explicitly declaring "no same-task reproduction, no numerical comparison", eliminating fabricated baselines); six references added (v1.0 had no reference list, a submission-blocking defect); §2 Baseline-2 error list given an executable traceability path.
