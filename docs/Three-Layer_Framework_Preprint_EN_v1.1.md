# A Three-Layer Division-of-Labor Framework for Nuclear Structure: Liquid-Drop Energy, Shell Gaps, and Spatial Correlations from a Four-Bond Saturated Network, with Application to ²³⁸U

**Preprint v1.1** · 2026-11-23

> Core thesis of the abstract: nuclear structure is not described by any single model replacing another, but by three physical scales each doing its own job in orthogonal complementarity — **the liquid-drop model provides the macroscopic energy skeleton, the shell model provides quantum gaps, and a four-bond saturated network provides few-body spatial correlations**. This paper gives quantitative evidence for this three-layer division of labor, and uses ²³⁸U to test the framework's reproduction accuracy and current boundaries in the heavy-nucleus region.

---

## Abstract

Nuclear-structure physics has long hosted two complementary yet mutually independent descriptions: the liquid-drop/macroscopic model, which characterizes mean binding energy and deformation as a continuous fluid, and the shell model, which characterizes magic-number gaps through single-particle orbitals. Between the two there is a missing layer: **few-body spatial correlations** (α clustering, quadrupole-moment sign, reaction chains), quantities that neither of the former possesses structurally. This paper proposes and systematically tests a three-layer division-of-labor framework in which that missing layer is supplied by a "four-bond saturated network."

The geometric layer contains only two local rules: a windowed linear bond well between nucleons, superposed with Pauli short-range repulsion; and a maximum of four bonds per nucleon, with bonds formed only by mutual consent (each endpoint must have a vacancy). This minimal rule set converges, in the infinite nuclear-matter limit, to a degenerate family of ideal tetrahedral four-coordinated networks, whose tetrahedral units are homologous to light-nucleus α clusters. Under a zero-per-nucleus-fitting constraint, liquid-drop coefficients emerge from the model energetics (a_v=16.060, a_s=19.010, a_c=0.7560, a_sym=19.50, identical across the two extraction routes — which share a common anchor, a self-consistency rather than an independence claim; §2.3), and linear differences of the same closed form yield four release-energy formulas — neutron separation, β⁻, α, and fission — each with a measured range of validity (S_n MAE 1.58 MeV; Q_β MAE 3.80 MeV; Q_α MAE 1.32 MeV; Q_fs 5/5 within band).

Quantitative evidence for the three-layer division is: the liquid drop achieves an RMSE of 3.33 MeV on binding energies with five global parameters (experimental-mass subset; 3.67 MeV on the full set of §3.3), versus 81.0 MeV for the geometric engine (a factor of 24), so energy belongs to the liquid drop; shell corrections of magic nuclei show a clear positive peak (⁴⁸Ca doubly magic, +5.80 MeV maximum) while geometric structural quantities correlate with them at ρ≈0, so shell gaps belong to the shell model; quadrupole-moment sign agreement of 69.2% and α-cluster geometric counts are quantities absent from every version of the liquid-drop model, so structure belongs to the geometric layer.

As a concrete application, the paper takes ²³⁸U (Z=92, N=146) as a heavy-nucleus case: under a global refit of liquid-drop + pairing + Myers–Swiatecki shell correction (six global parameters, 3487 nuclei, MAE 2.297 MeV), the ²³⁸U binding-energy deviation falls from +19.06 MeV to **+1.72 MeV (0.10%)**, the one-neutron separation energy from +1.51 to **+0.39 MeV (6.3%)**, and the Pd+Pd fission channel from +2.20 to **+0.35 MeV (0.18%)**; the α-decay energy falls from +2.71 to +1.21 MeV (28%), its residual deviation attributed to the absence of deformation energy in the spherical shell correction. A frozen holdout validation — all nuclei with Z≥88 removed before refitting the same six parameters — degrades these sample-internal numbers only to +5.44 MeV (0.30%) for the binding energy and −4.57 MeV (−2.40%) for the fission channel, with S_n and Q_α essentially unchanged (§3.4); the holdout numbers, not the in-sample ones, are the framework's first true predictions for ²³⁸U. All boundaries, failure lists, and the means of obtaining reproducible materials are disclosed in the text.

---

## §1 Introduction: Three Scales and One Missing Layer

### 1.1 Two Existing Pillars of Nuclear Structure

Nuclear structure spans at least three scales. At the macroscopic scale the nucleus approximates an incompressible liquid drop, with binding energy described by volume, surface, Coulomb, symmetry, and pairing terms (liquid-drop model [8]) at a precision of roughly 1–3 MeV, yet structurally devoid of any "few-body correlation" concept. At the quantum single-particle scale, nucleons fill orbitals in a mean potential, and magic numbers arise from the gaps in level filling (shell model), giving extra binding at magic numbers, but its computations (configuration mixing, quadrupole moments) are costly. The two are irreplaceable at their own scales, yet both fall silent on the intermediate scale of "spatial correlations."

### 1.2 The Missing Layer: Few-Body Spatial Correlations

Observables that are clearly present experimentally but that neither the liquid drop nor the shell model can directly provide structurally include at least:

1. **α clustering**: 2p2n correlations in N=Z even-even 4n nuclei (⁸Be, ¹²C, ¹⁶O, ²⁰Ne) [7]; the liquid drop describes these as a continuous fluid (no cluster concept), the shell model via single-particle orbitals (requiring complex configuration mixing to reveal).
2. **Quadrupole-moment sign**: the oblate/prolate shape of the proton distribution; the liquid-drop model has no such column, the shell model requires configuration-mixing calculations.
3. **Reaction chains**: the structural prerequisites of α-decay, fission, and other release-energy channels.

This paper proposes that this missing layer is supplied by a **four-bond saturated network model**, forming a three-layer orthogonal complementarity with the liquid-drop and shell models.

### 1.3 Core Thesis

> The liquid drop governs "mean energy," the shell governs "quantum gaps," the geometry governs "spatial correlations." The three are complementary **within their current validity domains** (liquid drop: all A; shell correction: A≥16; geometric layer: A≤40, §4.1). The geometric layer does not yet extend beyond A≤40; that is a boundary of the present version, not a division-of-labor conclusion.

This thesis is not a slogan but a quantitative proposition established with three independent datasets (§4), and anchored by the ²³⁸U application (§6).

---

## §2 The Geometric Layer: A Four-Bond Saturated Network Model

### 2.1 Minimal Rules

The model contains only two interaction rules (rule text finalized before any validation began, unmodified during validation):

- **Windowed linear bond well**: V_att(r) = −s·(r₂−r), r ∈ [r₁, r₂], zero outside; superposed with **Pauli short-range repulsion** V_rep(r) = C·max(0, r_eq−r)².
- **Bond-number saturation**: at most k_max = 4 bonds per nucleon; **mutual consent**: a bond forms only if both endpoints have a vacancy.

Mapping to empirical facts: the attractive window corresponds to the short range of the nuclear force (F3); saturation and mutual consent correspond to density saturation and constant binding energy (F1/F2); k_max=4 corresponds to the nucleon spin-isospin quartet and the light-nucleus magic number 4 (F5/F6). The model does not input shell structure, spin-orbit coupling, pairing interactions, density-dependent interactions, or any per-nucleus experimental data (see §7 boundaries).

### 2.2 The Diamond Limit and α-Cluster Homology

In the infinite nuclear-matter limit, k_max=4 and mutual consent jointly cap the coordination number at 4, reducing the problem to optimization of a four-coordinated network. All "ideal tetrahedral" networks (equal bond lengths, equal 109.47° bond angles) are **strictly degenerate** under this potential; cubic diamond is a representative member of the family, not the unique ground state. The nuclear-matter saturation point is calibrated to ρ0=0.160 fm⁻³, E/A=−15.8 MeV (the calibration targets are empirical values, so achievement is an acceptance criterion rather than an independent prediction); the compression modulus K=547 MeV has an analytic lower bound K_att = 4·k_n·s·c·ρ0^(−1/3) rooted in window freezing, and the framework does not predict giant-resonance energies.

The diamond tetrahedral unit is the geometry of a 4-nucleon α cluster. This homology between the "spin-isospin quartet k_max=4 + mutual consent" hypothesis and the empirical light-nucleus α clustering is the first cornerstone of the geometric layer: the same rules give a tetrahedral network in the infinite limit and a detectable tetrahedral sub-configuration in finite N=Z nuclei (§2.4).

### 2.3 Emergent Liquid-Drop Coefficients

Under a zero-per-nucleus-fitting constraint, the coefficients of the four-parameter closed form B(A,Z) = a_v·A − a_s·A^(2/3) − a_c·Z²/A^(1/3) − a_sym·(N−Z)²/A are derived from the model energetics via two mutually independent routes:

| Route | a_v | a_s | a_c | a_sym | Method |
|---|---|---|---|---|---|
| 1 (rescaled cluster sum) | 16.060 | 19.010 | 0.7560 | 19.50 | diamond-cluster sum + well-depth rescaling |
| 2 (model-grid calibration) | 13.851 | 1.592 | 0.7265 | 19.50 | four-bond network full energy + least squares |

**a_sym=19.50 is identical between the two routes.** We state candidly what this number is and is not: both routes extract the isovector coefficient from the *same* calibrated potential — route 1 via a cluster-count regression after the pre-registered well-depth rescaling (×2.61 to the saturation point), route 2 via a least-squares fit to the same model's grid energies. The shared anchor is the same potential and its same np/pp/nn bond-energy asymmetry, so the exact agreement demonstrates that the two extraction procedures are **self-consistent under a common anchor**, not that two independent models converge. What independence survives is limited to the extraction side: the routes differ in sampling (cluster sums vs. grid), in fitting procedure, and in the treatment of kinetic and Coulomb terms — and correspondingly, a_v and a_s do *not* agree exactly between them. The route-1 closed form gives an MAE of 12.6 MeV over 2811 nuclei (A=40–240) of the full AME2020 table.

Linear differences of the same closed form yield four release-energy formulas, with measured ranges of validity:

| Observable | Definition | Sample | MAE | Verdict |
|---|---|---|---|---|
| F1 neutron separation S_n | B(A,Z)−B(A−1,Z) | 85 nuclei (mostly A≥40) | 1.58 MeV | 59 pass/25 half/1 fail |
| F2 β⁻ decay Q_β | B(A,Z+1)−B(A,Z)−0.782 | 10-nucleus test set | 3.80 MeV (Lane form) | order-of-magnitude usable |
| F3 α decay Q_α | B(A−4,Z−2)+B_α−B(A,Z) | 14 nuclei (12 with A≥180) | 1.32 MeV (1.06 for A≥180) | 3 pass/9 half/2 fail |
| F4 fission Q_fs | B(f1)+B(f2)−B(parent) | 5 fission channels | 3.83 MeV | 5/5 within band |

### 2.4 Validation of Geometric Observables

Per-nucleus quantitative comparison over the full H→Ca set of 389 nuclei:

| Observable | Criterion | n | Result |
|---|---|---|---|
| Quadrupole-moment sign | sign agreement | 65 | 69.2% (45/65) |
| α-cluster count | ≥50% of ideal cluster count | 10 | 8/10 |
| Charge radius | ratio ∈ [0.7, 1.3] | 97 | 66.0% (median 0.794) |

**The quadrupole-moment sign and the α-cluster count are differential explanatory quantities that the liquid-drop and shell models lack structurally** — this is the irreplaceable position of the geometric layer (quantified in §4).

**Charge radii are a scale problem, not a shape problem.** The model radii are systematically too small: the median model/experiment ratio is 0.794, i.e. the model is too compact by ~21% (root cause: missing quantum zero-point swelling, most visible in the Z≤9 light nuclei, whose ratios fall to 0.35–0.70). Introducing a single radius scale factor f = 0.794 removes the systematic offset: the fraction of nuclei (n=95 of the 97 with both model and experimental values) within ±15% of experiment rises from 30.5% to 52.6%, and the median scaled ratio becomes 1.027. The residual scatter remains substantial (mean |scaled ratio − 1| = 0.29, dominated by the light nuclei and a few unconstrained basin-hopping configurations), so this version reports the scale diagnosis without claiming shape accuracy.

---

## §3 The Liquid-Drop + Shell Layer: A Macroscopic-Microscopic Mass Model

### 3.1 Why Energy Must Belong to the Liquid Drop

The geometric engine's binding-energy RMSE reaches 81.0 MeV (root cause: the engine lacks a surface term, systematically overbinding mid-to-heavy nuclei, with the α-particle binding energy 17% too low — the engine's α is 23.38 MeV versus 28.30 MeV experimental). The liquid drop, with five global parameters, achieves RMSE 3.33 MeV on the experimental-mass subset (2484 nuclei; on the full 3487-nucleus set of §3.3 the same refit gives 3.67 MeV). Handing energy to the liquid drop relieves the geometric model of its least-competent burden, letting it focus on its structural strengths.

### 3.2 Shell Correction: Myers–Swiatecki Form

The liquid-drop residual (B_exp − B_LDM) shows a clear positive peak at magic numbers (§4.2 table), the signature of single-particle level filling. This paper adopts the standard Myers–Swiatecki spherical shell correction [10]:

> δB_shell(Z,N) = C·[ c·A^(1/3) − (F(N)+F(Z))/(A/2)^(2/3) ],
> F(X) = (3/5)[ (M_i^(5/3)−M_{i−1}^(5/3))/(M_i−M_{i−1})·(X−M_{i−1}) − (X^(5/3)−M_{i−1}^(5/3)) ]

with the magic-number set {2, 8, 14, 28, 50, 82, 126, 184}, the fitted amplitude C, and a fixed shape constant **c = 0.25** (frozen; only C enters the fit).

### 3.3 Global Refit

The liquid-drop bulk, pairing, and shell correction must be fitted **jointly and globally** (naive superposition breaks the differences; §6.2 records this failed experiment). Least-squares fitting of six global parameters over all 3487 nuclei (A≥16) of AME2020:

> B(Z,N) = aV·A − aS·A^(2/3) − aC·Z(Z−1)/A^(1/3) − aA·(N−Z)²/A + aP·δ/√A + δB_shell

Fitted result: **aV=15.287, aS=16.307, aC=0.690, aA=22.251, aP=11.049, C_shell=3.226**, global **MAE = 2.297 MeV** (RMS 3.069 MeV). This is six global parameters over 3487 nuclei, **not per-nucleus fitting**.

**Precision comparison on the same nuclei set and the same error definition** (binding-energy deviation against AME2020, in MeV — a comparison the text previously owed to [3], which was cited as a benchmark but never compared against):

| Model | Fitting | Nuclei | MAE (MeV) | RMS (MeV) |
|---|---|---|---|---|
| This work (macro–micro, above) | 6 global parameters, zero per-nucleus fitting | 3487 (A≥16, incl. AME2020 estimates) | 2.297 | 3.069 |
| Pure liquid drop, refit on the identical set | 5 global parameters | 3487, same set | 2.756 | 3.668 |
| FRDM2012 [3] | 10 macroscopic constants + folded-Yukawa microscopic + per-nucleus deformation optimization | 3456 of the same 3487 present in the FRDM2012 table | 0.603 | 0.988 |

Notes: (i) FRDM2012's published rms of 0.56 MeV refers to its own adjustment set (2149 nuclei, 2012 mass database) [3]; evaluated against the same AME2020 set used here, its deviation is RMS 0.988 MeV — still about three times better than this work. (ii) The "RMSE 3.33 MeV" quoted for the pure liquid drop elsewhere in this paper is the experimental-only subset (2484 nuclei); on the identical 3487-nucleus set used in this table the refit gives 3.668 MeV, which is the number quoted here for set consistency.

### 3.4 Holdout Validation: ²³⁸U as a True Out-of-Sample Prediction

All results of §3.3 and §6 share one structural limitation: ²³⁸U (Z=92) is inside the 3487-nucleus fitting set, so its 0.10% binding-energy deviation is a **sample-internal reproduction** number. The following holdout protocol was therefore fixed before being run (rule frozen; the alternative rule "A≥220 removed" was not attempted): **remove all nuclei with Z≥88 from the fitting set, refit the same six global parameters on the remainder, and only then evaluate ²³⁸U.** The refit uses the same implementation as §3.3, validated by reproducing the in-sample fit to MAE 2.296 MeV with all six parameters agreeing to better than 0.5% (c fixed at 0.25 as in the original fit).

Refit on the remaining 2902 nuclei (585 removed): **aV=15.084, aS=15.821, aC=0.671, aA=21.602, aP=11.055, C_shell=3.036**; in-domain MAE 2.012 MeV (RMS 2.864 MeV).

Out-of-sample ²³⁸U deviations (²³⁸U and ²³⁴Th lie outside the holdout domain, Z≥88 removed; the ¹¹⁹Pd fragments lie inside it):

| Observable | In-sample (§3.3 fit) | Holdout (Z≥88 removed) |
|---|---|---|
| B (exp 1801.69 MeV) | +1.72 MeV (+0.10%) | **+5.44 MeV (+0.30%)** |
| S_n (exp 6.15 MeV) | +0.39 MeV (+6.3%) | **+0.47 MeV (+7.6%)** |
| Q_α (exp 4.270 MeV) | +1.21 MeV (+28%) | **+0.95 MeV (+22%)** |
| Q_fs(Pd+Pd) (exp 190.12 MeV) | +0.35 MeV (+0.18%) | **−4.57 MeV (−2.40%)** |

Reading: the binding energy degrades from 0.10% to 0.30% (+1.72 → +5.44 MeV) — the in-sample number should never be quoted without this one beside it. The separation energy is nearly unchanged (+6.3% → +7.6%); the α-decay energy even improves slightly (+28% → +22%), since the refit trades uranium-region binding against global fit quality; the fission channel stays within its differential-cancellation regime (0.18% → 2.40%, still inside the 1–7% band). No quality class collapses. These holdout numbers — not the in-sample ones — are the framework's first genuine predictions for ²³⁸U, and they define the quotable precision boundary: **binding energy at the 0.3% (≈5 MeV) level and fission Q-values at the 2–3% level, from a six-global-parameter model that never saw ²³⁸U.**

---

## §4 The Unified Three-Layer Evidence Table: Quantitative Evidence of Orthogonality

### 4.1 Three-Layer Division Summary

| Layer | Model | Observable | Key number | Current validity domain | Substitutability test |
|---|---|---|---|---|---|
| Macroscopic bulk | Liquid-drop model | Binding energy B | RMSE **3.33 MeV** (exp. subset; 3.67 on full set) | **all A** | Geometric engine 81.0 MeV (24× worse) → irreplaceable |
| Quantum single-particle | Shell model | Shell correction δ | Doubly-magic nuclei +2.4~+5.8 MeV | **A≥16** (fitting domain of §3.3) | Geometric-structure correlation ρ≈0 → irreplaceable |
| Few-body correlation | **Geometric model** | Q0 / α clusters | Q0 sign 69.2%, α clusters 8/10 | **A≤40** (engine limit) | Liquid drop **has no such concept** → irreplaceable |

The validity-domain column states capability boundaries as they are: the geometric layer is entirely absent above A≤40. That absence is a boundary of the present version, not an intended division of labor — the substitutability tests above are meaningful only inside these domains.

### 4.2 The Shell-Model Layer's Magic-Number Signature

| Nucleus | Identity | B_exp | B_LDM | Shell correction δ |
|---|---|---|---|---|
| ¹⁶O | doubly magic 8 | 127.62 | 125.21 | **+2.41** |
| ⁴⁰Ca | doubly magic 20 | 342.05 | 337.89 | **+4.16** |
| ⁴⁸Ca | doubly magic 20/28 | 416.00 | 410.20 | **+5.80** |
| ¹⁸O | Z magic 8 (N non-magic) | 139.81 | 142.64 | **−2.83** |

⁴He is excluded from this table: A=4 lies outside the liquid-drop fitting domain (A≥16, §3.3), so its droplet-based shell correction (δ=+4.28) has no physical meaning. The three doubly-magic rows and the negative ¹⁸O row carry the evidence.

The shell correction δ has a clear structure: doubly-magic nuclei peak highest, singly-magic nuclei lower, and a Z-magic nucleus with non-magic N is negative. This is the signature of the shell model (single-particle level filling), producible by neither the liquid drop nor the geometry.

### 4.3 The Geometric Layer's Exclusive Columns (Absent from the Liquid Drop)

| Nucleus | B_LDM | Shell δ | Q0_geometric | α clusters |
|---|---|---|---|---|
| ⁴He † | 24.02 | +4.28 | 1.41 | 1/1 |
| ¹⁶O | 125.21 | +2.41 | 11.46 | 2/4 |
| ²⁰Ne | 160.93 | −0.29 | 94.07 | 4/5 |
| ⁴⁰Ca | 337.89 | +4.16 | 35.17 | 5/10 |

† ⁴He: the liquid-drop and shell values are outside the fitting domain (A≥16) and are shown only for context; for this nucleus the meaningful entries are the geometric columns.

The `Q0_geometric` and `α clusters` columns do not exist in any version of the liquid-drop model — this is the irreplaceable position of the geometric model.

### 4.4 Conclusion on Three-Layer Orthogonality

Three independent datasets jointly establish: **energy, shell gaps, and spatial correlations are governed by different models and are mutually irreplaceable within the validity domains stated in §4.1** (liquid drop all A, shell correction A≥16, geometry A≤40). The previously passive narrative — "geometric model net contribution only 1.5pp, valves all failed, binding-energy deviation 81 MeV" — is redefined here as **a precise mapping of division-of-labor boundaries**: energy deviations belong to the liquid drop, shell gaps to the shell model, structure to the geometry. One caveat keeps this honest: above A=40 the framework presently has only two layers (liquid drop + shell correction); the geometric layer has not been extended there, and the ²³⁸U card of §6 tests only the two-layer energy side.

---

## §5 Method Overview (Auditability Statement)

The complete data flow of the three-layer framework is:

1. **Geometric layer**: structure-generation engine (diamond-lattice initialization + L-BFGS + basin-hopping global search [6]) → configuration → Q0, α clusters, bond network;
2. **Liquid-drop + shell layer**: AME2020 global refit → binding energy, shell correction, pairing;
3. **Difference layer**: linear differences of the closed form B → S_n / Q_β / Q_α / Q_fs.

The three layers share the same (Z,N) input, each computing its own part without interference. All key numbers correspond to archives and are traceable digit by digit. Complete reproducible code and datasets are described in the code and data availability statement of §8.

---

## §6 Heavy-Nucleus Application: The ²³⁸U Test Card

**Scope disclaimer (read first): this card tests the energy layer — the liquid-drop + shell model of §3.3 and the difference layer of §2.3 — not the geometric layer.** The geometric engine applies only to A≤40, so the structural quantities of ²³⁸U (Q0, α clusters) are outside this card's scope; this is stated at the top so it cannot be missed.

### 6.1 Method

²³⁸U (Z=92, N=146) is representative of heavy, strongly deformed, strongly shell-affected nuclei, and is a stress test of this framework. What is tested is the reproduction accuracy of the **liquid-drop + shell layer** in the heavy region and the release-energy precision of the **difference layer** (α decay, fission); with the holdout of §3.4, the same card also carries the framework's first true out-of-sample predictions for this nucleus.

### 6.2 Two Experiments: Naive Superposition Fails, Global Refit Succeeds

**Experiment 1 (naive superposition)**: superimposing the Myers–Swiatecki shell correction directly on the v5 closed-form liquid drop — B improves (+19.06→+11.15), but Q_α worsens (+2.71→+3.95, +63%→+92%) and Q_fs(Pd+Pd) worsens (+2.20→+6.77). **Cause of failure**: the v5 liquid-drop parameters are not AME2020-fitted values, so direct superposition breaks the differences. The shell correction must be jointly globally refit with the bulk parameters.

**Experiment 2 (global refit)**: the six-parameter joint fit of §3.3. Results below.

### 6.3 Results Compared (as-is)

| Observable | Experiment (AME2020) | v5 closed form (no shell) | Global refit (liquid-drop+shell+pairing) | Improvement |
|---|---|---|---|---|
| Binding energy B | 1801.69 | +19.06 (+1.06%) | **+1.72 (+0.10%)** | 11× |
| S_n one-neutron separation | 6.15 | +1.51 (+25%) | **+0.39 (+6.3%)** | 4× |
| Q_α decay energy | 4.27 | +2.71 (+63%) | **+1.21 (+28%)** | 2.2× |
| Q_fs(Xe+Sr) fission | 171.20 | −6.72 (−3.9%) | **−3.13 (−1.8%)** | 2× |
| Q_fs(Pd+Pd) fission | 190.12 | +2.20 (+1.15%) | **+0.35 (+0.18%)** | 6× |

Overall MAE over the 14 F3 α-decay nuclei: 1.06 MeV → **0.893 MeV**.

All values in this table are **sample-internal reproduction** (²³⁸U is inside the §3.3 fitting set); the out-of-sample numbers are those of §3.4.

### 6.4 Origin of the Residual Error

The residual Q_α deviation of +1.21 MeV is attributed to **deformation effects**: the uranium/thorium region is strongly deformed, while the implemented Myers–Swiatecki correction is spherical and single-constant, without deformation degrees of freedom or Strutinsky smoothing — an explicit work item for the next version. Quantified against the benchmark it cites: this work's global RMS of 3.07 MeV is about three times worse than FRDM2012 evaluated on the same nuclei set (RMS 0.99 MeV, §3.3 table; FRDM2012's own published rms is 0.56 MeV on its adjustment set [3]). The gap is the price of what this model buys with six global parameters, zero per-nucleus fitting, and a closed form whose differences are analytic; the trade is stated here rather than hidden. After the holdout of §3.4, the quotable precision for ²³⁸U is: **binding energy 0.30% (+5.44 MeV) and Pd+Pd fission channel 2.40% — out-of-sample numbers — not the in-sample 0.10% and 0.18%.**

---

## §7 Boundaries and Honest Disclosure

This paper deliberately retains complete failure records, because they are precise mappings of division-of-labor boundaries rather than defects:

1. **Shells do not emerge**: the k-distributions of magic nuclei are no more closed than those of neighbors; ¹⁶O is even among the most dangling-bond-rich on its isotone chain. Root cause: the structural engine relaxes toward "minimum total bond energy," whereas shell closure is a quantum single-particle gap — the two have different sources. Shell gaps therefore belong to the shell model, not the geometric layer (§4 already established ρ≈0).
2. **Dripline patches cannot be rescued**: the structural engine's S_n-sign judgment on the 12 dripline-patch nuclei is 0/12 correct, overbound by +31~+106 MeV. Root cause: the engine has no surface term, severely overbinding extremely neutron-rich nuclei. Dripline judgment depends on NUBASE patches, whose correctness does not come from the model.
3. **α-cluster ordering reverses in the mid-to-heavy region**: after multi-seed convergence, the α-cluster ordering gives Spearman r=−0.683 (p=0.042) — a **significant negative (reversed) correlation, not a null result**: the light region (A≤20) orders correctly, while the mid-to-heavy region (²⁴Mg–⁴⁰Ca) orders in reverse. Root cause: the purely geometric α-cluster criterion fails in large nuclei (accidentally close 2p2n groups are misidentified).
4. **Engine α binding energy 17% too low**: the engine's ⁴He binding energy is 23.38 MeV (28.30 MeV experimental), so a loose 2α configuration is energetically inferior to a compact cluster, and ⁸Be is crushed into a compact sphere. This is the manifestation of "a constant-slope bond-energy potential cannot simultaneously describe deuterium and the α particle" (an effective interaction lacking density dependence).
5. **Valve validation all four configurations failed**: replacing the decay-judgment layer's energy side from AME with the present closed form gives 33.9/57.2/66.7/66.7%, all below the 90% pre-registered red line. The "independent prediction" wording is accordingly abandoned.
6. **Half-life precision is order-of-magnitude**: log ft anchored at a single value, multivariate regression R²=0.138, half-lives only order-of-magnitude.

**Parameter-layering list (pre-registered)**: geometric engine C_REP=12, R_EQ=1.608 (four-light-nucleus calibration); closed-form four coefficients (cluster sum / grid calibration, zero per-nucleus fitting); macroscopic-microscopic six parameters (AME2020 global fit, not per-nucleus); B_α anchor 28.296 MeV (the only explicit experimental borrowing).

---

## §8 Conclusion and Collaboration Invitation

This paper establishes a three-layer division-of-labor framework: the liquid drop supplies macroscopic energy, the shell model supplies quantum gaps, and the four-bond saturated network supplies few-body spatial correlations — the three complementary **within their current validity domains** (liquid drop all A, shell correction A≥16, geometry A≤40; §4.1). The ²³⁸U application tests the framework's reproduction accuracy in the heavy region — in-sample binding energy 0.10%, and, under the frozen Z≥88 holdout of §3.4, **0.30%**, the latter being the first true out-of-sample prediction — and precisely locates the residual error (deformation effects).

**Code and data availability**: the complete reproducible package — the structure-generation engine, the global-refit and holdout scripts of §3.3–3.4, all per-nucleon datasets, and a `REPRODUCE.md` mapping every key number of this paper to its generating script and output file — is released under the **MIT license** at [repository URL, to be added at publication] and archived on Zenodo (DOI: 10.5281/zenodo.XXXXXXX, assigned upon publication). No access control applies; independent verification requires no correspondence with the authors.

---

## References

[1] M. Wang, W.J. Huang, F.G. Kondev, G. Audi, S. Naimi, Chin. Phys. C 45, 030001 (2021) (AME2020; source of all mass/Q-value/binding-energy data in this paper).

[2] F.G. Kondev, M. Wang, W.J. Huang, S. Naimi, G. Audi, Chin. Phys. C 45, 030002 (2021) (NUBASE2020; decay labels and half-lives).

[3] P. Möller, A.J. Sierk, T. Ichikawa, H. Sagawa, At. Data Nucl. Data Tables 109–110, 1 (2016) (FRDM2012; precision benchmark).

[4] I. Angeli, K.P. Marinova, At. Data Nucl. Data Tables 99, 69 (2013) (charge-radius experimental values).

[5] N.J. Stone, At. Data Nucl. Data Tables 111–112, 1 (2016) (quadrupole-moment experimental values).

[6] D.J. Wales, J.P.K. Doye, J. Phys. Chem. A 101, 5111 (1997) (basin-hopping algorithm).

[7] W. von Oertzen, M. Freer, Y. Kanada-En'yo, Phys. Rep. 432, 43 (2006) (α-clustering review).

[8] C.F. von Weizsäcker, Z. Phys. 96, 431 (1935); H.A. Bethe, R.F. Bacher, Rev. Mod. Phys. 8, 82 (1936) (classic liquid-drop model).

[9] R. Hofstadter, Rev. Mod. Phys. 28, 214 (1956) (electron scattering and nuclear charge distributions).

[10] W.D. Myers, W.J. Swiatecki, Nucl. Phys. 81, 1 (1966) (Myers–Swiatecki shell correction).
