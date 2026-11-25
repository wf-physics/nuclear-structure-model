# A Four-Bond Saturated Network Model of Nuclear Structure: From the Diamond Limit to Emergent Liquid-Drop Coefficients and a Difference-Formula Family

**Paper I (final manuscript, v2.2)** (revised 2026-11-24; merged Paper A + Paper B)
Rule system v5.0 (frozen). Throughout this paper, "zero fitting" means "zero *per-nucleus* fitting". All numbers correspond to the archives in `repo/results/` (mapping in Appendix D); one-command reproduction in Appendix E. The decay-verdict layer (99.0% on the full A≤40 set; four-protocol blind test) is reported separately in Paper C, whose net contribution of +1.5 percentage points is disclosed here. All results in the main text are **tests/reproductions** against known data; the only forecast in this paper is the pre-registered blind prediction of ³⁴Ca in Appendix F, kept strictly separate from the main text.

---

## Abstract

We propose and systematically test a phenomenological model of nuclear structure built on local bonding rules. The only nucleon–nucleon interaction is a windowed linear bond well (depth a calibration constant) plus Pauli short-range repulsion; each nucleon forms at most four bonds, and a bond requires mutual consent (both partners have a vacant slot). In the infinite nuclear-matter limit this minimal rule set converges to a **degenerate family of ideal tetrahedral four-coordinated networks** (cubic diamond is a representative member, not the unique ground state; intra-family energy differences ≤0.03 MeV/nucleon). The tetrahedral coordination unit is homologous to the α-cluster structure of light nuclei.

The model is validated quantitatively at three levels. (i) The nuclear-matter saturation point is *calibrated to* $\rho_0 = 0.160$ fm$^{-3}$ and $E/A = -15.8$ MeV (the calibration targets are the empirical values; meeting them is an acceptance criterion, not an independent prediction); the incompressibility $K = 547$ MeV possesses an analytic lower bound from window freezing, and the framework makes no prediction for the giant-monopole resonance energy. (ii) With zero per-nucleus fitting, the liquid-drop coefficients *emerge* from the model energetics: two sampling routes give $a_s/a_v = 1.18$ (empirical 1.14) and $a_{\mathrm{sym}} = 19.50$ MeV — the $a_{\mathrm{sym}}$ agreement is the *self-consistency of one energetics under two sampling schemes*, not a cross-validation of two independent models (mechanism in §4.3); the MAE over 2811 nuclides is 12.6 MeV. Linear differences of the same closed form produce four release-energy formulas — neutron separation $S_n$ (A≥40; 85 nuclides, MAE 1.58 MeV), $\beta^-$ decay energy $Q_{\beta^-}$ (10 nuclides, MAE 3.80 MeV), $\alpha$ decay energy $Q_\alpha$ (A≥180; 12 nuclides, MAE 1.06 MeV), and fission energy $Q_{fs}$ (5/5 in-band, MAE 3.83 MeV) — each with a domain of validity delimited by pre-registered criteria. (iii) Geometric observables (electric quadrupole sign, charge radius, α-cluster count) are quantified nucleus-by-nucleus over the full H→Ca set of 389 nuclides, among which the α-cluster count and the quadrupole sign are the only discriminators that the liquid-drop and shell models structurally lack.

The valve test, which substitutes the model's closed form for the AME on the energy side of the decay-verdict layer, fails under all four gate×calibration combinations (33.9/57.2/66.7/66.7% against a pre-registered 90% red line); we therefore abandon any "independent prediction" claim. The complete failure list, the pre-registration text, the falsification conditions, and the nucleus-by-nucleus data are published with this paper.

---

## §1 Introduction: From Empirical Facts to Minimal Assumptions

### 1.1 Six empirical facts

Any nuclear-structure model must confront six empirical facts that carry no model bias:

| # | Fact | Key numbers | Source |
|---|---|---|---|
| F1 | Nuclear-matter density saturation | $\rho_0 = 0.16 \pm 0.01$ fm$^{-3}$, almost independent of A | Electron scattering (Hofstadter series) [9] |
| F2 | Binding energy per nucleon ≈ constant | $B/A \approx 8$ MeV for A>20, variation <10% | AME2020 [1] |
| F3 | Short range + short-range hard core | Attraction range ~1–2 fm; strong repulsion for r<0.5 fm | NN phase-shift analyses |
| F4 | α clustering in light nuclei | 2p2n correlations in N=Z even-even 4n nuclei ($^8$Be, $^{12}$C, $^{16}$O, $^{20}$Ne) | α-cluster literature [7] |
| F5 | 4 is magic in light nuclei | $^4$He doubly magic, $B = 28.30$ MeV far above neighbours | AME2020 [1] |
| F6 | Spin-isospin quartet | 4 intrinsic states per nucleon (p↑ p↓ n↑ n↓) | elementary |

### 1.2 From facts to minimal assumptions: the mapping table

The central methodological stance of this paper is that **every model rule must be traceable to the empirical fact it answers**; otherwise it is a superfluous assumption. The complete mapping:

| Fact | Constraint it forces | Minimal assumption chosen | Why not a richer one |
|---|---|---|---|
| F3 short range + hard core | The potential must be nonzero in a finite window with short-range repulsion | Windowed linear well $V(r) = -s(r_2 - r)$ for $r\in[r_1,r_2]$ + Pauli repulsion | Morse/Reid forms introduce 5–10 shape parameters and forfeit "minimal"; the linear well is the fewest-parameter form satisfying F3 |
| F1 density saturation | A saturation mechanism must exist, else collapse | **No explicit density dependence**: $k_{\max}=4$ + mutual consent caps the bond count as neighbours fill up, giving saturation | An explicit $F(\rho)$ term bakes the answer in; "emergent saturation" is what *explains* F1 |
| F2 constant B/A | Bond energy per nucleon must approach a constant with A | Same: coordination ≤4 → $B/A \to (k_{\max}/2)\epsilon_{\text{bond}}$ (plus surface correction) | Same |
| F5+F6 quartet | Bond-count cap matches intrinsic-state count | $k_{\max} = 4$ | $k_{\max}=3$ or 5 have testable consequences (§3.3, falsification 2); 4 also matches the α-cluster identity (§3.3) |
| F4 α clusters | Tetrahedral units should arise naturally | Mutual-consent bonding | One-sided bonding would allow >4 coordination and destroy tetrahedral units |

### 1.3 First-principles derivation of the decay trigger: the monotonicity impasse and the global Q value

Beyond the bonding rules, the model must answer "what switches decay on?". This part is derived not from assumptions but **from counterexamples**:

**The minimal counterexample (monotonicity impasse)**. $^2$H is stable while $^3$H β-decays ($t_{1/2}=12.3$ yr). In $^3$H every neutron sits in a *stronger* local proton field (local Coulomb field, coordination environment) than the single neutron in $^2$H. If decay were triggered by "local stimulus above threshold" (a monotonic criterion), the $^3$H neutron should be *more* stable than $^2$H's — contradicting experiment. If instead one uses a non-monotonic "window" criterion, the windows cannot be self-consistently calibrated across nuclides. We actually tested five local trigger variants during model evolution (scalar-Coulomb escape, bond-environment stimulus, quark wave-function overlap, internal Coulomb threshold, nn-pairing enhancement) — all died on this counterexample (control-evaluation loss floor 36; the nn-pairing variant, although it fitted well, forced a negative nn coefficient and misjudged 5/5 heavy stable nuclei). **The data force the conclusion: any "local stimulus → decay rate" causal chain, monotonic or windowed, is incompatible with experiment.**

**The surviving criterion (global Q-value mass balance)**. Whether decay can proceed is determined by the initial–final mass difference: $Q \le 0$ strictly forbids the channel, $Q>0$ opens it. This criterion is global, not local, and dissolves the impasse by construction: $^2$H has $Q_{\beta^-}<0$ (forbidden) while $^3$H has $Q_{\beta^-}=+18.6$ keV (open) — the difference comes from the global nuclear mass difference, not local geometry. The criterion has zero free parameters and yields 99.7% stability verdicts over all A≤40 nuclides (a single exception) with **zero stable-nucleus misjudgment**; it is the only one of seven tested variants that satisfies the hard constraint "stable-nucleus misjudgment rate strictly zero" (details in Paper C).

This derivation is itself a falsifiable statement: if a nuclide with $Q>0$ were ever found to not decay at all (and not via a forbidden transition), the verdict framework is falsified.

### 1.4 What the model does NOT take as input (against "baked-in answer" charges)

| Model input | Model does NOT take as input |
|---|---|
| Windowed linear well form (depth a calibration constant) | Density-dependent interaction terms |
| Pauli short-range repulsion (two shape parameters, fixed by geometry) | Shell structure / single-particle orbits |
| $k_{\max} = 4$ | Spin–orbit coupling |
| Mutual-consent rule | Pairing interaction (odd-even term) |
| Q-value trigger (zero parameters) | α pre-formation probability |
| Phase-space + WKB rate layer (6 verdict parameters, fitted on training set only) | Any per-nucleus experimental data (energy layer closed form) |
| | Any mass-table values (energy layer closed form) |

### 1.5 Falsifiability statement

The following three experimental/computational results would overturn the model:

1. **Self-contradiction of rules**: if some stable nuclide's bond network admits no legal configuration under $k_{\max}=4$ (structure search finds no configuration with energy below the free-nucleon sum);
2. **Geometric-architecture failure**: if α-cluster recognition on the N=Z even-even 4n chain falls significantly below a random baseline (§5.3 reports 8/10 against a random control);
3. **Large-A limit failure**: if the emergent $a_s/a_v$ deviates from the empirical liquid-drop value by >30% (measured: +3.5%).

### 1.6 Relation to existing theory

We do not claim to replace the shell model or the liquid-drop model. The liquid-drop terms are the inevitable A→∞ limit of *any* short-range saturating two-body force (a theorem; derivation at the start of §4), so "the model automatically produces liquid-drop terms with matching values" is an **acceptance criterion**, not a bonus. In numerical precision the model currently trails the fitted liquid drop by a factor of 3–6 (§4.4, stated honestly). The only differential discriminators are the **α substructure and the geometric observables** — the shell model describes them by single-particle orbits, the liquid-drop model by a continuous fluid, and neither possesses an α-clustering mechanism structurally (§3.3, §5.3). Boundaries such as half-life precision and drip-line verdicts are covered in §10–11 and Paper C.

### 1.7 Main results table (all numbers in the text refer to this table)

| # | Result | Value | Sample domain | Source |
|---|---|---|---|---|
| R1 | Saturation point | $\rho_0=0.160$ fm$^{-3}$, $E/A=-15.8$ MeV (calibrated to) | infinite limit | §3.2 |
| R2 | Incompressibility | $K=547$ MeV (analytic lower bound exists; no GMR prediction) | infinite limit | §3.2/§7.4 |
| R3 | Route-1 closed form | $a_v=16.060$, $a_s=19.010$, $a_c=0.7560$, $a_{\mathrm{sym}}=19.50$ | 2811 nuclides, MAE 12.6 | §4.1 |
| R4 | Route-2 closed form | $a_v=13.851$, $a_s=1.592$, $a_c=0.7265$, $a_{\mathrm{sym}}=19.50$ | 99 grid points, RMS 128.7 (shell noise) | §4.2 |
| R5 | F1 $S_n$ | MAE 1.58 MeV (59 pass/25 half/1 fail) | 85 nuclides, A≥40 dominated | §8.1 |
| R6 | F2 $Q_{\beta^-}$ | MAE 3.80 MeV (Lane convention); residual fully absorbed by shell-correction difference | 10 test nuclides; 32 regression nuclides | §8.2 |
| R7 | F3 $Q_\alpha$ | MAE 1.32 MeV (3 pass/9 half/2 fail) | 14 nuclides (12 with A≥180 + 2 light controls) | §8.3 |
| R8 | F4 $Q_{fs}$ | MAE 3.83 MeV (5/5 in-band) | 5 fission channels | §8.4 |
| R9 | Valve test | 33.93/57.19/66.74/66.74%, all <90% | A≤40 blind protocol | §10 |
| R10 | Decay-verdict layer | full-set 99.0% (311/314); blind mean 97.4±1.1%; net +1.5 pp | A≤40 | Paper C |
| R11 | Geometric validation | B 19.8%/55.0% (not in the scorecard); $r_{ch}$ 66.0% (n=97); $Q_s$ 69.2% incl. half-pass (n=65); α 8/10 | H→Ca 389 nuclides | §5 |

---

## §2 The v5.0 Rule System

### 2.1 Rule 1 (bonding and saturation)

The only nucleon–nucleon interaction is the following pair potential (the rule text was frozen before any validation began; zero modification during validation):

- **Attractive well**: $V_{\mathrm{att}}(r) = -s(r_2 - r)$, active in the window $r \in [r_1, r_2]$, zero outside; numerically $s = 5$ MeV/fm, $r_1 = 0.8$ fm, $r_2 = 2.2$ fm (the depth $s$ is a formula-layer calibration constant, §2.3);
- **Pauli repulsion**: $V_{\mathrm{rep}}(r) = C\,\max(0, r_{eq} - r)^2$, added to (not replacing) the attraction; the two parameters are constrained by the pairing equilibrium $r_{\min} = r_{eq} - s/(2C) \in [1.15, 1.6]$ fm (a mathematical constraint, not a fit degree of freedom);
- **Hard core**: strong quadratic repulsion for $r < r_1$;
- **Bond-count saturation**: at most $k_{\max} = 4$ bonds per nucleon;
- **Mutual consent**: a bond requires both partners to have a vacant slot (bidirectional occupancy < $k_{\max}$).

Correspondence of the three sub-rules to facts: attractive window ← F3 (short range); saturation + consent ← F1/F2 (density and binding saturation); $k_{\max}=4$ ← F5/F6 (quartet).

### 2.2 Rule 2 (decay trigger)

A three-layer architecture (the result of the §1.3 derivation):

1. **Trigger layer (zero free parameters)**: all decay channels (β⁻/β⁺/EC/n/p/α) are opened or forbidden by the mass-balance Q value; Q is computed directly from AME2020 mass excesses;
2. **Rate layer**: β channels via the exact phase-space integral $f(Z_d, Q)$ (Fermi function $F(Z_d,E) = 2\pi\eta/(1-e^{-2\pi\eta})$, with the sign convention $Z_d = Z+1$ for β⁻ and $Z_d = -(Z-1)$ for β⁺ — missing this convention necessarily breaks EC-dominated nuclides) competing with WKB barrier penetration for particle channels; the fastest is the dominant mode;
3. **Matrix-element layer**: the log $ft$ baseline is anchored to the $^3$H literature value 3.049. **Reason for anchoring**: $^3$H is the only nuclide satisfying simultaneously "allowed transition + precise experiment + belongs to the calibration set"; once anchored, this parameter is never re-tuned over all A≤40. The cost is order-of-magnitude half-life precision (multivariate regression $R^2 = 0.138$ on log $ft$; strongly hindered transitions can be off by >3 orders of magnitude), listed honestly in §11.

### 2.3 Freeze discipline and pre-registration

The rule text was frozen (v5.0) before any validation began and was never modified during validation. **Pre-registration clause, quoted verbatim** (freeze-audit document §〇·二):

> "Conclusion: for the energy-formula layer to reach the target, **the well-depth constant must be rescaled as a 'formula-layer calibration constant' (the linear-well structure of the rule is unchanged)**, or one must introduce density-dependent bond-length contraction. The former changes only a constant and requires re-running the engine for geometric validation…… the latter touches Rule 2 §1.1 text and, per the freeze order, **requires explicit user approval**."

That is, the *functional form* of the rules is frozen; only the numerical scale is set by nuclear-matter calibration (the ×2.61 rescaling of §4.1 Route 1 follows this clause). Energy results that depend on the well depth (§4 Route 1, §7) use the rescaled constant; geometric validation (§5) uses the independent engine calibration ($C_{\mathrm{REP}}=12$, $r_{eq}=1.608$) without touching the well depth; the two are labelled separately.

### 2.4 Structure-search algorithm and stochasticity

Structure generation uses diamond-tetrahedral-lattice initialization + L-BFGS local minimization + basin-hopping global search (the scipy implementation of the Wales–Doye algorithm [6]). Convergence: energy-gradient norm <10⁻¹⁰ and two successive minimizations differing by <10⁻¹².

**Stochasticity treatment (multi-seed measurement, 6 nuclei × 5 seeds, `seed_dispersion.json`)**: energy observables are seed-robust (B/A scatter mean 3.4%, max 7.7%); **geometric shape observables are seed-sensitive** — $Q_0$ std up to 1154 ($^{20}$Ne), α-cluster count in the A≥20 region scatters over 1–5 across seeds, $r_{ch}$ std up to 2.9 fm; **the α-cluster count in the A≤14 region is fully seed-stable**. The main archive uses a fixed single-seed convention (seed=13, pre-registered before the dispersion measurement); geometric conclusions in the A≥20 region are therefore downgraded to qualitative evidence in §11. The six verdict-layer parameters ($\nu_{\mathrm{try}}$, $S_\alpha$, EC density factor, etc.) are fitted on the training set only within the blind-test protocol (Paper C).

---

## §3 The Diamond Limit

### 3.1 From rules to lattice: the degenerate-family statement and the exclusivity argument

Mutual consent together with $k_{\max}=4$ caps the network coordination at 4. In the infinite nuclear-matter limit the problem becomes energy minimization of a four-coordinated network under the given nearest-neighbour bond potential.

**Degenerate-family statement**. All "ideal tetrahedral" four-coordinated networks — equal bond lengths, all bond angles 109.47° — have second-neighbour distance $2\,r_{nn}\sin 54.74^\circ$. At $\rho_0$, $r_{nn}=1.595$ fm, so the second-neighbour distance is always 2.605 fm, beyond the window upper edge 2.2 fm; hence the energies of such networks are **strictly degenerate** under this potential: the family includes cubic diamond, hexagonal diamond (lonsdaleite), and amorphous tetrahedral networks (a-C type). **Cubic diamond is the representative member, not the unique ground state.**

**Exclusivity argument (why no lower-energy four-coordinated network exists)**. Fixed $\rho_0$ → fixed bond length → bond energy is determined only by bond count and length; the count is already at the cap $k_{\max}=4$, so ideal tetrahedral networks are strictly degenerate. Any non-ideal tetrahedral network that pulls second neighbours into the window (extra attraction) requires a bond angle <87.2° (a 22.3° deviation from ideal), at the cost of angular strain over all bonds; for an amorphous tetrahedral network with bond-angle spread σ≈10°, the fraction of bond pairs falling inside the window is ~1.3%, amounting to ~0.026 MeV/nucleon — not a lower-energy state. The infinite limit of this model is therefore the **degenerate family of ideal tetrahedral networks**; all quantitative diamond-configuration conclusions in this paper hold equally for any family member, and this indistinguishability is stated honestly.

**Why FCC necessarily fails, quantitatively**. FCC close packing has nearest-neighbour distance 2.07 fm at $\rho_0$, hugging the outer window edge 2.2 fm — the bond is open over only 9% of width, making the bond depth acutely density-sensitive (this is why all FCC solutions in the §3.2 table have $K\ge846$), and 12-fold coordination under $k_{\max}=4$ produces a large number of unsatisfied dangling bonds. Diamond at $\rho_0=0.16$ has nearest-neighbour distance 1.595 fm, comfortably inside the window (opening-width factor 2.6).

### 3.2 Systematic solution and pathological controls

Four independent parameterizations of the frozen potential (multiplicative density function $F(\rho)$, Morse control, density-dependent Pauli core, direct solution of the analytic equations) and two lattices are solved systematically (280-parameter grid, DE global optimization):

| Scheme | Valid solutions | $\rho_0$ (fm⁻³) | $E/A$ (MeV) | $K$ (MeV) | Pathology |
|---|---|---|---|---|---|
| FCC + $F(\rho)$ (7-param DE) | 1 | 0.30 | −15.0 | 846 | $\rho_0$ +88% off; K 2.5× empirical |
| FCC + density-dependent Pauli (5-param DE) | 1 | 0.30 | −11.4 | 1446 | same + E/A −28% off |
| FCC + direct analytic (280 grid) | 5 | — | — | ≥6.7×10⁴ | all K>350 (criterion); non-physical compressibility |
| **Diamond four-bond network + direct analytic** | **120** | **0.160** | **−15.8** | **547** | K misses 240±30 (see below) |

**Pathology criterion**: any of "$\rho_0$ off >25%", "E/A off >15%", "K>350" marks a pathological solution; the diamond solution misses only the K item. Saturation points and binding energies of the four schemes are shown in **Figure 1**.

**The analytic wall (K lower-bound theorem)**. $K=547$ misses the empirical $240\pm30$ because of an analytic constraint under window freezing, not a parameter problem: the second derivative of the attractive term gives $K_{\mathrm{att}} = 4\,k_n s\, c\,\rho_0^{-1/3}$ ($k_n=2$ bonds per nucleon, $c$ the diamond geometry constant; numerically 63.8 MeV at $s=5$, 166.5 MeV at the rescaled depth — analytic and numerical values agree exactly; **step-by-step derivation in §7.4**). All four parameterizations and both lattices hit the same wall — this is a **consequence of a formal theorem**, not a fitting failure. K enters only GMR-type compression observables; none of this paper's target observables depends on K, and the framework makes no GMR prediction (§11 limitation 1).

### 3.3 α-cluster identity

The diamond tetrahedral coordination unit (centre + 4 vertices, $r_{nn}=1.59$ fm) is the geometry of the 4-nucleon α cluster — **this is the homology point between the "quartet $k_{\max}=4$ + consent" assumption and the empirical fact of α clustering in light nuclei (F4)**: the same rules give a tetrahedral network in the infinite limit and detectable tetrahedral sub-configurations in finite N=Z nuclei. α-cluster detection is performed independently of nuclear-matter calibration (§5.3: 8/10 hits on the N=Z even-even 4n chain, fully seed-stable in A≤14), constituting a self-consistency check of the **configurational homology among nuclear-matter limit, finite-nucleus geometry, and α clustering**. Honest boundary: this is **configurational homology**, not "proof of α pre-formation" — the latter requires dynamical and spectroscopic evidence, beyond this framework.

---

## §4 Emergent Liquid-Drop Coefficients: Two Independent Routes

The four-parameter closed form $B(A,Z) = a_v A - a_s A^{2/3} - a_c Z^2/A^{1/3} - a_{\mathrm{sym}}(N-Z)^2/A$ has its coefficients derived from the model energetics by two mutually independent routes, both with zero per-nucleus fitting.

### 4.1 Route 1: rescaled cluster summation

**Analytic upper-bound derivation**: under four-bond saturation the bond count per nucleon → 2 (as A→∞, bonds = 2A); a single bond at $\rho_0$ has energy $-s(2.2-1.595) = -3.03$ MeV (original depth $s=5$), so the analytic upper bound is $E/A = -6.05$ MeV. Per the pre-registration clause (§2.3) the well depth is rescaled as a formula-layer calibration constant (×2.61) to $E/A(\rho_0) = -15.8$. Diamond spherical clusters (N=8–600) are summed and linearly regressed ($E/A$ vs $A^{-1/3}$) to extract **$a_v=16.06$, $a_s=19.01$, $a_s/a_v=1.18$** (empirical 1.14, +3.5% off). The closed form over the full AME2020 table A=40–240 (2811 nuclides): **MAE 12.6 MeV, median relative deviation −0.3%, all within a 15% band** (per-nucleus deviations with three mass-region colouring in **Figure 2**).

### 4.2 Route 2: model-grid calibration

On the discrete four-bond network with the full Pauli-core repulsion, two-component Fermi kinetic energy, and Coulomb (Metropolis optimal proton assignment), 99 grid points (A=40–240, N−Z∈[−24,24]) are sampled and least-squares fitted to extract **$a_v=13.851$, $a_s=1.592$, $a_c=0.7265$, $a_{\mathrm{sym}}=19.50$** — four coefficients, zero experimental borrowing. The fit RMS of 128.7 MeV arises from discrete-lattice shell noise (σ≈100–130 MeV, uniform across all |N−Z| bins) — this simultaneously explains why $a_s$ is unresolvable on this route (the surface signal ~$19A^{2/3}$ is drowned by shell noise in the residual).

### 4.3 Cross-validation and differences

- **$a_{\mathrm{sym}}=19.50$ agrees exactly between the two routes — mechanism and honest grading**: the two routes are not two mutually independent models — Route 1 (cluster summation) and Route 2 (grid calibration) sample **the same bond-type-selectivity mechanism of the same model energetics** (§7.7), with no experimental symmetry-energy data at the input. The exact agreement therefore **cannot be called "cross-validation of two independent routes"; it is the self-consistency of two sampling routes under the same anchor (the same model energetics)**. What it demonstrates is robustness of the number against the sampling scheme (analytic summation vs Metropolis grid + least squares), not evidence external to the model;
- $a_c$ agrees universally between routes but is **uninformative** ($0.72 \approx 0.6e^2/1.2$, an arithmetic result any Coulomb model gives) — kept as an honesty flag;
- $a_v$ agrees to within the kinetic-energy interpretation (16.06 vs 13.85+kinetic);
- $a_s$ is explicit only in Route 1 (19.01): Route 2's shell noise drowns the surface signal.

### 4.4 Precision comparison (disclosed in the same table)

| | This closed form | Fitted liquid drop | FRDM2012 |
|---|---|---|---|
| Metric | MAE 12.6 / RMSE ~15 MeV | RMSE 2–4 MeV | RMS 0.7 MeV |
| Per-nucleus fit parameters | **0** | 4–6 | ~30 |
| Coefficient source | geometric summation of the potential (emergent) | fit to experiment | fit to experiment |

**Same-set footnote**: FRDM2012's self-reported RMS 0.57 MeV is on its **own fit set (2149 nuclides)**; on the intersection with this paper's closed-form domain (A=40–240), n=2910, FRDM measures **MAE 0.45 / RMS 0.63 MeV** (archive `results/frdm_same_set.json`, script `src/frdm_same_set.py`). Same-set comparison: this paper MAE 12.6 / RMS ~15 vs FRDM 0.45 / 0.63 — **a 24–28× gap**, with a clear origin: this closed form has no shell-correction term, no deformation degrees of freedom, and no Strutinsky smoothing; this is an explicit work item for the next version. The price buys **zero per-nucleus fitting, coefficients fully emergent from the geometric summation of the potential, and an analytically differentiable closed form** (the four release-energy formulas of §6–8 all rely on this property). The 0.7 MeV in the table is the common literature shorthand; the same-set numbers are the precise comparison.

Achieving liquid-drop-grade approximation under zero per-nucleus fitting, with a 3–6× gap (against the fitted liquid drop) stated honestly; this paper does not claim liquid-drop or macroscopic-microscopic precision.

---

## §5 Geometric-Observable Validation (H→Ca, Full Set of 389 Nuclides)

All observables below follow a five-part template: verdict band (pre-registered text) → sample → three-tier counts with per-nucleus list (Appendix A) → structured attribution of failures → seed dispersion. Verdict bands and sample domains all correspond to the single archive `results/big_validation_table.csv` (one generation, multiple columns, no engine-convention switching).

### 5.1 Binding energy B (engine acceptance quantity, not in the scorecard)

- **Band** (pre-registered text): relative deviation $|B_{\mathrm{model}} - B_{\mathrm{exp}}|/B_{\mathrm{exp}} \le 15\%$ = pass, $\le 30\%$ = half-pass, >30% = fail;
- **Sample**: full H→Ca set of 389 nuclides (structure engine relaxed nucleus-by-nucleus; the $B_v$ column of the same archive);
- **Result**: **77/389 pass (19.8%), 137 half-pass, 175 fail; including half-pass 55.0% (214/389)** — this is the honest energy-side score of this engine, identical to the early-work archive (no re-measurement, no convention switch);
- **Structured attribution**: the structure engine has no surface term and no density dependence, so mid-to-heavy nuclei are systematically over-bound — precisely the reason the energy-formula layer (§6–7, explicit surface-term closed form) exists; the failure is *connected to* subsequent work, not an isolated defect. Distinction: this item is **engine-energy** acceptance (an honestly recorded non-pass), whereas §4's **closed form** (MAE 12.6) is the formula-layer result — the two are different sources and must not be conflated;
- **Statement**: binding energy is an engine acceptance quantity and is **not counted in the model scorecard** (scorecard = §5.2–5.3 + §8).

### 5.2 Charge radius $r_{ch}$

- **Band** (pre-registered text): model/experiment ratio $\in [0.7, 1.3]$ = pass, outside = fail (two tiers);
- **Sample**: **n=97** (those of the 389 nuclides with IAEA charge-radius measurements [4]; the other 292 have no experimental value and are not judged);
- **Result**: **64/97 = 66.0% in-band**, median ratio 0.794 (systematic −21% bias);
- **Attribution**: missing zero-point-motion correction — the 0.794 systematic undershoot is structural (zero-point motion inflates the charge distribution), a model-layer defect (§11 limitation 7), not a validation oversight; per-nucleus list in Appendix A.
- **Scaling diagnostic (separating a "scaling problem" from a "distribution problem")**: rescaling all model radii by the median ratio 0.794 as a global scale factor raises the fraction within the ±15% band from 30.5% to 52.6% (re-checked on the archived per-nucleus data, n=95, `results/radius_scaling_result.json`; differs from the main-text n=97 by 2 boundary nuclides, whose archived median 0.815 vs 0.794 comes from the same convention difference), with the remaining scatter concentrated in the light region. Restated conclusion: **the model's radius *shape* is right and the *scale* is off by an approximately constant factor** — absolute radii are unusable (scale off by ~21%), while relative trends after rescaling are usable.

### 5.3 Electric quadrupole moment and α-cluster count

**Electric quadrupole** ($Q_s = Q_0/(1-\lambda)$, classical-rotor conversion). Two sample domains, two bands, each reported honestly:

- **Main convention (pre-registered band text)**: "sign consistent and magnitude within 5× = pass; sign consistent but magnitude out of band = half-pass; sign inconsistent = fail". **Sample n=65** (those of the 389 with IAEA quadrupole measurements [5]): **34 pass + 11 half-pass = 69.2% including half-pass (45/65); the sign-consistency rate is also 45/65 = 69.2%** (all 20 "fails" are sign-inconsistent — the band structure makes every pass/half-pass necessarily sign-consistent);
- **Sub-sample control**: intersection with the Stone 2016 recommended-value table, 13 nuclides, under the single "sign-consistent" band: **9/13** (the recommended-value sub-sample is stricter: it drops low-precision IAEA measurements); the two numbers are reads of the *same* computation under different experimental compilations and bands, not two independent validations;
- Rank correlation of |Q|: $\rho = 0.005$ ($p = 0.986$, n=65) — **sign right, magnitude wrong**: magnitude requires a quantitative theory of moment of inertia and deformation, which this framework currently supports only at the sign level;
- Failures concentrate in strongly deformed regions (configurational distortion is seed-sensitive, §11 limitation 6).

**α-cluster count** (2p2n tetrahedral sub-configurations ≥ 50% of ideal cluster count = pass):

- Band (pre-registered text): number of recognized 2p2n tetrahedral sub-configurations ≥ $\lfloor A/4\rfloor \times 50\%$ = pass;
- Sample: N=Z even-even 4n chain, 10 nuclides ($^4$He–$^{40}$Ca);
- Result **8/10**; the $^4$He–$^{28}$Si chain passes throughout under the pre-registered band (light region A≤20 in the correct direction); **in the mid-to-heavy region (A≥32) the identification ratio shows a reversed ordering trend** — Spearman rank correlation of the identified/ideal cluster ratio against A is r=−0.603, p=0.065 (n=10, marginally significant, archive `results/alpha_cluster_spearman.json`), with the shortfall concentrated at ³²S–⁴⁰Ca (Coulomb stretching + under-bound distortion) — a graded trend, not random failure;
- Seed dispersion: **A≤14 fully seed-stable; A≥20 count scatters over 1–5, qualitative evidence only** (§2.4);
- This item and the quadrupole sign are the differential discriminators the liquid-drop/shell models structurally lack.

### 5.4 Honest summary of geometric validation

| Observable | Criterion (pre-registered) | n | Pass | Note |
|---|---|---|---|---|
| Binding energy | ≤15% pass / ≤30% half | 389 | 19.8% pass / 55.0% incl. half | engine acceptance, not in scorecard |
| Charge radius | ratio ∈ [0.7, 1.3] | 97 | 66.0%, median 0.794 | systematic undershoot is structural |
| Quadrupole | sign + magnitude ≤5× (main) | 65 | 69.2% incl. half (69.2% sign-consistent) | magnitude wrong (rank corr. ~0) |
| Quadrupole (sub-sample) | sign-consistent | 13 | 9/13 | Stone 2016 recommended intersection |
| α cluster | ≥50% sub-configurations | 10 | 8/10 | A≤14 seed-stable; A≥20 qualitative |

Pass rates by mass region for the four observables are shown in **Figure 3** (black dashed line = overall pass rate, "pass" tier only, excluding "half-pass").

---

## §6 From Closed Form to Differences: Why Linear Differences (Bridge Section)

§4 established the closed form $B(A,Z)$; §8 onward splits the same closed form into four release-energy formulas. The two are **two segments of one chain**; this section states the connection and its cost.

1. **Release energy = binding-energy difference — a physical definition, not a computational trick.** $S_n \equiv B(A,Z) - B(A-1,Z)$; $Q_\beta \equiv B(A,Z+1) - B(A,Z)$; $Q_\alpha \equiv B(A-4,Z-2) + B_\alpha - B(A,Z)$; $Q_{fs} \equiv B(f_1) + B(f_2) - B(\text{parent})$ — all four are linear combinations of initial and final binding energies, so *any* energy theory must pass through this channel to confront decay/fission data.

2. **Numerical value of differencing: common-term error cancellation.** Differencing the closed form term by term, the volume term $a_v A$ cancels directly in every "one-nucleon-less" difference (F1) and is strictly zero in A-invariant differences (F2); only the $a_c$, $a_s$, $a_{\mathrm{sym}}$ differences remain. The closed form's 12.6 MeV-level absolute-energy error contracts to the MeV level after differencing (1.58 in F1) — differencing is an *amplifier of closed-form precision*.

3. **Cost of differencing: what fails to cancel becomes the residual.** Structure the closed form lacks — shells, odd-even pairing, light-nucleus clustering — is *not* cancelled by differencing and is exposed as-is: this is the common root of F2's −3 to −9 MeV shell-effect negative bias (§8.2), F1's failure in the A<40 extrapolation region, and the §10 valve-test failure. Differencing cannot create physics the closed form does not have.

4. **Bridge.** §4's closed form is the *acceptance ruler* of the model energetics (coefficients emergent from the model, form the inevitable liquid-drop limit); §8 onward uses the same closed form as *building blocks*, split into four observables validated one by one. Which route feeds which observable is decided by the sensitivity structure of each difference (§9).

---

## §7 Foundations: From Potential to Closed Form (Derivation Chain)

### 7.1 From $V(r)$ to the bond-energy sum

The total model binding energy is the sum of the pair potential over the bond network (§2.1 Rule 1). For a given configuration $\{r_{ij}\}$:

> $B_{\mathrm{model}} = \tfrac{1}{2}\sum_{i\ne j} W_{ij}\,V(r_{ij}) + V_{\mathrm{Pauli}} + T_{\mathrm{Fermi}} + E_{\mathrm{Coul}}$

Three structural factors: the **½ factor** from double-counting in the pair sum (each bond counted twice); the **window cutoff** making pairs with $r>2.2$ fm contribute zero, with the bond count capped by mutual consent (§2.1); the **Pauli core** adding to the attraction in the $r<r_{eq}$ region to give the equilibrium bond length (pairing-equilibrium constraint $r_{\min}\in[1.15,1.6]$ fm).

### 7.2 Uniform limit: lattice sum → energy density

In the infinite uniform limit, bonds per nucleon → $k_{\max}/2 = 2$, and the nearest-neighbour distance is $d(\rho) = c\,\rho^{-1/3}$ (diamond geometry constant **$c = 0.866$ (dimensionless)**, obtained from $d = 1.595$ fm at $\rho_0 = 0.16$ fm⁻³: $c = d_0 \rho_0^{1/3} = 1.595 \times 0.543$). The energy density splits into the attractive sum and the Fermi kinetic energy:

> $e(\rho) = -k_n s\,(2.2 - c\,\rho^{-1/3}) + T(\rho)$,  $k_n = 2$ (bonds per nucleon)

with $T(\rho)$ the two-component Fermi-gas expression (numerically calibrated).

### 7.3 Saturation condition

Two calibration equations: $e(\rho_0) = -15.8$ MeV (the attractive term alone gives −6.05 MeV/A at the original depth; the ×2.61 rescaling of §4.1) and $\mathrm{d}e/\mathrm{d}\rho|_{\rho_0} = 0$ (attractive negative slope balanced by Pauli positive slope). $\rho_0=0.160$ and $E/A=-15.8$ are **calibration targets** (empirical values); meeting them is an acceptance criterion. The $e(\rho)$ curve is in A-Fig. 1.

### 7.4 Analytic constraint on the incompressibility (the structural origin of K)

Analytic second derivative of the attractive term (every step traceable; numerical verification agrees):

> $E_{\mathrm{att}}/A = -k_n s\,(2.2 - c\,\rho^{-1/3})$
> $\mathrm{d}^2(E_{\mathrm{att}}/A)/\mathrm{d}\rho^2 = k_n s\, c\, (4/9)\,\rho^{-7/3}$
> **$K_{\mathrm{att}} \equiv 9\rho_0^2 E''|_{\rho_0} = 4\,k_n s\, c\,\rho_0^{-1/3}$**

Numerically (diamond $c=0.866$, $\rho_0^{-1/3}=1.843$): original depth $s=5$ → **$K_{\mathrm{att}} = 63.8$ MeV**; rescaled depth $s=13.05$ → **$K_{\mathrm{att}} = 166.5$ MeV**. The analytic and numerical gradients agree exactly. The full-model $K=547$ MeV is the numerical solution of the attractive convexity + Pauli-core shape + Coulomb (the `nm_formula_v4` direct solution gives 546.993). **Structural inference**: lowering K requires halving the attractive convexity, equivalent to burying the equilibrium bond length deeper into the window (outer edge >2.6 fm) — which in turn breaks the surface-term magnitude (the K–surface trade-off, §3.2). K enters only GMR-type compression observables; this framework makes no GMR prediction.

> **Writing note**: an earlier working record contained a numerical approximation of the form "K ≳ 0.23×(71.6k−192)"; at the present paper-grade review its coefficients could not be traced step-by-step, so it has been rewritten as the traceable derivation above, and the old form is void (honest-correction record; the §3.2 analytic-wall statement has been synchronized with this form).

### 7.5 Cluster summation → surface term

Spherical diamond clusters (N=8–600) are summed and the energy decomposed into A and $A^{2/3}$: the volume term $a_v$ from the linear-regression intercept, the surface term $a_s$ from the $A^{-1/3}$ slope. The "surface too cheap" bias of the original depth is the un-rescaled well depth — after ×2.61 rescaling, $a_s = 19.01$ (+5.6% vs empirical 18), $a_s/a_v = 1.18$ (+3.5% vs 1.14). MAE 12.6 MeV over the full AME2020 A=40–240 table (2811 nuclides), median relative deviation −0.3%, all within 15%.

### 7.6 Coulomb term

The closed-form Coulomb term uses the uniform-sphere approximation $a_c Z^2/A^{1/3}$ ($=0.6 e^2 Z^2/R$, $R = r_0 A^{1/3}$). **Discrete lattice summation overestimates Coulomb** (point charges 0.93 vs continuous 0.72; the point-charge pair sum ignores the smearing of charge within a nucleon volume); the fact that Route 2's fitted $a_c = 0.7265$ matches the continuous value shows the least squares absorbs this systematic difference. $a_c$ agrees between routes but is **uninformative** — an arithmetic result of Coulomb's law, not model evidence (an honesty flag).

### 7.7 Symmetry energy

In the symmetry term $a_{\mathrm{sym}}(N-Z)^2/A$, the purely geometric contribution of bond-type selectivity (pp/nn/np weight differences) is at most ~11–14 MeV, short of the empirical 23 MeV — the deficit is taken up by the Lane term (an isovector correction with a physical source). Routes 1 and 2 sample the **same number $a_{\mathrm{sym}}=19.50$** (§4.3) — neither route has experimental symmetry-energy data at its input, so the agreement is the **self-consistency of one model energetics (the same bond-selectivity curve) under two sampling schemes**, not an independent cross-validation; the −3.5 deficit is one of F2's known sticking points (§8.2).

### 7.8 Closed-form summary

| Route | $a_v$ | $a_s$ | $a_c$ | $a_{\mathrm{sym}}$ | Borrowing | Use |
|---|---|---|---|---|---|---|
| 1 (rescaled cluster sum) | 16.060 | 19.010 | 0.7560 | 19.50 | well-depth rescaling only (pre-registered clause) | surface-sensitive differences (F4) |
| 2 (model grid) | 13.851 | 1.592 | 0.7265 | 19.50 | zero experimental borrowing | difference-response (F1/F2/F3) |

Physical basis for per-observable route selection in §9.

---

## §8 The Four Release-Energy Formulas F1–F4

Each section follows a six-part template: definition → derivation → criterion (pre-registered) → score table (by mass region, n/pass/half/fail/MAE) → failure mode with quantitative attribution → one-sentence domain of validity. **Sample convention**: F1's 85 nuclides are the official score convention (one near-stability-line nucleus per Z chain with experimental mass + the original 10 test nuclei + drip-line neighbours; the closed-form coefficients borrow nothing from AME, so there is no leakage); the 10-nucleus version is only an early test subset used during closed-form selection, and its correspondence is marked in the §9 selection table.

### 8.1 F1: neutron separation energy $S_n = B(A,Z) - B(A-1,Z)$

**① Definition**: $S_n^{\mathrm{model}} = B_{\mathrm{model}}(A,Z) - B_{\mathrm{model}}(A-1,Z)$ (Route-2 closed form).

**② Derivation**: term-by-term differencing — the $a_v$ term differs by 13.85; the $a_s$ term by $-1.59[A^{2/3} - (A-1)^{2/3}]$; the $a_c$ and $a_{\mathrm{sym}}$ terms are A- and Z-invariant under neutron removal → only the A and surface terms participate; odd-even block differences cancel entirely (same Z, same parity).

**③ Criterion** (pre-registered): $|dev| \le 2.0$ pass, $\le 5.0$ half, >5.0 fail.

**④ Score table** (85 nuclides):

| Mass region | n | pass | half | fail | MAE (MeV) |
|---|---|---|---|---|---|
| A<40 | 15 | 7 | 7 | 1 | 2.34 |
| 40≤A<140 | 53 | 41 | 12 | **0** | 1.45 |
| A≥140 | 17 | 11 | 6 | **0** | 1.30 |
| All | 85 | 59 | 25 | 1 | **1.58** |

Residual vs A and the |dev| histogram are in **Figure 4** (shaded bands: ±2 MeV pass, ±5 MeV half).

Route-1 closed form, same convention (added after the consistency revision): **MAE 1.54, 55 pass/30 half/0 fail** — the two closed forms are equivalent on the official 85-nucleus convention (0.04 MeV apart), and Route 2 has fewer extrapolation failures (zero fail); the 10-nucleus early subset comparison is in the §9 table note.

**⑤ Failure attribution**: the single "fail" is $^{18}$F (A<40 extrapolation region); the A<40 bias source is the absence of cluster/shell structure in the liquid-drop form in light nuclei (mutually corroborating the §10 valve failure).

**⑥ Domain**: **the $S_n$ closed form applies for A≥40**, MAE 1.58 MeV; A<40 is order-of-magnitude reference only.

### 8.2 F2: β⁻ decay energy $Q_{\beta^-} = B(A,Z+1) - B(A,Z) - 0.782$

**① Definition**: $Q_{\beta^-}^{\mathrm{model}} = B_{\mathrm{model}}(A,Z+1) - B_{\mathrm{model}}(A,Z) - 0.782$ (A invariant, Z+1).

**② Derivation**: the A term cancels; the $a_c$ term differs by $-a_c[(Z+1)^2 - Z^2]/A^{1/3}$; the $a_{\mathrm{sym}}$ term by $-a_{\mathrm{sym}}[(N-Z-2)^2 - (N-Z)^2]/A = +4a_{\mathrm{sym}}(N-Z)/A\cdot(-1)$…; term-by-term expansion in Appendix C. The Lane convention ($a_{\mathrm{sym}}=23$) is listed separately.

**③ Criterion** (pre-registered; same three tiers as F1).

**④ Score table** (**test set n=10**; four conventions compared. A separate **n=32 extended regression set** is used only for the shell-correction regression; the relation between the two is stated in ⑤):

| Convention | MAE (MeV) | pass/half/fail |
|---|---|---|
| Baseline ($a_{\mathrm{sym}}=19.5$) | 4.95 | 1/1/8 |
| **Lane ($a_{\mathrm{sym}}=23$)** | **3.80** | 1/1/8 (deviation systematically narrowed, 7/10 improved) |
| Baseline + pairing patch | 5.45 | rejected (see below) |
| Lane + pairing patch | 4.33 | rejected |

**⑤ Failure attribution (shell-correction regression)**. Regression of the F2 residual dev on the parent/daughter shell corrections ($B_{\mathrm{exp}} - B_{\mathrm{smooth}}$) (**regression set n=32**: the 10 test nuclei + extension to 32 nuclei with experimental masses; the test-set/regression-set relation is stated explicitly here):

> $dev = +1.028\,s_m - 1.041\,s_d - 2.512$,  **$R^2 = 0.978$, residual RMS = 0.41 MeV**

Coefficient-by-coefficient comparison with analytic expectation: **$a$ is 2.8% from +1, $b$ is 4.1% from −1** (both agree, showing the difference structure is indeed dominated by the shell-correction difference); **the intercept −2.512 is 1.73 MeV from the analytic −0.782** — an additive constant that shifts all dev uniformly without changing MAE, ranking, or per-nucleus verdicts; its source (Lane convention vs smooth-mass-formula constant convention) is not fully resolved, and this paper marks it honestly and **does not count it as model evidence**.

**A precise statement**: the closed form has no shell term → dev equals, *by construction*, the shell-correction difference plus a constant; hence the high $R^2$ proves "the residual is fully decomposable into the shell-correction difference", *not* the causal proposition "shell effects cause the residual" — the latter is guaranteed by construction, the former is what this regression empirically establishes. The defensible wording: **the residual is fully absorbed by the shell-correction difference ($R^2=0.978$), and beyond the shell correction no structure remains** (remaining RMS 0.41 MeV; correlation with magic-number distance $r=+0.08$, i.e. no residual shell structure).

Pairing term falsified on two routes (honest record): the closed-form layer's independent fit yields negative δ (−3.03, violating physical intuition) and is residual-contaminated (fit RMS 14.35 vs signal 2–3); the engine-level same-A differencing probe (198 points) has noise/signal ≈50×, unextractable; the quantitative evidence that the engine lacks the odd-even mechanism = the even-even minus odd-odd difference is systematically small by ~4 MeV (7/8 negative). **Conclusion: the pairing term cannot enter this framework with zero borrowing**, and is listed as engine-level future work. (The two routes of pairing falsification — closed-form-layer parity classification and engine-level same-A differencing — are in Figure 5.)

**⑥ Domain**: $Q_{\beta^-}$ usable for order-of-magnitude and trend judgment (MAE 3.80); not a precise decay-energy prediction.

### 8.3 F3: α decay energy $Q_\alpha = B(A-4,Z-2) + B_\alpha - B(A,Z)$

**① Definition**: $Q_\alpha^{\mathrm{model}} = B_{\mathrm{model}}(A-4,Z-2) + B_\alpha^{\mathrm{exp}} - B_{\mathrm{model}}(A,Z)$. **$B_\alpha = 28.296$ MeV is the experimental $^4$He binding energy, the only explicit experimental borrowing in the whole paper**, placed next to the formula rather than in a footnote: the difference formula needs the absolute binding energy of the α particle itself, and the model closed form's A=4 extrapolation precision is insufficient (light-nucleus failure, §8.1).

**② Derivation**: α differencing removes 2p2n simultaneously — the $a_c$ difference $\propto Z^2[1/A^{1/3} - 1/(A-4)^{1/3}]$ and the $a_{\mathrm{sym}}$ term largely cancel at constant (N−Z); the dominant term is the surface difference $a_s[A^{2/3} - (A-4)^{2/3} - 4^{2/3} + \dots]$.

**③ Criterion** (pre-registered; same as F1).

**④ Score table** (**14 nuclides** = 12 A≥180 Th–Cf α-chain nuclei + $^{147}$Sm/$^{152}$Gd two A<180 light-end controls, itemized in `q_formula_v4.json`): **MAE 1.32 MeV, 3 pass + 9 half + 2 fail**; restricted to the 12 A≥180 nuclei, **MAE 1.06 MeV**, the largest deviation $^{222}$Rn (1.78 MeV), the rest ≤1.31 MeV. **Sample-convention note**: the two light-end nuclei are *pre-registered as out-of-domain controls*, so their failure is the expected result, not a post-hoc removal — not removed, honestly counted into the 14-nucleus MAE.

**⑤ Failure attribution**: the two "fails" are the light-end Sm-147/Gd-152 — small-Q differences are surface-sensitive (sensitivity ~$\Delta a_s A^{2/3}$), and the sign can flip. The rescaled closed form ($a_s=19.01$) degrades F3 to 1.89 with a systematic +2.5–2.9 heavy-nucleus bias — **the surface term over-compensates $Q_\alpha$**, so F3 keeps the Route-2 closed form.

**⑥ Domain**: **the $Q_\alpha$ closed form applies for A≥180**.

### 8.4 F4: fission energy $Q_{fs} = B(f_1) + B(f_2) - B(\text{parent})$

**① Definition**: $Q_{fs} = B_{\mathrm{model}}(A_1,Z_1) + B_{\mathrm{model}}(A_2,Z_2) - B_{\mathrm{model}}(A,Z)$, fragments conserving complementarity ($A_1 + A_2 + A_n = A$, $Z_1 + Z_2 = Z$).

**② Derivation**: fragmentation turns one nucleus into two — the volume term is unchanged, the **surface term is paid twice** ($\Delta a_s \approx a_s[2A_f^{2/3} - A_p^{2/3}] \approx 110$–150 MeV), and the Coulomb term is released (inter-fragment Coulomb difference). **A closed form without a surface term must be systematically high on F4** — the failure mechanism of the fourth iteration.

**③ Criterion** (pre-registered): $|dev| \le 20\%\cdot Q \approx 35$ MeV band.

**④ Score table** (five fission channels; v4 surface-less vs v5 rescaled revival):

| Channel | dev v4 (MeV) | dev v5 (MeV) | Verdict |
|---|---|---|---|
| U-238 → Xe-140 + Sr-95 + 3n | −6.7 → (surface-less ~+177 high) | −6.7 | pass |
| U-238 → Pd-119 + Pd-119 | — | +2.2 | pass |
| U-235 → Ba-141 + Kr-89 + 3n | — | −9.0 | pass |
| Pu-240 → La-145 + Rb-92 + 3n | — | +0.3 | pass |
| Cf-252 → Ba-142 + Mo-107 + 3n | — | −1.0 | pass |
| **Summary** | MAE 184, 0/5 | **MAE 3.83, 5/5** | |

Five-channel v4 (surface-less) vs v5 (rescaled) deviation comparison (log axis + 20%·Q≈35 MeV band) in **Figure 6**.

**⑤ Falsify–correct narrative (the most complete one in this paper)**: round-4 F4 all-fail → attribution "model lacks the surface term; fission fragmentation pays no surface cost" (quantitative deficit 110–150 MeV, same order as measured) → round-5 adds the surface term ($a_s=19.01$) per the pre-registered well-depth-rescaling clause → 5/5 revival. This closed loop is a direct validation of the "attribution → pre-registered correction → re-test" discipline.

**⑥ Domain**: fission Q value usable at the formula layer (5/5 in the 20%·Q≈35 MeV band); no per-fragment refinement.

---

## §9 Closed-Form Selection and Unification

| Observable | Route 2 (13.851/1.592/0.7265/19.50) | Route 1 (16.060/19.010/0.7560/19.50) | Selection |
|---|---|---|---|
| F1 $S_n$ | **1.58 (85-nuclide official convention, 1 fail)** | 1.54 (85-nuclide same-convention re-test, 0 fail) | equivalent (either usable; F2/F3 same form, so Route 2 reported) |
| F2 $Q_{\beta}$ | **3.80 (+Lane, 10 nuclides)** | 5.31 | Route 2 |
| F3 $Q_\alpha$ | **1.32 (14 nuclides, incl. 2 light controls)** | 1.89 (heavy-nucleus bias) | Route 2 |
| F4 $Q_{fs}$ | 184 | **3.83 (5/5)** | Route 1 |

Table note: the two closed forms' 85-nucleus F1 scores are independent re-tests from the same batch (completed in this round's consistency revision; archive `results/f1_expanded_85.json` + the Route-1 re-test record); the early 10-nucleus subset comparison (Route 2 0.83 / Route 1 1.02) is kept as a selection-process record, not an official score. The MAE evolution over five iterations of the four channels (log axis) is in **Figure 7**.

**Physical basis (why not one-size-fits-all)**: the two observable families have different sensitivity structures to the closed-form terms — **difference-response type** (F1/F2/F3: differencing cancels the volume term, insensitive to the absolute scale of $a_v$, sensitive to missing shell/cluster) uses Route 2 (whose $a_v$ calibration is self-consistent with the Fermi kinetic energy); **surface-sensitive type** (F4: fragmentation pays a net extra surface) must use the explicit-surface Route 1. This is not parameter picking but a consequence of the difference structure; the two closed forms' coefficients are independent in origin (§4 dual routes), so there is no shared degree of freedom.

---

## §10 Valve Test: The Measured Boundary of Energy-Side Independence (Decisive Negative Result)

**Pre-registered criterion text**: replace the $S_n$ energy side of the decay-verdict layer (Paper C) from AME with this paper's closed form, run the full A≤40 blind protocol — ≥95% doubles the paper's value; 90–95% publishable (gate needs AME fallback); **<90% abandon the "independent prediction" wording**.

**Four-convention matrix** (30 seeds × 70/30, verdict pipeline reproduced verbatim, gate scale ERR ∈ {AMEerr, 0.3 MeV fixed} × calibration ∈ {raw, training-region linear calibration}):

| | gate=AMEerr | gate=0.3 MeV |
|---|---|---|
| raw | 33.93% | 57.19% |
| calibrated | 66.74% | 66.74% |

Four protocols and four conventions side by side in **Figure 8** (error bars = 1 std over 30 seeds; dashed line = pre-registered 90% red line; dotted = P4 permutation control ~30%).

**Verification that the two "calibrated" cells are bit-identical (not a copy error)**: per-seed verdict comparison across the two conventions (30 seeds all re-run, archive `results/t4_cal10_identical.json`) — **30/30 seeds have identical blind accuracy**. Reason: AMEerr's median is keV-scale, so thr = sig_k·ERR is far below the MeV-scale closed-form Q differences under both conventions, making the open-channel verdicts nucleus-by-nucleus identical; and the fitted random-number sequence is the same → per-seed accuracies are bit-equal. The two cells' coincidence is therefore a structural result of the verdicts, now annotated.

**All four cells < 90% pre-registered line**. The absolute value is gate-scale-sensitive (raw: 23 pp), but the **negative result is doubly robust to gate and calibration** — evidence of protocol strength, not protocol failure (the 23 pp must be reported; concealment is the fatal move).

**Failure mechanism**: in light nuclei the closed-form $S_n$ is an ~14 MeV near-constant (no shell/cluster structure), while drip-line verdicts need 0–5 MeV resolution; linear calibration does not generalize (training RMS 6.4 → test RMS 19.7 MeV).

**Conclusion**: the closed-form layer has no energy-independent closed loop to per-nucleus verdicts; A≤40 energy verdicts require explicit clustering; the verdict layer keeps the "mass-table pipeline + model structural gate" positioning (net +1.5 pp, Paper C). **This valve's negative result is a measured boundary under the pre-registered criterion, honestly replacing any "independent prediction" wording.** (The valve's three pillars vs the T5 pre-regression residual are in **Figure 9**.)

---

## §11 Discussion and Limitations (Eight Items; Each: Phenomenon → Quantitative Evidence → Scope)

1. **K=547 makes no GMR prediction.** Phenomenon: incompressibility exceeds the empirical 240±30; evidence: analytic lower bound $K_{\mathrm{att}} = 4k_n s c \rho_0^{-1/3}$ (166.5 MeV at rescaled depth; full-model numerical 547), four parameterizations and two lattices hitting the same wall (§3.2/§7.4); scope: compression-class observables only; this paper's target observables do not depend on K.

2. **A≤40 energy-side closed form fails.** Phenomenon: closed-form $S_n$ is an ~14 MeV near-constant in light nuclei; evidence: valve-test four conventions 33.9/57.2/66.7/66.7% all <90% pre-registered line (§10); scope: light-nucleus energy verdicts must take the explicit-cluster route; "independent prediction" wording abandoned.

3. **F2 systematic negative bias.** Phenomenon: −3 to −9 MeV; evidence: on the 32-nucleus regression set the residual is fully absorbed by the shell-correction difference ($R^2=0.978$, residual RMS 0.41 MeV), no structure beyond shell correction (§8.2); scope: $Q_\beta$ order-of-magnitude and trend usable, no precise prediction.

4. **F3 only for A≥180.** Phenomenon: light-end small-Q differences flip sign, both $^{147}$Sm/$^{152}$Gd out-of-domain controls "fail"; evidence: difference sensitivity to $a_s$, ~$\Delta a_s A^{2/3}$; scope: Th–Cf heavy region, 12 nuclei MAE 1.06 MeV, ≤1.31 MeV except $^{222}$Rn (1.78).

5. **Drip-line precision insufficient.** Phenomenon: drip-line verdicts need 0.1 MeV; evidence: closed-form MAE at best 0.83 (F1 early subset) / 1.58 (official); scope: drip-line verdicts rely on NUBASE patches (12 nuclei, list public), whose correctness does not come from the model.

6. **Annealing stochasticity.** Phenomenon: geometric shape quantities vary with seed; evidence: 6 nuclei × 5 seeds, $Q_0$ std up to 1154 ($^{20}$Ne), α count A≥20 scatter 1–5, B/A scatter 3.4% (§2.4); scope: energy conclusions robust, A≥20 geometric conclusions downgraded to qualitative.

7. **Charge radius systematically small.** Phenomenon: median ratio 0.794 (n=97); evidence: full per-nucleus table (Appendix A CSV); scope: absolute radius unusable, ratio trend usable.

8. **Half-life order-of-magnitude precision.** Phenomenon: single-value log ft anchoring; evidence: multivariate regression $R^2=0.138$ (Paper C); scope: half-life order-of-magnitude only; WKB α half-life systematically ~5 orders high (mode verdict unaffected).

**Parameter-tier list (pre-registered convention)**: v5.0 rule constants ($s$, window [0.8,2.2], $k_{\max}=4$, Pauli core) — frozen, non-adjustable; geometric engine $C_{\mathrm{REP}}=12$, $R_{\mathrm{EQ}}=1.608$ — four-light-nucleus calibration; nuclear-matter calibration $k=56.97$, $C=699.1$, $\delta=0.2$ — analytic-solution calibration; Route-1 closed-form four coefficients — cluster sum × rescale, zero per-nucleus fitting; Route-2 closed-form four coefficients — model-grid least squares, zero experimental borrowing; verdict-layer six parameters — training-set-only fitting within the blind protocol; $B_\alpha$ anchor 28.296 MeV — the only explicit experimental borrowing (§8.3).

---

## §12 Conclusion

Two local rules (a windowed linear bond well + four-bond mutual-consent saturation) generate, in the infinite limit, the degenerate family of ideal tetrahedral four-coordinated networks, with tetrahedral coordination homologous to α clustering; the liquid-drop coefficients emerge with zero per-nucleus fitting and the two sampling routes are self-consistent under one energetics (§4.3), and linear differences of the same closed form yield four release-energy formulas ($S_n/Q_\beta/Q_\alpha/Q_{fs}$) each with a measured domain of validity; geometric observables are quantified over the full H→Ca set, among which the α-cluster count and the quadrupole sign are discriminators the liquid-drop/shell models lack — the geometric layer's current domain is A≤40, which is this version's boundary rather than a division-of-labor conclusion, and extension to A>40 is explicit follow-up work; the energy-layer domains (A≥40 closed form, no GMR prediction, drip-line patch dependence, valve-test all-fail) are delimited by pre-registered negative results. Every result in the main text is a test/reproduction against known data; the only forecast in this paper — the ³⁴Ca blind prediction (Appendix F, pre-registered 2026-11-24) — is a falsifiable-in-whole output of the model for an unmeasured quantity. This paper's falsification conditions (§1.5), failure list (Appendix A), and full criterion text are published alongside; the model's range of applicability is defined by measured boundaries, not extrapolated from success cases.

---

## Figure Captions (Figures 1–9)

> All figures are generated by `src/make_figures.py`, output to `repo/figures/`. Captions numbered by order of first reference; repository filenames in parentheses.

**Figure 1** (`A_fig1_lattice_scan.png`) | **Saturation points and binding energies of four parameterizations × two lattices (pathology control).**
(a) Saturation density $\rho_0$; (b) binding energy per nucleon $E/A$. Bars are valid solutions per scheme; black dashed lines are empirical ($\rho_0=0.16$ fm⁻³, $E/A=-15.8$ MeV). The three FCC groups deviate by +88% ($\rho_0$), −28% ($E/A$), or have no valid solution at all (FCC direct analytic, all K>350 over the 280 grid); the diamond four-bond solution lands on the empirical values, missing only K. Criterion lines: $\rho_0$ off >25%, $E/A$ off >15%, or K>350 marks a pathology.

**Figure 2** (`A_fig2_mass_table.png`) | **Binding-energy residual vs mass number over 2811 nuclides (zero per-nucleus fitting).**
Scatter of $B_{\mathrm{model}} - B_{\mathrm{exp}}$, coloured by three mass regions: A=40–65 (MAE 6.5), 65–100 (MAE 9.3), 100–240 (MAE 14.3); overall MAE 12.6 MeV, median relative deviation −0.3%, all nuclides within 15% relative. Coefficients $a_v=16.060$, $a_s=19.010$, $a_c=0.756$, $a_{\mathrm{sym}}=19.50$ (Route 1, zero per-nucleus fitting).

**Figure 3** (`A_fig5_geometry_regions.png`) | **Pass rates of the four geometric observables by mass region (H→Ca, 389 nuclides).**
Four panels: binding energy, quadrupole sign, α-cluster count, charge radius. Bars are the "pass"-tier fraction in each mass region; black dashed line = overall pass rate. n values: B 389, $Q_s$ 65, α 10, $r_{ch}$ 97. **"Pass" tier only, excluding "half-pass"** — this figure differs from the "including half-pass" convention of the §5.4 table; do not conflate when reading.

**Figure 4** (`A_fig3_f1_expanded.png`) | **F1 neutron-separation-energy residual (85 nuclides).**
(a) $dev = S_n^{\mathrm{model}} - S_n^{\mathrm{exp}}$ vs A, coloured by three-tier verdict (pass n=59 / half n=25 / fail n=1, the single "fail" $^{18}$F); shaded bands ±2 MeV (pass) and ±5 MeV (half). (b) |dev| histogram, black dashed line = MAE 1.58 MeV. Sample: A<40 (15), 40≤A<140 (53), A≥140 (17).

**Figure 5** (`B_fig3_pairing_falsify.png`) | **Two-route falsification of the pairing term.**
(a) Closed-form layer: mean residual of even-even / odd-A / odd-odd classes, error bars ±1 s.e. (from class sizes); the three classes are indistinguishable within an RMS≈14 MeV residual; the independent fit's δ=−3.03 MeV (negative, violating physical intuition) is rejected. (b) Engine layer: same-A even-even minus odd-odd binding-energy difference, 7 of 8 groups negative, mean −3.9 MeV (red dashed line) — direct evidence that the engine lacks the odd-even mechanism.

**Figure 6** (`B_fig2_f4_revival.png`) | **The surface term revives the F4 fission energy.**
|Q_fs deviation| for five fission channels (log axis). Red bars: v4 (surface-less closed form) MAE 184 MeV, 0/5 pass; blue bars: v5 (rescaled, $a_s=19.01$) MAE 3.83 MeV, 5/5 pass. Black dashed line = pre-registered 20%·Q≈35 MeV band. Band text in §8.4③.

**Figure 7** (`B_fig1_evolution.png`) | **MAE evolution over five iterations of the difference-formula family.**
Log-MAE (MeV) vs v1→v5. Four curves for F1 ($S_n$), F2 ($Q_{\beta^-}$), F3 ($Q_\alpha$), F4 ($Q_{fs}$). Key turning points: F4 drops from 184 to 3.83 at v4→v5 (surface term added), F3 degrades from 1.32 to 1.89 (surface term over-compensates $Q_\alpha$, so F3 keeps the v4 closed form).

**Figure 8** (`A_fig4_blind_and_T4.png`) | **Four-protocol blind test and the four-convention valve matrix side by side.**
First four bars P1–P4: P1 random 70/30 (30 seeds) 97.4±1.1%, P2 element-block 97.1±1.3%, P3 extrapolation 96.75%, P4 training-label permutation 30.1%. **Error bars = 1 std over 30 seeds** (the matrix's four cells are fixed-seed, no error bars). Last four bars: T4 four conventions (raw/calibrated × gate=AMEerr/0.3 MeV). Black dashed = pre-registered 90% red line; grey dotted = P4 permutation control.

**Figure 9** (`B_fig4_valve_and_T5.png`) | **Valve test vs shell-correction attribution.**
(a) Three bars: raw closed form 33.9%, training-region linear calibration 66.7%, AME reference 97.4%; black dashed = 90% red line — all three below, the measured boundary of "energy side not independent". (b) T5 **pre-regression** residual (Q_β dev) vs distance-to-nearest-magic-number (n=10, Pearson r=−0.64, Spearman in archive): the negative correlation shows the residual is shell-dominated; **after extension to n=32 and subtraction of the shell-correction difference, the remaining correlation with magic-number distance drops to r=+0.08 (§8.2)** — the two panels' n and statistics differ and are not interchangeable.

---

## Appendix 0: Unified Symbol Table

| Symbol | Meaning | First appearance |
|---|---|---|
| $s$ | linear well depth (MeV/fm; nominal 5, rescaled 13.05) | §2.1 |
| $r_1, r_2$ | well window edges 0.8 / 2.2 fm | §2.1 |
| $k_{\max}$ | max bonds per nucleon = 4 | §2.1 |
| $k_n$ | bonds per nucleon in the infinite limit = $k_{\max}/2$ = 2 | §7.2 |
| $c$ | diamond geometry constant (dimensionless) = $d_0\rho_0^{1/3}$ = 0.866 | §7.2 |
| $C, r_{eq}$ | Pauli-core strength and equivalent parameter (engine $C_{\mathrm{REP}}=12$, $R_{\mathrm{EQ}}=1.608$) | §2.1/§2.4 |
| $a_v, a_s, a_c, a_{\mathrm{sym}}$ | four liquid-drop coefficients (two routes, §4) | §4 |
| $B_\alpha$ | $^4$He experimental binding energy 28.296 MeV (only explicit borrowing) | §8.3 |
| $s_m, s_d$ | parent/daughter shell corrections ($B_{\mathrm{exp}} - B_{\mathrm{smooth}}$) | §8.2 |
| $K_{\mathrm{att}}$ | attractive-term incompressibility = $4k_n s c \rho_0^{-1/3}$ | §7.4 |
| $Q_s, Q_0$ | electric quadrupole moment (experimental conversion / model eigenvalue) | §5.3 |

## Appendix A: H→Ca Reaction Chain and Nucleus-by-Nucleus Geometric Total Table

Representative in-text traces: $^4$He (4 nucleons, k distribution {3:4} — perfect tetrahedron), $^{16}$O (16 nucleons, k distribution {0:1, 1:1, 3:5, 4:9}), $^{48}$Ca (48 nucleons, k distribution {0:5, 1:2, 2:4, 3:16, 4:21}). **Full table** (389 nuclides × verdict type × slowest channel + 315 nuclides × $E_{\mathrm{model}}/E_{\mathrm{bond}}$/k distribution) in the supplement `appendixA_summary.xlsx` (315 rows × 12 columns) and `results/big_validation_table.csv` (389 rows × 22 columns). The failure-nucleus list is published with the table.

## Appendix B: Single-Nucleus Per-Nucleon Parameter Table

One row per nucleon: index / type (p/n) / coordinates (x,y,z) / main-bond count $k_{\mathrm{main}}$ / soft-bond count $k_{\mathrm{soft}}$ / neighbour indices / $r_{ij}$ / $V_{ij}$ / bond type (main/soft) / Pauli-push flag / hard-core flag / Coulomb pair energy. Three in-text tables (geometric-engine live relaxation, seed=13, table energies self-consistent with per-bond energies):

- **$^4$He**: 4 nucleons, 6 bonds (each nucleon k=3, tetrahedron), E = −23.09 MeV;
- **$^{16}$O**: 16 nucleons, 60 bonds (main+soft), E = −91.87 MeV (B/A 5.74);
- **$^{48}$Ca**: 48 nucleons, E = −250.73 MeV (B/A 5.22).

Full per-nucleon long table (203,662 bond rows, 315 nuclei, verdict-layer engine convention) in the supplement `appendixB_all_nucleons.csv`. The two-engine-convention difference is declared in §2.4 and the repo README. The Rule-2 parameter table (per-channel Q values, phase-space f, WKB penetration, branching ratios) is in Paper C's appendix.

## Appendix C: Formula List and Difference Expansions

C1 windowed linear well; C2 Pauli repulsion and pairing-equilibrium constraint; C3 four-bond saturation and mutual consent; C4 second-neighbour geometry (2.605 fm); C5 bond-angle threshold (87.2°); C6 analytic upper bound −6.05 MeV/A; C7 liquid-drop closed form (two routes); C8 Q_β trigger criterion; C9 phase-space integral and Fermi function; C10 WKB penetration; C11 $K_{\mathrm{att}} = 4k_n s c \rho_0^{-1/3}$; C12 F1 difference expansion ($a_v$ difference 13.851, surface difference $-1.592[A^{2/3}-(A-1)^{2/3}]$, $a_c/a_{\mathrm{sym}}$ cancel); C13 F2 difference expansion ($a_c$: $-a_c(2Z+1)/A^{1/3}$, $a_{\mathrm{sym}}$: $-4a_{\mathrm{sym}}(N-Z)/A$ approximation); C14 F3 difference expansion; C15 F4 surface-deficit formula $\Delta a_s[2A_f^{2/3}-A_p^{2/3}]$; C16 shell-correction regression (§8.2). Intermediate steps of all expansions in the `repo/src/q_formula_*.py` comments and the README mapping table.

## Appendix D: Code and Algorithms

| Module | Function | Path |
|---|---|---|
| nm_formula_v4 | nuclear-matter calibration (diamond-network direct solution) | `src/nm_formula_v4.py` |
| geometry_validation | geometric engine (L-BFGS + basin hopping) | `src/geometry_validation.py` |
| big_validation | H→Ca full-set four-observable verdicts | `src/big_validation.py` |
| blind_test / chain_predict | blind protocol and reaction-chain verdicts | `src/blind_test.py` etc. |
| q_formula_v4–v6 | difference-formula family, four iterations | `src/q_formula_*.py` |
| valve_test / t4_* / t5_* | valve and sensitivity experiments | `src/valve_test.py` etc. |
| appendixB_export | per-nucleon parameter export | `src/appendixB_export.py` |
| make_figures | all nine figures | `src/make_figures.py` |

Third-party dependencies: numpy/scipy (BSD), pandas (BSD), matplotlib (PSF); basin hopping is `scipy.optimize.basinhopping`, the Wales–Doye algorithm [6]. **Code license: MIT** (repo/LICENSE). **Paper-number → script → output mapping** in the repo README and `REPRODUCE.md` (every number traceable, including entries 19–22 for the Appendix-F prediction card).

## Appendix E: Reproduction

```bash
git clone <repo-url>      # repository URL and Zenodo DOI in the Data Availability statement (filled at submission)
pip install -r requirements.txt   # numpy / scipy / pandas / matplotlib / openpyxl
bash run_all.sh core    # core results ~10 min, output field-by-field identical to the paper tables
bash run_all.sh full    # full geometric/blind test, results archived in results/
bash run_all.sh figs    # redraw the 9 figures only (reads results/ archives, no recomputation)
```

Fixed seeds (seed=13 verdict/geometric main convention; 30 seeds for the blind test) guarantee exact reproduction. **The code is released with the paper under the MIT license**; the repository URL and Zenodo DOI are filled at submission (Data Availability statement); the traceability path of every number is in the repo README and REPRODUCE.md mapping tables.

## Appendix F: Pre-Registered Blind Prediction (³⁴Ca; pre-registration date 2026-11-24)

This appendix is the **only forecast content** in this paper, kept strictly separate from the main text's "test/reproduction": every result in the main text is a test against known data; this appendix is the model's output for a **never-measured quantity**, whose value lies in being confirmable or falsifiable as a whole by future experiments. The pre-registration precedes measurement: ³⁴Ca (Z=20, N=14) had no experimentally observed property as of 2026-11-24 (the RAON 2024 experimental proposal lists it as a target nucleus, and the quantity they propose to measure is exactly item P3 of this card).

**Relation of the prediction model to the main-text closed forms (must-read)**: this card uses a six-parameter macroscopic-microscopic closed form (five liquid-drop parameters + a spherical shell correction C_shell, a macroscopic-microscopic implementation of this paper's windowed linear well, with the smoothing constant fixed at c=0.25); its 6 parameters are least-squares fitted on **2484 measured nuclides** (A≥16, no AME2020 estimates) to extrapolate the model to unmeasured nuclides — **a different source** from the two zero-fitting emergent closed forms of §4 (which involve no per-nucleus fitting), though sharing the same potential form and rule system. The prediction error bar is the residual RMS over the prediction neighbourhood (Z=16–24, N=8–20; 42 measured nuclides): **1.55 MeV**.

| # | Prediction | Value | Notes |
|---|---|---|---|
| P1 | Fusion energy Q(⁹C + ²⁵Si → ³⁴Ca + γ) | **+18.03 ± 1.55 MeV** | both participants are measured, β-bound nuclei (NUBASE2020 ground states, no particle emission) |
| P1′ | Alternate channel Q(¹⁰C + ²⁴Si → ³⁴Ca + γ) | +11.74 ± 1.55 MeV | same |
| P3 | Two-proton decay energy E_T(³⁴Ca → ³²Ar + 2p) | **+2.33 ± 1.55 MeV** | ³²Ar measured; a **single-end prediction** (only one mass, B(³⁴Ca), is predicted); the quantity RAON proposes to measure |
| P5 | B(³⁴Ca) | 244.07 ± 1.55 MeV | no entry for this nuclide in AME2020 |

**Geometric layer** (single seed=13 in-situ relaxation; the seed-robustness limitation of §2.4/§11.6 applies; qualitative reference): α clusters 3/3 (the ideal maximum), Q₀ = +170.33 fm² (prolate), r_ch = 4.09 fm.

**Falsification conditions**: the measured E_T falling within 2.33 ± 3.10 MeV (2σ) passes the prediction; falling outside, or ³⁴Ca proving p-decay dominated, or a measured S₂p < 0 (unbound), each constitutes falsification — and upon falsification the result enters the failure list per the §10 practice.

**Rejection record (why ³⁴Ca)**: candidate selection was an exhaustive enumeration (measured + measured → unmeasured; 52 candidate C nuclides; bound-nucleus filter archived in `results/fusion_channels_bound.json`, `nubase_bound.json`). (i) The light neutron drip line is rejected: ²⁸O was measured in 2023; the model deviates by −18 MeV while the AME estimate is off by only 0.83 — the model's symmetry term fails systematically in the extremely neutron-rich regime (mutually corroborating with the §10 valve failure). (ii) The superheavy region is rejected: the α-chain differences are not self-consistent (²⁹³Lv→²⁸⁹Fl: model −4.4 vs AME +10.9 MeV). (iii) The Z=16–24 proton-drip-line neighbourhood is selected: the validation cohort (³⁴K/²³Si/²⁷S etc.) was measured in 2023–2024 and agrees with the AME estimates; the local error bar of 1.55 MeV is verified on 42 measured nuclides. Reproduction script `src/prediction_ca34.py` (one command, digit-exact), archives `results/prediction_ca34_numbers.json`, `ca34_structure_seed13.json`.

---

## References (W4 rule: cite only what is used; fewer rather than padded; verify before submission)

[1] AME2020: M. Wang, W.J. Huang, F.G. Kondev, G. Audi, S. Naimi, Chin. Phys. C 45, 030001 (2021) (all masses/Q values/binding energies in this paper).

[2] NUBASE2020: F.G. Kondev, M. Wang, W.J. Huang, S. Naimi, G. Audi, Chin. Phys. C 45, 030002 (2021) (decay labels and half-lives).

[3] FRDM2012: P. Möller, A.J. Sierk, T. Ichikawa, H. Sagawa, At. Data Nucl. Data Tables 109–110, 1 (2016) (precision-comparison ruler).

[4] Charge radii: I. Angeli, K.P. Marinova, At. Data Nucl. Data Tables 99, 69 (2013) (§5.2 experimental values).

[5] Quadrupole moments: N.J. Stone, At. Data Nucl. Data Tables 111–112, 1 (2016) (§5.3 experimental values and recommended-value sub-sample).

[6] D.J. Wales, J.P.K. Doye, J. Phys. Chem. A 101, 5111 (1997) (basin-hopping algorithm); scipy.optimize.basinhopping.

[7] α-cluster reviews: W. von Oertzen, M. Freer, Y. Kanada-En'yo, Phys. Rep. 432, 43 (2006); M. Freer, Rep. Prog. Phys. 70, 2149 (2007) (§1.1 F4, §3.3 controls).

[8] Classical liquid-drop model: C.F. von Weizsäcker, Z. Phys. 96, 431 (1935); H.A. Bethe, R.F. Bacher, Rev. Mod. Phys. 8, 82 (1936) (§1.6 historical statement of the liquid-drop limit).

[9] R. Hofstadter, Rev. Mod. Phys. 28, 214 (1956) (electron scattering and nuclear charge distribution; §1.1 F1 saturation-density empirical value).

> **Data Availability**: all raw data from public databases (AME2020/NUBASE2020/IAEA charge-radius and quadrupole-moment compilations); all per-nucleus tables, blind-test records, and code generated by this work are released with the paper under the **MIT license** via the GitHub repository (Zenodo DOI to be filled at submission; the two cross-reference each other). The paper-number → script → output mapping is in the repo README and REPRODUCE.md. The Appendix-A summary table is provided as an Excel workbook (`appendixA_summary.xlsx`, 315 rows × 12 columns); the full per-nucleon long table is `appendixB_all_nucleons.csv` (203,662 bond rows / 315 nuclei, UTF-8).

> Citation-principle statement: this is a phenomenological work of independent conception; indirectly borrowed existing-theory content (the liquid-drop form as the A→∞ limit ruler, the Fermi-gas kinetic-energy expression, the WKB penetration approximation, the classical-rotor conversion) is explicitly marked at the point of use; no shell-model orbits or density-functional theory inputs are used. References are kept few rather than padded.

---

> **Author archive (not for submission)**: this manuscript is the polished final version of v2.1. Principal revisions: F3 sample count corrected from 12 to 14 (the archive `q_formula_v4.json` actually holds 14 entries); the §8.2 regression-intercept vs analytic-value 1.73 MeV discrepancy disclosed honestly (the original "agreement to 3–4%" wording does not hold for the intercept); nine figure numbers and captions completed; references filled in (originally [7]–[9] were empty placeholders); AME2020 author attribution corrected (originally misattributed to Kondev; should be Wang et al.); §7.2 well-depth inversion intermediate number corrected ($0.16^{1/3}=0.543$; the original miswrote 1.842, which is $\rho_0^{-1/3}$). All numbers re-verified against the archives by actual runs.
>
> **v2.2 revision record (2026-11-24)**: (1) the $a_{\mathrm{sym}}=19.50$ two-route agreement regraded from "independent cross-validation" to "self-consistency of one energetics under two sampling schemes" (abstract/§4.3/§7.7/§12); (2) §4.4 adds a same-set FRDM2012 comparison footnote (n=2910, MAE 0.45/RMS 0.63; the 24–28× gap stated honestly; archive `frdm_same_set.json`); (3) §5.2 adds the radius scale-factor diagnostic (±15% in-band 30.5%→52.6%; "right shape, off by one constant scale"); (4) §5.3 α-cluster "fails" reworded to "reversed ordering trend r=−0.603, p=0.065" (computed on this paper's archive, `alpha_cluster_spearman.json`); (5) Appendices D/E and the Data Availability statement now state the MIT license and the REPRODUCE.md mapping; (6) new Appendix F (³⁴Ca blind prediction pre-registration, the only forecast in the paper); (7) wording sweep: no residual "predictive power" claims; "independent prediction" appears only in the valve-test pre-registration context.
