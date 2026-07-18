# Causal Identification of the VPC Under Differential Detection: Theory, a Jointly-Calibrated Estimator, and Sharp Partial-Identification Bounds

*Phase-3 research report and grand roadmap. No real data — synthetic simulation only, as throughout this repository. All numbers below come from code that was executed in this session (`scripts/validation/phase3_causal_demo.py` reproduces the headline table); none are hypothetical.*

---

## 0. Executive summary

This phase reframes the repository's detection-bias problem as a **causal / missing-data identification problem** and delivers three connected advances that, to our knowledge, have not been combined in the I-MAIHDA literature:

1. **Identification theory** (§2): exact conditions under which the true VPC is point-identified, set-identified, or unidentified from under-detected data — resolving *why* the Phase-2 robustness study found the breakdown regimes it did.
2. **A jointly-calibrated estimator** `joint_calibrated_vpc()` (§3): fixes the correction-composition failure documented in `METHODS_NOTE_ROBUSTNESS.md` §6. Against the population estimand in the joint detection+sparsity regime, mean bias improves from **−16.3pp (naive) and −6.7pp (detection-only) to +0.9pp**, with 80–95% CI coverage across regimes.
3. **Sharp partial-identification bounds** `vpc_partial_bounds()` (§4): globally exact min/max of the VPC over the identified set, valid **even in the outcome-dependent-detection regime where every point correction provably fails** — truth-in-bounds 45/45 across all tested regimes including that one.

A key conceptual discovery en route (§3.1): the Phase-2 "composition failure" had *three* distinct causes, one of which was an **estimand subtlety** — the finite-sample `y_true` refit used as "truth" in Phase 2 itself carries sparse-data shrinkage and sits ~6pp below the population parameter in the capstone regime. Once the estimand is fixed to the population VPC, the detection-only correction is revealed to still be materially biased in sparse designs (−6.7pp), and the jointly-calibrated estimator closes that gap.

---

## 1. The grand roadmap (siêu kế hoạch)

A three-paper research program, each stage already having its computational core implemented and validated in this repository:

**Paper A — Methods + simulation (materialized in Phases 1–2).** The two QBA tools (`correct_detection_bias`, `sparse_strata_vpc`), their self-validating recovery tests, and the systematic assumption-violation study with formal non-inferiority / breakdown criteria (`METHODS_NOTE_ROBUSTNESS.md`). Status: **complete**; needs only manuscript formatting.

**Paper B — Causal identification and partial identification for the VPC (this phase).** The identification propositions of §2, the sharp bounds of §4, and the estimand clarification of §3.1. Novel connections: missing-data identification theory and Manski-style partial identification have not previously been applied to variance-partition functionals in intersectional MAIHDA (our PubMed scans of the MAIHDA methods literature found no measurement-error treatment at all). Status: **theory and algorithms complete and validated**; literature-verification pass and manuscript writing remain.

**Paper C — The jointly-calibrated estimator as applied guidance.** §3's estimator, its correlation-aware simulation design, coverage properties, and practical guidance (when to report point correction vs. joint calibration vs. bounds — the "ladder of assumptions" of §5). Status: **estimator complete and validated in three regimes**; broader scenario grid and R mirror remain.

**Beyond (not started, honestly deferred):** K-generalization (other stratum designs), PCV-functional bounds (the PCV is a *ratio* of variances — its identified set requires joint bounds over two correlated functionals, a genuinely harder optimization), survey weights, and longitudinal extensions.

---

## 2. Identification theory

Setting: strata `j = 1..K` with true prevalences `p*_j = expit(η_j)`; the estimand is the population variance `σ² = Var(η)` (equivalently `VPC = σ²/(σ²+π²/3)`). Observation: `y = y* · D` where `D` is the detection indicator.

**Proposition 1 (point identification under covariate-only detection).**
If `D ⊥ y* | S` (detection conditionally independent of the true outcome given the stratum) and `P(D=1|S=j) = d_j > 0` is known up to the hypothesized `(δ, score, baseline)`, then

  `E[ȳ_j] = E[y*·D | S=j] = p*_j · d_j`,

so `p*_j = E[ȳ_j]/d_j` and every functional of `{p*_j}` — including the VPC — is point-identified. *Proof:* one line by conditional independence, as displayed. This is exactly the inversion `correct_detection_bias()` performs, and explains the Phase-2 finding that misspecifying *which covariates* enter the score degrades gracefully (the inversion is still of the right functional form), while misspecifying the score's *functional form* (quadratic arm) breaks non-inferiority (the inversion applies a wrong `d_j` pattern that no δ can repair).

**Proposition 2 (non-identification under outcome-dependent detection).**
If detection may depend on the true outcome itself, `P(D=1|S=j, y*=1) = d_j1` unknown, then the observed stratum prevalence is `q_j = p*_j·d_j1`, and for any candidate `p̃_j ∈ [q_j, 1]` the pair `(p̃_j, d̃_j1 = q_j/p̃_j)` is observationally equivalent to the truth. The identified set for `p*_j` is the full interval `[q_j, 1]` absent further assumptions; the VPC is therefore only **set-identified**. *Proof:* direct construction, as displayed; `d̃_j1 ∈ (0,1]` iff `p̃_j ≥ q_j`. This is precisely why Phase 2's M1-2 arm degraded monotonically and crossed into breakdown at `severity_weight ≈ 0.8–0.9`: the SES-score inversion corrects only the covariate-driven thinning factor, leaving a residual outcome-driven factor `expit(base−δs−w)/expit(base−δs) < 1` uncorrected in every stratum — no stratum-level score can represent it.

**Proposition 3 (sharp bounds; global optimality of the algorithms).**
Under monotone under-ascertainment with per-stratum relative detection no worse than `d_j(δ_max)`, the true stratum logits lie in a box `η_j ∈ [a_j, b_j]`. Since `vpc_latent` is strictly increasing in the variance, sharp VPC bounds are `[min, max]` of `Var(x)` over the box:

- **Maximum**: `Var` is convex, so the max over a compact convex polytope is attained at a **vertex** (`x_j ∈ {a_j, b_j}`). For fixed other coordinates, `Var` is a convex quadratic in `x_j`, maximized at an endpoint — so coordinate exchange strictly ascends through vertices and terminates; multi-start recovers the global vertex. *Verified exhaustively*: 300 random boxes at K=3–11, coordinate-exchange matched full 2^K enumeration in **300/300** cases (and 30/30 again in the fast test suite).
- **Minimum**: a convex program whose solution is the clipped-constant `x_j = clip(c, a_j, b_j)`, with the scalar `c` minimizing the convex `g(c) = Σ_j dist(c, [a_j,b_j])²` — solved exactly by ternary search. Verified against multi-start L-BFGS-B on the same 300 boxes (never worse).

**Proposition 4 (local-to-global identification of the calibration inverse).**
The calibration map `m(σ²) = E[T | σ²]` satisfies: (i) in the no-sampling-noise limit (stratum sizes → ∞) `m(σ²) = σ²` **exactly** — verified numerically to within Monte-Carlo error (max gap 0.058 over `σ² ∈ [0,1.5]` at 60 reps/point, consistent with the MC standard error); (ii) across the three tested regimes, the estimated joint-pipeline curve was strictly increasing in 44/45 grid steps, the single violation (−0.012) being within the MC standard error at 1,500 reps/point. Local identifiability everywhere plus (empirical) strict monotonicity gives **global invertibility** of the calibrated inverse; the implementation adds an isotonic projection so the inversion is well-defined under residual MC wiggle. *Status: exact in the limit, empirically established at finite n; a distribution-theoretic proof (e.g., via convex ordering of the corrected-logit spread in σ²) is left as a stated conjecture for Paper B.*

---

## 3. The jointly-calibrated estimator

### 3.1 Why the Phase-2 composition failed: three causes, now separated

1. **Estimand mismatch.** Phase 2 measured all arms against the finite-sample refit of VPC on `y_true` — but in sparse designs that refit itself is shrunk (capstone: 20.15% vs. population 26.11%). Calibration-type estimators target the *population* parameter, so they looked like they "overshot" when they were partly correcting real shrinkage the reference value still contained.
2. **Noise-model mismatch.** Detection-corrected counts `e/d` carry sampling variance inflated by `(1−pd)/(d(1−p))` relative to plain binomial; the sparse calibration curve assumed plain binomial and attributed the excess spread to true variance.
3. **Correlation blindness.** Stratum prevalence and detection are driven by the *same* SES covariates, hence negatively correlated. An independent-Gaussian effects simulator misses this. Decisive diagnostic (12 replicates, capstone regime): real pipeline statistic 0.863; simulator with the true correlated `p_j` **0.859** (match); independent-Gaussian simulator at the same total variance **0.727** (16% too low → inversion overshoots).

### 3.2 Design

`joint_calibrated_vpc(df, delta, ...)` calibrates the **entire pipeline** — detection-correct, then compute the precision-weighted null-variance statistic — in one simulation whose stratum effects are decomposed as

  `effect_j = m̂_j (WLS main-effects pattern estimated from the corrected logits, held fixed) + N(0, σ²_res)`,

so the simulator inherits the data's own prevalence–detection correlation; only the residual variance `σ²_res` sits on the calibration grid. The reported total is `Var(m̂) + σ̂²_res`; the CI pushes bootstrap replicates of the pipeline statistic back through the inverse curve. Saturation at the grid top is flagged (`saturated=True`).

### 3.3 Results (population estimand, n_rep=20 per regime)

| Regime | naive bias (RMSE) | detection-only bias (RMSE) | **joint bias (RMSE)** | joint 95%-CI coverage |
|---|---:|---:|---:|---:|
| Capstone D+E (n=3500, sparse) | −16.26 (17.02) | −6.72 (9.07) | **+0.87 (9.17)** | 16/20 |
| Harsher (n=2000, shift −3, sparse) | −17.58 (18.16) | −6.08 (10.37) | **+2.24 (10.19)** | 17/20 |
| Dense (n=12000) | −13.86 (14.32) | −7.62 (8.00) | **−1.44 (3.31)** | 19/20 |

The joint estimator cuts point bias by ~85–95% relative to naive and ~65–90% relative to detection-only in every regime, with the tightest gains (RMSE 3.31, coverage 95%) where data are dense. Residual RMSE in sparse regimes (~9–10pp) reflects genuine information scarcity — which the CI honestly expresses — not correctable bias.

## 4. Sharp bounds: results

`vpc_partial_bounds(df, delta_max)` with a generous envelope (`delta_max=1.6` = 2× the true generating strength):

| Regime | truth-in-bounds | mean width |
|---|---:|---:|
| Capstone D+E | 15/15 | 62.1pp |
| Dense detection | 15/15 | 76.5pp |
| **M1-2 outcome-dependent detection, severity_weight=1.0 (point-ID provably broken)** | **15/15** | 75.7pp |

Width is the honest price of weak assumptions and contracts sharply as the envelope tightens: at `delta_max = 0.8` (the true value) the dense-regime bound is **[8.6, 30.0]** (width 21.5pp, containing the ~26 population truth); `delta_max = 0` collapses the bounds to a point (identity check in the test suite).

## 5. The ladder of assumptions (practical guidance)

1. **Willing to assume covariate-only detection with a specific (δ, score)?** → point correction (`correct_detection_bias`); sweep δ (`vpc_detection_bounds`) as sensitivity analysis.
2. **Also facing sparse strata, and want the population parameter?** → `joint_calibrated_vpc` (this phase). Do **not** compose the two Phase-1/2 corrections sequentially (§3.1).
3. **Unwilling to rule out outcome-dependent detection?** → point identification is impossible (Prop. 2); report `vpc_partial_bounds` with a defensible `delta_max` envelope, plus the tipping-point (`detection_tipping_point`) in the spirit of E-values (VanderWeele & Ding, 2017).

## 6. Limitations

- The estimand-vs-refit distinction (§3.1) means Phase-2 tables (measured against the finite-sample refit) and Phase-3 tables (population estimand) are not directly comparable — both are internally consistent; the methods note's capstone conclusion ("do not naively compose") stands, now with its mechanism fully explained.
- `Var(m̂)` in the joint estimator is treated as fixed in the CI (WLS noise is second-order at K=36 with 7 parameters, but not zero). Coverage 80–95% observed; below-nominal cases are in the sparsest regimes.
- Monotonicity of the joint calibration curve is exact in the no-noise limit and empirically established at finite n; a full proof is deferred (Prop. 4).
- Bounds assume monotone under-ascertainment (no false positives) and a correct envelope `delta_max`; both are stated, checkable-in-principle assumptions, not verifiable from the observed data alone.
- All results use the repository's fixed 36-stratum design; K-generalization is deferred (roadmap).

## 7. Reproduction

```bash
cd python && pip install -r requirements.txt
python ../scripts/validation/phase3_causal_demo.py   # headline table + figure (~35s)
pytest -q tests/test_causal.py                        # 6 regression tests (~10s)
pytest -q                                             # full suite: 33 tests
```

New code: `python/imaihda_sim/causal.py` (estimator + bounds + globally-verified optimizers), `python/tests/test_causal.py`, `scripts/validation/phase3_causal_demo.py`, `figures/phase3_causal.{png,csv}`.

## 8. References

1. Evans CR, Williams DR, Onnela J-P, Subramanian SV. A multilevel approach to modeling health inequalities at the intersection of multiple social identities. *SSM - Population Health*. 2018;6:149–157. doi:10.1016/j.ssmph.2018.08.005
2. VanderWeele TJ, Ding P. Sensitivity Analysis in Observational Research: Introducing the E-Value. *Annals of Internal Medicine*. 2017;167(4):268–274. doi:10.7326/M16-2607

*The partial-identification framing follows the spirit of the bounds literature in econometrics (Manski) and the graphical missing-data identification literature (Mohan & Pearl); exact citations to be verified in a dedicated literature pass before manuscript submission — deliberately not cited by DOI here to avoid an unverified reference.*

---

*No real data. Synthetic simulation only — this document makes no empirical claim about any population.*
