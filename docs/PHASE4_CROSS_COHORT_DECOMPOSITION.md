# Attributing Cross-Cohort VPC/PCV Differences to Artefact Channels vs. Genuine Structure: A Shapley Decomposition

*Phase-4 methods report. No real data — synthetic simulation only. Every number below is reproduced by `scripts/validation/cross_cohort_demo.py`; none is hand-transcribed.*

---

## 0. The question this answers

The PhD exposé (WZB Doctoral Position 400, *"Do Intersectional Patterns of Late-Life Health Inequality Generalize Across High- and Middle-Income Aging Contexts?"*) states its central methodological worry precisely: a HIC–MIC comparison of the variance partition coefficient (VPC) or proportional change in variance (PCV) "could be informative, but only if it does not mistake changes in **event prevalence, sparse strata, selective attrition, or unequal diagnosis** for differences in intersectional structure," and names as a core contribution "a transparent reproducibility component documenting **when VPC and PCV can and cannot support cross-cohort comparison**."

Phases 1–3 of this repository built and stress-tested corrections for three of those four artefact channels (prevalence, sparse strata, differential detection) and the causal-identification theory behind them. **Two gaps remained**: (1) selective **attrition** — the fourth channel — had no generator or characterization; (2) there was no single tool that answers the exposé's actual question, which is not "correct one cohort's VPC" but "**given a raw HIC–MIC VPC gap, how much of it could each artefact channel account for, and is the residual structural difference distinguishable from zero?**" Phase 4 closes both.

## 1. The tool

`cross_cohort_decomposition(spec_a, spec_b, metric="vpc"|"pcv")` (`python/imaihda_sim/cohort.py`) takes two cohort specifications on the shared 36-stratum design — each carrying the four artefact-channel parameters (prevalence shift, sparse allocation, detection strength, attrition strength) plus the true between-stratum structure — and decomposes the observed gap `VPC_B − VPC_A` into **five additive contributions**: one per artefact channel plus a **structural residual**.

**Why Shapley.** Flipping the channels one at a time from cohort A's value to cohort B's traces a path from `VPC_A` to `VPC_B`, but the channels are non-orthogonal (sparsity and detection interact, prevalence modulates both), so a single sequential (Oaxaca–Blinder) path gives components that depend on the arbitrary ordering. The **Shapley value** (Shorrocks' 2013 unified framework for decomposition) is the unique attribution that is (i) *symmetric* across channels, (ii) *exactly additive* — the five shares sum to the observed gap by construction — and (iii) *order-independent* (it averages each channel's marginal contribution over all `5! = 120` orderings, computed exactly here via the `2⁵ = 32`-subset value cache). No equivalent exists in CRAN `MAIHDA`, which offers no cross-cohort comparison machinery at all.

**Inference.** A bootstrap over independent simulation streams yields a confidence interval on the **structural** share, and the boolean `structure_distinguishable_from_zero` — the direct answer to "could nuisance differences alone explain the observed gap?"

**Scope.** This is a *design/sensitivity* tool exactly as the exposé frames the repository ("to pre-specify diagnostics and sensitivity scenarios… not a new estimator or a substitute for GLMM-based inference"). It reasons over *hypothesized* cohort profiles to say how large a spurious VPC gap a given set of nuisance differences could manufacture — not to estimate a real cohort's VPC.

## 2. Self-validation (the true structural gap is known by construction)

Three scenarios, `n_rep=12`, `n_boot=40` (`scripts/validation/cross_cohort_demo.py`):

| Scenario | VPC_A | VPC_B | observed gap | prevalence | sparsity | detection | attrition | **structure** | structure 95% CI | distinguishable? | true structural gap |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|:--:|---:|
| **null_structure** | 11.10 | 0.09 | −11.01 | −0.44 | −3.72 | −5.39 | −1.45 | **0.00** | [0.00, 0.00] | **No** | 0.00 |
| **pure_structure** | 10.12 | 25.01 | +14.89 | 0.00 | 0.00 | 0.00 | 0.00 | **14.89** | [12.61, 16.04] | **Yes** | 14.24 |
| **masked_structure** | 10.12 | 4.90 | −5.23 | −0.76 | −7.70 | −6.13 | −0.90 | **10.27** | [8.50, 12.70] | **Yes** | 14.24 |

![Cross-cohort VPC gap decomposition](figures/cross_cohort_decomposition.png)

Three findings, each answering the exposé directly:

1. **Null structure (the primary safeguard).** Two cohorts with *identical* intersectional structure but different prevalence, sparsity, detection, and attrition show an −11.0pp raw VPC gap — a MIC cohort that, at face value, looks dramatically *more additive* (VPC ≈ 0). The decomposition attributes **exactly 0.00pp to structure** (CI [0,0], not distinguishable) and the entire gap to the four artefacts (detection −5.4, sparsity −3.7, attrition −1.5, prevalence −0.4). This is the concrete demonstration that the exposé's fear is real and detectable: *a raw HIC–MIC VPC difference of this magnitude can be 100% artefact.*

2. **Pure structure (recovery).** With a genuine structural difference and no artefact difference, the structural share (14.89) equals the gap, the artefacts sum to 0.00, and it recovers the analytically-known true gap (14.24) within 0.7pp.

3. **Masked structure (the scientifically dangerous case).** A cohort with genuinely *higher* intersectional structure (true gap +14.24) but also more sparsity/detection/attrition shows a **negative** raw gap (−5.23) — at face value the MIC cohort looks *less* intersectional, the opposite of the truth. The decomposition unmasks a **positive, distinguishable structural share (+10.27, CI [8.5, 12.7])**. A naïve VPC comparison here would not merely be imprecise — it would point the wrong way. This is exactly the "greater discriminatory accuracy and larger residual intersectional component" the exposé hypothesizes for the middle-income anchor, and shows why it cannot be read off raw VPCs.

## 3. The attrition channel (`simulate_selective_attrition`)

The fourth artefact, previously absent. Retention is SES-patterned and mildly outcome-dependent —
`logit(P(retained)) = 2.0 − attrition_strength·(edu + wealth + 0.4·rural) − 0.3·attrition_strength·y_true` — so disadvantaged strata, and true cases within them, are differentially lost. Unlike detection (which zeroes the recorded outcome but keeps the person in the denominator), attrition removes the individual entirely, shrinking both a stratum's event count and its size — hence its distinct signature in the decomposition. `attrition_strength=0` reproduces `simulate_intersectional_data` exactly (regression-checked). In the scenarios above attrition is the smallest of the four artefact channels here (−0.9 to −1.5pp), but it is non-negligible and, crucially, *separately attributable*.

## 4. PCV decomposition and its honest limit

The same machinery runs on the PCV (`metric="pcv"`). It is reliable in mild regimes (pure-structure: structural share −48.3, CI [−51.4, −42.5], distinguishable) but the tool **flags itself unreliable** (`reliable=False`, `distinguishable=None`) under heavy nuisance, because the PCV is genuinely undefined when the null variance collapses to zero — the exact PCV-collapse pathology characterized in `METHODS_NOTE_ROBUSTNESS.md` §5.5. The decomposition does not paper over this; it refuses to report a spurious PCV attribution where the PCV itself is not identified.

## 5. Limitations

- **Simulation-based, not an estimator.** The decomposition reasons over hypothesized cohort profiles; applying it to real HRS/ELSA/SHARE/CHARLS data requires plugging in each cohort's estimated nuisance parameters (prevalence, stratum-size distribution, and — harder — plausible detection and attrition profiles from validation sub-studies or design knowledge). Turning those inputs into estimated (rather than specified) channel profiles is the bridge to the empirical PhD chapters.
- **Structural underestimate under extreme nuisance.** In the masked scenario the structural share (10.27) sits below the true gap (14.24): when other channels are at their MIC-heavy values, the VPC estimator is itself attenuated (the K=36 truncation-floor effect of Phase 3), so the structure channel's averaged marginal effect is damped. Plugging the Phase-3 `joint_calibrated_vpc()` estimator into the value function (instead of the raw fast VPC) is the natural fix and a concrete next step.
- **Fixed 36-stratum design and a single structural form.** Both cohorts share the 2×3×3×2 design and the additive-plus-two-spike structure; generalization to other designs reuses `simulate_k_strata`'s machinery but was not run here.
- **Outcome-dependent attrition/detection remain only set-identified** (Phase-3, Prop. 2): the decomposition treats the specified attrition/detection profiles as known; if they are themselves uncertain, the honest report is the `vpc_partial_bounds`-style envelope, not a point decomposition.

## 6. Reproduction

```bash
cd python && pip install -r requirements.txt
python ../scripts/validation/cross_cohort_demo.py   # table + figure (~20s)
pytest -q tests/test_cohort.py                       # 8 tests
```

## 7. References

1. Evans CR, Williams DR, Onnela J-P, Subramanian SV. A multilevel approach to modeling health inequalities at the intersection of multiple social identities. *SSM - Population Health*. 2018;6:149–157. doi:10.1016/j.ssmph.2018.08.005
2. O'Sullivan JL, Alonso-Perez E, Färber F, et al. Onset of Type 2 diabetes in adults aged 50 and older in Europe: an intersectional multilevel analysis of individual heterogeneity and discriminatory accuracy. *Diabetology & Metabolic Syndrome*. 2024;16:293. doi:10.1186/s13098-024-01533-3

*The Shapley-decomposition-of-inequality framing follows Shorrocks' unified attribution approach; the exact bibliographic reference is to be verified in a dedicated literature pass before manuscript submission — deliberately not cited by DOI here to avoid an unverified reference.*

---

*No real data. Synthetic simulation only — this document makes no empirical claim about any population.*
