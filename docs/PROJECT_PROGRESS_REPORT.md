# Project Progress Report — `imaihda` Methodological Companion to the PhD

*Status report mapping the repository's methodological work onto the PhD exposé (WZB Doctoral Position 400). No real data anywhere in this repository — synthetic simulation only. Prepared 2026-07.*

---

## 1. What the PhD needs, methodologically

The exposé asks whether the additive-dominant intersectional pattern of late-life multimorbidity (high PCV, low VPC — as in the SHARE type-2-diabetes finding, VPC 4.3% → 0.3%, PCV ≈ 92%) **generalizes** across high- and middle-income aging cohorts (HRS, ELSA, SHARE, CHARLS, with CHARLS as the middle-income anchor). It commits explicitly to *not* mistaking, for genuine structural difference, any of four artefacts:

> **event prevalence · sparse strata · selective attrition · unequal diagnosis (differential detection)**

and to delivering "a transparent reproducibility component documenting **when VPC and PCV can and cannot support cross-cohort comparison**."

This repository is that reproducibility companion. It does **not** compete with CRAN `MAIHDA` (which supplies the GLMM estimation the empirical chapters will use); it supplies the **diagnostic and sensitivity methodology that CRAN `MAIHDA` does not address at all** — measurement error, small-sample bias in the fast diagnostics, identification under differential detection, and cross-cohort comparability.

## 2. Coverage of the four artefact channels — now complete

| Artefact channel (exposé) | Methodology delivered | Where | In CRAN MAIHDA? |
|---|---|---|---|
| **Differential detection / diagnosis** | `correct_detection_bias()`, sensitivity bounds, tipping-point (E-value analogue); causal identification theory (point-ID vs. set-ID); `joint_calibrated_vpc()`; sharp partial-ID bounds | `detection.py`, `causal.py` | **No** |
| **Sparse strata** | `sparse_strata_vpc()` (calibration bias-correction + CI); K-generalization {8,36,108}; truncation-floor CI fix (`ci_method="test_inversion"`) | `uncertainty.py`, `robustness.py` | **No** |
| **Event prevalence** | Rare-prevalence × sparsity robustness study; prevalence channel in the decomposition | `robustness.py`, `cohort.py` | **No** |
| **Selective attrition** | `simulate_selective_attrition()`; attrition channel in the decomposition | `robustness.py`, `cohort.py` | **No** |
| **All four jointly + "when can VPC/PCV support comparison"** | `cross_cohort_decomposition()` — Shapley attribution of a HIC–MIC gap to the four artefacts vs. genuine structure, with an inferential "distinguishable from zero?" verdict | `cohort.py` | **No** |

Every capability above was built, executed, and regression-tested; the **54-test** Python suite runs in CI on every push.

## 3. The methodological arc (four phases)

1. **Phase 1–2 — detection & sparsity corrections + robustness study.** Two quantitative-bias-analysis tools with self-validating recovery tests, then a systematic study of how each behaves when its assumptions are violated, with formal non-inferiority / breakdown criteria (`METHODS_NOTE_ROBUSTNESS.md`). Honest negative findings retained (e.g. a *wrong functional form* for the detection score, or outcome-dependent detection above a threshold, genuinely breaks the correction; sparsity inflates the PCV toward a spurious "purely additive" reading).

2. **Phase 3 — causal identification.** Reframes detection bias as a missing-data identification problem: **point-identified** under covariate-only detection, **provably only set-identified** under outcome-dependent detection. Delivers `joint_calibrated_vpc()` (fixes the failure of naively composing the detection and sparsity corrections) and `vpc_partial_bounds()` (sharp, globally-verified bounds valid *even where point correction provably fails*). (`PHASE3_CAUSAL_IDENTIFICATION.md`)

3. **Phase 4 — cross-cohort decomposition (the capstone).** The direct answer to the exposé's central question: a Shapley decomposition of a raw HIC–MIC VPC/PCV gap into the four artefact channels plus a structural residual, self-validated against a known-by-construction true structural gap. Demonstrates all three regimes the PhD must distinguish — a gap that is 100% artefact (structure share exactly 0), a gap that is genuine structure, and — most importantly — a genuine structural difference that a naïve VPC comparison points the *wrong way* on because artefacts mask it. (`PHASE4_CROSS_COHORT_DECOMPOSITION.md`)

## 4. How this maps onto the PhD's subprojects

- **Subproject 1 (incident multimorbidity across cohorts).** Before any HIC–MIC VPC/PCV comparison is interpreted, the decomposition pre-specifies how large a spurious gap the cohorts' known prevalence/sparsity/attrition/detection differences could manufacture — turning "is this gap real?" into a quantified, reportable sensitivity statement. The corrections (`sparse_strata_vpc`, `correct_detection_bias`, `joint_calibrated_vpc`) provide the adjusted diagnostics to accompany the primary GLMM results.
- **Subproject 2 (diagnostic inequality / access as misclassification).** The entire detection-bias line (identification theory + correction + bounds) *is* the formal treatment of "diagnostic access as a source of outcome misclassification and inequality" the exposé names as a contribution.
- **Subproject 3 (biological aging, conditional).** Out of scope for this companion, consistent with the exposé's conditional framing.
- **Reproducibility & analytic safeguards.** The exposé promises a harmonization/design audit over "stratum definitions, event prevalence, minimum and effective stratum size, missingness, attrition, and measurement pathways." Each of those now has a corresponding synthetic diagnostic here, so the audit can be *pre-registered against simulated ground truth* before touching participant data.

## 5. What remains (honestly deferred)

- **Bridge to real data:** the decomposition currently reasons over *specified* cohort profiles; estimating each cohort's channel profile (especially detection and attrition) from validation sub-studies or design knowledge is the concrete next step, and would let `joint_calibrated_vpc()` replace the raw VPC inside the decomposition's value function (removing the structural underestimate seen under extreme nuisance).
- **PCV partial-identification bounds** (a ratio of two correlated variances — a genuinely harder optimization than the VPC bounds already built).
- **Literature-verification pass** for the Shapley-decomposition and Manski/Mohan–Pearl partial-identification framings (flagged, not yet cited by DOI — no fabricated references).
- **R mirror** of the Phase-3/Phase-4 modules (`causal.py`, `cohort.py`), and a full monotonicity proof for the joint-calibration curve (currently a stated conjecture with exact-in-the-limit + empirical support).
- **Design generalization:** other stratum designs (K ≠ 36) and structural forms beyond additive-plus-spike.

## 6. Bottom line

The repository now covers **all four artefact channels the exposé names**, plus the identification theory that says when each is correctable, plus a single capstone tool that answers the PhD's literal central question — *when can a raw HIC–MIC VPC/PCV difference be trusted as genuine structure?* — with a self-validating, inferential decomposition. None of this duplicates CRAN `MAIHDA`; all of it addresses methodology the intersectional-MAIHDA literature (per the project's PubMed scans) has left unaddressed. **54/54 tests pass; every reported number is regenerated by scripts in `scripts/validation/`.**

---

*No real data. Synthetic simulation only — this document makes no empirical claim about any population.*
