# Results summary — In-hospital digoxin and 6-month readmission

*Target-trial emulation, doubly-robust AIPW with cross-fitted Super-Learner. PhysioNet Zigong HF cohort, n = 2,008.*

## Headline — a convergence on confounding

**The headline result is methodological: three independent lines of evidence converge on the same conclusion — an apparent digoxin "harm" signal for 6-month readmission is residual confounding by indication, not a causal effect.** The digoxin estimate is the demonstration vehicle that makes the confounding visible.

The estimate: after doubly-robust adjustment, in-hospital digoxin was associated with a **higher** risk of 6-month all-cause readmission — risk ratio **1.24 (95% CI 1.11–1.40)**, risk difference **+8.4 percentage points (95% CI +3.9 to +12.8)**; adjusted risk 42.7% (digoxin) vs. 34.4% (no digoxin).

The three converging signals:

1. **Adjustment barely moves the estimate** (naive 1.26 → adjusted 1.24): measured confounders explain almost none of the crude association.
2. **The direction reverses the randomized DIG trial** (which *reduced* HF hospitalization, RR ≈ 0.72): an observational estimate opposite to the RCT is a confounding flag.
3. **The E-value is small** (1.79 point / 1.45 CI bound): a single unmeasured confounder of modest strength erases the result — and **atrial fibrillation**, the dominant indication for digoxin and an independent readmission driver, is **absent from the dataset** (Step 0: not measurable → carried by the E-value).

**Numerical verification.** E-value(1.2429) = 1.2429 + √(1.2429×0.2429) = 1.79; E-value(1.1066) = 1.1066 + √(1.1066×0.1066) = 1.45. Machine- and hand-checked, because a wrong E-value would undermine the entire sensitivity argument.

## Primary estimate (subdistribution / CIF view)

| Quantity | Estimate | 95% CI |
|---|---|---|
| Risk, digoxin | 0.427 | — |
| Risk, no digoxin | 0.344 | — |
| **Risk difference** | **+0.084** | +0.039, +0.128 |
| **Risk ratio** | **1.24** | 1.11, 1.40 |

Estimator: 5-fold cross-fitted AIPW; propensity and outcome nuisances each a stacked ensemble (penalized logistic + gradient boosting + random forest); influence-function-based CIs.

## Diagnostics

- **Positivity is excellent.** Cross-fitted propensity scores span **[0.14, 0.81]**; **0%** of units fall outside [0.05, 0.95]; **no trimming was required**.
- **Balance after weighting is excellent.** Maximum absolute standardized mean difference across the 14 confounders fell to **0.046** (all < 0.10) after IPTW weighting (see `love_plot.png`, `balance_table.csv`).
- Overlap histogram (`ps_overlap.png`) shows substantial common support in both arms.

## Sensitivity analyses (pre-specified)

| Analysis | RR | 95% CI | Reading |
|---|---|---|---|
| Naive (unadjusted) | 1.26 | 1.13, 1.41 | Adjustment barely moves the estimate |
| **AIPW adjusted (primary)** | **1.24** | **1.11, 1.40** | Primary |
| Cause-specific (exclude 57 deaths) | 1.21 | 1.08, 1.36 | **Agrees** with subdistribution |
| Any cardiac glycoside | 1.21 | 1.06, 1.39 | Consistent |
| DIG benchmark (HF hosp., RCT) | 0.72 | 0.66, 0.79 | External; opposite direction |

**(i) E-value.** For the point estimate, **E-value = 1.79**; for the CI bound nearest the null, **1.45**. An unmeasured confounder would need to be associated with both digoxin use and readmission by a risk ratio of **≥ 1.79** (point) / **≥ 1.45** (CI) — *above and beyond* every measured covariate — to fully explain away the association. Atrial fibrillation is a plausible candidate of roughly this strength (it is both a strong indication for digoxin and an independent readmission driver), so the finding is **not robust** to realistic unmeasured confounding.

**(ii) Cause-specific vs. subdistribution.** The two views agree (RR 1.21 vs. 1.24; overlapping CIs). Because digoxin is mortality-neutral here, agreement is expected; the absence of divergence argues against differential mortality-driven selection.

**(iii) Any-glycoside.** Adding deslanoside leaves the estimate essentially unchanged (RR 1.21).

**(iv) Benchmark contrast.** The naive and adjusted estimates are nearly identical (1.26 → 1.24): the measured confounders explain almost none of the crude association. Combined with the direction reversal vs. the DIG RCT, this pattern is the signature of **confounding by indication** that the measured covariates do not capture — most plausibly atrial fibrillation and the lack of medication-timing/dose data. We report what the data show and do **not** tune toward the benchmark.

**(v) Mortality (descriptive only).** Crude 6-month mortality was 1.9% (digoxin) vs. 3.8% (no digoxin); 57 deaths total. Presented as a cumulative-incidence figure (`mortality_cif.png`) for context only — **no causal claim**.

## Limitations (stated plainly)

1. **Single-center** Chinese cohort (Zigong Fourth People's Hospital); external validity is limited.
2. **Atrial fibrillation is unmeasured** — no baseline rhythm field exists in the data. As the dominant indication for digoxin, AF is the most likely residual confounder; the E-value (1.79) is within the range a single AF-strength confounder could produce. This is the principal threat to the finding.
3. **LVEF is 68% missing** and was excluded; HF phenotype (HFrEF vs. HFpEF) is therefore only partly captured (NYHA, Killip, BNP).
4. **No medication timing, dose, or duration.** Exposure is an "ever during admission" binary; immortal-time and dose-response structure cannot be modeled.
5. **57 deaths** limit power for the competing-risk contrast; the mortality analysis is descriptive.
6. Single imputation (with missingness indicators) was used; **MICE** is the planned refinement.

## Bottom line

In this single-center cohort, in-hospital digoxin is associated with **~24% higher 6-month all-cause readmission**, but the association is **almost entirely unexplained by measured confounders, reverses the direction of the DIG RCT, and is fragile to a single unmeasured confounder of atrial-fibrillation strength (E-value 1.79)**. The most defensible reading is **residual confounding by indication**, not a causal harm. The result is best treated as a hypothesis-generating, methods-demonstration analysis pending AF-aware data and medication-timing detail.
