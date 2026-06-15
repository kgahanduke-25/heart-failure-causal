# Pre-registered protocol — Target-trial emulation of in-hospital digoxin and 6-month readmission

**Data source.** Zhang Z, Cao L, Chen R, et al. *Hospitalized patients with heart failure: integrating electronic healthcare records and external outcome data.* PhysioNet (Zigong Fourth People's Hospital, Sichuan, China), v1.3. n = 2,008 index hospitalizations. Restricted Health Data License — raw data not redistributed.

## Causal question
Among patients hospitalized for heart failure, what is the effect of **in-hospital digoxin** versus **no digoxin** on **all-cause readmission within 6 months**?

## Estimands
- **Risk difference (RD)** and **risk ratio (RR)** for 6-month all-cause readmission.
- Average treatment effect in the population (marginal), via doubly-robust standardization.

## Exposure
- **Primary A:** oral digoxin — `Drug_name == "Digoxin tablet"` in `dat_md.csv`, collapsed to a per-patient binary.
- **Secondary:** any cardiac glycoside — `Digoxin tablet` **OR** `Deslanoside injection`.

## Outcome
- **Y:** `re.admission.within.6.months` (1/0).
- **Primary (subdistribution / CIF) view:** the 57 patients who died within 6 months **without** readmission are retained and coded Y = 0 (death treated as a competing event that precludes the outcome).
- **Secondary (cause-specific) view:** those 57 are excluded and the estimate recomputed. Because digoxin is expected to be mortality-neutral, the two views should agree; divergence would flag residual confounding or differential censoring.

## Confounders W (admission-time only — no post-treatment variables)
age (decade midpoint from `ageCat`), sex, prior admissions (`visit.times`), HF severity (NYHA class, Killip grade, BNP), renal function (creatinine, eGFR), Charlson comorbidity index (`CCI.score`), admission vitals (heart rate `pulse`, SBP, DBP), and electrolytes (potassium, sodium). Implausible zero vitals were set to missing before imputation.

## Step 0 — Atrial-fibrillation confounder: AUTO-DECISION

**Verdict: AF EXCLUDED from W — not cleanly measurable at baseline; carried by the E-value.**

A programmatic scan of all 166 dictionary variables and all `dat.csv` column names for `fibrill*`, `arrhythm*`, `rhythm`, `atrial`, and `flutter` returned **no matching field**. The dataset contains no baseline atrial-fibrillation / arrhythmia / rhythm indicator (and no ECG rhythm variable). AF is therefore **absent**, not merely incompletely measured. Per the pre-specified rule, AF is excluded from W, and its potential confounding is addressed quantitatively through the **E-value** in the sensitivity analysis. This is the single most important unmeasured confounder: AF is both a classic indication for digoxin (rate control) and an independent driver of readmission, so an open backdoor through AF is expected and is explicitly carried.

**LVEF EXCLUDED** by pre-specification (68.4% missing — not informatively recoverable).

## Missing data
Median (continuous) / mode (binary-ordinal) single imputation, plus a per-variable **missingness-indicator** column for every imputed variable. All baseline confounders have < 4% missingness. **Multiple imputation (MICE)** is the planned refinement.

## Estimator (Step 2)
Doubly-robust **AIPW** with **cross-fitting** (5-fold). Both nuisance functions — the propensity score g(W) = P(A=1|W) and the outcome regression m(A,W) = E[Y|A,W] — are estimated with a **stacked Super-Learner ensemble**: penalized logistic regression + gradient boosting + random forest, combined by a logistic meta-learner. Inference uses the **efficient influence function** (robust SE; 95% CIs on the RD and on log-RR).

## Diagnostics (Step 3)
Propensity-score overlap histogram by arm; positivity report (min/max PS, % outside [0.05, 0.95], trimming rule [0.02, 0.98]); standardized-mean-difference love plot before vs. after IPTW weighting.

## Pre-specified sensitivity analyses (Step 4)
1. **E-value** for the RR point estimate and its CI bound nearest the null.
2. **Cause-specific vs. subdistribution** agreement.
3. **Any-glycoside** exposure rerun.
4. **Naive vs. adjusted vs. DIG-trial benchmark** (forest). The DIG benchmark is external context only; estimates are reported as the data show them — not tuned toward the benchmark.
5. **Mortality:** descriptive cumulative incidence only — no causal estimate.

## Analysis-set provenance
Raw data are read only from `data/raw/` (git-ignored). The derived analysis dataset is written to `data/derived/` (git-ignored). No patient-level data is committed to version control.
