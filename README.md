# heart-failure-causal

**Pre-registered target-trial emulation of in-hospital digoxin and 6-month readmission in hospitalized heart-failure patients.**

Doubly-robust causal analysis (cross-fitted AIPW with a stacked Super-Learner) on the PhysioNet Zigong heart-failure cohort (n = 2,008). Protocol, DAG, analysis code, diagnostics, pre-specified sensitivity analyses, and results.

> **Status: Analysis complete — preprint draft.**

## Headline finding — a convergence on confounding

**The finding is methodological, not pharmacological: three independent lines of evidence converge to show that an apparent digoxin "harm" signal is residual confounding by indication.** The digoxin–readmission estimate is the *demonstration vehicle*, not the claim.

The estimate itself: after doubly-robust adjustment, in-hospital digoxin is associated with **higher** 6-month all-cause readmission — RR **1.24 (95% CI 1.11–1.40)**, RD **+8.4 pp (95% CI +3.9 to +12.8)**; adjusted risk 42.7% (digoxin) vs. 34.4% (no digoxin).

Three signals say this is confounding, not causation:

1. **Adjustment moves almost nothing** — naive 1.26 -> adjusted 1.24. The measured confounders explain virtually none of the crude association, the fingerprint of an open backdoor they do not capture.
2. **It reverses the RCT** — the randomized DIG trial *reduced* HF hospitalization (RR ~ 0.72). An observational estimate pointing the opposite way is a red flag, not a discovery.
3. **It is fragile by exactly the right amount** — **E-value 1.79** (point) / **1.45** (CI bound). A single unmeasured confounder of that modest strength explains the result away, and **atrial fibrillation** — the dominant indication for digoxin and an independent readmission driver — is **not recorded anywhere in this dataset** (Step 0 verdict: not measurable -> carried by the E-value).

The most defensible reading is **residual confounding by indication**, not causal harm. This is a methods-demonstration / hypothesis-generating analysis. See [`docs/results_summary.md`](docs/results_summary.md).

> **Numerical checks (verified independently).** E-value(RR 1.2429) = 1.2429 + sqrt(1.2429 x 0.2429) = **1.79**; E-value(CI bound 1.1066) = 1.1066 + sqrt(1.1066 x 0.1066) = **1.45**. Both machine- and hand-computed — a wrong E-value would quietly undermine the entire sensitivity argument, so it is checked, not assumed.

## Understanding the causal approach

A plain-language walkthrough of how the causal-inference approach works in this project — the five moves from clinical question to defensible answer — is in **[`docs/causal_reasoning.md`](docs/causal_reasoning.md)** (with a live Mermaid diagram of the reasoning).

## Honest caveats

- **Single-center** Chinese cohort — limited external validity.
- **Atrial fibrillation unmeasured** (no baseline rhythm field exists) — the principal residual-confounding threat, carried explicitly by the E-value. See the Step 0 verdict in [`docs/protocol.md`](docs/protocol.md).
- **LVEF 68% missing** — excluded; HF phenotype only partly captured.
- **No medication timing/dose/duration** — exposure is an "ever-during-admission" binary.
- **57 deaths** limit the competing-risk contrast; mortality is reported descriptively only.
- Single imputation with missingness indicators; **MICE** is the planned refinement.

## Methods at a glance

| Element | Choice |
|---|---|
| Design | Target-trial emulation |
| Exposure | Oral digoxin (`Digoxin tablet`); secondary: any cardiac glycoside |
| Outcome | 6-month all-cause readmission (subdistribution primary; cause-specific secondary) |
| Confounders | Admission-time only: age, sex, NYHA, Killip, BNP, creatinine, eGFR, Charlson, vitals, K/Na, prior admissions |
| Estimator | Cross-fitted AIPW (5-fold); Super-Learner (penalized logistic + GBM + RF) for both nuisances; influence-function CIs |
| Diagnostics | PS overlap, positivity report, SMD love plot |
| Sensitivity | E-value; cause-specific vs. subdistribution; any-glycoside; DIG benchmark forest; descriptive mortality CIF |

## Reasoning & methods detail

Three artifacts document the thinking behind the analysis — read them in this order: the intuition, then the scoping logic, then the formal math.

- **[`docs/causal_reasoning.md`](docs/causal_reasoning.md)** — *how the causal logic works.* A plain-language walkthrough of the five causal-inference moves (target trial → DAG → doubly-robust estimation → diagnostics → stress-testing → verdict), with a live Mermaid diagram of the reasoning.
- **[`docs/decision_tree.mmd`](docs/decision_tree.mmd)** (renders natively on GitHub) and **[`output/figures/decision_tree.png`](output/figures/decision_tree.png)** — *how the project was scoped.* Every decision (study type → exposure → drug → time-to-event method → outcome → EF phenotype → estimation), with the data fact or principle that forced each choice and rejected branches greyed.
- **[`docs/methods_equations.md`](docs/methods_equations.md)** — *the formal math.* The five core equations (estimand Δ, g-formula identification, propensity score, AIPW estimating equation, E-value), each with a gloss and its assumption.

## ⚠️ Data access & licensing

The analysis uses PhysioNet's **"Hospitalized patients with heart failure (Zigong)"** dataset, governed by the **PhysioNet Restricted Health Data License 1.5.0**. **Raw data are NOT in this repository and must never be committed** (license Term 3). To reproduce:

1. Create a credentialed PhysioNet account and complete the required training (CITI "Data or Specimens Only Research").
2. Sign the Data Use Agreement: https://physionet.org/content/heart-failure-zigong/
3. Download the files and place them at `data/raw/physionet/` (this path is git-ignored).

This repo contains only **code and aggregate, non-identifying results**, consistent with license Term 7 (encouraging publication code to be shared openly).

## Reproduce

```bash
pip install -r requirements.txt
cd src
python 00_build_dataset.py        # builds data/derived/analysis_dataset.csv (git-ignored)
python fit_one.py primary         # cross-fitted AIPW -> /tmp/aipw_primary.pkl
python fit_one.py causespecific
python fit_one.py glycoside
python assemble.py                # diagnostics, sensitivity, figures, tables, results.json
python make_dag.py                # DAG + conceptual figures
```

## Repository layout

```
docs/      protocol.md (incl. Step 0 AF verdict), results_summary.md,
           causal_reasoning.md, methods_equations.md, decision_tree.mmd
src/       00_build_dataset.py, estimators.py, fit_one.py, assemble.py,
           make_dag.py, make_decision_tree.py
output/
  figures/ dag.png, conceptual_target_trial.png, ps_overlap.png, love_plot.png,
           forest_estimates.png, mortality_cif.png, decision_tree.png
  tables/  balance_table.csv, primary_estimates.csv, sensitivity.csv, positivity_report.csv
  results.json
data/      raw/ and derived/  (BOTH git-ignored — never committed)
```

## Citation

Dataset: Zhang Z, Cao L, Chen R, et al. *Hospitalized patients with heart failure: integrating electronic healthcare records and external outcome data.* PhysioNet (2021). Companion paper: *Scientific Data* 2021;8:46. PMID 33547290. doi:10.1038/s41597-021-00835-9.

External benchmark: The Digitalis Investigation Group. *The effect of digoxin on mortality and morbidity in patients with heart failure.* N Engl J Med 1997;336:525–533.

## License

Code: **MIT** (see `LICENSE`). Data: **PhysioNet Restricted Health Data License 1.5.0** — not redistributed here.
