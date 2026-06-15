# Methods — core equations

> Canonical forms of the five equations underpinning the analysis, each with a one-line
> plain-language gloss and the assumption it rests on. Notation: `A` = digoxin exposure (1/0),
> `Y` = 6-month readmission (1/0), `W` = admission-time confounders, `Y^a` = potential outcome
> under exposure `a`. These match the estimator implemented in `src/`.

## 1. Target estimand — risk difference (and risk ratio)

$$\Delta \;=\; \mathbb{E}\!\left[Y^{a=1}\right] - \mathbb{E}\!\left[Y^{a=0}\right],
\qquad
\mathrm{RR} \;=\; \frac{\mathbb{E}\!\left[Y^{a=1}\right]}{\mathbb{E}\!\left[Y^{a=0}\right]}.$$

**Gloss.** The average causal contrast in 6-month readmission risk had *everyone* versus *no one* received digoxin.
**Rests on.** Well-defined potential outcomes — consistency / SUTVA (one version of treatment; no interference).

## 2. Identification — the g-formula

$$\mathbb{E}\!\left[Y^{a}\right] \;=\; \mathbb{E}_{W}\!\left[\,\mathbb{E}\!\left[Y \mid A=a,\, W\right]\right]
\;=\; \int \mathbb{E}\!\left[Y \mid A=a,\, W=w\right]\, dF_{W}(w).$$

**Gloss.** Replaces the unobservable counterfactual mean with the observed outcome regression, averaged (standardized) over the covariate distribution.
**Rests on.** Conditional exchangeability ($Y^a \perp A \mid W$, i.e. no unmeasured confounding given $W$), positivity, and consistency.

## 3. Propensity score

$$e(W) \;=\; \Pr\!\left(A=1 \mid W\right).$$

**Gloss.** The probability of receiving digoxin given baseline covariates; a balancing score, so conditioning on $e(W)$ alone removes measured confounding.
**Rests on.** Positivity — $0 < e(W) < 1$ for all $W$ (every patient could plausibly have gone either way).

## 4. AIPW (augmented inverse-probability weighting) estimator

With outcome model $m_a(W)=\mathbb{E}[Y\mid A=a,W]$ and propensity $\pi_a(W)=\Pr(A=a\mid W)$, the efficient influence function and estimator are

$$\varphi_a(O) \;=\; \frac{\mathbb{1}(A=a)}{\pi_a(W)}\bigl(Y - m_a(W)\bigr) \;+\; m_a(W) \;-\; \mu_a,
\qquad
\hat{\mu}_a \;=\; \frac{1}{n}\sum_{i=1}^{n}\left[\frac{\mathbb{1}(A_i=a)}{\pi_a(W_i)}\bigl(Y_i - m_a(W_i)\bigr) + m_a(W_i)\right],$$

$$\hat{\Delta}_{\mathrm{AIPW}} \;=\; \hat{\mu}_1 - \hat{\mu}_0,
\qquad
\widehat{\mathrm{Var}}(\hat{\Delta}_{\mathrm{AIPW}}) \;=\; \frac{1}{n}\,\widehat{\mathrm{Var}}\!\left(\varphi_1(O) - \varphi_0(O)\right).$$

**Gloss.** Doubly robust: consistent if *either* the propensity *or* the outcome model is correct; the influence function $\varphi$ supplies model-agnostic standard errors and 95% CIs.
**Rests on.** Exchangeability + positivity, and at least one nuisance model correctly specified — here both $\pi_a$ and $m_a$ are cross-fitted Super Learners to relax that.

## 5. E-value

$$\mathrm{E\text{-}value} \;=\; \mathrm{RR} + \sqrt{\mathrm{RR}\,\bigl(\mathrm{RR}-1\bigr)}
\qquad (\text{for } \mathrm{RR}\ge 1;\ \text{apply to } 1/\mathrm{RR} \text{ if } \mathrm{RR}<1).$$

**Gloss.** The minimum strength of association (on the risk-ratio scale) that an unmeasured confounder would need with *both* digoxin use and readmission — beyond all measured covariates — to fully explain away the observed effect.
**Rests on.** No additional modeling assumptions; it *bounds* the impact of residual confounding (e.g. unmeasured atrial fibrillation) rather than testing for it.
