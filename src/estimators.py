"""
estimators.py — cross-fitted AIPW (DML-style) with a stacked Super-Learner ensemble.

Nuisance learners (both propensity g(W) and outcome m(A,W)):
  stacked ensemble of penalized logistic regression + gradient boosting + random forest
  (sklearn StackingClassifier, logistic meta-learner, internal CV).

Estimands: risk difference (RD) and risk ratio (RR) for binary Y.
Inference: efficient influence function (EIF) -> robust SE and 95% CIs.
Cross-fitting: K-fold; nuisance fit on training folds, evaluated out-of-fold.
"""
import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import (GradientBoostingClassifier, RandomForestClassifier,
                              StackingClassifier)
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline


def make_superlearner(seed=42):
    """Stacked ensemble: penalized logistic + GBM + RF, logistic meta-learner."""
    base = [
        ("penlogit", Pipeline([("sc", StandardScaler()),
                               ("lr", LogisticRegression(penalty="l2", C=1.0,
                                                         max_iter=2000, solver="lbfgs"))])),
        ("gbm", GradientBoostingClassifier(n_estimators=120, max_depth=3,
                                           learning_rate=0.05, subsample=0.9,
                                           random_state=seed)),
        ("rf",  RandomForestClassifier(n_estimators=250, max_depth=None,
                                       min_samples_leaf=10, n_jobs=2,
                                       random_state=seed)),
    ]
    return StackingClassifier(
        estimators=base,
        final_estimator=LogisticRegression(max_iter=2000),
        stack_method="predict_proba",
        cv=3, n_jobs=2,
    )


def _predict_proba1(model, X):
    p = model.predict_proba(X)
    # column for class "1"
    j = list(model.classes_).index(1)
    return p[:, j]


def crossfit_aipw(W, A, Y, n_splits=5, seed=42,
                  ps_trim=(0.02, 0.98), report_bounds=(0.05, 0.95)):
    """
    Returns dict with RD, RR, EIF-based SEs and 95% CIs, plus diagnostics
    (cross-fitted propensity scores, potential-outcome predictions, trim info).
    """
    W = np.asarray(W, float); A = np.asarray(A, int); Y = np.asarray(Y, int)
    n = len(Y)
    g = np.full(n, np.nan)        # P(A=1|W)
    mu1 = np.full(n, np.nan)      # E[Y|A=1,W]
    mu0 = np.full(n, np.nan)      # E[Y|A=0,W]

    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    for tr, te in skf.split(W, A):
        # propensity
        gm = make_superlearner(seed); gm.fit(W[tr], A[tr])
        g[te] = _predict_proba1(gm, W[te])
        # outcome model on [W, A]
        WA = np.column_stack([W, A])
        om = make_superlearner(seed); om.fit(WA[tr], Y[tr])
        WA1 = np.column_stack([W[te], np.ones(len(te))])
        WA0 = np.column_stack([W[te], np.zeros(len(te))])
        mu1[te] = _predict_proba1(om, WA1)
        mu0[te] = _predict_proba1(om, WA0)

    # positivity report (pre-trim)
    pos = {
        "ps_min": float(np.min(g)), "ps_max": float(np.max(g)),
        "pct_outside_report_bounds": float(np.mean((g < report_bounds[0]) |
                                                   (g > report_bounds[1])) * 100),
        "report_bounds": list(report_bounds),
        "trim_bounds": list(ps_trim),
    }
    keep = (g >= ps_trim[0]) & (g <= ps_trim[1])
    pos["n_trimmed"] = int((~keep).sum())
    pos["n_kept"] = int(keep.sum())

    gk, A_k, Y_k = g[keep], A[keep], Y[keep]
    mu1k, mu0k = mu1[keep], mu0[keep]
    nk = keep.sum()

    # AIPW (efficient influence function) pseudo-outcomes
    psi1 = mu1k + A_k * (Y_k - mu1k) / gk
    psi0 = mu0k + (1 - A_k) * (Y_k - mu0k) / (1 - gk)
    r1, r0 = psi1.mean(), psi0.mean()

    # RD
    rd = r1 - r0
    if_rd = (psi1 - psi0) - rd
    se_rd = np.sqrt(np.var(if_rd, ddof=1) / nk)

    # RR with delta-method on log scale (EIF)
    rr = r1 / r0
    if_logrr = (psi1 - r1) / r1 - (psi0 - r0) / r0
    se_logrr = np.sqrt(np.var(if_logrr, ddof=1) / nk)

    z = 1.959963985
    out = {
        "n": int(nk),
        "risk_treated": float(r1), "risk_control": float(r0),
        "RD": float(rd), "RD_se": float(se_rd),
        "RD_lcl": float(rd - z * se_rd), "RD_ucl": float(rd + z * se_rd),
        "RR": float(rr), "logRR_se": float(se_logrr),
        "RR_lcl": float(np.exp(np.log(rr) - z * se_logrr)),
        "RR_ucl": float(np.exp(np.log(rr) + z * se_logrr)),
        "positivity": pos,
        "ps": g, "mu1": mu1, "mu0": mu0, "keep": keep,
    }
    return out


def naive_rr_rd(A, Y):
    """Unadjusted association."""
    A = np.asarray(A, int); Y = np.asarray(Y, int)
    r1 = Y[A == 1].mean(); r0 = Y[A == 0].mean()
    n1 = (A == 1).sum(); n0 = (A == 0).sum()
    rd = r1 - r0
    se_rd = np.sqrt(r1*(1-r1)/n1 + r0*(1-r0)/n0)
    rr = r1 / r0
    se_logrr = np.sqrt((1-r1)/(r1*n1) + (1-r0)/(r0*n0))
    z = 1.959963985
    return {
        "risk_treated": float(r1), "risk_control": float(r0),
        "RD": float(rd), "RD_lcl": float(rd-z*se_rd), "RD_ucl": float(rd+z*se_rd),
        "RR": float(rr), "RR_lcl": float(np.exp(np.log(rr)-z*se_logrr)),
        "RR_ucl": float(np.exp(np.log(rr)+z*se_logrr)),
    }


def evalue_rr(rr, lcl=None, ucl=None):
    """VanderWeele & Ding E-value for a risk ratio and the CI bound nearest the null."""
    def _e(x):
        if x < 1:
            x = 1.0 / x
        return x + np.sqrt(x * (x - 1.0))
    out = {"evalue_point": float(_e(rr))}
    if lcl is not None and ucl is not None:
        if lcl <= 1 <= ucl:
            out["evalue_ci"] = 1.0  # CI crosses null
        else:
            bound = lcl if rr > 1 else ucl
            out["evalue_ci"] = float(_e(bound))
    return out
