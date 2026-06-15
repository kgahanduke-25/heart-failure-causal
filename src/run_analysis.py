"""
run_analysis.py — full pipeline: estimation (Step 2), diagnostics (Step 3),
sensitivity (Step 4). Writes output/tables/*.csv, output/figures/*.png, results.json.
Reads derived (gitignored) analysis dataset; reads raw only for mortality timing.
"""
import os, json, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estimators import crossfit_aipw, naive_rr_rd, evalue_rr, make_superlearner, _predict_proba1

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DERIV = os.path.join(ROOT, "data", "derived")
RAW   = os.path.join(ROOT, "data", "raw", "physionet")
FIG   = os.path.join(ROOT, "output", "figures")
TAB   = os.path.join(ROOT, "output", "tables")
os.makedirs(FIG, exist_ok=True); os.makedirs(TAB, exist_ok=True)
SEED = 42

df = pd.read_csv(os.path.join(DERIV, "analysis_dataset.csv"))
Wcols = ["age","male","visit_times","NYHA","Killip","BNP","creatinine","eGFR","CCI",
         "heart_rate","SBP","DBP","potassium","sodium",
         "BNP_miss","creatinine_miss","eGFR_miss","CCI_miss","heart_rate_miss",
         "SBP_miss","DBP_miss","potassium_miss","sodium_miss"]
LABELS = {"age":"Age","male":"Male","visit_times":"Prior admissions","NYHA":"NYHA class",
          "Killip":"Killip grade","BNP":"BNP","creatinine":"Creatinine","eGFR":"eGFR",
          "CCI":"Charlson index","heart_rate":"Heart rate","SBP":"SBP","DBP":"DBP",
          "potassium":"Potassium","sodium":"Sodium"}
CORE = [c for c in Wcols if not c.endswith("_miss")]

results = {}

# ---------------- Step 2: primary (digoxin, subdistribution view) ----------------
W = df[Wcols].values
print(">> primary AIPW (digoxin, subdistribution)...", flush=True)
prim = crossfit_aipw(W, df["A_digoxin"].values, df["Y"].values, n_splits=5, seed=SEED)
results["primary_digoxin_subdist"] = {k: v for k, v in prim.items()
                                      if k not in ("ps","mu1","mu0","keep")}

# ---------------- Step 4(ii): cause-specific (exclude 57 died-without-readmit) ----
print(">> cause-specific AIPW (exclude 57)...", flush=True)
cs = df[df["died_no_readmit"] == 0].reset_index(drop=True)
csres = crossfit_aipw(cs[Wcols].values, cs["A_digoxin"].values, cs["Y"].values,
                      n_splits=5, seed=SEED)
results["causespecific_digoxin"] = {k: v for k, v in csres.items()
                                    if k not in ("ps","mu1","mu0","keep")}

# ---------------- Step 4(iii): any cardiac glycoside ----------------
print(">> any-glycoside AIPW...", flush=True)
gly = crossfit_aipw(W, df["A_glycoside"].values, df["Y"].values, n_splits=5, seed=SEED)
results["secondary_anyglycoside"] = {k: v for k, v in gly.items()
                                     if k not in ("ps","mu1","mu0","keep")}

# ---------------- naive ----------------
results["naive_digoxin"] = naive_rr_rd(df["A_digoxin"].values, df["Y"].values)

# ---------------- E-values (Step 4 i) ----------------
results["evalue_primary"] = evalue_rr(prim["RR"], prim["RR_lcl"], prim["RR_ucl"])

# ---------------- DIG trial benchmark (external; published) ----------------
# DIG (NEJM 1997): all-cause mortality RR ~0.99 (neutral); HF-hospitalization RR 0.72 (0.66-0.79).
results["dig_benchmark"] = {
    "label": "DIG trial (NEJM 1997) — HF hospitalization",
    "RR": 0.72, "RR_lcl": 0.66, "RR_ucl": 0.79,
    "mortality_RR": 0.99, "note": "RCT, outpatients, reduced EF, sinus rhythm — different population/outcome window"
}

# ================= Step 3: DIAGNOSTICS =================
g = prim["ps"]; A = df["A_digoxin"].values
# (a) PS overlap
plt.figure(figsize=(7,4.5))
plt.hist(g[A==1], bins=30, alpha=0.6, label="Digoxin (treated)", color="#2b6cb0", density=True)
plt.hist(g[A==0], bins=30, alpha=0.6, label="No digoxin (control)", color="#c05621", density=True)
for b in prim["positivity"]["report_bounds"]:
    plt.axvline(b, ls="--", c="grey", lw=1)
plt.xlabel("Estimated propensity score  P(digoxin | W)"); plt.ylabel("Density")
plt.title("Propensity-score overlap by arm (cross-fitted)")
plt.legend(); plt.tight_layout(); plt.savefig(os.path.join(FIG,"ps_overlap.png"), dpi=150); plt.close()

# (b) SMD love plot: before vs after IPTW weighting (stabilized, trimmed to [0.02,0.98])
def smd(x, a, w=None):
    x=np.asarray(x,float); a=np.asarray(a,int)
    if w is None:
        m1,m0=x[a==1].mean(),x[a==0].mean()
        v1,v0=x[a==1].var(ddof=1),x[a==0].var(ddof=1)
    else:
        w1,w0=w[a==1],w[a==0]; x1,x0=x[a==1],x[a==0]
        m1=np.average(x1,weights=w1); m0=np.average(x0,weights=w0)
        v1=np.average((x1-m1)**2,weights=w1); v0=np.average((x0-m0)**2,weights=w0)
    sp=np.sqrt((v1+v0)/2); return 0.0 if sp==0 else (m1-m0)/sp

keep = prim["keep"]
gk=np.clip(g,0.02,0.98)
w_iptw = A/gk + (1-A)/(1-gk)
smd_before = {c: abs(smd(df[c].values, A)) for c in CORE}
smd_after  = {c: abs(smd(df[c].values[keep], A[keep], w_iptw[keep])) for c in CORE}
order = sorted(CORE, key=lambda c: smd_before[c])
y=np.arange(len(order))
plt.figure(figsize=(7,6))
plt.scatter([smd_before[c] for c in order], y, label="Unadjusted", color="#c05621", zorder=3)
plt.scatter([smd_after[c] for c in order], y, label="IPTW-weighted", color="#2b6cb0", zorder=3)
for i,c in enumerate(order):
    plt.plot([smd_before[c],smd_after[c]],[i,i],color="grey",lw=0.8,zorder=1)
plt.axvline(0.1, ls="--", c="red", lw=1, label="0.10 threshold")
plt.yticks(y,[LABELS[c] for c in order]); plt.xlabel("Absolute standardized mean difference")
plt.title("Covariate balance: before vs after weighting"); plt.legend(); plt.tight_layout()
plt.savefig(os.path.join(FIG,"love_plot.png"), dpi=150); plt.close()

# balance table
bal = pd.DataFrame({"covariate":[LABELS[c] for c in CORE],
                    "smd_unadjusted":[round(smd_before[c],4) for c in CORE],
                    "smd_weighted":[round(smd_after[c],4) for c in CORE]})
bal.to_csv(os.path.join(TAB,"balance_table.csv"), index=False)

# ================= Step 4(v): mortality descriptive CIF =================
raw = pd.read_csv(os.path.join(RAW,"dat.csv"), index_col=0)
raw["inpatient.number"]=raw["inpatient.number"].astype(str)
df["inpatient.number"]=df["inpatient.number"].astype(str)
mt = raw[["inpatient.number","time.of.death..days.from.admission.","death.within.6.months"]].copy()
mt.columns=["inpatient.number","tdeath","death6m"]
m = df[["inpatient.number","A_digoxin"]].merge(mt,on="inpatient.number",how="left")
def cif_death(sub, horizon=180):
    # 1 - KM for death; censor at horizon
    t = sub["tdeath"].values.astype(float); ev=(sub["death6m"]==1).values
    t = np.where(np.isnan(t), horizon, np.minimum(t,horizon))
    ev = ev & (t<=horizon)
    order=np.argsort(t); t=t[order]; ev=ev[order]
    n=len(t); surv=1.0; xs=[0]; ys=[0.0]; atrisk=n
    for i in range(n):
        if ev[i]:
            surv*= (1 - 1/atrisk)
            xs.append(t[i]); ys.append(1-surv)
        atrisk-=1
    xs.append(horizon); ys.append(ys[-1])
    return np.array(xs), np.array(ys)
plt.figure(figsize=(7,4.5))
for arm,lab,col in [(1,"Digoxin","#2b6cb0"),(0,"No digoxin","#c05621")]:
    xs,ys=cif_death(m[m["A_digoxin"]==arm]); plt.step(xs,ys,where="post",label=lab,color=col)
plt.xlabel("Days from admission"); plt.ylabel("Cumulative incidence of death")
plt.title("Descriptive 6-month mortality CIF by arm (NOT causal)")
plt.legend(); plt.tight_layout(); plt.savefig(os.path.join(FIG,"mortality_cif.png"),dpi=150); plt.close()
results["mortality_descriptive"] = {
    "death6m_digoxin": float(m.loc[m.A_digoxin==1,"death6m"].mean()),
    "death6m_control": float(m.loc[m.A_digoxin==0,"death6m"].mean()),
    "n_deaths_total": int(m["death6m"].sum())
}

# ================= Forest plot (Step 4 iv) =================
def row(label, rr, lcl, ucl): return (label, rr, lcl, ucl)
forest = [
    row("Naive (unadjusted)", results["naive_digoxin"]["RR"],
        results["naive_digoxin"]["RR_lcl"], results["naive_digoxin"]["RR_ucl"]),
    row("AIPW adjusted (primary)", prim["RR"], prim["RR_lcl"], prim["RR_ucl"]),
    row("AIPW cause-specific", csres["RR"], csres["RR_lcl"], csres["RR_ucl"]),
    row("AIPW any-glycoside", gly["RR"], gly["RR_lcl"], gly["RR_ucl"]),
    row("DIG benchmark (HF hosp.)", 0.72, 0.66, 0.79),
]
plt.figure(figsize=(7.5,4.2))
yy=np.arange(len(forest))[::-1]
for (lab,rr,lc,uc),yi in zip(forest,yy):
    plt.plot([lc,uc],[yi,yi],color="#333",lw=1.6)
    plt.plot(rr,yi,"s",color="#2b6cb0",ms=7)
plt.axvline(1.0,ls="--",c="red",lw=1)
plt.yticks(yy,[f[0] for f in forest]); plt.xscale("log")
plt.xlabel("Risk ratio for 6-month readmission (log scale)")
plt.title("Digoxin and 6-month readmission: estimates vs DIG benchmark")
plt.tight_layout(); plt.savefig(os.path.join(FIG,"forest_estimates.png"),dpi=150); plt.close()

# ================= tables =================
def est_row(name, d, scale="readmit"):
    return {"analysis":name,
            "risk_treated":round(d.get("risk_treated",np.nan),4),
            "risk_control":round(d.get("risk_control",np.nan),4),
            "RD":round(d["RD"],4),"RD_lcl":round(d["RD_lcl"],4),"RD_ucl":round(d["RD_ucl"],4),
            "RR":round(d["RR"],4),"RR_lcl":round(d["RR_lcl"],4),"RR_ucl":round(d["RR_ucl"],4),
            "n":d.get("n","")}
prim_tab = pd.DataFrame([
    est_row("Naive (unadjusted)", results["naive_digoxin"]),
    est_row("AIPW adjusted (primary, subdistribution)", prim),
])
prim_tab.to_csv(os.path.join(TAB,"primary_estimates.csv"), index=False)

sens_tab = pd.DataFrame([
    est_row("Primary: AIPW subdistribution (digoxin)", prim),
    est_row("Cause-specific (exclude 57 deaths)", csres),
    est_row("Any cardiac glycoside", gly),
])
sens_tab["evalue_point"]=[round(results["evalue_primary"]["evalue_point"],3),"",""]
sens_tab["evalue_ci"]=[round(results["evalue_primary"].get("evalue_ci",np.nan),3),"",""]
sens_tab.to_csv(os.path.join(TAB,"sensitivity.csv"), index=False)

# positivity report table
pos=prim["positivity"]
pd.DataFrame([pos]).to_csv(os.path.join(TAB,"positivity_report.csv"), index=False)

results["positivity"]=pos
with open(os.path.join(ROOT,"output","results.json"),"w") as f:
    json.dump(results,f,indent=2)

# console summary
print("\n================ SUMMARY ================")
print(f"Primary RD = {prim['RD']:+.4f} ({prim['RD_lcl']:+.4f}, {prim['RD_ucl']:+.4f})")
print(f"Primary RR = {prim['RR']:.3f} ({prim['RR_lcl']:.3f}, {prim['RR_ucl']:.3f})")
print(f"E-value point={results['evalue_primary']['evalue_point']:.2f} ci={results['evalue_primary'].get('evalue_ci')}")
print(f"Cause-specific RR = {csres['RR']:.3f} ({csres['RR_lcl']:.3f}, {csres['RR_ucl']:.3f})")
print(f"Any-glycoside RR = {gly['RR']:.3f} ({gly['RR_lcl']:.3f}, {gly['RR_ucl']:.3f})")
print(f"Naive RR = {results['naive_digoxin']['RR']:.3f}")
print(f"Positivity: PS[{pos['ps_min']:.3f},{pos['ps_max']:.3f}], {pos['pct_outside_report_bounds']:.1f}% outside [0.05,0.95], trimmed {pos['n_trimmed']}")
print("DONE_MARKER")
