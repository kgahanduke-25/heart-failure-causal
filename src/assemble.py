import os, json, pickle, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from estimators import naive_rr_rd, evalue_rr

ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DERIV=os.path.join(ROOT,"data","derived"); RAW=os.path.join(ROOT,"data","raw","physionet")
FIG=os.path.join(ROOT,"output","figures"); TAB=os.path.join(ROOT,"output","tables")
os.makedirs(FIG,exist_ok=True); os.makedirs(TAB,exist_ok=True)

df=pd.read_csv(os.path.join(DERIV,"analysis_dataset.csv"))
df["inpatient.number"]=df["inpatient.number"].astype(str)
prim=pickle.load(open("/tmp/aipw_primary.pkl","rb"))
cs  =pickle.load(open("/tmp/aipw_causespecific.pkl","rb"))
gly =pickle.load(open("/tmp/aipw_glycoside.pkl","rb"))

CORE=["age","male","visit_times","NYHA","Killip","BNP","creatinine","eGFR","CCI",
      "heart_rate","SBP","DBP","potassium","sodium"]
LAB={"age":"Age","male":"Male","visit_times":"Prior admissions","NYHA":"NYHA class",
     "Killip":"Killip grade","BNP":"BNP","creatinine":"Creatinine","eGFR":"eGFR",
     "CCI":"Charlson index","heart_rate":"Heart rate","SBP":"SBP","DBP":"DBP",
     "potassium":"Potassium","sodium":"Sodium"}

res={}
def clean(d): return {k:v for k,v in d.items() if k not in ("ps","mu1","mu0","keep","_idx")}
res["primary_digoxin_subdist"]=clean(prim)
res["causespecific_digoxin"]=clean(cs)
res["secondary_anyglycoside"]=clean(gly)
res["naive_digoxin"]=naive_rr_rd(df["A_digoxin"].values, df["Y"].values)
res["evalue_primary"]=evalue_rr(prim["RR"],prim["RR_lcl"],prim["RR_ucl"])
res["dig_benchmark"]={"label":"DIG trial (NEJM 1997) — HF hospitalization","RR":0.72,
    "RR_lcl":0.66,"RR_ucl":0.79,"mortality_RR":0.99,
    "note":"RCT, ambulatory, reduced EF, sinus rhythm — different population/outcome than 6-mo all-cause readmission"}

# ---- diagnostics ----
g=prim["ps"]; A=df["A_digoxin"].values
plt.figure(figsize=(7,4.5))
plt.hist(g[A==1],bins=30,alpha=0.6,label="Digoxin (treated)",color="#2b6cb0",density=True)
plt.hist(g[A==0],bins=30,alpha=0.6,label="No digoxin (control)",color="#c05621",density=True)
for b in (0.05,0.95): plt.axvline(b,ls="--",c="grey",lw=1)
plt.xlabel("Estimated propensity score  P(digoxin | W)"); plt.ylabel("Density")
plt.title("Propensity-score overlap by arm (cross-fitted)"); plt.legend()
plt.tight_layout(); plt.savefig(os.path.join(FIG,"ps_overlap.png"),dpi=150); plt.close()

def smd(x,a,w=None):
    x=np.asarray(x,float); a=np.asarray(a,int)
    if w is None:
        m1,m0=x[a==1].mean(),x[a==0].mean(); v1,v0=x[a==1].var(ddof=1),x[a==0].var(ddof=1)
    else:
        w1,w0=w[a==1],w[a==0]; x1,x0=x[a==1],x[a==0]
        m1=np.average(x1,weights=w1);m0=np.average(x0,weights=w0)
        v1=np.average((x1-m1)**2,weights=w1);v0=np.average((x0-m0)**2,weights=w0)
    sp=np.sqrt((v1+v0)/2); return 0.0 if sp==0 else (m1-m0)/sp
keep=prim["keep"]; gk=np.clip(g,0.02,0.98); w=A/gk+(1-A)/(1-gk)
sb={c:abs(smd(df[c].values,A)) for c in CORE}
sa={c:abs(smd(df[c].values[keep],A[keep],w[keep])) for c in CORE}
order=sorted(CORE,key=lambda c:sb[c]); y=np.arange(len(order))
plt.figure(figsize=(7,6))
plt.scatter([sb[c] for c in order],y,label="Unadjusted",color="#c05621",zorder=3)
plt.scatter([sa[c] for c in order],y,label="IPTW-weighted",color="#2b6cb0",zorder=3)
for i,c in enumerate(order): plt.plot([sb[c],sa[c]],[i,i],color="grey",lw=.8,zorder=1)
plt.axvline(0.1,ls="--",c="red",lw=1,label="0.10 threshold")
plt.yticks(y,[LAB[c] for c in order]); plt.xlabel("Absolute standardized mean difference")
plt.title("Covariate balance: before vs after weighting"); plt.legend()
plt.tight_layout(); plt.savefig(os.path.join(FIG,"love_plot.png"),dpi=150); plt.close()
pd.DataFrame({"covariate":[LAB[c] for c in CORE],
             "smd_unadjusted":[round(sb[c],4) for c in CORE],
             "smd_weighted":[round(sa[c],4) for c in CORE]}).to_csv(
             os.path.join(TAB,"balance_table.csv"),index=False)

# ---- mortality descriptive CIF ----
raw=pd.read_csv(os.path.join(RAW,"dat.csv"),index_col=0); raw["inpatient.number"]=raw["inpatient.number"].astype(str)
mt=raw[["inpatient.number","time.of.death..days.from.admission.","death.within.6.months"]].copy()
mt.columns=["inpatient.number","tdeath","death6m"]
m=df[["inpatient.number","A_digoxin"]].merge(mt,on="inpatient.number",how="left")
def cif(sub,h=180):
    t=sub["tdeath"].values.astype(float); ev=(sub["death6m"]==1).values
    t=np.where(np.isnan(t),h,np.minimum(t,h)); ev=ev&(t<=h)
    o=np.argsort(t); t=t[o]; ev=ev[o]; n=len(t); s=1.; xs=[0]; ys=[0.]; ar=n
    for i in range(n):
        if ev[i]: s*=(1-1/ar); xs.append(t[i]); ys.append(1-s)
        ar-=1
    xs.append(h); ys.append(ys[-1]); return np.array(xs),np.array(ys)
plt.figure(figsize=(7,4.5))
for arm,lab,col in [(1,"Digoxin","#2b6cb0"),(0,"No digoxin","#c05621")]:
    xs,ys=cif(m[m.A_digoxin==arm]); plt.step(xs,ys,where="post",label=lab,color=col)
plt.xlabel("Days from admission"); plt.ylabel("Cumulative incidence of death")
plt.title("Descriptive 6-month mortality CIF by arm (NOT a causal estimate)")
plt.legend(); plt.tight_layout(); plt.savefig(os.path.join(FIG,"mortality_cif.png"),dpi=150); plt.close()
res["mortality_descriptive"]={"death6m_digoxin":float(m.loc[m.A_digoxin==1,"death6m"].mean()),
   "death6m_control":float(m.loc[m.A_digoxin==0,"death6m"].mean()),"n_deaths_total":int(m["death6m"].sum())}

# ---- forest ----
forest=[("Naive (unadjusted)",res["naive_digoxin"]["RR"],res["naive_digoxin"]["RR_lcl"],res["naive_digoxin"]["RR_ucl"]),
        ("AIPW adjusted (primary)",prim["RR"],prim["RR_lcl"],prim["RR_ucl"]),
        ("AIPW cause-specific",cs["RR"],cs["RR_lcl"],cs["RR_ucl"]),
        ("AIPW any-glycoside",gly["RR"],gly["RR_lcl"],gly["RR_ucl"]),
        ("DIG benchmark (HF hosp.)",0.72,0.66,0.79)]
plt.figure(figsize=(7.5,4.2)); yy=np.arange(len(forest))[::-1]
for (lab,rr,lc,uc),yi in zip(forest,yy):
    plt.plot([lc,uc],[yi,yi],color="#333",lw=1.6); plt.plot(rr,yi,"s",color="#2b6cb0",ms=7)
plt.axvline(1.0,ls="--",c="red",lw=1); plt.yticks(yy,[f[0] for f in forest]); plt.xscale("log")
import matplotlib.ticker as mt2
plt.gca().xaxis.set_major_formatter(mt2.FuncFormatter(lambda v,_:f"{v:.2f}"))
plt.xticks([0.6,0.8,1.0,1.2,1.4])
plt.xlabel("Risk ratio for 6-month readmission (log scale)")
plt.title("Digoxin and 6-month readmission: estimates vs DIG benchmark")
plt.tight_layout(); plt.savefig(os.path.join(FIG,"forest_estimates.png"),dpi=150); plt.close()

# ---- tables ----
def er(name,d):
    return {"analysis":name,"risk_treated":round(d.get("risk_treated",np.nan),4),
            "risk_control":round(d.get("risk_control",np.nan),4),
            "RD":round(d["RD"],4),"RD_lcl":round(d["RD_lcl"],4),"RD_ucl":round(d["RD_ucl"],4),
            "RR":round(d["RR"],4),"RR_lcl":round(d["RR_lcl"],4),"RR_ucl":round(d["RR_ucl"],4),
            "n":d.get("n","")}
pd.DataFrame([er("Naive (unadjusted)",res["naive_digoxin"]),
              er("AIPW adjusted (primary, subdistribution)",prim)]).to_csv(
              os.path.join(TAB,"primary_estimates.csv"),index=False)
st=pd.DataFrame([er("Primary: AIPW subdistribution (digoxin)",prim),
                 er("Cause-specific (exclude 57 deaths)",cs),
                 er("Any cardiac glycoside",gly)])
st["evalue_point"]=[round(res["evalue_primary"]["evalue_point"],3),"",""]
st["evalue_ci"]=[round(res["evalue_primary"].get("evalue_ci",np.nan),3),"",""]
st.to_csv(os.path.join(TAB,"sensitivity.csv"),index=False)
pos=prim["positivity"]; pd.DataFrame([pos]).to_csv(os.path.join(TAB,"positivity_report.csv"),index=False)
res["positivity"]=pos
json.dump(res,open(os.path.join(ROOT,"output","results.json"),"w"),indent=2)

print("PRIMARY  RR=%.3f (%.3f,%.3f)  RD=%+.4f (%+.4f,%+.4f)"%(prim["RR"],prim["RR_lcl"],prim["RR_ucl"],prim["RD"],prim["RD_lcl"],prim["RD_ucl"]))
print("CAUSE-SP RR=%.3f (%.3f,%.3f)"%(cs["RR"],cs["RR_lcl"],cs["RR_ucl"]))
print("GLYCOS   RR=%.3f (%.3f,%.3f)"%(gly["RR"],gly["RR_lcl"],gly["RR_ucl"]))
print("NAIVE    RR=%.3f (%.3f,%.3f)"%(res["naive_digoxin"]["RR"],res["naive_digoxin"]["RR_lcl"],res["naive_digoxin"]["RR_ucl"]))
print("E-value point=%.2f  ci=%.2f"%(res["evalue_primary"]["evalue_point"],res["evalue_primary"]["evalue_ci"]))
print("Mortality: dig=%.3f ctrl=%.3f (descriptive)"%(res["mortality_descriptive"]["death6m_digoxin"],res["mortality_descriptive"]["death6m_control"]))
print("Positivity PS[%.3f,%.3f] %.1f%% outside[0.05,0.95] trimmed=%d kept=%d"%(pos["ps_min"],pos["ps_max"],pos["pct_outside_report_bounds"],pos["n_trimmed"],pos["n_kept"]))
print("Max |SMD| after weighting:", round(max(sa.values()),3))
print("ASSEMBLE_DONE")
