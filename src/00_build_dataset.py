"""
00_build_dataset.py
Target-trial emulation: in-hospital digoxin and 6-month readmission.
Builds the admission-time analysis dataset from PhysioNet Zigong HF data.

RAW DATA IS NEVER WRITTEN TO A TRACKED LOCATION. Derived data -> data/derived/ (gitignored).
"""
import os, json
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW  = os.path.join(ROOT, "data", "raw", "physionet")
DERIV= os.path.join(ROOT, "data", "derived")
os.makedirs(DERIV, exist_ok=True)

dat = pd.read_csv(os.path.join(RAW, "dat.csv"), index_col=0)
md  = pd.read_csv(os.path.join(RAW, "dat_md.csv"), index_col=0)
dat["inpatient.number"] = dat["inpatient.number"].astype(str)
md["inpatient.number"]  = md["inpatient.number"].astype(str)
md["Drug_name"] = md["Drug_name"].str.strip()

# ----- Exposures (per-patient binary from medication long table) -----
dig = set(md.loc[md["Drug_name"].str.lower() == "digoxin tablet", "inpatient.number"])
des = set(md.loc[md["Drug_name"].str.lower() == "deslanoside injection", "inpatient.number"])
dat["A_digoxin"]   = dat["inpatient.number"].isin(dig).astype(int)        # primary: oral digoxin
dat["A_glycoside"] = dat["inpatient.number"].isin(dig | des).astype(int)  # secondary: any cardiac glycoside

# ----- Outcome & competing event -----
dat["Y"] = dat["re.admission.within.6.months"].astype(int)               # 1 = readmitted within 6m
dat["death6m"] = dat["death.within.6.months"].astype(int)
# died-without-readmission set (subdistribution keeps Y=0; cause-specific excludes)
dat["died_no_readmit"] = ((dat["death6m"] == 1) & (dat["Y"] == 0)).astype(int)

# ----- Confounders W (ADMISSION-TIME ONLY; no post-treatment variables) -----
# AF/arrhythmia: NOT in dataset (see Step 0 verdict in docs/protocol.md) -> excluded, carried by E-value.
# LVEF: 68.4% missing -> excluded by pre-specification.

# age: only decade categories exist -> map to decade midpoints (documented limitation)
agemap = {"(21,29]":25,"(29,39]":34,"(39,49]":44,"(49,59]":54,
          "(59,69]":64,"(69,79]":74,"(79,89]":84,"(89,110]":99}
dat["age"] = dat["ageCat"].map(agemap)

# ordinal severity
dat["NYHA"]   = dat["NYHA.cardiac.function.classification"].map({"II":2,"III":3,"IV":4})
dat["Killip"] = dat["Killip.grade"].map({"I":1,"II":2,"III":3,"IV":4})
dat["male"]   = (dat["gender"] == "Male").astype(int)

# implausible zero vitals -> missing (then imputed)
for v in ["pulse","systolic.blood.pressure","diastolic.blood.pressure"]:
    dat.loc[dat[v] == 0, v] = np.nan

W = {
    "age":                         dat["age"],
    "male":                        dat["male"],
    "visit_times":                 dat["visit.times"],          # prior utilization
    "NYHA":                        dat["NYHA"],                 # HF severity
    "Killip":                      dat["Killip"],               # HF severity
    "BNP":                         dat["brain.natriuretic.peptide"],
    "creatinine":                  dat["creatinine.enzymatic.method"],
    "eGFR":                        dat["glomerular.filtration.rate"],
    "CCI":                         dat["CCI.score"],            # Charlson
    "heart_rate":                  dat["pulse"],
    "SBP":                         dat["systolic.blood.pressure"],
    "DBP":                         dat["diastolic.blood.pressure"],
    "potassium":                   dat["potassium"],
    "sodium":                      dat["sodium"],
}
Wdf = pd.DataFrame(W)

CONT = ["age","visit_times","BNP","creatinine","eGFR","CCI",
        "heart_rate","SBP","DBP","potassium","sodium"]
CAT  = ["male","NYHA","Killip"]   # binary/ordinal -> mode impute

# ----- Missingness indicators + median/mode imputation -----
miss_report = {}
for c in Wdf.columns:
    m = Wdf[c].isna()
    miss_report[c] = float(m.mean())
    if m.any():
        Wdf[c + "_miss"] = m.astype(int)
        if c in CONT:
            Wdf[c] = Wdf[c].fillna(Wdf[c].median())
        else:
            Wdf[c] = Wdf[c].fillna(Wdf[c].mode().iloc[0])

analysis = pd.concat(
    [dat[["inpatient.number","A_digoxin","A_glycoside","Y","death6m","died_no_readmit"]].reset_index(drop=True),
     Wdf.reset_index(drop=True)], axis=1)

analysis.to_csv(os.path.join(DERIV, "analysis_dataset.csv"), index=False)  # gitignored

meta = {
    "n_total": int(len(analysis)),
    "n_readmit": int(analysis["Y"].sum()),
    "readmit_rate": float(analysis["Y"].mean()),
    "n_digoxin": int(analysis["A_digoxin"].sum()),
    "n_glycoside": int(analysis["A_glycoside"].sum()),
    "n_died_no_readmit": int(analysis["died_no_readmit"].sum()),
    "confounders_baseline_missingness": miss_report,
    "W_columns": [c for c in Wdf.columns],
    "AF_in_W": False,
    "LVEF_in_W": False,
    "imputation": "median (continuous) / mode (binary-ordinal) + per-variable missingness indicator; MICE noted as planned refinement",
}
with open(os.path.join(DERIV, "build_meta.json"), "w") as f:
    json.dump(meta, f, indent=2)

print(json.dumps(meta, indent=2))
