import os, sys, pickle, time, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from estimators import crossfit_aipw
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
df=pd.read_csv(os.path.join(ROOT,"data","derived","analysis_dataset.csv"))
Wcols=["age","male","visit_times","NYHA","Killip","BNP","creatinine","eGFR","CCI",
       "heart_rate","SBP","DBP","potassium","sodium","BNP_miss","creatinine_miss",
       "eGFR_miss","CCI_miss","heart_rate_miss","SBP_miss","DBP_miss","potassium_miss","sodium_miss"]
which=sys.argv[1]; t0=time.time()
if which=="primary":
    sub=df; A=sub["A_digoxin"].values; Y=sub["Y"].values
elif which=="causespecific":
    sub=df[df["died_no_readmit"]==0].reset_index(drop=True); A=sub["A_digoxin"].values; Y=sub["Y"].values
elif which=="glycoside":
    sub=df; A=sub["A_glycoside"].values; Y=sub["Y"].values
res=crossfit_aipw(sub[Wcols].values, A, Y, n_splits=5, seed=42)
res["_idx"]=sub.index.values
pickle.dump(res, open(f"/tmp/aipw_{which}.pkl","wb"))
print(which,"done in %.1fs"%(time.time()-t0),
      "RR=%.3f (%.3f,%.3f)"%(res["RR"],res["RR_lcl"],res["RR_ucl"]),
      "RD=%.4f (%.4f,%.4f)"%(res["RD"],res["RD_lcl"],res["RD_ucl"]))
