import os, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Ellipse
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG=os.path.join(ROOT,"output","figures")

# ---------- DAG ----------
fig,ax=plt.subplots(figsize=(8.5,5.2)); ax.axis("off")
nodes={"A":(0.18,0.30,"Digoxin\n(in-hospital)"),
       "Y":(0.82,0.30,"6-month\nreadmission"),
       "W":(0.50,0.72,"Measured baseline W\n(age, NYHA, Killip, BNP,\nrenal, CCI, vitals, K/Na,\nprior admits)"),
       "U":(0.50,0.06,"Unmeasured U\n(atrial fibrillation,\nLVEF, med timing)")}
for k,(x,y,lab) in nodes.items():
    col="#fde2cc" if k=="U" else ("#dbeafe" if k in("A","Y") else "#e6f4ea")
    e=Ellipse((x,y),0.26,0.17,facecolor=col,edgecolor="#333",lw=1.3,zorder=2)
    ax.add_patch(e); ax.text(x,y,lab,ha="center",va="center",fontsize=8.5,zorder=3)
def arrow(p,q,style="-",col="#333"):
    a=FancyArrowPatch(p,q,arrowstyle="-|>",mutation_scale=14,lw=1.5,color=col,
                      shrinkA=20,shrinkB=20,linestyle=style,zorder=1); ax.add_patch(a)
arrow((0.31,0.30),(0.69,0.30))                     # A -> Y (effect of interest)
arrow((0.42,0.66),(0.24,0.37))                     # W -> A
arrow((0.58,0.66),(0.76,0.37))                     # W -> Y
arrow((0.42,0.12),(0.24,0.24),col="#c05621")       # U -> A (confounding, unmeasured)
arrow((0.58,0.12),(0.76,0.24),col="#c05621")       # U -> Y
ax.text(0.50,0.345,"causal effect (target)",ha="center",fontsize=8,style="italic",color="#2b6cb0")
ax.text(0.50,0.93,"Adjusted via doubly-robust AIPW (Super-Learner)",ha="center",fontsize=9,weight="bold")
ax.text(0.50,-0.02,"Open backdoor through U (esp. AF) → residual confounding, quantified by E-value",
        ha="center",fontsize=8,color="#c05621")
ax.set_xlim(0,1); ax.set_ylim(-0.05,1.0)
plt.tight_layout(); plt.savefig(os.path.join(FIG,"dag.png"),dpi=150,bbox_inches="tight"); plt.close()

# ---------- conceptual target-trial figure ----------
fig,ax=plt.subplots(figsize=(9,4.2)); ax.axis("off")
steps=[("Eligibility","Hospitalized HF\npatients (n=2008)"),
       ("Treatment\nstrategies","Digoxin vs.\nno digoxin"),
       ("Assignment","Emulated via\nbaseline W"),
       ("Outcome","6-mo all-cause\nreadmission"),
       ("Estimand","RD & RR\n(AIPW, cross-fit)"),
       ("Sensitivity","E-value · CIF\nglycoside · DIG")]
x=0.04
for i,(t,b) in enumerate(steps):
    ax.add_patch(plt.Rectangle((x,0.38),0.13,0.30,facecolor="#dbeafe",edgecolor="#333",lw=1.2))
    ax.text(x+0.065,0.60,t,ha="center",va="center",fontsize=8.5,weight="bold")
    ax.text(x+0.065,0.47,b,ha="center",va="center",fontsize=7.6)
    if i<len(steps)-1:
        ax.add_patch(FancyArrowPatch((x+0.13,0.53),(x+0.17,0.53),arrowstyle="-|>",
                     mutation_scale=12,lw=1.4,color="#333"))
    x+=0.16
ax.text(0.5,0.85,"Target-trial emulation: in-hospital digoxin and 6-month readmission",
        ha="center",fontsize=11,weight="bold")
ax.text(0.5,0.16,"Pre-registered protocol → DAG → doubly-robust analysis → pre-specified sensitivity",
        ha="center",fontsize=8.5,style="italic",color="#555")
ax.set_xlim(0,1.02); ax.set_ylim(0,1)
plt.tight_layout(); plt.savefig(os.path.join(FIG,"conceptual_target_trial.png"),dpi=150,bbox_inches="tight"); plt.close()
print("DAG + conceptual figures written")
