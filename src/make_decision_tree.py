import os, textwrap, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG=os.path.join(ROOT,"output","figures")

CHOSEN=dict(fc="#dbeafe",ec="#2b6cb0",tc="#13294b")
REJ=dict(fc="#eeeeee",ec="#bbbbbb",tc="#888888")
DEM=dict(fc="#f3efe0",ec="#cebf80",tc="#8a7b3a")

fig,ax=plt.subplots(figsize=(13.5,17.5)); ax.axis("off")
ax.set_xlim(0,13.5); ax.set_ylim(0,18.3)

def box(x,y,txt,style,w=2.7,h=0.95,fs=9.5,bold=False):
    p=FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle="round,pad=0.06,rounding_size=0.10",
                     fc=style["fc"],ec=style["ec"],lw=1.6,
                     ls=("--" if style is not CHOSEN else "-"),zorder=3)
    ax.add_patch(p)
    ax.text(x,y,txt,ha="center",va="center",fontsize=fs,color=style["tc"],
            zorder=4,weight=("bold" if bold else "normal"))

def arrow(p,q,chosen=True):
    col="#2b6cb0" if chosen else "#bbbbbb"
    ls="-" if chosen else (0,(4,3))
    a=FancyArrowPatch(p,q,arrowstyle="-|>",mutation_scale=15,lw=1.7,color=col,
                      ls=ls,shrinkA=2,shrinkB=2,zorder=2); ax.add_patch(a)

def lbl(x,y,txt,chosen=True,ha="left"):
    col="#1b4965" if chosen else "#999999"
    ax.text(x,y,txt,ha=ha,va="center",fontsize=7.6,color=col,style="italic",zorder=5,
            bbox=dict(boxstyle="round,pad=0.18",fc="white",ec="none",alpha=0.85))

xc=5.0; xL=1.7; xR=9.4
# chosen spine y-positions
yA,yC,yD,yDG,yCIF,yOUT,yPOOL,yEST,yDIAG,ySENS = 17.4,15.5,13.6,11.7,9.8,8.0,6.2,4.6,3.2,1.6

# title
ax.text(6.75,18.05,"Project decision logic: from dataset to estimand",
        ha="center",fontsize=15,weight="bold",color="#13294b")

# nodes (chosen)
box(xc,yA,"Dataset:\nZigong HF cohort (n=2,008)",CHOSEN,w=3.2,bold=True)
box(xc,yC,"Causal inference",CHOSEN)
box(xc,yD,"Drug exposure",CHOSEN)
box(xc,yDG,"Oral digoxin\n(primary exposure)",CHOSEN,bold=True)
box(xc,yCIF,"Fixed-horizon\ncompeting-risks CIF",CHOSEN)
box(xc,yOUT,"6-month readmission\n= causal endpoint",CHOSEN,bold=True)
box(xc,yPOOL,"Pooled estimate\n(EF missingness = limitation)",CHOSEN,w=3.2)
box(xc,yEST,"DR estimation:\nTMLE / AIPW + Super Learner",CHOSEN,w=3.4)
box(xc,yDIAG,"Diagnostics:\noverlap · balance · positivity",CHOSEN,w=3.4)
box(xc,ySENS,"Sensitivity: E-value · cause-specific vs.\nsubdistribution · any-glycoside ·\nnaive vs. adjusted vs. DIG",CHOSEN,w=4.6,h=1.25,fs=9)

# rejected / demoted side nodes
box(xL,yC,"Predictive ML",REJ,w=2.4)
box(xL,yD,"Admission route /\ndepartment",REJ,w=2.4)
box(xR,yD+0.55,"Spironolactone",REJ,w=2.4)
box(xR,yDG-0.15,"Any glycoside\n(demoted → sensitivity)",DEM,w=2.7)
box(xL,yCIF,"Fine–Gray /\ncause-specific Cox",REJ,w=2.4)
box(xR,yOUT,"6-month mortality\n(demoted → descriptive)",DEM,w=2.7)
box(xL,yPOOL,"Stratify by\nHFrEF / HFpEF",REJ,w=2.4)

# chosen arrows + labels
def seg(y1,y2): return (xc,y1-0.48),(xc,y2+0.48)
arrow(*seg(yA,yC));  lbl(xc+0.25,(yA+yC)/2,"no treatment-effect work exists\nhere; matches RWE goal")
arrow(*seg(yC,yD));  lbl(xc+0.25,(yC+yD)/2,"well-defined, manipulable;\nequipoise + RCT benchmark")
arrow(*seg(yD,yDG)); lbl(xc+0.25,(yD+yDG)/2,"998 vs 1,010 balanced; DIG\nbenchmark; oral = what DIG randomized")
arrow(*seg(yDG,yCIF));lbl(xc+0.25,(yDG+yCIF)/2,"validated 0/1 indicators, complete;\ncompeting event handled by definition")
arrow(*seg(yCIF,yOUT));lbl(xc+0.25,(yCIF+yOUT)/2,"773 events — well-powered")
arrow(*seg(yOUT,yPOOL));lbl(xc+0.25,(yOUT+yPOOL)/2,"LVEF 68.4% missing; type.of.heart.failure\nis anatomical — can't rescue it")
arrow(*seg(yPOOL,yEST));lbl(xc+0.25,(yPOOL+yEST)/2,"doubly-robust + flexible nuisances")
arrow(*seg(yEST,yDIAG));lbl(xc+0.25,(yEST+yDIAG)/2,"check identification holds")
arrow(*seg(yDIAG,ySENS));lbl(xc+0.25,(yDIAG+ySENS)/2,"probe robustness")

# rejected/demoted diagonal arrows + labels
def diag(xfrom,yfrom,xto,yto,chosen=False):
    arrow((xfrom,yfrom),(xto,yto),chosen=chosen)
diag(xc-1.6,yA-0.30,xL+1.2,yC+0.30); lbl(2.55,(yA+yC)/2+0.15,"literature saturated:\nreadmission XGBoost, HP\ncomparisons, prioritization\ntool already published",chosen=False,ha="center")
diag(xc-1.6,yC-0.10,xL+1.2,yD+0.20); lbl(2.7,yD+0.62,"not manipulable;\ncollider",chosen=False,ha="center")
diag(xc+1.6,yD+0.20,xR-1.2,yD+0.45); lbl(7.35,yD+0.95,"1,833/2,008 treated, only\n175 unexposed — thin\ncomparator; muddier benchmark",chosen=False,ha="center")
diag(xc+1.6,yD-0.25,xR-1.35,yDG-0.05); lbl(7.4,yDG+0.62,"pools IV deslanoside;\nbreaks DIG benchmark",chosen=False,ha="center")
diag(xc-1.6,yDG-0.30,xL+1.2,yCIF+0.30); lbl(2.6,(yDG+yCIF)/2,"16/57 deaths untimed; no clock\nfor 1,050/1,178 event-free;\n78 indicator–time discordant",chosen=False,ha="center")
diag(xc+1.7,yCIF-0.20,xR-1.35,yOUT+0.25); lbl(7.5,(yCIF+yOUT)/2+0.1,"57 events — underpowered\nfor doubly-robust estimation",chosen=False,ha="center")
diag(xc-1.7,yOUT-0.25,xL+1.2,yPOOL+0.25); lbl(2.6,(yOUT+yPOOL)/2,"cannot define HFrEF/HFpEF\nreliably",chosen=False,ha="center")

# legend
from matplotlib.patches import Patch
leg=[Patch(fc=CHOSEN["fc"],ec=CHOSEN["ec"],label="Chosen"),
     Patch(fc=REJ["fc"],ec=REJ["ec"],label="Rejected"),
     Patch(fc=DEM["fc"],ec=DEM["ec"],label="Demoted to sensitivity")]
ax.legend(handles=leg,loc="lower left",fontsize=9,frameon=True,bbox_to_anchor=(0.01,0.005))
ax.text(6.75,0.35,"Node = a decision · branch label = the data fact or principle that forced the choice",
        ha="center",fontsize=8.5,color="#555",style="italic")

plt.tight_layout()
plt.savefig(os.path.join(FIG,"decision_tree.png"),dpi=150,bbox_inches="tight")
print("decision_tree.png written")
