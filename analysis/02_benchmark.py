"""
02. The falsification test: can direction be predicted from sequence alone?

This is the result that decided the whole project, and it is negative. Models trained on the
only large set of directionally measured ion-channel variants (291 KCNQ2 variants, calibrated
against 57 benign controls) do well when the test variants are drawn at random from the same
structural regions as the training variants, and collapse to chance when a whole structural
region is held out, which is the situation that actually matters, because a new patient's
mutation is frequently somewhere nobody has measured.

Reads  data/kcnq2_signed_labels.csv
Writes figures/02_benchmark.png, results/02_benchmark_rerun.csv

Reference values from the original run are committed in data/signed_prediction_benchmark.csv.
This script is a from-scratch reimplementation: the held-out-region split is defined below
(the region carrying the most gain-of-function variants is held out), so values can differ in
the second decimal. The conclusion does not depend on which region is held out.

Topology from UniProt O43526 (KCNQ2_HUMAN), 872 aa, retrieved 2026.
"""
import pathlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline

ROOT = pathlib.Path(__file__).resolve().parents[1]
RNG = 0
LENGTH = 872

# UniProt feature table for O43526 -------------------------------------------------------------
REGIONS = [
    ("N-terminal cytoplasmic", 1, 90),
    ("S1", 91, 113), ("extracellular S1-S2", 114, 123),
    ("S2", 124, 145), ("cytoplasmic S2-S3", 146, 163),
    ("S3", 164, 183), ("extracellular S3-S4", 184, 196),
    ("S4 voltage sensor", 197, 215), ("cytoplasmic S4-S5", 216, 227),
    ("S5", 228, 253), ("extracellular S5-pore", 254, 263),
    ("pore helix H5", 264, 276), ("extracellular pore-S6", 277, 287),
    ("S6", 288, 314), ("C-terminal cytoplasmic", 315, 872),
]
REGION_NAMES = [r[0] for r in REGIONS]

# Kyte-Doolittle hydropathy, formal charge at pH 7, and residue volume (A^3).
KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2, I=4.5,
          L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2)
CHG = dict(D=-1, E=-1, K=1, R=1, H=0.1)
VOL = dict(A=88.6, R=173.4, N=114.1, D=111.1, C=108.5, Q=143.8, E=138.4, G=60.1, H=153.2,
           I=166.7, L=166.7, K=168.6, M=162.9, F=189.9, P=112.7, S=89.0, T=116.1, W=227.8,
           Y=193.6, V=140.0)


def region_of(pos):
    for name, lo, hi in REGIONS:
        if lo <= pos <= hi:
            return name
    return "unassigned"


def featurise(df):
    X = pd.DataFrame(index=df.index)
    X["rel_pos"] = df["pos"] / LENGTH
    X["hyd_ref"] = df["ref"].map(KD)
    X["hyd_alt"] = df["alt"].map(KD)
    X["d_hyd"] = X["hyd_alt"] - X["hyd_ref"]
    X["d_charge"] = df["alt"].map(lambda a: CHG.get(a, 0)) - df["ref"].map(lambda a: CHG.get(a, 0))
    X["d_volume"] = (df["alt"].map(VOL) - df["ref"].map(VOL)) / 100.0
    X["to_proline"] = (df["alt"] == "P").astype(int)
    X["from_glycine"] = (df["ref"] == "G").astype(int)
    for name in REGION_NAMES:
        X["reg_" + name.replace(" ", "_")] = (df["region"] == name).astype(int)
    return X.fillna(0.0)


def models():
    return {
        "logistic regression": make_pipeline(StandardScaler(),
                                             LogisticRegression(max_iter=4000,
                                                                class_weight="balanced")),
        "random forest": RandomForestClassifier(n_estimators=600, min_samples_leaf=2,
                                                class_weight="balanced", random_state=RNG),
        "gradient boosting": GradientBoostingClassifier(random_state=RNG),
    }


lab = pd.read_csv(ROOT / "data" / "kcnq2_signed_labels.csv")
lab["region"] = lab["pos"].map(region_of)

# The prediction task: among variants with a confident direction, is this one gain of function?
task = lab[lab["signed"].isin(["LoF", "GoF"])].reset_index(drop=True)
y = (task["signed"] == "GoF").astype(int).values
X = featurise(task)
print(f"{len(task)} variants with a confident direction: "
      f"{int(y.sum())} gain of function, {int((1-y).sum())} loss of function")
print(f"class balance (what a coin flip scores on precision): {y.mean():.4f}\n")

rows = []

# --- split A: random 5-fold, repeated ---------------------------------------------------------
for name, mdl in models().items():
    aucs, aps = [], []
    for seed in range(20):
        cv = StratifiedKFold(5, shuffle=True, random_state=seed)
        pred = np.zeros(len(y), float)
        for tr, te in cv.split(X, y):
            m = models()[name]
            m.fit(X.iloc[tr], y[tr])
            pred[te] = m.predict_proba(X.iloc[te])[:, 1]
        aucs.append(roc_auc_score(y, pred))
        aps.append(average_precision_score(y, pred))
    rows.append(dict(model=name, split="random 5-fold (20 seeds)",
                     auroc=np.mean(aucs), auroc_sd=np.std(aucs), auprc=np.mean(aps),
                     baseline_auprc=y.mean()))
    print(f"{name:<20} random 5-fold      AUROC {np.mean(aucs):.3f} +/- {np.std(aucs):.3f}")

# --- split B: hold out an entire structural region --------------------------------------------
gof_by_region = task.loc[y == 1, "region"].value_counts()
held = gof_by_region.index[0]
te = (task["region"] == held).values
tr = ~te
print(f"\nholding out the region with the most gain-of-function variants: {held} "
      f"(n={te.sum()}, {int(y[te].sum())} of them gain of function)")
for name in models():
    m = models()[name]
    m.fit(X[tr], y[tr])
    p = m.predict_proba(X[te])[:, 1]
    auc = roc_auc_score(y[te], p)
    rows.append(dict(model=name, split=f"held-out region ({held})", auroc=auc, auroc_sd=np.nan,
                     auprc=average_precision_score(y[te], p), baseline_auprc=y[te].mean()))
    print(f"{name:<20} held-out region    AUROC {auc:.3f}")

res = pd.DataFrame(rows)
res.to_csv(ROOT / "results" / "02_benchmark_rerun.csv", index=False)

# --- figure -----------------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.4, 4.4))
names = list(models())
w = 0.36
rand = [res[(res.model == n) & (res.split.str.startswith("random"))].auroc.iloc[0] for n in names]
sds = [res[(res.model == n) & (res.split.str.startswith("random"))].auroc_sd.iloc[0] for n in names]
hold = [res[(res.model == n) & (res.split.str.startswith("held"))].auroc.iloc[0] for n in names]
xs = np.arange(len(names))
ax.bar(xs - w/2, rand, w, yerr=sds, capsize=3, label="random split", color="#4C72B0")
ax.bar(xs + w/2, hold, w, label=f"held-out region ({held})", color="#C44E52")
ax.axhline(0.5, color="k", ls="--", lw=1)
ax.text(len(names) - 0.55, 0.512, "chance", fontsize=8.5)
ax.set_xticks(xs)
ax.set_xticklabels(names, fontsize=9)
ax.set_ylabel("AUROC (gain of function vs loss of function)")
ax.set_ylim(0, 1.0)
ax.set_title("Sequence features predict direction only inside\nstructural regions they were "
             "trained on", fontsize=11)
ax.legend(fontsize=9, loc="upper right")
fig.tight_layout()
fig.savefig(ROOT / "figures" / "02_benchmark.png", dpi=200)
print("\nwrote figures/02_benchmark.png and results/02_benchmark_rerun.csv")
