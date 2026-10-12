"""
03. Sizing the first experiment by simulation.

Two questions a wet-lab budget depends on, answered before spending anything:

  (a) If a pooled drug-by-mutation screen leaves most entries unmeasured, can the rest be
      recovered?  Simulated by measuring a fraction of a low-rank matrix and completing it.
  (b) How many variants does the panel need, in order to predict the whole drug-response
      profile of a variant that was never in the pool at all?  Simulated by fitting the latent
      position of a variant from its sequence-derived features.

The answer to (b) is the one that changed the plan: returns on panel size fall off sharply after
a few hundred variants, so a first screen should be small in variants and wide in drugs.

Writes results/03_matrix_completion.csv, results/03_coldstart.csv, figures/03_simulation.png
Runs in about a minute on a laptop. Self-contained: no data files, no network.
"""
import pathlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parents[1]

N_VARIANTS = 2000     # size of the simulated patient-variant panel
N_DRUGS = 26          # the anti-seizure formulary plus vehicle and dose variations
NOISE_SD = 0.25       # assay noise, in units of the response standard deviation
N_FEATURES = 64       # dimensionality of the sequence-derived representation


def truth(rank, seed, feature_signal=None, n_variants=N_VARIANTS):
    """A drug-by-variant response matrix with `rank` latent mechanism classes.

    If feature_signal is given, each variant's latent coordinates are partly a linear
    function of its sequence features, which is what makes prediction for an unassayed
    variant possible at all.
    """
    rs = np.random.RandomState(seed)
    W = rs.normal(size=(N_DRUGS, rank))
    if feature_signal is None:
        V = rs.normal(size=(n_variants, rank))
        F = None
    else:
        F = rs.normal(size=(n_variants, N_FEATURES))
        B = rs.normal(size=(N_FEATURES, rank)) / np.sqrt(N_FEATURES)
        V = feature_signal * (F @ B) + np.sqrt(max(1 - feature_signal**2, 0)) * \
            rs.normal(size=(n_variants, rank))
    M = V @ W.T
    M /= M.std()
    return M, V, W, F, rs


def soft_impute(Y, mask, rank, iters=160):
    """Low-rank completion by iterated truncated SVD on the filled matrix."""
    X = np.where(mask, Y, 0.0)
    for _ in range(iters):
        U, s, Vt = np.linalg.svd(X, full_matrices=False)
        s[rank:] = 0.0
        X_low = (U * s) @ Vt
        X = np.where(mask, Y, X_low)
    return X_low


def completion_trial(rank, rate, seed):
    M, *_ = truth(rank, seed)
    rs = np.random.RandomState(seed + 10_000)
    mask = rs.rand(*M.shape) < rate
    Y = M + rs.normal(scale=NOISE_SD, size=M.shape)
    hat = soft_impute(Y.T, mask.T, rank).T
    held = ~mask
    return np.corrcoef(M[held], hat[held])[0, 1]


def coldstart_trial(rank, panel, feature_signal, seed):
    """Measure `panel` variants, then predict the full profile of 400 held-out variants."""
    n = panel + 400
    M, V, W, F, _ = truth(rank, seed, feature_signal=feature_signal, n_variants=n)
    rs = np.random.RandomState(seed + 20_000)
    Y = M + rs.normal(scale=NOISE_SD, size=M.shape)
    tr, te = slice(0, panel), slice(panel, n)
    # recover latent coordinates of the measured panel, then regress them on features
    U, s, Vt = np.linalg.svd(Y[tr], full_matrices=False)
    Z = (U[:, :rank] * s[:rank])          # panel latent coordinates
    load = Vt[:rank]                       # drug loadings
    coef, *_ = np.linalg.lstsq(F[tr], Z, rcond=None)
    pred = (F[te] @ coef) @ load
    return np.corrcoef(M[te].ravel(), pred.ravel())[0, 1]


if __name__ == "__main__":
    rows = []
    for rank in (3, 5, 8):
        for rate in (0.05, 0.10, 0.15, 0.20, 0.30, 0.50):
            rs = [completion_trial(rank, rate, s) for s in range(5)]
            rows.append(dict(rank=rank, rate=rate, r=np.mean(rs), sd=np.std(rs)))
            print(f"completion  rank {rank}  {rate:>4.0%} measured   r = {np.mean(rs):.3f}")
    comp = pd.DataFrame(rows)
    comp.to_csv(ROOT / "results" / "03_matrix_completion.csv", index=False)

    rows = []
    for rank in (3, 5, 8):
        for panel in (50, 100, 250, 500, 1000):
            for sig in (0.3, 0.6, 0.9):
                rs = [coldstart_trial(rank, panel, sig, s) for s in range(4)]
                rows.append(dict(rank=rank, panel=panel, feature_signal=sig, r=np.mean(rs)))
        print(f"cold start  rank {rank}  done")
    cold = pd.DataFrame(rows)
    cold.to_csv(ROOT / "results" / "03_coldstart.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(12.4, 4.6))

    ax = axes[0]
    for sig, c in zip((0.3, 0.6, 0.9), ("#C44E52", "#DD8452", "#4C72B0")):
        sub = cold[(cold["rank"] == 5) & (cold.feature_signal == sig)]
        ax.plot(sub.panel, sub.r, "o-", color=c,
                label=f"feature signal {sig:.1f}")
    ax.set_xscale("log")
    ax.set_xlabel("variants measured in the panel")
    ax.set_ylabel("accuracy on a variant never assayed (r)")
    ax.set_title("Diminishing returns after a few hundred variants;\n"
                 "representation quality sets the ceiling", fontsize=10.5)
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=9, loc="lower right")
    ax.grid(alpha=0.3)

    ax = axes[1]
    for rank, c in zip((3, 5, 8), ("#4C72B0", "#55A868", "#8172B2")):
        sub = comp[comp["rank"] == rank]
        ax.errorbar(sub.rate * 100, sub.r, yerr=sub.sd, fmt="o-", color=c,
                    capsize=3, label=f"{rank} mechanism classes")
    ax.set_xlabel("% of the drug-by-variant matrix actually measured")
    ax.set_ylabel("accuracy on unmeasured entries (r)")
    ax.set_title("Most of the matrix can be left unmeasured\nand recovered", fontsize=10.5)
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=9, loc="lower right")
    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(ROOT / "figures" / "03_simulation.png", dpi=200)
    print("\nwrote figures/03_simulation.png and two tables in results/")
