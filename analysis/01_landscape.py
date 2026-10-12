"""
01. The landscape: how many genes, how much measured, how much treated.

Reads  data/dee_gene_table.csv, data/clinical_trial_landscape.csv,
       data/cross_family_signed_coverage.csv
Writes figures/01_landscape.png, results/01_landscape_summary.txt

No network access required. Source queries that built these tables are in docs/METHODS.md.
"""
import pathlib
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parents[1]
genes = pd.read_csv(ROOT / "data" / "dee_gene_table.csv")
trials = pd.read_csv(ROOT / "data" / "clinical_trial_landscape.csv")
families = pd.read_csv(ROOT / "data" / "cross_family_signed_coverage.csv")

# --- 1. mechanism classes -------------------------------------------------------------------
mech = genes["mech_class"].value_counts()

# A gene is "modality-ambiguous" when the direction of a new missense variant cannot be
# predicted from the gene identity alone: either both directions are documented (mixed) or the
# literature has not resolved it (missense-predominant, LoF-tolerated).
ambiguous = int(mech.get("mixed LoF + missense", 0)
                + mech.get("missense-predominant, LoF tolerated", 0))
# (the 6 unclassified genes are excluded from this count: unclassified means the literature
#  was not read far enough to decide, which is a different statement from unresolved.)

# --- 2. clinical coverage -------------------------------------------------------------------
interventional = trials[trials["study_type"].str.contains("INTERVENTIONAL", case=False, na=False)]
industry = interventional[interventional["sponsor_class"].str.contains("INDUSTRY", case=False, na=False)]
genes_with_industry = set(industry["gene"]) & set(genes["symbol"])
gene_directed = industry[industry["modality_class"].str.startswith("gene-directed")]
genes_gene_directed = sorted(set(gene_directed["gene"]) & set(genes["symbol"]))

# --- 3. measurement coverage by gene family -------------------------------------------------
fam = families.sort_values("plp_missense", ascending=False)

summary = [
    f"genes in the curated epilepsy set            : {len(genes)}",
    f"ClinVar pathogenic / likely-pathogenic       : {int(genes['clinvar_plp'].sum()):,}",
    f"  of which missense                          : {int(genes['clinvar_missense'].sum()):,}"
    f" ({genes['clinvar_missense'].sum()/genes['clinvar_plp'].sum():.0%})",
    "",
    "mechanism class (how the gene causes disease):",
    *[f"  {k:<42}: {v}" for k, v in mech.items()],
    f"  -> direction NOT inferable from gene name  : {ambiguous} of {len(genes)}"
    f"  (mixed, or missense-predominant and unresolved)",
    "",
    f"genes with any industry interventional trial : {len(genes_with_industry)} of {len(genes)}",
    f"genes with an industry gene-directed trial   : {len(genes_gene_directed)} "
    f"({', '.join(genes_gene_directed)})",
    "",
    "direction measured, by gene family:",
    *[f"  {r.family:<28}: {r.signed_measured:>5} of {r.plp_missense:>6} "
      f"({r.pct_with_signed_call:.2f}%)" for r in fam.itertuples()],
    f"  {'TOTAL':<28}: {fam.signed_measured.sum():>5} of {fam.plp_missense.sum():>6} "
    f"({100*fam.signed_measured.sum()/fam.plp_missense.sum():.2f}%)",
]
(ROOT / "results" / "01_landscape_summary.txt").write_text("\n".join(summary) + "\n")
print("\n".join(summary))

# --- figure ---------------------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

ax = axes[0]
lab = [m.replace("_", " ") for m in mech.index]
ax.barh(range(len(mech)), mech.values, color="#4C72B0")
ax.set_yticks(range(len(mech)))
ax.set_yticklabels(lab, fontsize=9)
ax.invert_yaxis()
ax.set_xlabel("genes")
ax.set_title(f"How each of {len(genes)} epilepsy genes fails\n"
             f"({ambiguous} do not determine a new variant's direction)", fontsize=10)
for i, v in enumerate(mech.values):
    ax.text(v + 0.6, i, str(v), va="center", fontsize=9)

ax = axes[1]
have, lack = len(genes_with_industry), len(genes) - len(genes_with_industry)
ax.bar(["any industry\ntrial", "none"], [have, lack], color=["#55A868", "#C44E52"])
ax.set_ylabel("genes")
ax.set_title("Industry clinical activity\nacross the same genes", fontsize=10)
for i, v in enumerate([have, lack]):
    ax.text(i, v + 2, str(v), ha="center", fontsize=10)
ax.annotate(f"industry gene-directed programmes\nexist for only "
            f"{len(genes_gene_directed)}: {', '.join(genes_gene_directed)}",
            xy=(0.5, 0.72), xycoords="axes fraction", ha="center", fontsize=8.5,
            bbox=dict(boxstyle="round,pad=0.4", fc="#FFF3CD", ec="#D9B300", lw=0.8))

ax = axes[2]
ax.barh(range(len(fam)), fam["pct_with_signed_call"], color="#8172B2")
ax.set_yticks(range(len(fam)))
ax.set_yticklabels([f.replace(" / ", "/\n") for f in fam["family"]], fontsize=8.5)
ax.invert_yaxis()
ax.set_xlabel("% of pathogenic missense with\ndirection ever measured")
ax.set_title("Direction has essentially never\nbeen measured anywhere", fontsize=10)
for i, r in enumerate(fam.itertuples()):
    ax.text(r.pct_with_signed_call + 0.3, i, f"{r.signed_measured}/{r.plp_missense}",
            va="center", fontsize=8.5)
ax.set_xlim(0, 24)

fig.tight_layout()
fig.savefig(ROOT / "figures" / "01_landscape.png", dpi=200)
print("\nwrote figures/01_landscape.png")
