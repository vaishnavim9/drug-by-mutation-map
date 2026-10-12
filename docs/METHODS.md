# Methods: where every table came from

All data are public. Queries were run in 2026. Nothing here involved a wet-lab experiment.

## Gene universe

**Open Targets Platform** GraphQL API (`api.platform.opentargets.org`), `disease.associatedTargets`,
filtered to `genetic_association` datatype score >= 0.6.

- The 142-gene set (`data/dee_gene_table.csv`) is the union over three disease terms:
  MONDO_0100062 (genetic developmental and epileptic encephalopathy), MONDO_0100620 (DEE),
  MONDO_0020071 (infantile epilepsy syndrome).
- The 298-gene set (`data/audited_gene_universe.csv`) re-runs the same query over seven terms,
  adding MONDO_0005027 (epilepsy), MONDO_0020072, MONDO_0100619 and HP_0001250 (Seizure). It
  surfaces 156 genes the narrower query missed, including SLC6A1, KCNQ3, DYRK1A, GRIN2A, SLC2A1,
  MECP2, TSC1, PURA, IQSEC2 and the neuronal ceroid lipofuscinosis family.
- **Consequence: every count taken from the 142-gene table is a floor, not a total.** Both tables
  are committed so the difference is inspectable.

Constraint (`loeuf`) comes from the Open Targets target annotation `constraint` field (gnomAD
LOEUF). Note for anyone re-running the query: that type exposes `oeUpper`, `oe` and `upperBin6`,
not `upperBound`.

## Variant spectrum

**ClinVar** via NCBI E-utilities, counting pathogenic and likely-pathogenic submissions per gene
by molecular consequence. Throttle to two concurrent requests or fewer; above that, E-utilities
silently drops requests rather than returning an error, which produces quietly wrong counts.

`mech_class` in `data/dee_gene_table.csv` is assigned from the ratio of truncating to missense
pathogenic variants together with the published mechanism for that gene:

| class | rule |
|---|---|
| LoF / haploinsufficiency | truncating fraction high, no established gain-of-function |
| missense-dominant (GoF / dominant-negative) | missense-dominated and a published gain-of-function or dominant-negative mechanism |
| mixed LoF + missense | both directions documented in the same gene |
| missense-predominant, LoF tolerated | missense-dominated but direction unresolved in the literature |
| unclassified | not enough read to decide |

The last two classes plus `mixed` are why the README says the gene name alone does not determine
direction for 53 of 142 genes. The 6 unclassified genes are excluded from that count on purpose:
"not yet read" is a different claim from "unresolved."

## Clinical trials

**ClinicalTrials.gov API v2**, queried per gene and per disease term. Three traps, all of which
bit during the original work:

1. The field is `Phase`, not `Phases`.
2. Condition search performs loose synonym expansion, so results must be re-filtered on key
   tokens from the gene or disease name.
3. The same trial is returned by several gene queries; deduplicate on NCT identifier.

`modality_class` and `sponsor_class` in `data/clinical_trial_landscape.csv` were assigned by
reading intervention text; `sponsor_class` is the registry's own `leadSponsorClass`.

## Measured-direction census

`data/cross_family_signed_coverage.csv` counts, for 36 genes across six disease families, how many
pathogenic missense variants have ever had a **signed** functional measurement, meaning gain versus
loss rather than merely "damaging."

Sources were found by a **PubMed abstract sweep** (E-utilities esearch/esummary/efetch) plus
**MaveDB** (`api.mavedb.org`). Cohort sizes were extracted from abstracts.

Two corrections that matter, because both were wrong in the first pass:

- An earlier version of this analysis counted 827 signed SCN1A/SCN2A-family variants from the
  funNCion paper (doi 10.1126/scitranslmed.aay6848). **Those labels are phenotype-inferred, not
  measured**; only 87 variants in that work were functionally tested. The corpus must never be
  described as ~800 electrophysiology labels.
- The largest single measured signed cohort found anywhere is **331 KCNQ2 variants** (PMID
  38142737). The abstract sweep alone surfaced a maximum of 69 per gene, so abstract-level
  extraction systematically undercounts; **non-epilepsy family figures in this table are lower
  bounds.**

MaveDB coverage for context: 2,823 score sets covering 654 distinct human proteins, about 3.3% of
the proteome. Three of the 142 epilepsy genes have any multiplexed data (ATP7B, KCNQ2, SPTAN1).

## Benchmark

`data/kcnq2_signed_labels.csv`: 291 KCNQ2 clinical variants from MaveDB score set
`urn:mavedb:00000664`, calibrated against 57 benign controls in the same experiment to set the
range within which a variant is called unchanged. Resulting labels: 214 loss of function, 16 gain
of function, 9 mixed, 52 within control range.

The prediction task in `analysis/02_benchmark.py` is gain versus loss among the 230 confidently
directional variants. Features are relative position, a one-hot of the UniProt topological region
(accession O43526), and changes in hydropathy, formal charge and residue volume. Two splits:
random 5-fold repeated over 20 seeds, and hold out one entire structural region.

The committed `data/signed_prediction_benchmark.csv` is the original run, which held out a
different region and included a slightly larger feature set; the script in this repository is a
from-scratch reimplementation and reports its own numbers to `results/`. Both show the same
collapse.

A paralog test was also run and was **negative**: adding KCNQ1/3/4/5 alignment features changed
AUROC by between -0.047 and +0.023, i.e. noise. The mechanistic reason is that within a single
gene, paralog conservation is a position property already carried by relative position. This
closes the within-gene version of the objection but does not test cross-gene transfer, which
needs measured variants in a second gene.

## Simulation

`analysis/03_simulation.py` is self-contained and documents its own assumptions in the module
docstring. Two settings matter and are both arguable: the response matrix is assumed low-rank
(that is the hypothesis the first real dataset tests), and assay noise is set at 0.25 response
standard deviations.

Sequencing cost, computed separately: 20,000 variants × 26 conditions × 500 reads per variant per
condition is 260 million reads, roughly an eighth of one NovaSeq lane. Variant count is therefore
nearly free; the cost driver is the number of drug conditions.

## Protein class census

**UniProt** REST, reviewed human proteins carrying keyword KW-0225 (Disease variant), classified
by keyword into enzyme, transporter/pump, transcription/chromatin, receptor, ion channel, or
scaffold/adaptor/structural. 3,994 proteins total; 2,503 in a class where a two-directional
selection is constructible; 1,491 scaffolds where it is not. `data/assayability_by_protein_class.csv`
holds the per-gene assignment for the 298-gene epilepsy universe (209 of 298 selectable).
