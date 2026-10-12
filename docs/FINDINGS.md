# Findings

Each number below is produced by a script in `analysis/` from a table in `data/`. The caveat
attached to each is not optional reading; `docs/LIMITATIONS.md` has the full list.

## The field

| number | meaning | caveat |
|---|---|---|
| **142** | genes with genetic association >= 0.6 for developmental and epileptic encephalopathy | a floor; seven-term query returns 298 |
| **26,142** | pathogenic / likely-pathogenic ClinVar variants in those genes | submission-level, not variant-curated |
| **58%** | of them missense, the class where direction is ambiguous | truncating variants are overwhelmingly loss of function |
| **132 of 142** | genes with no industry interventional trial of any kind | registry-visible trials only |
| **2** | genes with an industry programme aimed at the gene itself (SCN1A, STXBP1) | as of the 2026 query |
| **53 of 142** | genes where the gene name alone does not determine a new variant's direction | 18 with both directions documented, 35 unresolved |

## The measurement gap

| family | direction measured | pathogenic missense | % |
|---|---|---|---|
| epilepsy / neuronal channel | 530 | 4,401 | 12.04% |
| cardiac arrhythmia | 7 | 4,058 | 0.17% |
| metabolic / transporter | 2 | 2,475 | 0.08% |
| skeletal muscle | 3 | 2,001 | 0.15% |
| neurotransmitter receptor | 292 | 1,572 | 18.58% |
| inborn errors of immunity | 105 | 824 | 12.74% |
| **total** | **939** | **15,331** | **6.12%** |

Two thirds of the 939 come from two studies (331 KCNQ2 variants, PMID 38142737; 272 GABA-A
receptor variants, PMID 42546502). Non-epilepsy rows are lower bounds from an abstract sweep.

## The benchmark

Gain of function versus loss of function among 230 KCNQ2 variants with a confident measured
direction (16 gain, 214 loss).

| model | random 5-fold | one structural region held out |
|---|---|---|
| logistic regression | 0.65 | 0.64 |
| random forest | 0.78 | 0.54 |
| gradient boosting | 0.80 | 0.44 |

Values from `results/02_benchmark_rerun.csv`; the original run
(`data/signed_prediction_benchmark.csv`) held out a different region and reported 0.72 / 0.77 /
0.80 random against 0.61 / 0.42 / 0.50 held out. The conclusion is the same under both splits:
useful inside the structural context the model was trained on, chance outside it.

A paralog-transfer test was also negative: adding KCNQ1/3/4/5 alignment features moved AUROC by
between -0.047 and +0.023.

## The simulation

From `results/03_matrix_completion.csv` and `results/03_coldstart.csv` (five latent mechanism
classes):

- Recovering unmeasured matrix entries: r = 0.46 at 10% of entries measured, 0.68 at 15%, 0.80 at
  20%, 0.90 at 30%.
- Predicting the full drug-response profile of a variant never assayed, at the middle
  feature-quality setting: r = 0.15 at a panel of 50 variants, 0.31 at 100, 0.46 at 250, 0.54 at
  500, 0.56 at 1,000. At the high setting: 0.51 / 0.75 / 0.85 / 0.88 / 0.89.
- **The returns are strongly diminishing rather than flat.** Growing the panel from 50 to 250
  variants adds 0.31 to r at the middle setting and 0.34 at the high setting; growing it from 250
  to 1,000, which is four times the wet-lab cost, adds 0.11 and 0.04. The ceiling is set by
  representation quality, not panel size: at a poor representation no panel size rescues it.

Sequencing 20,000 variants across 26 conditions at 500 reads each is 260 million reads, about an
eighth of a NovaSeq lane.

## The reach

| scope | number |
|---|---|
| epilepsy genes also in a 492-gene neurodevelopmental set | 29 of 142 (60 of the wider 298) |
| industry autism trials in phase 2/3 using a physiological primary endpoint | 0 of 91 |
| human ion channel genes carrying disease-causing mutations | 148 |
| human disease genes in a two-directionally selectable protein class | 2,503 of 3,994 |
| of those, scaffolding proteins where it is probably not constructible | 1,491 |
