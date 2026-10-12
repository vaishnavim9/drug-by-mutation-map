# Limitations

In priority order. Anyone citing a number from this repository should read this page first.

## 1. No experiment of my own has been run

There is no measurement in this repository that I generated. Every number is computation on
public data or simulation. The central claim, that a pooled, two-sided selection can read the
direction of channel variants at scale, is **unproven**. The published precedent covers a
36-base region of the cardiac sodium channel (12 codons, PMID 31928070, patch-clamp agreement 8
of 9) and was never scaled in the six years since. Whether that was a technical wall or a change
of priorities is unknown to me and is the first question I intend to answer.

## 2. The gene universe is a floor, not a total

The 142-gene set comes from three disease terms. Anchoring on seven terms including the seizure
phenotype returns 298 genes. Every count derived from the 142, whether variant totals, mechanism
class counts or trial coverage, is therefore an underestimate of the true field.

## 3. The benchmark is one gene

Finding 2 rests on KCNQ2, because it is the only gene with a large measured signed set. A single
gene cannot establish that the collapse out of structural context is general, only that it
happens in the best-measured case. The 16 gain-of-function variants also make the positive class
very small, so the held-out-region numbers are noisy by construction.

## 4. Published methods do better than mine where data exist

Direction prediction is not impossible. Activating-versus-deactivating prediction in kinases
reaches AUC 0.94 (PMID 41152984) on two decades of oncology labels, UK Biobank machine learning
has assigned direction in ion channel genes (PMID 39781574), and dual-conformation structure
modelling has inferred direction for JAK1 (PMID 41136770). The honest claim is that signed
prediction works wherever labels exist and sits near chance in ion channels where only 939
variants have been measured.

## 5. The simulation assumes what it would need to discover

Low rank is the central assumption: if real drug responses across variants are not well
approximated by a small number of latent mechanism classes, both the completion result and the
cold-start result weaken. The first real dataset tests this directly and could refute it.
Feature-signal levels in the cold-start panel are also stipulated rather than measured, which is
why three levels are shown instead of one.

## 6. Measurement counts outside epilepsy are lower bounds

Non-epilepsy family figures in `data/cross_family_signed_coverage.csv` come from an abstract-level
PubMed sweep, which systematically undercounts: the sweep found at most 69 variants for a gene
whose true largest cohort is 331. Treat "7 of 4,058 in cardiac arrhythmia genes" as the floor of
what has been published, not a complete census.

## 7. Heterologous cells are not neurons

Any screen of this kind would run in a cell line, and a variant that looks like loss of function
in a kidney-derived line may behave differently in an interneuron. Cystic fibrosis hit exactly
this problem and anchored to patient-derived material. Trafficking-defective variants are a
specific worry: a current-based selection may score them as "no effect" rather than "complete
loss," which is worse than a missing value because it is confidently wrong.

## 8. Several anti-seizure drugs are not inert to cell growth

Valproate affects proliferation and chromatin; a pooled survival readout would confound drug
toxicity with channel rescue unless non-channel control variants ride in the same pool.

## 9. Comorbidity is modelled at the gene level, not the patient level

The statement that 29 of 142 epilepsy genes are also neurodevelopmental genes is a gene-set
overlap. It is not a patient-level comorbidity estimate, and no individual-level data were used.

## 10. Prevalence figures are partly modelled

The roughly 80,000 prevalent paediatric US and EU cases quoted in the pitch rest on per-gene
estimates of which only about half are observation-grounded: 15 genes from population studies, 4
from referral cohorts, 71 from a de novo mutation model, and 52 with no usable evidence at all.
The top five genes account for about a third of the total.
