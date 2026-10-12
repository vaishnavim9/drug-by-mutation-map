# A drug-by-mutation map for childhood genetic epilepsy

**Analysis repository for an early-stage venture idea.** Everything here is computation on public
data. No experiment of my own has been run. The three scripts in `analysis/` reproduce every
number quoted below from the committed tables in `data/`, in about a minute.

---

## The problem, in one paragraph

Children born with a mutation in an ion channel gene often develop severe epilepsy in the first
year of life. The drug that controls one child's seizures can make another child worse, and what
decides it is the *direction* of the mutation: whether it made the channel do too much (gain of
function) or too little (loss of function). Sodium-channel blockers such as carbamazepine are the
treatment of choice in gain of function and are contraindicated in loss of function. Across the
gene families where this difference decides treatment, **15,331 mutations are known to cause
disease and only 939 have ever had their direction measured in a laboratory, which is 6.1%**,
and two thirds of those come from just two studies. For every other child, direction is inferred from symptoms,
usually from the age at which seizures started, and the consequence is sequential medication
trials of weeks each while seizures continue.

## Three findings

**1. The gap is the whole field, not one gene.** Of 142 genes with strong genetic evidence for
developmental and epileptic encephalopathy, 132 have no industry interventional trial of any
kind, and exactly two of them, SCN1A and STXBP1, have an industry programme aimed at the gene itself.
For 53 of the 142, the gene name alone does not determine a new variant's direction: either both
directions are documented in that gene, or the literature has not resolved it.
→ `analysis/01_landscape.py`, `figures/01_landscape.png`

**2. The gap is not an algorithm problem, and I tested this on my own idea first.** Using the
largest directionally measured ion-channel variant set that exists (291 KCNQ2 variants against 57
benign controls), models trained on sequence and topology features reach AUROC 0.65 to 0.80 under a
random split and **fall to 0.44 to 0.64 when an entire structural region is held out**, which is
the situation that matters, because a new patient's mutation is often somewhere nobody has
measured. Related published work predicts activating-versus-deactivating direction in kinases at
AUC 0.94, where two decades of oncology produced thousands of labelled variants. The difference
between 0.94 and chance is the training data, not the method.
→ `analysis/02_benchmark.py`, `figures/02_benchmark.png`

**3. The experiment that closes the gap can be small.** Simulating a pooled screen (patient
variants in one culture, split across drug conditions, read out by sequencing) shows that most
of the drug-by-mutation matrix can be left unmeasured and recovered (r of 0.77 to 0.90 from 20 to 30% of
entries), and more importantly that **the panel needed to predict a never-assayed variant's whole
drug-response profile shows strongly diminishing returns after a few hundred variants.** Growing
the panel from 50 to 250 adds about 0.3 to the correlation; growing it from 250 to 1,000, which is four
times the wet-lab cost, adds 0.04 to 0.11 depending on how good the sequence representation is.
The ceiling is set by representation quality, not panel size. So the first screen should be small
in variants and wide in drugs.
→ `analysis/03_simulation.py`, `figures/03_simulation.png`

## What I would build

Put a few hundred real patient mutations into cells as one pooled library, screen the existing
anti-seizure formulary across it, and train a model on the resulting drug-by-mutation map that
places the next mutation without assaying it. Every drug in the first screen is already approved
and already prescribed to these children, so the result is usable immediately rather than after a
development decade.

Cystic fibrosis is the precedent and not the plan. CF established that matching drug to mutation
class is clinically real and commercially large, but it did so one variant per well, in a single
gene that only ever breaks in one direction, to widen the label of molecules the company owned.
What is different here: the assay has to be **two-sided** (a selection that kills the cell when
the channel does too much and a different one when it does too little), the variants have to be
**pooled** for the economics to work (this is the unproven step and therefore the first
experiment), the drugs are **already generic**, and the genes have **relatives**, so the real
prize is transfer across genes, which CF never had.

## How far it could reach

- **Neurodevelopmental conditions.** The same genes cause autism and intellectual disability as
  well as seizures, and in several of them the direction of the mutation decides which. 29 of
  these 142 genes sit in a 492-gene neurodevelopmental set (60 of a wider 298-gene epilepsy
  universe). That field has no measurable endpoint at all: of 91 industry autism trials in phase
  2 or 3, none used a physiological primary endpoint. The drug does not transfer; the measurement
  does.
- **Other channel diseases.** The same protein architecture carries disease in the heart, skeletal
  muscle and inner ear, where **148 human ion channel genes carry disease-causing mutations**, and
  almost nothing has been measured there either: 7 variants of 4,058 in cardiac arrhythmia genes.
- **The ceiling, stated as a ceiling.** Of 3,994 human genes where a mutation causes disease,
  **2,503 are enzymes, transporters, receptors, channels or transcription factors**, classes where
  a two-directional selection is at least constructible; 1,491 are scaffolding proteins where it
  probably is not. Each class needs its own assay designed and some will fail. That is a reach,
  not a roadmap. → `analysis/01_landscape.py` reports the per-family census;
  `data/assayability_by_protein_class.csv` holds the per-gene classification.

## Run it

```bash
git clone https://github.com/vaishnavim9/drug-by-mutation-map
cd drug-by-mutation-map
pip install -r requirements.txt
python analysis/01_landscape.py     # landscape and measurement census
python analysis/02_benchmark.py     # the falsification test
python analysis/03_simulation.py    # sizing the first experiment (~1 min)
```

Outputs are written to `figures/` and `results/`. `01` and `02` read only committed data; `03` is
self-contained.

## Repository map

| path | what it is |
|---|---|
| `data/dee_gene_table.csv` | 142 epilepsy genes × 41 fields: genetic association, constraint, ClinVar spectrum, mechanism class |
| `data/audited_gene_universe.csv` | 298 genes the same query returns when anchored on seven disease and phenotype terms rather than three |
| `data/cross_family_signed_coverage.csv` | the 939-of-15,331 census, by disease family |
| `data/kcnq2_signed_labels.csv` | 291 KCNQ2 variants with measured direction, the benchmark's input |
| `data/signed_prediction_benchmark.csv` | reference results of the original benchmark run |
| `data/clinical_trial_landscape.csv` | 213 registered trials across these genes, classified by modality and sponsor |
| `data/assayability_by_protein_class.csv` | per-gene protein class and whether a two-directional selection is constructible |
| `data/matrix_*.csv` | reference results of the original simulation run |
| `docs/METHODS.md` | where every table came from, query by query |
| `docs/FINDINGS.md` | the numbers, with the caveat attached to each |
| `docs/LIMITATIONS.md` | what is wrong with this analysis, in priority order |

## What this is not

Not a company, not a product, not a dataset of my own measurements, and not a validated
prediction method. The benchmark in finding 2 is a **negative** result about the approach I
started with, kept here because it is the reason the plan changed. Read `docs/LIMITATIONS.md`
before quoting anything from this repository.

Code MIT licensed. Derived tables are built from Open Targets, ClinVar, ClinicalTrials.gov,
UniProt, MaveDB and PubMed; see `docs/METHODS.md` for terms and accessions.
