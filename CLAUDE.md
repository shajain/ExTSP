# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Setup

```bash
# Install in editable mode from project root
pip install -e .

# Or using conda environment
conda env create -f environment.yml
conda activate exTSP
```

Data files are managed with DVC. The `data/` directory contains raw, filtered, and processed datasets tracked via `.dvc/`.

## Running Tests

There is no automated test suite. Validation is done via notebooks:
- `tutorial.ipynb` — primary workflow demonstration
- `test.ipynb` — validation including AUC calculations for pathogenic vs benign variant discrimination

## Architecture

ExTSP computes a tissue-specific pathogenicity score (`exTSP = IsoPath × normalizedLR`) for missense variants by integrating MutPred2 pathogenicity predictions with tissue-specific gene expression data (GTEx, BrainSpan).

### Primary User API

**`src/ExTSP/tripletSets_exTSP.py`** — The main interface. `get_tripletSets_with_exTSP()` returns a DataFrame of variant-transcript-tissue triplets with exTSP scores, filtered by disease and clinical significance (P, PLP, B, BLB, VUS) or ASD case/control.

### Score Computation Pipeline

1. **`exTSP_compute/variantAnnotations.py`** — Reads VEP output per chromosome, maps variants to amino acid substitutions → `variant_annotations.tsv`
2. **`exTSP_compute/mutPred2IsoPath.py`** — Processes MutPred2 scores and maps to isoforms (IsoPath)
3. **`exTSP_compute/compute_PTSE.py`** — Computes Protein Tissue Specific Expression via row-wise softmax on z-scored GTEx/BrainSpan expression matrices → normalizedLR
4. **`exTSP_compute/compute_exTSP.py`** — Merges MutPred2, PTSE, and MANE transcript data → `exTSP_gtex.tsv` and `exTSP_brainspan_cortical.tsv`

### Isoform Selection

- **`isoformSelection/IsoformSelection_Variant.py`** — Selects transcript with max exTSP per variant-tissue pair
- **`isoformSelection/IsoformSelection_Gene.py`** — Gene-level isoform selection, computes tissue enrichment ranks

### Phenotype Integration

**`phenotypes/`** — Queries OMIM, MONDO, and MedGen APIs to map diseases to phenotype IDs used to filter ClinVar variants. API keys are configured in `config.py`.

### Variant Extraction

**`extractVariants/extractVariants.py`** — Extracts missense variants from ClinVar filtered by disease phenotypes → `extractedVariants.xlsx` with per-disease sheets.

### Supported Diseases

CM (Cardiomyopathy), PKD (Polycystic Kidney Disease), PH (Pulmonary Hypertension), PCD (Primary Ciliary Dyskinesia), ASD (Autism Spectrum Disorder). All disease definitions and file paths are centralized in **`src/ExTSP/config.py`**.
