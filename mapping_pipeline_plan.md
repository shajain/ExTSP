# Tiered Ensembl-113 transcript remapping pipeline (ID history → structural match → SQANTI3 categorization)

## Context

The lab's transcript-level data spans three annotation eras:
- **GTEx**: GENCODE v26 (= Ensembl 89), genome build **GRCh38** — same primary assembly as current Ensembl, just an older annotation.
- **BrainSpan**: GENCODE v10 (= Ensembl 65), genome build **GRCh37/hg19** — an older assembly *and* an older annotation.
- **Target**: Ensembl 113 (GRCh38.p14).

The collaborator's existing file (`gtex_brainspan_transcripts_mapped_to_ensembl_v113.txt`, 162,222 rows) was built by matching on gene + exon-coordinate similarity (it carries `max_exon_change`/`exon_difference` diagnostic columns), and `ensembl_mapping.py` only *consumes* it for QC — the script that produced it isn't in this workspace and its logic can't be inspected or extended. The PI considers its data loss unacceptable and wants a specific, more principled tiered pipeline instead:

1. **Tier 1 — ID history**: versionless direct match, then Ensembl's ID History Converter for anything retired/merged/split. She expects this alone to recover ~80% of both sets.
2. **Tier 2 — structural match**: gffcompare class code `=` (identical internal splice junctions, ignoring TSS/TES drift) — the SQANTI3 "Full Splice Match" definition.
3. **Tier 3 — SQANTI3 categorization** of whatever remains (ISM/NIC/NNC), to document *why* a transcript couldn't be matched rather than just dropping it.

Environment check (2026-08-25, run locally on `/Users/karnamendonca-kamath/Desktop/exTSP2`) confirmed: conda is available, network access to `rest.ensembl.org` and `ftp.ensembl.org` works, 347GB free disk. No bioinformatics tools (gffcompare, minimap2, samtools, CrossMap, SQANTI3) are currently installed, and no GTFs/chain files exist locally — everything needs to be fetched/built fresh. This is a green-field build, not a modification of `ensembl_mapping.py`'s current logic.

**Open question before implementation starts:** whether to run this pipeline locally or on an SSH server (deferred — no install/run should happen until that's decided).

## Estimated cost

**Disk space: ~10–15 GB total**
- GENCODE v26 GTF (gzip ~40MB / uncompressed ~1.4GB)
- GENCODE v10 GTF (much smaller, ~30MB / ~200MB uncompressed — 2012-era annotation)
- Ensembl 113 GTF (~50MB / ~1.5GB uncompressed)
- GRCh38 primary-assembly genome FASTA (~850MB compressed / ~3.1GB uncompressed) — needed only for Tier 3 (SQANTI3)
- hg19→hg38 liftover chain (~40MB compressed)
- SQANTI3's own conda env + cDNA_Cupcake + R deps (~2–5GB)
- Intermediate/output files (lifted GTFs, gffcompare/.tmap files, SQANTI3 outputs): low hundreds of MB

**Time: roughly 1–2 hours end-to-end**, front-loaded into setup rather than the actual matching:
- Environment setup: ~10–20 min for the main conda env; SQANTI3 install is the biggest wildcard — it builds cDNA_Cupcake and pulls R packages, historically the most failure-prone step, could be 20–40 min or need a retry.
- Reference downloads: a few minutes (GTFs are small; the genome FASTA is the biggest single file at ~850MB).
- Tier 1 (ID history): fast — even in the worst case of ~50k unmatched IDs, that's ~50 batched REST calls, a few minutes.
- Tier 2 (liftover + gffcompare): a few minutes; gffcompare against a ~250k-transcript reference is not expensive.
- Tier 3 (SQANTI3): the heaviest per-transcript step, but it only runs on whatever's left after Tiers 1–2 — if that remainder is genuinely small (hundreds to low thousands, per the PI's 80% expectation), likely 5–20 min.

## Plan

All new work goes in a new `mapping_pipeline/` directory (sibling to `ensembl_mapping.py`), so the existing script and the collaborator's file are left untouched for comparison.

### 1. Environment setup
- `mapping_pipeline/environment.yml`: conda env `extsp-mapping` with bioconda+conda-forge channels: `gffcompare`, `gffread`, `minimap2`, `samtools`, `crossmap`, `gffutils`, `pysam`, `pandas`, `requests`.
- SQANTI3 installed separately per its own upstream instructions (github.com/ConesaLab/SQANTI3) into a dedicated env — it pulls in cDNA_Cupcake and R dependencies that don't mix cleanly with the main env. Only invoked in Tier 3.

### 2. Reference data acquisition (`mapping_pipeline/references/`, via a `00_download_references.sh`)
- GENCODE v26 GTF (GRCh38) — query annotation for GTEx, no liftover needed.
- GENCODE v10 GTF (GRCh37/hg19) — query annotation for BrainSpan.
- Ensembl 113 GTF (GRCh38.p14) — the target reference for gffcompare/SQANTI3.
- GRCh38 primary-assembly genome FASTA — required by SQANTI3 (Tier 3 only).
- UCSC `hg19ToHg38.over.chain.gz` — for CrossMap-ing the BrainSpan/v10 GTF from hg19 → GRCh38 before Tier 2/3 (GTEx skips this step entirely since v26 is already GRCh38).

### 3. Query transcript universe
Extract the source ID list directly from the collaborator's existing file's `transcript_id` + `dataset` columns (161,289 GTEx / 932 BrainSpan) — this is the full known universe of transcripts needing a v113 mapping, so there's no need to re-derive it from the `ExTSP` triplet-set data.

### 4. Tier 1 — `01_tier1_id_history.py`
- Build the versionless v113 transcript ID set from the Ensembl 113 GTF.
- Strip version suffixes from source IDs and match directly.
- For unmatched IDs: batch-query the Ensembl REST `/archive/id/` endpoint (POST, up to 1000 IDs/request — confirmed reachable and rate-limit-friendly at 54,994/hour remaining). Resolve `possible_replacement` entries (handles retirements/merges/splits) against the v113 ID set.
- Output: `tier1_matches.tsv` (source_id, v113_id, method), `tier1_unmatched_{gtex,brainspan}.txt`.

### 5. Tier 2 — structural match via gffcompare
- `02_liftover_brainspan.sh`: CrossMap the GENCODE v10 GTF from hg19 → GRCh38 (BrainSpan only).
- `03_tier2_gffcompare.sh`: filter each query GTF (v26 as-is / lifted v10) down to Tier-1-unmatched IDs (via `gffread`), then run `gffcompare -r ensembl113.gtf <query.gtf>` per dataset.
- `04_parse_tmap.py`: parse the `.tmap` output, keep rows with `class_code == '='`, map query_id → ref_id (v113 transcript).
- Output: `tier2_matches.tsv`, `tier2_unmatched_{gtex,brainspan}.txt`.

### 6. Tier 3 — SQANTI3 categorization
- `05_tier3_sqanti3.sh`: run `sqanti3_qc.py` on the still-unmatched query GTFs against the Ensembl 113 GTF + GRCh38 genome FASTA.
- `06_parse_sqanti_classification.py`: pull `structural_category` (FSM/ISM/NIC/NNC/genic/intergenic/...) and `associated_transcript` per transcript — ISM cases get a best-effort reference ID even without a `=` class code.
- Output: `tier3_categories.tsv`.

### 7. Merge + report
- `07_merge_and_report.py`: combine all three tiers into one new table, `mapping_pipeline/results/gtex_brainspan_mapped_v113_v2.tsv`, adding `match_tier` and `match_method`/`structural_category` columns on top of the schema `ensembl_mapping.py` already expects (`transcript_id`, `dataset`, `transcript_id_v113`, ...), so it's a drop-in alternative input.
- Reuse `summarize_sets`/`print_stats` from `ensembl_mapping.py` (confirmed present, lines 21-70) to compute the same before/after retention stats against the new table, broken down by tier and by disease category from the `ExTSP` triplet sets. Note: `ensembl_mapping.py`'s existing join (line 94-100) joins `df["Transcript"]` (versioned) against `mapping["transcript"]`, despite its own comment saying it should use the versioned `transcript_id` column — the new report script will join on the correct versioned `transcript_id` column explicitly rather than reproducing that mismatch.

## Verification
- After Tier 1: report % of each dataset resolved by direct match vs. ID-history match — should be in the ballpark of the PI's ~80% expectation.
- After Tier 2: report additional % resolved via `=` class code, spot-check a handful of `.tmap` `=` rows manually against Ensembl to confirm TSS/TES-only differences.
- After Tier 3: report the ISM/NIC/NNC breakdown for the remainder — no transcript should end up simply "unmatched" without a category.
- Final: run the existing `print_stats`/`summarize_sets` QC logic against `gtex_brainspan_mapped_v113_v2.tsv` for every disease category and confirm materially less data loss than the collaborator's file, in a form the PI can review directly.
