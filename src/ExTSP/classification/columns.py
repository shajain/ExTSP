"""
Canonical column names for variant–isoform–tissue exTSP tables.

Isoform is represented by ``Transcript_id`` (ENST); ``Transcript`` is the ID without
version where applicable.
"""

GENE = "Gene"
TRANSCRIPT = "Transcript"
TRANSCRIPT_ID = "Transcript_id"  # isoform identifier
VARIANT = "Variant"
MUTPRED2 = "MutPred2"
ISO_PATH = "IsoPath"
PTSE = "ptse"
EXTSP = "exTSP"
MEAN_TPM = "meanTPM"
TISSUE = "Tissue"
CLINVAR_ANNOTATION = "ClinVar_annotation"

# Default deduplication key for triplet rows (variant + isoform + tissue)
DEFAULT_TRIPLET_DEDUP = (VARIANT, TRANSCRIPT_ID, TISSUE)
