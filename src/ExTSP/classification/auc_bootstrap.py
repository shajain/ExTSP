"""
Bootstrap gene-wise AUC for pathogenic vs benign variant–isoform–tissue triplets.

Expected columns (defaults in :mod:`ExTSP.classification.columns`): ``Gene``,
``Transcript``, ``Transcript_id`` (isoform), ``Variant``, ``MutPred2``, ``IsoPath``,
``ptse``, ``exTSP``, ``meanTPM``, ``Tissue``.

Bootstrap draws genes with replacement from the **union** of genes across
pathogenic target, benign target, and non-target arms; out-of-bag genes are those in
that union never drawn. No extra overlap filter is applied afterward.

Positive / negative sets follow the user-specified construction: tissue is chosen
from pathogenic-target rows restricted to bootstrapped genes; isoforms are chosen
from rows whose gene is out-of-bag (pathogenic target, benign target, and non-target
arms as applicable); positives use pathogenic target at the chosen tissue; negatives
pool benign and non-target arms plus pathogenic-target rows at other tissues.

Tissue is :func:`~ExTSP.classification.pickers.pick_tissue_mean_exTSP` on bootstrapped
pathogenic target; isoforms use :func:`~ExTSP.isoformSelection.IsoformSelection_Variant.exTSP_selected_isoform`.
"""
from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

#from ExTSP.classification.clinvar_splits import split_four_tables
#from ExTSP.classification.triplet_data import load_triplet_dataframes
from ExTSP.tripletSets_exTSP import get_tripletSets_with_exTSP, variantCountsWithDuplicates, sample_tripletSets
from ExTSP.isoformSelection.IsoformSelection_Variant import exTSP_selected_isoform
from ExTSP.config import VariantSets

from . import columns as C
from ExTSP.tissueEnrichment.tissueEnrichment import select_disease_tissue


@dataclass
class BootstrapAUCResult:
    """One bootstrap draw: AUC from IsoPath and from exTSP (may be nan if undefined).

    ``n_oob_variants`` counts distinct variants in the OOB gene slice across pathogenic
    target, benign target, and non-target arms when present.
    """
    auc_isopath: float
    auc_extsp: float
    best_tissue: str
    # n_oob_variants_pt: int
    # n_oob_variants_bt: int
    # n_oob_variants_pn: int
    # n_oob_variants_bn: int
    n_oob_variants_pos: int
    n_oob_variants_neg: int
    n_oob_genes_pos: int
    n_oob_genes_neg: int
    



def _genes_unique_union(
    dfs: Sequence[pd.DataFrame | None],
    gene_col: str,
) -> np.ndarray:
    """Unique ``gene_col`` values across non-empty dataframes (order: first-seen)."""
    parts = [df[gene_col].unique().tolist() for df in dfs]
    #flatten the list
    parts = [item for sublist in parts for item in sublist]
    parts = list(set(parts))
    return parts


def _bootstrap_draw_and_oob(
    genes_boot: np.ndarray,
    genes_all: np.ndarray,
    rng: np.random.Generator | None = None,
    *,
    prop_gene_draws: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Sample ``n_draws`` genes with replacement from ``unique_genes``.
    In-bag genes are those that appear at least once in the draw; final OOB genes are
    those in ``unique_genes`` that never appear in the draw (no further overlap checks).
    """
    genes_boot = np.array(genes_boot, dtype=object)
    genes_all = np.array(genes_all, dtype=object)
    assert len(genes_boot) > 1, "bootstrap genes are empty or only one gene"
    prop_gene_draws = prop_gene_draws if prop_gene_draws is not None else 1.0
    while True:
        n_draws = np.ceil(prop_gene_draws * len(genes_boot)).astype(int)
        draw = rng.choice(genes_boot, size=n_draws, replace=True)
        seen = np.unique(draw)
        oob = genes_all[~np.isin(genes_all, seen)]
        if genes_boot[~np.isin(genes_boot, seen)].size > 0:
            break
    return draw, oob

def filter_by_gene(
    df: pd.DataFrame,
    genes: np.ndarray,
    gene_col: str,
) -> pd.DataFrame:
    """Rows whose ``gene_col`` is in ``oob_genes``."""
    return df.loc[df[gene_col].isin(genes)]

def build_positive_negative_triplets(
    best_tissue: Any,
    path_target: pd.DataFrame,
    ben_target: pd.DataFrame,
    path_nontarget: pd.DataFrame | None,
    ben_nontarget: pd.DataFrame | None,
    oob_genes: np.ndarray,
    gene_col: str,
    *,
    balanced_negative: bool = False,
    tissue_col: str = C.TISSUE,
    isoform_col: str = C.TRANSCRIPT_ID,
    variant_col: str = C.VARIANT,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Positive: pathogenic target rows with selected isoform at ``best_tissue``.

    Negative (depends on which arms are non-empty after ``split_four_tables``):

    * Benign target / benign non-target: selected isoform, **all tissues** (empty in
      setting 3, which has **no** benign variants).
    * Pathogenic non-target: selected isoform, **all tissues**.
    * Pathogenic target: selected isoform, **all tissues other than** ``best_tissue``.

    Setting 3: only PLP pathogenics; negatives are non-target pathogenics (all
    tissues) plus pathogenic target at non-selected tissues—no B/LB rows.
    """
    # pt_oob, bt_oob, pn_oob, bn_oob = create_oob(path_target, ben_target, path_nontarget, ben_nontarget, oob_genes, gene_col)
    # iso_pt, iso_bt, iso_pn, iso_bn = select_isoforms(pt_oob, bt_oob, pn_oob, bn_oob)
    iso_pt = exTSP_selected_isoform(filter_by_gene(path_target, oob_genes, gene_col), isoform_col)
    iso_bt = exTSP_selected_isoform(filter_by_gene(ben_target, oob_genes, gene_col), isoform_col)
    iso_pn = exTSP_selected_isoform(filter_by_gene(path_nontarget, oob_genes, gene_col), isoform_col)
    iso_bn = exTSP_selected_isoform(filter_by_gene(ben_nontarget, oob_genes, gene_col), isoform_col)
    pos = iso_pt[iso_pt[tissue_col] == best_tissue]
    neg_path_target = iso_pt[iso_pt[tissue_col] != best_tissue]
    neg_ben_target = iso_bt
    neg_path_nontarget = iso_pn
    neg_ben_nontarget = iso_bn
    neg_path = pd.concat([neg_path_target, neg_path_nontarget], axis=0, ignore_index=True)
    neg_path = neg_path.drop_duplicates(subset=[variant_col, isoform_col, tissue_col])
    neg_ben = pd.concat([neg_ben_target, neg_ben_nontarget], axis=0, ignore_index=True)
    neg_ben = neg_ben.drop_duplicates(subset=[variant_col, isoform_col, tissue_col])
    if balanced_negative and not neg_ben.empty:
        n_path_variants = neg_path[variant_col].nunique()
        n_ben_variants = neg_ben[variant_col].nunique()
        if n_path_variants != n_ben_variants:
            neg_ben = sample_tripletSets(neg_ben, n_path_variants)
    neg = pd.concat([neg_path, neg_ben], axis=0, ignore_index=True)
    return pos, neg


def _auc_binary(pos: pd.DataFrame, neg: pd.DataFrame, score_col: str) -> float:
    if pos.empty or neg.empty:
        return float("nan")
    y = np.concatenate([np.ones(len(pos)), np.zeros(len(neg))])
    s = pd.concat([pos[score_col], neg[score_col]], ignore_index=True).to_numpy(dtype=float)
    m = np.isfinite(s)
    if m.sum() == 0 or len(np.unique(y[m])) < 2:
        return float("nan")
    return float(roc_auc_score(y[m], s[m]))


def create_bootstrap_result(pos: pd.DataFrame, neg: pd.DataFrame, best_tissue: str, variant_col: str, isopath_col: str, extsp_col: str, gene_col: str) -> BootstrapAUCResult:
    auc_iso = _auc_binary(pos, neg, isopath_col)
    auc_ex = _auc_binary(pos, neg, extsp_col)
    return BootstrapAUCResult(
        auc_isopath=auc_iso,
        auc_extsp=auc_ex,
        best_tissue=best_tissue,
        n_oob_variants_pos=len(pos[variant_col].unique()),
        n_oob_variants_neg=len(neg[variant_col].unique()),
        n_oob_genes_pos=len(pos[gene_col].unique()),
        n_oob_genes_neg=len(neg[gene_col].unique()),
    )

# def create_oob(pt: pd.DataFrame, bt: pd.DataFrame, pn: pd.DataFrame, bn: pd.DataFrame, oob_genes: np.ndarray, gene_col: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
#     pt_oob = filter_by_gene(pt, oob_genes, gene_col)
#     bt_oob = filter_by_gene(bt, oob_genes, gene_col) 
#     pn_oob = filter_by_gene(pn, oob_genes, gene_col) 
#     bn_oob = filter_by_gene(bn, oob_genes, gene_col) 
#     return pt_oob, bt_oob, pn_oob, bn_oob

# def select_isoforms(pt: pd.DataFrame, bt: pd.DataFrame, pn: pd.DataFrame, bn: pd.DataFrame, isoform_col: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
#     iso_pt = exTSP_selected_isoform(pt, isoform_col) 
#     iso_bt = exTSP_selected_isoform(bt, isoform_col) 
#     iso_pn = exTSP_selected_isoform(pn, isoform_col) 
#     iso_bn = exTSP_selected_isoform(bn, isoform_col)
#     return iso_pt, iso_bt, iso_pn, iso_bn

def run_one_bootstrap(
    path_target: pd.DataFrame,
    ben_target: pd.DataFrame,
    path_nontarget: pd.DataFrame,
    ben_nontarget: pd.DataFrame,
    *,
    balanced_negative: bool = False,
    prop_gene_draws: float | None = None,
    rng: np.random.Generator | None = None,
    gene_col: str = C.GENE,
    variant_col: str = C.VARIANT,
    tissue_col: str = C.TISSUE,
    isoform_col: str = C.TRANSCRIPT_ID,
    isopath_col: str = C.ISO_PATH,
    extsp_col: str = C.EXTSP,
    
) -> BootstrapAUCResult:
    
    unique_genes_all = _genes_unique_union([path_target, ben_target, path_nontarget, ben_nontarget], gene_col)
    unique_genes_pos = _genes_unique_union([path_target], gene_col)
    draw, oob_genes = _bootstrap_draw_and_oob(unique_genes_pos, unique_genes_all, rng, prop_gene_draws=prop_gene_draws)

    pt_boot = path_target[path_target[gene_col].isin(draw)]
    best_tissue = select_disease_tissue(pt_boot)

    pos, neg = build_positive_negative_triplets(
        best_tissue, path_target, ben_target, path_nontarget, ben_nontarget, oob_genes, gene_col,
        balanced_negative=balanced_negative, tissue_col=tissue_col, isoform_col=isoform_col, variant_col=variant_col)

    return create_bootstrap_result(pos, neg, best_tissue, variant_col, isopath_col, extsp_col, gene_col)


def bootstrap_auc(
    df_target_p: pd.DataFrame,
    df_target_b: pd.DataFrame,
    df_nontarget_p: pd.DataFrame,
    df_nontarget_b: pd.DataFrame,
    *,
    balanced_negative: bool = False,
    n_bootstrap: int = 100,
    prop_gene_draws: float | None = None,
    random_state: int = 42,
    gene_col: str = C.GENE,
    variant_col: str = C.VARIANT,
    tissue_col: str = C.TISSUE,
    isoform_col: str = C.TRANSCRIPT_ID,
    isopath_col: str = C.ISO_PATH,
    extsp_col: str = C.EXTSP,
) -> tuple[list[BootstrapAUCResult], dict[str, float]]:
    """
    Run ``n_bootstrap`` gene-level bootstrap iterations and return per-iteration
    results plus mean AUCs (ignoring nan).

    Parameters
    ----------
    df_target, df_nontarget
        Tables using default columns from :mod:`ExTSP.classification.columns` (or pass
        overrides via ``*_col`` parameters).
    setting
        1–5 as described in ``split_four_tables`` / project spec.
    n_bootstrap
        Number of bootstrap replicates.
    """
    rng = np.random.default_rng(random_state)
    results: list[BootstrapAUCResult] = []
    for _ in range(n_bootstrap):
        results.append(
            run_one_bootstrap(
                df_target_p,
                df_target_b,
                df_nontarget_p,
                df_nontarget_b,
                prop_gene_draws=prop_gene_draws,
                rng=rng,
                balanced_negative=balanced_negative,
                gene_col=gene_col,
                variant_col=variant_col,
                tissue_col=tissue_col,
                isoform_col=isoform_col,
                isopath_col=isopath_col,
                extsp_col=extsp_col,
            )
        )

    def _mean(vals: Sequence[float]) -> float:
        a = np.asarray(vals, dtype=float)
        a = a[np.isfinite(a)]
        return float(a.mean()) if a.size else float("nan")

    summary = {
        "mean_auc_isopath": _mean([r.auc_isopath for r in results]),
        "mean_auc_extsp": _mean([r.auc_extsp for r in results]),
    }
    return results, summary

def filter_columns(df: pd.DataFrame, columns: list[str]=None) -> pd.DataFrame:
    if columns is None:
        columns = [C.GENE, C.VARIANT, C.TRANSCRIPT_ID, C.TISSUE, C.ISO_PATH, C.EXTSP]
    return df[columns]


def bootstrap_auc_for_disease(
    disease: str,
    setting: int,
    balanced_negative: bool = False,
    n_bootstrap: int = 100,
    prop_gene_draws: float | None = None,
    random_state: int = 42,
) -> tuple[list[BootstrapAUCResult], dict[str, float]]:
    """
    Same as :func:`bootstrap_auc`, but builds ``df_target`` / ``df_nontarget`` from
    ``disease`` and ``setting`` via :func:`~ExTSP.classification.triplet_data.load_triplet_dataframes`
    and :func:`~ExTSP.tripletSets_exTSP.get_tripletSets_with_exTSP`.

    Tissue and isoform selection match :func:`bootstrap_auc` (mean exTSP tissue on
    bootstrapped pathogenic target; ``exTSP_selected_isoform`` per variant).

    Parameters
    ----------
    disease
        Target disease sheet (e.g. ``\"PKD\"``), not ``\"nonTarget\"``.
    setting
        1–5; determines non-target inclusion and ``type`` (``PLP`` vs ``all``) for triplet loading.
    numVars
        Passed to ``get_tripletSets_with_exTSP`` (``brainspan`` is fixed to ``False``).
    **kwargs
        Forwarded to :func:`bootstrap_auc` (e.g. ``n_bootstrap``, ``random_state``).
    """
    assert setting in [1, 2, 3, 4, 5], "setting must be 1, 2, 3, 4, or 5"
    assert disease in VariantSets, f"disease {disease} must be in VariantSets: {VariantSets}"
    if setting in [1,2]:
        typeP = "P" if disease!="ASD" else "case"
        typeB = "B" if disease!="ASD" else "control"
    elif setting in [3, 4, 5]:
        typeP = "PLP"
        typeB = "BLB"
    df_target_p = get_tripletSets_with_exTSP(disease, type=typeP)
    df_target_b = df_target_p.head(0).copy()
    df_nontarget_p = df_target_p.head(0).copy()
    df_nontarget_b = df_target_p.head(0).copy()
    if setting != 3:
        df_target_b = get_tripletSets_with_exTSP(disease, type=typeB)
        if setting in [2,5]:
            df_nontarget_p = get_tripletSets_with_exTSP('nonTarget', type=typeP)
            df_nontarget_b = get_tripletSets_with_exTSP('nonTarget', type=typeB)
    else: 
        df_nontarget_p = get_tripletSets_with_exTSP('nonTarget', type=typeP)
    df_target_p = filter_columns(df_target_p)
    df_target_b = filter_columns(df_target_b)
    df_nontarget_p = filter_columns(df_nontarget_p)
    df_nontarget_b = filter_columns(df_nontarget_b)
    if len(df_target_p[C.GENE].unique()) == 1:
        return [], {"mean_auc_isopath": float("nan"), "mean_auc_extsp": float("nan")}
    return bootstrap_auc(df_target_p, df_target_b, df_nontarget_p, df_nontarget_b,
                            balanced_negative=balanced_negative, 
                            n_bootstrap=n_bootstrap, 
                            prop_gene_draws=prop_gene_draws,
                            random_state=random_state,
                            )
