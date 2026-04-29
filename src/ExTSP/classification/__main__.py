"""
Run bootstrap AUC for one disease and setting from the command line::

    python -m ExTSP.classification --disease CM --setting 2 --n-bootstrap 50
"""
from __future__ import annotations

import argparse
import sys
import matplotlib.pyplot as plt
import numpy as np

from ExTSP.classification.auc_bootstrap import bootstrap_auc_for_disease
from ExTSP.config import Diseases


_ALL_DISEASES = ["ASD", "CM", "PCD"]


def _setting_type(v: str) -> "int | str":
    if v == "all":
        return v
    return int(v)


def _build_parser() -> argparse.ArgumentParser:
    disease_choices = sorted(Diseases.keys()) + ["all"]
    p = argparse.ArgumentParser(
        description="Bootstrap AUC (IsoPath and exTSP) for variant–isoform–tissue triplets.",
    )
    p.add_argument(
        "--disease",
        required=True,
        choices=disease_choices,
        metavar="NAME",
        help=f"Disease sheet: {', '.join(disease_choices[:-1])}; or 'all' for {_ALL_DISEASES}.",
    )
    p.add_argument(
        "--setting",
        type=_setting_type,
        required=True,
        choices=[1, 2, 3, 4, 5, "all"],
        metavar="N",
        help="ClinVar / non-target configuration (1–5), or 'all' to run every setting.",
    )
    p.add_argument(
        "--n-bootstrap",
        type=int,
        default=100,
        metavar="N",
        help="Number of bootstrap replicates (default: 100).",
    )
    p.add_argument(
        "--random-state",
        type=int,
        default=None,
        metavar="SEED",
        help="RNG seed for reproducibility.",
    )
    p.add_argument(
        "--balanced-negative",
        action="store_true",
        help="Balance benign negative variants to match pathogenic negative count when possible.",
    )
    p.add_argument(
        "--multivar-genes-only",
        action="store_true",
        help="Restrict to genes with at least 3 distinct pathogenic variants before bootstrapping.",
    )
    p.add_argument(
        "--use-old-data",
        action="store_true",
        help="Load triplets from the legacy data source.",
    )
    p.add_argument(
        "--use-old-data-overlap",
        action="store_true",
        help=(
            "Load triplets from both the legacy and current data sources and keep only "
            "rows present in both (inner merge on Gene, Variant, Transcript_id, and Tissue)."
        ),
    )
    p.add_argument(
        "--non-overlap",
        action="store_true",
        help=(
            "Load triplets from both the legacy and current data sources and keep only "
            "rows exclusive to one source (symmetric difference on Gene, Variant, Transcript_id, and Tissue)."
        ),
    )
    p.add_argument(
        "--old-target-only",
        action="store_true",
        help=(
            "Load target triplets from the legacy data source and non-target triplets "
            "from the current data source."
        ),
    )
    p.add_argument(
        "--split-evaluation",
        action="store_true",
        help=(
            "In addition to the combined evaluation, run two separate evaluations: "
            "one using only pathogenic variants as negatives and one using only benign variants as negatives."
        ),
    )
    p.add_argument(
        "--print-tissue",
        action="store_true",
        help="Print value counts of best_tissue selected across bootstrap replicates.",
    )
    return p

def plot_auc_comparison(results):
    """
    Plot AUC comparison boxplots for different conditions and scenarios.
    
    Parameters:
    -----------
    results : dict
        Nested dictionary with structure:
        {
            'condition_name': {
                'scenario_number': {
                    'MutPred2': list of AUC values,
                    'exTSP': list of AUC values
                }
            }
        }
        Expected conditions: 'ASD', 'Cardiomyopathy', 'PCD'
    
    Returns:
    --------
    fig, axes : matplotlib figure and axes objects
    """
    
    # Create figure with 3 subplots
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # Define colors for the two methods
    colors = {'MutPred2': '#60B85F', 'exTSP': '#FAA43A'}
    
    # Panel labels
    panel_labels = ['A', 'B', 'C']
    conditions = ['ASD', 'Cardiomyopathy', 'PCD']
    full_names = ['Autism Spectrum Disorder (ASD)', 'Cardiomyopathy (CM)', 'Primary Ciliary Dyskinesis (PCD)']
    
    for idx, (condition, full_name) in enumerate(zip(conditions, full_names)):
        ax = axes[idx]
        
        # Get scenarios for this condition
        scenarios = sorted(results[condition].keys(), key=int)
        n_scenarios = len(scenarios)
        
        # Prepare data for boxplot
        positions = []
        data_to_plot = []
        colors_list = []
        
        for i, scenario in enumerate(scenarios):
            # MutPred2
            positions.append(i * 3)
            data_to_plot.append(results[condition][scenario]['MutPred2'])
            colors_list.append(colors['MutPred2'])
            
            # exTSP
            positions.append(i * 3 + 1)
            data_to_plot.append(results[condition][scenario]['exTSP'])
            colors_list.append(colors['exTSP'])
        
        # Create boxplot
        bp = ax.boxplot(data_to_plot, positions=positions, widths=0.6,
                         patch_artist=True, showfliers=False,
                         boxprops=dict(linewidth=1.5, edgecolor='grey'),
                         whiskerprops=dict(linewidth=1.5, color='grey'),
                         capprops=dict(linewidth=1.5, color='grey'),
                         medianprops=dict(linewidth=2, color='grey'))
        
        # Color the boxes
        for patch, color in zip(bp['boxes'], colors_list):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        
        # Add mean values as text above each boxplot
        for i, (pos, data) in enumerate(zip(positions, data_to_plot)):
            mean_val = np.mean(data)
            ax.text(pos, mean_val + 0.02, f'{mean_val:.2f}', 
                    ha='center', va='bottom', fontsize=8, fontweight='bold')
        
        # Set x-axis labels
        ax.set_xticks([i * 3 + 0.5 for i in range(n_scenarios)])
        ax.set_xticklabels([f'Scenario {s}' for s in scenarios])
        
        # Set y-axis
        ax.set_ylabel('AUC' if idx == 0 else '')
        ax.set_ylim(0.0, 1.0)
        ax.set_yticks(np.arange(0.0, 1.1, 0.2))
        
        # Add title with panel label
        ax.text(-0.15, 1.05, panel_labels[idx], transform=ax.transAxes,
                fontsize=16, fontweight='bold', va='top')
        ax.set_title(full_name, fontsize=12, fontweight='bold', pad=20)
        
        # Add grid
        ax.grid(axis='y', alpha=0.3, linestyle='--', linewidth=0.5)
        ax.set_axisbelow(True)
        
        # Add legend only to the first panel at the top
        if idx == 0:
            from matplotlib.patches import Patch
            legend_elements = [Patch(facecolor=colors['MutPred2'], alpha=0.7, label='MutPred2'),
                              Patch(facecolor=colors['exTSP'], alpha=0.7, label='exTSP')]
            ax.legend(handles=legend_elements, loc='upper center', frameon=True, ncol=2)
    
    # Add x-axis label for middle panel
    axes[1].set_xlabel('Classification Scenarios', fontsize=12)
    
    plt.tight_layout()
    return fig, axes


_DISEASE_TO_CONDITION = {'ASD': 'ASD', 'CM': 'Cardiomyopathy', 'PCD': 'PCD'}


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    diseases = _ALL_DISEASES if args.disease == "all" else [args.disease]
    settings = [1, 2, 3, 4, 5] if args.setting == "all" else [args.setting]

    generate_plot = args.disease == "all" and args.setting == "all"
    plot_results = {cond: {} for cond in _DISEASE_TO_CONDITION.values()} if generate_plot else None

    for disease in diseases:
        for setting in settings:
            results, summary = bootstrap_auc_for_disease(
                disease,
                setting,
                balanced_negative=args.balanced_negative,
                split_evaluation=args.split_evaluation,
                filter_min_pathogenic=args.multivar_genes_only,
                n_bootstrap=args.n_bootstrap,
                random_state=args.random_state,
                use_old_data=args.use_old_data,
                use_old_data_overlap=args.use_old_data_overlap,
                non_overlap=args.non_overlap,
                old_target_only=args.old_target_only,
            )
            print(f"disease={disease!r} setting={setting} n_bootstrap={args.n_bootstrap}")
            print(f"combined:           isopath={summary['mean_auc_isopath']:.4f}  extsp={summary['mean_auc_extsp']:.4f}")
            if args.split_evaluation:
                print(f"pathogenic-neg only: isopath={summary['mean_auc_isopath_pathneg']:.4f}  extsp={summary['mean_auc_extsp_pathneg']:.4f}")
                print(f"benign-neg only:     isopath={summary['mean_auc_isopath_benneg']:.4f}  extsp={summary['mean_auc_extsp_benneg']:.4f}")
            if args.print_tissue and results:
                import pandas as pd
                counts = pd.Series([r.best_tissue for r in results]).value_counts()
                tissue_str = ", ".join(f"{tissue}-{n}" for tissue, n in counts.items())
                print(f"best_tissue: {tissue_str}")
            print(f"completed {len(results)} replicates")

            if plot_results is not None and results:
                condition = _DISEASE_TO_CONDITION[disease]
                plot_results[condition][str(setting)] = {
                    'MutPred2': [r.auc_isopath for r in results if np.isfinite(r.auc_isopath)],
                    'exTSP': [r.auc_extsp for r in results if np.isfinite(r.auc_extsp)],
                }

    if generate_plot:
        fig, _ = plot_auc_comparison(plot_results)
        pdf_path = "auc_comparison--old-data-new-nontarget.pdf"
        fig.savefig(pdf_path, bbox_inches='tight')
        plt.close(fig)
        print(f"Saved AUC comparison plot to {pdf_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
