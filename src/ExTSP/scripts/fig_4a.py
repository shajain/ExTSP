# %%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# data_path = "/Users/karnamendonca-kamath/Desktop/ex-TSP/data/exTSP_file_expression_ensemble_annotation_ASD_Satt_Fu_BrainSpan.tsv"
data_path = "/Users/karnamendonca-kamath/Desktop/ex-TSP/data/exTSP_score_beta_1_Cortical_Chau11.tsv"
df = pd.read_csv(data_path, sep="\t")
df

# %%
def select_isoform(df):
    '''
    For each variant we consider the isoform with the highest value of Ptse(i, t) as the most 
    important isoform for the variant. We refer to this isoform as Ptse-select isoform.

    Returns a reduced version of the dataframe with the selected isoform for each variant
    '''
    # Step 1: Calculate the max ptse score for each Variant-Transcript_id combination
    transcript_scores = df.groupby(['Variant', 'Transcript_id'])['ptse'].max().reset_index()
    transcript_scores.columns = ['Variant', 'Transcript_id', 'max_ptse']

    # Step 2: For each Variant, find the Transcript_id with the highest max_ptse
    # If there's a tie, sample(1) will pick one randomly
    best_transcripts = (transcript_scores
                        .groupby('Variant')
                        .apply(lambda x: x.loc[x['max_ptse'] == x['max_ptse'].max()]
                            .sample(1))
                        .reset_index()
                        [['Variant', 'Transcript_id']])

    # Step 3: Filter the original dataframe to keep only rows with the selected Variant-Transcript_id pairs
    df_reduced = df.merge(best_transcripts, on=['Variant', 'Transcript_id'], how='inner')
    return df_reduced

df = select_isoform(df)
df

# %%
def plot_extsp_isopath(df, agg="avg"):
    """
    Plot mean exTSP (bars) by developmental period,
    separated by ASD missense and Control missense.
    Prenatal = avg/max over P1-P7, Postnatal = avg/max over P8-P12.

    agg: "avg" uses mean across periods, "max" uses max across periods
         (applies to the Prenatal and Postnatal aggregate bars)
    """
    assert agg in ["avg", "max"]

    individual_periods = ['P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10', 'P11', 'P12']
    all_tissues = df['Tissue'].unique()
    prenatal_periods = [p for p in ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7'] if p in all_tissues]
    postnatal_periods = [p for p in ['P8', 'P9', 'P10', 'P11', 'P12'] if p in all_tissues]
    all_labels = individual_periods + ['Prenatal', 'Postnatal']

    asd_df = df[df['ClinVar_annotation'] == 'ASD']
    ctrl_df = df[df['ClinVar_annotation'] == 'Control']

    def get_stats(subset_df, periods, aggregate=False):
        data = subset_df[subset_df['Tissue'].isin(periods if isinstance(periods, list) else [periods])]
        if aggregate and agg == "max":
            # Per variant, take the max exTSP across periods, then mean/sem over variants
            per_variant = data.groupby('Variant')['exTSP'].max()
            return per_variant.mean(), per_variant.sem()
        return data['exTSP'].mean(), data['exTSP'].sem()

    asd_stats = [get_stats(asd_df, [p]) for p in individual_periods] + \
                [get_stats(asd_df, prenatal_periods, aggregate=True),
                 get_stats(asd_df, postnatal_periods, aggregate=True)]
    ctrl_stats = [get_stats(ctrl_df, [p]) for p in individual_periods] + \
                 [get_stats(ctrl_df, prenatal_periods, aggregate=True),
                  get_stats(ctrl_df, postnatal_periods, aggregate=True)]

    asd_extsp_means, asd_extsp_sems = zip(*asd_stats)
    ctrl_extsp_means, ctrl_extsp_sems = zip(*ctrl_stats)

    x = np.arange(len(all_labels))
    width = 0.35
    asd_color = '#CC5500'
    ctrl_color = '#3A9E60'

    fig, ax = plt.subplots(figsize=(14, 6))

    ax.bar(x - width / 2, asd_extsp_means, width, color=asd_color,
           yerr=asd_extsp_sems, capsize=3, label='ASD missense')
    ax.bar(x + width / 2, ctrl_extsp_means, width, color=ctrl_color,
           yerr=ctrl_extsp_sems, capsize=3, label='Control missense')

    agg_label = "Mean" if agg == "avg" else "Max"
    ax.set_xlabel('Developmental Period')
    ax.set_ylabel('Mean exTSP')
    ax.set_title(f'Mean exTSP by Developmental Period\n+ Prenatal and Postnatal {agg_label}s')
    ax.set_xticks(x)
    ax.set_xticklabels(all_labels)
    ax.legend()

    plt.tight_layout()
    plt.show()
    return fig

fig_dir = "/Users/karnamendonca-kamath/Desktop/exTSP2/ExTSP/src/ExTSP/figures/"
for agg in ["avg", "max"]:
    fig = plot_extsp_isopath(df, agg=agg)
    fig.savefig(fig_dir + f'fig_4a_{agg}.pdf', bbox_inches='tight')

# %%
