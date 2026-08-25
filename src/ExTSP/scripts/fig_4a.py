# %%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# data_path = "/Users/karnamendonca-kamath/Desktop/ex-TSP/data/exTSP_file_expression_ensemble_annotation_ASD_Satt_Fu_BrainSpan.tsv"
# data_path = "/Users/karnamendonca-kamath/Desktop/ex-TSP/data/exTSP_score_beta_1_Cortical_Chau11.tsv"
data_dir = "/Users/karnamendonca-kamath/Desktop/exTSP2/ExTSP/src/ExTSP/data/brainspan/"
brainspan_path = data_dir + "ASD_brainspan_by_period_no_P13.tsv"
df = pd.read_csv(brainspan_path, sep="\t")
region = "Cortex"
df = df[df['Region'] == region] 
df

# %%
def select_isoform_brainspan(df):
    '''
    Select a single isoform for each variant within each developmental period.

    For every ['Gene', 'Variant', 'Transcript_id', 'Developmental_Period']
    combination, this reduces it to a single
    ['Gene', 'Variant', 'Developmental_Period'] combination by keeping the
    isoform (Transcript_id) whose row has the highest "exTSP" value.
    '''
    return df.loc[df.groupby(["Gene", "Variant", "Developmental_Period"])["exTSP"].idxmax()]

# Split into prenatal (P1-P7) and postnatal (P8-P12), then select the
# Ptse-select isoform independently within each.
prenatal_periods = ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7']
postnatal_periods = ['P8', 'P9', 'P10', 'P11', 'P12']

prenatal_df = df[df['Developmental_Period'].isin(prenatal_periods)]
postnatal_df = df[df['Developmental_Period'].isin(postnatal_periods)]

prenatal_df = select_isoform_brainspan(prenatal_df)
postnatal_df = select_isoform_brainspan(postnatal_df)

df = pd.concat([prenatal_df, postnatal_df], ignore_index=True)
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
    all_tissues = df['Developmental_Period'].unique()
    prenatal_periods = [p for p in ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7'] if p in all_tissues]
    postnatal_periods = [p for p in ['P8', 'P9', 'P10', 'P11', 'P12'] if p in all_tissues]
    all_labels = individual_periods + ['Prenatal', 'Postnatal']

    asd_df = df[df['Status'] == 'ASD']
    ctrl_df = df[df['Status'] == 'Control']

    def get_stats(subset_df, periods, aggregate=False):
        data = subset_df[subset_df['Developmental_Period'].isin(periods if isinstance(periods, list) else [periods])]
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
    ax.set_xlabel('Developmental Period', fontsize=22, labelpad=-8)
    ax.set_ylabel('Mean exTSP', fontsize=22)
    ax.set_title(f'Mean exTSP by Developmental Period', fontsize=26)
    ax.set_xticks(x)
    ax.set_xticklabels(all_labels, fontsize=18, rotation=45, ha='right')
    ax.tick_params(axis='y', labelsize=18)
    ax.legend(fontsize=18)

    plt.tight_layout()
    plt.show()
    return fig

fig_dir = "/Users/karnamendonca-kamath/Desktop/exTSP2/ExTSP/src/ExTSP/figures/"
fig = plot_extsp_isopath(df, agg='avg')
fig.savefig(fig_dir + f'fig_4a.pdf', bbox_inches='tight')

# %%
