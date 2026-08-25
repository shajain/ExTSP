# %%
### Import Statements
import pandas as pd
import matplotlib.pyplot as plt
from ExTSP.tripletSets_exTSP import get_tripletSets_with_exTSP

# %%
### Load file
data_path = "/Users/karnamendonca-kamath/Desktop/exTSP2/ExTSP/src/ExTSP/data/brainspan/ASD_brainspan_by_period_no_P13.tsv"
# data_path = "/Users/karnamendonca-kamath/Desktop/exTSP2/ExTSP/src/ExTSP/data/brainspan/ASD_brainspan_cortical_by_period.tsv"
df = pd.read_csv(data_path, sep='\t')
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


# Split into prenatal (P1-P7) and postnatal (P8-P12), then by Region, and select
# the Ptse-select isoform independently within each Region x stage group.
prenatal_periods = ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7']
postnatal_periods = ['P8', 'P9', 'P10', 'P11', 'P12']

df['Stage'] = df['Developmental_Period'].apply(
    lambda p: 'Prenatal' if p in prenatal_periods else 'Postnatal'
)

df = pd.concat(
    [select_isoform_brainspan(group_df) for _, group_df in df.groupby(["Region", "Stage"])],
    ignore_index=True,
)

### Plot average exTSP per region, split by stage (prenatal/postnatal) and
### case (ASD) vs control missense variants -> four bars per region.
def stage_status_means(stage, status):
    return df[(df["Stage"] == stage) & (df["Status"] == status)].groupby("Region")["exTSP"].mean()

prenatal_case = stage_status_means("Prenatal", "ASD")
prenatal_control = stage_status_means("Prenatal", "Control")
postnatal_case = stage_status_means("Postnatal", "ASD")
postnatal_control = stage_status_means("Postnatal", "Control")

# Order regions by overall (prenatal+postnatal) case mean, descending
region_order = df[df["Status"] == "ASD"].groupby("Region")["exTSP"].mean().sort_values(ascending=False).index

prenatal_case = prenatal_case.reindex(region_order)
prenatal_control = prenatal_control.reindex(region_order)
postnatal_case = postnatal_case.reindex(region_order)
postnatal_control = postnatal_control.reindex(region_order)

x = range(len(region_order))
width = 0.2
prenatal_control_color = "#88C9A1"   # light green
prenatal_case_color = "#E8A06A"      # light orange
postnatal_control_color = "#3A9E60"  # dark green
postnatal_case_color = "#CC5500"     # dark orange

fig, ax = plt.subplots(figsize=(14, 6))
ax.bar([i - 1.5 * width for i in x], prenatal_case, width=width,
       color=prenatal_case_color, label="Prenatal case")
ax.bar([i - 0.5 * width for i in x], postnatal_case, width=width,
       color=postnatal_case_color, label="Postnatal case")
ax.bar([i + 0.5 * width for i in x], prenatal_control, width=width,
       color=prenatal_control_color, label="Prenatal control")
ax.bar([i + 1.5 * width for i in x], postnatal_control, width=width,
       color=postnatal_control_color, label="Postnatal control")

ax.set_xticks(list(x))
ax.set_xticklabels(list(region_order), rotation=90)
ax.set_xlabel("Region")
ax.set_ylabel("Average exTSP")
ax.set_title(
    "Average exTSP of ASD Case vs. Control Missense Variants by Region\n"
    "(split by prenatal vs. postnatal)",
    fontsize=12,
)
ax.legend()
plt.tight_layout()
plt.show()

fig_dir = "/Users/karnamendonca-kamath/Desktop/exTSP2/ExTSP/src/ExTSP/figures/"
fig.savefig(fig_dir + 'fig_4b.pdf', bbox_inches='tight')


# %%
