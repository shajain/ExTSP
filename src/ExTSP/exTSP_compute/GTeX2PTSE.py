

import pandas as pd
import numpy as np
from scipy.stats import zscore
from scipy.stats import norm
from ExTSP.config import GTeX_file, PROCESSED_DIR, PTSE_file
from pathlib import Path




#PTSE_file = f'{PROCESSED_DIR}/GTEx_PTSE.csv'

def PTSE(df_input):
    """
    Copy the dataframe, apply row-wise softmax to the first 52 columns,
    then divide those columns by their row-wise maximum. Returns the
    modified copy.
    """
    df = df_input.copy()
    df = zScoreNormalize(df)
    # Select first 52 columns as a NumPy array
    x = df.iloc[:, 0:52].to_numpy(dtype=float)
    # Row-wise softmax: subtract row max for numerical stability
    x_max = np.max(x, axis=1, keepdims=True)
    exps = np.exp(x - x_max)
    # compute ptse
    ptse = exps / np.sum(exps, axis=1, keepdims=True)
    ptse = np.maximum(ptse, 1e-10)
    df.iloc[:, 0:52] = ptse
    return df

def ptse2normalizedLR(df_ptse):
    df = df_ptse.copy()
    ptse = df.iloc[:, 0:52].to_numpy(dtype=float)
     # compute lr
    lr = ptse / (1 - ptse)
    # Divide by row-wise max of the softmax values
    lr_max = np.max(lr, axis=1, keepdims=True)
    # Avoid division by zero
    lr_max[lr_max == 0] = 1.0
    # scale by row-wise max of the LR
    normalized_lr = lr / lr_max
    # Put back into dataframe
    df.iloc[:, 0:52] = normalized_lr
    return df

def zScoreNormalize(df_input,robust=False):
    df = df_input.copy()
    if robust:
        q1 = df.iloc[:,0:52].quantile(0.25, axis=1)
        q3 = df.iloc[:,0:52].quantile(0.75, axis=1)
        iqr = q3 - q1
        robust_std = iqr / 1.349
        robust_std[robust_std==0] = 10**-6
        mu = df.iloc[:,0:52].median(axis=1)
        std = robust_std
    else:
        std =  df.iloc[:,0:52].std(axis=1)
    mu = df.iloc[:,0:52].mean(axis=1)
    df.iloc[:,0:52] = df.iloc[:,0:52].sub(mu, axis=0).div(std, axis=0)
    df.iloc[:,0:52] = df.iloc[:,0:52].fillna(0)
    return df

def zScoreNormalizeSignificant(df_input,robust=False):
    df = zScoreNormalize(df_input,robust)
    df.iloc[:,0:52] = 1-norm.cdf(df.iloc[:,0:52])
    df.iloc[:,0:52] = (df.iloc[:,0:52] <=0.001).astype(float)
    df["numSignificant"] =  df.iloc[:,0:52].sum(axis=1)
    return df

def sigTissueCounts(df):
    nTissuesSig = df.iloc[:,0:52].apply(sum, axis=1)
    print("number of transcripts in GTEx:", df.shape[0])
    print("number of transcripts in GTEx where at aleast one tissue has significant expression:", (nTissuesSig>0).sum())
    print(nTissuesSig.value_counts())

def meltTissues_df(df, new_col_name):
     # First 52 columns are tissues; remaining columns are variant/isoform/gene/etc.
    tissue_cols = df.columns[:52]
    meta_cols = df.columns[52:]

    # Long format for PTSE values
    df_long = df.melt(
        id_vars=meta_cols,
        value_vars=tissue_cols,
        var_name="tissue",
        value_name=new_col_name
    )
    return df_long

def main():
    df_raw = pd.read_csv(GTeX_file, sep='\t')
    #Check if file exists

    #df_raw.iloc[:,0:52] = np.log2(df_raw.iloc[:,0:52]+1)
    # df_raw = df_raw.rename(columns={"gene_name": "GeneSymbol"})
    # df_raw["totalTranscripts"] = df_raw.groupby("GeneSymbol").transform("size")
    # df = df_raw.copy()
    # df1 = zScoreNormalizeSignificant(df) 
    # sigTissueCounts(df1)
    # df2 = df_raw.copy()
    # nTranscripts = []
    # expCols = df_raw.columns[0:52]
    # for gene, idx in df2.groupby("GeneSymbol").groups.items():
    #     mu = df2.loc[idx,expCols].mean(axis=0)
    #     df2.loc[idx,expCols] = df2.loc[idx,expCols] - mu
    #     nTranscripts.append(len(idx))
    # df3 = zScoreNormalizeSignificant(df2, robust=False)
    # sigTissueCounts(df3)
    #df1 = zScoreNormalize(df_raw)
    df_ptse = PTSE(df_raw)
    df = ptse2normalizedLR(df_ptse)
    df_raw_melted = meltTissues_df(df_raw, "meanTPM")
    df_ptse_melted = meltTissues_df(df_ptse, "ptse")
    df_scaledLR_melted = meltTissues_df(df, "normalizedLR")
    df_melted = pd.concat([df_raw_melted, df_ptse_melted, df_scaledLR_melted], axis=1)
    df_melted = df_melted.loc[:, ~df_melted.columns.duplicated()]
    #df_exTSP = df_vars_mp2.loc[:,["Gene", "Transcript_id", "Variant", "Tissue", "exTSP", "ptse", "IsoPath"]]
    df_melted.to_csv(PTSE_file, index=False, sep='\t')
    print(f"PTSE scores saved in {PTSE_file}")

if __name__ == "__main__":
    if Path(PTSE_file).exists():
        print(f"PTSE scores are already exist in {PTSE_file}")
        pass
    main()
       

