

import pandas as pd
import numpy as np
from scipy.stats import zscore
from scipy.stats import norm
from ExTSP.config import GTeX_file, BrainSpan_expression_file, PTSE_file_gtex, PTSE_file_brainspan
from pathlib import Path




#PTSE_file = f'{PROCESSED_DIR}/GTEx_PTSE.csv'


def ptse2normalizedLR(ptse, eps=1e-10):
    """
    ptse: 2D array (n_rows, n_cols), row-wise probabilities in (0, 1] (e.g. after PTSE).

    lr = ptse / (1 - ptse), then each row is divided by its max LR.
    Rows where max LR is 0 are left unchanged (or use lr_max=1).

    Returns array same shape as ptse.
    """
    ptse = np.asarray(ptse, dtype=float)
    if ptse.ndim != 2:
        raise ValueError("ptse2normalizedLR expects a 2D array")

    # Avoid ptse == 1 -> 0 denominator; clip upper tail slightly if needed
    denom = np.maximum(1.0 - ptse, eps)
    lr = ptse / denom

    lr_max = np.max(lr, axis=1, keepdims=True)
    assert any(lr_max <= 0), "lr_max is less than or equal to 0"
    lr_max = np.where(lr_max == 0, 1.0, lr_max)
    normalized_lr = lr / lr_max
    return normalized_lr

def zscore_rows(x, eps=1e-10):
    """Per-row z-score across columns. x shape (n_rows, n_cols)."""
    x = np.asarray(x, dtype=float)
    mu = x.mean(axis=1, keepdims=True)
    std = x.std(axis=1, keepdims=True)
    std = np.where(std < eps, 1.0, std)
    z = (x - mu) / std
    z = np.nan_to_num(z, nan=0.0, posinf=0.0, neginf=0.0)
    return z

def PTSE(expr):
    """
    x: 2D array, shape (n_rows, n_cols). All columns are treated as features (e.g. tissues).
    Returns: same shape, row-wise softmax after row-wise z-score; values floored at 1e-10.
    """
    expr = np.asarray(expr, dtype=float)
    if expr.ndim != 2:
        raise ValueError("PTSE_array expects a 2D array")

    z = zscore_rows(expr)

    x_max = np.max(z, axis=1, keepdims=True)
    exps = np.exp(z - x_max)
    ptse = exps / np.sum(exps, axis=1, keepdims=True)
    ptse = np.maximum(ptse, 1e-10)
    return ptse

# def PTSE(df_input):
#     """
#     Copy the dataframe, apply row-wise softmax to the first 52 columns,
#     then divide those columns by their row-wise maximum. Returns the
#     modified copy.
#     """
#     df = df_input.copy()
#     df = zScoreNormalize(df)
#     # Select first 52 columns as a NumPy array
#     x = df.iloc[:, 0:52].to_numpy(dtype=float)
#     # Row-wise softmax: subtract row max for numerical stability
#     x_max = np.max(x, axis=1, keepdims=True)
#     exps = np.exp(x - x_max)
#     # compute ptse
#     ptse = exps / np.sum(exps, axis=1, keepdims=True)
#     ptse = np.maximum(ptse, 1e-10)
#     df.iloc[:, 0:52] = ptse
#     return df

# def ptse2normalizedLR(df_ptse):
#     df = df_ptse.copy()
#     ptse = df.iloc[:, 0:52].to_numpy(dtype=float)
#      # compute lr
#     lr = ptse / (1 - ptse)
#     # Divide by row-wise max of the softmax values
#     lr_max = np.max(lr, axis=1, keepdims=True)
#     # Avoid division by zero
#     lr_max[lr_max == 0] = 1.0
#     # scale by row-wise max of the LR
#     normalized_lr = lr / lr_max
#     # Put back into dataframe
#     df.iloc[:, 0:52] = normalized_lr
#     return df

# def zScoreNormalize(df_input,robust=False):
#     df = df_input.copy()
#     if robust:
#         q1 = df.iloc[:,0:52].quantile(0.25, axis=1)
#         q3 = df.iloc[:,0:52].quantile(0.75, axis=1)
#         iqr = q3 - q1
#         robust_std = iqr / 1.349
#         robust_std[robust_std==0] = 10**-6
#         mu = df.iloc[:,0:52].median(axis=1)
#         std = robust_std
#     else:
#         std =  df.iloc[:,0:52].std(axis=1)
#     mu = df.iloc[:,0:52].mean(axis=1)
#     df.iloc[:,0:52] = df.iloc[:,0:52].sub(mu, axis=0).div(std, axis=0)
#     df.iloc[:,0:52] = df.iloc[:,0:52].fillna(0)
#     return df

# def zScoreNormalizeSignificant(df_input,robust=False):
#     df = zScoreNormalize(df_input,robust)
#     df.iloc[:,0:52] = 1-norm.cdf(df.iloc[:,0:52])
#     df.iloc[:,0:52] = (df.iloc[:,0:52] <=0.001).astype(float)
#     df["numSignificant"] =  df.iloc[:,0:52].sum(axis=1)
#     return df

def sigTissueCounts(df):
    nTissuesSig = df.iloc[:,0:52].apply(sum, axis=1)
    print("number of transcripts in GTEx:", df.shape[0])
    print("number of transcripts in GTEx where at aleast one tissue has significant expression:", (nTissuesSig>0).sum())
    print(nTissuesSig.value_counts())

def meltTissues_df(df, new_col_name, idx, idx_meta):
     # First 52 columns are tissues; remaining columns are variant/isoform/gene/etc.
    tissue_cols = df.columns[idx]
    meta_cols = df.columns[idx_meta]

    # Long format for PTSE values
    df_long = df.melt(
        id_vars=meta_cols,
        value_vars=tissue_cols,
        var_name="tissue",
        value_name=new_col_name
    )
    return df_long

def main(df_expr, idx, idx_meta):
    all_zero_rows = df_expr.iloc[:,idx].sum(axis=1)==0
    idx_zero = df_expr.iloc[:,idx] == 0
    print(f"number of rows with all zeros: {all_zero_rows.sum()}")
    #df_expr = df_expr[~all_zero_rows]
    expr = df_expr.iloc[:,idx].to_numpy(dtype=float)
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
    ptse = PTSE(expr)
    #ptse[all_zero_rows] = 0.0
    ptse[idx_zero] = 0.0
    df_ptse = df_expr.copy()
    df_ptse.iloc[:,idx] = ptse
    normalizedLR = ptse2normalizedLR(ptse)
    #normalizedLR[all_zero_rows] = 0.0
    normalizedLR[idx_zero] = 0.0
    df_normalizedLR = df_expr.copy()
    df_normalizedLR.iloc[:,idx] = normalizedLR 
    df_raw_melted = meltTissues_df(df_expr, "meanTPM", idx, idx_meta)
    df_ptse_melted = meltTissues_df(df_ptse, "ptse", idx, idx_meta)
    df_scaledLR_melted = meltTissues_df(df_normalizedLR, "normalizedLR", idx, idx_meta)
    df_melted = pd.concat([df_raw_melted, df_ptse_melted, df_scaledLR_melted], axis=1)
    df_melted = df_melted.loc[:, ~df_melted.columns.duplicated()]
    return df_melted

if __name__ == "__main__":
    #if Path(PTSE_file_gtex).exists():
    if False:
        print(f"PTSE scores are already exist in {PTSE_file_gtex}")
    else:
        df_expr = pd.read_csv(GTeX_file, sep='\t')
        df_expr = df_expr.fillna(0.0)
        idx = slice(0, 52)
        idx_meta = slice(52, df_expr.shape[1])
        df_melted = main(df_expr, idx=idx, idx_meta=idx_meta)
        df_melted.to_csv(PTSE_file_gtex, index=False, sep='\t')
        print(f"PTSE scores saved in {PTSE_file_gtex}")
    if Path(PTSE_file_brainspan).exists():
        print(f"PTSE scores are already exist in {PTSE_file_brainspan}")
    else:
        df_expr = pd.read_csv(BrainSpan_expression_file, sep='\t')
        df_expr = df_expr.fillna(0.0)
        idx = slice(3,15)
        idx_meta = list(range(0, 3)) + list(range(15, df_expr.shape[1]))
        df_melted = main(df_expr, idx=idx, idx_meta=idx_meta)
        df_melted.to_csv(PTSE_file_brainspan, index=False, sep='\t')
        print(f"PTSE scores saved in {PTSE_file_brainspan}")
    print("PTSE scores computation completed")
   
    
       

