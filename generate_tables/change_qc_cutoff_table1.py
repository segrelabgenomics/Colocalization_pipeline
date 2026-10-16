#!/usr/bin/env python3

import pandas as pd
from pandas import DataFrame

def add_filter_column(df: DataFrame) -> DataFrame:
    """
    column if passes filtering
    """
    df["QC_name_QTL_FDR05_P1E04_GWAS_P1E05"] = False

    df.loc[(df["QTL_gwas_p_value"] <= 0.00001) &
           ((df["QTL_variant_signif_esQTL"] == True) |  
               (df["QTL_pval_nominal"] <= 0.0001)),
           "QC_name_QTL_FDR05_P1E04_GWAS_P1E05"] = True

    return df

file = "POAG_EUR_table1_formatted_with_QC.tsv"
df = pd.read_csv(file, sep="\t")
print(df.head())
df = add_filter_column(df)
outname = "POAG_EUR_table1_formatted_with_QC_UPD.tsv"

df.to_csv(outname, sep="\t", header=True, index=None)




