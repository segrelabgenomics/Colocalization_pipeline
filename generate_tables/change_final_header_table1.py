#!/usr/bin/env python3

import pandas as pd
from pandas import DataFrame

def update_column_header(df: DataFrame) -> DataFrame:
    """
    Rename header August, 2024
    """
    
    df.columns = df.columns.str.replace(" ","_")
    df = df.rename(columns={'GWAS_variant_(aligned_to_1KG)': "GWAS_variant_aligned_to_1KGP",
                             'Effect_allele':'GWAS_effect_allele', 
                             'Other_Allele':'GWAS_other_allele', 
                             'rsID_GWAS':'GWAS_rsID', 
                             'GWAS_OR':'GWAS_Odds_Ratio', 
                             'qtl_Variant_id':'QTL_variant_ID', 
                             'rsID_qtl':'QTL_rsID', 
                             'gene_id':'Gene_ENS_ID',
                             'Gene_id_intron_cluster':'Gene_ID_intron_cluster', 
                             'QTL_gwas_p_value':'QTL_GWAS_p_value', 
                             'STDERR' : 'GWAS_se',
                             'GWAS_beta_relative_ALT':'GWAS_effect_relative_ALT', 
                             'Qtl_gwas_relative_effect':'QTL_GWAS_relative_direction_effect',
                             'QTL_variant_signif_esQTL':'QTL_variant_signif_esQTL_FDR05',
                             'QC_name_QTL_FDR05_P1E04_GWAS_P5E02':'Pass_QC_QTL_FDR05_P1E04_GWAS_P1E02',
                             'QC_name_QTL_FDR05_P1E04_GWAS_P1E05':'Pass_QC_QTL_FDR05_P1E04_GWAS_P1E05'})
    return df

def add_column(df: DataFrame) -> DataFrame:
    """
    Adding a missing column from AH pipeline 
    """
    df["QTL_variant_ALT_(effect_allele)"] = df['QTL_variant_ID'].str.split('_').str[3]
    return df

def order_columns_df(df: DataFrame) -> DataFrame:
    """
    Ordering the dataframe on desired column names
    """
    print(df.columns)
    desired_order = ["GWAS_Trait","GWAS_variant",
                    "GWAS_variant_aligned_to_1KGP",
                    "GWAS_effect_allele","GWAS_other_allele",
                    "GWAS_rsID","GWAS_Odds_Ratio","GWAS_effect",
                    "ABS(GWAS_effect)","GWAS_p_value",
                    "LD_interval_start_(b38)","LD_interval_end_(b38)",
                    "LD_interval_size_(kb)",
                    "Number_variants_tested_in_defined_LD_interval_(r2_0.1_±_50kb)",
                    "QTL_Type","QTL_variant_ID","QTL_variant_ALT_(effect_allele)",
                    "QTL_rsID","Gene_ENS_ID",
                    "Gene_symbol","Gene_ID_intron_cluster",
                    "Tissue","QTL_pval_nominal","QTL_GWAS_p_value",
                    "QTL_slope_relative_ALT","QTL_se","GWAS_effect_relative_ALT",
                    "QTL_GWAS_relative_direction_effect",
                    "eCAVIAR_prob_in_causal_set","CLPP",
                    "QTL_variant_signif_esQTL_FDR05",
                    "Pass_QC_QTL_FDR05_P1E04_GWAS_P1E02",
                    "Pass_QC_QTL_FDR05_P1E04_GWAS_P1E05"] 
    df = df[desired_order]
    return df

def parse_df(df: DataFrame) -> DataFrame:
    """
    Parse Table1
    """
    df = update_column_header(df)
    print("Header updated")
    df = add_column(df)
    print("Missing Effect Allele Column added")
    df = order_columns_df(df)
    print("Columns orderd")
    
    return df

file = "POAG_EUR_table1_formatted_with_QC_UPD.tsv"
df = pd.read_csv(file, sep="\t")
df = parse_df(df)
outname = "POAG_EUR_table1_formatted_with_QC_UPD_renamed.tsv"
df.to_csv(outname, sep="\t", header=True, index=None)

print("Done")



