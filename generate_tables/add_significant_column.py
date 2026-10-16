#!/usr/bin/env python3

import os
import os.path
import sys
import re
import pandas as pd

from typing import List, Dict, Tuple
from pandas import DataFrame

def retrieve_constants_gtex() -> Dict[str, str]:
    """
    constants for significant variant gene pair files
    from GTEx
    """
    basedir = "../data/GTEx/V8_release/eQTLs/"
    constants_dict = {"sig_eqtl_dir": os.path.join(basedir, "GTEx_Analysis_v8_eQTL_all_associations/"),
                      "eqtl_file_suffix": f".allpairs.txt.gz",
                      "sig_sqtl_dir": os.path.join(basedir, "GTEx_Analysis_v8_sQTL_all_associations/"),
                      "sqtl_file_suffix": f".v8.sqtl_allpairs.txt.gz"}
    return constants_dict


def retrieve_tissue_qtl(constants_dict, qtl_type):
    """
    returns tissue_qtl
    key in form: sig_tissue_qtl_dir
    """
    list_qtl_tissue=[]
    if qtl_type.lower() == "eqtl":
       for file_qtl in os.listdir(constants_dict["sig_eqtl_dir"]):
          if file_qtl.endswith(constants_dict["eqtl_file_suffix"]):
              temp_tissue = file_qtl.strip().split(".")[0]
              list_qtl_tissue.append((temp_tissue, "eQTL"))
    if qtl_type.lower() == "sqtl":
       for file_qtl in os.listdir(constants_dict["sig_sqtl_dir"]):
          if file_qtl.endswith(constants_dict["sqtl_file_suffix"]):
              temp_tissue = file_qtl.strip().split(".")[0]
              list_qtl_tissue.append((temp_tissue, "sQTL"))
    print('Found! Returning to Main function')
    return list_qtl_tissue


def retrieve_files(constants_dict, tissue,qtl_type):
    """
    returns GTEx files from directory given qtl type
    """
    if qtl_type.lower() == "eqtl":
        target = tissue + constants_dict["eqtl_file_suffix"]
        files = os.listdir(constants_dict["sig_eqtl_dir"])
        return [file for file in files if file == target]

    if qtl_type.lower() == "sqtl":
        target = tissue + constants_dict["sqtl_file_suffix"]
        files = os.listdir(constants_dict["sig_sqtl_dir"])
        return [file for file in files if file == target]

def parse_gtex_df(constants_dict: Dict[str, str],
                  list_file_sig: list,
                  tissue: str,
                  qtl_type: str) -> List[str]:
    """
    reads gtex file
    selects relevant columns
    adds tissue header
    """
    list_parse_gtex_df=[]
    if qtl_type.lower() == "eqtl":
       for file_sig in list_file_sig:
          df = pd.read_csv(constants_dict["sig_eqtl_dir"] + file_sig, sep="\t")
          df["id"] = df["variant_id"]+"_"+df["gene_id"]+"_"+tissue
          list_parse_gtex_df.extend(set(df["id"]))
    if qtl_type.lower() == "sqtl":
       for file_sig in list_file_sig:
          df = pd.read_csv(constants_dict["sig_sqtl_dir"] + file_sig, sep="\t")
          df["id"] = df["variant_id"]+"_"+df["gene_id"]+"_"+tissue
          list_parse_gtex_df.extend(set(df["id"]))
    print('checking list_parse_gtex_df before returning to main')
    return list_parse_gtex_df


def add_filter_column(df: DataFrame) -> DataFrame:
    """
    column if passes filtering
    """
    df["QC_name_QTL_FDR05_P1E04_GWAS_P5E02"] = False
    print(df)
    df.loc[(df["QTL_gwas_p_value"] <= 0.05) & ((df["QTL_variant_signif_esQTL"] == True) | (df["QTL_pval_nominal"] <= 0.0001)), "QC_name_QTL_FDR05_P1E04_GWAS_P5E02"] = True
    return df

def parse_cutoff_output(cutoff_output: str,
                        constants_dict: Dict[str, str]) -> DataFrame:
    """
    loops through gtex tissue and determines if in outputs
    """
    print('Start Parsing..')
    cutoff_df = pd.read_csv(cutoff_output, sep="\t")
    cutoff_df_new = pd.DataFrame()
    list_esqtl_tissue=[]
    for qtl_type in ['eqtl','sqtl']:
        if qtl_type == 'eqtl':
           list_eqtl_tissue = retrieve_tissue_qtl(constants_dict, qtl_type)
           print(list_eqtl_tissue)
           next
        if qtl_type == 'sqtl':
           list_sqtl_tissue = retrieve_tissue_qtl(constants_dict, qtl_type)
           print(list_sqtl_tissue)

    for element in list_sqtl_tissue:
       list_eqtl_tissue.append(element)
    list_esqtl_tissue=list_eqtl_tissue
    for tissue, qtl in list_esqtl_tissue:
        print(tissue)
        print(qtl)
        df = cutoff_df.loc[(cutoff_df["Tissue"] == tissue) & (cutoff_df["QTL Type"] == qtl)].copy()
        if df.empty:
            print('df is empty')
            next
        print('df is not empty')
        df["id"] = df["qtl_Variant_id"]+"_"+df["gene_id"]+"_"+tissue
        sig_file = retrieve_files(constants_dict,tissue,qtl)
        sig_sets = parse_gtex_df(constants_dict, sig_file, tissue, qtl)
        df["QTL_variant_signif_esQTL"] = df["id"].apply(lambda x: True if x in sig_sets else False)
        cutoff_df_new = pd.concat([cutoff_df_new, df], ignore_index=True)

    cutoff_df_new = cutoff_df_new.drop(["id"], axis=1)
    cutoff_df_new_filtered = add_filter_column(cutoff_df_new)

    return cutoff_df_new


def main():
    """
    main function
    """
    
    cutoff_01_file = sys.argv[1]
    output_dir = sys.argv[2]
    constants_dict = retrieve_constants_gtex()
    
    cutoff_df_added_col = parse_cutoff_output(cutoff_01_file, constants_dict)
    cutoff_filename = output_dir + "POAG_EUR_table1_formatted_with_QC.tsv"
    cutoff_df_added_col.to_csv(cutoff_filename, sep="\t", header=True, index=None, na_rep="NA")

if __name__ == "__main__":
    main()
