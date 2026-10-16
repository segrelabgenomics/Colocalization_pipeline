#!/usr/bin/env python3

import sys
import os
import os.path
import argparse
import subprocess
import re
import pandas as pd
import numpy as np

from typing import Dict, List, Tuple
from pandas import DataFrame
from argparse import ArgumentParser, Namespace

def parse_args() -> Namespace:
    """
    arguments
    """
    parser = ArgumentParser(description="Arguments to pass when running script")

    parser.add_argument("--ecaviar_table",
                        type=str, help="ecaviar table")
    parser.add_argument("--tissue",
                        type=str, help="tissue")
    parser.add_argument("--qtl_type",
                        type=str, help="qtl")
    parser.add_argument("--gene",
                        type=str, help="gene")
    parser.add_argument("--gwas_trait",
                        type=str, help="gwas trait")
    parser.add_argument("--output_directory",
                        type=str, help="output directory")

    return parser.parse_args()

def constants() -> Dict[str,str]:
    """
    dictionary of pathways
    """
    basedir = "../data/GTEx/V8_release/eQTLs/"
    constants_dict = {"sig_eqtl_dir": os.path.join(basedir, "GTEx_Analysis_v8_eQTL_all_associations/"),
                      "eqtl_file_suffix": f".allpairs.txt.gz",
                      "sig_sqtl_dir": os.path.join(basedir, "GTEx_Analysis_v8_sQTL_all_associations/"),
                      "sqtl_file_suffix": f".v8.sqtl_allpairs.txt.gz"}
    return constants_dict

def parse_ecaviar_table(ecaviar_table_file: str,
                        qtl_type: str,
                        tissue: str,
                        gene: str) -> DataFrame:
    """
    subsets table by gene-qtl-clpp cutoff
    """
    df = pd.read_csv(ecaviar_table_file, sep="\t")
    df = df.loc[(df["qtl"] == qtl_type) &
                (df["CLPP"] >= 0.01) &
                (df["tissue"] == tissue) &
                (df["gene_id"] == gene)].copy()
    return df

def retrieve_chrom(df: DataFrame) -> str:
    """
    returns chromosome
    """
    df = df[["loci"]].copy()
    df["chr"] = df["loci"].apply(lambda x: x.split(":")[0])
    if len(df["chr"].unique()) == 1:
        return df["chr"].unique()[0]
    else:
        raise Exception("Chromosome is not unique.")
        sys.exit(1)

def retrieve_tissue_file(tissue: str,
                         qtl: str) -> str:
    """
    returns tissue file relevant to qtl and chromosome
    """
    constants_dict = constants()
    if qtl == "eQTL":
        print('checking the gtex eqtl path:')
        return f"{constants_dict['sig_eqtl_dir']+tissue+constants_dict['eqtl_file_suffix']}"
    if qtl == "sQTL":
        print('checking the gtex sqtl path:')
        return f"{constants_dict['sig_sqtl_dir']+tissue+constants_dict['sqtl_file_suffix']}"

def retrieve_genes(df: DataFrame) -> List[str]:
    """
    returns list of intron column
    """
    return list(set(df["gene_id"]))

def fix_snpid(df: DataFrame) -> DataFrame:
    """
    chrN_pos_ref_alt_b38 -> N:pos:ref:alt
    """
    df["SNP_ID"] = df["variant_id"]
    df["SNP_ID"] = df["SNP_ID"].str.replace("_b38", "")
    df["SNP_ID"] = df["SNP_ID"].str.replace("chr", "")
    df["SNP_ID"] = df["SNP_ID"].str.replace("_", ":")
    return df

def parse_gtex_allpairs_file(retina_tissue_file: str,
                             gene: str,
                             chrom: str) -> DataFrame:
    """
    zgrep's intron from GTEx tissue 
    """
    command = ["zgrep", gene, retina_tissue_file]
    grep_output = subprocess.Popen(command,stdout=subprocess.PIPE)
    header = subprocess.check_output(["bash", "-c", f"zcat {retina_tissue_file} | head -n 1"]).decode().strip().split("\t")
    retina_tissue = pd.read_csv(grep_output.stdout,sep="\t", header=None,names=header)
    print(retina_tissue.head())
                              
    retina_tissue = retina_tissue[["gene_id", "variant_id",
                                   "pval_nominal","slope", "slope_se"]]
    retina_tissue = fix_snpid(retina_tissue)

    return retina_tissue

def extract_genes_tissue_ecaviar(ecaviar_table_retina_tissue: DataFrame,
                                 retina_tissue_file: str,
                                 chrom: str) -> DataFrame:
    """
    extracts genes from gtex
    """
    genes = retrieve_genes(ecaviar_table_retina_tissue)

    full_df = pd.DataFrame()

    for gene in genes:
        retina_tissue = parse_gtex_allpairs_file(retina_tissue_file, gene, chrom)
        tmp_df = ecaviar_table_retina_tissue.loc[ecaviar_table_retina_tissue["gene_id"] == gene]
        df = pd.merge(tmp_df, retina_tissue, on=["gene_id", "SNP_ID"])
        full_df = pd.concat([full_df, df], ignore_index=True)
    return full_df

def main():
    """
    main function
    """
    args = parse_args()
    # gwas trait
    gwas_trait = args.gwas_trait

    # filter table by CLPP
    print("parsing table")
    ecaviar_table_df = parse_ecaviar_table(args.ecaviar_table, 
                                           args.qtl_type, args.tissue, args.gene)
    #retrieve chromosome
    print("Chromosome")
    print(ecaviar_table_df.head())
    chrom = retrieve_chrom(ecaviar_table_df)
    print(chrom)
     # # # tissue file 
    tissue_file = retrieve_tissue_file(args.tissue, args.qtl_type)
    # 
    print(args.tissue)
    print(args.qtl_type)
    print(args.gene)
    print(ecaviar_table_df) 
    print("Tissue file")
    print(tissue_file)
    if not ecaviar_table_df.empty:
        print("extracting slope and se")
        ecaviar_table_tissue_completed = extract_genes_tissue_ecaviar(ecaviar_table_df,
                                                                      tissue_file,
                                                                      chrom)
        filename = f"{args.output_directory}{gwas_trait}_{args.tissue}_{args.qtl_type}_{args.gene}_full_qtl_info_ecaviar_table.tsv"
        ecaviar_table_tissue_completed.to_csv(filename,sep="\t", header=True, index=None)
    else:
        filename = f"{args.output_directory}{gwas_trait}_{args.tissue}_{args.qtl_type}_{args.gene}_empty_table.tsv"
        ecaviar_table_retina_tissue.to_csv(filename,sep="\t", header=True, index=None)

if __name__ == "__main__":
    main()

