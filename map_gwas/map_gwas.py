#!/usr/bin/env python3

import os
import os.path
import re
import argparse
import pandas as pd
import numpy as np

def check_chromsome_column(df):
    """
    if chr header is not format, chrN, convert to format
    """
    df["chr"] = df["chr"].astype(str)
    chr_list = list(df["chr"].unique())
    improper_format_chr = [chr for chr in chr_list if not chr.startswith("chr")]
    if improper_format_chr:
        df["chr"] = "chr"+df["chr"].astype(str)
        return(df)
    return(df)

def parse_gwas(gwas_file):
    """
    reads gwas
    position is int
    compute zscore
    """
    gwas = pd.read_csv(gwas_file,sep=" ")
    print(gwas.head())
    gwas["chr"] = gwas[args.chr_name].astype(str).apply(lambda x: x if x.startswith("chr") else "chr" + x)
    #gwas["pos"] = gwas[args.bp_name].astype(int)
    gwas["pos"] = pd.to_numeric(gwas[args.bp_name], errors='coerce').astype("Int64")
    gwas = check_chromsome_column(gwas)
    gwas['effect_allele'] =gwas[args.a1].str.upper()
    gwas['Other_allele'] = gwas[args.a2].str.upper()
    gwas["zscore"] = gwas[args.beta]/gwas[args.slope]
    return(gwas)

def parse_variants(variants_file,var_col):
    """
    reads file
    separates into chrom, pos, ref, alt
    """
    variants = pd.read_csv(variants_file)
    sep = ":" if args.ref == "1KG" else "_"
    variants["variant"] = variants[var_col]
    split_cols = variants["variant"].str.split(sep, expand=True)
    if args.ref == "1KG":
        variants["chr"]  = split_cols[0]
        variants["pos"]  = split_cols[1]
        variants["ref"]  = split_cols[2]
        variants["alt"]  = split_cols[3]
    else:
        variants["chr"]   = split_cols[0]
        variants["pos"]   = split_cols[1]
        variants["ref"]   = split_cols[2]
        variants["alt"]   = split_cols[3]
    
    variants["pos"] = variants["pos"].astype(int)
    return(variants)

def align_gwas(gwas,variants):
    """
    merge two dfs
    align allele
    """
    trait_aligned = pd.merge(variants,gwas,on=["chr","pos"])
    print(trait_aligned.head())
    conditions  = [ ((trait_aligned["alt"] == trait_aligned["effect_allele"]) & (trait_aligned["ref"] == trait_aligned["Other_allele"])),((trait_aligned["alt"] == trait_aligned["Other_allele"]) & (trait_aligned["ref"] == trait_aligned["effect_allele"]))]
    choices = [ 1,-1 ]

    trait_aligned["record_change"] = np.select(conditions, choices, default=0)

    conditions_effect = [(trait_aligned["record_change"] == 1),
                             (trait_aligned["record_change"] == -1),
                             (trait_aligned["record_change"] == 0)]
    choices_effect =  [trait_aligned["zscore"],(-1)*trait_aligned["zscore"],trait_aligned["zscore"]]
    trait_aligned["new_zscore"] = np.select(conditions_effect, choices_effect)
    trait_aligned = trait_aligned.loc[trait_aligned["record_change"] != 0]
    trait_aligned["build"] = args.build
    trait_aligned = trait_aligned.drop(["ref","alt","build","zscore","record_change"],axis=1)
    trait_aligned = trait_aligned.rename(columns={"new_zscore":"zscore"})
    return(trait_aligned)

def parse_args():
    parser = argparse.ArgumentParser(description="Match GWAS with variant file and output results.")

    parser.add_argument("--gwas_file", required=True, help="Path to the GWAS input file")
    parser.add_argument("--variants_file", required=True, help="Reference Panel, use 1KG or GTEx")
    parser.add_argument("--ref", required=True, help="Path to the file containing variants. 1KG if the variants are ':' seperated and '_' otherwise")
    parser.add_argument("--output_name", required=True, help="Path for the output file")
    parser.add_argument("--chr_name", required=True, help="Column name for chromosome in GWAS file without chr prefix")
    parser.add_argument("--bp_name", required=True, help="Column name for base pair position in GWAS file")
    parser.add_argument("--build", required=True, help="Build 38")
    parser.add_argument("--a1", required=True, help="Column name for effect allele")
    parser.add_argument("--a2", required=True, help="Column name for other allele")
    parser.add_argument("--beta", required=True, help="Column name for beta/effect size")
    parser.add_argument("--slope", required=True, help="Column name for standard error")
    parser.add_argument("--variant_col_in_reference", required=True, help="Column name for the variant identifier in the reference panel")
    return parser.parse_args()


if __name__ == "__main__":


    args = parse_args()

    print("parsing gwas")
    gwas = parse_gwas(args.gwas_file)

    print("parsing variants")
    variants = parse_variants(args.variants_file,args.variant_col_in_reference)
        
    print("aligning gwas")
    aligned_gwas = align_gwas(gwas,variants)
    aligned_gwas.to_csv(args.output_name,
                        sep="\t",index=None,header=True)
