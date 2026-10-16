#!/usr/bin/env python3

import os
import os.path
import sys
import pandas as pd
import numpy as np
import re
import argparse

from typing import List
from argparse import ArgumentParser, Namespace

def parse_args():
    """
    arguments to pass
    """
    parser = ArgumentParser(description="Arguments to pass when parsing through eCAVIAR results")

    parser.add_argument("--ecaviar_directory",
                        type=str, help="Directory consisting of eCAVIAR outputs")
    parser.add_argument("--gwas_trait",
                        type=str, help="GWAS trait of interest, e.g. poag_transancestry")
    parser.add_argument("--output_file",
                        type=str, help="Name of output file")
    parser.add_argument("--incorrect_num_file",
                        type=str, help="Name of output file consisting of loci with incorrect matching of col and z-score files")

    return(parser.parse_args())

def split_file(file: str) -> List:
    """
    splits file into list
    """
    return(file.split("_"))

def retrieve_second_variant(file: str,
                            qtl: str) -> str:
    """
    second variant
    4 for eqtl
    5 for sqtl
    """
    split = split_file(file)
    print(split)
    if qtl.lower() == "eqtl":
        return(split[4])# 4 in normal case
    return(split[5])# 5 in normal case, because leafcutter was updated

def retrieve_gene(file: str,
                  qtl: str) -> str:
    """
    gene name
    index 2 for eqtl
    2:4 for sqtl
    """
    split = split_file(file)
    print(file)
    if qtl.lower() == "eqtl":
        return(split[2]) #4
    return("_".join(split[2:4])) #[2:4]

def retrieve_tissue(file: str,
                    qtl: str) -> str:
    """
    tissue name
    5:-2 for eqtl
    6:-2 for sqtl
    """
    split = split_file(file)
    if qtl.lower() == "eqtl":
        return("_".join(split[5:-2]))
    return("_".join(split[6:-2]))

def retrieve_gwas_trait(file: str) -> str:
    """
    0:2 index
    """
    split = split_file(file)
    return("_".join(split[0:2]))#0:2

def retrieve_qtl(sub_directory: str) -> List:
    """
    4th index
    """
    return(sub_directory.split("_")[3])

def retrieve_loci(sub_directory: str) -> List:
    """
    2nd index
    """
    return(sub_directory.split("_")[2])

def main():

    args = parse_args()
    sub_directories = [sub for sub in os.listdir(args.ecaviar_directory) if args.gwas_trait in sub]

    full_df = pd.DataFrame()
    incorrect_num_df = pd.DataFrame()

    for sub_directory in sub_directories:
        new_directory = args.ecaviar_directory+sub_directory
        new_files = [file for file in os.listdir(new_directory) if file.endswith("_col")]

        if not new_files:
            continue

        for new_file in new_files:
            qtl = retrieve_qtl(sub_directory)
            df = pd.read_csv(new_directory+"/"+new_file,sep="\t")
            df["loci"] = retrieve_loci(sub_directory)
            df["qtl"] = qtl
            df["gwas_variant_aligned"] = retrieve_second_variant(new_file,qtl)
            df["tissue"] = retrieve_tissue(new_file,qtl)
            df["gene"] = retrieve_gene(new_file,qtl)
            df["trait"] = retrieve_gwas_trait(new_file)
            df["num_variants_tested"] = len(df)
            # z-score df
            z_score_df_file = re.sub("_eCAVIAR_col", "_qtl_z_score.txt", new_file)
            try:
                z_score_df = pd.read_csv(new_directory+"/"+z_score_df_file, sep="\t", header=None)
            except:
                #print(z_score_df_file, new_file)
                with open(args.gwas_trait + "_something_strange.txt", "a") as f:
                    f.write(new_file+"\n")
                continue
            df["num_variants_z_score_file"] = len(z_score_df)
            
            df_sig_variant = df.loc[df["CLPP"] >= 0.01].copy()
            #print(df_sig_variant)
            if not df_sig_variant.empty:
                if len(z_score_df) != len(df):
                    print(len(z_score_df))
                    print(len(df))
                    incorrect_num_df = pd.concat(
                        [incorrect_num_df, df_sig_variant],
                        ignore_index=True
                    )
                else:
                    full_df = pd.concat(
                        [full_df, df_sig_variant],
                        ignore_index=True
                    )
            else:
                df_max_variant = df.loc[df["CLPP"].idxmax()].to_frame().T

                if len(z_score_df) != len(df):
                    incorrect_num_df = pd.concat(
                        [incorrect_num_df, df_max_variant],
                        ignore_index=True
                    )
                else:
                    full_df = pd.concat(
                        [full_df, df_max_variant],
                        ignore_index=True
                    )

    full_df["gene_id"] = full_df["gene"]
    full_df["gene"] = full_df["gene_id"].apply(lambda x: x.split(":")[-1])

    full_df = full_df.sort_values(["CLPP"], ascending=False)
    print("writing output")
    full_df.to_csv(args.output_file,sep="\t",header=True,index=None)

    try:
        incorrect_num_df["gene_id"] = incorrect_num_df["gene"]
        incorrect_num_df["gene"] = incorrect_num_df["gene_id"].apply(lambda x: x.split(":")[-1])
        incorrect_num_df = incorrect_num_df.sort_values(["CLPP"], ascending=False)
        incorrect_num_df.to_csv(args.incorrect_num_file, sep="\t", header=True, index=None)
    except:
        pass

if __name__ == "__main__":
    main()