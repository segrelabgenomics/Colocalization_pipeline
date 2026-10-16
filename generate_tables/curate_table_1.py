#!/usr/bin/env python3

import os
import os.path
import re
import subprocess
import argparse
import pandas as pd
import numpy as np
from typing import Dict, List, Tuple
from pandas import DataFrame
from argparse import ArgumentParser, Namespace

def parse_args() -> Namespace:
    """
    arguments to run script
    """
    parser = ArgumentParser(description="parameters for generating table 1")

    parser.add_argument("--ecaviar_output_table",
                        type=str, help="ecaviar outputs", required=True)
    parser.add_argument("--gwas_file",
                        type=str, help="gwas file", required=True)
    parser.add_argument("--output_filename",
                        type=str, help="name of output file", required=True)

    return parser.parse_args()

def create_id_column(df: DataFrame) -> DataFrame:
    """
    changes chrN_pos_ref_alt_build to N:pos:ref:alt
    """
    df["SNP_ID"] = df["variant_id"]
    df["SNP_ID"] = df["SNP_ID"].str.replace("chr", "")
    return df

def rename_column(df: DataFrame) -> DataFrame:
    """
    changes SNP to MARKERNAME
    """
    df = df.rename(columns = {"Marker": "MARKERNAME",
                    "variant":"variant_id",
                    "P-value":"gwas_p_value",
                    "mtag_beta":"BETA",
                    "Effect":"BETA"})

    return df

def fix_effect(df: DataFrame) -> DataFrame:
    """
    ensures z-score and effect are same sign
    """
    conditions  = [ ((df["zscore"] < 0) & (df["BETA"] < 0)),
                    ((df["zscore"] < 0) & (df["BETA"] > 0)),
		    ((df["zscore"] > 0) & (df["BETA"] < 0)),
                    ((df["zscore"] > 0) & (df["BETA"] > 0))]

    choices = [ 1,-1,-1,1]
    df["record_change"] = np.select(conditions, choices, default=0)

    conditions_effect = [(df["record_change"] == 1),
                         (df["record_change"] == -1)]

    choices_effect =  [df["BETA"],(-1)*df["BETA"]]

    df["BETA"] = np.select(conditions_effect, choices_effect)
    df = df.drop(["record_change"], axis=1)
    
    return df

def parse_gwas(gwas_file: str) -> DataFrame:
    """
    fixes effect
    """
    df = pd.read_csv(gwas_file, sep="\t")
    df = rename_column(df)
    df = fix_effect(df)

    #df = create_id_column(df)

    return df

def merge_gwas(ecaviar_output_table: str,
               gwas_df: DataFrame) -> DataFrame:
    """
    merges gwas with ecaviar output
    """
    ecaviar_df = pd.read_csv(ecaviar_output_table, sep="\t")
    
    ecaviar_df["chr"] = ecaviar_df["variant_id"].apply(lambda x: x.split("_")[0]).astype(str)
    ecaviar_df["pos"] = ecaviar_df["variant_id"].apply(lambda x: x.split("_")[1]).astype(int)
    
    gwas_df["chr"] = gwas_df["CHR"].astype(str)
    ecaviar_df["chr"] = ecaviar_df["chr"].astype(str)

    ecaviar_gwas = pd.merge(ecaviar_df, gwas_df, on=["chr", "pos"])

    return ecaviar_gwas

def main():

    # retrieve args
    args = parse_args()

    gwas_df = parse_gwas(args.gwas_file)

    ecaviar_output = merge_gwas(args.ecaviar_output_table, gwas_df)
    ecaviar_output = ecaviar_output.rename(columns = {"variant_id_x": "variant_id"})
    ecaviar_output.to_csv(args.output_filename, sep="\t", header=True, index=None)

if __name__ == "__main__":
    main()

