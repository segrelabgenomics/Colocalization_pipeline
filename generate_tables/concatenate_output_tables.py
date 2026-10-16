#!/usr/bin/env python3

import os
import os.path
import re
import argparse
import pandas as pd

from pandas import DataFrame
from argparse import ArgumentParser, Namespace
from typing import Dict, List

def parse_args() -> Namespace:
    """
    arguments to pass when running scripts
    """
    parser = ArgumentParser(description="arguments to pass")

    parser.add_argument("--base_dir_eqtl",
                        type=str, help="Directory containing outputs for eqtls")
    parser.add_argument("--base_dir_sqtl",
                        type=str, help="Directory containing outputs for sqtls")
    parser.add_argument("--output_directory",
                        type=str, help="Directory to place parsed files")
    parser.add_argument("--output_filename",
                        type=str, help="name of output file")

    return parser.parse_args()

def retrieve_files(dir: str) -> List[str]:
    """
    returns list of files
    """
    return [dir+x for x in os.listdir(dir) if "full_qtl_info" in x]

def concatenate_ecaviar_outputs(ecaviar_tables: List[str]) -> DataFrame:
    """
    concatenates all tables
    """
    full_df = pd.DataFrame()

    for table in ecaviar_tables:
        df = pd.read_csv(table, sep="\t")
        full_df = full_df.append(df)
    return full_df

def parse_eqtl_table(df: DataFrame) -> DataFrame:
    """
    removes gene columns
    """
    #df = df.drop(["num_variants_z_score_file", "gene"], axis=1)
    return df

def parse_sqtl_table(df: DataFrame) -> DataFrame:
    """
    reconfigures intron cluster
    """
    df = df.rename(columns={"gene_id": "Gene_id_intron_cluster",
                            "gene": "gene_id"})
    return df

def concatenate_full_table(ecaviar_concatenated_df_eqtl: DataFrame,
                           ecaviar_concatenated_df_sqtl: DataFrame) -> DataFrame:
    """
    concatenates two tables
    """
    df_eqtl_parsed = parse_eqtl_table(ecaviar_concatenated_df_eqtl)
    df_sqtl_parsed = parse_sqtl_table(ecaviar_concatenated_df_sqtl)

    return df_eqtl_parsed.append(df_sqtl_parsed)

def main():
    """
    main function
    """
    args = parse_args() # retrieve arguments

    ecaviar_tables_eqtl = retrieve_files(args.base_dir_eqtl)
    ecaviar_tables_sqtl = retrieve_files(args.base_dir_sqtl)

    ecaviar_concatenated_df_eqtl = concatenate_ecaviar_outputs(ecaviar_tables_eqtl)
    ecaviar_concatenated_df_sqtl = concatenate_ecaviar_outputs(ecaviar_tables_sqtl)

    ecaviar_concatenated_df = concatenate_full_table(ecaviar_concatenated_df_eqtl, ecaviar_concatenated_df_sqtl)

    ecaviar_concatenated_df.to_csv(args.output_directory+args.output_filename,
                                   sep="\t", header=True, index=None, na_rep="NA")

if __name__ == "__main__":
    main()



