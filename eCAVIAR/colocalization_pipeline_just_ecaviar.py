#!/usr/bin/env python3

"""
Steps

1. Compute LD interval
2. liftOver GWAS
3. Find eGenes and sGenes 
4. Run eCAVIAR

"""
from __future__ import division
import os
import os.path
import re
import datetime
import argparse
import subprocess
import pandas as pd
import numpy as np
import ecaviar_section_single_combination
import sys

def constants():
    """
    dictionary of pathways for GTEx files
    """
    basedir = os.environ.get("GTEX_QTL_DIR", "data/GTEx/V8_release/eQTLs/")
    constants_dict = {"eqtl_allpairs_directory": os.path.join(basedir, "GTEx_Analysis_v8_eQTL_all_associations/"),
                      "eqtl_allpairs_filename": f".allpairs.txt.gz",
                      "sqtl_allpairs_directory": os.path.join(basedir, "GTEx_Analysis_v8_sQTL_all_associations/"),
                      "sqtl_allpairs_filename": f".v8.sqtl_allpairs.txt.gz"}

    return(constants_dict)


def retrieve_current_date():
    """
    Retrieves the date.
    Date is in the form of MonthDay_Year, where month is abbreviated
    and time is based on a 24-hour clock (e.g. Mar21_2019_13:45:20).

    directory date will not include time, only MonthDay_Year.
    """
    date = datetime.datetime.now().strftime("%b%d_%Y")
    return (date)

def make_directory(directory_name):
    """
    internal function to create directory

    equivalent to mkdir in linux
    makes directory if not already exists
    """
    try:
        os.mkdir(directory_name)
    except FileExistsError:
        pass
    return(directory_name)

def create_output_directory(output_directory,gwas_trait,gwas_variant,qtl_type):
    """
    output directory
    """
    date = retrieve_current_date()
    output_directory = output_directory+"/"+gwas_trait+"_"+gwas_variant+"_"+qtl_type+"_"+date
    output_directory = make_directory(output_directory)
    return(output_directory+"/")

def format_gwas(gwas_df_file, chrom):
    """
    renames columns
    creates variant id

    window size
    """
    print("GWAS")
    gwas_df = pd.read_csv(gwas_df_file, sep="\t")
    rename_map = {"P-value": "gwas_p_value",
                  "mtag_pval": "gwas_p_value",
                  "P" :  "gwas_p_value",
                  "chr": "chromosome_ignore"}  # ignore
    if "variant_id" in gwas_df.columns:
        gwas_df = gwas_df.drop(columns=["variant_id"])
    # Rename 'variant' to 'variant_id'
    if "variant" in gwas_df.columns:
        gwas_df = gwas_df.rename(columns={"variant": "variant_id"})

    gwas_df = gwas_df.rename(columns=rename_map)
    print(gwas_df.head())
    #filter on window size
    gwas_df["CHR"] = gwas_df["CHR"].astype(str)  # ensure it's string type
    gwas_df["CHR"] = gwas_df["CHR"].apply(lambda x: x if x.startswith("chr") else "chr" + x)

    gwas_df = gwas_df.loc[gwas_df["CHR"] == chrom]
    

    print(gwas_df.head())
    return(gwas_df)


def format_variant(variant):
    """
    replaces ':' with '_'
    """
    return(re.sub(':', '_', variant))

def parse_args():
    """
    arguments to input when running pipeline
    """
    parser = argparse.ArgumentParser(description="arguments to input when running pipeline")

    #related to gwas
    parser.add_argument("--gwas_variant",type=str,help="gwas variant of interest",required=True)
    parser.add_argument("--gwas_file",type=str,help="gwas file",required=True)
    parser.add_argument("--qtl_type",type=str,help="QTL Type",required=True)
    parser.add_argument("--gwas_trait", type = str, help = "gwas trait of interst", required = True)
    parser.add_argument("--tissue", type = str, help = "tissue", required = True)
    parser.add_argument("--gene", type = str, help = "gene", required = True)
    parser.add_argument("--output_directory",
                        type=str, help="Base output directory")
    parser.add_argument("--ld_directory",
                        type=str,
                        help="Directory containing LD intervals files")


    #plink files
    parser.add_argument("--plink_directory",type=str,help="Directory containing plink files",required=True)
    parser.add_argument("--plink_subset",type=str,help="Population subset",required=True)

    return(parser.parse_args())

if __name__ == "__main__":

    #retrieve arguments
    args = parse_args()

    #constants
    constants_dict = constants()

    #variant chromosome
    chrom = args.gwas_variant.split(":")[0]

    # variant format for writing
    variant_writing = format_variant(args.gwas_variant)

    #output directory
    output_directory = create_output_directory(args.output_directory,
                                               args.gwas_trait,
                                               args.gwas_variant,
                                               args.qtl_type)

    with open(args.ld_directory+"/"+args.gwas_trait+"_"+variant_writing+"_ld_interval.txt", "r") as f:
        lines = f.readlines()
        variant = str(lines[1].strip())
        ld_start = int(lines[2].strip())
        ld_end = int(lines[3].strip())

    gwas_df = format_gwas(args.gwas_file, chrom)
    # print(args.plink_directory)
    # print(args.plink_directory[-1])
    num_variants = ecaviar_section_single_combination.ecaviar_section(chrom, constants_dict, output_directory, args.tissue, args.gene, variant_writing, variant, gwas_df, args.gwas_trait, ld_start, ld_end, args.plink_directory, args.plink_subset, args.qtl_type)
