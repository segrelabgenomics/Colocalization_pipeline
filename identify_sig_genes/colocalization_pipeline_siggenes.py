#!/usr/bin/env python3

"""
Find significant genes in LD interval

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
from significant_genes import extract_significant_genes
import sys


def constants(args):
    """
    dictionary of pathways
    """
    gtex_basedir = args.base_dir_to_gtex_qtls

    constants_dict = {
        "sig_eqtl_dir": os.path.join(gtex_basedir, "GTEx_Analysis_v8_eQTL/"),
        "sig_sqtl_dir": os.path.join(gtex_basedir, "GTEx_Analysis_v8_sQTL/"),
        "eqtl_file_suffix": ".v8.signif_variant_gene_pairs.txt.gz",
        "sqtl_file_suffix": ".v8.sqtl_signifpairs.txt.gz"
    }
    return constants_dict


def format_dir(x: str) -> str:
    """
    adds '/' at end
    """
    if x[-1] != "/":
        return f"{x}/"
    return x


def format_variant(variant):
    """
    replaces ":" with "_"
    """
    return re.sub(":", "_", variant)


def parse_args():
    """
    arguments to input when running pipeline
    """
    parser = argparse.ArgumentParser(description="arguments to input when running pipeline")

    parser.add_argument("--gwas_variant",
                        type=str,help="gwas variant of interest")
    parser.add_argument("--qtl_type",
                        type=str,help="QTL Type")
    parser.add_argument("--gwas_trait",
                        type=str, help="gwas trait of interst")
    parser.add_argument("--ld_interval_directory",
                        type=str, help="LD interval directory")
    parser.add_argument("--sig_genes_directory",
                        type=str, help="Directory to place significant genes")
    parser.add_argument("--base_dir_to_gtex_qtls", type=str, required=True)
    parser.add_argument("--min_var_in_ld_window", type=str, help="Minimum number of variants present to extract significant genes", required=True,default = 5)


    return parser.parse_args()


def retrieve_ld_coordinates(ldDir, gwas_trait, variant):
    """
    extracts ld interval from file
    """
    with open(f"{ldDir}{gwas_trait}_{variant}_ld_interval.txt", "r") as f:
        lines = f.readlines()
        return int(lines[2].strip()), int(lines[3].strip())


def create_output_file(sig_genes_dir,
                       sig_genes_dict,
                       gwas_trait,
                       variant,
                       qtl_type,
                       option="sig"):
    """
    creates name of output file
    """
    output = f"{sig_genes_dir}{gwas_trait}_{variant}_{qtl_type}"
    if option == "sig":
        return f"{output}_significant_genes_list_ecaviar.txt"
    return f"{output}_empty_genes.txt"


def main():
    """
    main section
    """
    #retrieve arguments
    args = parse_args()

    #constants
    constants_dict = constants(args)

    # format directories
    ld_int_dir = format_dir(args.ld_interval_directory)
    sig_genes_dir = format_dir(args.sig_genes_directory)

    variant_writing = format_variant(args.gwas_variant)
    ld_start, ld_end = retrieve_ld_coordinates(ld_int_dir, 
                                               args.gwas_trait,
                                               variant_writing)
    chrom = variant_writing.split("_")[0]
    sig_genes_dict = extract_significant_genes(constants_dict,
                                               chrom,
                                               args.qtl_type,
                                               ld_start,
                                               ld_end,args.min_var_in_ld_window)
   
    if sig_genes_dict:
        sig_genes_output_file = create_output_file(sig_genes_dir,
                                                   sig_genes_dict,
                                                   args.gwas_trait,
                                                   variant_writing,
                                                   args.qtl_type)
        with open(sig_genes_output_file, "w") as output_file:
            for key,value in sig_genes_dict.items():
                output_file.write("%s %s\n" % (key,",".join(value)))
    else:
        print("not enough significant variants")
        sig_genes_output_file = create_output_file(sig_genes_dir,
                                                   sig_genes_dict,
                                                   args.gwas_trait,
                                                   variant_writing,
                                                   args.qtl_type, option="empty")
        with open(sig_genes_output_file, "w") as output_file:
            output_file.write(variant_writing+"\n")


if __name__ == "__main__":
    main()


