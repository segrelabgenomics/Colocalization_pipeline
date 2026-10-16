#!/usr/bin/env python3

"""
Steps
Computes LD interval for a given loci

"""

import sys
import os
import os.path
import re
import argparse
import subprocess
import pandas as pd
import numpy as np
from plink_ld_interval import compute_ld_interval

from typing import List
from argparse import Namespace, ArgumentParser
from pandas import DataFrame


def format_variant(x: str) -> str:
    """
    replaces ":" with "_"
    """
    return re.sub(":", "_", x)


def parse_args():
    """
    arguments to input when running pipeline
    """
    parser = ArgumentParser(description="Arguments for computing ld interval")

    #related to gwas
    parser.add_argument("--gwas_variant",
                        type=str,help="GWAS variant",
                        required=True)
    parser.add_argument("--gwas_trait",
                        type=str, help="GWAS trait",
                        required=True)
    parser.add_argument("--ld_interval_directory",
                        type=str, help="Directory to place LD intervals",
                        required=True)

    #plink files
    parser.add_argument("--plink_directory",
                        type=str, help="Directory containing plink files",
                        required=True)
    parser.add_argument("--plink_subset",
                        type=str,help="Population subset",
                        required=True)
    parser.add_argument("--second_plink_directory",
                        type=str, help="Second directory for plink files for ld proxy",
                        required = True)
    parser.add_argument("--second_plink_subset", 
                        type=str, help="Second population subset for ld proxy",
                        required=True)

    return parser.parse_args()


def format_dir(x: str) -> str:
    """
    adds '/' at end
    """
    if x[-1] != "/":
        return f"{x}/"
    return x


def main() -> None:
    """
    main section
    """
    #retrieve arguments
    args = parse_args()

    # fix directories
    ld_dir = format_dir(args.ld_interval_directory)
    plink_dir = format_dir(args.plink_directory)
    second_plink_dir = format_dir(args.second_plink_directory)

    #variant chromosome
    chrom = args.gwas_variant.split(":")[0]

    # variant format for writing
    variant_writing = format_variant(args.gwas_variant)

    ld_start,ld_end,variant = compute_ld_interval(ld_dir,
                                                  chrom,
                                                  args.gwas_variant,
                                                  args.gwas_trait,
                                                  plink_dir,
                                                  args.plink_subset,
                                                  second_plink_dir,
                                                  args.second_plink_subset)

    filename = f"{ld_dir}{args.gwas_trait}_{variant_writing}_ld_interval.txt"
    with open(filename, "w") as f:
        f.write(args.gwas_variant+"\n"+variant+"\n"+str(ld_start)+"\n"+str(ld_end))


if __name__ == "__main__":
    main()


