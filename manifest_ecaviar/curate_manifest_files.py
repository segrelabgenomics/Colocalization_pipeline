#!/usr/bin/env python3

import os
import os.path
import argparse
import re
import pandas as pd
import numpy as np

def format_directory(directory):
    """
    ensures there is a '/' at end of directory
    """
    if directory[-1] != "/":
        return(directory + "/")
    return(directory)

def retrieve_files(directory, gwas_trait):
    """
    returns list of files corresponding to gwas trait
    """
    return([directory + file for file in os.listdir(directory) if gwas_trait in file])

def parse_sig_genes_file(sig_genes_file, loci, qtl):
    """
    explodes dataframe
    adds header for loci
    """
    sig_gene_df = pd.read_csv(sig_genes_file, sep = " ", header = None)
    sig_gene_df.columns = ["tissue", "genes"]
    sig_gene_df_parsed = sig_gene_df.set_index(['tissue']).genes.str.split(',', expand = True).stack().reset_index().drop(['level_1'], axis = 1).rename(columns={0:'genes'})
    print(sig_gene_df_parsed)
    sig_gene_df_parsed["gwas_loci"] = loci
    sig_gene_df_parsed["qtl"] = qtl
    return(sig_gene_df_parsed)

def retrieve_loci_sig_genes(sig_genes_file):
    """
    indices 2-5 of basename of file
    """
    basename = os.path.basename(sig_genes_file)
    return(":".join(basename.split("_")[2:6]))

def retrieve_qtl_sig_genes(sig_genes_file):
    """
    index 6 of basename of file
    """
    basename = os.path.basename(sig_genes_file)
    print(basename)
    return(basename.split("_")[6])

def parse_sig_genes(sig_genes_files, plink_directory, plink_subset, gwas_trait, gwas_file):
    """
    returns df with headers tissue, sig_genes
    """
    sig_genes_df = pd.DataFrame()

    for sig_genes_file in sig_genes_files:
        if "empty_genes" in sig_genes_file:
            continue
        gwas_loci = retrieve_loci_sig_genes(sig_genes_file)
        qtl = retrieve_qtl_sig_genes(sig_genes_file)
        sig_genes_df_parsed = parse_sig_genes_file(sig_genes_file, gwas_loci, qtl)
        sig_genes_df = pd.concat([sig_genes_df, sig_genes_df_parsed], ignore_index=True)

    # necessary columns
    sig_genes_df["plink_directory"] = plink_directory
    sig_genes_df["plink_subset"] = plink_subset
    sig_genes_df["gwas_trait"] = gwas_trait
    sig_genes_df["gwas_file"] = gwas_file

    # reorder columns
    sig_genes_df = sig_genes_df[["gwas_loci", "gwas_file", "qtl", "gwas_trait", "tissue", "genes",
                                 "plink_directory", "plink_subset"]]

    return(sig_genes_df)

def chunk_manifest_file(sig_genes_df, gwas_trait):
    """
    chunk manifest files in chunks of 1000
    """
    length_df = len(sig_genes_df) + 100 # adding 100 to ensure no rows are left behind
    for num in np.arange(0, length_df, 1000):
        df = sig_genes_df[num:num+1000]
        if not df.empty:
            df.to_csv(args.manifest_chunks_dir + "MANIFEST_chunked_" + args.gwas_trait + "_" + str(num) + "_" + str(num+1000) + "_loci_tissue_gene.txt",
	              header = None, sep = ";", index = None)

def parse_args():
    """
    arguments to run
    """
    parser = argparse.ArgumentParser(description = "arguments to pass when running script")

    parser.add_argument("--gwas_file", type = str, help = "gwas file location", required = True)
    parser.add_argument("--sig_genes_directory", type = str, help = "directory containing significant genes", required = True)
    parser.add_argument("--gwas_trait", type = str, help = "gwas trait", required = True)
    parser.add_argument("--manifest_dir", type = str, help = "directory to manifest file", required = True)
    parser.add_argument("--manifest_chunks_dir", type = str, help = "directory to chunks of 1000 for manifest file", required = True)
    parser.add_argument("--plink_directory", type = str, help = "second plink directory", required = True)
    parser.add_argument("--plink_subset", type = str, help = "second plink subset", required = True)

    return(parser.parse_args())

if __name__ == "__main__":

    # retrieve arguments
    args = parse_args()

    #format directories
    sig_genes_directory = format_directory(args.sig_genes_directory)
    plink_directory = format_directory(args.plink_directory)

    # retrieve files from directories
    sig_genes_files = retrieve_files(sig_genes_directory, args.gwas_trait)
    
    sig_genes_df = parse_sig_genes(sig_genes_files, plink_directory, args.plink_subset, args.gwas_trait, args.gwas_file)
    print(sig_genes_df)
    # print whole manifest file for reference
    sig_genes_df.to_csv(args.manifest_dir +"MANIFEST_" + args.gwas_trait + "_loci_tissue_gene.txt", header = True, sep = "\t", index = None)

    # print manifest files in chunks of 1000
    chunk_manifest_file(sig_genes_df, args.gwas_trait)
