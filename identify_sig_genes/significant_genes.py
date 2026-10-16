import os
import os.path
import re
import pandas as pd


def identify_genes_ld_interval(filename,chrom,ld_start,ld_end,min_var_in_ld_window):
    """
    finds genes in ld interval
    """
    df = pd.read_csv(filename,sep="\t")
    df = df[df["variant_id"].notna()
                    & df["variant_id"].str.contains("_", na=False)
                    & ~df["variant_id"].str.contains(r"\s", na=False)]
    df["chr"] = df["variant_id"].apply(lambda x: x.split("_")[0])
    df["pos"] = df["variant_id"].apply(lambda x: x.split("_")[1]).astype(int)

    df = df.loc[df["chr"] == chrom].copy()
    print("QTL processed")
    #check if variants in window
    df = df.loc[(df["pos"] >= ld_start) & (df["pos"] <= ld_end)].copy()
    num = int(df["variant_id"].nunique())
    print("Number of variants in LD Window: " + str(num))
    if num > int(min_var_in_ld_window):
        return list(df["gene_id"].unique())
    return ""


def retrieve_files(constants_dict, qtl_type):
    """
    returns GTEx files from directory given qtl type
    """
    if qtl_type.lower() == "eqtl":
        files = os.listdir(constants_dict["sig_eqtl_dir"])
        return [file for file in files if file.endswith(constants_dict["eqtl_file_suffix"])]

    if qtl_type.lower() == "sqtl":
        files = os.listdir(constants_dict["sig_sqtl_dir"])
        return [file for file in files if file.endswith(constants_dict["sqtl_file_suffix"])]


def retrieve_tissue(constants_dict, gtex_file, qtl_type):
    """
    extracts tissue
    """
    if qtl_type.lower() == "eqtl":
        return re.sub(constants_dict["eqtl_file_suffix"], "", gtex_file)

    if qtl_type.lower() == "sqtl":
        return re.sub(constants_dict["sqtl_file_suffix"], "", gtex_file)


def create_filename(constants_dict, gtex_file, qtl_type):
    """
    creates filename
    """
    if qtl_type.lower() == "eqtl":
        return f"{constants_dict['sig_eqtl_dir']}{gtex_file}"

    if qtl_type.lower() == "sqtl":
        return f"{constants_dict['sig_sqtl_dir']}{gtex_file}"
 
def extract_significant_genes(constants_dict,chrom,qtl_type,ld_start,ld_end,min_var_in_ld_window):
    """
    determine list of significant genes
    """
    sig_tissues_dict = {}

    gtex_files = retrieve_files(constants_dict, qtl_type)

    for gtex_file in gtex_files:
        tissue = retrieve_tissue(constants_dict, gtex_file, qtl_type)
        filename = create_filename(constants_dict, gtex_file, qtl_type)
        print(tissue)
        significant_genes = identify_genes_ld_interval(filename,chrom,ld_start,ld_end,min_var_in_ld_window)
        if significant_genes:
            sig_tissues_dict[tissue] = significant_genes

    return sig_tissues_dict
