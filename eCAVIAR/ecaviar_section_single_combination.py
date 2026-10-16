import os
import os.path
import sys
import subprocess
import pandas as pd
import numpy as np

def retrieve_bfile_ecaviar(chrom,plink_directory):
    """
    identify if plink directory contains /
    add kgp.clean. + chromosome
    """
    if plink_directory[-1] != "/":
        plink_directory+="/"
    return(plink_directory+"kgp.clean."+chrom)

def create_plink_out_name_ecaviar(output_directory,trait_name,variant,tissue_name,qtl_type,gene):
    return(output_directory+trait_name+"_"+gene+"_"+qtl_type+"_"+variant+"_"+tissue_name)

def call_plink_ld_matrix(chrom, output_directory,qtl_type,trait_name,tissue_name,variant,plink_directory,subset_file,variant_list,
              gene, memory="16000", threads="1"):
    """
    run plink to calculate ld matrix
    """
    bfile = retrieve_bfile_ecaviar(chrom,plink_directory)

    out_name = create_plink_out_name_ecaviar(output_directory,trait_name,variant,tissue_name,qtl_type,gene)

    ld_matrix_list = ["module load module-autoload/1.0 && module load PLINK && plink",
    "--bfile",bfile,"--keep",subset_file,
    "--r","--matrix","--extract",variant_list,"--write-snplist",
    "--memory",memory,"--threads",threads,"--out",out_name]

    cmd_str = " ".join(ld_matrix_list)+ " && sleep 10 && sync" 
    subprocess.call(cmd_str, shell = True)
    subprocess.call("sleep 10 && sync", shell=True) 

    return(out_name)

def filter_matrix(plink_snplist_file_name,plink_ld_file_name):
    """
    from the 201 variants inputted into plink, some are removed
    reformats the list of variants
    """
    try:
        plink_snplist = pd.read_csv(plink_snplist_file_name,header=None)
    except FileNotFoundError:
        return("","")
    plink_snplist.columns = ["variant_id"]

    plink_ld_matrix = pd.read_csv(plink_ld_file_name,header=None,sep=" ")

    a = {}
    for index,i in enumerate(plink_ld_matrix.columns):
        temp = plink_ld_matrix.iloc[:,i]
        if pd.isnull(temp).all():
            a[i] = "null"

    snplist_matrix = plink_snplist.join(plink_ld_matrix)

    x = []

    for i in a.keys():
        if i == snplist_matrix.columns[-1]:
            snplist_matrix = snplist_matrix.drop(i,axis=1)
            x.append(i)

    for j in x:
        del a[j]

    for i in a.keys():
        snplist_matrix = snplist_matrix.drop(i,axis=0)
        snplist_matrix = snplist_matrix.drop(i,axis=1)

    snplist = snplist_matrix[["variant_id"]].copy()
    snplist_matrix = snplist_matrix.drop(["variant_id"],axis=1)
    snplist_matrix = snplist_matrix.fillna(0)

    if snplist_matrix.shape[0] == snplist_matrix.shape[1]:
        print("Matrix is correct")

    return(snplist,snplist_matrix)

def filter_z_score_df(snplist,z_score_df):
    """
    remove variants from z_score_df that are not included in
    ld matrix
    """
    snplist_z_score_df = pd.merge(snplist,z_score_df,on=["variant_id"])
    snplist_z_score_df["pos"] = snplist_z_score_df["variant_id"].apply(lambda x: x.split(":")[1]).astype(int)
    snplist_z_score_df = snplist_z_score_df.sort_values(["pos"],ascending=True)
    snplist_z_score_df = snplist_z_score_df.drop(["pos"],axis=1)
    return(snplist_z_score_df)

def extract_plink_ld_file_name(out_name):
    return(out_name+".ld")

def extract_plink_snplist_file_name(out_name):
    return(out_name+".snplist")

def reformat_ecaviar_input_files(out_name,eqtl_z_score,gwas_z_score,variant):

    plink_ld_file_name = extract_plink_ld_file_name(out_name)
    plink_snplist_file_name = extract_plink_snplist_file_name(out_name)

    snplist,snplist_matrix  = filter_matrix(plink_snplist_file_name,plink_ld_file_name)
    if isinstance(snplist,str):
        return("","","")
    snplist_matrix.to_csv(plink_ld_file_name,sep=" ",index=None,header=None)

    eqtl_z_score_snplist = filter_z_score_df(snplist,eqtl_z_score)
    eqtl_z_score_snplist["variant_id"] = eqtl_z_score_snplist["variant_id"].str.replace("_b38","")
    eqtl_z_score_snplist["variant_id"] = eqtl_z_score_snplist["variant_id"].str.replace("_",":")
    eqtl_z_score_snplist["variant_id"] = eqtl_z_score_snplist["variant_id"].str.replace("chr","")

    gwas_z_score_snplist = filter_z_score_df(snplist,gwas_z_score)
    gwas_z_score_snplist["variant_id"] = gwas_z_score_snplist["variant_id"].str.replace("_b38","")
    gwas_z_score_snplist["variant_id"] = gwas_z_score_snplist["variant_id"].str.replace("_",":")
    gwas_z_score_snplist["variant_id"] = gwas_z_score_snplist["variant_id"].str.replace("chr","")

    return(plink_ld_file_name,eqtl_z_score_snplist,gwas_z_score_snplist)

def parse_variant(variant):
    """
    extract chromosome and position from variant
    """
    if "_" not in variant:
        raise ValueError("{} not in correct format."
                         "Variant need be in form chrN_pos_ref_alt_genomebuild."
                         "e.g. chr10_123422_A_G_b38".format(variant))
    chrom,position = variant.split("_")[0:2]
    return (chrom,int(position))

def ensure_pos_dtype(gwas):
    """
    ensure position header consists of integers
    """
    if gwas["pos"].dtypes == int:
        return(gwas["pos"])
    return(gwas["pos"].astype(int))

def prepare_gwas(gwas):
    """
    formats gwas file
    """
    #format chr header
    gwas["pos"] = ensure_pos_dtype(gwas)

    gwas = gwas.loc[gwas["zscore"].notna()]
    gwas = gwas.sort_values(["pos"],ascending=True).reset_index(drop=True)
    gwas = gwas.drop_duplicates(["CHR","pos"])

    return gwas

def parse_gtex_file(tissue,gene_id):

    #tissue["gene_id"] = tissue["gene_id"].str.extract(r".*(ENSG\d+)")

    tissue["Z_eQTL"] = tissue["slope"]/tissue["slope_se"]
    tissue = tissue.loc[tissue["Z_eQTL"].notna()].copy()
    tissue["chr"] = tissue["variant_id"].apply(lambda x: x.split(":")[0])
    tissue["position"] = tissue["variant_id"].apply(lambda x: x.split(":")[1])
    tissue["position"] = tissue["position"].astype(int)
    return(tissue)

def merge_with_gwas(tissue_file,gene_id,gwas_file):

    tissue = parse_gtex_file(tissue_file,gene_id)

    print("prepare gwas")
    gwas = prepare_gwas(gwas_file)
    #gwas["variant_id"] = gwas["variant_id"].str.replace(":", "_") + "_b38"
    print(gwas.head())
    tissue_gwas = pd.merge(tissue,gwas,on=["variant_id"])
    tissue_gwas = tissue_gwas.drop_duplicates(["variant_id"])
    print("Merged Tissue-GWAS")
    print(tissue_gwas.head())
    #delete old variables to save memory
    del tissue,gwas

    return tissue_gwas

def harmonize_variant_ids(gwas, qtl, gwas_col="variant_id", qtl_col="variant_id"):
    print(gwas.columns)
    print(qtl.columns)

    print(gwas.head())

    print(qtl.head())

    gwas_unique = set(gwas[gwas_col].unique())
    qtl_unique = set(qtl[qtl_col].unique())

    n_gwas = len(gwas_unique)
    n_qtl = len(qtl_unique)
    intersection = gwas_unique & qtl_unique
    n_intersection = len(intersection)

    print(f"Unique variants in GWAS file: {n_gwas}")
    print(f"Unique variants in tissue (QTL) file: {n_qtl}")
    print(f"Intersection (shared variant IDs): {n_intersection}")

    gwas_split = gwas[gwas_col].str.split(':', expand=True)
    gwas_split.columns = ['chr', 'pos', 'ref', 'alt']
    gwas_split['index'] = gwas.index

    qtl_split = qtl[qtl_col].str.split(':', expand=True)
    qtl_split.columns = ['chr', 'pos', 'ref', 'alt']
    qtl_split['index_qtl'] = qtl.index

    # Normalize chromosomes and alleles
    for df in (gwas_split, qtl_split):
        df['ref'] = df['ref'].str.upper()
        df['alt'] = df['alt'].str.upper()

    # Merge on chr and pos
    merged = gwas_split.merge(
        qtl_split, on=['chr', 'pos'], how='inner', suffixes=('_gwas', '_qtl'))

    # Identify swapped alleles
    swapped_mask = (
        (merged['ref_gwas'] == merged['alt_qtl']) &
        (merged['alt_gwas'] == merged['ref_qtl']))
    
    print("\n=== DEBUG: Checking swapped_mask ===")
    print(f"Total merged variants: {merged.shape[0]}")
    print(f"Swapped TRUE count: {swapped_mask.sum()}")
    print(f"Swapped FALSE count: {(~swapped_mask).sum()}")
    
    debug_cols = ['chr', 'pos', 'ref_gwas', 'alt_gwas', 'ref_qtl', 'alt_qtl']
    print("\nFirst 10 merged variants:")
    print(merged[debug_cols].head(10))

    print("\nRows where swapped_mask == True (if any):")
    print(merged.loc[swapped_mask, debug_cols])
 
    # Build mapping of old and new IDs
    switched = merged.loc[swapped_mask, ['index', 'index_qtl']]
    switched = switched.assign(
        old_id=gwas.loc[switched['index'], gwas_col].values,
        new_id=qtl.loc[switched['index_qtl'], qtl_col].values
    )[['old_id', 'new_id']]

    # Apply replacements
    gwas_updated = gwas.copy()
    if not switched.empty:
        gwas_updated[gwas_col] = gwas_updated[gwas_col].replace(
            dict(zip(switched['old_id'], switched['new_id'])))

    print(f"GWAS-QTL harmonization complete. {switched.shape[0]} variant IDs were swapped.")
    if not switched.empty:
        print(switched.head())

    return gwas_updated

def parse_tissue_gwas(tissue_file,gene_id,gwas_file,ld_start,ld_end):
    gwas_harmonized = harmonize_variant_ids(gwas_file, tissue_file, 
                    gwas_col="variant_id", qtl_col="variant_id")

    print("Merging GWAS and QTL summary statistics")
    tissue_gwas = merge_with_gwas(tissue_file,gene_id,gwas_harmonized)
    tissue_gwas = tissue_gwas.sort_values(["position"]).reset_index(drop=True)

    tissue_gwas_ld = tissue_gwas.loc[(tissue_gwas["position"] >= ld_start) & (tissue_gwas["position"] <= ld_end)].copy()

    if len(tissue_gwas_ld) > 1000:
        tissue_gwas_ld = tissue_gwas_ld.loc[tissue_gwas_ld["gwas_p_value"] <= 0.05]
        if len(tissue_gwas_ld) > 3500:
            tissue_gwas_ld = tissue_gwas_ld.loc[tissue_gwas_ld["gwas_p_value"] <= 0.01]
            if len(tissue_gwas_ld) > 4500:
                tissue_gwas_ld = tissue_gwas_ld.loc[tissue_gwas_ld["gwas_p_value"] <= 0.005]
                if len(tissue_gwas_ld) > 4500:
                    tissue_gwas_ld = tissue_gwas_ld.loc[tissue_gwas_ld["gwas_p_value"] <= 0.001]
                    

    eqtl_only = tissue_gwas_ld[["variant_id","Z_eQTL"]].copy()
    gwas_only = tissue_gwas_ld[["variant_id","zscore"]].copy()

    variants_list = np.array(tissue_gwas_ld["variant_id"])

    #delete old variables to save memory
    del tissue_gwas

    return(eqtl_only,gwas_only,variants_list)

def write_temp_file(variants_list,gwas_trait,tissue_name,variant,output_directory):
    #Create temp file name tf =tempfile.NamedTemporaryfile()
    name = output_directory+gwas_trait+"_"+tissue_name+"_"+variant+"_temp"
    print("TEMP processed!")
    #Write variants into a list file with tf.name as the name of the file
    with open(name, mode="wt") as variant_file:
        variant_file.write("\n".join(variants_list))
    return(name)


def create_file(constants_dict, qtl_type, tissue):
    """
    creates allpairs filename
    """
    if qtl_type.lower() == "eqtl":
        tissue_dir = constants_dict["eqtl_allpairs_directory"]
        filename = constants_dict["eqtl_allpairs_filename"]
        print(f"{tissue_dir}{tissue}{filename}")
        return f"{tissue_dir}{tissue}{filename}"
    else:
        tissue_dir = constants_dict["sqtl_allpairs_directory"]
        filename = constants_dict["sqtl_allpairs_filename"]
        return f"{tissue_dir}{tissue}{filename}"


def ecaviar_section(chrom, constants_dict, output_directory, tissue, gene, gwas_variant, second_variant, gwas_df, gwas_trait, ld_start, ld_end, plink_directory, plink_subset, qtl_type):
    """
    main for ecaviar section of pipeline
    """
    filename = create_file(constants_dict, qtl_type, tissue)
    header = subprocess.check_output(["bash", "-c", f"zcat {filename} | head -n 1"]).decode().strip().split("\t")
    grep_output = subprocess.Popen(["zgrep", gene, filename], stdout=subprocess.PIPE)
    tissue_df = pd.read_csv(grep_output.stdout, sep="\t", header=None, names=header)
    print("Tissue df") 
    tissue_df = tissue_df[["gene_id", "variant_id",
                            "pval_nominal",
                            "slope", "slope_se"]]
    
    print(tissue_df.head())

    tissue_df["variant_id"] = tissue_df["variant_id"].str.replace("_", ":").str.replace(":b38", "")


    eqtl_only, gwas_only, variants_list = parse_tissue_gwas(tissue_df, gene, gwas_df, ld_start, ld_end)
    
    #plink
    subset_file = write_temp_file(variants_list,gwas_trait,tissue,second_variant,output_directory)
    out_name = call_plink_ld_matrix(chrom, output_directory,qtl_type,gwas_trait,tissue,second_variant,plink_directory,plink_subset,subset_file,gene)

    #ecaviar
    plink_ld_file_name,eqtl_z_score_snplist,gwas_z_score_snplist = reformat_ecaviar_input_files(out_name,eqtl_only,gwas_only,second_variant)
    if isinstance(eqtl_z_score_snplist,str):
        sys.exit(0)
    eqtl_z_score_snplist.to_csv(out_name+"_qtl_z_score.txt",index=None,header=None,sep="\t")
    gwas_z_score_snplist.to_csv(out_name+"_gwas_z_score.txt",index=None,header=None,sep="\t")

    print("running eCAVIAR")
    ecaviar_out_name = out_name+"_eCAVIAR"
    
    eCAVIAR_list = ["module load caviar && eCAVIAR",
                    "-o",ecaviar_out_name,
                    "-l",plink_ld_file_name,
                    "-l",plink_ld_file_name,
                    "-z",out_name+"_gwas_z_score.txt",
                    "-z",out_name+"_qtl_z_score.txt"]
    
    #subprocess.call(eCAVIAR_list)
    cmd_str = " ".join(eCAVIAR_list)+ " && sleep 10 && sync" 
    print(cmd_str)
    subprocess.call(cmd_str, shell = True)
    subprocess.call("sleep 10 && sync", shell=True,executable='/bin/bash') 

    # remove ld matrix
    os.remove(plink_ld_file_name)

    return(len(eqtl_z_score_snplist))
