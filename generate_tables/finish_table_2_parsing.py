#!/usr/bin/env python3
import os
import os.path
import pandas as pd
import numpy as np
import re
import gzip

class gene:
    def __init__(self, chr, gene_name, gene_id, gene_start, gene_finish):
        self.gene_name = gene_name
        self.gene_id = gene_id
        self.start = int(gene_start)
        self.end = int(gene_finish)
        self.chr = chr
        self.gene_range = range(self.start, self.end)

#Parses gtf file, returns 'gene' objects
def parse_gtf(gtf_file, compression=None):
    """
    gtf_file: path to a GTF file, compressed or not
    """
    gene_dict = {}
    if compression == None:
        gtf = open(gtf_file)
    elif compression == 'gz':
        gtf = gzip.open(gtf_file, 'rb')
    else:
        raise Exception("Compression of gtf file not recognized. Please used "
                        ".gz or no compression.")
    with gtf:
        for row in gtf:
            if isinstance(row, bytes):
                row = row.decode('utf-8')
            row = row.strip().split('\t')
            if row[0][0]=='#' or row[2]!='gene': continue
            attr = dict([i.split() for i in row[8].replace('"','').split(';') if i!=''])
            gene_dict[attr['gene_id'].split(".")[0]]=(gene(row[0],
                      attr['gene_name'],attr['gene_id'].split('.')[0], row[3], row[4]))
    return(gene_dict)

def parse_gencode_gtf(gtf_file):
    """
    manages parse_gtf function

    ARH 5 Mar 2021
    """
    if gtf_file.endswith('gz'):
        gene_dict = parse_gtf(gtf_file, 'gz')
        return(gene_dict)
    elif gtf_file.endswith('gtf'):
        gene_dict = parse_gtf(gtf_file)
        return(gene_dict)
    else:
        raise Exception('GTF file passed cannot be parsed. Please pass a .gtf or .gz.')

def extract_gene_name(gene_dict, ensembl_gene_id):
    """
    returns gene name from ensembl gene id
    """
    try:
        gene_name = gene_dict[ensembl_gene_id].gene_name
        return(gene_name)
    except:
        return('NA')

def extract_ensembl_gene(gene_dict, gene_name):
    """
    returns gene name from ensembl gene id
    """
    for i in gene_dict:
        if gene_dict[i].gene_name == gene_name:
            return i
    return "NA"

gencode_file = "gencode.v26.annotation.gtf"
gencode_dict = parse_gencode_gtf(gencode_file)
gencode_dict['ENSG00000223972'].gene_id

## Eye disease gwas table parsing
eye_disease_file = "../data/gwas.txt"
eye_disease_df = pd.read_csv(eye_disease_file, sep=" ")
eye_disease_df["Genome build"] = "b38"
display(eye_disease_df)

eye_disease_df["Allele1"] = eye_disease_df["EA"].str.upper()
eye_disease_df["Allele2"] = eye_disease_df["NEA"].str.upper()

for column in eye_disease_df.columns:
        if column in ("p_2gc","gwas_p_value", "mtag_pval","P-value","P"):
            eye_disease_df = eye_disease_df.rename(columns = {column: "GWAS_p_value"}) 
        if column in ("rsid","SNP","ID","MarkerName","rsID"):
            eye_disease_df = eye_disease_df.rename(columns = {column: "MARKERNAME"})
        # Uncomment in case no "BP"
        if column in ("POS"):
            eye_disease_df = eye_disease_df.rename(columns = {column: "BP"})

# Ensure numeric columns are numeric before using them
eye_disease_df["CHR"] = eye_disease_df["CHR"].astype(str)
eye_disease_df["BP"] = eye_disease_df["BP"].astype(str)

# Prefix chromosome and build variant ID
eye_disease_df["chr"] = "chr" + eye_disease_df["CHR"]
eye_disease_df["variant_id"] = (
    eye_disease_df["chr"] + ":" +
    eye_disease_df["BP"] + ":" +
    eye_disease_df["Allele2"].astype(str) + ":" +
    eye_disease_df["Allele1"].astype(str))

eye_disease_df = eye_disease_df.rename(columns={"MARKERNAME": "rsID_GWAS",
                                                "chromosome": "CHR",
                                                "BP": "BP_Build38",
                                                "GWAS_p_value":"P-VALUE",
                                                "variant_id":"GWAS variant",
                                                "Allele1":"Effect allele",
                                                "Allele2":"Other Allele",
                                                "BETA":"OR or BETA",
                                               })

eye_disease_df["P-VALUE"] = eye_disease_df["P-VALUE"].astype(float)
eye_disease_df.loc[eye_disease_df["P-VALUE"] == 0, "P-VALUE"] = 1e-177
eye_disease_df = eye_disease_df.reset_index(drop=True)
eye_disease_df = eye_disease_df[["GWAS variant", "rsID_GWAS", "Effect allele", "Other Allele",
                                 "P-VALUE", "OR or BETA"]]
eye_disease_df.columns = ["GWAS variant", "rsID_GWAS", "Effect allele", "Other Allele","GWAS_p_value", "OR or BETA"]
eye_disease_df["OR or BETA"] = eye_disease_df["OR or BETA"].astype(float)


if (eye_disease_df["OR or BETA"] >= 0).all():
    print("OR")
    eye_disease_df = eye_disease_df.rename(columns={"OR or BETA": "GWAS_OR"})
    eye_disease_df["GWAS_effect"] =  np.log(eye_disease_df["GWAS_OR"])
else:
    print("Beta")
    eye_disease_df = eye_disease_df.rename(columns={"OR or BETA": "GWAS_effect"})
    eye_disease_df["GWAS_OR"] = np.exp(eye_disease_df["GWAS_effect"]).astype(float)
eye_disease_df["GWAS variant"] = eye_disease_df["GWAS variant"].str.replace(":","_")+"_b38"
eye_disease_df.head()

## LD interval parsing
def parse_ld_interval_file(file):
    """
    extracts ld interval and loci
    """
    with open(file, "r") as f:
        lines = f.readlines()
        df = pd.DataFrame({"GWAS variant": [lines[0].strip()],
                           "GWAS variant (aligned to 1KG)": [lines[1].strip()],
                           "LD interval start (b38)": [int(float(lines[2].strip()))],
                           "LD interval end (b38)": [int(float(lines[3].strip()))]
                          })
        df["LD interval size (kb)"] = (df["LD interval end (b38)"]-df["LD interval start (b38)"])/1000
        df["GWAS Trait"] = "POAG_EUR"
    return df

# %% [code] cell 13
ldint_dir = "../ld_intervals/POAG_EUR/"
ldint_files = [ldint_dir+file for file in os.listdir(ldint_dir) if ".txt" in file]

ld_interval_df = pd.DataFrame()

for file in ldint_files:
    #print(file)
    df = parse_ld_interval_file(file)
    #print(df.head())
    ld_interval_df = ld_interval_df.append(df)
print(ld_interval_df.head())    
ld_interval_df = ld_interval_df[["GWAS Trait", "GWAS variant", 
                                 "GWAS variant (aligned to 1KG)",
                                 "LD interval start (b38)", "LD interval end (b38)",
                                "LD interval size (kb)"
                                ]]
ld_interval_df = ld_interval_df.reset_index(drop=True)
ld_interval_df = ld_interval_df.drop_duplicates(["GWAS variant"])
ld_interval_df = ld_interval_df.reset_index(drop=True)
ld_interval_df["CHR"] = ld_interval_df["GWAS variant"].apply(lambda x: x.split(":")[0])
ld_interval_df["BP_Build38"] = ld_interval_df["GWAS variant"].apply(lambda x: x.split(":")[1])
ld_interval_df["GWAS variant"]=ld_interval_df["GWAS variant"].str.replace(":", "_")
ld_interval_df["GWAS variant"]=ld_interval_df["GWAS variant"]+ "_b38"

ld_interval_df["ID1"] = ld_interval_df["GWAS variant"]
ld_interval_df["ID2"] = ld_interval_df["GWAS variant (aligned to 1KG)"].str.replace(":","_")
ld_interval_df["ID2"] = ld_interval_df["ID2"] + "_b38"

ld_long = ld_interval_df.melt(
    id_vars=[col for col in ld_interval_df.columns if col not in ['ID1', 'ID2']],
    value_vars=['ID1', 'ID2'],
    var_name='LD_ID',
    value_name='GWAS variant to merge')
ld_long = ld_long.drop_duplicates(subset=['GWAS variant to merge'] + [col for col in ld_interval_df.columns if col not in ['ID1', 'ID2']])
eye_disease_df['GWAS variant to merge'] = eye_disease_df['GWAS variant']
eye_disease_df = pd.merge(eye_disease_df, ld_interval_df, on='GWAS variant',how="right")

eye_disease_df = eye_disease_df.rename(columns={'GWAS variant_x': 'GWAS variant'})
eye_disease_df["GWAS variant"]=eye_disease_df["GWAS variant"].str.replace("_", ":")
eye_disease_df["GWAS variant"]=eye_disease_df["GWAS variant"].str.replace(":b38", "")
eye_disease_df

eye_disease_df = eye_disease_df.rename(columns={'CHR_x': 'CHR'})
eye_disease_df["CHR"] = eye_disease_df["CHR"].astype(str)
eye_disease_df['ABS GWAS effect'] = eye_disease_df['GWAS_effect'].astype(float).abs()
eye_disease_df = eye_disease_df.rename(columns={"GWAS variant": "GWAS Loci"})
eye_disease_df.head()

## eCAVIAR Parsing
ecaviar_file = "../Intermediate_files/results_concatenated/POAG_EUR_pretable1.tsv"
table1_file = "POAG_EUR_table1_formatted_with_QC_UPD_renamed.tsv"
output_filename = "POAG_EUR_table2_Pass_QC_QTL_FDR05_P1E04_GWAS_P5E02.tsv"

ecaviar_df = pd.read_csv(ecaviar_file, sep="\t")

ecaviar_df_sig = pd.read_csv(table1_file, sep="\t")
ecaviar_df_sig.columns = ecaviar_df_sig.columns.str.replace('.', ' ') 
for column in ecaviar_df_sig.columns:
        if column in ("QTL_Type","QTL.type","QTL_type","QTL.Type"):
            eye_disease_df = eye_disease_df.rename(columns = {column: 'QTL Type'}) 
ecaviar_df_sig.head()

ecaviar_df["gene"] = ecaviar_df["gene"].str.split('.').str[0]
ecaviar_df["gene_id"] = ecaviar_df["gene_id"].str.split('.').str[0]
ecaviar_df["loci"] = ecaviar_df["loci"].str.replace("_b38", "").str.replace("_", ":")
ecaviar_df_eqtl = ecaviar_df.loc[ecaviar_df["qtl"].str.lower() == "eqtl"].copy()
ecaviar_df_sqtl = ecaviar_df.loc[ecaviar_df["qtl"].str.lower() == "sqtl"].copy()

ecaviar_df_eqtl["gene_name"] = ecaviar_df_eqtl["gene"].apply(lambda x: extract_gene_name(gencode_dict, x))

ecaviar_df_sig = ecaviar_df_sig.rename(columns = {"QTL_Type": 'QTL Type'}) 

ecaviar_df_sig = ecaviar_df_sig.loc[(ecaviar_df_sig["Pass_QC_QTL_FDR05_P1E04_GWAS_P5E02"] == True)].copy()
ecaviar_df_sig_eqtl = ecaviar_df_sig.loc[(ecaviar_df_sig["QTL Type"] == "eQTL") &
                                         (ecaviar_df_sig["CLPP"] >= 0.01)].copy()
ecaviar_df_sig_sqtl = ecaviar_df_sig.loc[(ecaviar_df_sig["QTL Type"] == "sQTL") &
                                         (ecaviar_df_sig["CLPP"] >= 0.01)].copy()

# ## Number eQTL-tissue pairs tested
df = ecaviar_df_eqtl[["loci", "tissue", "gene"]]
df = df.drop_duplicates(["loci", "tissue", "gene"])
df_ecaviar_eqtl_tissue_pairs = df["loci"].value_counts().rename_axis("GWAS Loci").reset_index(name="# eQTL-tissue pairs (eCAVIAR)")
df_ecaviar_eqtl_tissue_pairs.head()

df_ecaviar_eqtl_tissue_pairs = (
    eye_disease_df.merge(df_ecaviar_eqtl_tissue_pairs, on=["GWAS Loci"], how="left")
)

## Number eGenes tested
df = ecaviar_df_eqtl[["loci", "gene"]]
df = df.drop_duplicates(["loci", "gene"])
df_ecaviar_eqtl_tested = df["loci"].value_counts().rename_axis("GWAS Loci").reset_index(name="# eGenes tested (eCAVIAR)")
df_ecaviar_eqtl_tested.head()

# ## Number significant eQTL-tissue pairs tested
df = ecaviar_df_sig_eqtl[["GWAS_variant", "Tissue", "Gene_ENS_ID"]]
df = df.drop_duplicates(["GWAS_variant", "Tissue", "Gene_ENS_ID"])
df["GWAS_variant"] = df["GWAS_variant"].str.replace("_b38", "").str.replace("_", ":")
df_ecaviar_sig_eqtl_tissue_pairs = df["GWAS_variant"].value_counts().rename_axis("GWAS Loci").reset_index(name="# significant eQTL-tissue pairs (eCAVIAR)")
df_ecaviar_sig_eqtl_tissue_pairs.head()

## Number significant eGenes
df = ecaviar_df_sig_eqtl[["GWAS_variant", "Gene_ENS_ID"]]
df = df.drop_duplicates(["GWAS_variant", "Gene_ENS_ID"])
df["GWAS_variant"] = df["GWAS_variant"].str.replace("_b38", "").str.replace("_", ":")
df_ecaviar_sig_egenes = df["GWAS_variant"].value_counts().rename_axis("GWAS Loci").reset_index(name="# significant eGenes (eCAVIAR)")
df_ecaviar_sig_egenes.head()

## List of colocalizing eGenes
df = ecaviar_df_sig_eqtl[["GWAS_variant", "Gene_symbol"]]
df = df.drop_duplicates(["GWAS_variant", "Gene_symbol"])
df["GWAS_variant"] = df["GWAS_variant"].str.replace("_b38", "").str.replace("_", ":")
df_ecaviar_sig_egenes_list = df.groupby(["GWAS_variant"])["Gene_symbol"].apply(list).rename_axis("GWAS Loci").reset_index(name="List of colocalizing eGenes (eCAVIAR)")
df_ecaviar_sig_egenes_list["List of colocalizing eGenes (eCAVIAR)"] = df_ecaviar_sig_egenes_list["List of colocalizing eGenes (eCAVIAR)"].apply(lambda x: ",".join(x))
df_ecaviar_sig_egenes_list.head()

df = ecaviar_df_sig_eqtl[["GWAS_variant", "Tissue", "Gene_symbol"]]
df = df.drop_duplicates(["GWAS_variant", "Tissue", "Gene_symbol"])
df["GWAS_variant"] = df["GWAS_variant"].str.replace("_b38", "").str.replace("_", ":")
df_ecaviar_sig_egenes_list = df.groupby(["GWAS_variant"])["Gene_symbol"].apply(list).rename_axis("GWAS Loci").reset_index(name="List_of_colocalizing_eGenes(eCAVIAR)")
df_ecaviar_sig_egenes_list["List_of_colocalizing_eGenes(eCAVIAR)"] = df_ecaviar_sig_egenes_list["List_of_colocalizing_eGenes(eCAVIAR)"].apply(lambda x: ",".join(set(x)))
df_ecaviar_sig_egenes_tissue_list = df.groupby(["GWAS_variant"])["Tissue"].apply(list).rename_axis("GWAS Loci").reset_index(name="List_of_colocalizing_Tissues_eGenes(eCAVIAR)")
df_ecaviar_sig_egenes_tissue_list["List_of_colocalizing_Tissues_eGenes(eCAVIAR)"] = df_ecaviar_sig_egenes_tissue_list["List_of_colocalizing_Tissues_eGenes(eCAVIAR)"].apply(lambda x: ",".join(sorted(set(x))))

## Number sQTL-tissue pairs tested
df = ecaviar_df_sqtl[["loci", "tissue", "gene"]]
df = df.drop_duplicates(["loci", "tissue", "gene"])
df_ecaviar_sqtl_tissue_pairs = df["loci"].value_counts().rename_axis("GWAS Loci").reset_index(name="# sQTL-tissue pairs (eCAVIAR)")
df_ecaviar_sqtl_tissue_pairs.head()

# ## Number sGenes tested
df = ecaviar_df_sqtl[["loci", "gene"]]
df = df.drop_duplicates(["loci", "gene"])
df_ecaviar_sqtl_tested = df["loci"].value_counts().rename_axis("GWAS Loci").reset_index(name="# sGenes tested (eCAVIAR)")
df_ecaviar_sqtl_tested.head()

# ## Number significant sQTLs (sGene-tissue pairs)
df = ecaviar_df_sig_sqtl[["GWAS_variant", "Tissue", "Gene_ENS_ID"]]
df = df.drop_duplicates(["GWAS_variant", "Tissue", "Gene_ENS_ID"])
df["GWAS_variant"] = df["GWAS_variant"].str.replace("_b38", "").str.replace("_", ":")
df_ecaviar_sig_sqtl_tissue_pairs = df["GWAS_variant"].value_counts().rename_axis("GWAS Loci").reset_index(name="# significant sQTL-tissue pairs (eCAVIAR)")
df_ecaviar_sig_sqtl_tissue_pairs.head()

# ## Number significant sGenes
df = ecaviar_df_sig_sqtl[["GWAS_variant", "Gene_ENS_ID"]].copy()
df = df.drop_duplicates(["GWAS_variant", "Gene_ENS_ID"])
df["GWAS_variant"] = df["GWAS_variant"].str.replace("_b38", "").str.replace("_", ":")
df_ecaviar_sig_sgenes_count = df["GWAS_variant"].value_counts().rename_axis("GWAS Loci").reset_index(name="# significant sGenes (eCAVIAR)")

# ## List of significant sGenes
df = ecaviar_df_sig_sqtl[["GWAS_variant", "Gene_symbol"]]
df = df.drop_duplicates(["GWAS_variant", "Gene_symbol"])
df["GWAS_variant"] = df["GWAS_variant"].str.replace("_b38", "").str.replace("_", ":")
df_ecaviar_sig_sgenes_list = df.groupby(["GWAS_variant"])["Gene_symbol"].apply(list).rename_axis("GWAS Loci").reset_index(name="List of colocalizing sGenes (eCAVIAR)")
df_ecaviar_sig_sgenes_list["List of colocalizing sGenes (eCAVIAR)"] = df_ecaviar_sig_sgenes_list["List of colocalizing sGenes (eCAVIAR)"].apply(lambda x: ",".join(x))
df_ecaviar_sig_sgenes_list.head()

df = ecaviar_df_sig_sqtl[["GWAS_variant", "Tissue", "Gene_symbol"]]
df = df.drop_duplicates(["GWAS_variant", "Tissue", "Gene_symbol"])
df["GWAS_variant"] = df["GWAS_variant"].str.replace("_b38", "").str.replace("_", ":")
df_ecaviar_sig_sgenes_list = df.groupby(["GWAS_variant"])["Gene_symbol"].apply(list).rename_axis("GWAS Loci").reset_index(name="List_of_colocalizing_sGenes(eCAVIAR)")
df_ecaviar_sig_sgenes_list["List_of_colocalizing_sGenes(eCAVIAR)"] = df_ecaviar_sig_sgenes_list["List_of_colocalizing_sGenes(eCAVIAR)"].apply(lambda x: ",".join(set(x)))
df_ecaviar_sig_sgenes_tissue_list = df.groupby(["GWAS_variant"])["Tissue"].apply(list).rename_axis("GWAS Loci").reset_index(name="List_of_colocalizing_Tissues_sGenes(eCAVIAR)")
df_ecaviar_sig_sgenes_tissue_list["List_of_colocalizing_Tissues_sGenes(eCAVIAR)"] = df_ecaviar_sig_sgenes_tissue_list["List_of_colocalizing_Tissues_sGenes(eCAVIAR)"].apply(lambda x: ",".join(sorted(set(x))))


# ## Number significant sGenes by intron cluster
df = ecaviar_df_sig_sqtl[["GWAS_variant", "Gene_ENS_ID", "Gene_ID_intron_cluster"]].copy()
df["GWAS_variant"] = df["GWAS_variant"].str.replace("_b38", "").str.replace("_", ":")
df["cluster"] = df["Gene_ID_intron_cluster"].apply(lambda x: re.sub(":clu","",x.split("_")[0]))
df = df.drop_duplicates(["GWAS_variant", "GWAS_variant", "cluster"])
df_ecaviar_sig_sgenes_intron_cluster = df["GWAS_variant"].value_counts().rename_axis("GWAS Loci").reset_index(name="# significant sGenes by intron cluster (eCAVIAR)")
df_ecaviar_sig_sgenes_intron_cluster.head()

# ## List significant sGenes by intron cluster
df = ecaviar_df_sig_sqtl[["GWAS_variant", "Gene_symbol", "Gene_ENS_ID","Gene_ID_intron_cluster"]].copy()
df["GWAS_variant"] = df["GWAS_variant"].str.replace("_b38", "").str.replace("_", ":")
df["cluster"] = df["Gene_ID_intron_cluster"].apply(lambda x: re.sub(":clu","",x.split("_")[0]))
df["cluster"] = df["Gene_symbol"]+"-"+df["cluster"]
df = df.drop_duplicates(["GWAS_variant", "cluster"])
df_ecaviar_sig_sgenes_list_intron_cluster = df.groupby(["GWAS_variant"])["cluster"].apply(list).rename_axis("GWAS Loci").reset_index(name="List of colocalizing sGenes by intron cluster (eCAVIAR)")
df_ecaviar_sig_sgenes_list_intron_cluster["List of colocalizing sGenes by intron cluster (eCAVIAR)"] = df_ecaviar_sig_sgenes_list_intron_cluster["List of colocalizing sGenes by intron cluster (eCAVIAR)"].apply(lambda x: ",".join(x))
df_ecaviar_sig_sgenes_list_intron_cluster.head()

## Number e/sGenes tested (unique union)
df = ecaviar_df_sqtl.append(ecaviar_df_eqtl)
df = df[["loci", "gene"]]
df = df.drop_duplicates(["loci", "gene"])
df_ecaviar_tested = df["loci"].value_counts().rename_axis("GWAS Loci").reset_index(name="# e/sGenes tested (eCAVIAR)")
df_ecaviar_tested.head()

# ## Number significantly colocalizing eGene and/or sGenes (unique union)
df = ecaviar_df_sig_sqtl.append(ecaviar_df_sig_eqtl)
df = df[["GWAS_variant", "Gene_symbol"]]
df["GWAS_variant"] = df["GWAS_variant"].str.replace("_b38", "").str.replace("_", ":")
df = df.drop_duplicates(["GWAS_variant", "Gene_symbol"])
df_ecaviar_sig_egenes_sgenes = df["GWAS_variant"].value_counts().rename_axis("GWAS Loci").reset_index(name="# significantly colocalizing e/sGenes (eCAVIAR)")
df_ecaviar_sig_egenes_sgenes.head()

# ## List significantly colocalizing eGene and/or sGenes (unique union)
df = ecaviar_df_sig_sqtl.append(ecaviar_df_sig_eqtl)
df = df[["GWAS_variant", "Gene_symbol", "QTL Type"]]
df["GWAS_variant"] = df["GWAS_variant"].str.replace("_b38", "").str.replace("_", ":")
df = df.drop_duplicates(["GWAS_variant", "Gene_symbol", "QTL Type"])
df["QTL Type"] = df["QTL Type"].apply(lambda x: "e" if x == "eQTL" else "s")
df_groupby = (df.groupby(["GWAS_variant", "Gene_symbol"])["QTL Type"].apply(list).
              rename_axis(["GWAS_variant", "Gene_symbol"]).
              reset_index(name="QTL Type")
             )
df_groupby["QTL Type"] = df_groupby["QTL Type"].sort_values().apply(lambda x: sorted(x))
df_groupby["QTL Type"] = df_groupby["QTL Type"].apply(lambda x: ",".join(x))
df_groupby["gene_name"] = df_groupby["Gene_symbol"]+" ("+df_groupby["QTL Type"]+")"
df_groupby = df_groupby.drop(["QTL Type"], axis=1)

df_ecaviar_sig_union_list = df_groupby.groupby(["GWAS_variant"])["gene_name"].apply(list).rename_axis("GWAS Loci").reset_index(name="List of significantly colocalizing e/sGenes (eCAVIAR)")
df_ecaviar_sig_union_list["List of significantly colocalizing e/sGenes (eCAVIAR)"] = df_ecaviar_sig_union_list["List of significantly colocalizing e/sGenes (eCAVIAR)"].apply(lambda x: ", ".join(x))
df_ecaviar_sig_union_list.head()

# ## List of significantly colocalizing genes with their Posterior probabilities
df = ecaviar_df_sig_sqtl.append(ecaviar_df_sig_eqtl)
df = df[["GWAS_variant", "Gene_symbol", "QTL Type", "CLPP"]]
df["GWAS_variant"] = df["GWAS_variant"].str.replace("_b38", "").str.replace("_", ":")
df = df.loc[df.groupby(["GWAS_variant", "Gene_symbol", "QTL Type"])["CLPP"].idxmax()]
df["QTL Type"] = df["QTL Type"].apply(lambda x: "e" if x == "eQTL" else "s")
df = df.sort_values(["CLPP"], ascending=False)
df["CLPP"] = df["CLPP"].round(2)
df["tmp"] = df["QTL Type"]+","+"max_CLPP="+df["CLPP"].astype(str)
df_groupby = (df.groupby(["GWAS_variant", "Gene_symbol"])["tmp"].apply(list).
              rename_axis(["GWAS_variant", "Gene_symbol"]).
              reset_index(name="tmp")
             )
df_groupby["tmp"] = df_groupby["tmp"].sort_values().apply(lambda x: sorted(x))


df_groupby["tmp"] = df_groupby["tmp"].apply(lambda x: ";".join(x))
df_groupby["gene_name"] = df_groupby["Gene_symbol"]+" ("+df_groupby["tmp"]+")"
df_groupby = df_groupby.drop(["tmp"], axis=1)

df_ecaviar_sig_union_list_clpp = df_groupby.groupby(["GWAS_variant"])["gene_name"].apply(list).rename_axis("GWAS Loci").reset_index(name="List of significantly colocalizing e/sGenes by CLPP (eCAVIAR)")
df_ecaviar_sig_union_list_clpp["List of significantly colocalizing e/sGenes by CLPP (eCAVIAR)"] = df_ecaviar_sig_union_list_clpp["List of significantly colocalizing e/sGenes by CLPP (eCAVIAR)"].apply(lambda x: ", ".join(x))
df_ecaviar_sig_union_list_clpp.head()

# ## List of significantly colocalizing genes with their direction of effects
def count_directions(x):
    """
    returns sum of directions
    e.g. [4+, 3-]
    """
    plus_count = 0
    negative_count = 0

    for i in x:
        if i == "+":
            plus_count += 1
        else:
            negative_count += 1

    if plus_count == 0:
        return("["+str(negative_count)+"-"+"]")
    elif negative_count == 0:
        return("["+str(plus_count)+"+"+"]")
    else:
        return("["+str(plus_count)+"+,"+str(negative_count)+"-"+"]")

df = ecaviar_df_sig[["GWAS_variant", "Gene_symbol", "QTL Type", "Tissue", "CLPP","QTL_GWAS_relative_direction_effect"]].copy()
df = df.reset_index(drop=True)
df["GWAS_variant"] = df["GWAS_variant"].str.replace("_b38", "").str.replace("_", ":")

df = df.loc[df.groupby(["GWAS_variant", "Gene_symbol", "QTL Type", "Tissue"])["CLPP"].idxmax()]

df = df.drop_duplicates(["GWAS_variant", "Gene_symbol", "QTL Type", "Tissue"])

df["new_direction"] = df["QTL_GWAS_relative_direction_effect"].apply(lambda x: "-" if x == "OPPOSITE" else "+")
df = df.drop(["QTL_GWAS_relative_direction_effect", "Tissue"], axis=1)
df["QTL Type"] = df["QTL Type"].apply(lambda x: "e" if x == "eQTL" else "s")
df = df.reset_index(drop=True)
df_groupby = (df.groupby(["GWAS_variant", "Gene_symbol", "QTL Type"])["new_direction"].apply(list).
              rename_axis(["GWAS_variant", "Gene_symbol", "QTL Type"]).
              reset_index(name="tmp")
             )
df_groupby["counts"] = df_groupby["tmp"].apply(lambda x: count_directions(x))
df_groupby["new_tmp"] = df_groupby["QTL Type"]+df_groupby["counts"]
df_groupby = df_groupby.drop(["tmp", "counts"], axis=1)

df_groupby = (df_groupby.groupby(["GWAS_variant", "Gene_symbol"])["new_tmp"].apply(list).
              rename_axis(["GWAS_variant", "Gene_symbol"]).
              reset_index(name="tmp")
             )
df_groupby["tmp"] = df_groupby["tmp"].apply(lambda x: ";".join(x))
df_groupby["gene_name"] = df_groupby["Gene_symbol"]+" ("+df_groupby["tmp"]+")"
df_groupby = df_groupby.drop(["tmp"], axis=1)

df_ecaviar_sig_union_list_effect = df_groupby.groupby(["GWAS_variant"])["gene_name"].apply(list).rename_axis("GWAS Loci").reset_index(name="List of significantly colocalizing e/sGenes by direction of effect (eCAVIAR)")
df_ecaviar_sig_union_list_effect["List of significantly colocalizing e/sGenes by direction of effect (eCAVIAR)"] = df_ecaviar_sig_union_list_effect["List of significantly colocalizing e/sGenes by direction of effect (eCAVIAR)"].apply(lambda x: ", ".join(x))
df_ecaviar_sig_union_list_effect.head()

# ## List of significantly colocalizing genes w/ pp and effect
def count_directions(x):
    """
    returns sum of directions
    e.g. [4+, 3-]
    """
    plus_count = 0
    negative_count = 0

    for i in x:
        if i == "+":
            plus_count += 1
        else:
            negative_count += 1

    if plus_count == 0:
        return("["+str(negative_count)+"-"+"]")
    elif negative_count == 0:
        return("["+str(plus_count)+"+"+"]")
    else:
        return("["+str(plus_count)+"+,"+str(negative_count)+"-"+"]")


## posterior probability

df_pp = ecaviar_df_sig[["GWAS_variant", "Gene_symbol", "QTL Type", "CLPP"]].copy()
df_pp["GWAS_variant"] = df_pp["GWAS_variant"].str.replace("_b38", "").str.replace("_", ":")
df_pp = df_pp.reset_index(drop=True)
df_pp = df_pp.loc[df_pp.groupby(["GWAS_variant", "Gene_symbol", "QTL Type"])["CLPP"].idxmax()]
df_pp["QTL Type"] = df_pp["QTL Type"].apply(lambda x: "e" if x == "eQTL" else "s")
df_pp = df_pp.sort_values(["CLPP"], ascending=False)
df_pp["CLPP"] = df_pp["CLPP"].round(2)
df_pp["tmp_pp"] = "max_CLPP="+df_pp["CLPP"].astype(str)
df_pp = df_pp[["GWAS_variant", "Gene_symbol", "QTL Type","tmp_pp"]]

## effect direction

df_effect = ecaviar_df_sig[["GWAS_variant", "Gene_symbol", "QTL Type", "Tissue", "CLPP","QTL_GWAS_relative_direction_effect"]].copy()
df_effect["GWAS_variant"] = df_effect["GWAS_variant"].str.replace("_b38", "").str.replace("_", ":")
df_effect = df_effect.reset_index(drop=True)
df_effect = df_effect.loc[df_effect.groupby(["GWAS_variant", "Gene_symbol", "QTL Type", "Tissue"])["CLPP"].idxmax()]
df_effect = df_effect.drop_duplicates(["GWAS_variant", "Gene_symbol", "QTL Type", "Tissue"])
df_effect["new_direction"] = df_effect["QTL_GWAS_relative_direction_effect"].apply(lambda x: "-" if x == "OPPOSITE" else "+")
df_effect = df_effect.drop(["QTL_GWAS_relative_direction_effect", "Tissue"], axis=1)

df_effect["QTL Type"] = df_effect["QTL Type"].apply(lambda x: "e" if x == "eQTL" else "s")
df_effect = df_effect.reset_index(drop=True)

df_effect_groupby = (df_effect.groupby(["GWAS_variant", "Gene_symbol", "QTL Type"])["new_direction"].apply(list).
                     rename_axis(["GWAS_variant", "Gene_symbol", "QTL Type"]).
                     reset_index(name="tmp")
)
df_effect_groupby["counts"] = df_effect_groupby["tmp"].apply(lambda x: count_directions(x))
df_effect_groupby["new_tmp"] = df_effect_groupby["QTL Type"]+df_effect_groupby["counts"]
df_effect_groupby = df_effect_groupby.drop(["tmp", "counts"], axis=1)

df_pp_effect = pd.merge(df_pp, df_effect_groupby, on=["GWAS_variant", "Gene_symbol", "QTL Type"])
df_pp_effect = df_pp_effect.drop_duplicates(["GWAS_variant", "Gene_symbol", "QTL Type", "tmp_pp", "new_tmp"])

df_pp_effect["last_tmp"] = df_pp_effect["new_tmp"]+","+df_pp_effect["tmp_pp"]
df_pp_effect = df_pp_effect.drop(["tmp_pp", "new_tmp", "QTL Type"], axis=1)

df_pp_effect_groupby = (df_pp_effect.groupby(["GWAS_variant", "Gene_symbol"])["last_tmp"].apply(list).
                        rename_axis(["GWAS_variant", "Gene_symbol"]).
                        reset_index(name="last_tmp")
)
df_pp_effect_groupby["last_tmp"] = df_pp_effect_groupby["last_tmp"].apply(lambda x: ";".join(x))
df_pp_effect_groupby["gene_name"] = df_pp_effect_groupby["Gene_symbol"]+" ("+df_pp_effect_groupby["last_tmp"]+")"
df_pp_effect_groupby = df_pp_effect_groupby.drop(["last_tmp"], axis=1)

df_ecaviar_sig_union_list_clpp_effect = df_pp_effect_groupby.groupby(["GWAS_variant"])["gene_name"].apply(list).rename_axis("GWAS Loci").reset_index(name="List of significantly colocalizing e/sGenes by CLPP and direction of effect (eCAVIAR)")
df_ecaviar_sig_union_list_clpp_effect["List of significantly colocalizing e/sGenes by CLPP and direction of effect (eCAVIAR)"] = df_ecaviar_sig_union_list_clpp_effect["List of significantly colocalizing e/sGenes by CLPP and direction of effect (eCAVIAR)"].apply(lambda x: ", ".join(x))
df_ecaviar_sig_union_list_clpp_effect.head()

df_ecaviar_sig_egenes_tissue_list = df_ecaviar_sig_egenes_tissue_list.rename(columns={"GWAS_Loci":"GWAS Loci" })
df_ecaviar_sig_sgenes_tissue_list = df_ecaviar_sig_sgenes_tissue_list.rename(columns={"GWAS_Loci":"GWAS Loci" })

# ## Merge
ecaviar_table_2_df = (
    df_ecaviar_eqtl_tissue_pairs.
    merge(df_ecaviar_eqtl_tested, on=["GWAS Loci"], how="outer").
    merge(df_ecaviar_sig_eqtl_tissue_pairs, on=["GWAS Loci"], how="outer").
    merge(df_ecaviar_sig_egenes, on=["GWAS Loci"], how="outer").
    merge(df_ecaviar_sig_egenes_list, on=["GWAS Loci"], how="outer").
    merge(df_ecaviar_sqtl_tissue_pairs, on=["GWAS Loci"], how="outer").
    merge(df_ecaviar_sqtl_tested, on=["GWAS Loci"], how="outer").
    merge(df_ecaviar_sig_sqtl_tissue_pairs, on=["GWAS Loci"], how="outer").
    merge(df_ecaviar_sig_sgenes_count, on=["GWAS Loci"], how="outer").
    merge(df_ecaviar_sig_sgenes_list, on=["GWAS Loci"], how="outer").
    merge(df_ecaviar_sig_sgenes_intron_cluster, on=["GWAS Loci"], how="outer").
    merge(df_ecaviar_sig_sgenes_list_intron_cluster, on=["GWAS Loci"], how="outer").
    merge(df_ecaviar_sig_egenes_tissue_list, on=["GWAS Loci"], how="outer").
    merge(df_ecaviar_sig_sgenes_tissue_list, on=["GWAS Loci"], how="outer").
    merge(df_ecaviar_tested, on=["GWAS Loci"], how="outer").
    merge(df_ecaviar_sig_egenes_sgenes, on=["GWAS Loci"], how="outer").
    merge(df_ecaviar_sig_union_list, on=["GWAS Loci"], how="outer").
    merge(df_ecaviar_sig_union_list_clpp, on=["GWAS Loci"], how="outer").
    merge(df_ecaviar_sig_union_list_effect, on=["GWAS Loci"], how="outer").
    merge(df_ecaviar_sig_union_list_clpp_effect, on=["GWAS Loci"], how="outer")
)
ecaviar_table_2_df.head()

for col in ecaviar_table_2_df.columns:
    if col.startswith("#"):
        ecaviar_table_2_df[col] = ecaviar_table_2_df[col].fillna(0).astype(int)
    if col.startswith("List"):
        ecaviar_table_2_df[col] = ecaviar_table_2_df[col].fillna('NA')
ecaviar_table_2_df.head()

table_2_sorted = ecaviar_table_2_df
table_2_sorted = table_2_sorted[table_2_sorted["GWAS Loci"].isin(eye_disease_df["GWAS Loci"])]
table_2_sorted.columns

table_2_sorted = table_2_sorted.rename(columns = {"BP_Build38_x": "BP_Build38"})
table_2_sorted = table_2_sorted.loc[:, ~table_2_sorted.columns.str.endswith('_y')]
print(output_filename)
table_2_sorted.to_csv(output_filename,
               header=True,
               sep="\t",
               index=None, 
               na_rep="NA"
              )

