#!/usr/bin/env python3
import os
import os.path
import numpy as np
import pandas as pd
import gzip
import re

# ## gencode parsing
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
            gene_dict[attr['gene_id'].split('.')[0]]=(gene(row[0],
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

def extract_gene_version(gene_dict, ensembl_gene_id):
    """
    returns ensembl with decimal version
    """
    for i in gene_dict.keys():
        if gene_dict[i].gene_id == ensembl_gene_id:
            return(i)
    return("NA")

gencode_file = "gencode.v26.annotation.gtf"
gencode_dict = parse_gencode_gtf(gencode_file)
gencode_dict['ENSG00000223972'].gene_id

# ## Eye disease gwas table parsing

eye_disease_file = "../data/gwas/gwas.txt"
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
            
eye_disease_df["chr"]="chr"+eye_disease_df["CHR"].astype(str)
eye_disease_df["variant_id"] = (
    eye_disease_df["chr"].astype(str) + ":" +
    eye_disease_df["BP"].astype(str) + ":" +
    eye_disease_df["Allele2"] + ":" +
    eye_disease_df["Allele1"]
)
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

# ## LD interval parsing
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

ldint_dir = "../ld_intervals/POAG_EUR/"
ldint_files = [ldint_dir + file for file in os.listdir(ldint_dir) if ".txt" in file and "POAGMVP_AFR" in file]

ld_interval_df = pd.DataFrame()

for file in ldint_files:
    #print(file)
    df = parse_ld_interval_file(file)
    #print(df.head())
    ld_interval_df = ld_interval_df.append(df)
    
ld_interval_df = ld_interval_df[["GWAS Trait", "GWAS variant", 
                                 "GWAS variant (aligned to 1KG)",
                                 "LD interval start (b38)", "LD interval end (b38)",
                                "LD interval size (kb)"
                                ]]
ld_interval_df = ld_interval_df.reset_index(drop=True)
ld_interval_df = ld_interval_df.drop_duplicates(["GWAS variant"])

ld_interval_df.head()

ld_interval_df = ld_interval_df[["GWAS Trait", "GWAS variant", 
                                 "GWAS variant (aligned to 1KG)",
                                 "LD interval start (b38)", "LD interval end (b38)",
                                "LD interval size (kb)"
                                ]]
ld_interval_df = ld_interval_df.reset_index(drop=True)
ld_interval_df = ld_interval_df.drop_duplicates(["GWAS variant"])

ld_interval_df.head()

ld_interval_df["GWAS variant"] = ld_interval_df["GWAS variant"].str.replace(":","_")
ld_interval_df["GWAS variant"] = ld_interval_df["GWAS variant"] + "_b38"

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
eye_disease_df['ABS GWAS effect'] = eye_disease_df['GWAS_effect'].astype(float).abs()

ecaviar_file = "POAG_EUR_table1.tsv" #Read table 1
ecaviar_df = pd.read_csv(ecaviar_file,sep ="\t")
ecaviar_df.columns

split_variants = ecaviar_df["SNP_ID"].str.split(':', expand=True)
chrom = split_variants[0].str.replace("chr", "")
pos = split_variants[1]
ecaviar_df["SNP_ID"] = chrom + ":" + pos

ecaviar_df = ecaviar_df.drop(["chr",
                              "pos"
                              ], axis=1)

ecaviar_df.columns
if column in ("effect_allele"):
            ecaviar_df = ecaviar_df.rename(columns = {column: "ALT"}) 
if column in ("variant_id"):
            ecaviar_df = ecaviar_df.rename(columns = {column: "qtl_variant"}) 
if column in ("Offect_allele"):
            ecaviar_df = ecaviar_df.rename(columns = {column: "REF"}) 

ecaviar_df = ecaviar_df.rename(columns={"loci": "GWAS variant",
                                        "gwas_variant_aligned": "GWAS variant (aligned to 1KG)",
                                        "num_variants_tested":"Number variants tested in defined LD interval (r2 0.1 +- 50kb)",
                                        "qtl": "QTL Type",
                                        "variant_id": "qtl_Variant_id",
                                        "SNP_ID": "rsID_qtl",
                                        "tissue": "Tissue",
                                        "slope": "QTL_slope_relative_ALT",
                                        "slope_se":"QTL_se",
                                        "Effect":"GWAS_beta_relative_ALT",
                                        "pval_nominal": "QTL_pval_nominal",
                                        "Prob_in_pCausalSet": "eCAVIAR_prob in causal set",
                                        "gwas_p_value": "QTL_gwas_p_value", 
                                        "SE": "GWAS_se",
                                        "BETA":"GWAS_beta_relative_ALT"
                                       })
ecaviar_df['GWAS variant'] = ecaviar_df['GWAS variant'].str.replace(":","_")
ecaviar_df['GWAS variant'] = ecaviar_df['GWAS variant'] + "_b38"
display(ecaviar_df)

ecaviar_df = pd.merge(ecaviar_df, eye_disease_df, on=["GWAS variant"])

variants_not_in_eye_disease = ecaviar_df.loc[~ecaviar_df["GWAS variant"].isin(eye_disease_df["GWAS variant"])]
variants_not_in_eye_disease

ecaviar_df["ABS(GWAS_effect)"] = ecaviar_df["ABS GWAS effect"].astype(float).abs()
ecaviar_df.head()


ecaviar_df = ecaviar_df.rename(columns={"Number variants tested in defined LD interval (r2 0.1 +- 50kb)":
                                       "Number variants tested in defined LD interval (r2 0.1 ± 50kb)"
                                       })

# Flip GWAS beta if effect allele != ALT
ecaviar_df["GWAS_beta_relative_ALT"] = np.where(
    ecaviar_df["gwas_effect_allele"] == ecaviar_df["ALT"],#Check column names for the effect allele in GWAS and QTL studies
    ecaviar_df["GWAS_beta_relative_ALT"],           # already relative to ALT, no flip needed
    -ecaviar_df["GWAS_beta_relative_ALT"]           # effect allele is REF -> flip sign
)


conditions = [(ecaviar_df["QTL_slope_relative_ALT"] * ecaviar_df["GWAS_beta_relative_ALT"] > 0),
              (ecaviar_df["QTL_slope_relative_ALT"] * ecaviar_df["GWAS_beta_relative_ALT"] < 0)]

choices = [ "SAME", "OPPOSITE" ]
np.select(conditions, choices, default=0)


ecaviar_df["Qtl_gwas_relative_effect"] = np.select(conditions, choices, default=0)
ecaviar_df["Qtl_gwas_relative_effect"].value_counts()

ecaviar_df.columns

ecaviar_df["QTL Type"] = ecaviar_df["QTL Type"].apply(lambda x: "eQTL" if x == "eQTL" else "sQTL")
ecaviar_df["QTL Type"].unique()

ecaviar_df["gene"] = ecaviar_df["gene_id"].apply(lambda x: x.split(":")[-1])

ecaviar_df["tmp"] = ecaviar_df["gene"].apply(lambda x: x.split(".")[0])
ecaviar_df["Gene symbol"] = ecaviar_df["tmp"].apply(lambda x: extract_gene_name(gencode_dict, x))
ecaviar_df.head()

ecaviar_df = ecaviar_df.rename(columns = {'GWAS variant_x': "GWAS variant"}) 

ecaviar_df = ecaviar_df.rename(columns = {'GWAS variant (aligned to 1KG)_x': 'GWAS variant (aligned to 1KG)'}) 

ecaviar_df = ecaviar_df[["GWAS Trait", "GWAS variant", "GWAS variant (aligned to 1KG)",
                         "Effect allele", "Other Allele",
                         "rsID_GWAS", 
                         "GWAS_OR", "GWAS_effect",
                         "ABS(GWAS_effect)",
                         "GWAS_p_value","LD interval start (b38)", "LD interval end (b38)", "LD interval size (kb)",
                         "Number variants tested in defined LD interval (r2 0.1 ± 50kb)", 
                         "QTL Type", "qtl_Variant_id", 
                         "rsID_qtl",
                         "gene_id","Gene_id_intron_cluster","Gene symbol", "Tissue", 
                         "QTL_pval_nominal",
                         "QTL_gwas_p_value",
                         "QTL_slope_relative_ALT", 
                         "QTL_se",
                         "GWAS_beta_relative_ALT",
                         "Qtl_gwas_relative_effect",
                         "eCAVIAR_prob in causal set", "CLPP"
                        ]].copy()
ecaviar_df["chr"] = ecaviar_df["GWAS variant"].apply(lambda x: re.sub("chr","",x.split("_")[0])).astype(int)
ecaviar_df["pos"] = ecaviar_df["GWAS variant"].apply(lambda x: x.split("_")[1]).astype(int)

ecaviar_df = ecaviar_df.sort_values(["chr", "pos"], ascending=[True, True])
ecaviar_df = ecaviar_df.reset_index(drop=True)
ecaviar_df = ecaviar_df.drop(["chr", "pos"], axis=1)

ecaviar_df = ecaviar_df.loc[ecaviar_df["CLPP"] >= 0.01]
ecaviar_df["Gene symbol"].unique()

ecaviar_df.to_csv("POAG_EUR_table1_formatted_without_QC.tsv", 
                  sep="\t", header=True, index=None, na_rep="NA")
