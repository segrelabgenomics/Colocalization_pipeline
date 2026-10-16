import os
import os.path
import subprocess
import glob
import pandas as pd
import re
import sys

def retrieve_bfile_gtex(chrom, plink_directory):
    """
    identify if plink directory contains /
    add kgp.clean. + chromosome
    """
    if plink_directory[-1] != "/":
        plink_directory = plink_directory + "/"
    return(plink_directory+"GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01."+chrom)

def retrieve_bfile_kgenomes(chrom,plink_directory):
    """
    identify if plink directory contains /
    add kgp.clean. + chromosome
    """
    if plink_directory[-1] != "/":
        plink_directory = plink_directory + "/"
    return(plink_directory+"kgp.clean."+chrom)

def create_plink_out_name(output_directory,trait_name,variant):
    if output_directory[-1] != "/":
        output_directory += "/"
    return(output_directory+trait_name+"_"+variant)

def switch_alleles(variant):
    """
    switches variant alleles
    e.g. 1:1234:A:G -> 1:1234:G:A
    """
    variant_split = variant.split("_")
    return(variant_split[0]+"_"+variant_split[1]+"_"+variant_split[3]+"_"+variant_split[2]+"_b38")

def switch_alleles_kgenomes(variant):
    """
    switches variant alleles
    e.g. 1:1234:A:G -> 1:1234:G:A
    """
    variant_split = variant.split(":")
    return(variant_split[0]+":"+variant_split[1]+":"+variant_split[3]+":"+variant_split[2])

def run_plink(output_directory, chrom, plink_directory, plink_subset, trait_name,variant,memory="10000",r2_cutoff="0.5",ld_window="99999",ld_window_kb="1000"):
    """
    run plink to determine ld interval
    """
    bfile = retrieve_bfile_gtex(chrom, plink_directory)
    out_name = create_plink_out_name(output_directory,trait_name,variant)

    plink_command_list = ["module load module-autoload/1.0 && module load PLINK && plink",
                    "--bfile",bfile,"--memory",memory,
                    "--keep",plink_subset,"--ld-snp",variant,
                    "--r2","--ld-window-r2",r2_cutoff,
                    "--ld-window",ld_window,
                    "--ld-window-kb",ld_window_kb,"--out",out_name]
    plink_command = " ".join(plink_command_list)+ " && sleep 5 && sync" 
    subprocess.call(plink_command, shell = True)

    return(out_name)

def run_plink_kgenomes(output_directory,chrom, k_plink_directory, k_subset, trait_name,variant,memory="10000",r2_cutoff="0.1",ld_window="99999",ld_window_kb="1000"):
    """
    run plink to determine ld interval
    """

    bfile = retrieve_bfile_kgenomes(chrom, k_plink_directory)
    out_name = create_plink_out_name(output_directory,trait_name,variant)
    #variant = variant.replace("chr","") 

    plink_command_list = ["module load module-autoload/1.0 && module load PLINK && plink",
                    "--bfile",bfile,"--memory",memory,
                    "--keep",k_subset,"--ld-snp",variant,
                    "--r2","--ld-window-r2",r2_cutoff,
                    "--ld-window",ld_window,
                    "--ld-window-kb",ld_window_kb,"--out",out_name]
    plink_command = " ".join(plink_command_list)+ " && sleep 5 && sync" 
    subprocess.call(plink_command, shell = True)

    return(out_name)

def finish_ld_parsing(df, variant):
    """
    get ld intervals
    """
    print("parsing ld file")
    df["pos"] = df["SNP_B"].apply(lambda x: x.split(":")[1]).astype(int)
    df = df.rename(columns={"SNP_B":"variant"})

    try:
        start = df.loc[df["pos"].idxmin()]["pos"]
        end = df.loc[df["pos"].idxmax()]["pos"]
    except ValueError:
        variant_pos = int(variant.split(":")[1])
        start = variant_pos - 250000
        end = variant_pos + 250000
        return(start, end)
    start = int(start) - 50000
    end = int(end) + 50000

    return(start,end)

def find_nearest_variant(variant, chrom, plink_directory):
    """
    finds nearest variant by position
    """
    k1_bim = retrieve_bfile_kgenomes(chrom, plink_directory)
    k1_bim_df = pd.read_csv(k1_bim + ".bim", sep = "\t", header = None)
    k1_bim_df = k1_bim_df[[0,1,3,4,5]]
    k1_bim_df.columns = ["chr", "variant_id", "pos", "alt", "ref"]
    pos = int(variant.split(":")[1])

    k1_bim_df["delta"] = k1_bim_df["pos"] - pos
    k1_bim_df["abs"] = k1_bim_df["delta"].abs()
    nearest_variant = k1_bim_df.loc[k1_bim_df["abs"].idxmin()]["variant_id"]
    return(nearest_variant)

def finish_ld_parsing_gtex(df, original_variant, gtex_variant, chrom, plink_directory):
    """
    get ld intervals
    """
    df["pos"] = df["SNP_B"].apply(lambda x: x.split("_")[1]).astype(int)
    df = df.rename(columns={"SNP_B":"variant","CHR_B":"chr"})

    df = df[["chr","pos","R2"]]
    df["chr"] = "chr"+df["chr"].astype(str)

    k1_variant_tmp = re.sub("_b38","",gtex_variant)
    k1_variant_tmp = re.sub("_",":",k1_variant_tmp)
    k1_bim = retrieve_bfile_kgenomes(chrom, plink_directory)
    k1_bim_df = pd.read_csv(k1_bim+".bim",sep="\t",header=None)
    k1_bim_df = k1_bim_df[[0,1,3,4,5]]
    k1_bim_df.columns = ["chr","variant_id","pos","alt","ref"]
    k1_bim_df["chr"] = "chr" + k1_bim_df["chr"].astype(str)

    tmp = pd.merge(k1_bim_df,df,on=["chr","pos"])
    tmp = tmp.reset_index(drop=True)

    try:
        ld_proxy = tmp.loc[tmp["R2"].idxmax()]["variant_id"]
        return(ld_proxy)
    except ValueError:
        nearest_variant = find_nearest_variant(original_variant, chrom, plink_directory)
        return(nearest_variant)

def compute_ld_interval(output_directory,chrom, variant, trait_name, plink_directory, plink_subset, second_plink_directory, second_plink_subset):
    """
    find ld interval at r2 cutoff
    """
    out_name = run_plink_kgenomes(output_directory, chrom, plink_directory, plink_subset, trait_name, variant)
    try:
        df = pd.read_csv(out_name+".ld",sep="\s+")
        start,end = finish_ld_parsing(df, variant)
        return(start,end,variant)
    except FileNotFoundError:
        print("\nswitching alleles\n")
        new_variant = switch_alleles_kgenomes(variant)
        print("\nrestarting plink out name\n")
        new_out_name = run_plink_kgenomes(output_directory,chrom, plink_directory, plink_subset, trait_name,new_variant)
        try:
            df = pd.read_csv(new_out_name+".ld",sep="\s+")
            start,end = finish_ld_parsing(df, variant)
            return(start,end,new_variant)
        except FileNotFoundError:
            print("\nnot found with switching alleles. resorting to GTEx\n")
            gtex_variant = re.sub(":","_",variant)+"_b38"
            new_out_name = run_plink(output_directory,chrom, second_plink_directory, second_plink_subset, trait_name,gtex_variant)
            try:
                df = pd.read_csv(new_out_name+".ld",sep="\s+")
                ld_proxy_variant = finish_ld_parsing_gtex(df,variant, gtex_variant, chrom, plink_directory)
                proxy_out_name = run_plink_kgenomes(output_directory, chrom, plink_directory, plink_subset, trait_name,ld_proxy_variant)
                df = pd.read_csv(proxy_out_name+".ld",sep="\s+")
                start,end = finish_ld_parsing(df, variant)

                return(start,end,ld_proxy_variant)
            except FileNotFoundError:
                print("\nswitching alleles again\n")
                new_gtex_variant = switch_alleles(gtex_variant)
                print("\n",new_gtex_variant,"\n")
                new_out_name = run_plink(output_directory, chrom, second_plink_directory, second_plink_subset, trait_name, new_gtex_variant)
                try:
                    df = pd.read_csv(new_out_name+".ld",sep="\s+")
                    ld_proxy_variant = finish_ld_parsing_gtex(df,variant, new_gtex_variant, chrom, plink_directory)

                    proxy_out_name = run_plink_kgenomes(output_directory, chrom, second_plink_directory, second_plink_subset, trait_name, ld_proxy_variant)
                    df = pd.read_csv(proxy_out_name+".ld",sep="\s+")
                    start,end = finish_ld_parsing(df, variant)
                    return(start,end,ld_proxy_variant)
                except FileNotFoundError:
                    print("\n variant not found in anything \n")
                    #add module that takes nearest variant
                    nearest_variant = find_nearest_variant(variant, chrom, plink_directory)
                    
                    final_out_name = run_plink_kgenomes(output_directory, chrom, plink_directory, plink_subset, trait_name, nearest_variant)
                    df = pd.read_csv(final_out_name+".ld",sep="\s+")
                    start,end = finish_ld_parsing(df, variant)
                    return(start,end,nearest_variant)




