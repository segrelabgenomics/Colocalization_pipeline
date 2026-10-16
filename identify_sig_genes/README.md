# Identify significant genes within the LD window of the lead GWAS variant

The following README contains directions for identifying significant genes within the LD window of the lead GWAS variant. 

## Dependencies
* The environment requires Python 3.8+

## Running Environment and Run Time

It takes few hours but small amount of memory to identify significant genes with the LD window.

## Repository structure

The base_dir/identify_sig_genes/ folder consists of scripts to identify significant genes with the LD window.

```
run_significant_genes.sh
```

# Instructions:

Below are the following arguments. Use this shell script:

```
run_significant_genes.sh
```

* ``--gwas_variant``: GWAS variant of interest.
* ``--gwas_trait``: GWAS trait
* ``--ld_interval_directory``: LD interval directory
* ``--qtl_type``: QTL Type
* ``--sig_genes_directory``: Directory to place significant genes
* ``--base_dir_to_gtex_qtls``: Directory to place the GTEx e/sQTLs
* ``--min_var_in_ld_window``: Minimum number of variants present to extract significant genes, Default 5

#Additional inputs:

* ``loci.txt``: File contains list of genome-wide significant variants in one column in 1000G format without any header

Files in the base directory need to be split by sub-directory for QTL type
GTEx_Analysis_v8_eQTL - File in this directory have the suffix: ".v8.signif_variant_gene_pairs.txt.gz"
GTEx_Analysis_v8_sQTL -  File in this directory have the suffix: ".v8.sqtl_signifpairs.txt.gz""

