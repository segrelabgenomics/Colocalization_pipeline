# Map GWAS variants to reference panel

The following README contains directions for generating the GWAS file with variants aligned to 1KG or GTEx reference panel. The alignmnet is only for the zscores which are required for the colocalization analysis.

## Dependencies
* The environment requires Python 3.8+
* The GWAS file should be in the following format:

## Running Environment and Run Time

Reserving 10GB is plenty and takes at most 10 minutes.

## Repository structure

The base_dir/map_gwas/ folder consists of scripts to generate the reference aligned GWAS.

```
To align GWAS to 1000 Genomes reference panel use --ref as 1KG else use GTEx:
run_map_gwas.sh

```

# Instructions on how to align the gwas to reference panel

Below are the following arguments. Use this shell script:


```
run_map_gwas_gtex.sh
```

* ``--gwas_file``: GWAS file with summary statistics.
* ``--variants_file``: Path to unique list of variants.
* ``--ref``: Name of the reference panel (eg, 1KG, GTEx, VCF)
* ``--chr_name``: chromosome 
* ``--bp_name``: variant position
* ``--build``: Build of the GWAS (eg, b38)
* ``--a1``: Effect_allele
* ``--a2``: Other_allele
* ``--beta``: variant effect size
* ``--slope``: variant effect size standard error
* ``--variant_col_in_reference``: column in the reference panel file listing the variant identifier (e.g., `variant`)

Output file name is formatted based on the reference used.
