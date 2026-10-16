# Curate manifest files 

The following README contains directions for generating manifest files to launch parallel jobs for eCAVIAR. 

## Dependencies

This requires Python 3.8+. 


## Running Environment and Run Time

One can also run it locally.

## Repository structure

The base_dir/manifest_ecaviar/ folder consists of scripts to generate the manifest files.

```
run_curate_manifest_files.sh
```

Below are the following arguments. Use this shell script:

```
run_curate_manifest_files.sh
```
* ``--gwas_file``: Path to GWAS file
* ``--sig_genes_directory``: Directory containing significant genes
* ``--gwas_trait``: GWAS trait of interest, e.g. POAG\_EUR.
* ``--plink_directory``: Plink directoy
* ``--plink_subset``: Subset file for plink directory
* ``--manifest_dir``: Directory to manifest file
* ``--manifest_chunks_dir``: Directory to chunks of 1000 runs for manifest file


#Example of a command from the manifest file:

chr6:134051012:G:C;repo_dir/map_gwas/gwas_aligned_UPD_1KG_NYGC_B38.txt;eqtl;POAG_EUR;Brain_Cerebellar_Hemisphere;ENSG00000028839.10;base_dir/data/1000Genomes/NYGC_HG38/Plink_Files/;base_dir/data/1000Genomes/NYGC_HG38_unrelated/sample_subset.txt



