# Colocalization analysis using eCAVIAR

The following README contains directions for performing colocalization using eCAVIAR. 

## Dependencies

Computing the LD interval requires Python 3.8+ and PLINK. 


## Running Environment and Run Time

Reserving 25GB is good for most loci but can take upto 150GB for the HLA region on chromosome 6.


## Repository structure

The base_dir/eCAVIAR/ folder consists of scripts to perform colocalization.

```
run_colocalization_ecaviar.sh
```

Below are the following arguments. Use this shell script:


* ``--gwas_variant``: GWAS variant of interest.
* ``--gwas_file``: GWAS file
* ``--qtl_type``: QTL type
* ``--gwas_trait``: GWAS trait of interest
* ``--tissue``: Tissue
* ``--gene``: Gene or Splice event in case of sQTLs
* ``--plink_directory``: Plink directory.
* ``--plink_subset``:Plink subest.
* ``--output_directory``: Output directory
* ``--ld_directory``: Directory with LD intervals of the lead GWAS variant.
