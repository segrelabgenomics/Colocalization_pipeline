#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name="S3_eCAV"
#SBATCH -e logs/%j.err
#SBATCH -o logs/%j.out
#SBATCH --export=ALL


export base_dir="$PWD"
trait="POAG_EUR"
./concatenate_output_tables.py --base_dir_eqtl $base_dir"generate_tables/tables/Intermediate_files/"$trait"/eqtl/" \
  --base_dir_sqtl $base_dir"generate_tables/tables/Intermediate_files/"$trait"/sqtl/" \
  --output_directory $base_dir"generate_tables/tables/Intermediate_files/results_concatenated/" \
  --output_filename $trait"_pretable1.tsv"

