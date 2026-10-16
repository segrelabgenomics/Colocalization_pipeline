#!/bin/bash
#SBATCH --job-name="S1_eCAV"
#SBATCH -e logs/%j.err
#SBATCH -o logs/%j.out
#SBATCH --export=ALL

export base_dir="$PWD"
#Directory with eCAVIAR results
gwas_trait="POAG_EUR"
ecaviar_dir=$base_dir"outputs/"$gwas_trait"/"
#Output file path and name
output_file=$base_dir"/generate_tables/tables/"$gwas_trait"_eCAVIAR_Results_Table.txt"
#TEMP file path and name
incorrect_num_file=$base_dir"/generate_tables/tables/"$gwas_trait"_eCAVIAR_Results_temp_file.txt"

echo $gwas_trait

python3 ecaviar_results_table.py --ecaviar_directory $ecaviar_dir \
 --gwas_trait $gwas_trait \
 --output_file $output_file \
 --incorrect_num_file $incorrect_num_file
