#!/bin/bash
#SBATCH --partition=long
#SBATCH --job-name="S2_eCAV"
#SBATCH -e logs/%j.err
#SBATCH -o logs/%j.out
#SBATCH --export=ALL

export base_dir="$PWD"
trait="POAG_EUR"
#Path to results table
table=$base_dir"/generate_tables/tables/"$trait"_eCAVIAR_Results_Table.txt"
#Path to intermediate files

output=$base_dir"/generate_tables/tables/Intermediate_files/"$trait/

mkdir -p "$output"

while IFS=';' read -r tissue qtl gene
do 
  python3 curate_ecaviar_table_pretable1.py --ecaviar_table "$table" \
  --tissue "$tissue" \
  --qtl_type "$qtl" \
  --gene "$gene" \
  --gwas_trait "$trait" \
  --output_directory "$output"
done < $base_dir"generate_tables/tables/MANIFEST_"$trait"_slope_se.txt"