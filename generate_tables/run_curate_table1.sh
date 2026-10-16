#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name="S4_eCAV"
#SBATCH -e logs/%j.err
#SBATCH -o logs/%j.out
#SBATCH --export=ALL

export base_dir="$PWD"
trait="POAG_EUR"
table=$base_dir"generate_tables/tables/Intermediate_files/results_concatenated/"$trait"_pretable1.tsv"
gwas="$base_dir/data/gwas/gwas_aligned_UPD_1KG_NYGC_B38.txt"
output=$base_dir"generate_tables/tables/"$trait"_table1.tsv"

./curate_table_1.py --ecaviar_output_table $table \
  --gwas_file $gwas \
  --output_filename $output

