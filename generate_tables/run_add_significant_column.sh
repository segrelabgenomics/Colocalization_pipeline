#!/bin/bash

#SBATCH --partition=normal
#SBATCG --mem=25G
#SBATCH -e %j.err
#SBATCH -o %j.out

export base_dir="$PWD"
trait="POAG_EUR"
cutoff_file=$base_dir"generate_tables/tables/"$trait"_table1_formatted_without_QC.tsv"
outputs_dir=$base_dir"generate_tables/tables/"

python3 add_significant_column.py $cutoff_file $outputs_dir


