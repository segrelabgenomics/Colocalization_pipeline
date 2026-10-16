#!/bin/bash
#SBATCH --job-name=sig_gene
#SBATCH --partition=normal
#SBATCH --mem=25G
#SBATCH --cpus-per-task=4

set -euo pipefail

module purge
if [[ -n "${VIRTUAL_ENV_ACTIVATE:-}" ]]; then
  source "$VIRTUAL_ENV_ACTIVATE"
fi

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
base_dir="$repo_dir"
trait="POAG_EUR"
ld_interval="$base_dir/ld_interval/ld_intervals/$trait/"
sig_genes_dir="$base_dir/identify_sig_genes/$trait/"
qtl_type="eQTL" # or sQTL, run seperately for each QTL type to make the process parallel
base_dir_to_gtex_qtls="$base_dir/data/GTEx/V8_release/eQTLs/"
loci_file="$base_dir/ld_interval/inputs/loci.txt"

mkdir -p $sig_genes_dir

while IFS=$';' read loci
do
   echo $loci
   ./colocalization_pipeline_siggenes.py --gwas_variant "$loci" \
     --gwas_trait "$trait" \
     --ld_interval_directory "$ld_interval" \
     --qtl_type "$qtl_type" \
     --sig_genes_directory "$sig_genes_dir" \
     --base_dir_to_gtex_qtls "$base_dir_to_gtex_qtls" \
     --min_var_in_ld_window 5
     
done<"$loci_file"
