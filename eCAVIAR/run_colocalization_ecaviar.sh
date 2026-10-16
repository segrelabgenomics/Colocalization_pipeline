#!/bin/bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
trait="POAG_EUR"
base_dir="$repo_dir"
output_directory="$base_dir/outputs/$trait/" # Each trait has a seperate output dir
ld_directory="$base_dir/ld_interval/ld_intervals/$trait"
manifest_file="$base_dir/manifest_files/chunked_manifest_files/MANIFEST_chunked_${trait}_0_0_loci_tissue_gene.txt" # Chunks are in batch of 1000 jobs. Run each chunk.
CustomJobID=0001 #Change this for each chunk for precice job tracking
mkdir -p "$output_directory"
while IFS=$';' read gwas_loci gwas_file qtl gwas_trait tissue gene plink_directory plink_subset
do
  CustomJobIDFormatted=$(printf "%04d" "$CustomJobID")
  sbatch submit_jobs_ecaviar.sh "$gwas_loci" "$gwas_file" "$qtl" "$gwas_trait" "$tissue" "$gene" "$plink_directory" "$plink_subset" "$output_directory" "$ld_directory" "$CustomJobIDFormatted"
  ((CustomJobID++))
done <"$manifest_file"
