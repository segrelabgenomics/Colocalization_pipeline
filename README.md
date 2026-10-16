# eCAVIAR Colocalization Pipeline

**Segre lab, Massachusetts Eye and Ear, Harvard Medical School**

This repository contains the scripts and instructions for running a colocalization pipeline that runs eCAVIAR, a Bayesian-based method that estimates the probability that two co-occurring association signals, such as a genome-wide association study (GWAS) locus and an expression quantitative trait locus (eQTL) are tagging the same causal variant or haplotype, while allowing for more than one causal variant in a locus as specified by the user. The pipeline includes pre-processing of input GWAS and QTL files into eCAVIAR and post-processing of output files and summary tables. The eCAVIAR tool was written by Farhad Hormozdiari. Additional details in the eCAVIAR method can be found in the following publication:

> Hormozdiari F, van de Bunt M, Segrè AV, Li X, Joo JWJ, Bilow M, Sul JH, Sankararaman S, Pasaniuc B, Eskin E. **Colocalization of GWAS and eQTL Signals Detects Target Genes.** *Am J Hum Genet.* 2016;99(6):1245–1260. [doi:10.1016/j.ajhg.2016.10.003](https://doi.org/10.1016/j.ajhg.2016.10.003)

**eCAVIAR source code:** https://github.com/fhormoz/caviar

---

**Authors:**
Andrew Hamel<sup>1,2,3</sup>, Puja A. Mehta<sup>1,2,3</sup>, Yuyang (Leo) Luo<sup>1,2,3</sup>, Ayellet V. Segre<sup>1,2,3</sup>
Segre lab, Massachusetts Eye and Ear, Harvard Medical School, Boston, MA

**Affiliations:**
<sup>1</sup> Ocular Genomics Institute, Massachusetts Eye and Ear, Boston, MA  
<sup>2</sup> Department of Ophthalmology, Harvard Medical School, Boston, MA  
<sup>3</sup> Broad Institute of MIT and Harvard, Cambridge, MA


For questions or comments regarding this pipeline, please contact Puja Mehta at `pamehta [at] meei [dot] harvard [dot] edu` or Ayellet Segre at `ayellet_segre [at] meei [dot] harvard [dot] edu`.

---


## Installation

### Python

- Python 3.9 or newer
- Python packages listed in `requirements.txt`

The scripts use only two third-party Python packages:

```text
numpy
pandas
```

All other Python imports are from the standard library and do not require separate installation.

Install Python dependencies:

```bash
python3 -m pip install -r requirements.txt
```

Or create a clean virtual environment first:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install -r requirements.txt
```

### External Tools

These are not installed by `pip` and must be available in your shell:

```text
PLINK 1.9  compute LD matrices and LD intervals
eCAVIAR    run colocalization; executable name must be eCAVIAR
```

Install PLINK. One option is Conda/Bioconda:

```bash
conda install -c bioconda plink
```

Install eCAVIAR from source:

```bash
git clone https://github.com/fhormoz/caviar.git
cd caviar/CAVIAR-C++
make
```

Then add the built eCAVIAR executable to your `PATH`, or create a symlink named `eCAVIAR` in a directory already on `PATH`:

```bash
export PATH="/path/to/caviar/CAVIAR-C++:$PATH"
which eCAVIAR
```

If running the provided cluster wrappers, you also need:

```text
SLURM              
sbatch job submission

Environment Modules module load commands used on HPC systems
```

On a non-HPC machine, remove or adapt the `module load` lines in the shell/Python command strings.

## Repository Layout

```text
ld_interval/          Compute LD intervals around lead GWAS variants
map_gwas/             Align GWAS variants to a reference variant list
identify_sig_genes/   Identify QTL genes overlapping each LD interval
manifest_ecaviar/     Build eCAVIAR manifest files
manifest_files/       Manifest files for running eCAVIAR
eCAVIAR/              Run eCAVIAR for manifest rows
generate_tables/      Parse eCAVIAR outputs and build Table 1 and Table 2
data/                 Local input data
outputs/              Generated outputs
```

## Configuration

Wrappers use repo-relative defaults. Override paths with environment variables:
```bash
export PROJECT_DIR="$PWD"
export GWAS_TRAIT="TRAIT_NAME"   # e.g., "POAG"
export GWAS_FILE="data/gwas/gwas_summary.tsv"
export PLINK_DIR="data/1000Genomes/Plink_Files"
export PLINK_SUBSET="data/1000Genomes/sample_subset.txt"
export SECOND_PLINK_DIR="data/GTEx/genotypes/PLINK_files"
export SECOND_PLINK_SUBSET="data/GTEx/genotypes/sample_subset.txt" # e.g., sample IDs of EUR samples subset
export GTEX_QTL_DIR="data/GTEx/V8_release/eQTLs"
```

Users can also modify the pipeline to use their choice of reference panel.

Optional:

```bash
export VIRTUAL_ENV_ACTIVATE=".venv/bin/activate"
```

## Input Formats

### Lead GWAS Loci

Single column, no header:

```text
chrN:POSITION:REF:ALT
chrN:POSITION:REF:ALT
```

### GWAS Summary Statistics

Raw GWAS files can use custom column names; map them with `map_gwas.py` arguments. After alignment, the eCAVIAR step expects a tab-delimited file with these code-used columns:

- `CHR`: chromosome, with or without `chr` prefix
- `pos`: position
- `variant`: reference-panel variant ID; this is renamed internally to `variant_id`
- `zscore`: GWAS z-score
- `P`, `P-value`, `mtag_pval`, or `gwas_p_value`: GWAS p-value

```text
CHR pos variant zscore P
chrN POSITION chrN:POSITION:REF:ALT ZSCORE PVALUE
chrN POSITION chrN:POSITION:REF:ALT ZSCORE PVALUE
```

### Reference Variant List

One variant ID per row with a header:

```text
variant
chrN:POSITION:REF:ALT
chrN:POSITION:REF:ALT
```

For GTEx-style variant IDs:

```text
variant
chrN_POSITION_REF_ALT_b38
chrN_POSITION_REF_ALT_b38
```

### PLINK Sample Subset

No header. Use the format by PLINK `--keep` to include ancestry specific samples or desired set of samples.

```text
FID1 IID1
FID2 IID2
```

If your PLINK files use a single sample ID, duplicate it in both columns.

### Significant QTL Pairs

Significant eQTL/sQTL variant-gene pair files in GTEx-format are used to identify candidate genes in each LD interval. The code requires `variant_id` and `gene_id`; additional GTEx summary-statistic columns are allowed. If a GTEx-like file uses another feature column name, such as `phenotype_id`, user will need to rename it to `gene_id` before running this step.

```text
variant_id gene_id tss_distance ma_samples ma_count maf pval_nominal slope slope_se
chrN_POSITION_REF_ALT_b38 ENSG_ID DIST N COUNT MAF PVALUE BETA SE
chrN_POSITION_REF_ALT_b38 ENSG_ID DIST N COUNT MAF PVALUE BETA SE
```

### All QTL Associations

The All-associations files per tissue in GTEx-format provide the per variant-gene summary statistics for variant-gene pairs tested that is used to create eCAVIAR QTL z-scores. These files must be gzipped, tab-delimited, and include a header row. The eCAVIAR code reads the header from the compressed file, and uses the following columns:

- `gene_id`
- `variant_id`
- `pval_nominal`
- `slope`
- `slope_se`

The QTL z-score is computed as `slope / slope_se`. The GTEx all-associations format:

```text
gene_id variant_id tss_distance ma_samples ma_count maf pval_nominal slope slope_se
ENSG_ID chrN_POSITION_REF_ALT_b38 DIST N COUNT MAF PVALUE BETA SE
ENSG_ID chrN_POSITION_REF_ALT_b38 DIST N COUNT MAF PVALUE BETA SE
```

### eCAVIAR Z-score Files

The pipeline writes these intermediate files without headers because eCAVIAR expects two columns only:

```text
VARIANT_ID ZSCORE
N1:POSITION1:REF1:ALT1 ZSCORE1
N2:POSITION2:REF2:ALT2 ZSCORE2
```

### eCAVIAR Manifest file

This file contains on each row the following information for each tissue and QTL dataset to be analyzed with eCAVIAR:


Semicolon-delimited, no header:

```text
gwas_loci;gwas_file;qtl;gwas_trait;tissue;gene;plink_directory;plink_subset
chrN:POSITION:REF:ALT;data/gwas/aligned.tsv;eQTL;TRAIT;TISSUE;ENSG_ID;data/ref/plink;data/ref/samples.txt
```

## Workflow

1. Compute LD intervals:

```bash
cd ld_interval

sbatch run_ld_interval.sh chrN:POSITION:REF:ALT
```

2. Align GWAS variants to the reference panel:

```bash
cd map_gwas
bash run_map_gwas.sh
```

3. Identify genes with significant QTL variants that fall in LD intervals:

```bash
cd identify_sig_genes
sbatch run_significant_genes.sh
```

4. Build eCAVIAR manifest file:

```bash
cd manifest_ecaviar
sbatch run_curate_manifest_files.sh
```

5. Run eCAVIAR:

```bash
cd eCAVIAR
bash run_colocalization_ecaviar.sh
```

6. Generate tables:

```bash
cd generate_tables
            bash run_ecaviar_results_table.sh
            Rscript Manifest_slope_se.R
            bash run_curate_ecaviar_table_pretable1.sh
            bash run_concatenate_output_tables.sh
            bash run_curate_table1.sh
            python3 finish_table_1_parsing.py
            bash run_add_significant_column.sh
            python3 change_qc_cutoff_table1.py
            python3 change_final_header_table1.py
            python3 finish_table_2_parsing.py
```

## Notes
- Each directory includes a README describing its step of the analysis.
- Keep large datasets in `data/`; keep generated results in `outputs/`.
- The scripts assume matching genome builds between the GWAS, QTL, and reference panel input datasets.
- Before running on a new system, confirm that `plink`, `eCAVIAR`, and Python are available in your shell or module environment.

## License

Our code is distributed under the terms of the BSD 3-Clause License. See LICENSE.txt file for more details.

## Citations
1. Hamel et al., "Integrating genetic regulation and single-cell expression with GWAS prioritizes causal genes and cell types for glaucoma", Nature Communications
Jan 9;15(1):396 (2024). doi: 10.1038/s41467-023-44380-y.
2. GTEx Consortium, "The GTEx Consortium atlas of genetic regulatory effects across human tissues", Science 369, 1318–1330 (2020).
3. Purcell S, et al. PLINK: A tool set for whole-genome association and population-based linkage analyses. American Journal of Human Genetics 81, 559–575 (2007). https://doi.org/10.1086/519795
4. Hormozdiari F, et al. Colocalization of GWAS and eQTL Signals Detects Target Genes. American Journal of Human Genetics 99, 1245–1260 (2016). https://doi.org/10.1016/j.ajhg.2016.10.003
5. Byrska-Bishop M, et al. High-coverage whole-genome sequencing of the expanded 1000 Genomes Project cohort including 602 trios. Cell 185, 3426–3440 (2022). https://doi.org/10.1016/j.cell.2022.08.004