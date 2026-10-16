# Program: bcftools (Tools for variant calling and manipulating VCFs and BCFs)
# License: GNU GPLv3+, due to use of the GNU Scientific Library
# Version: 1.11 (using htslib 1.11)


#From VCF
bcftools  query -f '%CHROM\_%POS\_%REF\_%ALT\_b38\n' input.vcf.gz > variants.txt
cat variants.txt | sort | uniq > VCF_unique_variants.txt
(echo "variant"; cat VCF_unique_variants.txt) > VCF_unique_variants_w_header.txt

#From bim files
awk 'BEGIN{OFS=""} {print $2}' 1000_genomes_hg38.bim | sort | uniq > 1000G_unique_variants.txt
(echo "variant"; cat 1000G_unique_variants.txt) > 1000G_unique_variants_w_header.txt