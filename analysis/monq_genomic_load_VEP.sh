#!/bin/bash
#SBATCH -A fnrdewoody
#SBATCH --job-name=VEP_monq
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 5-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL

module load biocontainers
module load ensembl-vep/108.2
module load bcftools
module load htslib

RAW_GFF="jaqu_plus_chicChW_renamed.gff"
SORTED_GFF="jaqu_plus_chicChW_renamed.sorted.gff.gz"

SNP_VCF="monq_mcSNP.final.sort.vcf"
SV_VCF="monq_mcSV.final.sort.vcf"
SV_OUTLIERS_VCF="monq_SVoutliers_final.vcf"

GENOME_FASTA="monq-out_refp_renamed.fa"
OUTDIR="VEP_results"
THREADS=64

mkdir -p "${OUTDIR}"
mkdir -p filtered_vcf

# prep gff
sed 's/\r$//' "$RAW_GFF" > tmp.clean.gff
grep '^#' tmp.clean.gff > header.txt
grep -v '^#' tmp.clean.gff | \
    LC_ALL=C sort -t $'\t' -k1,1 -k4,4n > body.sorted.gff

cat header.txt body.sorted.gff > jaqu_plus_chicChW_renamed.sorted.gff
bgzip -f jaqu_plus_chicChW_renamed.sorted.gff
tabix -f -p gff jaqu_plus_chicChW_renamed.sorted.gff.gz

rm -f tmp.clean.gff header.txt body.sorted.gff

# sort and index vcfs 
bcftools sort -Oz -o monq_mcSNP.sorted.vcf.gz "$SNP_VCF"
tabix -f -p vcf monq_mcSNP.sorted.vcf.gz

bcftools sort -Oz -o monq_mcSV.sorted.vcf.gz "$SV_VCF"
tabix -f -p vcf monq_mcSV.sorted.vcf.gz

bcftools sort -Oz -o monq_SVoutliers.sorted.vcf.gz "$SV_OUTLIERS_VCF"
tabix -f -p vcf monq_SVoutliers.sorted.vcf.gz

# run VEP
run_vep () {
    INPUT_VCF="$1"
    PREFIX="$2"
    vep \
        --input_file "${INPUT_VCF}" \
        --output_file "${OUTDIR}/${PREFIX}.vep.vcf" \
        --vcf \
        --force_overwrite \
        --fasta "${GENOME_FASTA}" \
        --gff "${SORTED_GFF}" \
        --everything \
        --variant_class \
        --species custom \
        --assembly custom \
        --fork ${THREADS} \
        --no_stats
}

run_vep "$SNP_VCF" "monq_SNP"
run_vep "$SV_VCF" "monq_SV"
run_vep "$SV_OUTLIERS_VCF" "monq_SVoutliers"

# filter results for deleterious variants
filter_deleterious () {
    IN="$1"
    OUT="$2"
    grep "^#" "${IN}" > "${OUT}"
    grep -v "^#" "${IN}" | grep -E "\|(HIGH|MODERATE)\|" >> "${OUT}" || true
}

filter_deleterious "${OUTDIR}/monq_SNP.vep.vcf" \
                   "${OUTDIR}/monq_SNP_deleterious.vcf"
filter_deleterious "${OUTDIR}/monq_SV.vep.vcf" \
                   "${OUTDIR}/monq_SV_deleterious.vcf"
filter_deleterious "${OUTDIR}/monq_SVoutliers.vep.vcf" \
                   "${OUTDIR}/monq_SVoutliers_deleterious.vcf"
