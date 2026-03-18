#!/bin/bash
#SBATCH --job-name=admix
#SBATCH -A fnrdewoody
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 14-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL

module load biocontainers
module load angsd

BAM_LIST="bam.list"

#run angsd to get beagle file for admixture
angsd -bam "$BAM_LIST" -GL 1 -doGlf 2 -out MONQ_beagle -doMaf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doIBS 2 -doCounts 1 -doCov 1 -doHWE 1

mkdir -p ADX

reference_genome="/scratch/negishi/allen715/MONQ/monq_ref/ref.fa"
fai_file="/scratch/negishi/allen715/MONQ/monq_ref/ref.fa.fai"
beagle="/scratch/negishi/allen715/MONQ/shortread_analysis/MONQ_beagle.beagle.gz" 

cd ADX

# loop over K values
for K in 1 2 3 4 5 6 7 8 9 10; do
  for i in {1..5}; do
    output_file="admix_K${K}_run${i}"
    NGSadmix -likes "$beagle" -K $K -o $output_file -P 64
  done
done
