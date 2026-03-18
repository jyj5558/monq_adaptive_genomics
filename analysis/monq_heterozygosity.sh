#!/bin/bash
#SBATCH --job-name=heterozygosity
#SBATCH -A fnrdewoody
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 14-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL

mkdir jobs_het
mkdir HET

while read -a line
do 
        echo "#!/bin/bash
#SBATCH -A highmem
#SBATCH -n 20
#SBATCH -t 05:00:00
#SBATCH --job-name=${line[0]}_het_stats
#SBATCH --error=${line[0]}_het_stats.e
#SBATCH --output=${line[0]}_het_stats.o

module load biocontainers
module load angsd

# move to the bams folder
cd /scratch/negishi/allen715/MONQ/monq_shortread/

angsd -i ${line[0]}.bam -ref /scratch/negishi/allen715/MONQ/monq_ref/ref.fa  -anc /scratch/negishi/allen715/MONQ/monq_ref/ref.fa  -dosaf 1 -minMap$

realSFS -P 10 -fold 1 /scratch/negishi/allen715/MONQ/shortread_analysis/HET/${line[0]}.saf.idx > /scratch/negishi/allen715/MONQ/shortread_analysis$

done < ./sample.list

#for i in `ls -1 *sh`; do  echo "sbatch $i" ; done > jobs ; source ./jobs

#Get individual heterozygosity, proportion of heterozygotes:
#cat ./*ml

#Use output from cat command for calculate prop heterozygote
