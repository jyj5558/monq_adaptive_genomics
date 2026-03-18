#!/bin/bash
#SBATCH --job-name=monq_kraken
#SBATCH -A fnrpupfish
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 14-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=jeon96@purdue.edu
#SBATCH --mail-type=END,FAIL

# paths
CLEAN_SRA=/scratch/bell/jeon96/monq/short-read/cleaned
FILT_SRA=/scratch/bell/jeon96/monq/short-read/filtered
FIXED_SRA=/scratch/bell/jeon96/monq/short-read/fixed

export PATH=$HOME/app/kraken2-2.1.5:$PATH

cd ${CLEAN_SRA}
mkdir -p ../filtered
mkdir -p ../fixed

for i in `ls -lrt ${CLEAN_SRA} | tr -s ' ' | cut -d ' ' -f 9 | cut -d 'R' -f 1 | sort | uniq`; do
    id=$(echo $i | cut -d '_' -f 1)
    echo "$id started"

    k2 classify --db $HOME/app/kraken2-2.1.5/pluspf_08gb --threads 64 --unclassified-out ${FILT_SRA}/${id}_filt_R1.fq --classified-out ${FILT_SRA}/${id}_cont_R1.fq --output ${FILT_SRA}/${id}_R1_out --report ${FILT_SRA}/${id}_R1_kraken --confidence 0.3 --minimum-hit-groups 2 --use-names ${i}R1_001_val_1.fq.gz
    
    k2 classify --db $HOME/app/kraken2-2.1.5/pluspf_08gb --threads 64 --unclassified-out ${FILT_SRA}/${id}_filt_R2.fq --classified-out ${FILT_SRA}/${id}_cont_R2.fq --output ${FILT_SRA}/${id}_R2_out --report ${FILT_SRA}/${id}_R2_kraken --confidence 0.3 --minimum-hit-groups 2 --use-names ${i}R2_001_val_2.fq.gz 

    echo "Finished processing $id"
done

# Fixing paired-end read orders
module load biocontainers
module load bbmap
module load htslib

cd ${FILT_SRA}

for i in `ls -lrt ${CLEAN_SRA} | tr -s ' ' | cut -d ' ' -f 9 | cut -d 'R' -f 1 | sort | uniq`; do
    id=$(echo $i | cut -d '_' -f 1)
    echo "$id started"
    
    repair.sh in1=${id}_filt_R1.fq in2=${id}_filt_R2.fq out1=${FIXED_SRA}/${id}_fixed_R1.fq.gz out2=${FIXED_SRA}/${id}_fixed_R2.fq.gz outs=${FIXED_SRA}/${id}_fixed_s.fq.gz

    echo "Finished processing $base"
done


# NM WGS samples from NCBI
cd ../
mkdir -p ./sra/filtered/
mkdir -p ./sra/fixed/

cd ./sra/cleaned/
for i in `ls -lt ./ | grep "fq" | tr -s ' ' | cut -d " " -f 9 | cut -d "_" -f 1 | sort | uniq`; do
    id=$(echo $i | cut -d '_' -f 1)
    echo "$id started"

    k2 classify --db $HOME/app/kraken2-2.1.5/pluspf_08gb --threads 64 --unclassified-out ../filtered/${id}_filt_R1.fq --classified-out ../filtered/${id}_cont_R1.fq --output ../filtered/${id}_R1_out --report ../filtered/${id}_R1_kraken --confidence 0.3 --minimum-hit-groups 2 --use-names ${i}_1_val_1.fq
    
    k2 classify --db $HOME/app/kraken2-2.1.5/pluspf_08gb --threads 64 --unclassified-out ../filtered/${id}_filt_R2.fq --classified-out ../filtered/${id}_cont_R2.fq --output ../filtered/${id}_R2_out --report ../filtered/${id}_R2_kraken --confidence 0.3 --minimum-hit-groups 2 --use-names ${i}_2_val_2.fq 

    echo "Finished processing $id"
done

# Fixing paired-end read orders
module load biocontainers
module load bbmap
module load htslib

cd ../filtered/

for i in `ls -lrt ./ | tr -s ' ' | cut -d ' ' -f 9 | cut -d '_' -f 1 | sort | uniq`; do
    id=$(echo $i | cut -d '_' -f 1)
    echo "$id started"
    
    repair.sh in1=${id}_filt_R1.fq in2=${id}_filt_R2.fq out1=../fixed/${id}_fixed_R1.fq.gz out2=../fixed/${id}_fixed_R2.fq.gz outs=../fixed/${id}_fixed_s.fq.gz

    echo "Finished processing $base"
done