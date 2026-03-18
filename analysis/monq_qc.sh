#!/bin/bash
#SBATCH --job-name=monq_qc
#SBATCH -A johnwayne
#SBATCH -t 14-00:00:00
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=jeon96@purdue.edu
#SBATCH --mail-type=END,FAIL


#########################
#DO NOT EDIT BELOW CODE #
#########################

n=20
cd /scratch/negishi/jeon96/monq/short-read/

mkdir -p ./qc/
mkdir -p ./cleaned/


# Data Quality Control
module load biocontainers
module load fastqc/0.12.1
module load trim-galore/0.6.10
module load cutadapt/3.7
module load kraken2/2.1.3

## initial quality check using Fastqc
cd ./raw/

for g in `ls -lt ./ | grep "fastq.gz" | tr -s ' ' | cut -d " " -f 9 | cut -d "R" -f 1 | sort | uniq`
  do
  fastqc ${g}R1_001.fastq.gz --extract --quiet -t 20 -o ../qc/
  fastqc ${g}R2_001.fastq.gz --extract --quiet -t 20 -o ../qc/
  
  cat ../qc/${g}R1_001_fastqc/summary.txt ../qc/${g}R2_001_fastqc/summary.txt > ../qc/${g}fastqc_summary.txt
  FAIL=$(grep "FAIL" ../qc/${g}fastqc_summary.txt)
  echo "raw"
  echo "$FAIL" 
done
#cd ../

## Adapter & low quality reads removal using Trim_galore

for g in `ls -lt ./ | grep "fastq.gz" | tr -s ' ' | cut -d " " -f 9 | cut -d "R" -f 1 | sort | uniq`
  do
  trim_galore --stringency 1 --length 30 --quality 20 --fastqc_args "--nogroup" -o ../cleaned/ --paired ${g}R1_001.fastq.gz ${g}R2_001.fastq.gz
  cat ../cleaned/${g}1.fastq_trimming_report.txt ../cleaned/${g}2.fastq_trimming_report.txt > ../cleaned/${g}fastqc_summary.txt
  FAIL=$(grep "FAIL" ../cleaned/${g}fastqc_summary.txt)
  echo "cleaned"
  echo "$FAIL"
done
#cd ../


# NM WGS samples from NCBI
cd ../
mkdir -p ./sra/raw/
mkdir -p ./sra/qc/
mkdir -p ./sra/cleaned/

module load biocontainers
module load sra-tools
N=64

cd ./sra/raw/
for g in SRR11514058 SRR11514059 SRR11514060 SRR11514061 SRR11514063 SRR11514064 SRR11514065
  do 
  mkdir ${g} 
  cd ${g} 
  prefetch ${g} --max-size 500GB -O ./
  fasterq-dump ${g} -e ${N} --progress; find . -name '*.fastq' -exec mv {} ../ \;
  cd ../
  rm -r ${g}
done

for g in `ls -lt ./ | grep "fastq" | tr -s ' ' | cut -d " " -f 9 | cut -d "_" -f 1 | sort | uniq`
do
  fastqc ${g}_1.fastq --extract --quiet -t 20 -o ../qc/
  fastqc ${g}_2.fastq --extract --quiet -t 20 -o ../qc/
  cat ../qc/${g}_2_fastqc/summary.txt ../qc/${g}_1_fastqc/summary.txt > ../qc/${g}_fastqc_summary.txt
  FAIL=$(grep "FAIL" ../qc/${g}_fastqc_summary.txt)
  echo "raw"
  echo "$FAIL"
done

cd ../cleaned
for g in `ls -lt ./ | grep "fastq" | tr -s ' ' | cut -d " " -f 9 | cut -d "_" -f 1 | sort | uniq`
  do
  trim_galore --stringency 1 --length 30 --quality 20 --fastqc_args "--nogroup" -o ./ --paired ../raw/${g}_1.fastq ../raw/${g}_2.fastq
  #cat ../cleaned/${g}_1.fastq_trimming_report.txt ../cleaned/${g}_2.fastq_trimming_report.txt > ../cleaned/${g}_fastqc_summary.txt
  unzip ${g}_1_val_1_fastqc.zip 
  unzip ${g}_2_val_2_fastqc.zip 
  cat ${g}_1_val_1_fastqc/summary.txt ${g}_2_val_2_fastqc/summary.txt > ${g}_fastqc_summary.txt
  FAIL=$(grep "FAIL" ${g}_fastqc_summary.txt)
  echo "cleaned"
  echo "$FAIL"
done