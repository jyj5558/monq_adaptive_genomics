#!/bin/bash
#SBATCH -A fnrpupfish
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 14-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --job-name=monq_mapping2
#SBATCH --mail-user=jeon96@purdue.edu
#SBATCH --mail-type=END,FAIL

module load biocontainers
module load bwa
module load samtools
module load picard
module load gatk
module load bamtools
module load qualimap

#Define variables to shorten commands
REF=/scratch/bell/jeon96/monq/ref/ref.fa
DICT=/scratch/bell/jeon96/monq/ref/ref.dict
FILT=/scratch/bell/jeon96/monq/ref/ok.bed
N=64
CLEAN_SRA=/scratch/bell/jeon96/monq/short-read/cleaned
FILT_SRA=/scratch/bell/jeon96/monq/short-read/filtered
FIXED_SRA=/scratch/bell/jeon96/monq/short-read/fixed

#Move to fastq containing directory
cd /scratch/bell/jeon96/monq/short-read/


bwa index $REF
samtools faidx $REF
PicardCommandLine CreateSequenceDictionary reference=$REF output=$DICT

mkdir -p ./aligned/
mkdir -p ./final_bams/

cd ./aligned/

for i in `ls -lrt ${CLEAN_SRA} | tr -s ' ' | cut -d ' ' -f 9 | cut -d 'R' -f 1 | sort | uniq`
    do
    id=$(echo $i | cut -d '_' -f 1)
    bwa mem -t ${N} -M -R "@RG\tID:group1\tSM:${i}\tPL:illumina\tLB:lib1\tPU:unit1" $REF ${FIXED_SRA}/${id}_fixed_R1.fq.gz ${FIXED_SRA}/${id}_fixed_R2.fq.gz > ${id}.sam
done

cd ../final_bams/

for i in `ls -lrt ${CLEAN_SRA} | tr -s ' ' | cut -d ' ' -f 9 | cut -d 'R' -f 1 | sort | uniq`
    do
    id=$(echo $i | cut -d '_' -f 1)
    echo "processing ${id}.sam"
    
    #Move to aligned directory
    cd ../aligned/
    
    #Validate sam file
    PicardCommandLine ValidateSamFile I=${id}.sam MODE=SUMMARY O=${id}.sam.txt

    #Sort validated sam file by read coordinate
    PicardCommandLine SortSam INPUT=${id}.sam OUTPUT=sorted_${id}.bam SORT_ORDER=coordinate

    #Get summary stats on initial alignments:
    samtools flagstat sorted_${id}.bam > ${id}_mapped.txt

    samtools depth -a sorted_${id}.bam | awk '{c++;s+=$3}END{print s/c}' > ${id}.pre.meandepth.txt 

    samtools depth -a sorted_${id}.bam | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' > ${id}.pre.1xbreadth.txt

    #Mark duplicates
    PicardCommandLine MarkDuplicates INPUT=sorted_${id}.bam OUTPUT=dedup_${id}.bam METRICS_FILE=metrics_${id}.bam.txt
    #Index in prep for realignment
    PicardCommandLine BuildBamIndex INPUT=dedup_${id}.bam
    #Local realignment of reads
    gatk3 -T RealignerTargetCreator -nt ${N} -R $REF -I dedup_${id}.bam -o ${id}_forIndelRealigner.intervals
    #Realign with established intervals
    export _JAVA_OPTIONS="-Xmx120G"
    gatk3 -T IndelRealigner -R $REF -I dedup_${id}.bam -targetIntervals ${id}_forIndelRealigner.intervals -o realgined_${id}.bam
    #Make new directory
    #Fix mate info
    PicardCommandLine FixMateInformation INPUT=realgined_${id}.bam OUTPUT=fixmate_${id}.bam SO=coordinate CREATE_INDEX=true
    #   Remove unmapped (4), secondary (256), QC failed (512), duplicate (1024), and
    #   supplementary (2048) reads from indel-realigned BAMs, and keep only reads
    #   mapped in a proper pair (2) to regions in a BED file (produced from QC_reference.sh)
    samtools view -@ ${N} -q 30 -b -F 3844 -f 2 -L $FILT fixmate_${id}.bam > ../final_bams/${id}_filt.bam 

    #Move into the final directory
    cd ../final_bams/
    #Index bam file
    PicardCommandLine BuildBamIndex INPUT=${id}_filt.bam

    samtools depth -a ${id}_filt.bam | awk '{c++;s+=$3}END{print s/c}' > ${id}.post.meandepth.txt

    samtools depth -a ${id}_filt.bam | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' > ${id}.post.1xbreadth.txt
done

for i in `ls -lrt ${CLEAN_SRA} | tr -s ' ' | cut -d ' ' -f 9 | cut -d 'R' -f 1 | sort | uniq`
  do
  id=$(echo $i | cut -d '_' -f 1)
  bam=${id}_filt.bam 
  sample=$(basename "$bam" _filt.bam)
  mkdir -p qualimap_reports/"$sample"
  qualimap bamqc -bam "$bam" --java-mem-size=120G -outdir qualimap_reports/"$sample" -outformat HTML
done