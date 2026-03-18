#!/bin/bash 
#SBATCH --job-name=snp_call_5
#SBATCH -A fnrpupfish 
#SBATCH -t 4:00:00 
#SBATCH -p cpu 
#SBATCH -q standby 
#SBATCH -N 1 
#SBATCH -n 10 
#SBATCH -e %x_%j.err 
#SBATCH -o %x_%j.out 
#SBATCH --mail-user=jeon96@purdue.edu 
#SBATCH --mail-type=END,FAIL 

i=5

module load biocontainers
module load picard
module load samtools

MCPAN=/scratch/bell/jeon96/monq/pan 
N=128 

cd ${MCPAN}
cd ./aligned 
rm -r jobs 
mkdir -p jobs 
mkdir -p ${MCPAN}/called/snp/chr${i}

ls -lrt /scratch/bell/jeon96/monq/pan/aligned/chr${i}/ | tr -s ' ' | cut -d ' ' -f 9 | grep "bam" | cut -d '_' -f 1 | sort | uniq > sample.list

#PicardCommandLine CreateSequenceDictionary --REFERENCE ${MCPAN}/monq-out_refp.fa --OUTPUT ${MCPAN}/monq-out_refp.dict
cd ./chr${i}/
for bam in *sorted.marked.bam; do
  samtools index -@ ${N} ${bam}  
done
cd ../

while read -a line; do
    id=${line[0]}
  cat <<EOF > ./jobs/${id}_snp.sh
#!/bin/bash
#SBATCH --job-name=${id}_snpcall_${i}
#SBATCH -A fnrpupfish
#SBATCH -t 4:00:00
#SBATCH -p cpu
#SBATCH -q standby
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out

i=${i}

MCPAN=${MCPAN}
N=${N} 

cd \${MCPAN} 
cd ./aligned

module load biocontainers 
module load gatk4 
      
gatk --java-options "-Xmx100G -XX:ParallelGCThreads=4" HaplotypeCaller -R \${MCPAN}/monq-out_refp.fa -I \${MCPAN}/aligned/chr\${i}/${id}_aug_chr\${i}_surject.sorted.marked.bam  -O \${MCPAN}/called/snp/chr\${i}/${id}_chr\${i}.haplotypecaller.vcf.gz --tmp-dir . --native-pair-hmm-threads \${N}

echo "[DONE] ${id}"
EOF
done < ./sample.list 

for i in `ls -1 ./jobs/*sh`; do echo "sbatch $i" ; done > slurmm_jobs ; source ./slurmm_jobs