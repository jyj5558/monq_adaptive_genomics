#!/bin/bash
#SBATCH --job-name=monq_vg
#SBATCH -A fnrpupfish
#SBATCH -t 14-00:00:00
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=jeon96@purdue.edu
#SBATCH --mail-type=END,FAIL

REF=/scratch/bell/jeon96/monq/ref/
MCPAN=/scratch/bell/jeon96/monq/pan
APP=/home/jeon96/app
N=64
CLEAN_SRA=/scratch/bell/jeon96/monq/short-read/cleaned
FIXED_SRA=/scratch/bell/jeon96/monq/short-read/fixed
FILT_SRA=/scratch/bell/jeon96/monq/short-read/filtered
CALLED=/scratch/bell/jeon96/monq/pan/called/sv/vg/chr1/

# Preprocessing the pangenome file
cd ${MCPAN}/
mkdir -p ./aligned/tmp

gunzip monq-out.gfa.gz 
${APP}/vg convert -g monq-out.gfa -v > monq-out.vg

## Modifying the pangenome to flip nodes' strands to reduce the number of times paths change strands
mkdir -p ./tmp/
${APP}/vg mod -t ${N} -O monq-out.vg > monq-out_mod.vg
#${APP}/vg index -p -t ${N} -x monq-out_mod.xg monq-out_mod.vg

## Chopping nodes in the graph so they are not more than 256 bp and index the graph
${APP}/vg mod -t ${N} -X 256 monq-out_mod.vg > monq-out_mod_chopped.vg
#${APP}/vg index -t ${N} -x monq-out_mod_chopped.xg monq-out_mod_chopped.vg

## Pruning the graph with kmer size 45 and index the graph
#${APP}/vg prune -t ${N} -k 45 monq-out_mod_chopped.vg > monq-out_mod_chopped_pruned.vg
#${APP}/vg index -t ${N} -b ./tmp -p -g monq-out_mod_chopped_pruned.gcsa monq-out_mod_chopped_pruned.vg

## Indexing the graph with -L option
${APP}/vg index -t ${N} -L -b ./tmp monq-out_mod_chopped.vg -x monq-out_mod_chopped_new_L.xg #-L: preserve alt paths in the xg

## Indexing for vg giraffe 
${APP}/vg convert -t ${N} -f monq-out_mod_chopped_new_L.xg > monq-out_mod_chopped.gfa 
${APP}/vg gbwt --progress --num-threads ${N} -d ./tmp --gbz-format -g monq-out_mod_chopped.gbz -G monq-out_mod_chopped.gfa 
${APP}/vg index --progress -t ${N} -b ./tmp -j monq-out_mod_chopped.dist monq-out_mod_chopped.gfa 
${APP}/vg minimizer --progress -t ${N} -d monq-out_mod_chopped.dist -o monq-out_mod_chopped.min monq-out_mod_chopped.gbz 

 
# Mapping short reads
cd ./aligned/ 
for i in `ls -lrt ${CLEAN_SRA} | tr -s ' ' | cut -d ' ' -f 9 | cut -d 'R' -f 1 | sort | uniq`; do
  id=$(echo $i | cut -d '_' -f 1)
  
## Aligning the individuals
  ${APP}/vg giraffe -p -t ${N} -Z ${MCPAN}/monq-out_mod_chopped.gbz -d ${MCPAN}/monq-out_mod_chopped.dist -m ${MCPAN}/monq-out_mod_chopped.min -x ${MCPAN}/monq-out_mod_chopped_new_L.xg -f ${FIXED_SRA}/${id}_fixed_R1.fq.gz -f  ${FIXED_SRA}/${id}_fixed_R2.fq.gz > ${id}.gam
  #${APP}/vg map -t ${N} -f ${FIXED_SRA}/${id}_fixed_R1.fq.gz -f ${FIXED_SRA}/${id}_fixed_R2.fq.gz -x ../monq-out_mod_chopped.xg -g ../monq-out_mod_chopped_pruned.gcsa > ${i}aln.gam

## Filtering secondary and ambiguous read mappings out of the gam (for SV detection for now) -> filtering suppressed for comparisons with linear assemblies
  ${APP}/vg filter -t ${N} ${id}.gam -r 0.90 -fu -m 1 -q 15 -D 999 -x ${MCPAN}/monq-out_mod_chopped_new_L.xg > ${id}_filt.gam
  #${APP}/vg filter -t ${N} ${i}aln.gam -r 0.90 -fu -m 1 -q 15 -D 999 -x ../monq-out_mod_chopped.xg > ${i}aln.filtered.gam
done

export TMPDIR=/scratch/bell/jeon96/monq/pan/aligned/tmp 


# Sorting individual gams
for i in `ls -lrt ${CLEAN_SRA} | tr -s ' ' | cut -d ' ' -f 9 | cut -d 'R' -f 1 | sort | uniq`; do
  id=$(echo $i | cut -d '_' -f 1)
  echo "processing $id"
  ${APP}/vg gamsort -t ${N} -p ${id}_filt.gam -i ${id}_filt_sorted.gam.gai > ${id}_filt_sorted.gam
  echo "processed $id"
done
#cat $(for i in $(ls -lrt ${CLEAN_SRA} | tr -s ' ' | cut -d ' ' -f 9 | cut -d 'R' -f 1 | sort | uniq | cut -d '_' -f 1); do echo ${i}_filt_sorted.gam; done) > combined_filt_sorted.gam
#cat $(for i in $(ls -lrt ${CLEAN_SRA} | tr -s ' ' | cut -d ' ' -f 9 | cut -d 'R' -f 1 | sort | uniq | cut -d '_' -f 1); do echo ${i}_filt.gam; done) > combined_filt.gam
#${APP}/vg gamsort -t ${N} -p combined_filt.gam -i combined_filt_sorted.gam.gai > combined_filt_sorted.gam 

${APP}/vg paths -M -x ${MCPAN}/monq-out_mod_chopped_new_L.xg | grep "ref_p" | cut -f 1 | sort -V > chr_list.txt


# Chunking individual gams for each chromosome
for i in `ls -lrt ${CLEAN_SRA} | tr -s ' ' | cut -d ' ' -f 9 | cut -d 'R' -f 1 | sort | uniq`; do
  id=$(echo $i | cut -d '_' -f 1)
  mkdir -p ./${id}_chunk
  cd ./${id}_chunk
  echo "chuncking $id at $(date +\"%T\")"
  ${APP}/vg chunk -t ${N} -P ../chr_list.txt -c 9999999 -O pg -g -x ${MCPAN}/monq-out_mod_chopped_new_L.xg -a ../${id}_filt_sorted.gam
  echo "chuncked $id at $(date +\"%T\")"
  cd ../
done
#${APP}/vg chunk -t ${N} -P chr_list.txt -c 9999999 -O pg -g -x ${MCPAN}/monq-out_mod_chopped_new_L.xg -a combined_filt_sorted.gam


# Augmenting each chromosome graph using an individual-combined gam
cd ../
mkdir -p ./augmented/tmp/
cd ./augmented/

for i in {1..30}; do
  mkdir -p ./chr${i}
  cp ../aligned/F1958_chunk/*chr${i}_*.vg ./chr${i}.vg
  cat $(for j in $(ls -lrt ../aligned | grep "chunk" | tr -s ' ' | cut -d ' ' -f 9 | sort | uniq | cut -d '_' -f 1); do echo ../aligned/${j}_chunk/*chr${i}_*.gam; done) > chr${i}_combined_filt_sorted.gam
  ${APP}/vg augment -t ${N} chr${i}.vg chr${i}_combined_filt_sorted.gam -s -m 3 -q 5 -Q 5 -p -A ./chr${i}/monq_chr${i}_aug.gam > ./chr${i}/monq_chr${i}_aug.pg # chr1-5 need high memory
done

#${APP}/vg combine -p $(for i in {1..30}; do echo ./chr${i}/monq_chr${i}_aug.pg; done) > monq_merged_aug.pg # need high memory; downstream processing steps does not work due to oom
#${APP}/vg ids -j $(for i in {1..30}; do echo ./chr${i}/monq_chr${i}_aug.pg; done) # vg combine do this on the fly
#cat $(for i in $(seq 1 22; echo X; echo Y); do echo $i.vg; done) > monq

## Augmenting the graph with all variation from the GAM -> not working due to too large files combined
#${APP}/vg convert -t ${N} ../monq-out.xg -p > ../monq-out.pg
#${APP}/vg augment -t ${N} ../monq-out.pg combined_aln.filtered.gam -s -m 3 -q 5 -Q 5 -A ../monq-out_aug.gam > ../monq-out_aug.pg #-c 5 will filter out variants if their read depth < 5

## Augmenting the graph iteratively -> not efficient
#${APP}/vg convert -t ${N} ${MCPAN}/monq-out_mod_chopped_new_L.xg -p > ${MCPAN}/monq-out_mod_chopped.pg
#cp ../monq-out.pg running_aug.pg

#for i in `ls -lrt ${CLEAN_SRA} | tr -s ' ' | cut -d ' ' -f 9 | cut -d 'R' -f 1 | sort | uniq`; do
#    id=$(echo $i | cut -d '_' -f 1)
#    gam=${id}_aln.filtered.gam
#    echo "Augmenting with ${gam}"
#    ${APP}/vg augment -t ${N} running_aug.pg ${gam} -s -m 3 -q 5 -Q 5 > tmp_aug.pg
#    cp tmp_aug.pg ${id}_aug.pg
#    mv tmp_aug.pg running_aug.pg
#done


# Indexing the augmented graph 
for i in {1..30}; do
  cd ./chr${i}
  mkdir -p ./tmp/
  
  ${APP}/vg mod -t ${N} -X 256 monq_chr${i}_aug.pg > monq_chr${i}_aug_chopped.pg
  #${APP}/vg index -t ${N} -p -x monq_chr${i}_aug_chopped.xg monq_chr${i}_aug_chopped.pg
  ${APP}/vg index -t ${N} -L -b ./tmp -p monq_chr${i}_aug_chopped.pg -x monq_chr${i}_aug_chopped_L.xg #-L: preserve alt paths in the xg
  ${APP}/vg convert -t ${N} -f monq_chr${i}_aug_chopped_L.xg > monq_chr${i}_aug_chopped.gfa 
  #${APP}/vg gbwt --num-threads ${N} -d ./tmp/ -p -G monq_chr${i}_aug_chopped.gfa -o monq_chr${i}_aug_chopped.gbwt # for vg map
  #${APP}/vg ids -m monq_chr${i}_mapping monq_chr${i}_aug_chopped.pg # for vg map
  #cp monq_chr${i}_mapping monq_chr${i}_mapping.backup # for vg map
  #${APP}/vg prune -t ${N} -p -k 45 -u -g monq_chr${i}_aug_chopped.gbwt -a -m monq_chr${i}_mapping monq_chr${i}_aug_chopped.pg > monq_chr${i}_aug_chopped_pruned_new.pg # for vg map
  #${APP}/vg index -t ${N} -b ./tmp -p -g monq_chr${i}_aug_chopped_pruned.gcsa -f monq_chr${i}_mapping monq_chr${i}_aug_chopped_pruned.pg # for vg map
  #${APP}/vg convert -t ${N} monq_chr${i}_aug_chopped_L.xg -p > monq_chr${i}_aug_chopped_L.pg # for vg map
  ${APP}/vg gbwt --progress --num-threads ${N} -d ./tmp --gbz-format -g monq_chr${i}_aug_chopped.gbz -G monq_chr${i}_aug_chopped.gfa 
  ${APP}/vg index --progress -b ./tmp/ -t ${N} --snarl-limit 100000 -j monq_chr${i}_aug_chopped.dist monq_chr${i}_aug_chopped.gbz
  ${APP}/vg minimizer --progress -t ${N} -d monq_chr${i}_aug_chopped.dist -z monq_chr${i}_aug_chopped.zip -o monq_chr${i}_aug_chopped.min monq_chr${i}_aug_chopped.gbz  
  
  cd ../
done


# Computing pangenome stats
for i in {1..30}; do
  cd ./chr${i}
  ${APP}/vg stats --threads ${N} --progress --node-count --edge-count -lz monq_chr${i}_aug.pg > monq_chr${i}.stats
  cd ../
done


# Mapping short reads to the augmented graphs - vg map takes too long; goes with vg giraffe and per-sample, per-chromosome version in stadby queue
cd ../aligned/ 
for i in {1..30}; do
  for j in `ls -lrt ${CLEAN_SRA} | tr -s ' ' | cut -d ' ' -f 9 | cut -d 'R' -f 1 | sort | uniq`; do
    id=$(echo $j | cut -d '_' -f 1)
  
## Aligning the individuals
  ${APP}/vg giraffe -p -t ${N} -Z ../augmented/chr${i}/monq_chr${i}_aug_chopped.gbz -d ../augmented/chr${i}/monq_chr${i}_aug_chopped.dist -m ../augmented/chr${i}/monq_chr${i}_aug_chopped.min -z ../augmented/chr${i}/monq_chr${i}_aug_chopped.zip -x ../augmented/chr${i}/monq_chr${i}_aug_chopped_L.xg -f ${FIXED_SRA}/${id}_fixed_R1.fq.gz -f ${FIXED_SRA}/${id}_fixed_R2.fq.gz > ./chr${i}/${id}_aln_aug_chr${i}.gam
    #${APP}/vg map -t ${N} -f ${FIX_SRA}/${id}_fixed_R1.fq.gz -f ${FIX_SRA}/${id}_fixed_R2.fq.gz -x ../augmented/chr${i}/monq_chr${i}_aug_chopped_L.xg -g ../augmented/chr${i}/monq_chr${i}_aug_chopped_pruned.gcsa > ./chr${i}/${id}_aln_aug_chr${i}.gam

## Filtering secondary and ambiguous read mappings out of the gam (for SV detection for now) 
  ${APP}/vg filter -t ${N} ./chr${i}/${id}_aln_aug_chr${i}.gam -r 0.90 -fu -m 1 -q 30 -D 999 --proper-pairs -x ../augmented/chr${i}/monq_chr${i}_aug_chopped_L.xg > ./chr${i}/${id}_aln_aug_chr${i}.filt.gam
  done
done


# Indexing the augmented graph (whole-genome version; not working due to oom)
#${APP}/vg mod -t ${N} -X 256 monq_merged_aug.pg > monq_merged_aug_chopped.pg
#${APP}/vg index -t ${N} -b ./tmp -x monq_merged_aug_chopped.xg monq_merged_aug_chopped.pg
#${APP}/vg prune -t ${N} -k 45 monq_merged_aug_chopped.pg > monq_merged_aug_chopped_pruned.pg
#${APP}/vg index -t ${N} -b ./tmp -p -g monq_merged_aug_chopped_pruned.gcsa monq_merged_aug_chopped_pruned.pg

# Indexing the augmented graph with -L option
#${APP}/vg index -t ${N} -L -b ./tmp monq_merged_aug_chopped.pg -x monq_merged_aug_chopped_new_L.xg #-L: preserve alt paths in the xg
#${APP}/vg convert -t ${N} monq_merged_aug_chopped_new_L.xg -p > monq_merged_aug_new.pg
#${APP}/vg convert -t ${N} -f monq_merged_aug_chopped_new_L.xg > monq_aug_chopped.gfa 
#${APP}/vg gbwt --progress --num-threads ${N} -d ./tmp --gbz-format -g monq_aug_chopped.gbz -G monq_aug_chopped.gfa 
#${APP}/vg index --progress -t ${N} -b ./tmp -j monq_aug_chopped.dist monq_aug_chopped.gfa 
#${APP}/vg minimizer --progress -t ${N} -d monq_aug_chopped.dist -o monq_aug_chopped.min monq_aug_chopped.gbz 


# Variant calling (for SVs) for each individual using vg
mkdir -p ../called/sv/vg/

## Computing the snarls
for i in {1..30}; do
  ${APP}/vg snarls -t ${N} ../augmented/chr${i}/monq_chr${i}_aug_chopped.gbz > ../augmented/chr${i}/monq_chr${i}_aug.snarls

## Computing the support
  for j in `ls -lrt ${CLEAN_SRA} | tr -s ' ' | cut -d ' ' -f 9 | cut -d 'R' -f 1 | sort | uniq`; do
    id=$(echo $j | cut -d '_' -f 1)
    #${APP}/vg pack -t ${N} -x ./chr${i}/monq_chr${i}_aug.pg -g ./chr${i}/monq_chr${i}_aug.gam -Q 5 -o ./chr${i}/monq_chr${i}_aug.pack #-Q 5: ignore mapping and base qualitiy < 5
    ${APP}/vg pack -t ${N} -x ../augmented/chr${i}/monq_chr${i}_aug_chopped.gbz -g ./chr${i}/${id}_aln_aug_chr${i}.filt.gam -Q 5 -o ./chr${i}/${id}_chr${i}_aug.pack #-Q 5: ignore mapping and base qualitiy < 5
  done
done

## Calling variants; run this step using highmem queue (otherwise it can't be finished)
#for i in `ls -lrt ${FILT_SRA} | tr -s ' ' | cut -d ' ' -f 9 | cut -d 'R' -f 1 | sort | uniq`; do
#  ${APP}/vg call -t ${N} monq-out_augSV_new.pg -r monq-out_augSV.snarls -k ${i}augSV.pack -s ${i} -a -A -c 50 -C 100000 --progress > ../called/sv/vg/${i}mcSV.vcf #-a: calling every snarl using the same coordinates and including reference calls; -c: minimum length of variants to be called; -C: maximum length of variants to be called
#done

## Calling variants
for i in {1..30}; do #chr1-3, 29 need high memory
  mkdir -p ../called/sv/vg/chr${i}/
  for j in `ls -lrt ${FIXED_SRA} | tr -s ' ' | cut -d ' ' -f 9 | grep "s.fq.gz" | sort | uniq`; do
    id=$(echo $j | cut -d '_' -f 1)
    ${APP}/vg call -t ${N} ../augmented/chr${i}/monq_chr${i}_aug_chopped.gbz -r ../augmented/chr${i}/monq_chr${i}_aug.snarls -k ./chr${i}/${id}_chr${i}_aug.pack -s ${id} -a -A -c 50 -C 100000 --progress > ../called/sv/vg/chr${i}/${id}_chr${i}_mcSV.vcf
  done
done

## Concatenate variants of individual chromosomes per sample
module load biocontainers
module load bcftools
module load htslib

cd ../called/sv/vg/

for i in `ls -lrt ${FIXED_SRA} | tr -s ' ' | cut -d ' ' -f 9 | grep "s.fq.gz" | sort | uniq`; do
  id=$(echo $i | cut -d '_' -f 1)
  for j in {1..30}; do 
    bgzip ./chr${j}/${id}_chr${j}_mcSV.vcf
    tabix -f ./chr${j}/${id}_chr${j}_mcSV.vcf.gz
    bcftools sort ./chr${j}/${id}_chr${j}_mcSV.vcf.gz -Oz -o ./chr${j}/${id}_chr${j}_sorted_mcSV.vcf.gz
  done
  
  samples=""
  for j in {1..30}; do 
    samples+="./chr${j}/${id}_chr${j}_sorted_mcSV.vcf.gz "
  done
  
  bcftools concat ${samples} --threads ${N} -Oz -o ${id}_concat_mcSV.vcf.gz
  bcftools index ${id}_concat_mcSV.vcf.gz --threads ${N}
done


# Variant calling (for SNPs) for each individual using vg
module load biocontainers
module load samtools
module load picard
module load gatk4

mkdir -p ${MCPAN}/called/snp/
cd ${MCPAN}/aligned/ 

for i in {1..30}; do
## Creating a list of reference paths
  ${APP}/vg paths -x ../augmented/chr${i}/monq_chr${i}_aug.pg -S ref_p -L > ./chr${i}/refp_chr${i}_paths.txt # use .full graph just to make comprehensive path lists
done

for i in {1..30}; do
  for j in `ls -lrt ${CALLED} | tr -s ' ' | cut -d ' ' -f 9 | grep "vcf.gz.tbi" | sort | uniq`; do
    id=$(echo $j | cut -d '_' -f 1)
  
## Filtering secondary and ambiguous read mappings out of the gam (for SNP detection)
    ${APP}/vg filter -t ${N} ./chr${i}/${id}_aln_aug_chr${i}.gam -r 0.90 -fu -m 1 -q 15 -D 999 -i --proper-pairs -x ../augmented/chr${i}/monq_chr${i}_aug_chopped_L.xg > ./chr${i}/${id}_aln_aug_intlvd_chr${i}.filt.gam
    ${APP}/vg surject -x ../augmented/chr${i}/monq_chr${i}_aug_chopped_L.xg ./chr${i}/${id}_aln_aug_intlvd_chr${i}.filt.gam --threads ${N} --prune-low-cplx --interleaved -F ./chr${i}/refp_chr${i}_paths.txt -b -N ${id} -R "ID:1 LB:lib1 SM:${id} PL:illumina PU:unit1" > ./chr${i}/${id}_aug_chr${i}_surject.bam
  done
done

## Processing bams 
### Indexing bam files
for i in {1..30}; do
  for j in `ls -lrt ${CALLED} | tr -s ' ' | cut -d ' ' -f 9 | grep "vcf.gz.tbi" | sort | uniq`; do
    id=$(echo $j | cut -d '_' -f 1)
    samtools sort -@ ${N} -o ./chr${i}/${id}_aug_chr${i}_surject.sorted.bam ./chr${i}/${id}_aug_chr${i}_surject.bam 
    samtools index -@ ${N} -b ./chr${i}/${id}_aug_chr${i}_surject.sorted.bam  
  done
done

### Marking duplicates
for i in {1..30}; do
  for j in `ls -lrt ${CALLED} | tr -s ' ' | cut -d ' ' -f 9 | grep "vcf.gz.tbi" | sort | uniq`; do
    id=$(echo $j | cut -d '_' -f 1)
    PicardCommandLine MarkDuplicates --INPUT ./chr${i}/${id}_aug_chr${i}_surject.sorted.bam --OUTPUT ./chr${i}/${id}_aug_chr${i}_surject.sorted.marked.bam --METRICS_FILE ./chr${i}/${id}_aug_chr${i}_surject.metrics
    PicardCommandLine BuildBamIndex --INPUT ./chr${i}/${id}_aug_chr${i}_surject.sorted.marked.bam --OUTPUT ./chr${i}/${id}_aug_chr${i}_surject.sorted.marked.bam.bai
  done
done

### Merging per-chromosome bam files per individual
for j in `ls -lrt ${CALLED} | tr -s ' ' | cut -d ' ' -f 9 | grep "vcf.gz.tbi" | sort | uniq`; do
  id=$(echo $j | cut -d '_' -f 1)
  files=()
  for i in {1..30}; do
    f="./chr${i}/${id}_aug_chr${i}_surject.sorted.marked.bam"
    [[ -s "${f}" ]] && files+=("${f}") || echo "[$id] missing: ${f}" >&2
  done
  samtools merge -c -p -@ ${N} -o ${id}_aug_surject.merged.bam ${files[@]} # apply -f to overwrite
  samtools sort -@ ${N} -o ${id}_aug_surject.merged.sorted.bam ${id}_aug_surject.merged.bam 
  samtools index -@ ${N} -b ${id}_aug_surject.merged.sorted.bam
done

for id in E8947 F1856; do # reheader some bam files that did not map to some chromosomes
  mv ${id}_aug_surject.merged.sorted.bam ${id}_temp.bam
  samtools reheader -P F1958_aug_surject.merged.sorted.bam ${id}_temp.bam > ${id}_aug_surject.merged.sorted.bam
  samtools index -@ ${N} -b ${id}_aug_surject.merged.sorted.bam
done

## Calling SNPs using GATK4 haplotype caller
PicardCommandLine CreateSequenceDictionary --REFERENCE ${MCPAN}/monq-out_refp.fa --OUTPUT ${MCPAN}/monq-out_refp.dict

for i in {1..30}; do # run this per-chromosome; see the script "monq_snpcall_byID.sh"
  cd ./chr${i}/
  for bam in *sorted.marked.bam; do
    samtools index -@ ${N} ${bam}  
  done
  cd ../
done

for i in {1..30}; do # run this per-chromosome; see the script "monq_snpcall_byID.sh"
  mkdir -p ${MCPAN}/called/snp/chr${i}

  for j in `ls -lrt ${CALLED} | tr -s ' ' | cut -d ' ' -f 9 | grep "vcf.gz.tbi" | sort | uniq`; do
    id=$(echo $j | cut -d '_' -f 1)
    gatk --java-options "-Xmx100G -XX:ParallelGCThreads=4" HaplotypeCaller -R ${MCPAN}/monq-out_refp.fa -I ${MCPAN}/aligned/chr${i}/${id}_aug_chr${i}_surject.sorted.marked.bam  -O ${MCPAN}/called/snp/chr${i}/${id}_chr${i}.haplotypecaller.vcf.gz --tmp-dir . --native-pair-hmm-threads ${N}
  done
done

#for i in `ls -lrt ${CALLED} | tr -s ' ' | cut -d ' ' -f 9 | grep "vcf.gz.tbi" | sort | uniq`; do
# id=$(echo $i | cut -d '_' -f 1)
# echo "processing ${id}"
# for j in {1..30}; do 
#   vcf=./chr${j}/${id}_chr${j}.haplotypecaller.vcf.gz
#   [[ ! -s "$vcf" ]] && continue
#   echo "${id}" > samples_${id}_chr${j}.txt

#	  bcftools reheader \
#      	--samples samples_${id}_chr${j}.txt --threads ${N} \
#      	-o ${vcf%.vcf.gz}.SMfixed.vcf.gz \
#      	"$vcf"

#	  mv ${vcf%.vcf.gz}.SMfixed.vcf.gz "$vcf"
#   rm samples_${id}_chr${j}.txt
#	  echo "processing ${id} done"
#   done
#done

## Concatenate variants of individual chromosomes per sample
module load biocontainers
module load bcftools
module load htslib
module load vcftools

cd ../called/snp/

for i in `ls -lrt ${CALLED} | tr -s ' ' | cut -d ' ' -f 9 | grep "vcf.gz.tbi" | sort | uniq`; do
  id=$(echo $i | cut -d '_' -f 1)
  echo "processing ${id}"
  for j in {1..30}; do 
    bcftools sort ./chr${j}/${id}_chr${j}.haplotypecaller.vcf.gz -Oz -o ./chr${j}/${id}_chr${j}_sorted_mcSNP.vcf.gz
    echo "sorting chr${j} done"
  done
  
  files=()
  for j in {1..30}; do 
    vcf=./chr${j}/${id}_chr${j}_sorted_mcSNP.vcf.gz
    [[ -s "$vcf" ]] && files+=("$vcf")
  done
  
  bcftools concat ${files[@]} --threads ${N} -Oz -o ${id}_concat_mcSNP.vcf.gz
  bcftools index ${id}_concat_mcSNP.vcf.gz --threads ${N}
done

## Compressing and indexing each vcf file first
for i in `ls -lrt ${CALLED} | tr -s ' ' | cut -d ' ' -f 9 | grep "vcf.gz.tbi" | sort | uniq`; do
  id=$(echo $i | cut -d '_' -f 1)
  gunzip ${id}_concat_mcSNP.vcf.gz 
  sed 's/ref_p#0#//g' ${id}_concat_mcSNP.vcf > ${id}_mcSNP_temp.vcf # polish the header line to be compatible with bcftools
  bcftools sort ${id}_mcSNP_temp.vcf -Oz -o ${id}_mcSNP.sorted.vcf.gz
  bcftools index ${id}_mcSNP.sorted.vcf.gz --threads ${N}
  bcftools view ${id}_mcSNP.sorted.vcf.gz --threads ${N} | grep -v "##" | wc -l # number of variants
  rm ${id}_mcSNP_temp.vcf
done

## Combining separately called SNP vcf files
bcftools merge -m all -Oz -o monq_mcSNP.merged.vcf.gz --threads ${N} *_mcSNP.sorted.vcf.gz 


# Filtering with population-level parameters
vcftools --gzvcf monq_mcSNP.merged.vcf.gz --missing-indv 
vcftools --gzvcf monq_mcSNP.merged.vcf.gz --missing-site
vcftools --gzvcf monq_mcSNP.merged.vcf.gz --depth
vcftools --gzvcf monq_mcSNP.merged.vcf.gz --site-mean-depth
vcftools --gzvcf monq_mcSNP.merged.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.2' | cut -f1,2 >> badloci
bcftools stats monq_mcSNP.merged.vcf.gz > vcf-stats_snp.txt

# In R, plotting and summarizing vcf stats
module load r

R
library(tidyverse)
library(ggplot2)
var_qual <- read_delim("./out.lqual", delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
a + theme_light() + xlim(0,500)
summary(var_qual$qual)
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    30.0   108.6   182.6   229.2   293.0 99318.1

var_depth <- read_delim("./out.ldepth.mean", delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()   
a + theme_light() + xlim(0,100)        
summary(var_depth$mean_depth)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#   1.000    5.000    6.828    7.340    8.667 1389.000
quantile(var_depth$mean_depth, probs = seq(0, 1, 1/10))
#        0%        10%        20%        30%        40%        50%        60% 
#   1.00000    3.25000    4.33333    5.00000    6.00000    6.82812    7.50000
#       70%        80%        90%       100%
#   8.06897    9.00000   11.00000 1389.00000
quantile(var_depth$mean_depth, probs = c(0.01, 0.05, 0.95, 0.99))
#     1%      5%     95%     99%
# 2.0000  2.7500 12.6667 20.00007
 
var_miss <- read_delim("./out.lmiss", delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(var_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(var_miss$fmiss)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.01099 0.91209 0.97802 0.88939 0.98901 0.98901
                   
ind_depth <- read_delim("./out.idepth", delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)
a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(ind_depth$depth)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  1.988   7.013   7.480   7.576   7.952  11.343

ind_miss  <- read_delim("./out.imiss", delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(ind_miss$fmiss)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.8537  0.8782  0.8874  0.8894  0.8970  1.0000

quit()

# Filtering vcf file based on the determined conservative cut-off values from above results                                                                      
vcftools --gzvcf monq_mcSNP.merged.vcf.gz --out monq_mcSNP.filtered --minQ 20 --mac 2 --min-meanDP 6 --max-missing 0.8 --exclude-positions badloci --max-meanDP 15 --recode --recode-INFO-all #--remove-filtered lowad  
bcftools sort monq_mcSNP.filtered.recode.vcf -Oz -o monq_mcSNP.filtered.vcf.gz

bcftools norm -m -any --threads ${N} monq_mcSNP.filtered.vcf.gz -Ov -o monq_mcSNP.filt.decomp.vcf

vcftools --vcf monq_mcSNP.filt.decomp.vcf --missing-site --out out2
cat out2.lmiss | awk '!/CHR/' | awk '$6 > 0.2' | cut -f1,2 >> badloci2
vcftools --vcf monq_mcSNP.filt.decomp.vcf  --out monq_mcSNP.final --minQ 20 --mac 2 --max-missing 0.8 --exclude-positions badloci2 --recode --recode-INFO-all
bcftools sort monq_mcSNP.final.recode.vcf > monq_mcSNP.final.sort.vcf



# New Mexico samples from NCBI
cd /scratch/bell/jeon96/monq/pan/aligned/
mkdir -p ./new_mexico/
cd ./new_mexico/

for i in {1..30}; do 
  for j in SRR11514058 SRR11514059 SRR11514060 SRR11514061 SRR11514063 SRR11514064 SRR11514065; do 
    mkdir -p chr${i}; 
    ${APP}/vg giraffe -p -t ${N} -Z ../../augmented/chr${i}/monq_chr${i}_aug_chopped.gbz -d ../../augmented/chr${i}/monq_chr${i}_aug_chopped.dist -m ../../augmented/chr${i}/monq_chr${i}_aug_chopped.min -z ../../augmented/chr${i}/monq_chr${i}_aug_chopped.zip -x ../../augmented/chr${i}/monq_chr${i}_aug_chopped_L.xg -f /scratch/bell/jeon96/monq/short-read/sra/fixed/${j}_fixed_R1.fq.gz -f /scratch/bell/jeon96/monq/short-read/sra/fixed/${j}_fixed_R2.fq.gz > ./chr${i}/${j}_aln_aug_chr${i}.gam
    ${APP}/vg filter -t ${N} ./chr${i}/${j}_aln_aug_chr${i}.gam -r 0.90 -fu -m 1 -q 30 -D 999 --proper-pairs -x ../../augmented/chr${i}/monq_chr${i}_aug_chopped_L.xg > ./chr${i}/${j}_aln_aug_chr${i}.filt.gam
  done
done

## Computing the support
for i in {1..30}; do 
  for j in SRR11514058 SRR11514059 SRR11514060 SRR11514061 SRR11514063 SRR11514064 SRR11514065; do 
    ${APP}/vg pack -t ${N} -x ../../augmented/chr${i}/monq_chr${i}_aug_chopped.gbz -g ./chr${i}/${j}_aln_aug_chr${i}.filt.gam -Q 5 -o ./chr${i}/${j}_chr${i}_aug.pack #-Q 5: ignore mapping and base qualitiy < 5
  done
done

## Calling variants
for i in {1..30}; do #chr1-3, 29 need high memory
  mkdir -p ../../called/sv/vg/new_mexico/chr${i}/
  for id in SRR11514058 SRR11514059 SRR11514060 SRR11514061 SRR11514063 SRR11514064 SRR11514065; do 
    ${APP}/vg call -t ${N} ../../augmented/chr${i}/monq_chr${i}_aug_chopped.gbz -r ../../augmented/chr${i}/monq_chr${i}_aug.snarls -k ./chr${i}/${id}_chr${i}_aug.pack -s ${id} -a -A -c 50 -C 100000 --progress > ../../called/sv/vg/new_mexico/chr${i}/${j}_chr${i}_mcSV.vcf
  done
done

## Concatenate variants of individual chromosomes per sample
module load biocontainers
module load bcftools
module load htslib

cd ../../called/sv/vg/new_mexico/

for i in SRR11514058 SRR11514059 SRR11514060 SRR11514061 SRR11514063 SRR11514064 SRR11514065; do
  #id=$(echo $i | cut -d '_' -f 1)
  for j in {1..30}; do 
    bgzip ./chr${j}/${i}_chr${j}_mcSV.vcf
    tabix -f ./chr${j}/${i}_chr${j}_mcSV.vcf.gz
    bcftools sort ./chr${j}/${i}_chr${j}_mcSV.vcf.gz -Oz -o ./chr${j}/${i}_chr${j}_sorted_mcSV.vcf.gz
  done
done

for i in SRR11514058 SRR11514059 SRR11514060 SRR11514061 SRR11514063 SRR11514064 SRR11514065; do
  files=()
  for j in {1..30}; do 
    vcf=./chr${j}/${i}_chr${j}_sorted_mcSV.vcf.gz
    [[ -s "$vcf" ]] && files+=("$vcf")
  done
  
  bcftools concat ${files[@]} --threads ${N} -Oz -o ${i}_concat_mcSV.vcf.gz
  bcftools index -f ${i}_concat_mcSV.vcf.gz --threads ${N}
done

cd ${MCPAN}/aligned/ 

for i in {1..30}; do
## Creating a list of reference paths
  ${APP}/vg paths -x ../augmented/chr${i}/monq_chr${i}_aug.pg -S ref_p -L > ./chr${i}/refp_chr${i}_paths.txt # use .full graph just to make comprehensive path lists
done

cd ./new_mexico

for i in {1..30}; do
  for id in SRR11514058 SRR11514059 SRR11514060 SRR11514061 SRR11514063 SRR11514064 SRR11514065; do
  
## Filtering secondary and ambiguous read mappings out of the gam (for SNP detection)
    ${APP}/vg filter -t ${N} ./chr${i}/${id}_aln_aug_chr${i}.gam -r 0.90 -fu -m 1 -q 15 -D 999 -i --proper-pairs -x ../../augmented/chr${i}/monq_chr${i}_aug_chopped_L.xg > ./chr${i}/${id}_aln_aug_intlvd_chr${i}.filt.gam
    ${APP}/vg surject -x ../../augmented/chr${i}/monq_chr${i}_aug_chopped_L.xg ./chr${i}/${id}_aln_aug_intlvd_chr${i}.filt.gam --threads ${N} --prune-low-cplx --interleaved -F ../chr${i}/refp_chr${i}_paths.txt -b -N ${id} -R "ID:1 LB:lib1 SM:${id} PL:illumina PU:unit1" > ./chr${i}/${id}_aug_chr${i}_surject.bam
  done
done

## Processing bams 
### Indexing bam files
for i in {1..30}; do
  for id in SRR11514058 SRR11514059 SRR11514060 SRR11514061 SRR11514063 SRR11514064 SRR11514065; do
    samtools sort -@ ${N} -o ./chr${i}/${id}_aug_chr${i}_surject.sorted.bam ./chr${i}/${id}_aug_chr${i}_surject.bam 
    samtools index -@ ${N} -b ./chr${i}/${id}_aug_chr${i}_surject.sorted.bam  
  done
done

### Marking duplicates
for i in {1..30}; do
  for id in SRR11514058 SRR11514059 SRR11514060 SRR11514061 SRR11514063 SRR11514064 SRR11514065; do
    PicardCommandLine MarkDuplicates --INPUT ./chr${i}/${id}_aug_chr${i}_surject.sorted.bam --OUTPUT ./chr${i}/${id}_aug_chr${i}_surject.sorted.marked.bam --METRICS_FILE ./chr${i}/${id}_aug_chr${i}_surject.metrics
    PicardCommandLine BuildBamIndex --INPUT ./chr${i}/${id}_aug_chr${i}_surject.sorted.marked.bam --OUTPUT ./chr${i}/${id}_aug_chr${i}_surject.sorted.marked.bam.bai
  done
done

### Merging per-chromosome bam files per individual
for id in SRR11514058 SRR11514059 SRR11514060 SRR11514061 SRR11514063 SRR11514064 SRR11514065; do
  files=()
  for i in {1..30}; do
    f="./chr${i}/${id}_aug_chr${i}_surject.sorted.marked.bam"
    [[ -s "${f}" ]] && files+=("${f}") || echo "[$id] missing: ${f}" >&2
  done
  samtools merge -c -p -@ ${N} -o ${id}_aug_surject.merged.bam ${files[@]} # apply -f to overwrite
  samtools sort -@ ${N} -o ${id}_aug_surject.merged.sorted.bam ${id}_aug_surject.merged.bam 
  samtools index -@ ${N} -b ${id}_aug_surject.merged.sorted.bam
done

## Calling SNPs using GATK4 haplotype caller
PicardCommandLine CreateSequenceDictionary --REFERENCE ${MCPAN}/monq-out_refp.fa --OUTPUT ${MCPAN}/monq-out_refp.dict

for i in {1..30}; do # run this per-chromosome; see the script "monq_snpcall_byID.sh"
  cd ./chr${i}/
  for bam in *sorted.marked.bam; do
    samtools index -@ ${N} ${bam}  
  done
  cd ../
done

for i in {1..30}; do # run this per-chromosome; see the script "monq_snpcall_byID.sh"
  mkdir -p ${MCPAN}/called/snp/new_mexico/chr${i}

  for id in SRR11514058 SRR11514059 SRR11514060 SRR11514061 SRR11514063 SRR11514064 SRR11514065; do
    gatk --java-options "-Xmx100G -XX:ParallelGCThreads=4" HaplotypeCaller -R ${MCPAN}/monq-out_refp.fa -I ${MCPAN}/aligned/new_mexico/chr${i}/${id}_aug_chr${i}_surject.sorted.marked.bam  -O ${MCPAN}/called/snp/new_mexico/chr${i}/${id}_chr${i}.haplotypecaller.vcf.gz --tmp-dir . --native-pair-hmm-threads ${N}
  done
done


## Concatenate variants of individual chromosomes per sample (for SLiM simulations)
module load biocontainers
module load bcftools
module load htslib
module load vcftools

cd ${MCPAN}/called/snp/new_mexico

for id in SRR11514058 SRR11514059 SRR11514060 SRR11514061 SRR11514063 SRR11514064 SRR11514065; do
  echo "processing ${id}"
  for j in {1..30}; do 
    bcftools sort ./chr${j}/${id}_chr${j}.haplotypecaller.vcf.gz -Oz -o ./chr${j}/${id}_chr${j}_sorted_mcSNP.vcf.gz
    echo "sorting chr${j} done"
  done
  
  files=()
  for j in {1..30}; do 
    vcf=./chr${j}/${id}_chr${j}_sorted_mcSNP.vcf.gz
    [[ -s "$vcf" ]] && files+=("$vcf")
  done
  
  bcftools concat ${files[@]} --threads ${N} -Oz -o ${id}_concat_mcSNP.vcf.gz
  bcftools index ${id}_concat_mcSNP.vcf.gz --threads ${N}
done

## Compressing and indexing each vcf file first
for id in SRR11514058 SRR11514059 SRR11514060 SRR11514061 SRR11514063 SRR11514064 SRR11514065; do
  gunzip ${id}_concat_mcSNP.vcf.gz 
  sed 's/ref_p#0#//g' ${id}_concat_mcSNP.vcf > ${id}_mcSNP_temp.vcf # polish the header line to be compatible with bcftools
  bcftools sort ${id}_mcSNP_temp.vcf -Oz -o ${id}_mcSNP.sorted.vcf.gz
  bcftools index ${id}_mcSNP.sorted.vcf.gz --threads ${N}
  bcftools view ${id}_mcSNP.sorted.vcf.gz --threads ${N} | grep -v "##" | wc -l # number of variants
  rm ${id}_mcSNP_temp.vcf
done

## Combining separately called SNP vcf files of NM only
bcftools merge -m all -Oz -o monqNM_mcSNP.merged.vcf.gz --threads ${N} SRR*_mcSNP.sorted.vcf.gz 

## Filtering with population-level parameters
vcftools --gzvcf monqNM_mcSNP.merged.vcf.gz --missing-indv 
vcftools --gzvcf monqNM_mcSNP.merged.vcf.gz --missing-site
vcftools --gzvcf monqNM_mcSNP.merged.vcf.gz --depth
vcftools --gzvcf monqNM_mcSNP.merged.vcf.gz --site-mean-depth
vcftools --gzvcf monqNM_mcSNP.merged.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.2' | cut -f1,2 >> badloci
bcftools stats monqNM_mcSNP.merged.vcf.gz > vcf-stats_snp.txt

# In R, plotting and summarizing vcf stats
module load r

R
library(tidyverse)
library(ggplot2)
var_qual <- read_delim("./out.lqual", delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
a + theme_light() + xlim(0,500)
summary(var_qual$qual)
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    30.0   108.6   182.6   229.2   293.0 99318.1

var_depth <- read_delim("./out.ldepth.mean", delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()   
a + theme_light() + xlim(0,100)        
summary(var_depth$mean_depth)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#   1.000    5.000    6.828    7.340    8.667 1389.000
quantile(var_depth$mean_depth, probs = seq(0, 1, 1/10))
#        0%        10%        20%        30%        40%        50%        60% 
#   1.00000    3.25000    4.33333    5.00000    6.00000    6.82812    7.50000
#       70%        80%        90%       100%
#   8.06897    9.00000   11.00000 1389.00000
quantile(var_depth$mean_depth, probs = c(0.01, 0.05, 0.95, 0.99))
#     1%      5%     95%     99%
# 2.0000  2.7500 12.6667 20.00007
 
var_miss <- read_delim("./out.lmiss", delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(var_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(var_miss$fmiss)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.01099 0.91209 0.97802 0.88939 0.98901 0.98901
                   
ind_depth <- read_delim("./out.idepth", delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)
a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(ind_depth$depth)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  1.988   7.013   7.480   7.576   7.952  11.343

ind_miss  <- read_delim("./out.imiss", delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(ind_miss$fmiss)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.8537  0.8782  0.8874  0.8894  0.8970  1.0000

quit()

# Filtering vcf file based on the determined conservative cut-off values from above results   
bcftools query -f '%CHROM\t%POS\n' ${CALLED}/monq_mcSNP.final.sort.vcf > vcf_regions.txt # querying the sites of the main vcf file
bcftools norm -m -any --threads ${N} monqNM_mcSNP.merged.vcf.gz -Ov -o monqNM_mcSNP.decomp.vcf
bcftools sort monqNM_mcSNP.decomp.vcf -Oz -o monqNM_mcSNP.decomp.sort.vcf.gz
bcftools index monqNM_mcSNP.decomp.sort.vcf.gz
bcftools view -R vcf_regions.txt monqNM_mcSNP.decomp.sort.vcf.gz -Oz -o monqNM_mcSNP.subset.vcf.gz # make compatible with the main vcf 

# Combining separately called SNP vcf files again for all individuals including NM and of this study
#bcftools merge -m all -Oz -o monqNM_mcSV.merged.vcf.gz --threads ${N} *_mcSV.sorted.vcf.gz 
bcftools merge -m all -Oz -o monq_mcSNP_wNM.merged.vcf.gz --threads ${N} *_mcSNP.sorted.vcf.gz 

# Filtering vcf file based on the determined conservative cut-off values from above results       
bcftools query -f '%CHROM\t%POS\n' ${CALLED}/monq_mcSNP.final.sort.vcf > vcf_regions.txt # querying the sites of the main vcf file
bcftools norm -m -any --threads ${N} monq_mcSNP_wNM.merged.vcf.gz -Ov -o monq_mcSNP_wNM.decomp.vcf
bcftools sort monq_mcSNP_wNM.decomp.vcf -Oz -o monq_mcSNP_wNM.decomp.sort.vcf.gz
bcftools index monq_mcSNP_wNM.decomp.sort.vcf.gz
bcftools view -R vcf_regions.txt monq_mcSNP_wNM.decomp.sort.vcf.gz -Oz -o monq_mcSNP_wNM.subset.vcf.gz 



# Preparation for genetic rescue simulations
ml biocontainers
ml bcftools
ml htslib
ml ensembl-vep

cd ${MCPAN}/called/

## Generating dummy vcf and fasta files
#bcftools view -Ov -o monq_mcSV_wNM.subset.vcf monq_mcSV_wNM.subset.vcf.gz
tabix -p vcf monq_mcSNP_wNM.subset.vcf.gz # make sure the files were moved in the directory
tabix -p vcf monq_mcSV_wNM.subset.vcf.gz # make sure the files were moved in the directory
bcftools query -l monq_mcSNP_wNM.subset.vcf.gz > samples.snp
bcftools query -l monq_mcSV_wNM.subset.vcf.gz  > samples.sv
diff samples.snp samples.sv # check sample order

bcftools view -h monq_mcSNP_wNM.subset.vcf.gz | grep "^##contig" > contig.snp
bcftools view -h monq_mcSV_wNM.subset.vcf.gz  | grep "^##contig" > contig.sv
diff contig.snp contig.sv 3 check chromosome order

bcftools isec -n=2 monq_mcSNP_wNM.subset.vcf.gz monq_mcSV_wNM.subset.vcf.gz # check overlaps -> there are some
bcftools query -f '%CHROM\t%POS\n' monq_mcSV_wNM.subset.vcf.gz | sort -n -u > sv_sites.txt
bcftools view -T ^sv_sites.txt -Oz -o monq_mcSNP_wNM_noisect.vcf.gz monq_mcSNP_wNM.subset.vcf.gz
tabix -p vcf monq_mcSNP_wNM_noisect.vcf.gz

bcftools concat -a -Oz -o monq_mcSNPSV_wNM.subset.vcf.gz monq_mcSNP_wNM_noisect.vcf.gz monq_mcSV_wNM.subset.vcf.gz # concat SNP and SV vcf files
tabix -p vcf monq_mcSNPSV_wNM.subset.vcf.gz

bcftools view -Ov -o monq_mcSNPSV_wNM.subset.vcf monq_mcSNPSV_wNM.subset.vcf.gz
python3 ${APP}/recode_gt_ge2_to1.py monq_mcSNPSV_wNM.subset.vcf > monq_mcSNPSV_wNM_dummybi.vcf # SLiM does not accept multi allelic loci (replace genotypes of 2 to 1)
awk 'BEGIN{OFS="\t"} /^#/ {print; next} {$4 = "A"; $5 = "C"; print}' monq_mcSNPSV_wNM_dummybi.vcf > monq_mcSNPSV_wNM_dummy1.vcf # make all REF alleles to "A" and ALT alleles to "C" 
awk 'BEGIN{OFS="\t"} FNR==NR {if ($1!~/^#/) {key = $1":"$2; keep[key] = 1} next} /^#/ {print; next} {key = $1":"$2; if(key in keep) {$4 = "G"}; print}' monq_SVoutliers_final.vcf monq_mcSNPSV_wNM_dummy1.vcf > monq_mcSNPSV_wNM_dummy2.vcf # make REF alleles of outlier SVs to "G"
bcftools +setGT monq_mcSNPSV_wNM_dummy2.vcf -Ov -o monq_mcSNPSV_wNM_dummy3.vcf -- -t . -n 0 # SLiM does not accept missing genotypes
bcftools view -H monq_mcSNPSV_wNM_dummy3.vcf | head -n 5 # check manually
bcftools view -H monq_mcSNPSV_wNM_dummy2.vcf | head -n 5 # check manually
#awk 'BEGIN{OFS="\t"; pos=0} /^#/ {print; next} {pos++; $1 = "concat"; $2 = pos; print}' monq_mcSV_wNM_dummy3.vcf > monq_mcSV_wNM_forSlim.vcf # reorder position to be 1, 2, 3, ...

#python3 ${APP}/concat_fasta_fromVCF.py monq_mcSV_wNM.subset.vcf monq_SVoutliers_final.vcf concat_forSlim # generate a dummy fasta for SLiM
python3 ${APP}/make_dummy_fasta_fromVCF_forSLiM.py ${MCPAN}/monq-out_refp_renamed.fa.fai monq_mcSNPSV_wNM_dummy3.vcf dummy_fasta_forSLiM # generate a dummy fasta for SLiM
python3 ${APP}/add_T_at_segment_ends.py dummy_fasta_forSLiM.fa dummy_fasta_forSLiM_segments.tsv dummy_fasta_withT_forSlim.fa
cat dummy_fasta_withT_forSlim.fa | grep -v ">" | grep "A" | wc -l
cat dummy_fasta_withT_forSlim.fa | grep -v ">" | grep "G" | wc -l

## Subsetting by samples to use in SLiM
bgzip monq_mcSNPSV_wNM_dummy3.vcf
tabix monq_mcSNPSV_wNM_dummy3.vcf.gz
bcftools view -s E9067,E9570,F1306,F1953,F1930,F1933,F1934,F1955,F1957,F1958 -Oz -o monqCTX_mcSNPSV_forSlim.vcf.gz monq_mcSNPSV_wNM_dummy3.vcf.gz # C TX
bcftools view -s F1876,F1878,F1879,F1882,F1916,F1819,F1825,F1827,F1828,F1832,F1833,F1835,F1910 -Oz -o monqSEMX_mcSNPSV_forSlim.vcf.gz monq_mcSNPSV_wNM_dummy3.vcf.gz # SE MX
bcftools view -s SRR11514058,SRR11514059,SRR11514060,SRR11514061,SRR11514063,SRR11514064,SRR11514065 -Oz -o monqNM_mcSNPSV_forSlim.vcf.gz monq_mcSNPSV_wNM_dummy3.vcf.gz # NM
#bcftools index monqCTX_mcSNPSV_forSlim.vcf.gz 
#bcftools index monqNM_mcSNPSV_forSlim.vcf.gz 
#bcftools index monqSEMX_mcSNPSV_forSlim.vcf.gz 
for f in monqCTX_mcSNPSV_forSlim.vcf.gz monqNM_mcSNPSV_forSlim.vcf.gz monqSEMX_mcSNPSV_forSlim.vcf.gz; do
  base=${f%.vcf.gz}
  bcftools annotate -x INFO/AT,FORMAT/AD,FORMAT/GL -Oz -o ${base}.cleaned.vcf.gz $f
  bcftools index ${base}.cleaned.vcf.gz
done
bcftools merge -m all monqCTX_mcSNPSV_forSlim.cleaned.vcf.gz monqNM_mcSNPSV_forSlim.cleaned.vcf.gz monqSEMX_mcSNPSV_forSlim.cleaned.vcf.gz -Ov -o monq_mcSNPSV_forSlim.cleaned.vcf # merge again to load in SLiM with a known order of populations
#bcftools merge -m all monqCTX_mcSNPSV_forSlim.vcf.gz monqNM_mcSNPSV_forSlim.vcf.gz monqSEMX_mcSNPSV_forSlim.vcf.gz -Ov -o monq_mcSNPSV_forSlim.vcf # merge again to load in SLiM with a known order of populations

cp monq*_mcSNPSV_forSlim.vcf ../../slim/
cp dummy_fasta_withT_forSlim.fa ../../slim/

## Subsetting by samples to estimate genetic load to use in SLiM
bcftools sort -Oz -o monq_mcSNPSV_wNM.subset.sorted.vcf.gz monq_mcSNPSV_wNM.subset.vcf.gz
tabix -f -p vcf monq_mcSNPSV_wNM.subset.sorted.vcf.gz

#bcftools view -s E9067,E9570,F1306,F1953,F1930,F1933,F1934,F1955,F1957,F1958 -Ov -o monqCTX_mcSNPSV_forVEP.vcf monq_mcSNPSV_wNM.subset.sorted.vcf.gz
#bcftools view -s F1876,F1878,F1879,F1882,F1916,F1819,F1825,F1827,F1828,F1832,F1833,F1835,F1910 -Ov -o monqSEMX_mcSNPSV_forVEP.vcf monq_mcSNPSV_wNM.subset.sorted.vcf.gz
#bcftools view -s SRR11514058,SRR11514059,SRR11514060,SRR11514061,SRR11514063,SRR11514064,SRR11514065 -Ov -o monqNM_mcSNPSV_forVEP.vcf monq_mcSNPSV_wNM.subset.sorted.vcf.gz

## prep gff
sed 's/\r$//' ${MCPAN}/jaqu_plus_chicChW_renamed.gff > tmp.clean.gff
grep '^#' tmp.clean.gff > header.txt
grep -v '^#' tmp.clean.gff | \
    LC_ALL=C sort -t $'\t' -k1,1 -k4,4n > body.sorted.gff
cat header.txt body.sorted.gff > jaqu_plus_chicChW_renamed.sorted.gff
bgzip -f jaqu_plus_chicChW_renamed.sorted.gff
tabix -f -p gff jaqu_plus_chicChW_renamed.sorted.gff.gz

rm -f tmp.clean.gff header.txt body.sorted.gff

#for INPUT_VCF in monqCTX_mcSNPSV_forVEP.vcf monqSEMX_mcSNPSV_forVEP.vcf monqNM_mcSNPSV_forVEP.vcf; do
#  PREFIX=$(basename "$INPUT_VCF" _mcSNPSV_forVEP.vcf)
#  vep --input_file ${INPUT_VCF} --output_file ${PREFIX}.vep.sift --force_overwrite --fasta ${MCPAN}/monq-out_refp_renamed.fa --gff jaqu_plus_chicChW_renamed.sorted.gff.gz \
#      --everything --variant_class --species custom --assembly custom --fork ${N} --allow_non_variant
#done
vep --input_file monq_mcSNPSV_wNM.subset.sorted.vcf.gz --vcf --output_file monq_mcSNPSV_wNM.vep.vcf --force_overwrite --fasta ${MCPAN}/monq-out_refp_renamed.fa --gff jaqu_plus_chicChW_renamed.sorted.gff.gz \
    --everything --variant_class --species custom --assembly custom --fork ${N} --allow_non_variant
vep --input_file monq_SVoutliers_wNM.subset.vcf.gz --vcf --output_file monq_SVoutliers_wNM.vep.vcf --force_overwrite --fasta ${MCPAN}/monq-out_refp_renamed.fa --gff jaqu_plus_chicChW_renamed.sorted.gff.gz \
    --everything --variant_class --species custom --assembly custom --fork ${N} --allow_non_variant

#bcftools sort monq_mcSNP_wNM.subset.vcf.gz -Oz -o monq_mcSNP_wNM.subset.sorted.vcf.gz
#bcftools sort monq_mcSV_wNM.subset.vcf.gz -Oz -o monq_mcSV_wNM.subset.sorted.vcf.gz
#vep --input_file monq_mcSNP_wNM.subset.sorted.vcf.gz --vcf --output_file monq_mcSNP_wNM.vep.vcf --force_overwrite --fasta ${MCPAN}/monq-out_refp_renamed.fa --gff jaqu_plus_chicChW_renamed.sorted.gff.gz \
#    --everything --variant_class --species custom --assembly custom --fork ${N} --allow_non_variant
#vep --input_file monq_mcSV_wNM.subset.sorted.vcf.gz --vcf --output_file monq_mcSV_wNM.vep.vcf --force_overwrite --fasta ${MCPAN}/monq-out_refp_renamed.fa --gff jaqu_plus_chicChW_renamed.sorted.gff.gz \
#    --everything --variant_class --species custom --assembly custom --fork ${N} --allow_non_variant


# per-individual allele counts by worst IMPACT
VCF=monq_mcSNPSV_wNM.vep.vcf

## split multi-allelic sites so ALT is a single allele per record
bcftools query -f '%ALT\t%INFO/CSQ[\t%GT]\n' "$VCF" \
| awk -F'\t' -v samples="$(bcftools query -l "$VCF" | paste -sd, -)" '
BEGIN{
  n=split(samples,S,",");
  sev["MODIFIER"]=1; sev["LOW"]=2; sev["MODERATE"]=3; sev["HIGH"]=4;
  print "Sample\tHIGH_impact\tMODERATE_impact\tLOW_impact\tMODIFIER";
}
function worst_impact_for_alt(alt, csq,   m,j,f,worst,wsev){
  worst=""; wsev=0
  m=split(csq,A,",")
  for(j=1;j<=m;j++){
    split(A[j],f,"|")
    if(f[1]==alt && (f[3] in sev) && sev[f[3]]>wsev){ wsev=sev[f[3]]; worst=f[3] }
  }
  return worst
}
function worst_impact_for_gt(gt, alts, csq,   a,k,m,ai,alt,imp,worst,wsev){
  if(gt=="./."||gt==".|.") return ""
  split(alts,ALT,",")
  m=split(gt,a,/[\/|]/)
  worst=""; wsev=0
  for(k=1;k<=m;k++){
    ai=a[k]; if(ai==""||ai=="."||ai==0) continue
    alt=ALT[ai]
    imp=worst_impact_for_alt(alt,csq)
    if(imp!="" && sev[imp]>wsev){ wsev=sev[imp]; worst=imp }
  }
  return worst
}
{
  alts=$1; csq=$2
  for(i=3;i<=NF;i++){
    imp=worst_impact_for_gt($i,alts,csq)
    if(imp!="") C[i,imp]++
  }
}
END{
  for(i=1;i<=n;i++){
    col=i+2
    printf "%s\t%d\t%d\t%d\t%d\n", S[i], C[col,"HIGH"]+0, C[col,"MODERATE"]+0, C[col,"LOW"]+0, C[col,"MODIFIER"]+0
  }
}'

VCF=monq_SVoutliers_wNM.vep.vcf

## split multi-allelic sites so ALT is a single allele per record
bcftools query -f '%ALT\t%INFO/CSQ[\t%GT]\n' "$VCF" \
| awk -F'\t' -v samples="$(bcftools query -l "$VCF" | paste -sd, -)" '
BEGIN{
  n=split(samples,S,",");
  sev["MODIFIER"]=1; sev["LOW"]=2; sev["MODERATE"]=3; sev["HIGH"]=4;
  print "Sample\tHIGH_impact\tMODERATE_impact\tLOW_impact\tMODIFIER";
}
function worst_impact_for_alt(alt, csq,   m,j,f,worst,wsev){
  worst=""; wsev=0
  m=split(csq,A,",")
  for(j=1;j<=m;j++){
    split(A[j],f,"|")
    if(f[1]==alt && (f[3] in sev) && sev[f[3]]>wsev){ wsev=sev[f[3]]; worst=f[3] }
  }
  return worst
}
function worst_impact_for_gt(gt, alts, csq,   a,k,m,ai,alt,imp,worst,wsev){
  if(gt=="./."||gt==".|.") return ""
  split(alts,ALT,",")
  m=split(gt,a,/[\/|]/)
  worst=""; wsev=0
  for(k=1;k<=m;k++){
    ai=a[k]; if(ai==""||ai=="."||ai==0) continue
    alt=ALT[ai]
    imp=worst_impact_for_alt(alt,csq)
    if(imp!="" && sev[imp]>wsev){ wsev=sev[imp]; worst=imp }
  }
  return worst
}
{
  alts=$1; csq=$2
  for(i=3;i<=NF;i++){
    imp=worst_impact_for_gt($i,alts,csq)
    if(imp!="") C[i,imp]++
  }
}
END{
  for(i=1;i<=n;i++){
    col=i+2
    printf "%s\t%d\t%d\t%d\t%d\n", S[i], C[col,"HIGH"]+0, C[col,"MODERATE"]+0, C[col,"LOW"]+0, C[col,"MODIFIER"]+0
  }
}'


# Subsetting SNPs by samples, too, for Ne estimation
bcftools view -m2 -M2 -v snps -Oz -o monq_mcSNP.biallelic.vcf.gz monq_mcSNP.final.sort.vcf
bcftools view -m2 -M2 -v snps -Oz -o monqNM_mcSNP.biallelic.vcf.gz monqNM_mcSNP.subset.vcf.gz
bcftools view -s E9067,E9570,F1306,F1953,F1930,F1933,F1934,F1955,F1957,F1958 -Ov -o monqCTX_mcSNP.vcf monq_mcSNP.biallelic.vcf.gz # C TX
bcftools view -s F1876,F1878,F1879,F1882,F1916,F1819,F1825,F1827,F1828,F1832,F1833,F1835,F1910 -Ov -o monqSEMX_mcSNP.vcf monq_mcSNP.biallelic.vcf.gz # SE MX
bcftools view -Ov -o monqNM_mcSNP.vcf monqNM_mcSNP.biallelic.vcf.gz # NM

bcftools view -s E9067,E9570,F1306,F1953,F1930,F1933,F1934,F1955,F1957,F1958,E8954,E9030,E9031,E9032,E9035,E9037,E9040,E9042,E9044,E9046,E9047,E9054,E9055,E9056,E9058 -Ov -o monqTX_mcSNP.vcf monq_mcSNP.biallelic.vcf.gz # TX
bcftools view -s F1788,F1789,F1792,F1793,F1800,F1801,F1840,F1852,F1853,F1857,F1858,F1876,F1878,F1879,F1882,F1916,F1819,F1825,F1827,F1828,F1832,F1833,F1835,F1910 -Ov -o monqMXCT_mcSNP.vcf monq_mcSNP.biallelic.vcf.gz
bcftools view -s F1811,F1814,F1815,F1872,F1873,F1874,F1888,F1889,F1890,F1891,F1895,F1897,F1900,F1901,F1902,F1903,F1904,F1911,F1913,F1914,F1915,F1940,F1941,F1942,F1944,F1946,F1947,F1949,F1952 -Ov -o monqSN_mcSNP.vcf monq_mcSNP.biallelic.vcf.gz # Sonora

bgzip monqSN_mcSNP.vcf
bgzip monqNM_mcSNP.vcf
tabix monqSN_mcSNP.vcf.gz
tabix monqNM_mcSNP.vcf.gz
bcftools merge -m all -Ov -o monqNMSN_mcSNP.vcf --threads ${N} monqSN_mcSNP.vcf.gz monqNM_mcSNP.vcf.gz


# Estimate Ne
${APP}/currentNe2/currentne2 -t ${N} monqCTX_mcSNP.vcf 30 # C TX -> not converge; 
${APP}/currentNe2/currentne2 -t ${N} monqTX_mcSNP.vcf 30 # TX -> 132.40 (89.55 - 195.77; 90% CI); half for CTX? 660
${APP}/currentNe2/currentne2 -t ${N} monqSEMX_mcSNP.vcf 30 # SE MX -> 433.57 (155.81 - 1206.46; 90% CI)
${APP}/currentNe2/currentne2 -t ${N} monqMXCT_mcSNP.vcf 30 # MX City -> 225.55 (140.30 - 362.60; 90% CI)
${APP}/currentNe2/currentne2 -t ${N} monqNM_mcSNP.vcf 30 # NM -> not converge; 36,417 from Mathur & DeWoody 2021
${APP}/currentNe2/currentne2 -t ${N} monqSN_mcSNP.vcf 30 # Sonora -> 314.88 (196.35 - 504.98; 90% CI)
${APP}/currentNe2/currentne2 -t ${N} monqNMSN_mcSNP.vcf 30 # NM-SN -> 336.53 (217.82 - 519.92; 90% CI)
