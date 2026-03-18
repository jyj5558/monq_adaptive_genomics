#!/bin/bash
#SBATCH --job-name=pangenome
#SBATCH -A fnrdewoody
#SBATCH -t 5-00:00:00
#SBATCH -n 128
#SBATCH -N 1
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL

# long read assemblies were first renamed and masked following https://github.com/nataliemallen/Pangenome/tree/main
# then pangenome assembly with minigraph cactus:

module purge
module load anaconda
conda activate cactus_env
source /home/allen715/cactus/cactus-bin-v2.7.0/venv-cactus-v2.7.0/bin/activate

REF=refp
OUT=/scratch/negishi/allen715/MONQ/longread_assem/pangenome/out
REFDIR=/scratch/negishi/allen715/MONQ/longread_assem/pangenome/ref
GENOME=reference_primary_renamed

cactus-minigraph ./jobstore seqfile.txt monq_pan.gfa --reference $REF --binariesMode local

cactus-graphmap ./jobstore seqfile.txt monq_pan.gfa monq_pan.paf --outputFasta monq_pan.gfa.fa --reference $REF --binariesMode local --delFilter 10000000

keep_contigs=$(awk '{print $1}' ${REFDIR}/${GENOME}.fasta.fai)

cactus-graphmap-split ./jobstore seqfile.txt monq_pan.gfa monq_pan.paf --reference $REF --outDir $OUT --binariesMode local --otherContig contigOther --refContigs $(for i in $keep_contigs; do echo ${i}; done)

cactus-align ./jobstore $OUT/chromfile.txt $OUT/chrom-alignments --batch --pangenome --reference $REF --outVG --maxLen 10000 --binariesMode local

cactus-graphmap-join ./jobstore --vg $OUT/chrom-alignments/*.vg --hal $OUT/chrom-alignments/*.hal --outDir $OUT --outName monq_pan-out --reference $REF --binariesMode local --gbz clip full --gfa --filter 2 --clip 10000 --giraffe clip --vcf --chrom-vg --chrom-og --viz --xg --draw
