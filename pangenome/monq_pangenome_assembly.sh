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

### step 1: prep long read assemblies 

# find all files ending in *filt.fa.gz from the current directory
find . -type f -name "*filt.fa.gz" | while read -r file_path; do
    # check if the file exists (redundant, but safe)
    if [ -e "$file_path" ]; then
        # create a temporary file to store the modified content
        temp_file="${file_path}.tmp"

        # counter for contig renaming
        count=1

        # unzip, modify, and write to temp
        gunzip -c "$file_path" | while IFS= read -r line; do
            if [[ $line == ">"* ]]; then
                echo ">A1_pri$count" >> "$temp_file"
                ((count++))
            else
                echo "$line" >> "$temp_file"
            fi
        done

        # gzip the temp file and overwrite the original
        gzip -c "$temp_file" > "$file_path"

        # remove the temp file
        rm "$temp_file"

        echo "Contigs in $file_path renamed successfully."
    else
        echo "Fasta file $file_path not found."
    fi
done

### step 2: prep reference genome

#renames fasta headings of reference files. primary should be "chr", alt should be "contig" (must change in script)

cd /scratch/negishi/allen715/MONQ/longread_assem/pangenome/

##primary
# Define input file
genome_file="monq_ref_p.fa" 
# Define output file
genome_output="monq_ref_p_renamed.fasta" 

# Function to rename chromosomes in a genome file
rename_chromosomes() {
    input_file=$1
    output_file=$2

    # Use awk to rename chromosomes
    awk '/^>/{print ">chr" ++i; next}{print}' < "$input_file" > "$output_file"
}

# Rename chromosomes in the genome
rename_chromosomes "$genome_file" "$genome_output"

##alt
# Define input file
genome_file2="monq_ref_a.fa" 
# Define output file
genome_output2="monq_ref_a_renamed.fasta" 

# Function to rename chromosomes in a genome file
rename_chromosomes() {
    input_file=$1
    output_file=$2

    # Use awk to rename chromosomes
    awk '/^>/{print ">contig" ++i; next}{print}' < "$input_file" > "$output_file"
}

# Rename chromosomes in the genome
rename_chromosomes "$genome_file2" "$genome_output2"

### step 3: index reference genome

module load biocontainers
module load samtools

samtools faidx monq_ref_p_renamed.fasta

samtools faidx monq_ref_a_renamed.fasta

### step 4: mask long read assemblies

module purge
module load biocontainers
module load repeatmasker
module load bedtools
export PATH=$PATH:/home/allen715/winmasker/ncbi_cxx--25_2_0/GCC850-ReleaseMT64/bin/

# loop through all *_filt.fa.gz files
for file in *_filt.fa.gz; do
    # get base name without path and extension
    base=$(basename "$file" .fa.gz)     # e.g., F_1789_p_filt
    prefix=${base%_filt}                # e.g., F_1789_p

    echo "Processing $file"

    # make output directory
    mkdir -p "$prefix"
    cd "$prefix" || exit 1

    # decompress and copy fasta into the working folder
    cp ../"$file" .
    gunzip "$file"   # now have $base.fa

    # run windowmasker
    windowmasker -mk_counts -in "$base.fa" -out "stage1_${prefix}.counts"
    windowmasker -ustat "stage1_${prefix}.counts" -in "$base.fa" -outfmt fasta -out "${prefix}_winmask.fasta"

    # run RepeatMasker (writes output in current dir)
    RepeatMasker -pa 32 -xsmall -species aves "$base.fa"

    # convert .out file to BED format
    awk 'BEGIN{OFS="\t"} {if($1 !~ /^#/) print $5, $6-1, $7, $11, $0}' "${base.fa}.out" | tail -n +4 > output.bed

    # apply final soft masking
    bedtools maskfasta -soft -fi "${prefix}_winmask.fasta" -bed output.bed -fo "${prefix}_masked.fasta"

    # return to starting directory
    cd ..
done

### step 5: pangenome assembly

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
