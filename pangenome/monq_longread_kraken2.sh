#!/bin/bash
#SBATCH --job-name=F1952
#SBATCH -A fnrblack
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 14-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=jeon96@purdue.edu
#SBATCH --mail-type=END,FAIL

# paths
BASE_DIR=/scratch/bell/jeon96/MONQ/lr-assembly
WORK_DIR=/scratch/bell/jeon96/MONQ/lr-assembly/filt

export PATH=$HOME/app/kraken2-2.1.5:$PATH

cd "$WORK_DIR"

for fa_file in F_1847.asm.p_ctg.purged.fa F_1847.asm.a_ctg.purged.fa; do
    echo "$fa_file started"
    base=$(basename "$fa_file" .purged.fa)
    
    # extract ID 
    id=$(echo "$base" | cut -d. -f1)
    
    k2 classify --db $HOME/app/kraken2-2.1.5/pluspf_08gb --threads 64 --unclassified-out "${base}_filt.fa" --classified-out "${base}_cont.fa" --output "${base}_cseq" --report "${base}_kraken" --confidence 0.3 --minimum-hit-groups 2 --use-names "${BASE_DIR}/$fa_file"

    echo "Finished processing $base"
done