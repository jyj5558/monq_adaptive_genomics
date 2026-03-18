#!/bin/bash
#SBATCH --job-name=hifiasm
#SBATCH -A highmem
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 1-00:00:00
#SBATCH --error=hifiasm.err
#SBATCH --output=hifiasm.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=allen715@purdue.edu

#### Step 1: first assemble genomes in hifiasm 

module --force purge
module load biocontainers
module load hifiasm

cd /scratch/negishi/allen715/

hifiasm -o E_9056.asm --primary -t 64 E_9056.fastq.gz
hifiasm -o F_1789.asm --primary -t 64 F_1789.tfastq.gz
hifiasm -o F_1803.asm --primary -t 64 F_1803.fastq.gz
hifiasm -o F_1847.asm --primary -t 64 F_1847.fastq.gz
hifiasm -o F_1865.asm --primary -t 64 F_1865.fastq.gz
hifiasm -o F_1867.asm --primary -t 64 F_1867.fastq.gz
hifiasm -o F_1915.asm --primary -t 64 F_1915.fastq.gz
hifiasm -o F_1927.asm --primary -t 64 F_1927.fastq.gz
hifiasm -o F_1934.asm --primary -t 64 F_1934.fastq.gz
hifiasm -o F_1952.asm --primary -t 64 F_1952.tfastq.gz

#### Step 2: run jellyfish

module load biocontainers
module load kmer-jellyfish

cd /scratch/negishi/allen715/

jellyfish count -m 21 -s 100M -t 64 -C -o E_9056_p.jf E_9056_p.fa && jellyfish histo -t 64 E_9056_p.jf > E_9056_p_histo.txt
jellyfish count -m 21 -s 100M -t 64 -C -o E_9056_a.jf E_9056_a.fa && jellyfish histo -t 64 E_9056_a.jf > E_9056_a_histo.txt
jellyfish count -m 21 -s 100M -t 64 -C -o F_1789_p.jf F_1789_p.fa && jellyfish histo -t 64 F_1789_p.jf > F_1789_p_histo.txt
jellyfish count -m 21 -s 100M -t 64 -C -o F_1789_a.jf F_1789_a.fa && jellyfish histo -t 64 F_1789_a.jf > F_1789_a_histo.txt
jellyfish count -m 21 -s 100M -t 64 -C -o F_1803_p.jf F_1803_p.fa && jellyfish histo -t 64 F_1803_p.jf > F_1803_p_histo.txt
jellyfish count -m 21 -s 100M -t 64 -C -o F_1803_a.jf F_1803_a.fa && jellyfish histo -t 64 F_1803_a.jf > F_1803_a_histo.txt
jellyfish count -m 21 -s 100M -t 64 -C -o F_1847_p.jf F_1847_p.fa && jellyfish histo -t 64 F_1847_p.jf > F_1847_p_histo.txt
jellyfish count -m 21 -s 100M -t 64 -C -o F_1847_a.jf F_1847_a.fa && jellyfish histo -t 64 F_1847_a.jf > F_1847_a_histo.txt
jellyfish count -m 21 -s 100M -t 64 -C -o F_1865_p.jf F_1865_p.fa && jellyfish histo -t 64 F_1865_p.jf > F_1865_p_histo.txt
jellyfish count -m 21 -s 100M -t 64 -C -o F_1865_a.jf F_1865_a.fa && jellyfish histo -t 64 F_1865_a.jf > F_1865_a_histo.txt
jellyfish count -m 21 -s 100M -t 64 -C -o F_1867_p.jf F_1867_p.fa && jellyfish histo -t 64 F_1867_p.jf > F_1867_p_histo.txt
jellyfish count -m 21 -s 100M -t 64 -C -o F_1867_a.jf F_1867_a.fa && jellyfish histo -t 64 F_1867_a.jf > F_1867_a_histo.txt
jellyfish count -m 21 -s 100M -t 64 -C -o F_1915_p.jf F_1915_p.fa && jellyfish histo -t 64 F_1915_p.jf > F_1915_p_histo.txt
jellyfish count -m 21 -s 100M -t 64 -C -o F_1915_a.jf F_1915_a.fa && jellyfish histo -t 64 F_1915_a.jf > F_1915_a_histo.txt
jellyfish count -m 21 -s 100M -t 64 -C -o F_1927_p.jf F_1927_p.fa && jellyfish histo -t 64 F_1927_p.jf > F_1927_p_histo.txt
jellyfish count -m 21 -s 100M -t 64 -C -o F_1927_a.jf F_1927_a.fa && jellyfish histo -t 64 F_1927_a.jf > F_1927_a_histo.txt
jellyfish count -m 21 -s 100M -t 64 -C -o F_1934_p.jf F_1934_p.fa && jellyfish histo -t 64 F_1934_p.jf > F_1934_p_histo.txt
jellyfish count -m 21 -s 100M -t 64 -C -o F_1934_a.jf F_1934_a.fa && jellyfish histo -t 64 F_1934_a.jf > F_1934_a_histo.txt
jellyfish count -m 21 -s 100M -t 64 -C -o F_1952_p.jf F_1952_p.fa && jellyfish histo -t 64 F_1952_p.jf > F_1952_p_histo.txt
jellyfish count -m 21 -s 100M -t 64 -C -o F_1952_a.jf F_1952_a.fa && jellyfish histo -t 64 F_1952_a.jf > F_1952_a_histo.txt

# Upload *_histo.txt files to GenomeScope2 (web)

#### Step 3: run jellyfish

module purge
module load biocontainers
module load purge_dups
module load minimap2/2.26

cd /scratch/negishi/allen715/

# generate config files first for each assembly:
# pd_config.py E_9056_p.fa pb.fofn && mv config.json E_9056_p_config.json

run_purge_dups.py -p bash E_9056_p_config.json /usr/local/bin E_9056_p
run_purge_dups.py -p bash E_9056_a_config.json /usr/local/bin E_9056_a
run_purge_dups.py -p bash F_1789_p_config.json /usr/local/bin F_1789_p
run_purge_dups.py -p bash F_1789_a_config.json /usr/local/bin F_1789_a
run_purge_dups.py -p bash F_1803_p_config.json /usr/local/bin F_1803_p
run_purge_dups.py -p bash F_1803_a_config.json /usr/local/bin F_1803_a
run_purge_dups.py -p bash F_1847_p_config.json /usr/local/bin F_1847_p
run_purge_dups.py -p bash F_1847_a_config.json /usr/local/bin F_1847_a
run_purge_dups.py -p bash F_1865_p_config.json /usr/local/bin F_1865_p
run_purge_dups.py -p bash F_1865_a_config.json /usr/local/bin F_1865_a
run_purge_dups.py -p bash F_1867_p_config.json /usr/local/bin F_1867_p
run_purge_dups.py -p bash F_1867_a_config.json /usr/local/bin F_1867_a
run_purge_dups.py -p bash F_1915_p_config.json /usr/local/bin F_1915_p
run_purge_dups.py -p bash F_1915_a_config.json /usr/local/bin F_1915_a
run_purge_dups.py -p bash F_1927_p_config.json /usr/local/bin F_1927_p
run_purge_dups.py -p bash F_1927_a_config.json /usr/local/bin F_1927_a
run_purge_dups.py -p bash F_1934_p_config.json /usr/local/bin F_1934_p
run_purge_dups.py -p bash F_1934_a_config.json /usr/local/bin F_1934_a
run_purge_dups.py -p bash F_1952_p_config.json /usr/local/bin F_1952_p
run_purge_dups.py -p bash F_1952_a_config.json /usr/local/bin F_1952_a
