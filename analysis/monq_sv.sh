#!/bin/bash
#SBATCH --job-name=monq_vg
#SBATCH -A fnrtowhee
#SBATCH -t 14-00:00:00
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=jeon96@purdue.edu
#SBATCH --mail-type=END,FAIL

MCPAN=/scratch/bell/jeon96/monq/pan
APP=/home/jeon96/app
N=64
CLEAN_SRA=/scratch/bell/jeon96/monq/short-read/cleaned
FIXED_SRA=/scratch/bell/jeon96/monq/short-read/fixed
FILT_SRA=/scratch/bell/jeon96/monq/short-read/filtered
FIXED_SRA=/scratch/bell/jeon96/monq/short-read/fixed

module load biocontainers
module load bcftools
module load vcftools
module load htslib

cd ${MCPAN}/called/sv/vg
CALLED=/scratch/bell/jeon96/monq/pan/called/sv/vg/

# Preprocessing
## Compressing and indexing each vcf file first
for id in `ls -lrt ${CALLED}/chr1 | tr -s ' ' | cut -d ' ' -f 9 | cut -d '_' -f 1 | sort | uniq`; do
  #id=$(echo $i | cut -d '_' -f 1)
  gunzip ${id}_concat_mcSV.vcf.gz 
  sed 's/ref_p#0#//g' ${id}_concat_mcSV.vcf > ${id}_mcSV_temp.vcf # polish the header line to be compatible with bcftools
  bcftools sort ${id}_mcSV_temp.vcf -Oz -o ${id}_mcSV.sorted.vcf.gz
  bcftools index ${id}_mcSV.sorted.vcf.gz --threads ${N}
  bcftools view ${id}_mcSV.sorted.vcf.gz --threads ${N} | grep -v "##" | wc -l # number of variants
  rm ${id}_mcSV_temp.vcf
done

## Combining separately called SV vcf files
bcftools merge -m all -Oz -o monq_mcSV.merged.vcf.gz --threads ${N} *_mcSV.sorted.vcf.gz 

## Filtering with population-level parameters
vcftools --gzvcf monq_mcSV.merged.vcf.gz --missing-indv 
vcftools --gzvcf monq_mcSV.merged.vcf.gz --missing-site
vcftools --gzvcf monq_mcSV.merged.vcf.gz --depth
vcftools --gzvcf monq_mcSV.merged.vcf.gz --site-mean-depth
vcftools --gzvcf monq_mcSV.merged.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.2' | cut -f1,2 >> badloci
bcftools stats monq_mcSV.merged.vcf.gz > vcf-stats.txt

## In R, plotting and summarizing vcf stats
module load r

R
library(tidyverse)
library(ggplot2)
var_qual <- read_delim("./out.lqual", delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
a + theme_light() + xlim(0,500)
summary(var_qual$qual)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#    -4.79     15.94     46.18     98.50    113.14 124270.00
   
var_depth <- read_delim("./out.ldepth.mean", delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()   
a + theme_light() + xlim(0,100)        
summary(var_depth$mean_depth)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#   0.0000    0.1209    3.4396    4.3707    7.3736 1978.0400
quantile(var_depth$mean_depth, probs = seq(0, 1, 1/10))
#          0%          10%          20%          30%          40%          50%         60%          70%          80%          90%         100%
#   0.0000000    0.0000000    0.0374074    0.3296700    1.4317960    3.4395600   5.2645060    6.6923100    8.0659300    9.5934100 1978.0400000
quantile(var_depth$mean_depth, probs = c(0.01, 0.05, 0.95, 0.99))
#     1%      5%     95%     99%
# 0.0000  0.0000 10.4648 12.6667
 
var_miss <- read_delim("./out.lmiss", delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(var_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(var_miss$fmiss)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.00000 0.02198 0.07692 0.34406 0.78022 1.00000
                       
ind_depth <- read_delim("./out.idepth", delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)
a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(ind_depth$depth)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.000   3.943   4.302   4.374   4.754   7.110

ind_miss  <- read_delim("./out.imiss", delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(ind_miss$fmiss)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.3142  0.3297  0.3348  0.3441  0.3382  0.9995

quit()

# Filtering vcf file based on the determined conservative cut-off values from above results                                                                      
vcftools --gzvcf monq_mcSV.merged.vcf.gz --out monq_mcSV.filtered --minQ 20 --mac 2 --min-meanDP 3 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all #--remove-filtered lowad --max-meanDP 15 
bcftools sort monq_mcSV.filtered.recode.vcf -Oz -o monq_mcSV.filtered.vcf.gz
${APP}/vcfbub --input monq_mcSV.filtered.vcf.gz --max-allele-length 100000 --max-level 0 > monq_mcSV.filtered.popped.vcf # remove large (>100 kb) alleles in graphs. This removes all non-top-level bubbles from the VCF unless they were nested inside a top-level bubble with a reference length exceeding 100 kb; that is, top-level bubbles longer than that are replaced by their child nodes in the snarl tree.

singularity run ${APP}/vcflib-1.09.simg
N=64
vcfwave -I 1000 -t ${N} monq_mcSV.filtered.popped.vcf > monq_mcSV.filt.pop.decomp.vcf # decompose complex alleles into primitive ones
exit

vcftools --gzvcf monq_mcSV.filt.pop.decomp.vcf --missing-site --out out2
cat out2.lmiss | awk '!/CHR/' | awk '$6 > 0.2' | cut -f1,2 >> badloci2
vcftools --vcf monq_mcSV.filt.pop.decomp.vcf  --out monq_mcSV.final --minQ 20 --mac 2 --max-missing 0.8 --exclude-positions badloci2 --recode --recode-INFO-all


# Profile SVs
GEAGO=/scratch/negishi/jeon96/monq/gea-go
MCPAN=/scratch/negishi/jeon96/monq/pan
APP=/home/jeon96/app
N=64

module load anaconda
module load use.own
module load conda-env/mypackages-py3.9.13 # Negishi
module load biocontainers
module load bcftools
module load vcftools

cd ${GEAGO}

truvari anno svinfo -m 0 -o monq_mcSV.final.annot.vcf monq_mcSV.final.recode.vcf # annotate SV info with the minimum allele size of 50
echo "Type summary" > all_sv_summary.txt
bcftools query -f '%SVTYPE\n' monq_mcSV.final.annot.vcf | sort | uniq -c >> all_sv_summary.txt # sv profiling
bcftools query -f '%SVTYPE\t%SVLEN\n' monq_mcSV.final.annot.vcf > all_sv_length.txt # sv length distribution
echo "Length average" >> all_sv_summary.txt
cat all_sv_length.txt | awk '{sum[$1]+=$2; count[$1]++} END {for (i in sum) print i, sum[i]/count[i]}' >> all_sv_summary.txt 
#bcftools query -f '%CHROM\t%POS\t%SVTYPE\n' monq_mcSV.final.annot.vcf > all_sv_regions.tsv
bcftools query -f '%CHROM\t%SVTYPE\n' monq_mcSV.final.annot.vcf > all_sv_counts.tsv

bcftools sort monq_mcSV.final.annot.vcf > monq_mcSV.final.sort.vcf

# After RDA (see RDA script),
grep "^#" monq_mcSV.final.sort.vcf > monq_SVoutliers_final.vcf
awk 'BEGIN {FS=OFS="\t"} NR==FNR {key[$1 FS $2 FS $3 FS $4 FS $5]; next} ($1 FS $2 FS $3 FS $4 FS $5) in key' GEAoutliers_full_lines_final.tsv monq_mcSV.final.sort.vcf >> monq_SVoutliers_final.vcf

truvari anno svinfo -m 0 -o monq_SVoutliers.annot.vcf monq_SVoutliers_final.vcf # annotate SV info with the minimum allele size of 0
echo "Type summary" > outlier_sv_summary.txt
bcftools query -f '%SVTYPE\n' monq_SVoutliers.annot.vcf | sort | uniq -c >> outlier_sv_summary.txt # sv profiling
bcftools query -f '%SVTYPE\t%SVLEN\n' monq_SVoutliers.annot.vcf > outlier_sv_length.txt # sv length distribution
echo "Length average" >> outlier_sv_summary.txt
cat outlier_sv_length.txt | awk '{sum[$1]+=$2; count[$1]++} END {for (i in sum) print i, sum[i]/count[i]}' >> outlier_sv_summary.txt 
#bcftools query -f '%CHROM\t%POS\t%SVTYPE\n' monq_SVoutliers.annot.vcf >  outlier_sv_regions.tsv
bcftools query -f '%CHROM\t%SVTYPE\n' monq_SVoutliers.annot.vcf > outlier_sv_counts.tsv

R
library(ggplot2)
library(dplyr)
library(gtools)

#regions <- read.table("all_sv_regions.tsv", header = FALSE, col.names = c("CHROM", "POS", "SVTYPE"))
#table(regions$SVTYPE)
#hist(regions$END - regions$START, breaks=100, main="All SV length distribution")

# Read the 'all' data
counts <- read.table("all_sv_counts.tsv", header = FALSE, col.names = c("Chromosome", "SVTYPE"))
counts <- counts[counts$SVTYPE != "SNP",]
counts$Chromosome <- gsub("chr", "", counts$Chromosome)
genome <- read.table("/scratch/negishi/jeon96/monq/pan/monq-out_refp.genome", header = FALSE, col.names = c("Chromosome", "Length"))
genome$Chromosome <- gsub("chr", "", genome$Chromosome)
sv_summary <- counts %>% group_by(Chromosome, SVTYPE) %>% summarise(Count = n()) %>% ungroup()
sv_summary$Chromosome <- factor(sv_summary$Chromosome, levels = mixedsort(unique(sv_summary$Chromosome)))
sv_genome <- left_join(sv_summary, genome, by = "Chromosome") 
sv_genome$StandardizedCount <- sv_genome$Count / sv_genome$Length * 1000000 # number of SVs per Mb
sv_genome$Chromosome <- factor(sv_genome$Chromosome, levels = mixedsort(unique(sv_genome$Chromosome)))

p <- ggplot(sv_summary, aes(x = Chromosome, y = Count, fill = SVTYPE)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Distribution of SV Types per Chromosome", x = "Chromosome", y = "Count") +
  theme_minimal(base_size = 12) +
  theme(axis.ticks = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
  plot.title = element_text(hjust = 0.5))
ggsave("all_sv_distribution.svg", plot = p, width = 8, height = 8, units = "in", device = "svg")
ggsave("all_sv_distribution.pdf", plot = p, width = 8, height = 8, units = "in", device = cairo_pdf)

p <- ggplot(sv_genome, aes(x = Chromosome, y = StandardizedCount, fill = SVTYPE)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Distribution of SV Types per Chromosome", x = "Chromosome", y = "Length-corrected Count") +
  theme_minimal(base_size = 12) +
  theme(axis.ticks = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
  plot.title = element_text(hjust = 0.5))
ggsave("all_sv_standardized_distribution.svg", plot = p, width = 8, height = 8, units = "in", device = "svg")
ggsave("all_sv_standardized_distribution.pdf", plot = p, width = 8, height = 8, units = "in", device = cairo_pdf)

lengths <- read.table("all_sv_length.txt", header = FALSE, col.names = c("SVTYPE", "SVLEN"))
lengths <- lengths[lengths$SVTYPE != "SNP",]
lengths$SVLEN <- abs(lengths$SVLEN)
lengths$LengthBin <- cut(
  lengths$SVLEN,
  breaks = c(0, 50, 1000, 10000, 100000),
  labels = c("0-50", "50-1k", "1k-10k", "10k-100k"),
  right = FALSE
)
sv_length_summary <- lengths %>% group_by(SVTYPE, LengthBin) %>% summarise(Count = n(), .groups = "drop")

p <- ggplot(sv_length_summary, aes(x = LengthBin, y = Count, fill = SVTYPE)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = Count), position = position_dodge(width = 0.9), vjust = -0.3, size = 3) +
  labs(
    title = "SV Length Distribution by Type",
    x = "SV Length Bin (bp)",
    y = "Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggsave("all_sv_lengths.svg", plot = p, width = 8, height = 8, units = "in", device = "svg")
ggsave("all_sv_lengths.pdf", plot = p, width = 8, height = 8, units = "in", device = cairo_pdf)

# re-plot (raw distribution of SVs)
library(VariantAnnotation)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggridges)

vcf_file <- "monq_mcSV.final.sort.vcf"  
vcf <- readVcf(vcf_file)

# contig length from Seqinfo
si <- seqinfo(vcf)
chr_df <- tibble(chrom = seqnames(si) |> as.character(), length = seqlengths(si) |> as.numeric()) %>%
  filter(!is.na(length))

rr <- rowRanges(vcf)

sv_df <- tibble(chrom = as.character(seqnames(rr)), pos = start(rr), svtype = info(vcf)$SVTYPE |> as.character()) %>%
  filter(svtype %in% c("INS", "DEL"))

names(info(vcf))

chr_levels <- c(paste0("chr", 1:30))

chr_df <- chr_df %>%
  mutate(chrom = factor(chrom, levels = chr_levels)) %>%
  arrange(chrom) %>%
  mutate(chrom2 = str_replace(chrom, "^chr", ""))
  
sv_df <- sv_df %>%
  mutate(chrom2 = str_replace(chrom, "^chr", ""))

chr_df2 <- chr_df %>%
  mutate(chrom = factor(chrom, levels = chr_levels)) %>%
  filter(!is.na(chrom)) %>%
  arrange(chrom) %>% 
  mutate(y = rev(seq_len(n())), ymin = y - 0.25, ymax = y + 0.25, xmin = 0, xmax = length)

bias_amt   <- 0.20 
jitter_amt <- 0.10
limit      <- 0.40

sv_df2 <- sv_df %>%
  mutate(chrom = factor(chrom, levels = chr_levels)) %>%
  left_join(chr_df2 %>% select(chrom, y), by = "chrom") %>%
  filter(!is.na(y))  %>%
  mutate(bias = case_when(svtype == "INS" ~  bias_amt, svtype == "DEL" ~ -bias_amt, TRUE ~  0), y_jit = y + bias + runif(n(), -jitter_amt, jitter_amt), y_jit = pmin(pmax(y_jit, y - limit), y + limit))

col_chr_fill   <- "#FFDD55"
col_chr_stroke <- "#806600"

col_del_fill   <- "#DE8787"
col_del_stroke <- "#800000"

col_ins_fill   <- "#5599FF"
col_ins_stroke <- "#000080"

p <- ggplot() +
  # Chromosome bars
  geom_segment(data = chr_df2, aes(x = 0, xend = length, y = y, yend = y), color = col_chr_stroke, linewidth = 6, lineend = "round") +
  geom_segment(data = chr_df2, aes(x = 0, xend = length, y = y, yend = y), color = col_chr_fill, linewidth = 5, lineend = "round") +
  # SV points
  geom_point(data = sv_df2, aes(x = pos, y = y_jit, fill = svtype, color = svtype), shape = 21, size = 0.5, stroke = 0.2, alpha = 0.5) +
  scale_fill_manual(values = c(DEL = col_del_fill, INS = col_ins_fill)) +
  scale_color_manual(values = c(DEL = col_del_stroke, INS = col_ins_stroke)) +
  scale_y_continuous(breaks = chr_df2$y, labels = ifelse(chr_df2$chrom == "chr29", "Z", ifelse(chr_df2$chrom == "chr30", "W", chr_df2$chrom)), name = "Chromosome") +
  scale_x_continuous(name = "Position (bp)", labels = scales::label_number(scale_cut = scales::cut_si(""), accuracy = 1)) +
  theme_classic(base_size = 12) +
  theme(legend.title = element_blank(), axis.ticks.y = element_blank()) +
  guides(color = "none", fill = guide_legend(override.aes = list(size = 3, stroke = 0.6)))
ggsave("all_sv_distribution_chromosome.svg", plot = p, width = 8, height = 8, units = "in", device = "svg")
ggsave("all_sv_distribution_chromosome.pdf", plot = p, width = 8, height = 8, units = "in", device = cairo_pdf)

# Read the 'outlier' data
counts <- read.table("outlier_sv_counts.tsv", header = FALSE, col.names = c("Chromosome", "SVTYPE"))
counts <- counts[counts$SVTYPE != "SNP",]
counts$Chromosome <- gsub("chr", "", counts$Chromosome)
genome <- read.table("/scratch/negishi/jeon96/monq/pan/monq-out_refp.genome", header = FALSE, col.names = c("Chromosome", "Length"))
genome$Chromosome <- gsub("chr", "", genome$Chromosome)
sv_summary <- counts %>% group_by(Chromosome, SVTYPE) %>% summarise(Count = n()) %>% ungroup()
sv_summary$Chromosome <- factor(sv_summary$Chromosome, levels = mixedsort(unique(sv_summary$Chromosome)))
sv_genome <- left_join(sv_summary, genome, by = "Chromosome") 
sv_genome$StandardizedCount <- sv_genome$Count / sv_genome$Length * 1000000 # number of SVs per Mb
sv_genome$Chromosome <- factor(sv_genome$Chromosome, levels = mixedsort(unique(sv_genome$Chromosome)))

p <- ggplot(sv_summary, aes(x = Chromosome, y = Count, fill = SVTYPE)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Distribution of SV Types per Chromosome", x = "Chromosome", y = "Count") +
  theme_minimal(base_size = 12) +
  theme(axis.ticks = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
  plot.title = element_text(hjust = 0.5))
ggsave("outlier_sv_distribution.svg", plot = p, width = 8, height = 8, units = "in", device = "svg")
ggsave("outlier_sv_distribution.pdf", plot = p, width = 8, height = 8, units = "in", device = cairo_pdf)

p <- ggplot(sv_genome, aes(x = Chromosome, y = StandardizedCount, fill = SVTYPE)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Distribution of SV Types per Chromosome", x = "Chromosome", y = "Length-corrected Count") +
  theme_minimal(base_size = 12) +
  theme(axis.ticks = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
  plot.title = element_text(hjust = 0.5))
ggsave("outlier_sv_standardized_distribution.svg", plot = p, width = 8, height = 8, units = "in", device = "svg")
ggsave("outlier_sv_standardized_distribution.pdf", plot = p, width = 8, height = 8, units = "in", device = cairo_pdf)

lengths <- read.table("outlier_sv_length.txt", header = FALSE, col.names = c("SVTYPE", "SVLEN"))
lengths <- lengths[lengths$SVTYPE != "SNP",]
lengths$SVLEN <- abs(lengths$SVLEN)
lengths$LengthBin <- cut(
  lengths$SVLEN,
  breaks = c(0, 50, 1000, 10000, 100000),
  labels = c("0-50", "50-1k", "1k-10k", "10k-100k"),
  right = FALSE
)
sv_length_summary <- lengths %>% group_by(SVTYPE, LengthBin) %>% summarise(Count = n(), .groups = "drop")

p <- ggplot(sv_length_summary, aes(x = LengthBin, y = Count, fill = SVTYPE)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = Count), position = position_dodge(width = 0.9), vjust = -0.3, size = 3) +
  labs(
    title = "SV Length Distribution by Type",
    x = "SV Length Bin (bp)",
    y = "Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggsave("outlier_sv_lengths.svg", plot = p, width = 8, height = 8, units = "in", device = "svg")
ggsave("outlier_sv_lengths.pdf", plot = p, width = 8, height = 8, units = "in", device = cairo_pdf)

# re-plot (raw distribution of SVs)
library(VariantAnnotation)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggridges)

vcf_file <- "monq_SVoutliers_final.vcf"  
vcf <- readVcf(vcf_file)

# contig length from Seqinfo
si <- seqinfo(vcf)
chr_df <- tibble(chrom = seqnames(si) |> as.character(), length = seqlengths(si) |> as.numeric()) %>%
  filter(!is.na(length))

rr <- rowRanges(vcf)

sv_df <- tibble(chrom = as.character(seqnames(rr)), pos = start(rr), svtype = info(vcf)$SVTYPE |> as.character()) %>%
  filter(svtype %in% c("INS", "DEL"))

names(info(vcf))

chr_levels <- c(paste0("chr", 1:30))

chr_df <- chr_df %>%
  mutate(chrom = factor(chrom, levels = chr_levels)) %>%
  arrange(chrom) %>%
  mutate(chrom2 = str_replace(chrom, "^chr", ""))
  
sv_df <- sv_df %>%
  mutate(chrom2 = str_replace(chrom, "^chr", ""))

chr_df2 <- chr_df %>%
  mutate(chrom = factor(chrom, levels = chr_levels)) %>%
  filter(!is.na(chrom)) %>%
  arrange(chrom) %>% 
  mutate(y = rev(seq_len(n())), ymin = y - 0.25, ymax = y + 0.25, xmin = 0, xmax = length)

bias_amt   <- 0.20 
jitter_amt <- 0.10
limit      <- 0.40

sv_df2 <- sv_df %>%
  mutate(chrom = factor(chrom, levels = chr_levels)) %>%
  left_join(chr_df2 %>% select(chrom, y), by = "chrom") %>%
  filter(!is.na(y))  %>%
  mutate(bias = case_when(svtype == "INS" ~  bias_amt, svtype == "DEL" ~ -bias_amt, TRUE ~  0), y_jit = y + bias + runif(n(), -jitter_amt, jitter_amt), y_jit = pmin(pmax(y_jit, y - limit), y + limit))

col_chr_fill   <- "#FFDD55"
col_chr_stroke <- "#806600"

col_del_fill   <- "#DE8787"
col_del_stroke <- "#800000"

col_ins_fill   <- "#5599FF"
col_ins_stroke <- "#000080"

p <- ggplot() +
  # Chromosome bars
  geom_segment(data = chr_df2, aes(x = 0, xend = length, y = y, yend = y), color = col_chr_stroke, linewidth = 6, lineend = "round") +
  geom_segment(data = chr_df2, aes(x = 0, xend = length, y = y, yend = y), color = col_chr_fill, linewidth = 5, lineend = "round") +
  # SV points
  geom_point(data = sv_df2, aes(x = pos, y = y_jit, fill = svtype, color = svtype), shape = 21, size = 1.5, stroke = 0.6, alpha = 0.75) +
  scale_fill_manual(values = c(DEL = col_del_fill, INS = col_ins_fill)) +
  scale_color_manual(values = c(DEL = col_del_stroke, INS = col_ins_stroke)) +
  scale_y_continuous(breaks = chr_df2$y, labels = ifelse(chr_df2$chrom == "chr29", "Z", ifelse(chr_df2$chrom == "chr30", "W", chr_df2$chrom)), name = "Chromosome") +
  scale_x_continuous(name = "Position (bp)", labels = scales::label_number(scale_cut = scales::cut_si(""), accuracy = 1)) +
  theme_classic(base_size = 12) +
  theme(legend.title = element_blank(), axis.ticks.y = element_blank()) +
  guides(color = "none", fill = guide_legend(override.aes = list(size = 3, stroke = 0.6)))
ggsave("outlier_sv_distribution_chromosome.svg", plot = p, width = 8, height = 8, units = "in", device = "svg")
ggsave("outlier_sv_distribution_chromosome.pdf", plot = p, width = 8, height = 8, units = "in", device = cairo_pdf)



# Gene ontology and KEGG pathway enrichment
module load biocontainers
module load bcftools
module load htslib
module load samtools
module load bedtools
module load bedops
module load anaconda
module load use.own
module load conda-env/mypackages-py3.9.13 # Negishi

## Normalizing outlier SVs
sed 's/ref_p#0#//g' ${MCPAN}/monq-out_refp.fa > ${MCPAN}/monq-out_refp_renamed.fa
samtools faidx ${MCPAN}/monq-out_refp_renamed.fa > ${MCPAN}/monq-out_refp_renamed.fa.fai

## Processing gene annotation files
cd ${MCPAN}
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/577/835/GCF_001577835.2_Coturnix_japonica_2.1/GCF_001577835.2_Coturnix_japonica_2.1_genomic.gtf.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gtf.gz 
gunzip GCF_001577835.2_Coturnix_japonica_2.1_genomic.gtf.gz
gunzip GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gtf.gz
cat GCF_001577835.2_Coturnix_japonica_2.1_genomic.gtf | grep "NC" | grep -E -v 'NC_029544.1|NC_029545.1|NC_003408.1' > jaqu_chrs.gtf # japanese quail chromosomes ~chr Z
cat GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gtf | grep "NC_052571.1" > chic_chW.gtf # chicken chromosome W
cat jaqu_chrs.gtf chic_chW.gtf > jaqu_plus_chicChW.gtf
awk 'BEGIN{FS=OFS="\t"}
     NR==FNR {                       # read mapping first
       if($1 !~ /^#/ && NF>=2) map[$7]=$1
       next
     }
     /^#/ {print; next}              # keep GTF comment lines
     {
       key=$1
       # (optional) strip version for matching, e.g., NC_000001.11 -> NC_000001
       # tmp=key; sub(/\.[0-9]+$/,"",tmp)
       # if(tmp in map) key=tmp
       if(key in map) $1=map[key]
       print
     }' chr_naming_pair.txt jaqu_plus_chicChW.gtf > jaqu_plus_chicChW_renamed.gtf # rename chromosomes to match with vcf files
#awk 'BEGIN{OFS="\t"}{print $1,0,$2}' monq-out_refp.fa.fai > genome.bed # whole genome as BED (0-based, half-open)
#sed 's/ref_p#0#//g' genome.bed > genome_renamed.bed
#awk '$3=="gene"' jaqu_plus_chicChW_renamed.gtf | gtf2bed | sort -k1,1 -k2,2n | bedtools merge -i - > functional.bed # define functional region (remove whole gene spans)
#bedtools subtract -a genome_renamed.bed -b functional.bed | sort -k1,1 -k2,2n | bedtools merge -i - > neutral.bed # make neutral region bed file
#bedtools intersect -a neutral.bed -b ok.bed | sort -k1,1 -k2,2n | bedtools merge -i - > neutral.ok.bed # generate 'new' analysis region

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/577/835/GCF_001577835.2_Coturnix_japonica_2.1/GCF_001577835.2_Coturnix_japonica_2.1_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gff.gz 
gunzip GCF_001577835.2_Coturnix_japonica_2.1_genomic.gff.gz
gunzip GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gff.gz
cat GCF_001577835.2_Coturnix_japonica_2.1_genomic.gff | grep "NC" | grep -E -v 'NC_029544.1|NC_029545.1|NC_003408.1' > jaqu_chrs.gff # japanese quail chromosomes ~chr Z
cat GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gff | grep "NC_052571.1" > chic_chW.gff # chicken chromosome W
cat jaqu_chrs.gff chic_chW.gff > jaqu_plus_chicChW.gff
awk 'BEGIN{FS=OFS="\t"}
     NR==FNR {                       # read mapping first
       if($1 !~ /^#/ && NF>=2) map[$7]=$1
       next
     }
     /^#/ {print; next}              # keep GFF comment lines
     {
       key=$1
       # (optional) strip version for matching, e.g., NC_000001.11 -> NC_000001
       # tmp=key; sub(/\.[0-9]+$/,"",tmp)
       # if(tmp in map) key=tmp
       if(key in map) $1=map[key]
       print
     }' chr_naming_pair.txt jaqu_plus_chicChW.gff > jaqu_plus_chicChW_renamed.gff # rename chromosomes to match with vcf files

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/577/835/GCF_001577835.2_Coturnix_japonica_2.1/GCF_001577835.2_Coturnix_japonica_2.1_gene_ontology.gaf.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_gene_ontology.gaf.gz 
gunzip GCF_001577835.2_Coturnix_japonica_2.1_gene_ontology.gaf.gz
gunzip GCF_016699485.2_gene_ontology.gaf.gz 

## Preparation for GO enrichment analysis
cd ${GEAGO}
bcftools norm -f ${MCPAN}/monq-out_refp_renamed.fa -Oz -o monq_SVoutliers.norm.vcf.gz monq_SVoutliers_final.vcf
tabix -p vcf monq_SVoutliers.norm.vcf.gz
bcftools query -f '%CHROM\t%POS0\t%END\t%SVTYPE\t%SVLEN\t%ID\n' monq_SVoutliers.norm.vcf.gz > monq_SVoutliers.bed # convert to bed
sort -k1,1 -k2,2n monq_SVoutliers.bed > monq_SVoutliers.sorted.bed

bcftools norm -f ${MCPAN}/monq-out_refp_renamed.fa -Oz -o monq_mcSV.final.norm.vcf.gz monq_mcSV.final.sort.vcf # process all SVs just like outliers
tabix -p vcf monq_mcSV.final.norm.vcf.gz 
bcftools query -f '%CHROM\t%POS0\t%END\t%SVTYPE\t%SVLEN\t%ID\n' monq_mcSV.final.norm.vcf.gz  > monq_SVall.bed # convert to bed
sort -k1,1 -k2,2n monq_SVall.bed > monq_SVall.sorted.bed

#awk -v OFS='\t' {'print $1,$2'} ${MCPAN}/monq-out_refp_renamed.fa.fai > ${MCPAN}/monq-out_refp.genome # make a genome file for bedtools

#bcftools view --no-header monq_SVoutliers.norm.vcf.gz | cut -f 8 | awk -F'SVLEN=' '/SVLEN=/ { split($2,a,";"); print a[1] }' | awk 'NR==1 || $1>max { max=$1 } END { print max }' # 7213
#bedtools slop -i monq_SVoutliers.sorted.bed -g ${MCPAN}/monq-out_refp.genome -b 7213 > monq_SVoutliers.sorted.padded.bed # add padding of the maximum SV size found
#awk -v OFS='\t' '
#  NR==FNR { gs[$1]=$2; next }                        # 1st: genome sizes
#  {
#    chrom=$1; start=$2; end=$3; svt=$4; svlen=$5; id=$6;
#    # for INS with empty end
#    if (end < start) end = start + 1;
#
#    pad=0;
#    if (svt=="INS") {
#      if (svlen=="." || svlen=="" ) svlen=0;
#      if (svlen < 0) svlen = -svlen;                 
#      pad = svlen;
#    }

#    newS = start - pad; if (newS < 0) newS = 0;
#    newE = end   + pad;
#    if (chrom in gs && newE > gs[chrom]) newE = gs[chrom];

#    print chrom, newS, newE, svt, svlen, id;
#  }
#' ${MCPAN}/monq-out_refp.genome monq_SVoutliers.sorted.bed > monq_SVoutliers.sorted.padded.bed # padding only INS

#bcftools view --no-header monq_mcSV.final.norm.vcf.gz | cut -f 8 | awk -F'SVLEN=' '/SVLEN=/ { split($2,a,";"); print a[1] }' | awk 'NR==1 || $1>max { max=$1 } END { print max }' # 14177
#bedtools slop -i monq_SVall.sorted.bed -g ${MCPAN}/monq-out_refp.genome -b 14177 > monq_SVall.sorted.padded.bed

bedtools sort -i monq_SVoutliers.sorted.bed | bedtools merge -i - -c 4,4 -o distinct,count > SVoutliers_sorted.merged.bed
bedtools sort -i monq_SVall.sorted.bed | bedtools merge -i - -c 4,4 -o distinct,count > SVall_sorted.merged.bed

#awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} $3=="gene"' ${MCPAN}/jaqu_plus_chicChW_renamed.gtf > functional.gtf
#bedtools sort -i functional.gtf > functional_sorted.gtf   # bedtools can sort GTF tabularly
#bedtools intersect -sorted -wa -wb -a SV_sorted.merged.bed -b functional_sorted.gtf > outlier_functional_hits.tsv 

GFF=${MCPAN}/jaqu_plus_chicChW_renamed.gff
GAF_QUAIL=${MCPAN}/GCF_001577835.2_Coturnix_japonica_2.1_gene_ontology.gaf
GAF_CHICKEN=${MCPAN}/GCF_016699485.2_gene_ontology.gaf

awk -F'\t' ' 
$0 ~ /^#/ { next }
$3=="gene" {
  attr=$9
  gid=""
  name=""
  # extract GeneID from Dbxref
  if (match(attr, /GeneID:([0-9]+)/, m)) gid=m[1]
  # extract Name= (gene symbol) if present
  if (match(attr, /(^|;)Name=([^;]+)/, n)) name=n[2]
  if (gid!="") {
    start=$4-1; if (start<0) start=0
    print $1, start, $5, gid, name, $7
  }
}' OFS='\t' "$GFF" \
| bedtools sort -i - > genes.bed # Build a genes BED from the merged GFF (chrom, start0, end, GeneID, gene_name, strand)

cat "$GAF_QUAIL" "$GAF_CHICKEN" \
  | awk -F'\t' '!/^!/' \
  | awk -F'\t' '
    {
      id=$2; go=$5
      # keep only proper GO terms
      if (go !~ /^GO:/) next
      # normalize common ID patterns; keep the numeric GeneID if present
      gsub(/^GeneID:/, "", id)
      # if your GAF uses Ensembl or UniProt IDs, see the optional mapping step below
      print id "\t" go
    }' \
  | sort -u > gene2go.raw.tsv # Merge the two GAFs and extract a 2-column gene→GO table
cp gene2go.raw.tsv gene2go.tsv

### collapse multiple GO terms per GeneID into a single semicolon-separated list
awk '{go[$1]=(go[$1]?go[$1]";":"")$2} END{for(k in go) print k"\t"go[k]}' \
  gene2go.tsv | sort -u > geneid_to_golist.tsv

### append GO list to the genes.bed (as column 7)
awk 'BEGIN{FS=OFS="\t"} NR==FNR{golist[$1]=$2; next}
     { print $0, ( ($4 in golist) ? golist[$4] : "." ) }' \
     geneid_to_golist.tsv genes.bed \
| bedtools sort -i - > genes_with_GO.bed

### append gene length
awk 'BEGIN{FS=OFS="\t"}
  NR==FNR {
    L = $3 - $2;                           # end - start
    if (!( $4 in len ) || L > len[$4])     # maximum for redundant gene IDs
      len[$4] = (L > 0 ? L : 1);           # safeguard
    next
  }
  {
    id = $1; go = $2;
    L  = (id in len ? len[id] : ".");
    print id, go, L
  }' genes_with_GO.bed geneid_to_golist.tsv > geneid_to_golist.with_len.tsv

bedtools intersect -wa -wb -a SVoutliers_sorted.merged.bed -b genes_with_GO.bed > SVoutlier_gene_hits.tsv 
bedtools intersect -wa -wb -a SVall_sorted.merged.bed -b genes_with_GO.bed > SVall_gene_hits.tsv

## Preparation for KEGG pathway enrichiment analysis
curl -s 'https://rest.kegg.jp/link/pathway/cjo'  > cjo_gene2path.tsv
curl -s 'https://rest.kegg.jp/list/pathway/cjo'  > cjo_pathnames.tsv
curl -s 'https://rest.kegg.jp/link/pathway/gga'  > gga_gene2path.tsv
curl -s 'https://rest.kegg.jp/list/pathway/gga'  > gga_pathnames.tsv

cat cjo_pathnames.tsv gga_pathnames.tsv > ko_pathnames.tsv

awk 'BEGIN{FS=OFS="\t"} $1=="chr30"{print $4}' genes.bed | sort -u > chickenW.ids
awk 'BEGIN{FS=OFS="\t"} $1!="chr30"{print $4}' genes.bed | sort -u > quail.ids

### Build quail gene ➜ KEGG pathway
awk 'NR==FNR{Q[$1]=1; next}
     {split($1,a,":"); if(a[1]=="cjo" && a[2] in Q) print a[2]"\t"$2}' \
    quail.ids cjo_gene2path.tsv \
  | sort -u > gene2path_quail.tsv

### Build chicken-W gene ➜ KEGG pathway
awk 'NR==FNR{W[$1]=1; next}
     {split($1,a,":"); if(a[1]=="gga" && a[2] in W) print a[2]"\t"$2}' \
    chickenW.ids gga_gene2path.tsv \
  | sort -u > gene2path_chickenW.tsv

cat gene2path_quail.tsv gene2path_chickenW.tsv | sort -u > gene2path_merged.tsv


# Gene ontology enrichment
R

#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#
#BiocManager::install("goseq")
#install.packages(c("data.table","dplyr"))

library(goseq)
library(data.table)
library(dplyr)
library(splitstackshape)

gene_info = read.table("geneid_to_golist.with_len.tsv", header=F, stringsAsFactors=F, sep="\t")
colnames(gene_info) <- c("GeneID", "GeneGo", "length")
head(gene_info)
row.names(gene_info) <- gene_info$GeneID

## split the GO terms
go_split = cSplit(gene_info, "GeneGo", sep = ";", type.convert=FALSE)
#go_split$contig <- gene_info$GeneID
head(go_split[,1:10])

## linearize the matrix
terms = colnames(select(go_split, contains("GeneGo")))
go_long = melt(go_split, measure.vars = terms, id.vars = "GeneID", na.rm = TRUE)
head(go_long)
go_ready = as.data.frame(go_long[, c(1, 3)])
head(go_ready)

## upload genes intersecting with all SVs
all_genes <- read.table("SVall_gene_hits.tsv", header = F)
colnames(all_genes)[9] <- "GeneID"
dim(all_genes) # how many?
head(all_genes)

## add size info
all_genes <- left_join(all_genes, gene_info[, c(1,3)])
dim(all_genes)
head(all_genes)

## make unique
all_genes_unique <- all_genes %>% distinct(GeneID, .keep_all = TRUE)
dim(all_genes_unique) # see the matrix has reduced...
# and we need to name the rows b gene names for goseq...
row.names(all_genes_unique) <- all_genes_unique$GeneID
head(all_genes_unique)

## genes in outliers
outliers_genes <- read.table("SVoutlier_gene_hits.tsv", header=F)
colnames(outliers_genes)[9] <- "GeneID"
head(outliers_genes)
dim(outliers_genes)

## add a column to the matrix with all the genes indicating whether this gene is found or not inside the outliers list. this will be a 0/1 vector
all_genes_unique$outliers_genes <- as.numeric(all_genes_unique$GeneID %in% outliers_genes$GeneID)
head(all_genes_unique)
all_genes_unique_valid <- subset(all_genes_unique, !is.na(length)) 

## prepare the data and integrate length bias for nullp. It requires vectors so here we go transforming data
measured_genes_v = as.vector(all_genes_unique_valid$GeneID)
outliers_genes_v = setNames(as.integer(all_genes_unique_valid$outliers_genes), all_genes_unique_valid$GeneID)
length_v = setNames(as.numeric(all_genes_unique_valid$length), all_genes_unique_valid$GeneID)
pwf_outliers = nullp(outliers_genes_v, bias.data = length_v)
row.names(pwf_outliers) <- row.names(all_genes_unique_valid)
head(pwf_outliers)

## test enrichment using our database
enrich_outliers = goseq(pwf_outliers, gene2cat = go_ready, use_genes_without_cat = TRUE)
head(enrich_outliers)

## correct for multiple testing with Benjamini & Hochberg correction
enrich_outliers$over_represented_padjust <- p.adjust(enrich_outliers$over_represented_pvalue, method = "BH")
head(enrich_outliers)
write.table(enrich_outliers, "GO_enrichment_SVoutliers.txt", sep="\t")

## see only the the significant enrichments
enrich_outliers[which(enrich_outliers$over_represented_padjust < 0.10), ] # nothing

subset(enrich_outliers, over_represented_pvalue < 0.01)
#       category over_represented_pvalue under_represented_pvalue numDEInCat
#334  GO:0003906            0.0009836821                1.0000000          2
#3578 GO:0071364            0.0043953933                1.0000000          2
#848  GO:0005929            0.0075412979                0.9990298          4
#     numInCat                                                      term
#334         2 DNA-(apurinic or apyrimidinic site) endonuclease activity
#3578        2     cellular response to epidermal growth factor stimulus
#848        20                                                    cilium
#     ontology over_represented_padjust
#334        MF                        1
#3578       BP                        1
#848        CC                        1

subset(enrich_outliers, over_represented_pvalue < 0.05)

hits <- subset(enrich_outliers, over_represented_pvalue < 0.01)$category
cat2genes <- split(go_ready$GeneID, go_ready$value)
go_genes <- names(outliers_genes_v)[outliers_genes_v == 1]
lapply(hits, function(go) intersect(cat2genes[[go]], go_genes))
#[[1]]
#[1] "107310660" "107313323" # APLF, NEIL3

#[[2]]
#[1] "107310110" "107318728" # GAREM1, PDE8A

#[[3]]
#[1] "107318197" "107313813" "107322881" "107314654" # TULP3, WDR19, IFT52, ABCC8


# KEGG pathway enrichment analysis
## universe = the exact gene IDs used in your GOseq vectors
universe_ids <- names(length_v)

## read merged mapping (two columns: gene_id, path_id)
g2p <- fread("gene2path_merged.tsv", col.names = c("gene_id","path"))

## keep only universe genes
g2p_universe <- g2p[gene_id %in% universe_ids]

## drop duplicates
kegg_ready <- unique(g2p_universe)

fwrite(kegg_ready, "gene2path_merged.universe.tsv", sep = "\t", col.names = FALSE)

enrich_kegg <- goseq(pwf_outliers, gene2cat = as.data.frame(kegg_ready),
                  use_genes_without_cat = TRUE)

enrich_kegg$over_represented_padjust <- p.adjust(enrich_kegg$over_represented_pvalue, "BH")

head(enrich_kegg)
write.table(enrich_kegg, "KEGG_enrichment_SVoutliers.txt", sep="\t")

## see only the the significant enrichments
enrich_kegg[which(enrich_kegg$over_represented_padjust < 0.10), ] # nothing

## Optional: pretty names for pathways (KO-based or gga-based both OK)
paths <- fread("ko_pathnames.tsv", col.names=c("path","term")) 
enrich_kegg$category <- sub("^path:", "", enrich_kegg$category)
enrich_kegg <- merge(enrich_kegg, paths, by.x="category", by.y="path", all.x = TRUE)
enrich_kegg[order(enrich_kegg$over_represented_padjust), ][1:20,
         c("category","term","numDEInCat","numInCat","over_represented_pvalue","over_represented_padjust")] # nothing significant
enrich_kegg[order(enrich_kegg$over_represented_pvalue), ][1:20,
         c("category","term","numDEInCat","numInCat","over_represented_pvalue","over_represented_padjust")] 

subset(enrich_kegg, over_represented_pvalue < 0.01) 
#   category over_represented_pvalue under_represented_pvalue numDEInCat
#16 cjo00230             0.006183038                0.9984079          8
#   numInCat over_represented_padjust
#16       43                        1
#                                                     term
#16 Purine metabolism - Coturnix japonica (Japanese quail)

subset(enrich_kegg, over_represented_pvalue < 0.05)
#    category over_represented_pvalue under_represented_pvalue numDEInCat
#16  cjo00230             0.006183038                0.9984079          8
#120 cjo04010             0.025263966                0.9891333         12
#    numInCat over_represented_padjust
#16        43                        1
#120      112                        1
#                                                           term
#16       Purine metabolism - Coturnix japonica (Japanese quail)
#120 MAPK signaling pathway - Coturnix japonica (Japanese quail)

hits <- subset(enrich_kegg, over_represented_pvalue < 0.05)$category
cat2genes <- split(kegg_ready$gene_id, sub("^path:", "", kegg_ready$path))
go_genes <- names(outliers_genes_v)[outliers_genes_v == 1]
lapply(hits, function(go) intersect(cat2genes[[go]], go_genes))
#[[1]]
#[1] "107311729" "107315810" "107316592" "107316682" "107316827" "107318728" # ENPP3, ADK, PDE1A, PDE3A, ADCY5, PDE8A
#[7] "107319927" "107319940" # NME7, FHIT

#[[2]]
# [1] "107306834" "107308786" "107308787" "107309296" "107309802" "107310899" # LOC107306834, MYD88, LOC107308787, TGFBR2, EGFR, SOS1
# [7] "107312792" "107313060" "107313520" "107315612" "107316660" "107320005" # LOC107312792, PDGFC, EREG, RET, ATF2, CACNA2D3


# SV PCA
R
library(vcfR)
library(tidyverse)
library(ggplot2)
library(pcadapt)
library(vegan)
library(ggrepel)

## All SVs
obj.vcfR <- read.vcfR("monq_mcSV.final.sort.vcf")
geno <- extract.gt(obj.vcfR) # extract and format genotype matrix

G <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno)) # an empty matrix, (9 stands for missing data)
G[geno %in% c("0/0", "0|0")] <- 0
G[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2
colnames(G) <- colnames(geno)
missing_ind <- colMeans(is.na(G))
missing_loci <- rowMeans(is.na(G))
G_clean <- G[missing_loci < 0.1, missing_ind < 0.2] # exclude "F1856", "F1913"

### PCAdapt -> missing data handled on the fly
G_pca <- read.pcadapt(G_clean, type = "pcadapt")
G_pca[is.na(G_pca)] <- 9
x <- pcadapt(input = G_pca, K = 20, pca.only = TRUE) 
pca_result <- x$u
eigenvalues <- x$d^2
prop_var <- eigenvalues / sum(eigenvalues) # PC1: 15.89%, PC2: 11.76%, PC3: 7.37%

meta <- as.matrix(read.table("/scratch/negishi/jeon96/monq/gea-go/monq_latlong.txt", header=TRUE))
meta <- as.data.frame(meta) %>% filter(!str_detect(sample_name, "F1913|F1856")) 

pca_df <- data.frame(pca_result)
pca_df <- cbind(meta, pca_df)
pca_df$latitude <- as.numeric(pca_df$latitude)
pca_df$longitude <- as.numeric(pca_df$longitude)
pca_df$lat_long_color <- (scale(pca_df$latitude) + scale(pca_df$longitude)) / 2
pca_df$label <- ifelse(pca_df$sample_name %in% c("E9056", "F1306", "F1789", "F1803", "F1847", "F1865", "F1867", "F1915", "F1927", "F1934", "F1952", "F1856", "F1889", "F1876"), pca_df$sample_name, NA)

svg("monq_pan_sv_pca1-2.svg", width=6, height=5)
ggplot(pca_df, aes(x = X1, y = X2, color = lat_long_color)) + 
  geom_point(size = 2.5, alpha = 0.5) + 
  geom_label_repel(aes(label = label), color = "black", max.overlaps = 14, box.padding = 0.75, alpha = 0.5) + 
  labs(x = "Principal Component 1 (15.89%)", y = "Principal Component 2 (11.76%)") + 
  theme_classic(base_size=11) + scale_color_gradientn(colors = c("red", "blue")) +
  theme(legend.position = "none")
dev.off()

svg("monq_pan_sv_pca2-3.svg", width=6, height=5)
ggplot(pca_df, aes(x = X2, y = X3, color = lat_long_color)) + 
  geom_point(size = 2.5, alpha = 0.5) + 
  geom_label_repel(aes(label = label), color = "black", max.overlaps = 14, box.padding = 0.75, alpha = 0.5) + 
  labs(x = "Principal Component 2 (11.76%)", y = "Principal Component 3 (7.37%)") + 
  theme_classic(base_size=11) + scale_color_gradientn(colors = c("red", "blue")) +
  theme(legend.position = "none")
dev.off()

## outlier SVs
obj.vcfR <- read.vcfR("monq_SVoutliers_final.vcf")
geno <- extract.gt(obj.vcfR) # extract and format genotype matrix

G <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno)) # an empty matrix, (9 stands for missing data)
G[geno %in% c("0/0", "0|0")] <- 0
G[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2
colnames(G) <- colnames(geno)
missing_ind <- colMeans(is.na(G))
missing_loci <- rowMeans(is.na(G))
G_clean <- G[missing_loci < 0.1, missing_ind < 0.2] # exclude "F1856", "F1913"

### PCAdapt -> missing data handled on the fly
G_pca <- read.pcadapt(G_clean, type = "pcadapt")
G_pca[is.na(G_pca)] <- 9
x <- pcadapt(input = G_pca, K = 20, pca.only = TRUE) 
pca_result <- x$u
eigenvalues <- x$d^2
prop_var <- eigenvalues / sum(eigenvalues) # PC1: 29.37%, PC2: 13.52%, PC3: 7.37%

meta <- as.matrix(read.table("/scratch/negishi/jeon96/monq/gea-go/monq_latlong.txt", header=TRUE))
meta <- as.data.frame(meta) %>% filter(!str_detect(sample_name, "F1913|F1856")) 

pca_df <- data.frame(pca_result)
pca_df <- cbind(meta, pca_df)
pca_df$latitude <- as.numeric(pca_df$latitude)
pca_df$longitude <- as.numeric(pca_df$longitude)
pca_df$lat_long_color <- (scale(pca_df$latitude) + scale(pca_df$longitude)) / 2
pca_df$label <- ifelse(pca_df$sample_name %in% c("E9056", "F1306", "F1789", "F1803", "F1847", "F1865", "F1867", "F1915", "F1927", "F1934", "F1952", "F1856", "F1889", "F1876"), pca_df$sample_name, NA)

svg("monq_pan_svoutlier_pca1-2.svg", width=6, height=5)
ggplot(pca_df, aes(x = X1, y = X2, color = lat_long_color)) + 
  geom_point(size = 2.5, alpha = 0.5) + 
  geom_label_repel(aes(label = label), color = "black", max.overlaps = 14, box.padding = 0.75, alpha = 0.5) + 
  labs(x = "Principal Component 1 (29.37%)", y = "Principal Component 2 (13.52%)") + 
  theme_classic(base_size=11) + scale_color_gradientn(colors = c("red", "blue")) +
  theme(legend.position = "none")
dev.off()

svg("monq_pan_svoutlier_pca2-3.svg", width=6, height=5)
ggplot(pca_df, aes(x = X2, y = X3, color = lat_long_color)) + 
  geom_point(size = 2.5, alpha = 0.5) + 
  geom_label_repel(aes(label = label), color = "black", max.overlaps = 14, box.padding = 0.75, alpha = 0.5) + 
  labs(x = "Principal Component 2 (13.52%)", y = "Principal Component 3 (7.37%)") + 
  theme_classic(base_size=11) + scale_color_gradientn(colors = c("red", "blue")) +
  theme(legend.position = "none")
dev.off()


# Site frequency spectrum
vcftools --vcf monq_mcSV.final.sort.vcf --freq2 --out SVall_freq
vcftools --vcf monq_SVoutliers_final.vcf --freq2 --out SVoutlier_freq

R
library(data.table)
library(tidyr)
library(dplyr)
library(stringr)

freq <- fread("SVall_freq.frq", sep = "\t", fill = TRUE, quote = "", header = TRUE, data.table = FALSE, check.names = FALSE)
pair_start <- 5L # firt column that pair starts

freq_long <- data.table(
  CHROM = freq$CHROM,
  POS   = freq$POS,
  N_CHR = freq$N_CHR,
  pairs = apply(freq[, pair_start:ncol(freq)], 1, function(z) paste(z[!is.na(z)], collapse=" "))
) |>
  separate_rows(pairs, sep="\\s+") |>
  mutate(pairs = str_replace_all(pairs, "[{}]", "")) |>
  separate(pairs, c("ALLELE","FREQ"), ":", convert=TRUE)

by_site <- freq_long %>%
  group_by(CHROM, POS) %>%
  summarise(
    N_CHR = first(N_CHR),
    f_max = max(FREQ, na.rm=TRUE),
    .groups="drop"
  ) |>
  mutate(
    mac = round(pmax(0, N_CHR * (1 - f_max)))   # fold to minor count
  )

by_site <- by_site |>
  filter(mac >= 1, mac <= floor(N_CHR/2))

norm <- function(x) x/sum(x) # function to normalize 
sfs <- tabulate(by_site$mac)
sfs <- norm(sfs)  # the variable categories of the sfs 
svg("monq_SVall_sfs_plot.svg", width=5, height=5)
barplot(sfs, xlab = "Allele frequency", names = 1:length(sfs), ylab = "Proportions", main = "Site Frequency Spectrum plot", col = 'blue')
dev.off()

freq <- fread("SVoutlier_freq.frq", sep = "\t", fill = TRUE, quote = "", header = TRUE, data.table = FALSE, check.names = FALSE)
pair_start <- 5L # firt column that pair starts

freq_long <- data.table(
  CHROM = freq$CHROM,
  POS   = freq$POS,
  N_CHR = freq$N_CHR,
  pairs = apply(freq[, pair_start:ncol(freq)], 1, function(z) paste(z[!is.na(z)], collapse=" "))
) |>
  separate_rows(pairs, sep="\\s+") |>
  mutate(pairs = str_replace_all(pairs, "[{}]", "")) |>
  separate(pairs, c("ALLELE","FREQ"), ":", convert=TRUE)

by_site <- freq_long %>%
  group_by(CHROM, POS) %>%
  summarise(
    N_CHR = first(N_CHR),
    f_max = max(FREQ, na.rm=TRUE),
    .groups="drop"
  ) |>
  mutate(
    mac = round(pmax(0, N_CHR * (1 - f_max)))   # fold to minor count
  )

by_site <- by_site |>
  filter(mac >= 1, mac <= floor(N_CHR/2))

norm <- function(x) x/sum(x) # function to normalize 
sfs <- tabulate(by_site$mac)
sfs <- norm(sfs)  # the variable categories of the sfs 
svg("monq_SVoutlier_sfs_plot.svg", width=5, height=5)
barplot(sfs, xlab = "Allele frequency", names = 1:length(sfs), ylab = "Proportions", main = "Site Frequency Spectrum plot", col = 'blue')
dev.off()



# New Mexico samples from NCBI
cd ${MCPAN}/called/sv/vg/new_mexico

## Compressing and indexing each vcf file first
for i in `ls -lrt /scratch/bell/jeon96/monq/short-read/sra/fixed/ | tr -s ' ' | cut -d ' ' -f 9 | cut -d '_' -f 1 | sort | uniq`; do
  #id=$(echo $i | cut -d '_' -f 1)
  gunzip ${i}_concat_mcSV.vcf.gz 
  sed 's/ref_p#0#//g' ${i}_concat_mcSV.vcf > ${i}_mcSV_temp.vcf # polish the header line to be compatible with bcftools
  bcftools sort ${i}_mcSV_temp.vcf -Oz -o ${i}_mcSV.sorted.vcf.gz
  bcftools index ${i}_mcSV.sorted.vcf.gz --threads ${N}
  bcftools view ${i}_mcSV.sorted.vcf.gz --threads ${N} | grep -v "##" | wc -l # number of variants
  rm ${i}_mcSV_temp.vcf
done

## Combining separately called SV vcf files
#bcftools merge -m all -Oz -o monqNM_mcSV.merged.vcf.gz --threads ${N} *_mcSV.sorted.vcf.gz 
bcftools merge -m all -Oz -o monq_mcSV_wNM.merged.vcf.gz --threads ${N} *_mcSV.sorted.vcf.gz 

## Filtering vcf file based on the determined conservative cut-off values from above results                                                                      
vcftools --gzvcf monq_mcSV_wNM.merged.vcf.gz --out monq_mcSV_wNM.filtered --minQ 20 --mac 2 --recode --recode-INFO-all #--remove-filtered lowad --max-meanDP 15 
bcftools sort monq_mcSV_wNM.filtered.recode.vcf -Oz -o monq_mcSV_wNM.filtered.sorted.vcf.gz
${APP}/vcfbub --input monq_mcSV_wNM.filtered.sorted.vcf.gz --max-allele-length 100000 --max-level 0 > monq_mcSV_wNM.filt.sort.popped.vcf.gz # remove large (>100 kb) alleles in graphs. This removes all non-top-level bubbles from the VCF unless they were nested inside a top-level bubble with a reference length exceeding 100 kb; that is, top-level bubbles longer than that are replaced by their child nodes in the snarl tree.

bgzip -c monq_mcSV_wNM.filt.sort.popped.vcf.gz > monq_mcSV_wNM.filt.sort.popped.vcf.bgz
tabix monq_mcSV_wNM.filt.sort.popped.vcf.bgz
singularity run ${APP}/vcflib-1.09.simg
N=64
vcfwave -I 1000 -t ${N} monq_mcSV_wNM.filt.sort.popped.vcf.bgz > monq_mcSV_wNM.filt.popped.decomp.vcf.gz # decompose complex alleles into primitive ones
exit

bcftools query -f '%CHROM\t%POS\n' ${CALLED}/monq_mcSV.final.sort.vcf > vcf_regions.txt # querying the sites of the main vcf file
bcftools sort monq_mcSV_wNM.filt.popped.decomp.vcf.gz -Oz -o monq_mcSV_wNM.filt.popped.decomp.sort.vcf.gz
bcftools index monq_mcSV_wNM.filt.popped.decomp.sort.vcf.gz
bcftools view -R vcf_regions.txt monq_mcSV_wNM.filt.popped.decomp.sort.vcf.gz -Oz -o monq_mcSV_wNM.subset.vcf.gz # make compatible with the main vcf 

bcftools query -f '%CHROM\t%POS\n' ${CALLED}/monq_SVoutliers_final.vcf > vcf_regions.txt # querying the sites of the main vcf file
bcftools view -R vcf_regions.txt monq_mcSV_wNM.subset.vcf.gz -Oz -o monq_SVoutliers_wNM.subset.vcf.gz
