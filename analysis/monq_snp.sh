#!/bin/bash
#SBATCH --job-name=monq_snp
#SBATCH -A fnrpupfish
#SBATCH -p cpu
#SBATCH -t 14-00:00:00
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=jeon96@purdue.edu
#SBATCH --mail-type=END,FAIL

# Population genetics of SNPs

## Preprocessing the reference genome and bam files
module purge
module load biocontainers
module load angsd
module load picard
module load samtools
module load bcftools
module load boost
module load genmap
module load repeatmasker
module load bedtools
module load bbmap
module load bedops

cd ${MCPAN}
${APP}/vg index -t ${N} -p monq-out_mod.vg -x monq-out.mod.xg
${APP}/vg paths -t ${N} -x monq-out.mod.xg -L | head -n 50 # check path names
${APP}/vg paths -t ${N} -x monq-out.mod.xg -F -S ref_p > monq-out_refp.fa
samtools faidx monq-out_refp.fa # reference assembly which was used in pangenome building
cp monq-out_refp.fa original.fa
reformat.sh in=original.fa out=new.fa trd=t -Xmx20g overwrite=T # reduce fasta header length
sortbyname.sh in=new.fa out=ref.fa -Xmx20g length descending overwrite=T # sort by lengthv
reformat.sh minlength=100000 in=ref.fa out=ref_100kb.fa overwrite=T # remove sequences smaller that 100kb prior to any repeatmasking
samtools faidx ref_100kb.fa

RepeatMasker -pa 64 -a -qq -species aves -dir . ref_100kb.fa # repeat mask
samtools faidx ref_100kb.fa
cat ref_100kb.fa.out  | tail -n +4 | awk '{print $5,$6,$7,$11}' | sed 's/ /\t/g' > repeats.bed 

genmap index -F ref_100kb.fa -I index -S 50 # assess mappability of reference
mkdir mappability
genmap map -K 100 -E 2 -T 64 -I index -O mappability -t -w -bg    

sortBed -i repeats.bed > repeats_sorted.bed # sort bed 
awk 'BEGIN {FS="\t"}; {print $1 FS $2}' ref_100kb.fa.fai > ref.genome
awk '{print $1, $2, $2}' ref.genome > ref2.genome
sed -i 's/ /\t/g' ref2.genome
sortBed -i ref2.genome > ref3.genome
awk '{print $1, $2 }' ref3.genome > ref_sorted.genome
sed -i 's/ /\t/g' ref_sorted.genome
bedtools complement -i repeats_sorted.bed -g ref_sorted.genome > nonrepeat.bed # find nonrepeat regions
awk '$4 == 1' mappability/ref_100kb.genmap.bedgraph > mappability/map.bed # clean mappability file to remove sites with <1 mappability                                                  
awk 'BEGIN {FS="\t"}; {print $1 FS $2 FS $3}' mappability/map.bed > mappability/mappability.bed
sortBed -i mappability/mappability.bed > mappability/mappability_sorted.bed
sed -i 's/ /\t/g' mappability/mappability_sorted.bed
bedtools subtract -a mappability/mappability_sorted.bed -b repeats_sorted.bed > mappability/map_nonreapeat.bed
bedtools sort -i mappability/map_nonreapeat.bed > mappability/filter_sorted.bed
bedtools merge -i mappability/filter_sorted.bed > mappability/merged.bed

awk '{ print $1, $2, $2 }' ref_100kb.fa.fai > ref_100kb.info # make bed file with the 100k and merged.bed (no repeats, mappability =1 sites)
awk '$2="0"' ref_100kb.info > ref_100kb.bed
sed -i 's/ /\t/g' ref_100kb.bed
bedtools intersect -a ref_100kb.bed -b ./mappability/merged.bed > ok.bed # only include scaffolds in merged.bed if they are in ref_100kb.bed	 
cut -f 1 ref_100kb.bed > chrs.txt # make chrs.txt

rm new.fa
rm ref.genome
rm ref2.genome
rm ref3.genome

# ANGSD for PCA
module purge
module load biocontainers
module load angsd
module load ngsld
module load r
module load pcangsd
module load bcftools

mkdir -p ./popgen
cd ./popgen/

ls /scratch/bell/jeon96/monq/pan/aligned/*merged.sorted.bam | sort > ./bam.filelist

awk '{print $1"\t"$2+1"\t"$3}' ../ok.bed > ./angsd.file

MIND=$((`wc -l < ./bam.filelist` / 2))
angsd sites index ./angsd.file
angsd -bam ./bam.filelist -anc ../ref_100kb.fa -ref ../ref_100kb.fa -rf ../chrs.txt -sites ./angsd.file -GL 1 -doglf 2 -minmaf 0.05 -snp_pval 1e-6 -dopost 1 -minQ 30 -minInd ${MIND} -docounts 1 -domajorminor 5 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out pan_snp_pca # generate SNP beagle

## LD pruning
zcat pan_snp_pca.beagle.gz | tail -n +2 | awk 'NR % 10 == 0' | cut -f 4- | gzip  > pan_snp_pca_subsampled.geno.beagle.gz # prepare a geno file by subsampling one SNP in every 10 SNPs in the beagle file
zcat pan_snp_pca.mafs.gz | tail -n +2 | cut -f 1,2 | awk 'NR % 10 == 0' | sed 's/:/_/g'| gzip > pan_snp_pca_subsampled.pos.gz # prepare a pos file by subsampling one SNP in every 10 SNPs in the beagle file 
NSITES=$(zcat pan_snp_pca_subsampled.pos.gz | wc -l)

ngsLD --geno pan_snp_pca_subsampled.geno.beagle.gz --pos pan_snp_pca_subsampled.pos.gz --probs --n_ind 91 --n_sites ${NSITES} --max_kb_dist 1000 --n_threads ${N} --out out_pan_snp.ngsld 

START_SITE=$(cat out_pan_snp.ngsld | grep "chr1" | cut -f1 | cut -f2 -d ":" | sort -n | head -n 1)
END_SITE=$(cat out_pan_snp.ngsld | grep "chr1" | cut -f2 | cut -f2 -d ":" | sort -n | tail -n 1)
#bash "$HOME/LD_blocks_modified.sh" "ref_p#0#chr1" "$START_SITE" "$END_SITE" < out_pan_snp.ngsld
cat out_pan_snp.ngsld | grep "chr1" > out_pan_snp_1stchr.ngsld
#cat out_pan_snp.ngsld | grep "chr1" > out_pan_snp_1mbchr.ngsld # this is the chromosome with its length of 1Mb
bash "$HOME/LD_blocks_modified.sh" "ref_p#0#chr1" "$START_SITE" "$END_SITE" < out_pan_snp_1stchr.ngsld
bash "${HOME}/LD_blocks_modified.sh" "ref_p#0#chr1" 5000 1005000 < out_pan_snp_1stchr.ngsld
echo -e "pos1\tpos2\tdist\tr2_ExpG\tD\tDp\tr2" > decay.header
cat decay.header out_pan_snp_1stchr.ngsld > out_pan_snp_1stchr_header.ngsld
echo "${PWD}/out_pan_snp_1stchr_header.ngsld" > decay_snp.ld
Rscript --vanilla --slave ${HOME}/fit_LDdecay.R --ld_files decay_snp.ld --max_kb_dist 1000 --out fit_SNP_LDdecay.pdf --n_ind 91 --fit_level 100 --plot_x_lim 100 --plot_size 5,6 --plot_axis_scale free

${APP}/prune_graph/target/release/prune_graph --in out_pan_snp.ngsld --weight-field "column_7" --weight-filter "column_3 <= 100000 && column_7 >= 0.5" --out out_pan_snp_unlinked.ngsld
sed 's/:/_/g' out_pan_snp_unlinked.ngsld > out_pan_snp_unlinked

zcat pan_snp_pca.beagle.gz | head -n 1 > pan_snp_pca_ld.beagleheader
zcat pan_snp_pca.beagle.gz | grep -Fwf out_pan_snp_unlinked > pan_snp_pca_ld.beagle # keep only unlinked loci
cat pan_snp_pca_ld.beagleheader pan_snp_pca_ld.beagle | gzip > pan_snp_pca_ld.beagle.gz
zcat pan_snp_pca_ld.beagle.gz | head
cat out_pan_snp_unlinked | wc -l # check the number of loci if matching
zcat pan_snp_pca_ld.beagle.gz | tail -n +2 | wc -l # check the number of loci if matching

pcangsd -b pan_snp_pca_ld.beagle.gz -o monq_pan_snp_pca -t ${N}

ml sqlite
module load geos
module load proj
module load gdal

R
library(ggplot2)
library(ggrepel)

# PCAngsd
meta <- as.matrix(read.table("/scratch/bell/jeon96/monq/pan/popgen/monq_latlong.txt", header=TRUE))
C <- as.matrix(read.table("/scratch/bell/jeon96/monq/pan/popgen/monq_pan_snp_pca.cov"))
e <- eigen(C)
pca_result <- prcomp(C, scale = TRUE)
summary(pca_result) # PC1: 21.43%, PC2: 13.18%, PC3: 1.75%
pca_df <- data.frame(pca_result$x)
pca_df <- cbind(meta, pca_df)
pca_df$latitude <- as.numeric(pca_df$latitude)
pca_df$longitude <- as.numeric(pca_df$longitude)
pca_df$lat_long_color <- (scale(pca_df$latitude) + scale(pca_df$longitude)) / 2
pca_df$label <- ifelse(pca_df$sample_name %in% c("E9056", "F1306", "F1789", "F1803", "F1847", "F1865", "F1867", "F1915", "F1927", "F1934", "F1952", "F1856", "F1889", "F1876"), pca_df$sample_name, NA)

svg("monq_pan_snp_pca1-2.svg", width=6, height=5)
ggplot(pca_df, aes(x = PC1, y = PC2, color = lat_long_color)) + 
  geom_point(size = 2.5, alpha = 0.5) + 
  geom_label_repel(aes(label = label), color = "black", max.overlaps = 10, box.padding = 0.75, alpha = 0.5) + 
  labs(x = "Principal Component 1 (21.43%)", y = "Principal Component 2 (13.18%)") + 
  theme_classic(base_size=11) + scale_color_gradientn(colors = c("red", "blue")) +
  theme(legend.position = "none")
dev.off()

svg("monq_snp_pca2-3.svg", width=6, height=5)
ggplot(pca_df, aes(x = PC2, y = PC3, color = lat_long_color)) + 
  geom_point(size = 2.5, alpha = 0.5) + 
  geom_label_repel(aes(label = label), color = "black", max.overlaps = 10, box.padding = 0.75, alpha = 0.5) + 
  labs(x = "Principal Component 2 (13.18%)", y = "Principal Component 3 (1.75%)") + 
  theme_classic(base_size=11) + scale_color_gradientn(colors = c("red", "blue")) +
  theme(legend.position = "none")
dev.off()

library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(maptiles)
library(ggspatial)

points_sf <- st_as_sf(pca_df, coords = c("longitude", "latitude"), crs = 4326)
points_sf_3857 <- st_transform(points_sf, 3857)
#world <- ne_countries(scale = "medium", returnclass = "sf")
#bbox <- sf::st_bbox(points_sf)

label_df <- pca_df[!is.na(pca_df$label), ]
label_sf <- st_as_sf(label_df, coords = c("longitude", "latitude"), crs = 4326)
label_sf_3857 <- st_transform(label_sf, 3857)
label_coords <- st_coordinates(label_sf_3857)
label_df$X <- label_coords[, 1]
label_df$Y <- label_coords[, 2]

x_range <- range(pca_df$longitude, na.rm = TRUE)
y_range <- range(pca_df$latitude, na.rm = TRUE)
x_buffer <- diff(x_range) * 0.05
y_buffer <- diff(y_range) * 0.05
xlim_expanded <- c(x_range[1] - x_buffer, x_range[2] + x_buffer)
ylim_expanded <- c(y_range[1] - y_buffer, y_range[2] + y_buffer)

bbox_wgs84 <- st_as_sfc(st_bbox(c(xmin = xlim_expanded[1], xmax = xlim_expanded[2], ymin = ylim_expanded[1], ymax = ylim_expanded[2]), crs = 4326))
bbox_3857 <- st_transform(bbox_wgs84, 3857)
bbox_coords <- st_bbox(bbox_3857)
xlim_3857 <- c(bbox_coords["xmin"], bbox_coords["xmax"])
ylim_3857 <- c(bbox_coords["ymin"], bbox_coords["ymax"])

#basemap <- get_tiles(bbox_3857, provider = "Stadia.StamenTerrainBackground", apikey="391f7013-698d-4c56-9b6f-f725ebd8d19a", crop = TRUE, zoom = 8)
basemap <- get_tiles(bbox_3857, provider = "Thunderforest.Landscape", crop = TRUE, zoom = 8)
basemap <- terra::flip(basemap, direction = "vertical")
world <- ne_countries(scale = "medium", returnclass = "sf") |> st_transform(3857) |> st_crop(bbox_3857) 
states <- ne_states(returnclass = "sf") |> st_transform(3857) |> st_crop(bbox_3857)

png("monq_sample_pca.png", width=6, height=5, units = "in", res = 300)
ggplot() +
  #layer_spatial(basemap) + 
  geom_sf(data = world, fill = NA, color = "black", linewidth = 0.4) +   
  geom_sf(data = states, fill = NA, color = "gray40", linewidth = 0.2) + 
  #geom_sf(data = world, fill = "gray95", color = "gray70") +
  geom_sf(data = points_sf_3857, aes(color = lat_long_color), size = 2.5, alpha = 0.5) +
  geom_label_repel(data = label_df, aes(x = X, y = Y, label = label), force = 2.5, box.padding = 0.5, point.padding = 0.5, max.overlaps = Inf, min.segment.length = 0, alpha = 0.75) +
  #geom_label_repel(data = label_df, aes(x = X, y = Y, label = label), max.overlaps = 10, box.padding = 0.75, alpha = 0.5) +
  scale_color_gradientn(colors = c("red", "blue")) +
  theme_classic(base_size=11) +
  #coord_sf(xlim = xlim_expanded, ylim = ylim_expanded, expand = FALSE) +
  #coord_sf(xlim = xlim_3857, ylim = ylim_3857, expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme(legend.position = "none")
dev.off()

#svg("monq_sample_pca.svg", width=6, height=5)
p <- ggplot() +
  #layer_spatial(basemap) + 
  geom_sf(data = world, fill = NA, color = "black", linewidth = 0.6) +   
  geom_sf(data = states, fill = NA, color = "gray40", linewidth = 0.3) + 
  #geom_sf(data = world, fill = "gray95", color = "gray70") +
  geom_sf(data = points_sf_3857, aes(color = lat_long_color), size = 2.5, alpha = 0.5) +
  geom_label_repel(data = label_df, aes(x = X, y = Y, label = label), force = 2.5, box.padding = 0.5, point.padding = 0.5, max.overlaps = Inf, min.segment.length = 0, alpha = 0.75) +
  #geom_label_repel(data = label_df, aes(x = X, y = Y, label = label), max.overlaps = 10, box.padding = 0.75, alpha = 0.5) +
  scale_color_gradientn(colors = c("red", "blue")) +
  theme_classic(base_size=11) +
  #coord_sf(xlim = xlim_expanded, ylim = ylim_expanded, expand = FALSE) +
  #coord_sf(xlim = xlim_3857, ylim = ylim_3857, expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme(legend.position = "none")
#dev.off()
ggsave("monq_sample_pca.svg", p, width = 6, height = 5, bg = "transparent")

quit()

# Effective population size to use in SLiM
## Subset individuals of each population from a vcf
module load biocontainers
module load bcftools

bcftools view -s SampleA,SampleC -Ov -o monqCT_mcSNP.vcf /scratch/bell/jeon96/monq/pan/called/snp/monq_mcSNP.final.sort.vcf