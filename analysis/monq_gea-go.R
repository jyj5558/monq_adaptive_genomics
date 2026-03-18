module load geos
module load gdal
module load proj
module load r
module load sqlite

cd /scratch/negishi/jeon96/monq/gea-go/

R
# Install packages
#packages = c("pegas", "terra", "raster", "ggplot2", "rgda", "rnaturalearth", "rnaturalearthdata",
#             "RColorBrewer", "ggpubr", "vegan", "robust", "ggVennDiagram",
#             "cowplot","corrplot","rgeos","sf","dplyr","corrplot", "ggnewscale")
#if (!require(packages, quietly = TRUE))
#  install.packages(packages)
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("LEA")
#BiocManager::install("qvalue")
#BiocManager::install("WMDB")
#install.packages(c("MatrixModels", "quantreg", "car", "rstatix", "ggpubr"), lib = "/home/jeon96/R/negishi/4.2.2-gcc-12.2.0-524zdta")

# Change below path to your directory
#install.packages("/home/jeon96/app/WMDB_1.0.tar.gz",
#                 repos=NULL, type='source')
#if (!("bigstatsr" %in% installed.packages())){install.packages("bigstatsr")}
#if (!("bigsnpr" %in% installed.packages())){devtools::install_github("privefl/bigsnpr")}

# Load libraries
library(vcfR)
library(tidyverse)
library(pegas)
library(ggplot2)
library(raster)
#library(rgdal)
library(LEA)
library(rnaturalearth)
library(rnaturalearthdata)
library(RColorBrewer)
library(Matrix, lib.loc = "/home/jeon96/R/negishi/4.2.2-gcc-12.2.0-524zdta")
library(MatrixModels, lib.loc="/home/jeon96/R/negishi/4.2.2-gcc-12.2.0-524zdta")
library(ggpubr)
library(vegan)
library(qvalue)
library(robust)
library(WMDB)
library(ggVennDiagram)
library(cowplot)
library(corrplot)
#library(rgeos)
library(terra)
library(sf)
library(dplyr)
library(ggnewscale)
library(missMDA)
library(svglite) 
library(foreach)
library(gtools)
library(OutFLANK)
library(bigstatsr)
library(bigsnpr)   # package for SNP trimming
library(pcadapt)

# Fst outlier test - no avail
obj.vcfR <- read.vcfR("monq_mcSV.final.sort.vcf")
chr_pos <- data.frame(
  id = getID(obj.vcfR),
  chr = getCHROM(obj.vcfR),
  pos = as.numeric(getPOS(obj.vcfR)),
  stringsAsFactors = FALSE
)
rownames(chr_pos) <- chr_pos$id

#position <- getPOS(obj.vcfR) # Positions in bp
#chromosome <- getCHROM(obj.vcfR) # Chromosome information
#id_sv <- getID(obj.vcfR) # ID of the SNP
#chr_pos <- as.data.frame(cbind(id_sv, chromosome, position)) 
#chr_pos$position <- as.numeric(as.character(chr_pos$position))
write.table(chr_pos, "monq_sv_pos.txt", sep="\t", quote=F, row.names=F)
geno <- extract.gt(obj.vcfR) # extract and format genotype matrix

G <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno)) # an empty matrix, (9 stands for missing data)
G[geno %in% c("0/0", "0|0")] <- 0
G[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2
colnames(G) <- colnames(geno)

mean(is.na(G)) # 0.02195726 # proportion of missing values in the whole matrix

missing_ind <- colMeans(is.na(G)) # proportion of missing per individual
summary(missing_ind)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.003448 0.004956 0.006000 0.021957 0.007824 0.999834
hist(missing_ind, main = "Missing rate per individual", xlab = "Proportion missing")

missing_loci <- rowMeans(is.na(G)) # proportion of missing per locus
summary(missing_loci)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.00000 0.01099 0.01099 0.02196 0.02198 0.97802
hist(missing_loci, main = "Missing rate per locus", xlab = "Proportion missing")

G_clean <- G[missing_loci < 0.1, missing_ind < 0.2]  # keep individuals with <20% missing (excluding "F1856", "F1913")
sum(is.na(G_clean))               # total number of NAs
rowSums(is.na(G_clean))           # missing count per locus
colSums(is.na(G_clean))           # missing count per individual

chr_pos_clean <- chr_pos[missing_loci < 0.1, ]

chr_pos_filtered <- chr_pos_clean[!zero_var, ]
infos.chr <- as.integer(gsub("chr", "", chr_pos_filtered$chr))
infos.pos <- chr_pos_filtered$pos

keep_ind  <- which(missing_ind < 0.2)
kept_ids  <- colnames(G)[keep_ind]  

X <- t(G_clean)
rownames(X) <- kept_ids

col_uniq <- apply(X, 2, function(v) length(unique(v[!is.na(v)])))
table(col_uniq)
col_uniq

zero_var <- col_uniq <= 1
sum(zero_var)  # should be 9
X_filtered <- X[, !zero_var, drop = FALSE]
dim(X_filtered) 

nb <- missMDA::estim_ncpPCA(X_filtered, method = "Regularized", scale = TRUE) # impute missing data
G_imp <- missMDA::imputePCA(X_filtered, ncp = nb$ncp) # impute missing data
G_imp <- G_imp$completeObs
G_imp[1:10,1:10] # overview of our data and its first 10 rows/10 columns

#for (i in 1:ncol(G))
#{
#  G[which(is.na(G[,i])),i] <- median(G[-which(is.na(G[,i])),i], na.rm=TRUE) # imputation
#}
#table(as.vector(G))
#dim(G)
#G[1:10,1:10] # overview of our data and its first 10 rows/10 columns
write.table(G_imp, "monq_sv_geno_imp_matrix.txt", sep="\t", col.names=F, row.names=F) # will be useful later

#### variant trimming ####
G2 <- add_code256(big_copy(G_imp,type = "raw"),code=bigsnpr:::CODE_012)
n <- nrow(G2)
alt_counts <- bigstatsr::big_colstats(G2)$sum
alleles    <- 2 * n
p_alt      <- alt_counts / alleles
maf        <- pmin(p_alt, 1 - p_alt)
mac        <- maf * alleles
maf_keep <- which(!is.na(maf) & maf >= 0.02 & mac >= 10) # Compute MAF & MAC on your imputed 0/1/2 data
all_idx   <- seq_len(ncol(G2))
to_exclude <- setdiff(all_idx, maf_keep)
#keep_clump <- snp_clumping(G = G2, infos.chr = as.integer(gsub("chr", "", chr_pos$chromosome)), infos.pos = chr_pos$position, exclude = to_exclude, thr.r2 = 0.2, size = 250, ncores = 64) # in kb (typical LD window) # LD clumping on the filtered set (same-ish to autoSVD defaults)
keep_clump <- snp_clumping(G = G2, infos.chr = infos.chr, infos.pos = infos.pos, exclude = to_exclude, thr.r2 = 0.2, size = 250, ncores = 64) # in kb (typical LD window) # LD clumping on the filtered set (same-ish to autoSVD defaults)
#tab_clump <- table(chromosome[keep_clump])
tab_clump <- table(infos.chr[keep_clump])
good_chr  <- names(tab_clump)[tab_clump >= 3]
#final_keep <- keep_clump[chromosome[keep_clump] %in% good_chr] # Drop chromosomes that would have < 3 variants AFTER clumping
final_keep <- keep_clump[infos.chr[keep_clump] %in% good_chr] # Drop chromosomes that would have < 3 variants AFTER clumping
#ord <- order(chromosome[final_keep], chr_pos$position[final_keep])
ord <- order(infos.chr[final_keep], infos.pos[final_keep])
ind_final <- final_keep[ord]
roll_size_safe <- 3L
size_safe      <- 51L
#newpc <- snp_autoSVD(G = G2, infos.chr = as.integer(gsub("chr", "", chr_pos$chromosome)), infos.pos = chr_pos$position, ind.col = ind_final, ncores = 64, size = size_safe, roll.size = roll_size_safe)
newpc <- snp_autoSVD(G = G2, infos.chr = infos.chr, infos.pos = infos.pos, ind.col = ind_final, ncores = 64, size = size_safe, roll.size = roll_size_safe)
which_pruned <- attr(newpc, which="subset") # Indexes of remaining SNPS after pruning
length(which_pruned)
#### variant trimming ####

G_hard <- floor(G_imp + 0.5)   # 0.49->0, 0.50–1.49->1, 1.50+->2
G_hard[G_hard < 0] <- 0
G_hard[G_hard > 2] <- 2
storage.mode(G_hard) <- "integer"  # ensure integer, not double

info_samples <- read.csv("monq_sample_info.csv", header=T)
head(info_samples)
match_idx  <- match(rownames(G_hard), info_samples$id)
pop_vector <- info_samples$k3[match_idx]
#my_fst <- MakeDiploidFSTMat(t(G), locusNames = id_sv, popNames = pop_vector) # FST matrix with OutFLANK
my_fst <- MakeDiploidFSTMat(G_hard, locusNames = chr_pos_filtered$id, popNames = pop_vector) # FST matrix with OutFLANK
out_trim <- OutFLANK(my_fst[which_pruned,], NumberOfSamples=3, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
str(out_trim)
OutFLANKResultsPlotter(out_trim, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =FALSE, RightZoomFraction = 0.05, titletext = NULL) # have a look at the results
hist(out_trim$results$pvaluesRightTail)
P1 <- pOutlierFinderChiSqNoCorr(my_fst, Fstbar = out_trim$FSTNoCorrbar, dfInferred = out_trim$dfInferred, qthreshold = 0.05, Hmin=0.1) # run OutFLANK on all SNPs of the data, corrected by the trim dataset
head(P1)

#P1_pos<-left_join(P1, chr_pos, by=c("LocusName"="id_sv"))
P1_pos<-left_join(P1, chr_pos_filtered, by=c("LocusName"="id"))
dim(P1_pos)
P1_pos[P1_pos$OutlierFlag==TRUE,] # no outlier

ggplot(P1_pos, aes(x=He, y=FST, colour=OutlierFlag)) + geom_point() + theme_classic() # look at the FSt as a function of heterozygosity to understand which snps have been evaluated, which one appear true or false outliers 
ggplot(P1_pos, aes(x=position/1000000, y=FST, colour=OutlierFlag))+ geom_point() + theme_classic() + facet_grid(cols = vars(chromosome), scales = "free_x", space="free_x") + labs(x = "position (in MB)") 
write.table(P1_pos, "sv_outflank_fst_outliers.txt", sep="\t", row.names=F, quote=F)

# PCAdapt -> missing data handled on the fly
G_pca <- read.pcadapt(G, type = "pcadapt")
G_pca[is.na(G_pca)] <- 9
#G_pca <- read.pcadapt(t(G_hard), type = "pcadapt")
x <- pcadapt(input = G_pca, K = 20) 
plot(x, option = "screeplot", scale = none) # K = 6 seems optiomal
#plot(x, option = "scores", Pop = info_samples$k2)
S <- scores(x)
df <- data.frame(PC1 = S[,1], PC2 = S[,2], pop = info_samples$k2)
#df <- data.frame(PC1 = S[,1], PC2 = S[,2], pop = info_samples$k2[match_idx])
ggplot(df, aes(PC1, PC2, color = pop)) + geom_point(size = 1.5) + coord_equal() + labs(x = "PC1", y = "PC2") + theme_minimal()
#plot(x, option = "scores", i = 3, j = 4, Pop = info_samples$k3)
df <- data.frame(PC3 = S[,3], PC4 = S[,4], pop = info_samples$k3)
#df <- data.frame(PC3 = S[,3], PC4 = S[,4], pop = info_samples$k3[match_idx])
ggplot(df, aes(PC3, PC4, color = pop)) + geom_point(size = 1.5) + coord_equal() + labs(x = "PC3", y = "PC4") + theme_minimal()
#plot(x, option = "scores", i = 4, j = 5, Pop = info_samples$k5)
df <- data.frame(PC4 = S[,4], PC5 = S[,5], pop = info_samples$k5)
#df <- data.frame(PC4 = S[,4], PC5 = S[,5], pop = info_samples$k5[match_idx])
ggplot(df, aes(PC4, PC5, color = pop)) + geom_point(size = 1.5) + coord_equal() + labs(x = "PC4", y = "PC5") + theme_minimal()

x <- pcadapt(G_pca, K = 6, min.maf = 0.01) # do using pruned SV set, too; but see below (LD affected minimally)
summary(x)
plot(x, option = "manhattan")
plot(x, option = "qqplot")
hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
plot(x, option = "stat.distribution")

par(mfrow = c(2, 2)) # evaluate if LD affected PCs - little in PC3
for (i in 1:4){
  plot(x$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))
}

qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.05
pca_outliers <- which(qval < alpha)
length(pca_outliers)

G_pca_ld <- read.pcadapt(G[which_pruned,], type = "pcadapt")
#G_pca_ld <- read.pcadapt(t(G_hard[,which_pruned]), type = "pcadapt")
x_ld <- pcadapt(input = G_pca_ld, K = 20) 
plot(x_ld, option = "screeplot") # K = 5 was optiomal
#plot(x_ld, option = "scores", Pop = info_samples$k2)
S <- scores(x_ld)
df <- data.frame(PC1 = S[,1], PC2 = S[,2], pop = info_samples$k2)
#df <- data.frame(PC1 = S[,1], PC2 = S[,2], pop = info_samples$k2[match_idx])
ggplot(df, aes(PC1, PC2, color = pop)) + geom_point(size = 1.5) + coord_equal() + labs(x = "PC1", y = "PC2") + theme_minimal()
#plot(x_ld, option = "scores", i = 3, j = 4, Pop = info_samples$k3)
df <- data.frame(PC3 = S[,3], PC4 = S[,4], pop = info_samples$k3)
#df <- data.frame(PC3 = S[,3], PC4 = S[,4], pop = info_samples$k3[match_idx])
ggplot(df, aes(PC3, PC4, color = pop)) + geom_point(size = 1.5) + coord_equal() + labs(x = "PC3", y = "PC4") + theme_minimal()
#plot(x_ld, option = "scores", i = 4, j = 5, Pop = info_samples$k5)
df <- data.frame(PC4 = S[,4], PC5 = S[,5], pop = info_samples$k5)
#df <- data.frame(PC4 = S[,4], PC5 = S[,5], pop = info_samples$k5[match_idx])
ggplot(df, aes(PC4, PC5, color = pop)) + geom_point(size = 1.5) + coord_equal() + labs(x = "PC4", y = "PC5") + theme_minimal()

x_ld <- pcadapt(G_pca_ld, K = 5, min.maf = 0.01) # pruned SV set looks clearer
summary(x_ld)
plot(x_ld, option = "manhattan")
plot(x_ld, option = "qqplot")
hist(x_ld$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
plot(x_ld, option = "stat.distribution")

qval_ld <- qvalue(x_ld$pvalues)$qvalues
alpha <- 0.05
pca_ld_outliers <- which(qval_ld < alpha)
length(pca_ld_outliers)

sv_pc <- get.pc(x, pca_outliers) # Association between PCs and outliers


# 0. Format genetic dataset
vcf <- read.vcfR("monq_mcSV.final.sort.vcf")  # Read VCF file
gt <- extract.gt(vcf, element = "GT", as.numeric = FALSE) # Extract genotype matrix

count_minor <- function(gt) { # Function to count minor alleles per genotype
  if (is.na(gt)) return(NA)
  alleles <- unlist(strsplit(gt, "[/|]"))
  return(sum(alleles == "1"))  # count the number of ALT alleles
}
allele_counts <- apply(gt, c(1, 2), count_minor) # Apply to matrix
allele_df <- as.data.frame(t(allele_counts))  # transpose to samples-as-rows
rownames(allele_df) <- colnames(gt)  # sample names
chrom <- vcf@fix[, "CHROM"] # Extract fields
pos   <- vcf@fix[, "POS"] # Extract fields
ref   <- vcf@fix[, "REF"] # Extract fields
alt   <- vcf@fix[, "ALT"] # Extract fields
site_key <- paste(chrom, pos, sep = "_") # Count how many variants share the same CHROM-POS-REF
site_counter <- ave(site_key, site_key, FUN = function(x) seq_along(x)) # Use a counter per duplicated site_key
variant_ids <- paste0(site_key, "_", site_counter) # Build the final unique ID
colnames(allele_df) <- variant_ids  # # Apply variant IDs to  allele count matrix
head(allele_df) # View example
write.csv(allele_df, file = "monq_MCsv.csv", quote = FALSE)
pca_outliers_adf <- allele_df[,pca_outliers]
write.csv(pca_outliers_adf, file = "monq_pcadapt_sv_final.csv", quote = FALSE)


# 1. Loading and formatting data 
# 1.1 Genetic dataset 
# Loading genetic dataset
#Genotypes <- read.csv("monq_MCsv.csv", header = T, row.names = 1)
Genotypes <- read.csv("monq_pcadapt_sv_final.csv", header = T, row.names = 1)

# Loading samples metadata
InfoInd <- read.csv("monq_sample_info.csv", header = T)
Genotypes <- Genotypes[match(InfoInd$id, row.names(Genotypes), nomatch = 0),]

# Estimating population allele frequencies
#AllFreq <- aggregate(Genotypes, by = list(InfoInd$k3), function(x) mean(x, na.rm = T)/2)
#row.names(AllFreq) <- as.character(AllFreq$Group.1)

# Imputing missing data
mean(is.na(Genotypes)) # 0.02047958 # proportion of missing values in the whole matrix

missing_ind <- rowMeans(is.na(Genotypes)) # proportion of missing per individual
summary(missing_ind)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.004814 0.008769 0.011348 0.031187 0.013755 1.000000
hist(missing_ind, main = "Missing rate per individual", xlab = "Proportion missing")

missing_loci <- colMeans(is.na(Genotypes)) # proportion of missing per locus
summary(missing_loci)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.01099 0.01099 0.02198 0.03119 0.04396 0.19780
hist(missing_loci, main = "Missing rate per locus", xlab = "Proportion missing")

Genotypes_clean <- Genotypes[missing_ind < 0.2, missing_loci < 0.1]  # keep individuals with <20% missing
sum(is.na(Genotypes_clean))               # total number of NAs
colSums(is.na(Genotypes_clean))           # missing count per locus
rowSums(is.na(Genotypes_clean))           # missing count per individual

X <- Genotypes_clean
zero_var <- vapply(X, function(v) {
  vals <- unique(v[!is.na(v)])
  length(vals) <= 1
}, logical(1))

sum(zero_var)              # how many columns
names_zero <- names(X)[zero_var]
head(names_zero)           # a peek at their names

X_zero <- X[, zero_var, drop = FALSE]
dim(X_zero)
head(X_zero) 

val_counts <- lapply(X_zero, function(v) {
  sort(table(v, useNA = "ifany"), decreasing = TRUE)
})
val_counts

# zero_var is TRUE for columns you want to remove
X_filtered <- X[, !zero_var, drop = FALSE]

# check dimensions before and after
dim(X)
dim(X_filtered)

nb <- missMDA::estim_ncpPCA(X_filtered, method = "Regularized", scale = TRUE) # impute missing data
Geno_imp <- missMDA::imputePCA(X_filtered, ncp = nb$ncp) # impute missing data

# Ordering loci based on their scaffold
Geno_imp <- Geno_imp$completeObs[,order(colnames(Geno_imp$completeObs))]
#Genotypes <- Genotypes[,order(colnames(Genotypes))]

# 1.2. Environmental variables
# Loading the environmental rasters
ras_present <- stack(list.files("./present", pattern = ".tif$", full.names = T))
names(ras_present) <- tools::file_path_sans_ext(basename(list.files("./present", pattern = ".tif$", full.names = T)))

ras_ssp126_2011_2040 <- stack(list.files("./ssp126_2011-2040", pattern = ".tif$", full.names = T))
names(ras_ssp126_2011_2040) <- tools::file_path_sans_ext(basename(list.files("./ssp126_2011-2040", pattern = ".tif$", full.names = T)))
ras_ssp126_2041_2070 <- stack(list.files("./ssp126_2041-2070", pattern = ".tif$", full.names = T))
names(ras_ssp126_2041_2070) <- tools::file_path_sans_ext(basename(list.files("./ssp126_2041-2070", pattern = ".tif$", full.names = T)))
ras_ssp126_2071_2100 <- stack(list.files("./ssp126_2071-2100", pattern = ".tif$", full.names = T))
names(ras_ssp126_2071_2100) <- tools::file_path_sans_ext(basename(list.files("./ssp126_2071-2100", pattern = ".tif$", full.names = T)))

ras_ssp585_2011_2040 <- stack(list.files("./ssp585_2011-2040", pattern = ".tif$", full.names = T))
names(ras_ssp585_2011_2040) <- tools::file_path_sans_ext(basename(list.files("./ssp585_2011-2040", pattern = ".tif$", full.names = T)))
ras_ssp585_2041_2070 <- stack(list.files("./ssp585_2041-2070", pattern = ".tif$", full.names = T))
names(ras_ssp585_2041_2070) <- tools::file_path_sans_ext(basename(list.files("./ssp585_2041-2070", pattern = ".tif$", full.names = T)))
ras_ssp585_2071_2100 <- stack(list.files("./ssp585_2071-2100", pattern = ".tif$", full.names = T))
names(ras_ssp585_2071_2100) <- tools::file_path_sans_ext(basename(list.files("./ssp585_2071-2100", pattern = ".tif$", full.names = T)))

# A perfect match between non missing pixels in the different climatic rasters is required to facilitate some of the following analyses. 
# To homogenize the different time periods we use the remove.NAs.stack function.
remove.NAs.stack<-function(rast.stack){
  nom<-names(rast.stack)
  test1<-calc(rast.stack, fun=sum)
  test1[!is.na(test1)]<-1
  test2<-rast.stack*test1
  test2<-stack(test2)
  names(test2)<-nom
  return(test2)
}
ras_present <- remove.NAs.stack(ras_present)

ras_ssp126_2011_2040 <- remove.NAs.stack(ras_ssp126_2011_2040)
ras_ssp126_2041_2070 <- remove.NAs.stack(ras_ssp126_2041_2070)
ras_ssp126_2071_2100 <- remove.NAs.stack(ras_ssp126_2071_2100)
ras_ssp585_2011_2040 <- remove.NAs.stack(ras_ssp585_2011_2040)
ras_ssp585_2041_2070 <- remove.NAs.stack(ras_ssp585_2041_2070)
ras_ssp585_2071_2100 <- remove.NAs.stack(ras_ssp585_2071_2100)

# Extracting environmental values for each source population
Env <- data.frame(extract(ras_present, InfoInd[,2:3])) # 2:3 = x,y

# Standardization of the variables
Env <- scale(Env, center=TRUE, scale=TRUE) # center=TRUE, scale=TRUE are the defaults for scale()

# Recovering scaling coefficients
scale_env <- attr(Env, 'scaled:scale')
center_env <- attr(Env, 'scaled:center')

# Climatic table
Env <- as.data.frame(Env)
row.names(Env) <- c(InfoInd$id)

Env_filt <- Env[rownames(Env) %in% rownames(Geno_imp), ]

# 1.3. Phenotypic traits 

# 1.4. Inferring population structure -> update later using neutral SNP only
# Loading PCAngsd result
C <- as.matrix(read.table("monq_pan_snp_pca.cov"))
e <- eigen(C)
pca <- prcomp(C, scale = TRUE)
screeplot(pca, type = "barplot", npcs=10, main="PCA Eigenvalues")

# Population structure table
PCs <- scores(pca, choices=c(1:2), display="sites", scaling=0) # Based on the screeplot, two PCs would be a reasonable set to retain as a proxy for neutral population structure in downstream analyses
PopStruct <- data.frame(Individual = InfoInd$id, PCs)
colnames(PopStruct) <- c("Individual", "PC1", "PC2")
PopStruct_filt <- PopStruct[PopStruct$Individual %in% rownames(Geno_imp), ] 

# 1.5. Merging all variables into a single table + loading mapping features 
# Table gathering all variables
InfoInd_filt <- InfoInd[InfoInd$id %in% rownames(Geno_imp), ]
Variables <- data.frame(InfoInd_filt[,1:3], PopStruct_filt[,-1], Env_filt)
#write.csv(Variables, file = "merged_variables_final.csv", quote = FALSE)
#Variables <- read.csv("merged_variables_final.csv", header = T, row.names = 1)

# Species range shapefile (download information at beginning of tutorial)
range <- st_read("monq_450km_90y_onLand.shp") 
admin <- ne_countries(scale = "medium", returnclass = "sf")

# 2. Variable selection
# Null model
RDA0 <- rda(Geno_imp ~ 1,  Variables) 

# Full model
RDAfull <- rda(Geno_imp ~ bio2 + bio5 + bio7 + bio8 + bio9 + bio14 + bio15 + bio18 + aspect + slope, Variables)

# Stepwise procedure with ordiR2step function
mod <- ordiR2step(RDA0, RDAfull, Pin = 0.05, R2permutations = 1000, R2scope = T) # stopping criteria: variable significance of p < 0.01 using 1000 permutations, and the adjusted R2 of the global model
## Result (pca outliers)
#Step: R2.adj= 0.1656966
#Call: Geno_imp ~ bio15 + bio7 + bio9 + bio8


# 3. Variance partitioning: disentangling the drivers of genetic variation 
# Full model
pRDAfull <- rda(Geno_imp ~ bio2 + bio5 + bio7 + bio8 + bio9 + bio14 + bio15 + bio18 + aspect + slope + PC1 + PC2 + x + y, Variables)
#pRDAfull_subset <- rda(Geno_imp ~ bio2 + bio7 + bio8 + bio9 + PC1 + PC2 + x + y, Variables) # selected variables from above for pca outliers
pRDAfull_subset <- rda(Geno_imp ~ bio7 + bio8 + bio9 + bio15 + PC1 + PC2 + x + y, Variables) # selected variables from above for pca outliers
RsquareAdj(pRDAfull)
#$r.squared
#[1] 0.3242274

#$adj.r.squared
#[1] 0.1946271

RsquareAdj(pRDAfull_subset)
#$r.squared
#[1] 0.2624335

#$adj.r.squared
#[1] 0.1877433

anova(pRDAfull)
#Permutation test for rda under reduced model
#Permutation: free
#Number of permutations: 999

#Model: rda(formula = Geno_imp ~ bio2 + bio5 + bio7 + bio8 + bio9 + bio14 + bio15 + bio18 + aspect + slope + PC1 + PC2 + x + y, data = Variables)
#         Df Variance      F Pr(>F)
#Model    14   92.457 2.5017  0.001 ***
#Residual 73  192.705
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(pRDAfull_subset)
#Permutation test for rda under reduced model
#Permutation: free
#Number of permutations: 999

#Model: rda(formula = Geno_imp ~ bio7 + bio8 + bio9 + bio15 + PC1 + PC2 + x + y, data = Variables)
#         Df Variance      F Pr(>F)
#Model     8   74.836 3.5136  0.001 ***
#Residual 79  210.326
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Pure environment model
pRDAenv <- rda(Geno_imp ~ bio2 + bio5 + bio7 + bio8 + bio9 + bio14 + bio15 + bio18 + aspect + slope + Condition(PC1 + PC2 + x + y),  Variables)
#pRDAenv_subset <- rda(Geno_imp ~ bio2 + bio7 + bio8 + bio9 + Condition(PC1 + PC2 + x + y),  Variables)
pRDAenv_subset <- rda(Geno_imp ~ bio7 + bio8 + bio9 + bio15 + Condition(PC1 + PC2 + x + y),  Variables)
RsquareAdj(pRDAenv)
#$r.squared
#[1] 0.1143315

#$adj.r.squared
#[1] 0.02280854

RsquareAdj(pRDAenv_subset)
#$r.squared
#[1] 0.05253763

#$adj.r.squared
#[1] 0.01592466

anova(pRDAenv)
#Permutation test for rda under reduced model
#Permutation: free
#Number of permutations: 999

#Model: rda(formula = Geno_imp ~ bio2 + bio5 + bio7 + bio8 + bio9 + bio14 + bio15 + bio18 + aspect + slope + Condition(PC1 + PC2 + x + y), data = Variables)
#         Df Variance      F Pr(>F)
#Model    10   32.603 1.2351  0.001 ***
#Residual 73  192.705
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(pRDAenv_subset)
#Permutation test for rda under reduced model
#Permutation: free
#Number of permutations: 999

#Model: rda(formula = Geno_imp ~ bio7 + bio8 + bio9 + bio15 + Condition(PC1 + PC2 + x + y), data = Variables)
#         Df Variance      F Pr(>F)
#Model     4   14.982 1.4068  0.001 ***
#Residual 79  210.326
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Pure climate model
pRDAclim <- rda(Geno_imp ~ bio2 + bio5 + bio7 + bio8 + bio9 + bio14 + bio15 + bio18 + Condition(aspect + slope + PC1 + PC2 + x + y),  Variables)
#pRDAclim_subset <- rda(Geno_imp ~ bio2 + bio7 + bio8 + bio9 + Condition(PC1 + PC2 + x + y),  Variables)
pRDAclim_subset <- rda(Geno_imp ~ bio7 + bio8 + bio9 + bio15 + Condition(PC1 + PC2 + x + y),  Variables)
RsquareAdj(pRDAclim)
#$r.squared
#[1] 0.09558288

#$adj.r.squared
#[1] 0.02312009 # the highest among the examined -> account for non-bioclim, genetic ancestry and recent IBD especially for the scattered sampling scheme

RsquareAdj(pRDAclim_subset)
#$r.squared
#[1] 0.05253763

#$adj.r.squared
#[1] 0.01592466

anova(pRDAclim)
#Permutation test for rda under reduced model
#Permutation: free
#Number of permutations: 999

#Model: rda(formula = Geno_imp ~ bio2 + bio5 + bio7 + bio8 + bio9 + bio14 + bio15 + bio18 + Condition(aspect + slope + PC1 + PC2 + x + y), data = Variables)
#         Df Variance      F Pr(>F)
#Model     8   27.257 1.2907  0.001 ***
#Residual 73  192.705
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(pRDAclim_subset)
#Permutation test for rda under reduced model
#Permutation: free
#Number of permutations: 999

#Model: rda(formula = Geno_imp ~ bio7 + bio8 + bio9 + bio15 + Condition(PC1 + PC2 + x + y), data = Variables)
#         Df Variance      F Pr(>F)
#Model     4   14.982 1.4068  0.001 ***
#Residual 79  210.326
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

pRDAclim2 <- rda(Geno_imp ~ bio2 + bio5 + bio7 + bio8 + bio9 + bio14 + bio15 + bio18 + Condition(PC1 + PC2 + x + y),  Variables)
RsquareAdj(pRDAclim2)
#$r.squared
#[1] 0.09534688

#$adj.r.squared
#[1] 0.02226303 

anova(pRDAclim2)
#Permutation test for rda under reduced model
#Permutation: free
#Number of permutations: 999

#Model: rda(formula = Geno_imp ~ bio2 + bio5 + bio7 + bio8 + bio9 + bio14 + bio15 + bio18 + Condition(PC1 + PC2 + x + y), data = Variables)
#         Df Variance      F Pr(>F)
#Model     8   27.189 1.2866  0.001 ***
#Residual 75  198.118
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

pRDAclim3 <- rda(Geno_imp ~ bio2 + bio5 + bio7 + bio8 + bio9 + bio14 + bio15 + bio18 + Condition(PC1 + PC2),  Variables)
RsquareAdj(pRDAclim3)
#$r.squared
#[1] 0.09502168

#$adj.r.squared
#[1] 0.02115863

anova(pRDAclim3)
#Permutation test for rda under reduced model
#Permutation: free
#Number of permutations: 999

#Model: rda(formula = Geno_imp ~ bio2 + bio5 + bio7 + bio8 + bio9 + bio14 + bio15 + bio18 + Condition(PC1 + PC2), data = Variables)
#         Df Variance     F Pr(>F)
#Model     8   27.097 1.278  0.001 ***
#Residual 77  204.066
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

pRDAclim4 <- rda(Geno_imp ~ bio2 + bio5 + bio7 + bio8 + bio9 + bio14 + bio15 + bio18 + Condition(x + y),  Variables)
RsquareAdj(pRDAclim4)
#$r.squared
#[1] 0.1102026

#$adj.r.squared
#[1] 0.03563286 # higher, but conditioning geocoordinates cannot account for genetic ancestry

anova(pRDAclim4)
#Permutation test for rda under reduced model
#Permutation: free
#Number of permutations: 999

#Model: rda(formula = Geno_imp ~ bio2 + bio5 + bio7 + bio8 + bio9 + bio14 + bio15 + bio18 + Condition(x + y), data = Variables)
#         Df Variance      F Pr(>F)
#Model     8   31.426 1.4618  0.001 ***
#Residual 77  206.919
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
pRDAclim5 <- rda(Geno_imp ~ bio2 + bio5 + bio7 + bio8 + bio9 + bio14 + bio15 + bio18 + Condition(aspect + slope),  Variables)
RsquareAdj(pRDAclim5)
#$r.squared
#[1] 0.2431062

#$adj.r.squared
#[1] 0.1710922

anova(pRDAclim5)
#Permutation test for rda under reduced model
#Permutation: free
#Number of permutations: 999

#Model: rda(formula = Geno_imp ~ bio2 + bio5 + bio7 + bio8 + bio9 + bio14 + bio15 + bio18 + Condition(aspect + slope), data = Variables)
#         Df Variance     F Pr(>F)
#Model     8   69.325 3.201  0.001 ***
#Residual 77  208.451
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Pure neutral population structure model 
pRDAstruct <- rda(Geno_imp ~ PC1 + PC2 + Condition(x + y + bio2 + bio5 + bio7 + bio8 + bio9 + bio14 + bio15 + bio18 + aspect + slope),  Variables)
pRDAstruct_subset <- rda(Geno_imp ~ PC1 + PC2 + Condition(x + y + bio7 + bio8 + bio9 + bio15),  Variables)
RsquareAdj(pRDAstruct)
#$r.squared
#[1] 0.0312126

#$adj.r.squared
#[1] 0.01473

RsquareAdj(pRDAstruct_subset)
#$r.squared
#[1] 0.03463372

#$adj.r.squared
#[1] 0.01714346

anova(pRDAstruct)
#Permutation test for rda under reduced model
#Permutation: free
#Number of permutations: 999

#Model: rda(formula = Geno_imp ~ PC1 + PC2 + Condition(x + y + bio2 + bio5 + bio7 + bio8 + bio9 + bio14 + bio15 + bio18 + aspect + slope), data = Variables)
#         Df Variance      F Pr(>F)
#Model     2    8.901 1.6859  0.001 ***
#Residual 73  192.705
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(pRDAstruct_subset)
#Permutation test for rda under reduced model
#Permutation: free
#Number of permutations: 999

#Model: rda(formula = Geno_imp ~ PC1 + PC2 + Condition(x + y + bio7 + bio8 + bio9 + bio15), data = Variables)
#         Df Variance      F Pr(>F)
#Model     2    9.876 1.8548  0.001 ***
#Residual 79  210.326
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Pure geography model 
pRDAgeog <- rda(Geno_imp ~ x + y + Condition(bio2 + bio5 + bio7 + bio8 + bio9 + bio14 + bio15 + bio18 + aspect + slope + PC1 + PC2),  Variables)
pRDAgeog_subset <- rda(Geno_imp ~ x + y + Condition(bio7 + bio8 + bio9 + bio15 + PC1 + PC2),  Variables)
RsquareAdj(pRDAgeog)
#$r.squared
#[1] 0.02050491

#$adj.r.squared
#[1] 0.002309083

RsquareAdj(pRDAgeog_subset)
#$r.squared
#[1] 0.0212079

#$adj.r.squared
#[1] 0.00272313

anova(pRDAgeog)
#Permutation test for rda under reduced model
#Permutation: free
#Number of permutations: 999

#Model: rda(formula = Geno_imp ~ x + y + Condition(bio2 + bio5 + bio7 + bio8 + bio9 + bio14 + bio15 + bio18 + aspect + slope + PC1 + PC2), data = Variables)
#         Df Variance      F Pr(>F)
#Model     2    5.847 1.1075  0.032 *
#Residual 73  192.705
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(pRDAgeog_subset)
#Permutation test for rda under reduced model
#Permutation: free
#Number of permutations: 999

#Model: rda(formula = Geno_imp ~ x + y + Condition(bio7 + bio8 + bio9 + bio15 + PC1 + PC2), data = Variables)
#         Df Variance      F Pr(>F)
#Model     2    6.048 1.1358   0.02 *
#Residual 79  210.326
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Pure topography model
pRDAtopo <- rda(Geno_imp ~ aspect + slope + Condition(bio2 + bio5 + bio7 + bio8 + bio9 + bio14 + bio15 + bio18 + PC1 + PC2 + x + y),  Variables)
pRDAtopo_subset <- rda(Geno_imp ~ aspect + slope + Condition(bio7 + bio8 + bio9 + bio15 + PC1 + PC2 + x + y),  Variables)
RsquareAdj(pRDAtopo)
#$r.squared
#[1] 0.01898459

#$adj.r.squared
#[1] 0.0005455126

RsquareAdj(pRDAtopo_subset)
#$r.squared
#[1] 0.01868501

#$adj.r.squared
#[1] 1.405688e-05

anova(pRDAtopo)
#Permutation test for rda under reduced model
#Permutation: free
#Number of permutations: 999

#Model: rda(formula = Geno_imp ~ aspect + slope + Condition(bio2 + bio5 + bio7 + bio8 + bio9 + bio14 + bio15 + bio18 + PC1 + PC2 + x + y), data = Variables)
#         Df Variance      F Pr(>F)
#Model     2    5.414 1.0254  0.332
#Residual 73  192.705

anova(pRDAtopo_subset)
#Permutation test for rda under reduced model
#Permutation: free
#Number of permutations: 999

#Model: rda(formula = Geno_imp ~ aspect + slope + Condition(bio7 + bio8 + bio9 + bio15 + PC1 + PC2 + x + y), data = Variables)
#         Df Variance      F Pr(>F)
#Model     2    5.328 1.0007  0.471
#Residual 77  204.998


# 4. Genotype-Environment Associations: identifying loci under selection 
# 4.1. Conducting the genome scan using pRDA 
RDA_env <- rda(Geno_imp ~ bio2 + bio5 + bio7 + bio8 + bio9 + bio14 + bio15 + bio18 + aspect + slope + Condition(PC1 + PC2 + x + y), Variables)
RDA_env_subset <- rda(Geno_imp ~ bio7 + bio8 + bio9 + bio15 + Condition(PC1 + PC2 + x + y), Variables)
screeplot(RDA_env, main="Eigenvalues of constrained axes")
screeplot(RDA_env_subset, main="Eigenvalues of constrained axes")
summary(eigenvals(RDA_env, model = "constrained"))
#Importance of components:
#                        RDA1   RDA2   RDA3    RDA4    RDA5    RDA6    RDA7
#Eigenvalue            7.4905 4.6416 3.5459 3.03967 2.73958 2.63826 2.42425
#Proportion Explained  0.2297 0.1424 0.1088 0.09323 0.08403 0.08092 0.07436
#Cumulative Proportion 0.2297 0.3721 0.4809 0.57411 0.65814 0.73906 0.81342
#                         RDA8    RDA9   RDA10
#Eigenvalue            2.32668 1.98510 1.77139
#Proportion Explained  0.07136 0.06089 0.05433
#Cumulative Proportion 0.88478 0.94567 1.00000

summary(eigenvals(RDA_env_subset, model = "constrained"))
#Importance of components:
#                        RDA1  RDA2   RDA3   RDA4
#Eigenvalue            5.6940 3.430 3.0128 2.8445
#Proportion Explained  0.3801 0.229 0.2011 0.1899
#Cumulative Proportion 0.3801 0.609 0.8101 1.0000

signif.full <- anova.cca(RDA_env, parallel=getOption("mc.cores")) # default is permutation=999
signif.full
#Permutation test for rda under reduced model
#Permutation: free
#Number of permutations: 999

#Model: rda(formula = Geno_imp ~ bio2 + bio5 + bio7 + bio8 + bio9 + bio14 + bio15 + bio18 + aspect + slope + Condition(PC1 + PC2 + x + y), data = Variables)
#         Df Variance      F Pr(>F)
#Model    10   32.603 1.2351  0.001 ***
#Residual 73  192.705
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

signif.axis <- anova.cca(RDA_env, by="axis", parallel=getOption("mc.cores"))
signif.axis
#Permutation test for rda under reduced model
#Forward tests for axes
#Permutation: free
#Number of permutations: 999

#Model: rda(formula = Geno_imp ~ bio2 + bio5 + bio7 + bio8 + bio9 + bio14 + bio15 + bio18 + aspect + slope + Condition(PC1 + PC2 + x + y), data = Variables)
#         Df Variance      F Pr(>F)
#RDA1      1    7.491 2.8375  0.001 ***
#RDA2      1    4.642 1.7824  0.003 **
#RDA3      1    3.546 1.3801  0.352
#RDA4      1    3.040 1.1988  0.878
#RDA5      1    2.740 1.0947  0.984
#RDA6      1    2.638 1.0679  0.984
#RDA7      1    2.424 0.9938  0.996
#RDA8      1    2.327 0.9659  0.996
#RDA9      1    1.985 0.8344  0.998
#RDA10     1    1.771 0.7538  0.998
#Residual 73  192.705
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

RsquareAdj(RDA_env)
#$r.squared
#[1] 0.1143315

#$adj.r.squared
#[1] 0.02280854

# Function rdadapt
source("rdadapt.R")

# Running the function with K = 2
rdadapt_env <- rdadapt(RDA_env, 2)
rdadapt_env_subset <- rdadapt(RDA_env_subset, 2)

# P-values threshold after Bonferroni correction
#thres_env <- 0.01/length(rdadapt_env$p.values)
thres_env <- 0.05

# Identifying the loci that are below the p-value threshold
outliers <- data.frame(Loci = colnames(Geno_imp)[which(rdadapt_env$q.values<thres_env)], q.value = rdadapt_env$q.values[which(rdadapt_env$q.values<thres_env)], contig = unlist(lapply(strsplit(colnames(Geno_imp)[which(rdadapt_env$q.values<thres_env)], split = "_"), function(x) x[1])))
outliers_subset <- data.frame(Loci = colnames(Geno_imp)[which(rdadapt_env_subset$q.values<thres_env)], q.value = rdadapt_env_subset$q.values[which(rdadapt_env_subset$q.values<thres_env)], contig = unlist(lapply(strsplit(colnames(Geno_imp)[which(rdadapt_env_subset$q.values<thres_env)], split = "_"), function(x) x[1])))
intersects <- intersect(as.character(unlist(outliers[1])), as.character(unlist(outliers_subset[1])))
full_only <- setdiff(outliers[1], outliers_subset[1])
subset_only <- setdiff(outliers_subset[1], outliers[1])
nU <- length(Geno_imp); nF <- nrow(outliers); nS <- nrow(outliers_subset); nI <- length(intersects)
jaccard <- nI / length(union(full_only, subset_only)) # 146
precision_full  <- nI / nF           # 0.5309091; of full’s outliers, how many replicate?
recall_full     <- nI / nS           # 0.8156425; symmetry check relative to subset
mat <- matrix(c(nI, nF - nI, nS - nI, nU - nF - nS + nI), nrow = 2)
fisher_res <- fisher.test(mat, alternative = "greater")  # enrichment of overlap
fisher_res$p.value # 0; Significant overlap → strong concordance; proceed with full set

outliers <- outliers[order(outliers$contig, outliers$q.value),]
#saveRDS(outliers, file = "RDA_outliers_final.rds")

foo <- matrix(nrow=(length(outliers$Loci)), ncol=10)  # 10 columns for 10 predictors
colnames(foo) <- c("aspect","bio14", "bio15", "bio18", "bio2","bio5","bio7","bio8","bio9","slope")

for (i in 1:length(outliers$Loci)) {
  nam <- outliers$Loci[i]
  sv.gen <- Geno_imp[,nam]
  foo[i,] <- apply(Env_filt,2,function(x) cor(x,sv.gen))
}

outliers_env <- cbind.data.frame(outliers,foo)  
head(outliers_env)

outliers_env_uniq <- outliers_env[!duplicated(outliers_env$Loci),] # the same

for (i in 1:length(outliers_env_uniq$Loci)) {
  bar <- outliers_env_uniq[i,]
  outliers_env_uniq[i,14] <- names(which.max(abs(bar[4:13]))) # gives the name of the variable with the maximum correlation (out of 6) 
  outliers_env_uniq[i,15] <- max(abs(bar[4:13]))              # gives the maximum correlation value
}

colnames(outliers_env_uniq)[14] <- "predictor"
colnames(outliers_env_uniq)[15] <- "correlation"

#saveRDS(outliers_env_uniq, file = "RDA_outliers_wPredictor_final.rds")

# List of Top hit outlier names per contig
outliers_rdadapt_env <- as.character(outliers$Loci[!duplicated(outliers$contig)])

# Plotting
sel <- outliers_env_uniq$Loci
env <- outliers_env_uniq$predictor
env[env=="aspect"] <- '#fb9a99'
env[env=="bio14"] <- '#33a02c' 
env[env=="bio15"] <- '#1f78b4'
env[env=="bio18"] <- '#a6cee3'  
env[env=="bio2"] <- '#6a3d9a'
env[env=="bio5"] <- '#ffff33'
env[env=="bio7"] <- '#b2df8a'
env[env=="bio8"] <- '#e31a1c'
env[env=="bio9"] <- '#bf5b17'
env[env=="slope"] <- '#CC79A7'

# color by predictor:
col.pred <- rownames(RDA_env$CCA$v) # pull the variant names

for (i in 1:length(sel)) {           # color code candidate variants
  foo <- match(sel[i], col.pred)
  col.pred[foo] <- env[i]
}

col.pred[grep("chr", col.pred)] <- '#f1eef6' # non-candidate variants
empty <- col.pred
empty[grep("#f1eef6", empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty == "#00FF0000", "#00FF0000", "gray32")
bg <- c('#fb9a99', '#33a02c', '#1f78b4', '#a6cee3', '#6a3d9a', '#ffff33', '#b2df8a', '#e31a1c', '#bf5b17', '#CC79A7')

# axes 1 & 2
svg("rda_plot_wPredictor.svg", width = 8, height = 8)
plot(RDA_env, type = "n", scaling = 3, xlim = c(-1,1), ylim = c(-1,1))
points(RDA_env, display = "species", pch = 21, cex = 1, col = "gray32", bg = col.pred, scaling = 3)
points(RDA_env, display = "species", pch = 21, cex = 1, col = empty.outline, bg = empty, scaling = 3)
text(RDA_env, scaling = 3, display = "bp", col = "black", cex = 1)
legend("bottomright", legend = c("aspect", "bio14", "bio15", "bio18", "bio2", "bio5", "bio7", "bio8", "bio9", "slope"), bty = "n", col = "gray32", pch = 21, cex = 1, pt.bg = bg)
dev.off()

pdf("rda_plot_wPredictor.pdf", width = 8, height = 8)
plot(RDA_env, type = "n", scaling = 3, xlim = c(-1,1), ylim = c(-1,1))
points(RDA_env, display = "species", pch = 21, cex = 1, col = "gray32", bg = col.pred, scaling = 3)
points(RDA_env, display = "species", pch = 21, cex = 1, col = empty.outline, bg = empty, scaling = 3)
text(RDA_env, scaling = 3, display = "bp", col = "black", cex = 1)
legend("bottomright", legend = c("aspect", "bio14", "bio15", "bio18", "bio2", "bio5", "bio7", "bio8", "bio9", "slope"), bty = "n", col = "gray32", pch = 21, cex = 1, pt.bg = bg)
dev.off()

# Formatting table for ggplot 
locus_scores <- scores(RDA_env, choices = c(1:2), display = "species", scaling = "none") # vegan references "species", here these are the loci
TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores)
TAB_loci$type <- "Neutral"
TAB_loci$type[TAB_loci$names%in%outliers$Loci] <- "All outliers"
TAB_loci$type[TAB_loci$names%in%outliers_rdadapt_env] <- "Top outliers"
TAB_loci$type <- factor(TAB_loci$type, levels = c("Neutral", "All outliers", "Top outliers"))
TAB_loci <- TAB_loci[order(TAB_loci$type),]
TAB_var <- as.data.frame(scores(RDA_env, choices = c(1,2), display = "bp")) # pull the biplot scores

top_outlier_names <- sort(outliers_rdadapt_env) 
TAB_loci$shape_id <- NA
TAB_loci$shade <- NA
n <- length(top_outlier_names) # 16
shape_ids <- rep(0:15, length.out = n)
shades <- rep(c("wine", "purple"), length.out = n)
shape_lookup <- setNames(shape_ids, top_outlier_names)
shade_lookup <- setNames(shades, top_outlier_names)

TAB_loci$shape_id <- as.character(shape_lookup[TAB_loci$names])
TAB_loci$shade <- as.character(shade_lookup[TAB_loci$names])
TAB_loci$shape_id <- as.factor(TAB_loci$shape_id)
TAB_loci$shade <- factor(TAB_loci$shade, levels = c("wine", "purple"))
TAB_loci$color_group <- as.character(TAB_loci$type)
TAB_loci$color_group[!is.na(TAB_loci$shape_id)] <- paste("Top", TAB_loci$shade[!is.na(TAB_loci$shape_id)])
TAB_loci$color_group <- factor(TAB_loci$color_group, levels = c("Neutral", "All outliers", "Top wine", "Top purple"))

# Biplot of RDA loci and variables scores
p <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray80", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray80", linewidth = 0.5) +
  geom_point(data = TAB_loci[is.na(TAB_loci$shape_id), ], aes(x = RDA1 * 20, y = RDA2 * 20, colour = color_group), size = 1.2, alpha = 0.85) +
  geom_point(data = TAB_loci[!is.na(TAB_loci$shape_id), ], aes(x = RDA1 * 20, y = RDA2 * 20, shape = shape_id, color = color_group), size = 2.5, stroke = 1, alpha = 0.85, show.legend = c(shape = FALSE, color = TRUE)) +
  scale_color_manual(values = c("Neutral" = "gray90", "All outliers" = "#F9A242", "Top wine" = "#7F264B", "Top purple" = "#6B4596"), breaks = c("Neutral", "All outliers", "Top wine", "Top purple"), labels = c("Neutral", "All outliers", "Top chromosomal outliers", "Top chromosomal outliers"), name = "Locus type") +
  scale_shape_manual(values = 0:15, na.translate = FALSE) +
  geom_segment(data = TAB_var, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), colour="black", linewidth = 0.3, linetype = 1, arrow = arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x = 1.1 * RDA1, y = 1.1 * RDA2, label = row.names(TAB_var)), size = 3, family = "Helvetica") +
  xlab("RDA 1") + ylab("RDA 2") +
  facet_wrap(~"RDA space") +
  guides(color=guide_legend(title = "Locus type")) +
  guides(shape = "none") +
  theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(panel.grid = element_blank(), axis.ticks = element_blank(), legend.key = element_blank(), panel.background = element_blank(), legend.background = element_blank(), legend.title = element_text(size = 10), legend.text = element_text(size = 9), 
  plot.background = element_blank(), strip.text = element_text(size = 11, face = "bold"), plot.margin = margin(10, 10, 10, 10), aspect.ratio = 1)
ggsave("rda_plot.svg", plot = p, width = 8, height = 8, units = "in", device = "svg")
ggsave("rda_plot.pdf", plot = p, width = 8, height = 8, units = "in", device = cairo_pdf)

# Manhattan plot
Outliers <- rep("Neutral", length(colnames(Geno_imp)))
Outliers[colnames(Geno_imp)%in%outliers$Loci] <- "All outliers"
Outliers[colnames(Geno_imp)%in%outliers_rdadapt_env] <- "Top outliers"
Outliers <- factor(Outliers, levels = c("Neutral", "All outliers", "Top outliers"))

chr_names <- gsub("chr", "", sapply(strsplit(colnames(Geno_imp), "_"), `[`, 1))

TAB_manhatan <- data.frame(#pos = 1:length(colnames(Geno_imp)), 
                           Chrom = chr_names,
                           qvalues = rdadapt_env$q.values, 
                           Outliers = Outliers,
                           stringsAsFactors = FALSE)
TAB_manhatan$Chrom <- factor(TAB_manhatan$Chrom, levels = mixedsort(unique(TAB_manhatan$Chrom)))                           
top_outlier_names <- colnames(Geno_imp)[colnames(Geno_imp) %in% outliers_rdadapt_env]
TAB_manhatan$shape_id <- NA
TAB_manhatan$shade <- NA
shape_ids <- rep(0:15, length.out = n) # n= 16
shades <- rep(c("wine", "purple"), length.out = n)
TAB_manhatan$shape_id <- as.character(shape_lookup[colnames(Geno_imp)])
TAB_manhatan$shade <- as.character(shade_lookup[colnames(Geno_imp)])
TAB_manhatan$shape_id <- as.factor(TAB_manhatan$shape_id)
TAB_manhatan$shade <- factor(TAB_manhatan$shade, levels = c("wine", "purple"))
TAB_manhatan$color_group <- as.character(TAB_manhatan$Outliers)
chr_levels <- levels(TAB_manhatan$Chrom)
for (i in seq_along(chr_levels)) {
  chrom <- chr_levels[i]
  grey_color <- ifelse(i %% 2 == 0, "gray85", "gray30")
  idx <- TAB_manhatan$Chrom == chrom & TAB_manhatan$Outliers == "Neutral"
  TAB_manhatan$color_group[idx] <- grey_color
}
TAB_manhatan$color_group[!is.na(TAB_manhatan$shape_id)] <- paste("Top", TAB_manhatan$shade[!is.na(TAB_manhatan$shape_id)])
TAB_manhatan$color_group <- factor(TAB_manhatan$color_group, levels = c("gray85", "gray30", "All outliers", "Top wine", "Top purple"))
TAB_manhatan <- TAB_manhatan[order(TAB_manhatan$Chrom), ]
TAB_manhatan$pos <- 1:nrow(TAB_manhatan)

keep_chroms <- c(as.character(1:10), "29", "30")
chr_ticks <- aggregate(pos ~ Chrom, data = TAB_manhatan, FUN = function(x) mean(range(x)))
chr_ticks <- chr_ticks[chr_ticks$Chrom %in% keep_chroms, ]
x_ticks <- chr_ticks$pos
x_labels <- chr_ticks$Chrom

p <- ggplot(data = TAB_manhatan) +
  geom_point(data = TAB_manhatan[is.na(TAB_manhatan$shape_id), ], aes(x = pos, y = -log10(qvalues), color = color_group), size = 1.1, alpha = 0.85) +
  geom_point(data = TAB_manhatan[!is.na(TAB_manhatan$shape_id), ], aes(x = pos, y = -log10(qvalues), color = color_group, shape = shape_id), size = 2.2, stroke = 1, alpha = 0.85, show.legend = c(color = TRUE, shape = FALSE)) +
  scale_color_manual(values = c("gray85" = "gray85", "gray30" = "gray30", "All outliers" = "#F9A242", "Top wine" = "#7F264B", "Top purple" = "#6B4596"), breaks = c("gray85", "All outliers", "Top wine", "Top purple"), labels = c("Neutral", "All outliers", "Top outliers", "Top outliers"), name = "Locus type") +
  scale_shape_manual(values = 0:15, na.translate = FALSE) +
  scale_x_continuous(breaks = x_ticks, labels = x_labels) +
  xlab("Chromosome") + ylab(expression(-log[10](q))) +
  geom_hline(yintercept = -log10(thres_env), linetype = "dashed", color = "black", linewidth = 0.5) +
  facet_wrap(~"Manhattan plot", nrow = 3, strip.position = "top") +
  guides(color = guide_legend(title = "Locus type")) +
  guides(shape = "none") +
  theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(panel.grid = element_blank(), axis.ticks = element_blank(), strip.text = element_text(size = 11, face = "bold"), legend.position = "right", legend.background = element_blank(), legend.key = element_blank(), legend.box.background = element_blank(), 
  legend.title = element_text(size = 10), legend.text = element_text(size = 9), plot.background = element_blank(), panel.background = element_blank(), panel.spacing = unit(1, "lines"), plot.margin = margin(10, 10, 10, 10), aspect.ratio = 0.4)
ggsave("manhattan_plot.pdf", plot = p, width = 10, height = 8, units = "in", device = cairo_pdf)
ggsave("manhattan_plot.svg", plot = p, width = 10, height = 8, units = "in", device = "svg")

# 4.3. Not accounting for population structure 
# Running a simple RDA model
RDA_env_unconstrained <- rda(Geno_imp ~ bio2 + bio5 + bio7 + bio8 + bio9 + bio14 + bio15 + bio18,  Variables)
#RDA_env_unconstrained <- rda(Geno_imp ~ bio2 + bio7 + bio8 + bio9 + bio15,  Variables)

# Running the rdadapt function
rdadapt_env_unconstrained <- rdadapt(RDA_env_unconstrained, 2)

# Setting the p-value threshold 
thres_env <- 0.05

# Identifying the outliers for the simple RDA
outliers_unconstrained <- data.frame(Loci = colnames(Geno_imp)[which(rdadapt_env_unconstrained$q.values<thres_env)], q.value = rdadapt_env_unconstrained$q.values[which(rdadapt_env_unconstrained$q.values<thres_env)], contig = unlist(lapply(strsplit(colnames(Geno_imp)[which(rdadapt_env_unconstrained$q.values<thres_env)], split = "_"), function(x) x[1])))
outliers_unconstrained <- outliers_unconstrained[order(outliers_unconstrained$contig, outliers_unconstrained$q.value),]
outliers_rdadapt_env_unconstrained <- as.character(outliers_unconstrained$Loci[!duplicated(outliers_unconstrained$contig)])

# For all the outliers
list_outliers_RDA_all <- list(RDA_constrained = as.character(outliers$Loci), RDA_unconstrained = as.character(outliers_unconstrained$Loci))
ggVennDiagram(list_outliers_RDA_all, category.names = c("partial RDA", "simple RDA"), lty="solid", size=0.25) + 
  scale_fill_gradient(low = "white", high = "gray60") + scale_color_manual(values = rep("gray70", 4)) + guides(fill = "none") + theme_void(base_family = "Helvetica") +
  theme(text = element_text(size = 14, family = "Helvetica"), legend.position = "none", plot.margin = margin(10, 10, 10, 10))

# Only for the top hit locus per contig
list_outliers_RDA_top <- list(RDA_constrained = outliers_rdadapt_env, RDA_unconstrained = outliers_rdadapt_env_unconstrained)
ggVennDiagram(list_outliers_RDA_top, category.names = c("partial RDA", "simple RDA"), lty="solid", size=0.25) + 
  scale_fill_gradient(low = "white", high = "gray60") + scale_color_manual(values = rep("gray70", 4)) + guides(fill = "none") + theme_void(base_family = "Helvetica") +
  theme(text = element_text(size = 14, family = "Helvetica"), legend.position = "none", plot.margin = margin(10, 10, 10, 10))

#common_outliers_RDA_top <- Reduce(intersect, list_outliers_RDA_top)


# 5. Adaptive landscape: projecting adaptive gradient(s) across space 
# 5.1. Adaptively enriched genetic space 
# Adaptively enriched RDA
outliers_RDA_top <- outliers[outliers$q.value < 0.01,] # n = 550 -> n = 450
#RDA_outliers <- rda(Geno_imp[,outliers_RDA_top$Loci] ~ bio2 + bio5 + bio7 + bio8 + bio9 + bio14 + bio15 + bio18,  Variables)
RDA_outliers <- rda(Geno_imp[,outliers$Loci] ~ bio2 + bio5 + bio7 + bio8 + bio9 + bio14 + bio15 + bio18 + aspect + slope, Variables)
RDA_enriched <- rda(Geno_imp[,outliers_RDA_top$Loci] ~ bio2 + bio5 + bio7 + bio8 + bio9 + bio14 + bio15 + bio18 + aspect + slope, Variables)
summary(eigenvals(RDA_outliers, model = "constrained"))
#Importance of components:
#                         RDA1   RDA2    RDA3    RDA4    RDA5    RDA6    RDA7
#Eigenvalue            27.2641 6.0385 3.68558 3.05380 1.27790 1.18060 1.03124
#Proportion Explained   0.5932 0.1314 0.08019 0.06645 0.02781 0.02569 0.02244
#Cumulative Proportion  0.5932 0.7246 0.80481 0.87125 0.89906 0.92475 0.94718
#                         RDA8    RDA9   RDA10
#Eigenvalue            0.93162 0.80851 0.68724
#Proportion Explained  0.02027 0.01759 0.01495
#Cumulative Proportion 0.96745 0.98505 1.00000

summary(eigenvals(RDA_enriched, model = "constrained"))
#Importance of components:
#                         RDA1   RDA2    RDA3    RDA4    RDA5    RDA6    RDA7
#Eigenvalue            23.9527 5.2782 3.17223 2.55803 1.15760 1.03535 0.93659
#Proportion Explained   0.5957 0.1313 0.07889 0.06361 0.02879 0.02575 0.02329
#Cumulative Proportion  0.5957 0.7269 0.80582 0.86944 0.89823 0.92397 0.94726
#                         RDA8    RDA9   RDA10
#Eigenvalue            0.82536 0.69738 0.59781
#Proportion Explained  0.02053 0.01734 0.01487
#Cumulative Proportion 0.96779 0.98513 1.00000

# RDA biplot
TAB_loci <- as.data.frame(scores(RDA_enriched, choices=c(1:2), display="species", scaling="none"))
TAB_var <- as.data.frame(scores(RDA_enriched, choices=c(1:2), display="bp"))
p <- ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = "gray80", linewidth=0.5) +
  geom_vline(xintercept=0, linetype="dashed", color = "gray80", linewidth=0.5) +
  geom_point(data = TAB_loci, aes(x=RDA1*20, y=RDA2*20), colour = "#EB8055FF", size = 2, alpha = 0.75) + #"#F9A242FF"
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", linewidth=0.3, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 4, family = "Helvetica") +
  xlab("RDA 1 (24%)") + ylab("RDA 2 (5%)") +
  facet_wrap(~"Adaptively enriched RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(panel.grid = element_blank(), axis.ticks = element_blank(), legend.key = element_blank(), panel.background = element_blank(), legend.background = element_blank(), legend.title = element_text(size = 10), legend.text = element_text(size = 9), 
  plot.background = element_blank(), strip.text = element_text(size=11, face = "bold"), plot.margin = margin(10, 10, 10, 10), aspect.ratio = 1)
ggsave("rda_outlier_plot.svg", plot = p, width = 6, height = 6, units = "in", device = "svg")
ggsave("rda_outlier_plot.pdf", plot = p, width = 6, height = 6, units = "in", device = cairo_pdf)

# 5.2. Adaptive index across the landscape  
# The adaptive index provides an estimate of adaptive genetic similarity or difference of all pixels on the landscape as a function of the values of the environmental predictors at that location. 
# When projected on a map it allows visualizing the different adaptive gradients across a species range.

# Function to predict the adaptive index across the landscape
source("adaptive_index.R")

# Running the function for all the climatic pixels of distribution range
#drop_idx <- grep("aspect|slope|bio5|bio14|bio18", names(ras_present), ignore.case = TRUE)
#ras_present_bioclim <- ras_present[[ -drop_idx ]]
ras_present_bioclim <- ras_present
res_RDAoutliers_proj_current <- adaptive_index(RDA = RDA_outliers, K = 2, env = ras_present_bioclim, method = "loadings", scale_env = scale_env, center_env = center_env) #, env_mask = range, scale_env = scale_env, center_env = center_env
res_RDAenriched_proj_current <- adaptive_index(RDA = RDA_enriched, K = 2, env = ras_present_bioclim, method = "loadings", scale_env = scale_env, center_env = center_env) #, env_mask = range, scale_env = scale_env, center_env = center_env

# Vectorization of the climatic rasters for ggplot
RDAoutliers_proj <- list(res_RDAoutliers_proj_current$RDA1, res_RDAoutliers_proj_current$RDA2)
names(RDAoutliers_proj[[1]]) <- "val_RDA1"
names(RDAoutliers_proj[[2]]) <- "val_RDA2"
RDAoutliers_pts <- lapply(RDAoutliers_proj, function(x) rasterToPoints(x))
for(i in 1:length(RDAoutliers_pts)){
  RDAoutliers_pts[[i]][,3] <- (RDAoutliers_pts[[i]][,3]-min(RDAoutliers_pts[[i]][,3]))/(max(RDAoutliers_pts[[i]][,3])-min(RDAoutliers_pts[[i]][,3]))
}

RDAenriched_proj <- list(res_RDAenriched_proj_current$RDA1, res_RDAenriched_proj_current$RDA2)
names(RDAenriched_proj[[1]]) <- "val_RDA1"
names(RDAenriched_proj[[2]]) <- "val_RDA2"
RDAenriched_pts <- lapply(RDAenriched_proj, function(x) rasterToPoints(x))
for(i in 1:length(RDAenriched_pts)){
  RDAenriched_pts[[i]][,3] <- (RDAenriched_pts[[i]][,3]-min(RDAenriched_pts[[i]][,3]))/(max(RDAenriched_pts[[i]][,3])-min(RDAenriched_pts[[i]][,3]))
}

# Adaptive genetic turnover projected across species range for RDA1 and RDA2 indexes
TAB_RDAoutliers <- as.data.frame(do.call(rbind, RDAoutliers_pts[1:2]))
colnames(TAB_RDAoutliers)[3] <- "value"
TAB_RDAoutliers$variable <- factor(c(rep("RDA1", nrow(RDAoutliers_pts[[1]])), rep("RDA2", nrow(RDAoutliers_pts[[2]]))), levels = c("RDA1","RDA2"))

TAB_RDAenriched <- as.data.frame(do.call(rbind, RDAenriched_pts[1:2]))
colnames(TAB_RDAenriched)[3] <- "value"
TAB_RDAenriched$variable <- factor(c(rep("RDA1", nrow(RDAenriched_pts[[1]])), rep("RDA2", nrow(RDAenriched_pts[[2]]))), levels = c("RDA1","RDA2"))

bbox <- st_bbox(range)
p <- ggplot(data = TAB_RDAoutliers) + 
  geom_sf(data = admin, fill=gray(0.95), color = NA) +
  geom_raster(aes(x = x, y = y, fill = cut(value, breaks=seq(0, 1, length.out=10), include.lowest = T))) + 
  scale_fill_viridis_d(alpha = 0.85, direction = -1, option = "A", labels = c("Negative scores", rep("", 3), "Intermediate scores", rep("", 3), "Positive scores")) +
  geom_sf(data = admin, fill=NA, color = "black", size = 0.2) +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"]), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude") +
  guides(fill=guide_legend(title="Adaptive index")) +
  facet_grid(~ variable) +
  theme_minimal(base_size = 12, base_family = "Helvetica")  +
  theme(panel.grid = element_blank(), axis.ticks = element_blank(), legend.title = element_text(size = 10), legend.text = element_text(size = 9), 
  plot.margin = margin(10, 10, 10, 10), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size = 11, face = "bold"))
ggsave("adaptive_index_outliers.svg", plot = p, width = 10, height = 10, units = "in", device = "svg")
ggsave("adaptive_index_outliers.pdf", plot = p, width = 10, height = 10, units = "in", device = cairo_pdf)

p <- ggplot(data = TAB_RDAenriched) + 
  geom_sf(data = admin, fill=gray(0.95), color = NA) +
  geom_raster(aes(x = x, y = y, fill = cut(value, breaks=seq(0, 1, length.out=10), include.lowest = T))) + 
  scale_fill_viridis_d(alpha = 0.85, direction = -1, option = "A", labels = c("Negative scores", rep("", 3), "Intermediate scores", rep("", 3), "Positive scores")) +
  geom_sf(data = admin, fill=NA, color = "black", size = 0.2) +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"]), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude") +
  guides(fill=guide_legend(title="Adaptive index")) +
  facet_grid(~ variable) +
  theme_minimal(base_size = 12, base_family = "Helvetica")  +
  theme(panel.grid = element_blank(), axis.ticks = element_blank(), legend.title = element_text(size = 10), legend.text = element_text(size = 9), 
  plot.margin = margin(10, 10, 10, 10), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size = 11, face = "bold"))
ggsave("adaptive_index_enriched.svg", plot = p, width = 10, height = 10, units = "in", device = "svg")
ggsave("adaptive_index_enriched.pdf", plot = p, width = 10, height = 10, units = "in", device = cairo_pdf)

# which one to save between RDAoutliers vs RDAenchriched - landscape patterns are almost the same, so go with the normal RDAoutliers -> clip with 5km_1y layer in ArcGIS later
writeRaster(res_RDAoutliers_proj_current$RDA1, filename = "RDA1_present_final.tif", filetype = "GTiff", overwrite = TRUE)
writeRaster(res_RDAoutliers_proj_current$RDA2, filename = "RDA2_present_final.tif", filetype = "GTiff", overwrite = TRUE)


# 6. Predicting local maladaptation: genomic offset 
# 6.1. Predicting future genomic offset using RDA
# Once the genetic ~ environment relationship is characterized, it can be extrapolated to future environments to predict a potential shift in adaptive optimum induced by climate change
# RDA can be used to predict the optimal adaptive genetic composition for each environmental pixel under consideration (Adaptive Index in Section 5), using both current and future environmental conditions.
# The difference between the two predictions provides an estimate of the change in genetic composition that would be required to track climate change.

# Function to predict genomic offset from a RDA model
source("genomic_offset.R") 
#genomic_offset <- function(RDA, K, env_pres, env_fut, env_mask = NULL, method = "loadings")
#{
#  ## CHECKS -------------------------------------------------------------------
#  if (!inherits(RDA, "rda")) { stop("\n RDA must be a 'rda' object") }
#  if (!("CCA" %in% names(RDA)) || !("eig" %in% names(RDA$CCA))) {
#    stop("\n RDA$CCA$eig seems not to exist")
#  }
#  if (K > length(RDA$CCA$eig)) {
#    K <- length(RDA$CCA$eig)
#    warning("\n Not enough RDA axis available, K is set to length(RDA$CCA$eig)")
#  }

#  ## FUNCTION -----------------------------------------------------------------

#  ## Make predictions
#  AI_pres <- adaptive_index(RDA = RDA
#                            , K = K
#                            , env = env_pres
#                            , env_mask = env_mask
#                            , method = method)
#  AI_fut <- adaptive_index(RDA = RDA
#                           , K = K
#                           , env = env_fut
#                           , env_mask = env_mask
#                           , method = method)

#  ## Single axis genetic offset -----------------------------------------------
#  offset <- foreach(i = 1:K) %do%
#    {
#      ras <- abs(AI_pres[[i]] - AI_fut[[i]])
#      names(ras) <- paste0("RDA", i)
#      return(ras)
#    }
#  offset <- rast(offset)

#  ## Weight current and future adaptive indices based on eigen values of associated axes
#  weights <- RDA$CCA$eig / sum(RDA$CCA$eig)
#  weights <- weights[1:K]
#  AI_pres_w <- AI_pres * weights
#  AI_fut_w <- AI_fut * weights

#  ## Predict a global genetic offset
#  diff_stack <- AI_pres_w - AI_fut_w
#  squared_diff <- diff_stack ^ 2
#  offset_global <- sqrt(sum(squared_diff))
#  #offset_global <- offset[[1]]
#  #offset_global[!is.na(offset_global)] <- sapply(1:ncell(AI_pres_w), function(x) {
#  #  tmp_mat <- rbind(AI_pres_w[x], AI_fut_w[x])
#  #  return(dist(tmp_mat, method = "euclidean"))
#  #})
#  names(offset_global) <- "Global_offset"

#  return(list(Proj_pres = AI_pres,
#              Proj_fut = AI_fut,
#              Proj_offset = offset,
#              Proj_offset_global = offset_global,
#              weights = weights))
#}

# Running the function for future using outliers
#drop_idx <- grep("aspect|slope|bio5|bio14|bio18", names(ras_present), ignore.case = TRUE)
#ras_present_bioclim <- ras_present[[ -drop_idx ]]
#ras_ssp126_2011_2040_bioclim <- ras_ssp126_2011_2040[[ -drop_idx ]]
#ras_ssp126_2041_2070_bioclim <- ras_ssp126_2041_2070[[ -drop_idx ]]
#ras_ssp126_2071_2100_bioclim <- ras_ssp126_2071_2100[[ -drop_idx ]]
#ras_ssp585_2011_2040_bioclim <- ras_ssp585_2011_2040[[ -drop_idx ]]
#ras_ssp585_2041_2070_bioclim <- ras_ssp585_2041_2070[[ -drop_idx ]]
#ras_ssp585_2071_2100_bioclim <- ras_ssp585_2071_2100[[ -drop_idx ]]

#res_RDA_proj2040 <- genomic_offset(RDA_env, K = 2, env_pres = rast_present_bioclim, env_fut = rast_2011_2040_bioclim, method = "loadings") #, range = range,  scale_env = scale_env, center_env = center_env)
#res_RDA_proj2070 <- genomic_offset(RDA_env, K = 2, env_pres = rast_present_bioclim, env_fut = rast_2041_2070_bioclim, method = "loadings") #, range = range,  scale_env = scale_env, center_env = center_env)
#res_RDA_proj2100 <- genomic_offset(RDA_env, K = 2, env_pres = rast_present_bioclim, env_fut = rast_2071_2100_bioclim, method = "loadings") #, range = range,  scale_env = scale_env, center_env = center_env)

res_RDA_ssp126_proj2040 <- genomic_offset(RDA_outliers, K = 2, env_pres = ras_present, env_fut = ras_ssp126_2011_2040, method = "loadings", range = range, scale_env = scale_env, center_env = center_env)
res_RDA_ssp126_proj2070 <- genomic_offset(RDA_outliers, K = 2, env_pres = ras_present, env_fut = ras_ssp126_2041_2070, method = "loadings", range = range, scale_env = scale_env, center_env = center_env)
res_RDA_ssp126_proj2100 <- genomic_offset(RDA_outliers, K = 2, env_pres = ras_present, env_fut = ras_ssp126_2071_2100, method = "loadings", range = range, scale_env = scale_env, center_env = center_env)
res_RDA_ssp585_proj2040 <- genomic_offset(RDA_outliers, K = 2, env_pres = ras_present, env_fut = ras_ssp585_2011_2040, method = "loadings", range = range, scale_env = scale_env, center_env = center_env)
res_RDA_ssp585_proj2070 <- genomic_offset(RDA_outliers, K = 2, env_pres = ras_present, env_fut = ras_ssp585_2041_2070, method = "loadings", range = range, scale_env = scale_env, center_env = center_env)
res_RDA_ssp585_proj2100 <- genomic_offset(RDA_outliers, K = 2, env_pres = ras_present, env_fut = ras_ssp585_2071_2100, method = "loadings", range = range, scale_env = scale_env, center_env = center_env)

# Table global genetic offset predicted for 2050 and 2080
RDA_ssp126_proj_offset <- data.frame(rbind(rasterToPoints(res_RDA_ssp126_proj2040$Proj_offset_global), rasterToPoints(res_RDA_ssp126_proj2070$Proj_offset_global), rasterToPoints(res_RDA_ssp126_proj2100$Proj_offset_global)), Date = c(rep("2040", nrow(rasterToPoints(res_RDA_ssp126_proj2040$Proj_offset_global))), rep("2070", nrow(rasterToPoints(res_RDA_ssp126_proj2070$Proj_offset_global))), rep("2100", nrow(rasterToPoints(res_RDA_ssp126_proj2100$Proj_offset_global)))))
RDA_ssp585_proj_offset <- data.frame(rbind(rasterToPoints(res_RDA_ssp585_proj2040$Proj_offset_global), rasterToPoints(res_RDA_ssp585_proj2070$Proj_offset_global), rasterToPoints(res_RDA_ssp585_proj2100$Proj_offset_global)), Date = c(rep("2040", nrow(rasterToPoints(res_RDA_ssp585_proj2040$Proj_offset_global))), rep("2070", nrow(rasterToPoints(res_RDA_ssp585_proj2070$Proj_offset_global))), rep("2100", nrow(rasterToPoints(res_RDA_ssp585_proj2100$Proj_offset_global)))))

# Projecting genomic offset on a map
range_vals <- range(RDA_ssp126_proj_offset$Global_offset, na.rm = TRUE)
#breaks <- seq(range_vals[1], range_vals[2], length.out = 9)  # 8 bins = 9 breakpoints
quantile_breaks <- quantile(RDA_ssp126_proj_offset$Global_offset, probs = seq(0, 1, length.out = 9), na.rm = TRUE)
#bin_labels <- paste0(format(round(breaks[-length(breaks)], 2), nsmall = 2), "–", format(round(breaks[-1], 2), nsmall = 2))
bin_labels <- paste0(format(round(quantile_breaks[-length(quantile_breaks)], 2), nsmall = 2), "–", format(round(quantile_breaks[-1], 2), nsmall = 2))
#fill = cut(RDA_proj_offset$Global_offset, breaks = breaks, include.lowest = TRUE)
fill = cut(RDA_ssp126_proj_offset$Global_offset, breaks = quantile_breaks, include.lowest = TRUE, labels = bin_labels)
colors <- colorRampPalette(brewer.pal(11, "Spectral")[c(6:1)])(8)
p <- ggplot(data = RDA_ssp126_proj_offset) + 
  geom_sf(data = admin, fill=gray(.9), size=0) +
  geom_raster(aes(x = x, y = y, fill = fill), alpha = 1) + 
  scale_fill_manual(values = colors, labels = bin_labels, guide = guide_legend(title="Genomic offset", title.position = "top", title.hjust = 0.5, ncol = 1, label.position="right"), na.translate = F) +
  geom_sf(data = admin, fill=NA, color = "black", size = 0.2) +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"])) + #, expand = FALSE
  xlab("Longitude") + ylab("Latitude") +
  facet_grid(~ Date) +
  theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11, face = "bold"))
ggsave("RDA_ssp126_genomic_offset.svg", plot = p, width = 15, height = 10, units = "in", device = "svg")
ggsave("RDA_ssp126_genomic_offset.pdf", plot = p, width = 15, height = 10, units = "in", device = cairo_pdf)

range_vals <- range(RDA_ssp585_proj_offset$Global_offset, na.rm = TRUE)
#breaks <- seq(range_vals[1], range_vals[2], length.out = 9)  # 8 bins = 9 breakpoints
quantile_breaks <- quantile(RDA_ssp585_proj_offset$Global_offset, probs = seq(0, 1, length.out = 9), na.rm = TRUE)
#bin_labels <- paste0(format(round(breaks[-length(breaks)], 2), nsmall = 2), "–", format(round(breaks[-1], 2), nsmall = 2))
bin_labels <- paste0(format(round(quantile_breaks[-length(quantile_breaks)], 2), nsmall = 2), "–", format(round(quantile_breaks[-1], 2), nsmall = 2))
#fill = cut(RDAenriched_proj_offset$Global_offset, breaks = breaks, include.lowest = TRUE)
fill = cut(RDA_ssp585_proj_offset$Global_offset, breaks = quantile_breaks, include.lowest = TRUE, labels = bin_labels)
colors <- colorRampPalette(brewer.pal(11, "Spectral")[c(6:1)])(8)
p <- ggplot(data = RDA_ssp585_proj_offset) + 
  geom_sf(data = admin, fill=gray(.9), size=0) +
  geom_raster(aes(x = x, y = y, fill = fill), alpha = 1) + 
  scale_fill_manual(values = colors, labels = bin_labels, guide = guide_legend(title="Genomic offset", title.position = "top", title.hjust = 0.5, ncol = 1, label.position="right"), na.translate = F) +
  geom_sf(data = admin, fill=NA, color = "black", size = 0.2) +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"])) + #, expand = FALSE
  xlab("Longitude") + ylab("Latitude") +
  facet_grid(~ Date) +
  theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11, face = "bold"))
ggsave("RDA_ssp585_genomic_offset.svg", plot = p, width = 15, height = 10, units = "in", device = "svg")
ggsave("RDA_ssp585_genomic_offset.pdf", plot = p, width = 15, height = 10, units = "in", device = cairo_pdf)

writeRaster(res_RDA_ssp126_proj2040$Proj_offset_global, filename = "RDA_ssp126_offset_2040_final.tif", filetype = "GTiff", overwrite = TRUE)
writeRaster(res_RDA_ssp126_proj2070$Proj_offset_global, filename = "RDA_ssp126_offset_2070_final.tif", filetype = "GTiff", overwrite = TRUE)
writeRaster(res_RDA_ssp126_proj2100$Proj_offset_global, filename = "RDA_ssp126_offset_2100_final.tif", filetype = "GTiff", overwrite = TRUE)
writeRaster(res_RDA_ssp585_proj2040$Proj_offset_global, filename = "RDA_ssp585_offset_2040_final.tif", filetype = "GTiff", overwrite = TRUE)
writeRaster(res_RDA_ssp585_proj2070$Proj_offset_global, filename = "RDA_ssp585_offset_2070_final.tif", filetype = "GTiff", overwrite = TRUE)
writeRaster(res_RDA_ssp585_proj2100$Proj_offset_global, filename = "RDA_ssp585_offset_2100_final.tif", filetype = "GTiff", overwrite = TRUE)


# Postprocessing
# Identify IDs of outliers in the original vcf
vcf <- read.table("monq_mcSV.final.sort.vcf", comment.char = "#", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(vcf)[1:5] <- c("CHROM", "POS", "ID", "REF", "ALT")
vcf$base_id <- paste(vcf$CHROM, vcf$POS, sep = "_")
vcf$allele_index <- ave(vcf$base_id, vcf$base_id, FUN = seq_along)
vcf$variant_id <- paste0(vcf$base_id, "_", vcf$allele_index)

matched <- vcf[vcf$variant_id %in% outliers$Loci, ]
write.table(matched, "GEAoutliers_full_lines_final.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


# Final plotting -> come back after ArcGIS if I need to clip
# Vectorization of the climatic rasters for ggplot
RDAoutliers_proj_clip <- list(rast("RDA1_present_clip.tif"), rast("RDA2_present_clip.tif"))
RDAoutliers_proj_clip <- lapply(RDAoutliers_proj_clip, function(x) rasterToPoints(raster(x)))
for(i in 1:length(RDAoutliers_proj_clip)){
  RDAoutliers_proj_clip[[i]][,3] <- (RDAoutliers_proj_clip[[i]][,3]-min(RDAoutliers_proj_clip[[i]][,3]))/(max(RDAoutliers_proj_clip[[i]][,3])-min(RDAoutliers_proj_clip[[i]][,3]))
}

# Adaptive genetic turnover projected across species range for RDA1 and RDA2 indexes
TAB_RDAoutliers_clip <- as.data.frame(do.call(rbind, RDAoutliers_proj_clip[1:2]))
colnames(TAB_RDAoutliers_clip)[3] <- "value"
TAB_RDAoutliers_clip$variable <- factor(c(rep("RDA1", nrow(RDAoutliers_proj_clip[[1]])), rep("RDA2", nrow(RDAoutliers_proj_clip[[2]]))), levels = c("RDA1","RDA2"))

p <- ggplot(data = TAB_RDAoutliers_clip) + 
  geom_sf(data = admin, fill=gray(0.95), color = NA) +
  geom_raster(aes(x = x, y = y, fill = cut(value, breaks=seq(0, 1, length.out=10), include.lowest = T))) + 
  scale_fill_viridis_d(alpha = 0.85, direction = -1, option = "A", labels = c("Negative scores", rep("", 3), "Intermediate scores", rep("", 3), "Positive scores")) +
  geom_sf(data = admin, fill=NA, color = "black", size = 0.2) +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"]), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude") +
  guides(fill=guide_legend(title="Adaptive index")) +
  facet_grid(~ variable) +
  theme_minimal(base_size = 12, base_family = "Helvetica")  +
  theme(panel.grid = element_blank(), axis.ticks = element_blank(), legend.title = element_text(size = 10), legend.text = element_text(size = 9), 
  plot.margin = margin(10, 10, 10, 10), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size = 11, face = "bold"))
ggsave("adaptive_index_outliers_clip.svg", plot = p, width = 10, height = 10, units = "in", device = "svg")
ggsave("adaptive_index_outliers_clip.pdf", plot = p, width = 10, height = 10, units = "in", device = cairo_pdf)


# Table global genetic offset predicted for 2050 and 2080
RDA_ssp126_proj_offset_clip <- data.frame(rbind(rasterToPoints(raster("RDA_ssp126_offset_2040_clip.tif")), rasterToPoints(raster("RDA_ssp126_offset_2070_clip.tif")), rasterToPoints(raster("RDA_ssp126_offset_2100.tif"))), Date = c(rep("2040", nrow(rasterToPoints(raster("RDA_ssp126_offset_2040_clip.tif")))), rep("2070", nrow(rasterToPoints(raster("RDA_ssp126_offset_2070_clip.tif")))), rep("2100", nrow(rasterToPoints(raster("RDA_ssp126_offset_2100.tif"))))))
RDA_ssp585_proj_offset_clip <- data.frame(rbind(rasterToPoints(raster("RDA_ssp585_offset_2040_clip.tif")), rasterToPoints(raster("RDA_ssp585_offset_2070_clip.tif")), rasterToPoints(raster("RDA_ssp585_offset_2100.tif"))), Date = c(rep("2040", nrow(rasterToPoints(raster("RDA_ssp585_offset_2040_clip.tif")))), rep("2070", nrow(rasterToPoints(raster("RDA_ssp585_offset_2070_clip.tif")))), rep("2100", nrow(rasterToPoints(raster("RDA_ssp585_offset_2100.tif"))))))

# Projecting genomic offset on a map
range_vals <- range(RDA_ssp126_proj_offset_clip$RDA_ssp126_offset_2040_clip, na.rm = TRUE)
quantile_breaks <- quantile(RDA_ssp126_proj_offset_clip$RDA_ssp126_offset_2040_clip, probs = seq(0, 1, length.out = 9), na.rm = TRUE)
bin_labels <- paste0(format(round(quantile_breaks[-length(quantile_breaks)], 2), nsmall = 2), "–", format(round(quantile_breaks[-1], 2), nsmall = 2))
fill = cut(RDA_ssp126_proj_offset_clip$RDA_ssp126_offset_2040_clip, breaks = quantile_breaks, include.lowest = TRUE, labels = bin_labels)
colors <- colorRampPalette(brewer.pal(11, "Spectral")[c(6:1)])(8)
p <- ggplot(data = RDA_ssp126_proj_offset_clip) + 
  geom_sf(data = admin, fill=gray(.9), size=0) +
  geom_raster(aes(x = x, y = y, fill = fill), alpha = 1) + 
  scale_fill_manual(values = colors, labels = bin_labels, guide = guide_legend(title="Genomic offset", title.position = "top", title.hjust = 0.5, ncol = 1, label.position="right"), na.translate = F) +
  geom_sf(data = admin, fill=NA, color = "black", size = 0.2) +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"])) + #, expand = FALSE
  xlab("Longitude") + ylab("Latitude") +
  facet_grid(~ Date) +
  theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11, face = "bold"))
ggsave("RDA_ssp126_genomic_offset_clip.svg", plot = p, width = 15, height = 10, units = "in", device = "svg")
ggsave("RDA_ssp126_genomic_offset_clip.pdf", plot = p, width = 15, height = 10, units = "in", device = cairo_pdf)

range_vals <- range(RDA_ssp585_proj_offset_clip$RDA_ssp585_offset_2040_clip, na.rm = TRUE)
quantile_breaks <- quantile(RDA_ssp585_proj_offset_clip$RDA_ssp585_offset_2040_clip, probs = seq(0, 1, length.out = 9), na.rm = TRUE)
bin_labels <- paste0(format(round(quantile_breaks[-length(quantile_breaks)], 2), nsmall = 2), "–", format(round(quantile_breaks[-1], 2), nsmall = 2))
fill = cut(RDA_ssp585_proj_offset_clip$RDA_ssp585_offset_2040_clip, breaks = quantile_breaks, include.lowest = TRUE, labels = bin_labels)
colors <- colorRampPalette(brewer.pal(11, "Spectral")[c(6:1)])(8)
p <- ggplot(data = RDA_ssp585_proj_offset_clip) + 
  geom_sf(data = admin, fill=gray(.9), size=0) +
  geom_raster(aes(x = x, y = y, fill = fill), alpha = 1) + 
  scale_fill_manual(values = colors, labels = bin_labels, guide = guide_legend(title="Genomic offset", title.position = "top", title.hjust = 0.5, ncol = 1, label.position="right"), na.translate = F) +
  geom_sf(data = admin, fill=NA, color = "black", size = 0.2) +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"])) + #, expand = FALSE
  xlab("Longitude") + ylab("Latitude") +
  facet_grid(~ Date) +
  theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11, face = "bold"))
ggsave("RDA_ssp585_genomic_offset_clip.svg", plot = p, width = 15, height = 10, units = "in", device = "svg")
ggsave("RDA_ssp585_genomic_offset_clip.pdf", plot = p, width = 15, height = 10, units = "in", device = cairo_pdf)


# 6.2 Predicting future genomic offset using Generalized Dissimilarity Modeling (modifed from Gougherty et al. (2021) version)
# Settings ----------------------------------------------------------------
nCores<-50
nBreak<-100

# Loading required R packages ------------------------------------------------
Sys.setenv(R_LIBS_USER="/home/jeon96/R/negishi/4.2.2-gcc-12.2.0-524zdta")
.libPaths(Sys.getenv("R_LIBS_USER"))

library(adehabitatLT)
library(fields)
library(reshape)
library(parallel)
library(doParallel)
library(geosphere)
library(gdm)
library(foreach)
library(raster)
library(tidyverse)
library(terra)
library(sf)
library(missMDA)

# Function to remove intercept from gdm predictions
removeIntercept <- function(mod,pred){
  adjust <- 0 - log(1-pred) - mod$intercept
  adjustDissim <- 1-exp(0-adjust)
  return(adjustDissim)
}

rast.aggr_present <- aggregate(rast(ras_present), fact = 10, fun = 'mean', na.rm = TRUE, cores = nCores)

rast.aggr_ssp126_2011_2040 <- aggregate(rast(ras_ssp126_2011_2040), fact = 10, fun = 'mean', na.rm = TRUE, cores = nCores)
rast.aggr_ssp126_2041_2070 <- aggregate(rast(ras_ssp126_2041_2070), fact = 10, fun = 'mean', na.rm = TRUE, cores = nCores)
rast.aggr_ssp126_2071_2100 <- aggregate(rast(ras_ssp126_2071_2100), fact = 10, fun = 'mean', na.rm = TRUE, cores = nCores)
rast.aggr_ssp585_2011_2040 <- aggregate(rast(ras_ssp585_2011_2040), fact = 10, fun = 'mean', na.rm = TRUE, cores = nCores)
rast.aggr_ssp585_2041_2070 <- aggregate(rast(ras_ssp585_2041_2070), fact = 10, fun = 'mean', na.rm = TRUE, cores = nCores)
rast.aggr_ssp585_2071_2100 <- aggregate(rast(ras_ssp585_2071_2100), fact = 10, fun = 'mean', na.rm = TRUE, cores = nCores)

rast.aggr_scen.fut <- list(
  ssp126_2011_2040 = rast.aggr_ssp126_2011_2040,
  ssp126_2041_2070 = rast.aggr_ssp126_2041_2070,
  ssp126_2071_2100 = rast.aggr_ssp126_2071_2100,
  ssp585_2011_2040 = rast.aggr_ssp585_2011_2040,
  ssp585_2041_2070 = rast.aggr_ssp585_2041_2070,
  ssp585_2071_2100 = rast.aggr_ssp585_2071_2100
)

# Species range shapefile (download information at beginning of tutorial)
shp <- shapefile("monq_450km_90y_onLand.shp")

outliers <- readRDS("RDA_outliers.rds")

##############
#Load and prep PCA-based genetic distance data
##############
# Converting discrete genotype values to PCA-based continuous values of candidate SNPs 
Geno_outlier <- Geno_imp[,colnames(Geno_imp) %in% outliers$Loci]

Geno_outlier_pca <- prcomp(Geno_outlier, center = TRUE, scale. = TRUE)
eigvals <- Geno_outlier_pca$sdev^2
var_expl <- eigvals / sum(eigvals)
cumvar   <- cumsum(var_expl)
target <- 1.00 # change as needed
rank_pc <- length(eigvals)
if (all(cumvar < target)) {
  n_pc <- rank_pc
} else {
  n_pc <- min(which(cumvar >= target)[1], rank_pc)
}
n_pc # 87

pc_n <- Geno_outlier_pca$x[, seq_len(n_pc), drop = FALSE]
D <- as.matrix(dist(pc_n))                   # n x n
diag(D) <- 0
mx <- max(D[upper.tri(D)], na.rm = TRUE)
D01 <- D / mx

colnames(D01) <- rownames(Geno_outlier)
dist <- as.data.frame(D01)
dist <- cbind(inds = rownames(Geno_outlier), dist)

# Reading in population locations
inds <- InfoInd_filt[,1:3]
ids  <- as.character(inds[, 1])

# Confirming the two are in the same order
all(dist$inds %in% ids)
all(dist$inds == ids)

##############
# Load and prep shapefile and climate data
##############
# Loading shapefile; optional

# Choosing predictors
predNames <- c("bio2", "bio5", "bio7", "bio8", "bio9", "bio14", "bio15", "bio18")
presClims <- rast.aggr_present #ras_present
presClims <- presClims[[predNames]]
presClims_r <- stack(presClims) 

# Creating pred data for gdm (cols = population name, long, lat, climate data)
pred <- data.frame(site = inds$id, long = inds$x, lat = inds$y, extract(presClims_r, y = inds[,c("x","y")]),
                   stringsAsFactors = FALSE)

######################
#GDM model
######################
# Creating site pair table
dist$site <- as.integer(sub("^[^0-9]+", "", dist$inds)) # needs numeric site IDs
old_dist_cols <- setdiff(names(dist), c("inds","site"))
new_dist_cols <- sub("^[^0-9]+", "", old_dist_cols)
names(dist)[match(old_dist_cols, names(dist))] <- new_dist_cols
dist_cols <- as.character(dist$site)
dist <- cbind(site = dist$site, dist[, dist_cols, drop = FALSE])

sitePair <- formatsitepair(bioDat = dist, bioFormat = 3, siteColumn = "site", XColumn = "long", YColumn = "lat", predData = pred)

# Creating gdm
gdmMod <- gdm(na.omit(sitePair), geo = FALSE)
gdmMod_full <- gdm(na.omit(sitePair), geo = TRUE)

summary(gdmMod)
# GDM Modelling Summary
# Creation Date:  Mon Oct 27 20:51:24 2025
#
# Name:  gdmMod
#
# Data:  na.omit  Data:  sitePair
#
# Samples:  3828
#
# Geographical distance used in model fitting?  FALSE
#
# NULL Deviance:  479.147
# GDM Deviance:  183.899
# Percent Deviance Explained:  61.62
#
# Intercept:  0.427

summary(gdmMod_full)
#  GDM Modelling Summary
# Creation Date:  Mon Oct 27 20:51:24 2025
#
# Name:  gdmMod_full
#
# Data:  na.omit  Data:  sitePair
#
# Samples:  3828
#
# Geographical distance used in model fitting?  TRUE
#
# NULL Deviance:  479.147
# GDM Deviance:  183.899
# Percent Deviance Explained:  61.62
#
# Intercept:  0.427


# Loading future climate data
#futClims <- rast.aggr_ssp126_2011_2040 #ras_ssp126_2011_2040 # future climate layers
#futClimDat <- as.data.frame(futClims, xy=TRUE, na.rm=TRUE)

# Getting all coordinates in range map in current climate
#popDat <- na.omit(as.data.frame(mask(raster(presClims), shp), xy=TRUE))
#popDat <- data.frame(distance=1, weight=1, popDat)
#popDat <- split(popDat, seq(nrow(popDat)))

###############
#Reverse offset calculation
##############
#Getting all coordinates in the range in current climate
popDat <- na.omit(as.data.frame(mask(presClims_r, shp), xy = TRUE))
popDat <- data.frame(popDat)

for (scen_fut in names(rast.aggr_scen.fut)) {
  #Gets climate data from the range in future climate
  futClims <- rast.aggr_scen.fut[[scen_fut]] #ras_ssp126_2011_2040 # future climate layers
  futClims <- futClims[[predNames]]
  futClims_r <- stack(futClims)

  futClimDat <- na.omit(as.data.frame(mask(futClims_r, shp), xy = TRUE))

  # Setting up for prediction
  futClimDat <- data.frame(distance = 1, weight = 1, futClimDat)


  cl <- makeCluster(15)
  registerDoParallel(cl)
  parallel::clusterEvalQ(cl, {
    .libPaths("/home/jeon96/R/negishi/4.2.2-gcc-12.2.0-524zdta")
    library(gdm)
    library(fields)
    library(geosphere)
  })

  reverseOffsetGDM <- foreach(i = 1:nrow(futClimDat), .packages = c("fields","gdm","geosphere")) %dopar%{
  
    #get the focal population
    onePop <- futClimDat[i, ]
  
    #set up a dataframe where the first site is the focal population, and the second population are sites across the range
    setUp <- cbind(onePop,popDat)
    colnames(setUp) <- c("distance", "weights",
                         "s1.xCoord", "s1.yCoord", paste("s1.", predNames, sep=""), 
                         "s2.xCoord", "s2.yCoord", paste("s2.", predNames, sep=""))
  
    #rearrange the colums for the gdm prediction
    dat <- setUp[,c("distance","weights","s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord",
                    paste("s1.", predNames, sep = ""), 
                    paste("s2.", predNames, sep = ""))]
  
    #do the prediction and set up a dataframe with second sites x/y and predicted genetic distance
    combinedDat <- predict(object = gdmMod, dat, time = FALSE)
    combinedDat <- data.frame(dat[, c("s2.xCoord","s2.yCoord")], predGenDist = removeIntercept(gdmMod, combinedDat))
  
    ##Get metrics for the focal population
    #coordinate of focal population
    coord <- onePop[, c("x","y")]
  
    #choose the pixels with the minimum offset
    minCoords <- combinedDat[which(combinedDat$predGenDist == min(combinedDat$predGenDist)), ]
  
    #calculate the distance to the sites with minimum offset, and selct the one with the shortest distance
    minCoords["dists"] <- distGeo(p1 = coord, p2 = minCoords[, 1:2])
    minCoords <- minCoords[which(minCoords$dists == min(minCoords$dists)), ]

    if (nrow(minCoords) > 1) {
       cat(sprintf(
       "Ties at focal (%.6f, %.6f): %s\n",
       coord$x, coord$y,
       paste(sprintf("(%.6f, %.6f)", minCoords$x, minCoords$y), collapse = "; ")
       ))
    }
    
    #if multiple sites have the same offset, and same distance, one is randomly chosen
    minCoords <- minCoords[sample(1:nrow(minCoords), 1), ]
  
    #get local offset
    offset <- combinedDat[which(combinedDat$s2.xCoord == coord$x & combinedDat$s2.yCoord == coord$y), "predGenDist"]
  
    #get the minimum predicted offset - reverse offset in this case
    minVal <- minCoords$predGenDist
  
    #get distance and coordinates of site that minimizes offset
    toGo <- minCoords$dists
    minPt <- minCoords[, c("s2.xCoord", "s2.yCoord")]
  
    #get bearing to the site that minimizes offset
    bear <- bearing(coord, minPt)
  
    #write out
    out <- c(x1=coord[[1]], y1=coord[[2]], localOffset=offset, reverseOffset=minVal, predGeoDist=toGo, bearing=bear, x2=minPt[[1]], y2=minPt[[2]])
  
  }

  stopCluster(cl)

  #in this resultant dataframe the columns are:
  #x1/y1: focal coordinates
  #localOffset: local offset
  #reverseOffset: reverse offset in genetic distance (e.g., Fst or PCA distance)
  #predGeoDist: distance to site of reverse offset
  #bearing: bearing to site of reverse offset
  #x2/y2: coordinate of site of reverse offset
  
  reverseOffsetGDM <- do.call(rbind, reverseOffsetGDM)

  write.csv(reverseOffsetGDM,paste0("./reverseOffsetGDM_", scen_fut, ".csv"), row.names = FALSE)
}

GDM_ssp126_proj2040 <- read.csv("reverseOffsetGDM_ssp126_2011_2040.csv", header = T)
GDM_ssp126_proj2070 <- read.csv("reverseOffsetGDM_ssp126_2041_2070.csv", header = T)
GDM_ssp126_proj2100 <- read.csv("reverseOffsetGDM_ssp126_2071_2100.csv", header = T)

GDM_ssp126_proj_offset <- data.frame(rbind(GDM_ssp126_proj2040, GDM_ssp126_proj2070, GDM_ssp126_proj2100), Date = c(rep("2040", nrow(GDM_ssp126_proj2040)), rep("2070", nrow(GDM_ssp126_proj2070)), rep("2100", nrow(GDM_ssp126_proj2100))))

GDM_ssp585_proj2040 <- read.csv("reverseOffsetGDM_ssp585_2011_2040.csv", header = T)
GDM_ssp585_proj2070 <- read.csv("reverseOffsetGDM_ssp585_2041_2070.csv", header = T)
GDM_ssp585_proj2100 <- read.csv("reverseOffsetGDM_ssp585_2071_2100.csv", header = T)

GDM_ssp585_proj_offset <- data.frame(rbind(GDM_ssp585_proj2040, GDM_ssp585_proj2070, GDM_ssp585_proj2100), Date = c(rep("2040", nrow(GDM_ssp585_proj2040)), rep("2070", nrow(GDM_ssp585_proj2070)), rep("2100", nrow(GDM_ssp585_proj2100))))

# Projecting genomic offset on a map
range_vals <- range(GDM_ssp126_proj_offset$localOffset, na.rm = TRUE)
#breaks <- seq(range_vals[1], range_vals[2], length.out = 9)  # 8 bins = 9 breakpoints
quantile_breaks <- quantile(GDM_ssp126_proj_offset$localOffset, probs = seq(0, 1, length.out = 9), na.rm = TRUE)
#bin_labels <- paste0(format(round(breaks[-length(breaks)], 2), nsmall = 2), "–", format(round(breaks[-1], 2), nsmall = 2))
bin_labels <- paste0(format(round(quantile_breaks[-length(quantile_breaks)], 2), nsmall = 2), "–", format(round(quantile_breaks[-1], 2), nsmall = 2))
fill = cut(GDM_ssp126_proj_offset$localOffset, breaks = quantile_breaks, include.lowest = TRUE, labels = bin_labels)
colors <- colorRampPalette(brewer.pal(11, "Spectral")[c(6:1)])(8)
p <- ggplot(data = GDM_ssp126_proj_offset) + 
  geom_sf(data = admin, fill=gray(.9), size=0) +
  geom_raster(aes(x = x1, y = y1, fill = fill), alpha = 1) + 
  scale_fill_manual(values = colors, labels = bin_labels, guide = guide_legend(title="Local genomic offset", title.position = "top", title.hjust = 0.5, ncol = 1, label.position="right"), na.translate = F) +
  geom_sf(data = admin, fill=NA, color = "black", size = 0.2) +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"])) + #, expand = FALSE
  xlab("Longitude") + ylab("Latitude") +
  facet_grid(~ Date) +
  theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11, face = "bold"))
ggsave("GDM_ssp126_genomic_local_offset.svg", plot = p, width = 15, height = 10, units = "in", device = "svg")
ggsave("GDM_ssp126_genomic_local_offset.pdf", plot = p, width = 15, height = 10, units = "in", device = cairo_pdf)

range_vals <- range(GDM_ssp126_proj_offset$reverseOffset, na.rm = TRUE)
#breaks <- seq(range_vals[1], range_vals[2], length.out = 9)  # 8 bins = 9 breakpoints
quantile_breaks <- quantile(GDM_ssp126_proj_offset$reverseOffset, probs = seq(0, 1, length.out = 9), na.rm = TRUE)
#bin_labels <- paste0(format(round(breaks[-length(breaks)], 2), nsmall = 2), "–", format(round(breaks[-1], 2), nsmall = 2))
bin_labels <- paste0(format(round(quantile_breaks[-length(quantile_breaks)], 5), nsmall = 2), "–", format(round(quantile_breaks[-1], 5), nsmall = 2))
fill = cut(GDM_ssp126_proj_offset$reverseOffset, breaks = quantile_breaks, include.lowest = TRUE, labels = bin_labels)
colors <- colorRampPalette(brewer.pal(11, "Spectral")[c(6:1)])(8)
p <- ggplot(data = GDM_ssp126_proj_offset) + 
  geom_sf(data = admin, fill=gray(.9), size=0) +
  geom_raster(aes(x = x1, y = y1, fill = fill), alpha = 1) + 
  scale_fill_manual(values = colors, labels = bin_labels, guide = guide_legend(title="Reverse genomic offset", title.position = "top", title.hjust = 0.5, ncol = 1, label.position="right"), na.translate = F) +
  geom_sf(data = admin, fill=NA, color = "black", size = 0.2) +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"])) + #, expand = FALSE
  xlab("Longitude") + ylab("Latitude") +
  facet_grid(~ Date) +
  theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11, face = "bold"))
ggsave("GDM_ssp126_genomic_reverse_offset.svg", plot = p, width = 15, height = 10, units = "in", device = "svg")
ggsave("GDM_ssp126_genomic_reverse_offset.pdf", plot = p, width = 15, height = 10, units = "in", device = cairo_pdf)

range_vals <- range(GDM_ssp585_proj_offset$localOffset, na.rm = TRUE)
#breaks <- seq(range_vals[1], range_vals[2], length.out = 9)  # 8 bins = 9 breakpoints
quantile_breaks <- quantile(GDM_ssp585_proj_offset$localOffset, probs = seq(0, 1, length.out = 9), na.rm = TRUE)
#bin_labels <- paste0(format(round(breaks[-length(breaks)], 2), nsmall = 2), "–", format(round(breaks[-1], 2), nsmall = 2))
bin_labels <- paste0(format(round(quantile_breaks[-length(quantile_breaks)], 2), nsmall = 2), "–", format(round(quantile_breaks[-1], 2), nsmall = 2))
fill = cut(GDM_ssp585_proj_offset$localOffset, breaks = quantile_breaks, include.lowest = TRUE, labels = bin_labels)
colors <- colorRampPalette(brewer.pal(11, "Spectral")[c(6:1)])(8)
p <- ggplot(data = GDM_ssp585_proj_offset) + 
  geom_sf(data = admin, fill=gray(.9), size=0) +
  geom_raster(aes(x = x1, y = y1, fill = fill), alpha = 1) + 
  scale_fill_manual(values = colors, labels = bin_labels, guide = guide_legend(title="Local genomic offset", title.position = "top", title.hjust = 0.5, ncol = 1, label.position="right"), na.translate = F) +
  geom_sf(data = admin, fill=NA, color = "black", size = 0.2) +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"])) + #, expand = FALSE
  xlab("Longitude") + ylab("Latitude") +
  facet_grid(~ Date) +
  theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11, face = "bold"))
ggsave("GDM_ssp585_genomic_local_offset.svg", plot = p, width = 15, height = 10, units = "in", device = "svg")
ggsave("GDM_ssp585_genomic_local_offset.pdf", plot = p, width = 15, height = 10, units = "in", device = cairo_pdf)

range_vals <- range(GDM_ssp585_proj_offset$reverseOffset, na.rm = TRUE)
#breaks <- seq(range_vals[1], range_vals[2], length.out = 9)  # 8 bins = 9 breakpoints
quantile_breaks <- quantile(GDM_ssp585_proj_offset$reverseOffset, probs = seq(0, 1, length.out = 9), na.rm = TRUE)
#bin_labels <- paste0(format(round(breaks[-length(breaks)], 2), nsmall = 2), "–", format(round(breaks[-1], 2), nsmall = 2))
bin_labels <- paste0(format(round(quantile_breaks[-length(quantile_breaks)], 5), nsmall = 2), "–", format(round(quantile_breaks[-1], 5), nsmall = 2))
fill = cut(GDM_ssp585_proj_offset$reverseOffset, breaks = quantile_breaks, include.lowest = TRUE, labels = bin_labels)
colors <- colorRampPalette(brewer.pal(11, "Spectral")[c(6:1)])(8)
p <- ggplot(data = GDM_ssp585_proj_offset) + 
  geom_sf(data = admin, fill=gray(.9), size=0) +
  geom_raster(aes(x = x1, y = y1, fill = fill), alpha = 1) + 
  scale_fill_manual(values = colors, labels = bin_labels, guide = guide_legend(title="Reverse genomic offset", title.position = "top", title.hjust = 0.5, ncol = 1, label.position="right"), na.translate = F) +
  geom_sf(data = admin, fill=NA, color = "black", size = 0.2) +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"])) + #, expand = FALSE
  xlab("Longitude") + ylab("Latitude") +
  facet_grid(~ Date) +
  theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11, face = "bold"))
ggsave("GDM_ssp585_genomic_reverse_offset.svg", plot = p, width = 15, height = 10, units = "in", device = "svg")
ggsave("GDM_ssp585_genomic_reverse_offset.pdf", plot = p, width = 15, height = 10, units = "in", device = cairo_pdf)

GDM_ssp126_proj2040_localOffset_df <- GDM_ssp126_proj_offset %>% dplyr::filter(Date == "2040") %>% dplyr::select(x = x1, y = y1, z = localOffset) 
GDM_ssp126_proj2040_localOffset_r <- terra::rast(GDM_ssp126_proj2040_localOffset_df, type = "xyz", crs = "EPSG:4326")
terra::writeRaster(GDM_ssp126_proj2040_localOffset_r, filename = "GDM_ssp126_local_offset_2040.tif", filetype = "GTiff", overwrite = TRUE)
GDM_ssp126_proj2040_reverseOffset_df <- GDM_ssp126_proj_offset %>% dplyr::filter(Date == "2040") %>% dplyr::select(x = x1, y = y1, z = reverseOffset) 
GDM_ssp126_proj2040_reverseOffset_r <- terra::rast(GDM_ssp126_proj2040_reverseOffset_df, type = "xyz", crs = "EPSG:4326")
terra::writeRaster(GDM_ssp126_proj2040_reverseOffset_r, filename = "GDM_ssp126_reverse_offset_2040.tif", filetype = "GTiff", overwrite = TRUE)
GDM_ssp126_proj2070_localOffset_df <- GDM_ssp126_proj_offset %>% dplyr::filter(Date == "2070") %>% dplyr::select(x = x1, y = y1, z = localOffset) 
GDM_ssp126_proj2070_localOffset_r <- terra::rast(GDM_ssp126_proj2070_localOffset_df, type = "xyz", crs = "EPSG:4326")
terra::writeRaster(GDM_ssp126_proj2070_localOffset_r, filename = "GDM_ssp126_local_offset_2070.tif", filetype = "GTiff", overwrite = TRUE)
GDM_ssp126_proj2070_reverseOffset_df <- GDM_ssp126_proj_offset %>% dplyr::filter(Date == "2070") %>% dplyr::select(x = x1, y = y1, z = reverseOffset) 
GDM_ssp126_proj2070_reverseOffset_r <- terra::rast(GDM_ssp126_proj2070_reverseOffset_df, type = "xyz", crs = "EPSG:4326")
terra::writeRaster(GDM_ssp126_proj2070_reverseOffset_r, filename = "GDM_ssp126_reverse_offset_2070.tif", filetype = "GTiff", overwrite = TRUE)
GDM_ssp126_proj2100_localOffset_df <- GDM_ssp126_proj_offset %>% dplyr::filter(Date == "2100") %>% dplyr::select(x = x1, y = y1, z = localOffset) 
GDM_ssp126_proj2100_localOffset_r <- terra::rast(GDM_ssp126_proj2100_localOffset_df, type = "xyz", crs = "EPSG:4326")
terra::writeRaster(GDM_ssp126_proj2100_localOffset_r, filename = "GDM_ssp126_local_offset_2100.tif", filetype = "GTiff", overwrite = TRUE)
GDM_ssp126_proj2100_reverseOffset_df <- GDM_ssp126_proj_offset %>% dplyr::filter(Date == "2100") %>% dplyr::select(x = x1, y = y1, z = reverseOffset) 
GDM_ssp126_proj2100_reverseOffset_r <- terra::rast(GDM_ssp126_proj2100_reverseOffset_df, type = "xyz", crs = "EPSG:4326")
terra::writeRaster(GDM_ssp126_proj2100_reverseOffset_r, filename = "GDM_ssp126_reverse_offset_2100.tif", filetype = "GTiff", overwrite = TRUE)

GDM_ssp585_proj2040_localOffset_df <- GDM_ssp585_proj_offset %>% dplyr::filter(Date == "2040") %>% dplyr::select(x = x1, y = y1, z = localOffset) 
GDM_ssp585_proj2040_localOffset_r <- terra::rast(GDM_ssp585_proj2040_localOffset_df, type = "xyz", crs = "EPSG:4326")
terra::writeRaster(GDM_ssp585_proj2040_localOffset_r, filename = "GDM_ssp585_local_offset_2040.tif", filetype = "GTiff", overwrite = TRUE)
GDM_ssp585_proj2040_reverseOffset_df <- GDM_ssp585_proj_offset %>% dplyr::filter(Date == "2040") %>% dplyr::select(x = x1, y = y1, z = reverseOffset) 
GDM_ssp585_proj2040_reverseOffset_r <- terra::rast(GDM_ssp585_proj2040_reverseOffset_df, type = "xyz", crs = "EPSG:4326")
terra::writeRaster(GDM_ssp585_proj2040_reverseOffset_r, filename = "GDM_ssp585_reverse_offset_2040.tif", filetype = "GTiff", overwrite = TRUE)
GDM_ssp585_proj2070_localOffset_df <- GDM_ssp585_proj_offset %>% dplyr::filter(Date == "2070") %>% dplyr::select(x = x1, y = y1, z = localOffset) 
GDM_ssp585_proj2070_localOffset_r <- terra::rast(GDM_ssp585_proj2070_localOffset_df, type = "xyz", crs = "EPSG:4326")
terra::writeRaster(GDM_ssp585_proj2070_localOffset_r, filename = "GDM_ssp585_local_offset_2070.tif", filetype = "GTiff", overwrite = TRUE)
GDM_ssp585_proj2070_reverseOffset_df <- GDM_ssp585_proj_offset %>% dplyr::filter(Date == "2070") %>% dplyr::select(x = x1, y = y1, z = reverseOffset) 
GDM_ssp585_proj2070_reverseOffset_r <- terra::rast(GDM_ssp585_proj2070_reverseOffset_df, type = "xyz", crs = "EPSG:4326")
terra::writeRaster(GDM_ssp585_proj2070_reverseOffset_r, filename = "GDM_ssp585_reverse_offset_2070.tif", filetype = "GTiff", overwrite = TRUE)
GDM_ssp585_proj2100_localOffset_df <- GDM_ssp585_proj_offset %>% dplyr::filter(Date == "2100") %>% dplyr::select(x = x1, y = y1, z = localOffset) 
GDM_ssp585_proj2100_localOffset_r <- terra::rast(GDM_ssp585_proj2100_localOffset_df, type = "xyz", crs = "EPSG:4326")
terra::writeRaster(GDM_ssp585_proj2100_localOffset_r, filename = "GDM_ssp585_local_offset_2100.tif", filetype = "GTiff", overwrite = TRUE)
GDM_ssp585_proj2100_reverseOffset_df <- GDM_ssp585_proj_offset %>% dplyr::filter(Date == "2100") %>% dplyr::select(x = x1, y = y1, z = reverseOffset) 
GDM_ssp585_proj2100_reverseOffset_r <- terra::rast(GDM_ssp585_proj2100_reverseOffset_df, type = "xyz", crs = "EPSG:4326")
terra::writeRaster(GDM_ssp585_proj2100_reverseOffset_r, filename = "GDM_ssp585_reverse_offset_2100.tif", filetype = "GTiff", overwrite = TRUE)


# 6.3. Predicting future genomic offset using gradientForest (modifed from Lachmuth et al. (2023) and Gougherty et al. (2021) version)
# Settings ----------------------------------------------------------------
nCores<-50
nBreak<-100

# Loading required R packages ------------------------------------------------
Sys.setenv(R_LIBS_USER="/home/jeon96/R/negishi/4.2.2-gcc-12.2.0-524zdta")
.libPaths(Sys.getenv("R_LIBS_USER"))

library(adehabitatLT)
library(fields)
library(reshape)
library(parallel)
library(doParallel)
library(geosphere)
library(gradientForest)
library(foreach)
library(raster)
library(tidyverse)
library(terra)
library(sf)
library(missMDA)

rast.aggr_present <- aggregate(rast(ras_present), fact = 10, fun = 'mean', na.rm = TRUE, cores = nCores)

rast.aggr_ssp126_2011_2040 <- aggregate(rast(ras_ssp126_2011_2040), fact = 10, fun = 'mean', na.rm = TRUE, cores = nCores)
rast.aggr_ssp126_2041_2070 <- aggregate(rast(ras_ssp126_2041_2070), fact = 10, fun = 'mean', na.rm = TRUE, cores = nCores)
rast.aggr_ssp126_2071_2100 <- aggregate(rast(ras_ssp126_2071_2100), fact = 10, fun = 'mean', na.rm = TRUE, cores = nCores)
rast.aggr_ssp585_2011_2040 <- aggregate(rast(ras_ssp585_2011_2040), fact = 10, fun = 'mean', na.rm = TRUE, cores = nCores)
rast.aggr_ssp585_2041_2070 <- aggregate(rast(ras_ssp585_2041_2070), fact = 10, fun = 'mean', na.rm = TRUE, cores = nCores)
rast.aggr_ssp585_2071_2100 <- aggregate(rast(ras_ssp585_2071_2100), fact = 10, fun = 'mean', na.rm = TRUE, cores = nCores)

rast.aggr_scen.fut <- list(
  ssp126_2011_2040 = rast.aggr_ssp126_2011_2040,
  ssp126_2041_2070 = rast.aggr_ssp126_2041_2070,
  ssp126_2071_2100 = rast.aggr_ssp126_2071_2100,
  ssp585_2011_2040 = rast.aggr_ssp585_2011_2040,
  ssp585_2041_2070 = rast.aggr_ssp585_2041_2070,
  ssp585_2071_2100 = rast.aggr_ssp585_2071_2100
)

# Species range shapefile (download information at beginning of tutorial)
shp <- shapefile("monq_450km_90y_onLand.shp")

outliers <- readRDS("RDA_outliers.rds")

##############
#Load and prep PCA-based genetic distance data
##############
# Converting discrete genotype values to PCA-based continuous values of candidate SNPs 
Geno_outlier <- Geno_imp[,colnames(Geno_imp) %in% outliers$Loci]

Geno_outlier_pca <- prcomp(Geno_outlier, center = TRUE, scale. = TRUE)
eigvals <- Geno_outlier_pca$sdev^2
var_expl <- eigvals / sum(eigvals)
cumvar   <- cumsum(var_expl)
target <- 1.00 # change as needed
rank_pc <- length(eigvals)
if (all(cumvar < target)) {
  n_pc <- rank_pc
} else {
  n_pc <- min(which(cumvar >= target)[1], rank_pc)
}
n_pc

pc_n <- Geno_outlier_pca$x[, seq_len(n_pc), drop = FALSE]

# Reading in population locations
inds <- InfoInd_filt[,1:3]
ids  <- as.character(inds[, 1])

# Confirming the two are in the same order
all(rownames(pc_n) %in% ids)
all(rownames(pc_n) == ids)

##############
# Load and prep shapefile and climate data
##############
# Loading shapefile; optional

# Choosing predictors
predNames <- c("bio2", "bio5", "bio7", "bio8", "bio9", "bio14", "bio15", "bio18")
presClims <- rast.aggr_present #ras_present
presClims <- presClims[[predNames]]
presClims_r <- stack(presClims) 

# Creating pred data for gdm (cols = population name, long, lat, climate data)
pred <- data.frame(site = inds$id, long = inds$x, lat = inds$y, extract(presClims_r, y = inds[,c("x","y")]),
                   stringsAsFactors = FALSE)

######################
#GF model
######################
# Reading in genetic data
svs <- cbind(rownames(pc_n), pc_n)
colnames(svs) <- c("inds", unlist(colnames(pc_n)))

# Getting snp names
svPCs <- colnames(svs)[-1]

# Merging climate data and maf
svs <- merge(pred, svs, by.x = "site", by.y = "inds", all.x = TRUE)

# Creating gdm
gfMod <- gradientForest(data = svs, predictor.vars = predNames, response.vars = svPCs,
                        ntree = 500, 
                        maxLevel = log2(0.368*nrow(svs)/2), trace = T, 
                        corr.threshold = 0.70)

# Loading future climate data
#futClims <- rast.aggr_ssp126_2011_2040 #ras_ssp126_2011_2040 # future climate layers
#futClimDat <- as.data.frame(futClims, xy=TRUE, na.rm=TRUE)

# Getting all coordinates in range map in current climate
#popDat <- na.omit(as.data.frame(mask(raster(presClims), shp), xy=TRUE))
#popDat <- data.frame(distance=1, weight=1, popDat)
#popDat <- split(popDat, seq(nrow(popDat)))

###############
#Reverse offset calculation
##############
#Getting all coordinates in the range in current climate
popDat <- na.omit(as.data.frame(mask(presClims_r, shp), xy = TRUE))
popDat <- data.frame(popDat[, c("x","y")], predict(gfMod, popDat[,predNames]))

for (scen_fut in names(rast.aggr_scen.fut)) {
  #Gets climate data from the range in future climate
  futClims <- rast.aggr_scen.fut[[scen_fut]] #ras_ssp126_2011_2040 # future climate layers
  futClims <- futClims[[predNames]]
  futClims_r <- stack(futClims)

  futClimDat <- na.omit(as.data.frame(mask(futClims_r, shp), xy = TRUE))

  # Setting up for prediction
  futClimDat <- data.frame(futClimDat[, c("x","y")], predict(gfMod, futClimDat[, predNames]))


  cl <- makeCluster(15)
  registerDoParallel(cl)
  parallel::clusterEvalQ(cl, {
    .libPaths("/home/jeon96/R/negishi/4.2.2-gcc-12.2.0-524zdta")
    library(gradientForest)
    library(fields)
    library(geosphere)
  })

  reverseOffsetGF <- foreach(i = 1:nrow(futClimDat), .packages = c("fields","gradientForest","geosphere")) %dopar%{
  
    #get the focal population
    onePop <- futClimDat[i, ]
  
    combinedDat <- popDat[, c("x","y")]
    combinedDat["gfOffset"] <- c(rdist(onePop[, predNames], popDat[, predNames]))
  
    ##Get metrics for the focal population
    #coordinate of focal population
    coord <- onePop[, c("x","y")]
  
    #choose the pixels with the minimum offset
    minCoords <- combinedDat[which(combinedDat$gfOffset == min(combinedDat$gfOffset)), ]
  
    #calculate the distance to the sites with minimum offset, and selct the one with the shortest distance
    minCoords["dists"] <- distGeo(p1 = coord, p2 = minCoords[, 1:2])
    minCoords <- minCoords[which(minCoords$dists == min(minCoords$dists)), ]
    
    if (nrow(minCoords) > 1) {
       cat(sprintf(
       "Ties at focal (%.6f, %.6f): %s\n",
       coord$x, coord$y,
       paste(sprintf("(%.6f, %.6f)", minCoords$x, minCoords$y), collapse = "; ")
       ))
    }

    #if multiple sites have the same offset, and same distance, one is randomly chosen
    minCoords <- minCoords[sample(1:nrow(minCoords), 1), ]
  
    #get local offset
    offset <- combinedDat[which(combinedDat$x == coord$x & combinedDat$y == coord$y), "gfOffset"]
  
    #get the minimum predicted offset - reverse offset in this case
    minVal <- minCoords$gfOffset
  
    #get distance and coordinates of site that minimizes offset
    toGo <- minCoords$dists
    minPt <- minCoords[, c("x", "y")]
  
    #get bearing to the site that minimizes offset
    bear <- bearing(coord, minPt)
  
    #write out
    out <- c(x1=coord[[1]], y1=coord[[2]], localOffset=offset, reverseOffset=minVal, predGeoDist=toGo, bearing=bear, x2=minPt[[1]], y2=minPt[[2]])
  
  }

  stopCluster(cl)

  #in this resultant dataframe the columns are:
  #x1/y1: focal coordinates
  #localOffset: local offset
  #reverseOffset: reverse offset based on gradient forest
  #predGeoDist: distance to site of reverse offset
  #bearing: bearing to site of reverse offset
  #x2/y2: coordinate of site of reverse offset
  
  reverseOffsetGF <- do.call(rbind, reverseOffsetGF)

  write.csv(reverseOffsetGF, paste0("./reverseOffsetGF_", scen_fut, ".csv"), row.names = FALSE)
}

GF_ssp126_proj2040 <- read.csv("reverseOffsetGF_ssp126_2011_2040.csv", header = T)
GF_ssp126_proj2070 <- read.csv("reverseOffsetGF_ssp126_2041_2070.csv", header = T)
GF_ssp126_proj2100 <- read.csv("reverseOffsetGF_ssp126_2071_2100.csv", header = T)

GF_ssp126_proj_offset <- data.frame(rbind(GF_ssp126_proj2040, GF_ssp126_proj2070, GF_ssp126_proj2100), Date = c(rep("2040", nrow(GF_ssp126_proj2040)), rep("2070", nrow(GF_ssp126_proj2070)), rep("2100", nrow(GF_ssp126_proj2100))))

GF_ssp585_proj2040 <- read.csv("reverseOffsetGF_ssp585_2011_2040.csv", header = T)
GF_ssp585_proj2070 <- read.csv("reverseOffsetGF_ssp585_2041_2070.csv", header = T)
GF_ssp585_proj2100 <- read.csv("reverseOffsetGF_ssp585_2071_2100.csv", header = T)

GF_ssp585_proj_offset <- data.frame(rbind(GF_ssp585_proj2040, GF_ssp585_proj2070, GF_ssp585_proj2100), Date = c(rep("2040", nrow(GF_ssp585_proj2040)), rep("2070", nrow(GF_ssp585_proj2070)), rep("2100", nrow(GF_ssp585_proj2100))))

# Projecting genomic offset on a map
range_vals <- range(GF_ssp126_proj_offset$localOffset, na.rm = TRUE)
#breaks <- seq(range_vals[1], range_vals[2], length.out = 9)  # 8 bins = 9 breakpoints
quantile_breaks <- quantile(GF_ssp126_proj_offset$localOffset, probs = seq(0, 1, length.out = 9), na.rm = TRUE)
#bin_labels <- paste0(format(round(breaks[-length(breaks)], 2), nsmall = 2), "–", format(round(breaks[-1], 2), nsmall = 2))
bin_labels <- paste0(format(round(quantile_breaks[-length(quantile_breaks)], 2), nsmall = 2), "–", format(round(quantile_breaks[-1], 2), nsmall = 2))
fill = cut(GF_ssp126_proj_offset$localOffset, breaks = quantile_breaks, include.lowest = TRUE, labels = bin_labels)
colors <- colorRampPalette(brewer.pal(11, "Spectral")[c(6:1)])(8)
p <- ggplot(data = GF_ssp126_proj_offset) + 
  geom_sf(data = admin, fill=gray(.9), size=0) +
  geom_raster(aes(x = x1, y = y1, fill = fill), alpha = 1) + 
  scale_fill_manual(values = colors, labels = bin_labels, guide = guide_legend(title="Local genomic offset", title.position = "top", title.hjust = 0.5, ncol = 1, label.position="right"), na.translate = F) +
  geom_sf(data = admin, fill=NA, color = "black", size = 0.2) +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"])) + #, expand = FALSE
  xlab("Longitude") + ylab("Latitude") +
  facet_grid(~ Date) +
  theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11, face = "bold"))
ggsave("GF_ssp126_genomic_local_offset.svg", plot = p, width = 15, height = 10, units = "in", device = "svg")
ggsave("GF_ssp126_genomic_local_offset.pdf", plot = p, width = 15, height = 10, units = "in", device = cairo_pdf)

range_vals <- range(GF_ssp126_proj_offset$reverseOffset, na.rm = TRUE)
#breaks <- seq(range_vals[1], range_vals[2], length.out = 9)  # 8 bins = 9 breakpoints
quantile_breaks <- quantile(GF_ssp126_proj_offset$reverseOffset, probs = seq(0, 1, length.out = 9), na.rm = TRUE)
#bin_labels <- paste0(format(round(breaks[-length(breaks)], 2), nsmall = 2), "–", format(round(breaks[-1], 2), nsmall = 2))
bin_labels <- paste0(format(round(quantile_breaks[-length(quantile_breaks)], 5), nsmall = 2), "–", format(round(quantile_breaks[-1], 5), nsmall = 2))
fill = cut(GF_ssp126_proj_offset$reverseOffset, breaks = quantile_breaks, include.lowest = TRUE, labels = bin_labels)
colors <- colorRampPalette(brewer.pal(11, "Spectral")[c(6:1)])(8)
p <- ggplot(data = GF_ssp126_proj_offset) + 
  geom_sf(data = admin, fill=gray(.9), size=0) +
  geom_raster(aes(x = x1, y = y1, fill = fill), alpha = 1) + 
  scale_fill_manual(values = colors, labels = bin_labels, guide = guide_legend(title="Reverse genomic offset", title.position = "top", title.hjust = 0.5, ncol = 1, label.position="right"), na.translate = F) +
  geom_sf(data = admin, fill=NA, color = "black", size = 0.2) +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"])) + #, expand = FALSE
  xlab("Longitude") + ylab("Latitude") +
  facet_grid(~ Date) +
  theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11, face = "bold"))
ggsave("GF_ssp126_genomic_reverse_offset.svg", plot = p, width = 15, height = 10, units = "in", device = "svg")
ggsave("GF_ssp126_genomic_reverse_offset.pdf", plot = p, width = 15, height = 10, units = "in", device = cairo_pdf)

range_vals <- range(GF_ssp585_proj_offset$localOffset, na.rm = TRUE)
#breaks <- seq(range_vals[1], range_vals[2], length.out = 9)  # 8 bins = 9 breakpoints
quantile_breaks <- quantile(GF_ssp585_proj_offset$localOffset, probs = seq(0, 1, length.out = 9), na.rm = TRUE)
#bin_labels <- paste0(format(round(breaks[-length(breaks)], 2), nsmall = 2), "–", format(round(breaks[-1], 2), nsmall = 2))
bin_labels <- paste0(format(round(quantile_breaks[-length(quantile_breaks)], 2), nsmall = 2), "–", format(round(quantile_breaks[-1], 2), nsmall = 2))
fill = cut(GF_ssp585_proj_offset$localOffset, breaks = quantile_breaks, include.lowest = TRUE, labels = bin_labels)
colors <- colorRampPalette(brewer.pal(11, "Spectral")[c(6:1)])(8)
p <- ggplot(data = GF_ssp585_proj_offset) + 
  geom_sf(data = admin, fill=gray(.9), size=0) +
  geom_raster(aes(x = x1, y = y1, fill = fill), alpha = 1) + 
  scale_fill_manual(values = colors, labels = bin_labels, guide = guide_legend(title="Local genomic offset", title.position = "top", title.hjust = 0.5, ncol = 1, label.position="right"), na.translate = F) +
  geom_sf(data = admin, fill=NA, color = "black", size = 0.2) +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"])) + #, expand = FALSE
  xlab("Longitude") + ylab("Latitude") +
  facet_grid(~ Date) +
  theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11, face = "bold"))
ggsave("GF_ssp585_genomic_local_offset.svg", plot = p, width = 15, height = 10, units = "in", device = "svg")
ggsave("GF_ssp585_genomic_local_offset.pdf", plot = p, width = 15, height = 10, units = "in", device = cairo_pdf)

range_vals <- range(GF_ssp585_proj_offset$reverseOffset, na.rm = TRUE)
#breaks <- seq(range_vals[1], range_vals[2], length.out = 9)  # 8 bins = 9 breakpoints
quantile_breaks <- quantile(GF_ssp585_proj_offset$reverseOffset, probs = seq(0, 1, length.out = 9), na.rm = TRUE)
#bin_labels <- paste0(format(round(breaks[-length(breaks)], 2), nsmall = 2), "–", format(round(breaks[-1], 2), nsmall = 2))
bin_labels <- paste0(format(round(quantile_breaks[-length(quantile_breaks)], 5), nsmall = 2), "–", format(round(quantile_breaks[-1], 5), nsmall = 2))
fill = cut(GF_ssp585_proj_offset$reverseOffset, breaks = quantile_breaks, include.lowest = TRUE, labels = bin_labels)
colors <- colorRampPalette(brewer.pal(11, "Spectral")[c(6:1)])(8)
p <- ggplot(data = GF_ssp585_proj_offset) + 
  geom_sf(data = admin, fill=gray(.9), size=0) +
  geom_raster(aes(x = x1, y = y1, fill = fill), alpha = 1) + 
  scale_fill_manual(values = colors, labels = bin_labels, guide = guide_legend(title="Reverse genomic offset", title.position = "top", title.hjust = 0.5, ncol = 1, label.position="right"), na.translate = F) +
  geom_sf(data = admin, fill=NA, color = "black", size = 0.2) +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"])) + #, expand = FALSE
  xlab("Longitude") + ylab("Latitude") +
  facet_grid(~ Date) +
  theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11, face = "bold"))
ggsave("GF_ssp585_genomic_reverse_offset.svg", plot = p, width = 15, height = 10, units = "in", device = "svg")
ggsave("GF_ssp585_genomic_reverse_offset.pdf", plot = p, width = 15, height = 10, units = "in", device = cairo_pdf)

GF_ssp126_proj2040_localOffset_df <- GF_ssp126_proj_offset %>% dplyr::filter(Date == "2040") %>% dplyr::select(x = x1, y = y1, z = localOffset) 
GF_ssp126_proj2040_localOffset_r <- terra::rast(GF_ssp126_proj2040_localOffset_df, type = "xyz", crs = "EPSG:4326")
terra::writeRaster(GF_ssp126_proj2040_localOffset_r, filename = "GF_ssp126_local_offset_2040.tif", filetype = "GTiff", overwrite = TRUE)
GF_ssp126_proj2040_reverseOffset_df <- GF_ssp126_proj_offset %>% dplyr::filter(Date == "2040") %>% dplyr::select(x = x1, y = y1, z = reverseOffset) 
GF_ssp126_proj2040_reverseOffset_r <- terra::rast(GF_ssp126_proj2040_reverseOffset_df, type = "xyz", crs = "EPSG:4326")
terra::writeRaster(GF_ssp126_proj2040_reverseOffset_r, filename = "GF_ssp126_reverse_offset_2040.tif", filetype = "GTiff", overwrite = TRUE)
GF_ssp126_proj2070_localOffset_df <- GF_ssp126_proj_offset %>% dplyr::filter(Date == "2070") %>% dplyr::select(x = x1, y = y1, z = localOffset) 
GF_ssp126_proj2070_localOffset_r <- terra::rast(GF_ssp126_proj2070_localOffset_df, type = "xyz", crs = "EPSG:4326")
terra::writeRaster(GF_ssp126_proj2070_localOffset_r, filename = "GF_ssp126_local_offset_2070.tif", filetype = "GTiff", overwrite = TRUE)
GF_ssp126_proj2070_reverseOffset_df <- GF_ssp126_proj_offset %>% dplyr::filter(Date == "2070") %>% dplyr::select(x = x1, y = y1, z = reverseOffset) 
GF_ssp126_proj2070_reverseOffset_r <- terra::rast(GF_ssp126_proj2070_reverseOffset_df, type = "xyz", crs = "EPSG:4326")
terra::writeRaster(GF_ssp126_proj2070_reverseOffset_r, filename = "GF_ssp126_reverse_offset_2070.tif", filetype = "GTiff", overwrite = TRUE)
GF_ssp126_proj2100_localOffset_df <- GF_ssp126_proj_offset %>% dplyr::filter(Date == "2100") %>% dplyr::select(x = x1, y = y1, z = localOffset) 
GF_ssp126_proj2100_localOffset_r <- terra::rast(GF_ssp126_proj2100_localOffset_df, type = "xyz", crs = "EPSG:4326")
terra::writeRaster(GF_ssp126_proj2100_localOffset_r, filename = "GF_ssp126_local_offset_2100.tif", filetype = "GTiff", overwrite = TRUE)
GF_ssp126_proj2100_reverseOffset_df <- GF_ssp126_proj_offset %>% dplyr::filter(Date == "2100") %>% dplyr::select(x = x1, y = y1, z = reverseOffset) 
GF_ssp126_proj2100_reverseOffset_r <- terra::rast(GF_ssp126_proj2100_reverseOffset_df, type = "xyz", crs = "EPSG:4326")
terra::writeRaster(GF_ssp126_proj2100_reverseOffset_r, filename = "GF_ssp126_reverse_offset_2100.tif", filetype = "GTiff", overwrite = TRUE)

GF_ssp585_proj2040_localOffset_df <- GF_ssp585_proj_offset %>% dplyr::filter(Date == "2040") %>% dplyr::select(x = x1, y = y1, z = localOffset) 
GF_ssp585_proj2040_localOffset_r <- terra::rast(GF_ssp585_proj2040_localOffset_df, type = "xyz", crs = "EPSG:4326")
terra::writeRaster(GF_ssp585_proj2040_localOffset_r, filename = "GF_ssp585_local_offset_2040.tif", filetype = "GTiff", overwrite = TRUE)
GF_ssp585_proj2040_reverseOffset_df <- GF_ssp585_proj_offset %>% dplyr::filter(Date == "2040") %>% dplyr::select(x = x1, y = y1, z = reverseOffset) 
GF_ssp585_proj2040_reverseOffset_r <- terra::rast(GF_ssp585_proj2040_reverseOffset_df, type = "xyz", crs = "EPSG:4326")
terra::writeRaster(GF_ssp585_proj2040_reverseOffset_r, filename = "GF_ssp585_reverse_offset_2040.tif", filetype = "GTiff", overwrite = TRUE)
GF_ssp585_proj2070_localOffset_df <- GF_ssp585_proj_offset %>% dplyr::filter(Date == "2070") %>% dplyr::select(x = x1, y = y1, z = localOffset) 
GF_ssp585_proj2070_localOffset_r <- terra::rast(GF_ssp585_proj2070_localOffset_df, type = "xyz", crs = "EPSG:4326")
terra::writeRaster(GF_ssp585_proj2070_localOffset_r, filename = "GF_ssp585_local_offset_2070.tif", filetype = "GTiff", overwrite = TRUE)
GF_ssp585_proj2070_reverseOffset_df <- GF_ssp585_proj_offset %>% dplyr::filter(Date == "2070") %>% dplyr::select(x = x1, y = y1, z = reverseOffset) 
GF_ssp585_proj2070_reverseOffset_r <- terra::rast(GF_ssp585_proj2070_reverseOffset_df, type = "xyz", crs = "EPSG:4326")
terra::writeRaster(GF_ssp585_proj2070_reverseOffset_r, filename = "GF_ssp585_reverse_offset_2070.tif", filetype = "GTiff", overwrite = TRUE)
GF_ssp585_proj2100_localOffset_df <- GF_ssp585_proj_offset %>% dplyr::filter(Date == "2100") %>% dplyr::select(x = x1, y = y1, z = localOffset) 
GF_ssp585_proj2100_localOffset_r <- terra::rast(GF_ssp585_proj2100_localOffset_df, type = "xyz", crs = "EPSG:4326")
terra::writeRaster(GF_ssp585_proj2100_localOffset_r, filename = "GF_ssp585_local_offset_2100.tif", filetype = "GTiff", overwrite = TRUE)
GF_ssp585_proj2100_reverseOffset_df <- GF_ssp585_proj_offset %>% dplyr::filter(Date == "2100") %>% dplyr::select(x = x1, y = y1, z = reverseOffset) 
GF_ssp585_proj2100_reverseOffset_r <- terra::rast(GF_ssp585_proj2100_reverseOffset_df, type = "xyz", crs = "EPSG:4326")
terra::writeRaster(GF_ssp585_proj2100_reverseOffset_r, filename = "GF_ssp585_reverse_offset_2100.tif", filetype = "GTiff", overwrite = TRUE)



# Phenotype-Genotype regression
# Load data
Geno_imp <- readRDS("Geno_imp_final.rds")
outliers <- readRDS("RDA_outliers_final.rds")
Geno_outlier <- Geno_imp[,colnames(Geno_imp) %in% outliers$Loci]

Variables <- read.csv("merged_variables_final.csv", header = T, row.names = 1)
bio <- Variables[, c("id", grep("^bio", colnames(Variables), value = TRUE))]
rownames(bio) <- bio$id
bio$id <- NULL
pop <- Variables[, c("id", grep("^PC", colnames(Variables), value = TRUE))]
rownames(pop) <- pop$id
pop$id <- NULL

Phenotypes <- read.csv("monq_phenotypes.csv", header = T) # mean and sd of L, a, b, mass

Y_raw <- Phenotypes[Phenotypes$ID != "F1913", c("ID", "avg.L", "avg.a", "avg.b", "mass", "age")]
rownames(Y_raw) <- Y_raw$ID
Y_raw$ID <- NULL
Y <- scale(Y_raw[, c("avg.L", "avg.a", "avg.b")])

# Converting scaled bioclim variables to PCA axes
bio_pca <- prcomp(bio, center=TRUE, scale.=TRUE)
eigvals <- bio_pca$sdev^2
var_expl <- eigvals / sum(eigvals)
cumvar   <- cumsum(var_expl)
target <- 1.00 # change as needed
rank_pc <- length(eigvals)
if (all(cumvar < target)) {
  n_pc <- rank_pc
} else {
  n_pc <- min(which(cumvar >= target)[1], rank_pc)
}
n_pc # 8

bio_pheno_pcx <- bio_pca$x[rownames(bio_pca$x) %in% Phenotypes$ID,]
pop_pheno <- pop[rownames(pop) %in% Phenotypes$ID,]

# Converting discrete genotype values to PCA-based continuous values of candidate SNPs 
Geno_outlier_pca <- prcomp(Geno_outlier, center = TRUE, scale. = TRUE)

eigvals <- Geno_outlier_pca$sdev^2
var_expl <- eigvals / sum(eigvals)
cumvar   <- cumsum(var_expl)
target <- 0.95 # change as needed
rank_pc <- length(eigvals)
if (all(cumvar < target)) {
  n_pc <- rank_pc
} else {
  n_pc <- min(which(cumvar >= target)[1], rank_pc)
}
n_pc # 66

Geno_outlier_pheno_pcx <- Geno_outlier_pca$x[rownames(Geno_outlier_pca$x) %in% Phenotypes$ID,]

X <- as.data.frame(Geno_outlier_pheno_pcx[, seq_len(n_pc)]) 
colnames(X) <- paste0("genoPC", 1:n_pc)

all(rownames(bio_pheno_pcx) == rownames(X))
all(rownames(pop_pheno) == rownames(X))
all(rownames(Y) == rownames(X))

COV <- data.frame(strPC1 = pop_pheno[, 1], strPC2 = pop_pheno[, 2],
                  bioPC1_2      = bio_pheno_pcx[, 1:2])
colnames(COV) <- c("strPC1", "strPC2", paste0("bioPC", 1:2))

# iterate for K = 1:12
Kmax <- 12

rda_results <- lapply(1:Kmax, function(k) {

  Xk <- as.data.frame(Geno_outlier_pheno_pcx[, 1:k, drop = FALSE])
  colnames(Xk) <- paste0("genoPC", 1:k)

  pheno_rda <- rda(Y ~ ., data = Xk, Condition = COV)

  a <- anova(pheno_rda, permutations = 999)

  data.frame(
    k = k,
    p_value = a$`Pr(>F)`[1],
    R2 = RsquareAdj(pheno_rda)$r.squared,
    R2_adj = RsquareAdj(pheno_rda)$adj.r.squared,
    constrained_rank = pheno_rda$CCA$rank,
    residual_rank = pheno_rda$CA$rank
  )
})

rda_results <- do.call(rbind, rda_results)
rda_results
#    k p_value        R2    R2_adj constrained_rank residual_rank
#1   1   0.003 0.2156534 0.1854862                1             3
#2   2   0.008 0.2545620 0.1949269                2             3
#3   3   0.017 0.2672986 0.1757109                3             3
#4   4   0.043 0.2702014 0.1432799                3             3
#5   5   0.084 0.2876153 0.1257097                3             3
#6   6   0.088 0.3361536 0.1464832                3             3
#7   7   0.090 0.3719791 0.1521718                3             3
#8   8   0.050 0.4392720 0.2031760                3             3
#9   9   0.034 0.4995532 0.2493298                3             3
#10 10   0.056 0.5048179 0.2135344                3             3
#11 11   0.061 0.5443074 0.2310188                3             3
#12 12   0.027 0.6119104 0.3014386                3             3


# p-value vs. k
svg("rda_pvalue_k.svg")
plot(rda_results$k, rda_results$p_value, type = "b",
     ylab = "Permutation p-value",
     xlab = "Number of genoPCs included",
     ylim = c(0, 1))
abline(h = 0.05, col = "red", lty = 2)
dev.off()

# R2_adj vs. k
svg("rda_adjR2_k.svg")
plot(rda_results$k, rda_results$R2_adj, type = "b",
     ylab = "Adjusted R²",
     xlab = "Number of genoPCs included")
dev.off()

# Partial RDA: Y ~ X | covariates
X3 <- as.data.frame(Geno_outlier_pheno_pcx[, 1:3, drop = FALSE])
colnames(X3) <- paste0("genoPC", 1:3)

pheno_rda_fit <- rda(Y ~ ., data=X3, Condition = COV) # Condition(structure PC1 + structure PC2 + bioclim PC1 + bioclim PC2)

summary(pheno_rda_fit)
#Call:
#rda(formula = Y ~ genoPC1 + genoPC2 + genoPC3, data = X3, Condition = COV)

#Partitioning of variance:
#              Inertia Proportion
#Total          3.0000     1.0000
#Constrained    0.8019     0.2673
#Unconstrained  2.1981     0.7327

#Eigenvalues, and their contribution to the variance

#Importance of components:
#                        RDA1    RDA2      RDA3    PC1    PC2     PC3
#Eigenvalue            0.7176 0.08409 2.024e-04 1.4703 0.4565 0.27129
#Proportion Explained  0.2392 0.02803 6.746e-05 0.4901 0.1522 0.09043
#Cumulative Proportion 0.2392 0.26723 2.673e-01 0.7574 0.9096 1.00000

#Accumulated constrained eigenvalues
#Importance of components:
#                        RDA1    RDA2      RDA3
#Eigenvalue            0.7176 0.08409 0.0002024
#Proportion Explained  0.8949 0.10486 0.0002524
#Cumulative Proportion 0.8949 0.99975 1.0000000

# Permutation test for overall model, constrained by covariates
anova(pheno_rda_fit, permutations=999)
#Permutation test for rda under reduced model
#Permutation: free
#Number of permutations: 999

#Model: rda(formula = Y ~ genoPC1 + genoPC2 + genoPC3, data = X3, Condition = COV)
#         Df Variance      F Pr(>F)
#Model     3   0.8019 2.9185   0.01 **
#Residual 24   2.1981
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Test marginal terms (genoPC axes)
anova(pheno_rda_fit, by="terms", permutations=999)
#Permutation test for rda under reduced model
#Terms added sequentially (first to last)
#Permutation: free
#Number of permutations: 999

#Model: rda(formula = Y ~ genoPC1 + genoPC2 + genoPC3, data = X3, Condition = COV)
#         Df Variance      F Pr(>F)
#genoPC1   1  0.64696 7.0638  0.003 **
#genoPC2   1  0.11673 1.2745  0.279
#genoPC3   1  0.03821 0.4172  0.676
#Residual 24  2.19810
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Test RDA axes
anova(pheno_rda_fit, by="axis", permutations=999)
#Permutation test for rda under reduced model
#Forward tests for axes
#Permutation: free
#Number of permutations: 999

#Model: rda(formula = Y ~ genoPC1 + genoPC2 + genoPC3, data = X3, Condition = COV)
#         Df Variance      F Pr(>F)
#RDA1      1  0.71761 7.8352  0.011 *
#RDA2      1  0.08409 0.9564  0.709
#RDA3      1  0.00020 0.0024  1.000
#Residual 24  2.19810
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Effect size
RsquareAdj(pheno_rda_fit)
#$r.squared
#[1] 0.2672986

#$adj.r.squared
#[1] 0.1757109


# Univariate multiple regression
library(summarytools)
library(sjstats)
library(car)

df <- data.frame(
  L = Y_raw$avg.L,
  a = Y_raw$avg.a,
  b = Y_raw$avg.b,
  mass = Y_raw$mass,
  age = Y_raw$age,
  genoPC1 = X$genoPC1, genoPC2 = X$genoPC2, genoPC3 = X$genoPC3,
  strPC1 = COV$strPC1, strPC2 = COV$strPC2,
  bioPC1 = COV$bioPC1, bioPC2 = COV$bioPC2
)

qqnorm(df$L) 
qqline(df$L)
descr(df$L, stats = c("Skewness", "Kurtosis"), style="rmarkdown") # Skewness = 0.38, Kurtosis = -0.88

L_fit <- lm(L ~ genoPC1 + genoPC2 + strPC1 + strPC2 + bioPC1 + bioPC2, data=df)
summary(L_fit)
#Call:
#lm(formula = L ~ genoPC1 + genoPC2 + strPC1 + strPC2 + bioPC1 +
#    bioPC2, data = df)

#Residuals:
#    Min      1Q  Median      3Q     Max
#-2.8899 -1.3757 -0.1665  1.2947  3.5924

#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)
#(Intercept) 20.73945    0.96264  21.544 8.43e-16 ***
#genoPC1      0.34133    0.43539   0.784    0.442
#genoPC2      0.23593    0.38895   0.607    0.551
#strPC1      -0.32683    0.86590  -0.377    0.710
#strPC2       0.66285    0.59233   1.119    0.276
#bioPC1       0.01409    0.51152   0.028    0.978
#bioPC2       0.90616    0.67896   1.335    0.196
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 2.001 on 21 degrees of freedom
#Multiple R-squared:  0.7626,    Adjusted R-squared:  0.6948
#F-statistic: 11.24 on 6 and 21 DF,  p-value: 1.222e-05

Anova(L_fit, type = 3)
#Anova Table (Type III tests)

#Response: L
#             Sum Sq Df  F value    Pr(>F)
#(Intercept) 1859.11  1 464.1563 8.433e-16 ***
#genoPC1        2.46  1   0.6146    0.4418
#genoPC2        1.47  1   0.3679    0.5506
#strPC1         0.57  1   0.1425    0.7096
#strPC2         5.02  1   1.2523    0.2758
#bioPC1         0.00  1   0.0008    0.9783
#bioPC2         7.13  1   1.7812    0.1963
#Residuals     84.11 21
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova_stats(L_fit)
#etasq | partial.etasq | omegasq | partial.omegasq | epsilonsq | cohens.f
#------------------------------------------------------------------------
#0.397 |         0.626 |   0.382 |           0.549 |     0.386 |    1.294
#0.002 |         0.008 |  -0.009 |          -0.031 |    -0.009 |    0.090
#0.309 |         0.565 |   0.294 |           0.485 |     0.298 |    1.141
#0.034 |         0.125 |   0.022 |           0.066 |     0.022 |    0.377
#0.001 |         0.003 |  -0.011 |          -0.035 |    -0.011 |    0.053
#0.020 |         0.078 |   0.009 |           0.027 |     0.009 |    0.291
#      |               |         |                 |           |
#
#etasq |      term |   sumsq | df |  meansq | statistic | p.value | power
#------------------------------------------------------------------------
#0.397 |   genoPC1 | 140.744 |  1 | 140.744 |    35.139 |  < .001 | 1.000
#0.002 |   genoPC2 |   0.679 |  1 |   0.679 |     0.170 |   0.685 | 0.070
#0.309 |    strPC1 | 109.457 |  1 | 109.457 |    27.328 |  < .001 | 0.999
#0.034 |    strPC2 |  11.969 |  1 |  11.969 |     2.988 |   0.099 | 0.408
#0.001 |    bioPC1 |   0.233 |  1 |   0.233 |     0.058 |   0.812 | 0.057
#0.020 |    bioPC2 |   7.134 |  1 |   7.134 |     1.781 |   0.196 | 0.266
#      | Residuals |  84.113 | 21 |   4.005 |           |         |

qqnorm(df$a) 
qqline(df$a)
descr(df$a, stats = c("Skewness", "Kurtosis"), style="rmarkdown") # Skewness = 0.19, Kurtosis = -0.90

a_fit <- lm(a ~ genoPC1 + genoPC2 + strPC1 + strPC2 + bioPC1 + bioPC2, data=df)
summary(a_fit)
#Call:
#lm(formula = a ~ genoPC1 + genoPC2 + strPC1 + strPC2 + bioPC1 +
#    bioPC2, data = df)

#Residuals:
#    Min      1Q  Median      3Q     Max
#-1.9928 -0.5476 -0.1176  0.5329  1.7782

#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)
#(Intercept)  0.21862    0.49925   0.438    0.666
#genoPC1      0.02300    0.22581   0.102    0.920
#genoPC2     -0.30622    0.20172  -1.518    0.144
#strPC1      -0.06063    0.44908  -0.135    0.894
#strPC2       0.11435    0.30720   0.372    0.713
#bioPC1       0.36191    0.26529   1.364    0.187
#bioPC2       0.37191    0.35213   1.056    0.303

#Residual standard error: 1.038 on 21 degrees of freedom
#Multiple R-squared:  0.3196,    Adjusted R-squared:  0.1252
#F-statistic: 1.644 on 6 and 21 DF,  p-value: 0.1846

Anova(a_fit, type = 3)
#Anova Table (Type III tests)

#Response: a
#             Sum Sq Df F value Pr(>F)
#(Intercept)  0.2066  1  0.1917 0.6659
#genoPC1      0.0112  1  0.0104 0.9198
#genoPC2      2.4827  1  2.3045 0.1439
#strPC1       0.0196  1  0.0182 0.8939
#strPC2       0.1493  1  0.1386 0.7134
#bioPC1       2.0049  1  1.8610 0.1870
#bioPC2       1.2018  1  1.1155 0.3029
#Residuals   22.6242 21

anova_stats(a_fit)
#etasq | partial.etasq | omegasq | partial.omegasq | epsilonsq | cohens.f
#------------------------------------------------------------------------
#0.052 |         0.071 |   0.019 |           0.022 |     0.020 |    0.277
#0.078 |         0.103 |   0.044 |           0.048 |     0.045 |    0.338
#0.046 |         0.063 |   0.013 |           0.015 |     0.013 |    0.259
#0.030 |         0.042 |  -0.003 |          -0.003 |    -0.003 |    0.209
#0.078 |         0.103 |   0.044 |           0.048 |     0.045 |    0.338
#0.036 |         0.050 |   0.004 |           0.004 |     0.004 |    0.230
#      |               |         |                 |           |

#etasq |      term |  sumsq | df | meansq | statistic | p.value | power
#----------------------------------------------------------------------
#0.052 |   genoPC1 |  1.741 |  1 |  1.741 |     1.616 |   0.218 | 0.246
#0.078 |   genoPC2 |  2.585 |  1 |  2.585 |     2.400 |   0.136 | 0.340
#0.046 |    strPC1 |  1.522 |  1 |  1.522 |     1.413 |   0.248 | 0.221
#0.030 |    strPC2 |  0.988 |  1 |  0.988 |     0.917 |   0.349 | 0.160
#0.078 |    bioPC1 |  2.589 |  1 |  2.589 |     2.403 |   0.136 | 0.341
#0.036 |    bioPC2 |  1.202 |  1 |  1.202 |     1.115 |   0.303 | 0.184
#      | Residuals | 22.624 | 21 |  1.077 |           |         |

qqnorm(df$b) 
qqline(df$b)
descr(df$b, stats = c("Skewness", "Kurtosis"), style="rmarkdown") # Skewness = -0.10, Kurtosis = -1.16

b_fit <- lm(b ~ genoPC1 + genoPC2 + strPC1 + strPC2 + bioPC1 + bioPC2, data=df)
summary(b_fit)
#Call:
#lm(formula = b ~ genoPC1 + genoPC2 + strPC1 + strPC2 + bioPC1 +
#    bioPC2, data = df)

#Residuals:
#     Min       1Q   Median       3Q      Max
#-0.92296 -0.41186 -0.09974  0.47524  1.28258

#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)
#(Intercept)  0.56184    0.30405   1.848   0.0788 .
#genoPC1     -0.13174    0.13752  -0.958   0.3490
#genoPC2      0.21005    0.12285   1.710   0.1020
#strPC1       0.30638    0.27350   1.120   0.2753
#strPC2      -0.15433    0.18709  -0.825   0.4187
#bioPC1      -0.11112    0.16157  -0.688   0.4991
#bioPC2      -0.04306    0.21445  -0.201   0.8428
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.6321 on 21 degrees of freedom
#Multiple R-squared:  0.3371,    Adjusted R-squared:  0.1477
#F-statistic:  1.78 on 6 and 21 DF,  p-value: 0.1521

Anova(b_fit, type = 3)
#Anova Table (Type III tests)

#Response: b
#            Sum Sq Df F value  Pr(>F)
#(Intercept) 1.3644  1  3.4146 0.07876 .
#genoPC1     0.3667  1  0.9177 0.34897
#genoPC2     1.1681  1  2.9234 0.10204
#strPC1      0.5014  1  1.2549 0.27527
#strPC2      0.2719  1  0.6805 0.41869
#bioPC1      0.1890  1  0.4730 0.49913
#bioPC2      0.0161  1  0.0403 0.84280
#Residuals   8.3913 21
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova_stats(b_fit)
#etasq | partial.etasq | omegasq | partial.omegasq | epsilonsq | cohens.f
#------------------------------------------------------------------------
#0.197 |         0.229 |   0.161 |           0.158 |     0.166 |    0.546
#0.037 |         0.053 |   0.005 |           0.006 |     0.005 |    0.236
#0.077 |         0.105 |   0.044 |           0.049 |     0.046 |    0.342
#0.007 |         0.011 |  -0.024 |          -0.028 |    -0.024 |    0.104
#0.017 |         0.025 |  -0.014 |          -0.017 |    -0.015 |    0.159
#0.001 |         0.002 |  -0.029 |          -0.035 |    -0.030 |    0.044
#      |               |         |                 |           |

#etasq |      term | sumsq | df | meansq | statistic | p.value | power
#---------------------------------------------------------------------
#0.197 |   genoPC1 | 2.498 |  1 |  2.498 |     6.253 |   0.021 | 0.704
#0.037 |   genoPC2 | 0.469 |  1 |  0.469 |     1.174 |   0.291 | 0.191
#0.077 |    strPC1 | 0.980 |  1 |  0.980 |     2.452 |   0.132 | 0.346
#0.007 |    strPC2 | 0.091 |  1 |  0.091 |     0.229 |   0.638 | 0.077
#0.017 |    bioPC1 | 0.212 |  1 |  0.212 |     0.531 |   0.474 | 0.113
#0.001 |    bioPC2 | 0.016 |  1 |  0.016 |     0.040 |   0.843 | 0.055
#      | Residuals | 8.391 | 21 |  0.400 |           |         |

df$age <- factor(df$age, levels = c("Juvenile", "Adult"))

qqnorm(df$mass) 
qqline(df$mass)
descr(df$mass, stats = c("Skewness", "Kurtosis"), style="rmarkdown") # Skewness = -0.12, Kurtosis = -0.62

mass_fit <- lm(mass ~ genoPC1 + genoPC2 + strPC1 + strPC2 + bioPC1 + bioPC2 + age, data=df)
summary(mass_fit)
#Call:
#lm(formula = mass ~ genoPC1 + genoPC2 + strPC1 + strPC2 + bioPC1 +
#    bioPC2 + age, data = df)

#Residuals:
#    Min      1Q  Median      3Q     Max
#-21.245  -6.219   1.376   8.626  14.553

#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)
#(Intercept) 182.8618     7.4709  24.477 7.87e-16 ***
#genoPC1       1.1584     2.6991   0.429   0.6726
#genoPC2       1.9135     2.6411   0.724   0.4776
#strPC1       -1.3996     5.3998  -0.259   0.7983
#strPC2        4.2808     3.6545   1.171   0.2559
#bioPC1        1.8287     3.2853   0.557   0.5843
#bioPC2        0.9124     4.1993   0.217   0.8303
#ageAdult     13.6326     4.8498   2.811   0.0112 *
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 12.34 on 19 degrees of freedom
#  (1 observation deleted due to missingness)
#Multiple R-squared:  0.5925,    Adjusted R-squared:  0.4424
#F-statistic: 3.947 on 7 and 19 DF,  p-value: 0.008008

Anova(mass_fit, type = 3)
#Anova Table (Type III tests)
#
#Response: mass
#            Sum Sq Df  F value    Pr(>F)
#(Intercept)  91160  1 599.1060 7.871e-16 ***
#genoPC1         28  1   0.1842   0.67262
#genoPC2         80  1   0.5249   0.47759
#strPC1          10  1   0.0672   0.79828
#strPC2         209  1   1.3721   0.25593
#bioPC1          47  1   0.3098   0.58428
#bioPC2           7  1   0.0472   0.83031
#age           1202  1   7.9015   0.01115 *
#Residuals     2891 19
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova_stats(mass_fit)
#etasq | partial.etasq | omegasq | partial.omegasq | epsilonsq | cohens.f
#------------------------------------------------------------------------
#0.065 |         0.137 |   0.042 |           0.070 |     0.043 |    0.399
#0.010 |         0.023 |  -0.012 |          -0.021 |    -0.012 |    0.154
#0.320 |         0.440 |   0.293 |           0.340 |     0.299 |    0.887
#0.013 |         0.031 |  -0.008 |          -0.015 |    -0.008 |    0.179
#0.014 |         0.034 |  -0.007 |          -0.012 |    -0.007 |    0.189
#0.001 |         0.002 |  -0.020 |          -0.037 |    -0.021 |    0.044
#0.169 |         0.294 |   0.145 |           0.204 |     0.148 |    0.645
#      |               |         |                 |           |

#etasq |      term |    sumsq | df |   meansq | statistic | p.value | power
#--------------------------------------------------------------------------
#0.065 |   genoPC1 |  459.439 |  1 |  459.439 |     3.019 |   0.098 | 0.411
#0.010 |   genoPC2 |   68.540 |  1 |   68.540 |     0.450 |   0.510 | 0.103
#0.320 |    strPC1 | 2272.362 |  1 | 2272.362 |    14.934 |   0.001 | 0.971
#0.013 |    strPC2 |   93.013 |  1 |   93.013 |     0.611 |   0.444 | 0.122
#0.014 |    bioPC1 |  102.840 |  1 |  102.840 |     0.676 |   0.421 | 0.130
#0.001 |    bioPC2 |    5.660 |  1 |    5.660 |     0.037 |   0.849 | 0.054
#0.169 |       age | 1202.289 |  1 | 1202.289 |     7.901 |   0.011 | 0.800
#      | Residuals | 2891.042 | 19 |  152.160 |           |         |


# PLS
library(mixOmics)

# Tune PLS (leave-one-out validation to determine ncomp)
Ymat <- as.matrix(Y[, c("avg.L", "avg.a", "avg.b")])
storage.mode(Ymat) <- "double"
Xmat <- as.matrix(X[, c("genoPC1","genoPC2","genoPC3")])
storage.mode(Xmat) <- "double"
all(rownames(Ymat) == rownames(Xmat))

pls_fit <- pls(Xmat, Ymat, ncomp = 3, mode = 'regression')
Q2.pls <- perf(pls_fit, validation = 'Mfold', folds = 5, nrepeat = 100)
Q2.pls$measures$Q2
#$summary
#  feature comp       mean         sd
#1   avg.L    1  0.3390614 0.02689963
#2   avg.L    2 -0.2512991 0.11641088
#3   avg.L    3 -0.4432416 0.21683320
#4   avg.a    1 -0.1602426 0.10552970
#5   avg.a    2 -0.5064958 0.39411224
#6   avg.a    3 -0.7344124 0.35090249
#7   avg.b    1  0.1007850 0.07459838
#8   avg.b    2 -0.4433203 0.26042525
#9   avg.b    3 -0.6330671 0.32936694

svg("Q2.pls.svg")
plot(Q2.pls, criterion = 'Q2') # ncomp = 1; only L 
dev.off()

#list.keepX = c(1:3)
#pls_tune <- tune.spls(X, Y, ncomp = 1, test.keepX = list.keepX, validation = 'Mfold', folds = 5, nrepeat = 100, mode = 'regression', measure = 'cor')
#plot(pls_tune) 

# full model
pheno_pls <- pls(Xmat, Ymat, ncomp = 1, mode = "regression")
pheno_pls_perf <- perf(pheno_pls, validation = 'Mfold', folds = 5, nrepeat = 100)
pheno_pls_perf $measures$Q2
#$summary
#  feature comp       mean         sd
#1   avg.L    1  0.3375866 0.03434003
#2   avg.a    1 -0.1671632 0.11214345
#3   avg.b    1  0.1080429 0.07132048

# summary
#summary(pheno_pls) # error
pheno_pls$prop_expl_var$X   # X explained variance (comp1); 0.911977
pheno_pls$prop_expl_var$Y   # Y explained variance (comp1); 0.6814907
head(pheno_pls$variates$X)  # genotype latent score (t)
head(pheno_pls$variates$Y)  # phenotype latent score (u)
pheno_pls$variates$X[,1] # = genotype PLS1 score
pheno_pls$variates$Y[,1] # = phenotype PLS1 score

pheno_pls$loadings$X # genotype's contribution; genoPC1 -0.5760671, genoPC2 -0.5606078, genoPC3  0.5948659
pheno_pls$loadings$Y # phenotype's contribution; avg.L -0.7512626, avg.a  0.3343098, avg.b -0.5690707

#plotIndiv(pheno_pls, comp = 1)
#plotVar(pheno_pls, comp = 1)
plotLoadings(pheno_pls, comp = 1, contrib = "max")

t1 <- pheno_pls$variates$X[,1]   # genotype score
u1 <- pheno_pls$variates$Y[,1]   # phenotype score

svg("full_pls_corr.svg")
plot(t1, u1, xlab="PLS1 score (X: genotype)", ylab="PLS1 score (Y: phenotype)")
abline(lm(u1 ~ t1), lty=2)
dev.off()

cor(t1, u1, use="complete.obs") # 0.6034859
cor(t1, Ymat[, "avg.L"]) # -0.6330237
cor(t1, Ymat[, "avg.a"]) # 0.2816938
cor(t1, Ymat[, "avg.b"]) # -0.4795064

# permutation (optional)
set.seed(123)

obs_cor <- cor(
  pheno_pls$variates$X[,1],
  pheno_pls$variates$Y[,1]
)

nperm <- 999
perm_cor <- numeric(nperm)

for (i in 1:nperm) {
  Yperm <- Ymat[sample(nrow(Ymat)), ]  # phenotype permute
  pls_perm <- pls(Xmat, Yperm, ncomp = 1, mode = "regression")
  perm_cor[i] <- cor(pls_perm$variates$X[,1],
                     pls_perm$variates$Y[,1])
}

p_perm <- mean(abs(perm_cor) >= abs(obs_cor))
p_perm # 0.002002002

# residual model - str + bio
Y_COV <- cbind(Y, COV)
L_res <- residuals(lm(avg.L ~ strPC1 + strPC2 + bioPC1 + bioPC2, data=Y_COV))
a_res <- residuals(lm(avg.a ~ strPC1 + strPC2 + bioPC1 + bioPC2, data=Y_COV))
b_res <- residuals(lm(avg.b ~ strPC1 + strPC2 + bioPC1 + bioPC2, data=Y_COV))

Y_res <- cbind(L_res, a_res, b_res)
rownames(Y_res) <- rownames(Y_COV)

pls_residual_fit <- pls(Xmat, Y_res, ncomp = 3, mode = 'regression')
Q2.pls_residual <- perf(pls_residual_fit, validation = 'Mfold', folds = 5, nrepeat = 100)
Q2.pls_residual$measures$Q2
#$summary
#  feature comp       mean        sd
#1   L_res    1 -0.1641537 0.0896763
#2   L_res    2 -0.2684662 0.1466805
#3   L_res    3 -0.4973559 0.2554287
#4   a_res    1 -0.2996816 0.2814675
#5   a_res    2 -0.3012360 0.2500467
#6   a_res    3 -0.7515372 0.5484193
#7   b_res    1 -0.1820657 0.1157741
#8   b_res    2 -0.2740282 0.1647069
#9   b_res    3 -0.5377760 0.3085555

svg("Q2.pls_residual.svg")
plot(Q2.pls_residual, criterion = 'Q2')
dev.off()

pheno_residual_pls <- pls(Xmat, Y_res, ncomp = 1, mode = "regression")
pheno_residual_pls_perf <- perf(pheno_residual_pls, validation = 'Mfold', folds = 5, nrepeat = 100)
pheno_residual_pls_perf$measures$Q2
#$summary
#  feature comp       mean         sd
#1   L_res    1 -0.1751511 0.08912059
#2   a_res    1 -0.2953531 0.22558956
#3   b_res    1 -0.1845460 0.11383084

# summary
#summary(pheno_residual_pls) # error
pheno_residual_pls$prop_expl_var$X   # X explained variance (comp1); 0.8605407
pheno_residual_pls$prop_expl_var$Y   # Y explained variance (comp1); 0.6584153
head(pheno_residual_pls$variates$X)  # genotype latent score (t)
head(pheno_residual_pls$variates$Y)  # phenotype latent score (u)
pheno_residual_pls$variates$X[,1] # = genotype PLS1 score
pheno_residual_pls$variates$Y[,1] # = phenotype PLS1 score

pheno_residual_pls$loadings$X # genotype's contribution; genoPC1 -0.1447108, genoPC2 -0.9585190, genoPC3  0.2455607
pheno_residual_pls$loadings$Y # phenotype's contribution; L_res -0.5939927, a_res  0.6527572, b_res -0.4701921

#plotIndiv(pheno_residual_pls, comp = 1)
#plotVar(pheno_residual_pls, comp = 1)
plotLoadings(pheno_residual_pls, comp = 1, contrib = "max")

t1 <- pheno_residual_pls$variates$X[,1]   # genotype score
u1 <- pheno_residual_pls$variates$Y[,1]   # phenotype score

svg("residual_pls_corr.svg")
plot(t1, u1, xlab="PLS1 score (X: genotype)", ylab="PLS1 score (Y: phenotype)")
abline(lm(u1 ~ t1), lty=2)
dev.off()

cor(t1, u1, use="complete.obs") # 0.2400144
cor(t1, Y_res[, "L_res"]) # -0.197829
cor(t1, Y_res[, "a_res"]) # 0.2174005
cor(t1, Y_res[, "b_res"]) # -0.1565973

# permutation (optional)
set.seed(123)

obs_cor <- cor(
  pheno_residual_pls$variates$X[,1],
  pheno_residual_pls$variates$Y[,1]
)

nperm <- 999
perm_cor <- numeric(nperm)

for (i in 1:nperm) {
  Yperm <- Y_res[sample(nrow(Y_res)), ]  # phenotype permute
  pls_residual_perm <- pls(Xmat, Yperm, ncomp = 1, mode = "regression")
  perm_cor[i] <- cor(pls_residual_perm$variates$X[,1],
                     pls_residual_perm$variates$Y[,1])
}

p_perm <- mean(abs(perm_cor) >= abs(obs_cor))
p_perm # 0.6886887

# residual model - str
L_res2 <- residuals(lm(avg.L ~ strPC1 + strPC2, data=Y_COV))
a_res2 <- residuals(lm(avg.a ~ strPC1 + strPC2, data=Y_COV))
b_res2 <- residuals(lm(avg.b ~ strPC1 + strPC2, data=Y_COV))

Y_res2 <- cbind(L_res2, a_res2, b_res2)
rownames(Y_res2) <- rownames(Y_COV)

pls_residual2_fit <- pls(Xmat, Y_res2, ncomp = 3, mode = 'regression')
Q2.pls_residual2 <- perf(pls_residual2_fit, validation = 'Mfold', folds = 5, nrepeat = 100)
Q2.pls_residual2$measures$Q2
#$summary
#  feature comp       mean         sd
#1  L_res2    1 -0.1734770 0.09800471
#2  L_res2    2 -0.2553314 0.16277747
#3  L_res2    3 -0.5008410 0.26043612
#4  a_res2    1 -0.2552016 0.16350350
#5  a_res2    2 -0.3324968 0.19745971
#6  a_res2    3 -0.6501455 0.36402184
#7  b_res2    1 -0.1843096 0.14473286
#8  b_res2    2 -0.2921409 0.17421414
#9  b_res2    3 -0.5462792 0.29962237

svg("Q2.pls_residual2.svg")
plot(Q2.pls_residual2, criterion = 'Q2')
dev.off()

pheno_residual2_pls <- pls(Xmat, Y_res2, ncomp = 1, mode = "regression")
pheno_residual2_pls_perf <- perf(pheno_residual2_pls, validation = 'Mfold', folds = 5, nrepeat = 100)
pheno_residual2_pls_perf$measures$Q2
#$summary
#  feature comp       mean        sd
#1  L_res2    1 -0.1694605 0.1052035
#2  a_res2    1 -0.2360956 0.1478834
#3  b_res2    1 -0.1806802 0.1097305

# summary
pheno_residual2_pls$prop_expl_var$X   # X explained variance (comp1); 0.8709373
pheno_residual2_pls$prop_expl_var$Y   # Y explained variance (comp1); 0.6305261
head(pheno_residual2_pls$variates$X)  # genotype latent score (t)
head(pheno_residual2_pls$variates$Y)  # phenotype latent score (u)
pheno_residual2_pls$variates$X[,1] # = genotype PLS1 score
pheno_residual2_pls$variates$Y[,1] # = phenotype PLS1 score

pheno_residual2_pls$loadings$X # genotype's contribution; genoPC1 -0.1704718, genoPC2 -0.9387417, genoPC3  0.2995054
pheno_residual2_pls$loadings$Y # phenotype's contribution; L_res2 -0.6630313, a_res2  0.5868384, b_res2 -0.4647690

plotLoadings(pheno_residual2_pls, comp = 1, contrib = "max")

t1 <- pheno_residual2_pls$variates$X[,1]   # genotype score
u1 <- pheno_residual2_pls$variates$Y[,1]   # phenotype score

svg("residual2_pls_corr.svg")
plot(t1, u1, xlab="PLS1 score (X: genotype)", ylab="PLS1 score (Y: phenotype)")
abline(lm(u1 ~ t1), lty=2)
dev.off()

cor(t1, u1, use="complete.obs") # 0.2387201
cor(t1, Y_res2[, "L_res2"]) # -0.2140298
cor(t1, Y_res2[, "a_res2"]) # 0.1894344
cor(t1, Y_res2[, "b_res2"]) # -0.1500298

# permutation (optional)
set.seed(123)

obs_cor <- cor(
  pheno_residual2_pls$variates$X[,1],
  pheno_residual2_pls$variates$Y[,1]
)

nperm <- 999
perm_cor <- numeric(nperm)

for (i in 1:nperm) {
  Yperm <- Y_res2[sample(nrow(Y_res2)), ]  # phenotype permute
  pls_residual2_perm <- pls(Xmat, Yperm, ncomp = 1, mode = "regression")
  perm_cor[i] <- cor(pls_residual2_perm$variates$X[,1],
                     pls_residual2_perm$variates$Y[,1])
}

p_perm <- mean(abs(perm_cor) >= abs(obs_cor))
p_perm # 0.7067067

# residual model - bio
L_res3 <- residuals(lm(avg.L ~ bioPC1 + bioPC2, data=Y_COV))
a_res3 <- residuals(lm(avg.a ~ bioPC1 + bioPC2, data=Y_COV))
b_res3 <- residuals(lm(avg.b ~ bioPC1 + bioPC2, data=Y_COV))

Y_res3 <- cbind(L_res3, a_res3, b_res3)
rownames(Y_res3) <- rownames(Y_COV)

pls_residual3_fit <- pls(Xmat, Y_res3, ncomp = 3, mode = 'regression')
Q2.pls_residual3 <- perf(pls_residual3_fit, validation = 'Mfold', folds = 5, nrepeat = 100)
Q2.pls_residual3$measures$Q2
#$summary
#  feature comp        mean        sd
#1  L_res3    1 -0.06621388 0.1214581
#2  L_res3    2 -0.52536599 1.5784420
#3  L_res3    3 -0.60023508 1.0904253
#4  a_res3    1 -0.22953286 0.4249979
#5  a_res3    2 -0.52319804 0.5982791
#6  a_res3    3 -0.86014028 1.1618897
#7  b_res3    1 -0.16286068 0.3111791
#8  b_res3    2 -0.43171105 0.7225713
#9  b_res3    3 -0.74934078 1.4643028

svg("Q2.pls_residual3.svg")
plot(Q2.pls_residual3, criterion = 'Q2')
dev.off()

pheno_residual3_pls <- pls(Xmat, Y_res3, ncomp = 1, mode = "regression")
pheno_residual3_pls_perf <- perf(pheno_residual3_pls, validation = 'Mfold', folds = 5, nrepeat = 100)
pheno_residual3_pls_perf$measures$Q2
#$summary
#  feature comp        mean         sd
#1  L_res3    1 -0.07198936 0.07382271
#2  a_res3    1 -0.16779300 0.17715736
#3  b_res3    1 -0.10412108 0.11227407

# summary
pheno_residual3_pls$prop_expl_var$X   # X explained variance (comp1); 0.9081394
pheno_residual3_pls$prop_expl_var$Y   # Y explained variance (comp1); 0.6808874
head(pheno_residual3_pls$variates$X)  # genotype latent score (t)
head(pheno_residual3_pls$variates$Y)  # phenotype latent score (u)
pheno_residual3_pls$variates$X[,1] # = genotype PLS1 score
pheno_residual3_pls$variates$Y[,1] # = phenotype PLS1 score

pheno_residual3_pls$loadings$X # genotype's contribution; genoPC1 -0.4692666, genoPC2 -0.7215115, genoPC3  0.5091268
pheno_residual3_pls$loadings$Y # phenotype's contribution; L_res2 -0.6141025, a_res2  0.6106372, b_res2 -0.5000003

plotLoadings(pheno_residual3_pls, comp = 1, contrib = "max")

t1 <- pheno_residual3_pls$variates$X[,1]   # genotype score
u1 <- pheno_residual3_pls$variates$Y[,1]   # phenotype score

svg("residual3_pls_corr.svg")
plot(t1, u1, xlab="PLS1 score (X: genotype)", ylab="PLS1 score (Y: phenotype)")
abline(lm(u1 ~ t1), lty=2)
dev.off()

cor(t1, u1, use="complete.obs") # 0.3403544
cor(t1, Y_res3[, "L_res3"]) # -0.2960478
cor(t1, Y_res3[, "a_res3"]) # 0.2943772
cor(t1, Y_res3[, "b_res3"]) # -0.2410411

# permutation (optional)
set.seed(123)

obs_cor <- cor(
  pheno_residual3_pls$variates$X[,1],
  pheno_residual3_pls$variates$Y[,1]
)

nperm <- 999
perm_cor <- numeric(nperm)

for (i in 1:nperm) {
  Yperm <- Y_res3[sample(nrow(Y_res3)), ]  # phenotype permute
  pls_residual3_perm <- pls(Xmat, Yperm, ncomp = 1, mode = "regression")
  perm_cor[i] <- cor(pls_residual3_perm$variates$X[,1],
                     pls_residual3_perm$variates$Y[,1])
}

p_perm <- mean(abs(perm_cor) >= abs(obs_cor))
p_perm # 0.2522523


# Further inspection

genoPC1_loadings <- Geno_outlier_pca$rotation[, 1]
pc1_loadings <- data.frame(variant = names(genoPC1_loadings), loading = genoPC1_loadings)
pc1_loadings <- pc1_loadings %>% mutate(abs_loading = abs(loading)) %>% arrange(desc(abs_loading)) %>% mutate(rank = row_number())

pc1_loadings[grepl("chr1", pc1_loadings$variant), ]
pc1_loadings["chr1_59728592_2", ] # chr1_59728592_2 chr1_59728592_2 0.04854993  0.04854993  199
pc1_loadings["chr1_68382845_1", ] # chr1_68382845_1 chr1_68382845_1  0.01478846  0.01478846  467
pc1_loadings["chr1_74973144_2", ] # chr1_74973144_2 chr1_74973144_2 -0.03301838  0.03301838  337
pc1_loadings[grepl("chr2", pc1_loadings$variant), ]
pc1_loadings["chr2_96188458_1", ] # chr2_96188458_1 chr2_96188458_1 0.03920401  0.03920401  303
pc1_loadings["chr2_96188458_2", ] # chr2_96188458_2 chr2_96188458_2 0.03920401  0.03920401  304
pc1_loadings[grepl("chr3", pc1_loadings$variant), ]
pc1_loadings["chr3_285215_1", ] # chr3_285215_1 chr3_285215_1 -0.05866289  0.05866289   73
pc1_loadings["chr3_51891085_1", ] # chr3_51891085_1   chr3_51891085_1  0.063185758 0.063185758   40
pc1_loadings[grepl("chr4", pc1_loadings$variant), ]
pc1_loadings["chr4_37517378_1", ] # chr4_37517378_1 chr4_37517378_1 0.03734278  0.03734278  317
pc1_loadings["chr4_62626603_1", ] # chr4_62626603_1 chr4_62626603_1 -0.05480331  0.05480331  117
pc1_loadings[grepl("chr5", pc1_loadings$variant), ]
pc1_loadings["chr5_16029134_1", ] # chr5_16029134_1 chr5_16029134_1 0.01327634  0.01327634  484
pc1_loadings[grepl("chr6", pc1_loadings$variant), ]
pc1_loadings["chr6_12744614_1", ] # chr6_12744614_1 chr6_12744614_1 0.01163656  0.01163656  496
pc1_loadings["chr6_12744614_2", ] # chr6_12744614_2 chr6_12744614_2 0.01163656  0.01163656  497
pc1_loadings[grepl("chr7", pc1_loadings$variant), ]
pc1_loadings["chr7_12391327_1", ] # chr7_12391327_1 chr7_12391327_1  0.06636067  0.06636067   25
pc1_loadings["chr7_19896871_1", ] # chr7_19896871_1 chr7_19896871_1  0.05471959  0.05471959  119
pc1_loadings[grepl("chr10", pc1_loadings$variant), ]
pc1_loadings["chr10_9831264_1", ] # chr10_9831264_1 chr10_9831264_1 0.03583839  0.03583839  327
pc1_loadings[grepl("chr12", pc1_loadings$variant), ]
pc1_loadings["chr12_9730381_1", ] # chr12_9730381_1 chr12_9730381_1 0.01950704  0.01950704  426
pc1_loadings[grepl("chr20", pc1_loadings$variant), ]
pc1_loadings["chr20_3117406_1", ] # chr20_3117406_1 chr20_3117406_1 0.02267192  0.02267192  407

scores(pheno_rda_fit, choices = 1, display = "bp") # predictor loadings on RDA1; genoPC1  0.9477311, genoPC2  0.9221732, genoPC3 -0.9786687
scores(pheno_rda_fit, choices = 1, display = "sp") # response loadings on RDA1; avg.L  1.1018312, avg.a -0.4884759, avg.b  0.8367664