#module load proj
#module load geos
#module load gdal
#module load r

R
library(terra, lib.loc = "/home/jeon96/R/negishi/4.2.2-gcc-12.2.0-524zdta")
set.seed(123)

# Load raster (should be 259 x 271)
r <- rast("/scratch/negishi/jeon96/monq/slim/hsm/monq_hsm.tif")

# Get raster dimensions
nrow_r <- nrow(r)
ncol_r <- ncol(r)

# Extract values and find valid cells
values_all <- values(r)
valid_cells <- which(values_all != 0)
valid_values <- values_all[valid_cells]

# Normalize for probability sampling
prob_weights <- valid_values / sum(valid_values)

# Number of samples
n_samples <- 11200 # 112000

# Sample 1D cell indices
sampled_indices <- sample(valid_cells, size = n_samples, replace = TRUE, prob = prob_weights)

# Convert to (row, col)
rows <- rowFromCell(r, sampled_indices)
cols <- colFromCell(r, sampled_indices)

# Add uniform random jitter in [0,1) to each row/col
rows_jittered <- rows - runif(n_samples)
rows_flipped <- nrow_r - rows_jittered #In raster maps (like terra::rast), row 1 is at the top, and row indices increase downward. In Cartesian coordinate systems (like in SLiM and most standard XY plots), y = 0 is at the bottom, and values increase upward.
cols_jittered <- cols - runif(n_samples)

# Combine into data frame
sampled_locations <- data.frame(col_float = cols_jittered, row_float = rows_flipped)

# Save to CSV
write.table(sampled_locations,"monq_prelocations_11200.csv", sep=',', row.names = FALSE, col.names = FALSE)