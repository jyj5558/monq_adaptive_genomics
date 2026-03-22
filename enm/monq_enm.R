library(ENMeval) 
library(dismo) 
library(raster)
library(parallel)
options(cores = detectCores() - 5) 
library(dplyr)
library(biomod2)
library(lattice)
library(tidyterra)
library(ggtext)

# Run in Negishi
#library(terra, lib.loc = "/home/jeon96/R/negishi/4.2.2-gcc-12.2.0-524zdta")
#library(ENMeval, lib.loc = "/home/jeon96/R/negishi/4.2.2-gcc-12.2.0-524zdta")
#library(sp, lib.loc = "/home/jeon96/R/negishi/4.2.2-gcc-12.2.0-524zdta")
#library(raster, lib.loc = "/home/jeon96/R/negishi/4.2.2-gcc-12.2.0-524zdta")
#library(dismo, lib.loc = "/home/jeon96/R/negishi/4.2.2-gcc-12.2.0-524zdta") 
#library(parallel)
#options(cores = detectCores() - 5) 
#library(tibble, lib.loc = "/home/jeon96/R/negishi/4.2.2-gcc-12.2.0-524zdta")
#library(dplyr)
#library(ggplot2)
#library(lattice, lib.loc = "/home/jeon96/R/negishi/4.2.2-gcc-12.2.0-524zdta")
#library(tidyterra, lib.loc = "/home/jeon96/R/negishi/4.2.2-gcc-12.2.0-524zdta")
#library(class)
#library(nlme, lib.loc = "/home/jeon96/R/negishi/4.2.2-gcc-12.2.0-524zdta")
#library(caret, lib.loc = "/home/jeon96/R/negishi/4.2.2-gcc-12.2.0-524zdta")
#library(randomForest, lib.loc = "/home/jeon96/R/negishi/4.2.2-gcc-12.2.0-524zdta")
#library(mda, lib.loc = "/home/jeon96/R/negishi/4.2.2-gcc-12.2.0-524zdta")
#library(gam, lib.loc = "/home/jeon96/R/negishi/4.2.2-gcc-12.2.0-524zdta")
#library(gbm, lib.loc = "/home/jeon96/R/negishi/4.2.2-gcc-12.2.0-524zdta")
#library(nnet, lib.loc = "/home/jeon96/R/negishi/4.2.2-gcc-12.2.0-524zdta")
#library(Formula, lib.loc = "/home/jeon96/R/negishi/4.2.2-gcc-12.2.0-524zdta")
#library(plotrix, lib.loc = "/home/jeon96/R/negishi/4.2.2-gcc-12.2.0-524zdta")
#library(plotmo, lib.loc = "/home/jeon96/R/negishi/4.2.2-gcc-12.2.0-524zdta")
#library(earth, lib.loc = "/home/jeon96/R/negishi/4.2.2-gcc-12.2.0-524zdta")
#library(xgboost, lib.loc = "/home/jeon96/R/negishi/4.2.2-gcc-12.2.0-524zdta")
#library(rpart, lib.loc = "/home/jeon96/R/negishi/4.2.2-gcc-12.2.0-524zdta")
#library(biomod2, lib.loc = "/home/jeon96/R/negishi/4.2.2-gcc-12.2.0-524zdta")
#library(ggtext, lib.loc = "/home/jeon96/R/negishi/4.2.2-gcc-12.2.0-524zdta")
#library(maxnet, lib.loc = "/home/jeon96/R/negishi/4.2.2-gcc-12.2.0-524zdta")
#detectCores()

setwd("D:/Research_Data_Backup/Montezuma_quail/enm")
#setwd("/scratch/negishi/jeon96/monq/enm")

# Read data
monq_envs.files <- list.files(path = "bioclim/present/processed", pattern = ".tif$", full.names=TRUE)
#monq_envs.files <- list.files(path = "present_processed", pattern = ".tif$", full.names=TRUE)
monq_envs <- stack(monq_envs.files)
proj_wgs84 <- crs(monq_envs)

monq_points <- read.csv("monq_enm_gis/monq_rarefied_points.csv", header = T) # there are few weirdly formated rows. manually check and correct them.
#monq_points <- read.csv("monq_rarefied_points.csv", header = T)
monq_coordinates <- cbind(as.numeric(monq_points$x), as.numeric(monq_points$y))
colnames(monq_coordinates ) <- c("Lon", "Lat")
monq_occs <- as.data.frame(monq_coordinates)
monq_occs.z <- cbind(monq_occs, raster::extract(monq_envs, monq_occs)) # extract raster values at points

monq_bias <- raster("monq_enm_gis/monq_sampling_bias.asc")
#monq_bias <- raster("monq_sampling_bias.asc")

# Plot first raster in the stack, the mean annual temperature
plot(monq_bias)

# Add points for all the occurrence points onto the raster
points(monq_occs, col = "red")

# Randomly sample 10,000 background points from one background extent raster (n=10000 is recommended for MaxEnt by Barbet-Massin et al. 2012 in Methods in Ecology and Evolution)
monq_bgpoints <- dismo::randomPoints(monq_bias, n = 10000, prob = T) %>% as.data.frame()
#saveRDS(monq_bgpoints, "monq_bg.rds")
#monq_bgpoints <- readRDS("monq_bg.rds")
colnames(monq_bgpoints) <- colnames(monq_coordinates)
monq_bg.z <- cbind(monq_bgpoints, raster::extract(monq_envs, monq_bgpoints)) # extract raster values at points

plot(monq_bias)
points(monq_bgpoints, pch = 20, cex = 0.2)

# Run ENMeval
e.mx.monq <- ENMevaluate(occs = monq_occs.z, bg = monq_bg.z, 
                        algorithm = 'maxent.jar', partitions = 'block', doClamp = TRUE,
                        tune.args = list(fc = c("L","LQ","LQH","LQHP","H"), rm = seq(0.5,5,0.5)), 
                        parallel = TRUE, numCores = 4)
e.mx.monq
str(e.mx.monq, max.level=2)

# Results table with summary statistics for cross validation on test data.
eval.results(e.mx.monq) %>% head()

# Results table with cross validation statistics for each test partition.
eval.results.partitions(e.mx.monq) %>% head()

# Visualize tuning results
evalplot.stats(e = e.mx.monq, stats = c("or.mtp", "auc.val"), color = "fc", x.var = "rm", error.bars = FALSE)

# Overall results
res.monq <- eval.results(e.mx.monq)

# Select the model with the lowest AICc score
res.monq[order(res.monq$delta.AICc),] %>% head()
opt.aicc.monq <- res.monq %>% filter(delta.AICc == 0)
opt.aicc.monq # fc: H, rm: 1

# the sequential criteria
res.monq[order(res.monq$or.10p.avg, -res.monq$auc.val.avg),] %>% head()
opt.seq.monq <- res.monq %>% 
  filter(or.10p.avg == min(or.10p.avg)) %>% 
  filter(auc.val.avg == max(auc.val.avg))
opt.seq.monq # fc: LQ, rm: 0.5

# result model
monq_model <- e.mx.monq@models$fc.H_rm.1

svg("monq_responseCurves.svg", width = 10, height = 10) # response curve
response(monq_model)
dev.off()

svg("monq_variableContribution.svg", width = 7, height = 7) # variable contribution
plot(monq_model)
dev.off()

# model evaluation
var.import <- eval.variable.importance(e.mx.monq) # variable importance (variable contribution, permutation importance)
monq_var.import <- var.import$fc.H_rm.1
monq_var.import
write.csv(monq_var.import, "var_importance.csv")

monq_eval <- eval.models(e.mx.monq)
monq_eval <- monq_eval$fc.H_rm.1
monq_eval@results # auc training = 0.7941 (0.8486143 from e.mx.monq@results)
e.mx.monq@results # auc validation average = 0.8311749

# Dismo
#pred.present <- dismo::predict(monq_envs, monq_model)
#pred.present # values     : 2.031947e-05, 0.9317202  (min, max)
#plot(pred.present)

#terra::writeRaster(pred.present, "monq_present.tif", overwrite=TRUE, NAflag=0)
# then in ArcGIS Pro, convert this to 8-bit unsigned png file with checking "Scale Pixel Value"

# 2011-2040 model
#monq_envs_2040.files <- list.files(path = "bioclim/future/processed/2011-2040", pattern = ".tif$", full.names=TRUE)
#monq_envs_2040 <- stack(monq_envs_2040.files)
#proj_wgs84 <- crs(monq_envs_2040)

#pred.2040 <- dismo::predict(monq_envs_2040, monq_model)


# Biomod
monq_pre <- monq_occs
monq_pre$presence <- 1
monq_abs <- monq_bgpoints
monq_abs$presence <- NA
monq_res <- rbind(monq_pre, monq_abs)

# Format the data
monq_data <- BIOMOD_FormatingData(
  resp.var = monq_res['presence'],
  resp.xy = monq_res[, c('Lon', 'Lat')],
  expl.var = monq_envs,
  resp.name = "monq",
  PA.nb.rep = 1,
  PA.nb.absences = 10000,
  PA.strategy = 'random'
)

plot(monq_data)
summary(monq_data)
head(monq_data@PA.table)

# Run model
monq_model <- BIOMOD_Modeling(
  bm.format = monq_data,
  models = c("MAXENT", "RF", "GLM", "ANN", "FDA", "GAM", "GBM", "MARS", "SRE", "XGBOOST", "CTA"),
  OPT.strategy = 'tuned',
  CV.strategy = 'block',
  #CV.nb.rep = 3,
  CV.perc = 0.7,
  var.import = 3,
  #nb.cpu = 1,
  metric.eval = c('TSS','ROC'),
  modeling.id = "monq_sdm"
)

# Run model algo-by-algo (in Negishi cluster)
monq_cta <- BIOMOD_Modeling(
  bm.format = monq_data,
  models = c("CTA"),
  OPT.strategy = 'tuned',
  CV.strategy = 'random',
  CV.nb.rep = 3,
  CV.perc = 0.7,
  nb.cpu = 120,
  var.import = 3,
  metric.eval = c('TSS','ROC'),
  modeling.id = "monq_cta",
  seed.val = 123
)

monq_cta_scores <- get_evaluations(monq_cta) # Get model evaluation scores
dim(monq_cta_scores)
dimnames(monq_cta_scores)
(monq_cta_eval.scor_mean <- aggregate(data = monq_cta_scores, calibration ~ metric.eval, FUN = mean))
  # TSS=0.6016, ROC=0.8394
  
monq_xgboost <- BIOMOD_Modeling(
  bm.format = monq_data,
  models = c("XGBOOST"),
  OPT.strategy = 'tuned',
  CV.strategy = 'random',
  CV.nb.rep = 3,
  CV.perc = 0.7,
  nb.cpu = 120,
  var.import = 3,
  metric.eval = c('TSS','ROC'),
  modeling.id = "monq_xgboost",
  seed.val = 123
)

monq_xgboost_scores <- get_evaluations(monq_xgboost) # Get model evaluation scores
dim(monq_xgboost_scores)
dimnames(monq_xgboost_scores)
(monq_xgboost_eval.scor_mean <- aggregate(data = monq_xgboost_scores, calibration ~ metric.eval, FUN = mean))
  # TSS=0.5770, ROC=0.8486

monq_sre <- BIOMOD_Modeling(
  bm.format = monq_data,
  models = c("SRE"),
  OPT.strategy = 'tuned',
  CV.strategy = 'random',
  CV.nb.rep = 3,
  CV.perc = 0.7,
  nb.cpu = 120,
  var.import = 3,
  metric.eval = c('TSS','ROC'),
  modeling.id = "monq_sre",
  seed.val = 123
)

monq_sre_scores <- get_evaluations(monq_sre) # Get model evaluation scores
dim(monq_sre_scores)
dimnames(monq_sre_scores)
(monq_sre_eval.scor_mean <- aggregate(data = monq_sre_scores, calibration ~ metric.eval, FUN = mean))
  # TSS=0.4028, ROC=0.7012

monq_mars <- BIOMOD_Modeling(
  bm.format = monq_data,
  models = c("MARS"),
  OPT.strategy = 'tuned',
  CV.strategy = 'random',
  CV.nb.rep = 3,
  CV.perc = 0.7,
  nb.cpu = 120,
  var.import = 3,
  metric.eval = c('TSS','ROC'),
  modeling.id = "monq_mars",
  seed.val = 123
)

monq_mars_scores <- get_evaluations(monq_mars) # Get model evaluation scores
dim(monq_mars_scores)
dimnames(monq_mars_scores)
(monq_mars_eval.scor_mean <- aggregate(data = monq_mars_scores, calibration ~ metric.eval, FUN = mean))
  # TSS=0.57625, ROC=0.84550

monq_gbm <- BIOMOD_Modeling(
  bm.format = monq_data,
  models = c("GBM"),
  OPT.strategy = 'tuned',
  CV.strategy = 'random',
  CV.nb.rep = 3,
  CV.perc = 0.7,
  nb.cpu = 120,
  var.import = 3,
  metric.eval = c('TSS','ROC'),
  modeling.id = "monq_gbm",
  seed.val = 123
)

monq_gbm_scores <- get_evaluations(monq_gbm) # Get model evaluation scores
dim(monq_gbm_scores)
dimnames(monq_gbm_scores)
(monq_gbm_eval.scor_mean <- aggregate(data = monq_gbm_scores, calibration ~ metric.eval, FUN = mean))
  # TSS=0.8726, ROC=0.9778

monq_gam <- BIOMOD_Modeling(
  bm.format = monq_data,
  models = c("GAM"),
  OPT.strategy = 'tuned',
  CV.strategy = 'random',
  CV.nb.rep = 3,
  CV.perc = 0.7,
  nb.cpu = 120,
  var.import = 3,
  metric.eval = c('TSS','ROC'),
  modeling.id = "monq_gam",
  seed.val = 123
)

monq_gam_scores <- get_evaluations(monq_gam) # Get model evaluation scores
dim(monq_gam_scores)
dimnames(monq_gam_scores)
(monq_gam_eval.scor_mean <- aggregate(data = monq_gam_scores, calibration ~ metric.eval, FUN = mean))
  # TSS=0.5164, ROC=0.8130

monq_fda <- BIOMOD_Modeling(
  bm.format = monq_data,
  models = c("FDA"),
  OPT.strategy = 'tuned',
  CV.strategy = 'random',
  CV.nb.rep = 3,
  CV.perc = 0.7,
  nb.cpu = 120,
  var.import = 3,
  metric.eval = c('TSS','ROC'),
  modeling.id = "monq_fda",
  seed.val = 123
)

monq_fda_scores <- get_evaluations(monq_fda) # Get model evaluation scores
dim(monq_fda_scores)
dimnames(monq_fda_scores)
(monq_fda_eval.scor_mean <- aggregate(data = monq_fda_scores, calibration ~ metric.eval, FUN = mean))
  # TSS=0.5814, ROC=0.8432

monq_ann <- BIOMOD_Modeling(
  bm.format = monq_data,
  models = c("ANN"),
  OPT.strategy = 'bigboss',
  CV.strategy = 'random',
  CV.nb.rep = 3,
  CV.perc = 0.7,
  nb.cpu = 120,
  var.import = 3,
  metric.eval = c('TSS','ROC'),
  modeling.id = "monq_ann",
  seed.val = 123
)

monq_ann_scores <- get_evaluations(monq_ann) # Get model evaluation scores
dim(monq_ann_scores)
dimnames(monq_ann_scores)
(monq_ann_eval.scor_mean <- aggregate(data = monq_ann_scores, calibration ~ metric.eval, FUN = mean))
  # TSS=0.5682, ROC=0.8288

monq_glm <- BIOMOD_Modeling(
  bm.format = monq_data,
  models = c("GLM"),
  OPT.strategy = 'tuned',
  CV.strategy = 'random',
  CV.nb.rep = 3,
  CV.perc = 0.7,
  nb.cpu = 120,
  var.import = 3,
  metric.eval = c('TSS','ROC'),
  modeling.id = "monq_glm",
  seed.val = 123
)

monq_glm_scores <- get_evaluations(monq_glm) # Get model evaluation scores
dim(monq_glm_scores)
dimnames(monq_glm_scores)
(monq_glm_eval.scor_mean <- aggregate(data = monq_glm_scores, calibration ~ metric.eval, FUN = mean))
  # TSS=0.5870, ROC=0.8538

monq_rf <- BIOMOD_Modeling(
  bm.format = monq_data,
  models = c("RF"),
  OPT.strategy = 'tuned',
  CV.strategy = 'random',
  CV.nb.rep = 3,
  CV.perc = 0.7,
  nb.cpu = 120,
  var.import = 3,
  metric.eval = c('TSS','ROC'),
  modeling.id = "monq_rf",
  seed.val = 123
)
#saveRDS(monq_rf, file = "monq_rf.rds")
monq_rf_scores <- get_evaluations(monq_rf) # Get model evaluation scores
dim(monq_rf_scores)
dimnames(monq_rf_scores)
(monq_rf_eval.scor_mean <- aggregate(data = monq_rf_scores, calibration ~ metric.eval, FUN = mean))
  # TSS=0.99, ROC=1.00

monq_max <- BIOMOD_Modeling(
  bm.format = monq_data,
  models = c("MAXNET"),
  OPT.strategy = 'tuned',
  CV.strategy = 'random',
  CV.nb.rep = 3,
  CV.perc = 0.7,
  nb.cpu = 120,
  var.import = 3,
  metric.eval = c('TSS','ROC'),
  modeling.id = "monq_max",
  seed.val = 123
)

monq_max_scores <- get_evaluations(monq_max) # Get model evaluation scores
dim(monq_max_scores)
dimnames(monq_max_scores)
(monq_max_eval.scor_mean <- aggregate(data = monq_max_scores, calibration ~ metric.eval, FUN = mean))
  # TSS=0.5728, ROC=0.8504
## The best model = random forest

# Plot model evaluation scores
bm_PlotEvalMean(monq_rf)

# Check variable importance
(monq_rf_var_import <- get_variables_importance(monq_rf))

# Make the mean of variable importance by algorithm
(monq_rf_var.imp_mean <- aggregate(data = monq_rf_var_import, var.imp ~ expl.var, FUN = mean))
  # aspect=0.0664688, slope=0.0656562, bio2=0.0721172, bio5=0.0918364, bio7=0.3383939, bio8=0.0417030, bio9=0.1459687, bio14=0.1118603, bio15=0.2460208, bio18=0.2271341

# Model response plots
monq_rf_eval_plot <- 
  bm_PlotResponseCurves(
    bm.out  = monq_rf,
    new.env = get_formal_data(monq_rf,'expl.var'), 
    show.variables= get_formal_data(monq_rf,'expl.var.names'),
    do.bivariate = FALSE,
    fixed.var = 'median',
  )

# Run the ensemble model
monq_rf_ensemble <- 
  BIOMOD_EnsembleModeling(
    bm.mod = monq_rf,
    models.chosen = 'all',
    em.by = 'all',
    em.algo = 'EMwmean',
    metric.eval = c('TSS','ROC'),
    metric.select = "TSS",
    var.import = 3,
    nb.cpu = 120,
    seed.val = 123
  )
#saveRDS(monq_rf_ensemble, file = "monq_rf_ensemble.rds")

# Assess ensemble models quality ----
(monq_rf_ensemble_scores <- get_evaluations(monq_rf_ensemble))
  # TSS=0.928, ROC=0.994

# Check variable importance
(monq_rf_ensemble_var_import <- get_variables_importance(monq_rf_ensemble))
(monq_rf_ensemble_var.imp_mean <- aggregate(data = monq_rf_ensemble_var_import, var.imp ~ expl.var, FUN = mean))
  # aspect=0.05776400, slope=0.05595733, bio2=0.06643600, bio5=0.08808467, bio7=0.35282933, bio8=0.03599067, bio9=0.13687967, bio14=0.09517300, bio15=0.25168367, bio18=0.21459833

# Model response plots
monq_rf_ensemble_eval_plot <- 
  bm_PlotResponseCurves(
    bm.out  = monq_rf_ensemble,
    new.env = get_formal_data(monq_rf_ensemble,'expl.var'), 
    show.variables= get_formal_data(monq_rf_ensemble,'expl.var.names'),
    do.bivariate = FALSE,
    fixed.var = 'median',
  )

# Present environment projections
monq_rf_proj_present <- 
  BIOMOD_Projection(
    bm.mod = monq_rf,
    models.chosen = "all",
    new.env = monq_envs,
    proj.name = "monq_present",
    metric.binary = "all",
    output.format = ".tif",
    do.stack = FALSE
  )

monq_rf_ensemble_proj_present <- 
  BIOMOD_EnsembleForecasting(
    bm.em = monq_rf_ensemble,
    models.chosen = "all",
    bm.proj = monq_rf_proj_present,
    metric.binary = "all",
    output.format = ".tif",
    do.stack = FALSE,
    nb.cpu = 120
  )
plot(monq_rf_ensemble_proj_present)

monq_rf_ensemble_present.pred <- get_predictions(monq_rf_ensemble_proj_present)
monq_rf_ensemble_present.pred.bin <- get_predictions(monq_rf_ensemble_proj_present, metric.binary = "TSS") #This is the binary representation
plot(monq_rf_ensemble_present.pred)
plot(monq_rf_ensemble_present.pred.bin)

for (year in c("2040","2070","2100")){
  for (scen in c("ssp126")) {
    envs_files.name <- paste(year,scen,"envs", sep="_")
    assign(envs_files.name, list.files(path = paste0(year,"_processed"), pattern = ".tif$", full.names=TRUE))
    envs_stack.name <- paste0(year,"_",scen,"_envs",".stack")
    assign(envs_stack.name, stack(get(envs_files.name)))
    envs_stack <- get(envs_stack.name)
    envs_stack <- setMinMax(envs_stack)
    envs_stack
    new_env_names <- c("aspect", "bio14", "bio15", "bio18", "bio2", "bio5", "bio7", "bio8", "bio9", "slope") # alphabetically ordered
    print(envs_stack)
    names(envs_stack) <- new_env_names
    print(envs_stack)
    sngl.name <- paste0("monq_rf_proj_",year,"_",scen)
    assign(sngl.name, BIOMOD_Projection(
      bm.mod = monq_rf,
      models.chosen = "all",
      new.env = envs_stack,
      proj.name = paste0("monq_",year,"_",scen),
      metric.binary = "all",
      output.format = ".tif",
      do.stack = FALSE
    ))
    esmbl.name <- paste0("monq_rf_ensemble_proj_",year,"_",scen)
    assign(esmbl.name, BIOMOD_EnsembleForecasting(
      bm.em = monq_rf_ensemble,
      models.chosen = "all",
      bm.proj = get(sngl.name),
      metric.binary = "all",
      output.format = ".tif",
      do.stack = FALSE,
      nb.cpu = 120
    ))
    pred_file.name <- paste0(paste("monq_rf_ensemble",year,scen, sep="_"),".pred")
    assign(pred_file.name, get_predictions(get(esmbl.name)))
    pred.bin_file.name <- paste0(paste("monq_rf_ensemble",year,scen, sep="_"),".pred.bin")
    assign(pred.bin_file.name, get_predictions(get(esmbl.name), metric.binary = "TSS"))
  }
} 

for (year in c("minus0","minus50","minus100","minus150","minus200")){
  for (scen in c("paleo")) {
    envs_files.name <- paste(year,scen,"envs", sep="_")
    assign(envs_files.name, list.files(path = paste0(year,"_processed"), pattern = ".tif$", full.names=TRUE))
    envs_stack.name <- paste0(year,"_",scen,"_envs",".stack")
    assign(envs_stack.name, stack(get(envs_files.name)))
    envs_stack <- get(envs_stack.name)
    envs_stack <- setMinMax(envs_stack)
    envs_stack
    new_env_names <- c("aspect", "bio14", "bio15", "bio18", "bio2", "bio5", "bio7", "bio8", "bio9", "slope") # alphabetically ordered
    print(envs_stack)
    names(envs_stack) <- new_env_names
    print(envs_stack)
    sngl.name <- paste0("monq_rf_proj_",year,"_",scen)
    assign(sngl.name, BIOMOD_Projection(
      bm.mod = monq_rf,
      models.chosen = "all",
      new.env = envs_stack,
      proj.name = paste0("monq_",year,"_",scen),
      metric.binary = "all",
      output.format = ".tif",
      do.stack = FALSE
    ))
    esmbl.name <- paste0("monq_rf_ensemble_proj_",year,"_",scen)
    assign(esmbl.name, BIOMOD_EnsembleForecasting(
      bm.em = monq_rf_ensemble,
      models.chosen = "all",
      bm.proj = get(sngl.name),
      metric.binary = "all",
      output.format = ".tif",
      do.stack = FALSE,
      nb.cpu = 120
    ))
    pred_file.name <- paste0(paste("monq_rf_ensemble",year,scen, sep="_"),".pred")
    assign(pred_file.name, get_predictions(get(esmbl.name)))
    pred.bin_file.name <- paste0(paste("monq_rf_ensemble",year,scen, sep="_"),".pred.bin")
    assign(pred.bin_file.name, get_predictions(get(esmbl.name), metric.binary = "TSS"))
  }
} 

# Plotting -> in ArcGIS
custom_palette <- colorRampPalette(c("#FFD966", "#00A600"))
zlim <- c(0, 1)

monq_rf_ensemble_present.pred <- raster("./monq/proj_monq_present/individual_projections/monq_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo.tif")
raster_outline <- rasterToPolygons(monq_rf_ensemble_present.pred, dissolve = TRUE)

tiff("./monq/monq_present.tiff", width=664, height=664)
plot(monq_rf_ensemble_present.pred, col = "black", legend=FALSE, asp=1)
plot(monq_rf_ensemble_present.pred, col = custom_palette(100), breaks = seq(zlim[1], zlim[2], length.out = 100), zlim = zlim, legend=FALSE, add=TRUE)
dev.off()


# Interpolating -> in ArcGIS
library(raster)
library(gstat)
library(sp)
library(sf)

setwd("D:/Research_Data_Backup/Montezuma_quail/enm/monq_enm_R")
#suit_present <- raster("monq/models/present/individual_projections/monq_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo.tif")
#suit_present <- raster("monq/models/present/suit_present_merged_ext.tif")
extent(suit_present) 
#xmin       : -111.6835 
#xmax       : -97.98347 
#ymin       : 18.80819 
#ymax       : 34.89986
#suit_2040 <- raster("monq/models/2011-2040/individual_projections/monq_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo.tif")
#suit_2040 <- raster("monq/models/2011-2040/suit_2040_merged.tif")
extent(suit_2040)
#xmin       : -113.2168 
#xmax       : -96.60014 
#ymin       : 17.49986 
#ymax       : 36.20819  
#suit_2070 <- raster("monq/models/2041-2070/individual_projections/monq_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo.tif")
#suit_2070 <- raster("monq/models/2041-2070/suit_2070_merged.tif")
extent(suit_2070)
#xmin       : -114.8251 
#xmax       : -95.25847 
#ymin       : 16.34986 
#ymax       : 37.55819 
#suit_2100 <- raster("monq/models/2071-2100/individual_projections/monq_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo.tif")
#suit_2100 <- raster("monq/models/2071-2100/suit_2100_merged.tif")
extent(suit_2100)
#xmin       : -116.4668 
#xmax       : -93.90014 
#ymin       : 15.65819 
#ymax       : 38.91653 
#suit_minus0 <- raster("monq/models/minus0/individual_projections/monq_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo.tif")
extent(suit_minus0)
#xmin       : -117.0085 
#xmax       : -93.37514 
#ymin       : 15.63319 
#ymax       : 39.36653
#suit_minus50 <- raster("monq/models/minus50/individual_projections/monq_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo.tif")
extent(suit_minus50)
#xmin       : -117.0085 
#xmax       : -93.37514 
#ymin       : 15.63319 
#ymax       : 39.36653
#suit_minus100 <- raster("monq/models/minus100/individual_projections/monq_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo.tif")
extent(suit_minus100)
#xmin       : -117.0085 
#xmax       : -93.31681 
#ymin       : 15.63319 
#ymax       : 39.36653 
#suit_minus150 <- raster("monq/models/minus150/individual_projections/monq_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo.tif")
extent(suit_minus150)
#xmin       : -117.0085 
#xmax       : -93.30847 
#ymin       : 15.63319 
#ymax       : 39.36653
#suit_minus200 <- raster("monq/models/minus200/individual_projections/monq_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo.tif")
extent(suit_minus200)
#xmin       : -117.0085 
#xmax       : -93.30847 
#ymin       : 15.63319 
#ymax       : 39.36653 
# min(xmin): -117.0085, max(xmax): -93.30847, min(ymin): 15.63319, max(ymax): 39.36653

# below was not run (interpolation was done in ArcGIS)
for (year in c("present", "2040", "2070", "2100")){
  original_raster <- get(paste0("suit_",year))
  filled_raster <- extend(original_raster, rect_extent, value = 0)
  raster_points <- rasterToPoints(original_raster, fun = function(x) !is.na(x))
  raster_values <- as.data.frame(raster_points)
  colnames(raster_values) <- c("x", "y", "value")
  #raster_values <- raster_values[!is.na(raster_values$layer), ]  # Remove NA points
  coordinates(raster_values) <- ~x + y
  proj4string(raster_values) <- crs(original_raster)
  idw_model <- gstat::gstat(formula = value ~ 1, data = raster_values, set = list(idp = 2))
  pred_grd <- raster(original_raster)
  idw_result <- interpolate(pred_grd, idw_model)
  interpolated_raster <- raster(idw_result)
  #colnames(raster_points) <- c("x", "y", "value")
  #coordinates(raster_points) <- ~x + y
  #proj4string(raster_points) <- crs(original_raster)
  #vgm_model <- variogram(value ~ 1, data = raster_points)
  #fit_model <- fit.variogram(vgm_model, model = vgm("Sph"))
  #pred_grid <- as(original_raster, "SpatialGrid")
  #kriging_result <- krige(formula = value ~ 1, locations = raster_points, 
  #                      newdata = pred_grid, model = fit_model)
  #interpolated_raster <- raster(kriging_result)
  filling_raster <- mask(interpolated_raster, original_raster)
  final_raster <- merge(filled_raster, filling_raster)
  writeRaster(final_raster, paste0("MonqSuit_",year,"_final.tif"), format = "GTiff", overwrite = TRUE)
  assign(paste0("suit_",year,"_final"), final_raster)
} 

setwd("D:/Research_Data_Backup/Montezuma_quail/enm/monq_enm_gis")
# Padding and scaling
#for (year in c("present", "2040", "2070", "2100")){
#  original_raster <- get(paste0("suit_",year))
#  no_na_raster <- calc(original_raster, fun = function(x) { ifelse(is.na(x), 0, x) })
#  filled_raster <- extend(no_na_raster, rect_extent, value = 0)
#  scaled_raster <- filled_raster / 1000
#  plot(original_raster)
#  plot(no_na_raster)
#  plot(filled_raster)
#  plot(scaled_raster)
#  writeRaster(scaled_raster, paste0("MonqSuit_",year,"_final.tif"), format = "GTiff", overwrite = TRUE)
#}

#for (year in c("minus0", "minus50", "minus100", "minus150", "minus200")){
#  original_raster <- get(paste0("suit_",year))
#  no_na_raster <- calc(original_raster, fun = function(x) { ifelse(is.na(x), 0, x) })
#  filled_raster <- extend(no_na_raster, rect_extent, value = 0)
#  scaled_raster <- filled_raster / 1000
#  plot(original_raster)
#  plot(no_na_raster)
#  plot(filled_raster)
#  plot(scaled_raster)
#  writeRaster(scaled_raster, paste0("MonqSuit_",year,"_final.tif"), format = "GTiff", overwrite = TRUE)
#}

#for (year in c("present", "2040", "2070", "2100", "minus0", "minus50", "minus100", "minus150", "minus200")){
#  original_raster <- raster(paste0("MonqSuit_",year,"_final.tif"))
#  onezero_raster <- calc(original_raster, fun = function(x) { ifelse(x > 0, 1, 0) })
#  land_raster <- extend(onezero_raster, rect_extent, value = 0)
#  plot(original_raster)
#  plot(onezero_raster)
#  plot(land_raster)
#  writeRaster(land_raster, paste0("MonqLand_",year,"_final.tif"), format = "GTiff", overwrite = TRUE)
#}

for (year in c("present", "2040", "2070", "2100", "minus0", "minus50", "minus100", "minus150", "minus200")){
  original_raster <- raster(paste0("suit_",year,"_res10.tif"))
  print(extent(original_raster))
}

rect_extent <- extent(-117.0085, -93.34181, 15.63319, 39.36653)

for (year in c("minus0", "minus50", "minus100", "minus150", "minus200")){
  original_raster <- raster(paste0("suit_",year,"_res10.tif"))
  no_na_raster <- calc(original_raster, fun = function(x) { ifelse(is.na(x), 0, x) })
  aligned_extent <- alignExtent(rect_extent, no_na_raster, snap = "out") # this should be the same for these layers
  filled_raster <- extend(no_na_raster, aligned_extent, value = 0)
  scaled_raster <- filled_raster / 1000
  print(scaled_raster)
  plot(original_raster)
  plot(no_na_raster)
  plot(filled_raster)
  plot(scaled_raster)
  writeRaster(scaled_raster, paste0("D:/Research_Data_Backup/Montezuma_quail/enm/monq_enm_R/","MonqSuit_",year,"_final_res10.tif"), format = "GTiff", overwrite = TRUE)
}

for (year in c("present", "2040", "2070", "2100")){
  original_raster <- raster(paste0("suit_",year,"_res10.tif"))
  no_na_raster <- calc(original_raster, fun = function(x) { ifelse(is.na(x), 0, x) })
  filled_raster <- extend(no_na_raster, aligned_extent, value = 0)
  scaled_raster <- filled_raster / 1000
  print(scaled_raster)
  plot(original_raster)
  plot(no_na_raster)
  plot(filled_raster)
  plot(scaled_raster)
  writeRaster(scaled_raster, paste0("D:/Research_Data_Backup/Montezuma_quail/enm/monq_enm_R/","MonqSuit_",year,"_final_res10.tif"), format = "GTiff", overwrite = TRUE)
}

for (year in c("present", "2040", "2070", "2100", "minus0", "minus50", "minus100", "minus150", "minus200")){
  original_raster <- raster(paste0("suit_",year,"_res10.tif"))
  onezero_raster <- calc(original_raster, fun = function(x) { ifelse(is.na(x), 0, 1) })
  land_raster <- extend(onezero_raster, rect_extent, value = 0)
  plot(original_raster)
  plot(onezero_raster)
  plot(land_raster)
  writeRaster(land_raster, paste0("D:/Research_Data_Backup/Montezuma_quail/enm/monq_enm_R/","MonqLand_",year,"_final_res10.tif"), format = "GTiff", overwrite = TRUE)
}
