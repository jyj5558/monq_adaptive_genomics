#ml geos
#ml gdal
#ml proj
#module load r

library(ggplot2)
library(tidyr)
library(dplyr)
library(raster)      # or use terra
library(viridis)     # for nice color scales
library(ggspatial)   # optional: for scale bar / compass
library(terra)
library(sf)
library(scales)
library(ggnewscale)
library(scico)
library(rnaturalearth)
library(rnaturalearthdata)

# Geonomix
setwd("/scratch/negishi/jeon96/monq/gnx/GNX_mod-monq_gnx_params_ssp585/")
ncol_gnx <- 258
nrow_gnx <- 291

hsm_0 <- raster("/scratch/negishi/jeon96/monq/gnx/monq_hsm.tif")
hsm_0_df <- as.data.frame(hsm_0, xy = TRUE, na.rm = FALSE)
colnames(hsm_0_df) <- c("x", "y", "suitability")

e <- ext(hsm_0)
x_ras <- c(e$xmin, e$xmax)
y_ras <- c(e$ymin, e$ymax)

#hsm_0_df2 <- hsm_0_df[!is.na(hsm_0_df$suitability) & hsm_0_df$suitability > 0, ]
#ext <- extent(hsm_0_df2)
#x_ras <- c(ext@xmin, ext@xmax)
#y_ras <- c(ext@ymin, ext@ymax)

admin <- ne_countries(scale = "medium", returnclass = "sf")
admin_proj <- st_transform(admin, crs = st_crs(hsm_0))
range <- st_read("/scratch/negishi/jeon96/monq/gea-go/monq_450km_90y_onLand.shp") 
range_proj <- st_transform(range, crs = st_crs(hsm_0))
bbox <- st_bbox(range_proj)

gens <- c(0, 30, 40, 50, 60)

for (gen in gens) {
  loc_file <- sprintf("/scratch/negishi/jeon96/monq/gnx/GNX_mod-monq_gnx_params_ssp585/it--1/spp-MONQ/mod-monq_gnx_params_it--1_t-%d_spp-MONQ.csv", gen)
  location_0 <- read.csv(loc_file, header = TRUE)
  location_0 <- location_0 %>%
    mutate(
      x_norm = (x - 0) / (ncol_gnx - 1),
      y_norm = (y - 0) / (nrow_gnx - 1),
      x_projected = x_norm * diff(x_ras) + x_ras[1],
      y_projected = (1 - y_norm) * diff(y_ras) + y_ras[1]
    )
  
  p <- ggplot() +
    # Basemap
    geom_sf(data = admin_proj, fill=gray(0.95), color = NA) +
    geom_sf(data = admin_proj, fill=NA, color = "black", size = 0.2) +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"]), expand = TRUE) +
  
    # Dummy for legend
    geom_point(data = data.frame(x = -10*e$xmax, y = -10*e$ymax, level = c(0.25, 0.5, 0.75, 0.95)), aes(x = x, y = y, fill = factor(level)), shape = 22, size = 5, alpha = 1, show.legend = TRUE) +
    scale_fill_manual(values = c("0.25" = "#FFFFFF00", "0.5" = "#D8D2A8", "0.75" = "#B8AD70", "0.95" = "#7B6F2A"), name = "Relative density", labels = c("25%", "50%", "75%", "95%")) +

    # Reset the fill scale
    ggnewscale::new_scale_fill() +

    # Background raster layer
    geom_raster(data = dplyr::filter(hsm_0_df, !is.na(suitability) & suitability > 0), aes(x = x, y = y, fill = suitability), alpha = 0.6) +
    scale_fill_scico(name = "Suitability", palette = "vik", direction = -1, na.value = "transparent", limits = c(0, 1)) +

    # Reset the fill scale
    ggnewscale::new_scale_fill() +

    # 2D density layer over simulation points
    stat_density_2d(data = location_0, aes(x = x_projected, y = y_projected, fill = after_stat(level), alpha = after_stat(level) * 0.6), geom = "polygon", color = NA, linewidth = 0.15, show.legend = c(fill = TRUE, alpha = FALSE)) +
    scale_fill_gradientn(colours = c("#FFFFFF00", "#D8D2A8", "#B8AD70", "#7B6F2A"), name = "Relative density", breaks = c(0.25, 0.5, 0.75, 0.95), guide = guide_colorbar(order = 1), labels = c("25%", "50%", "75%", "95%")) +
    scale_alpha_continuous(range = c(0, 0.5), guide = "none") +

    #coord_fixed() +
    theme_minimal(base_size = 14, base_family = "Helvetica") +
    labs(title = "Population Density at Generation 0") +
    theme(legend.position = "right", legend.box = "vertical", plot.title = element_text(hjust = 0.5))
  
  out_svg <- sprintf("hsm_density_t-%d.svg", gen)
  out_pdf <- sprintf("hsm_density_t-%d.pdf", gen)

  ggsave(out_svg, plot = p, width = 10, height = 10, units = "in", device = "svg")
  ggsave(out_pdf, plot = p, width = 10, height = 10, units = "in", device = cairo_pdf)
}

#location_0 <- read.csv("/scratch/negishi/jeon96/monq/gnx/GNX_mod-monq_gnx_params_ssp585/it--1/spp-MONQ/mod-monq_gnx_params_it--1_t-0_spp-MONQ.csv", header = T)
#head(location_0)

#location_0 <- location_0 %>%
#  mutate(
#    x_norm = (x - 0) / (ncol_gnx - 1),
#    y_norm = (y - 0) / (nrow_gnx - 1),
#    x_projected = x_norm * diff(x_ras) + x_ras[1],
#    y_projected = (1 - y_norm) * diff(y_ras) + y_ras[1]  
#  )

#x_loc <- range(location_0$x)
#y_loc <- range(location_0$y)
#location_0$x_norm <- (location_0$x - x_loc[1]) / diff(x_loc)
#location_0$y_norm <- (location_0$y - y_loc[1]) / diff(y_loc)
#location_0$x_projected <- location_0$x_norm * diff(x_ras) + x_ras[1]
#location_0$y_projected <- (1 - location_0$y_norm) * diff(y_ras) + y_ras[1]

#location_0$x <- rescale(location_0$x, to = x_range)
#location_0$y <- rescale(location_0$y, to = rev(y_range))  # flip Y *during* rescaling
#location_0$x <- scales::rescale(location_0$x, to = c(bbox["xmin"], bbox["xmax"])) 
#location_0$y <- scales::rescale(location_0$y, to = c(bbox["ymin"], bbox["ymax"]))
#location_0$y <- max(location_0$y) - location_0$y + min(location_0$y)


# SSP585
hsm_60 <- raster("/scratch/negishi/jeon96/monq/gnx/hsm_ssp585/60_monq_hsm.tif")
hsm_60_df <- as.data.frame(hsm_60, xy = TRUE, na.rm = FALSE)
colnames(hsm_60_df) <- c("x", "y", "suitability")

e <- ext(hsm_60)
x_ras <- c(e$xmin, e$xmax)
y_ras <- c(e$ymin, e$ymax)

#hsm_30_df2 <- hsm_30_df[!is.na(hsm_30_df$suitability) & hsm_30_df$suitability > 0, ]
#ext <- extent(hsm_30_df2)
#x_ras <- c(ext@xmin, ext@xmax)
#y_ras <- c(ext@ymin, ext@ymax)

gens <- c(70, 80, 90)

for (gen in gens) {
  loc_file <- sprintf("/scratch/negishi/jeon96/monq/gnx/GNX_mod-monq_gnx_params_ssp585/it--1/spp-MONQ/mod-monq_gnx_params_it--1_t-%d_spp-MONQ.csv", gen)
  location_60 <- read.csv(loc_file, header = TRUE)
  location_60 <- location_60 %>%
    mutate(
      x_norm = (x - 0) / (ncol_gnx - 1),
      y_norm = (y - 0) / (nrow_gnx - 1),
      x_projected = x_norm * diff(x_ras) + x_ras[1],
      y_projected = (1 - y_norm) * diff(y_ras) + y_ras[1]
    )
  
  p <- ggplot() +
    # Basemap
    geom_sf(data = admin_proj, fill=gray(0.95), color = NA) +
    geom_sf(data = admin_proj, fill=NA, color = "black", size = 0.2) +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"]), expand = TRUE) +
  
    # Dummy for legend
    geom_point(data = data.frame(x = -10*e$xmax, y = -10*e$ymax, level = c(0.25, 0.5, 0.75, 0.95)), aes(x = x, y = y, fill = factor(level)), shape = 22, size = 5, alpha = 1, show.legend = TRUE) +
    scale_fill_manual(values = c("0.25" = "#FFFFFF00", "0.5" = "#D8D2A8", "0.75" = "#B8AD70", "0.95" = "#7B6F2A"), name = "Relative density", labels = c("25%", "50%", "75%", "95%")) +

    # Reset the fill scale
    ggnewscale::new_scale_fill() +

    # Background raster layer
    geom_raster(data = dplyr::filter(hsm_60_df, !is.na(suitability) & suitability > 0), aes(x = x, y = y, fill = suitability), alpha = 0.6) +
    scale_fill_scico(name = "Suitability", palette = "vik", direction = -1, na.value = "transparent", limits = c(0, 1)) +

    # Reset the fill scale
    ggnewscale::new_scale_fill() +

    # 2D density layer over simulation points
    stat_density_2d(data = location_60, aes(x = x_projected, y = y_projected, fill = after_stat(level), alpha = after_stat(level) * 0.6), geom = "polygon", color = NA, linewidth = 0.15, show.legend = c(fill = TRUE, alpha = FALSE)) +
    scale_fill_gradientn(colours = c("#FFFFFF00", "#D8D2A8", "#B8AD70", "#7B6F2A"), name = "Relative density", breaks = c(0.25, 0.5, 0.75, 0.95), guide = guide_colorbar(order = 1), labels = c("25%", "50%", "75%", "95%")) +
    scale_alpha_continuous(range = c(0, 0.5), guide = "none") +

    #coord_fixed() +
    theme_minimal(base_size = 14, base_family = "Helvetica") +
    labs(title = "Population Density at Generation 0") +
    theme(legend.position = "right", legend.box = "vertical", plot.title = element_text(hjust = 0.5))
  
  out_svg <- sprintf("hsm_density_t-%d.svg", gen)
  out_pdf <- sprintf("hsm_density_t-%d.pdf", gen)

  ggsave(out_svg, plot = p, width = 10, height = 10, units = "in", device = "svg")
  ggsave(out_pdf, plot = p, width = 10, height = 10, units = "in", device = cairo_pdf)
}


hsm_90 <- raster("/scratch/negishi/jeon96/monq/gnx/hsm_ssp585/90_monq_hsm.tif")
hsm_90_df <- as.data.frame(hsm_90, xy = TRUE, na.rm = FALSE)
colnames(hsm_90_df) <- c("x", "y", "suitability")

e <- ext(hsm_90)
x_ras <- c(e$xmin, e$xmax)
y_ras <- c(e$ymin, e$ymax)

gens <- c(100, 110, 120)

for (gen in gens) {
  loc_file <- sprintf("/scratch/negishi/jeon96/monq/gnx/GNX_mod-monq_gnx_params_ssp585/it--1/spp-MONQ/mod-monq_gnx_params_it--1_t-%d_spp-MONQ.csv", gen)
  location_90 <- read.csv(loc_file, header = TRUE)
  location_90 <- location_90 %>%
    mutate(
      x_norm = (x - 0) / (ncol_gnx - 1),
      y_norm = (y - 0) / (nrow_gnx - 1),
      x_projected = x_norm * diff(x_ras) + x_ras[1],
      y_projected = (1 - y_norm) * diff(y_ras) + y_ras[1]
    )
  
  p <- ggplot() +
    # Basemap
    geom_sf(data = admin_proj, fill=gray(0.95), color = NA) +
    geom_sf(data = admin_proj, fill=NA, color = "black", size = 0.2) +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"]), expand = TRUE) +
  
    # Dummy for legend
    geom_point(data = data.frame(x = -10*e$xmax, y = -10*e$ymax, level = c(0.25, 0.5, 0.75, 0.95)), aes(x = x, y = y, fill = factor(level)), shape = 22, size = 5, alpha = 1, show.legend = TRUE) +
    scale_fill_manual(values = c("0.25" = "#FFFFFF00", "0.5" = "#D8D2A8", "0.75" = "#B8AD70", "0.95" = "#7B6F2A"), name = "Relative density", labels = c("25%", "50%", "75%", "95%")) +

    # Reset the fill scale
    ggnewscale::new_scale_fill() +

    # Background raster layer
    geom_raster(data = dplyr::filter(hsm_90_df, !is.na(suitability) & suitability > 0), aes(x = x, y = y, fill = suitability), alpha = 0.6) +
    scale_fill_scico(name = "Suitability", palette = "vik", direction = -1, na.value = "transparent", limits = c(0, 1)) +

    # Reset the fill scale
    ggnewscale::new_scale_fill() +

    # 2D density layer over simulation points
    stat_density_2d(data = location_90, aes(x = x_projected, y = y_projected, fill = after_stat(level), alpha = after_stat(level) * 0.6), geom = "polygon", color = NA, linewidth = 0.15, show.legend = c(fill = TRUE, alpha = FALSE)) +
    scale_fill_gradientn(colours = c("#FFFFFF00", "#D8D2A8", "#B8AD70", "#7B6F2A"), name = "Relative density", breaks = c(0.25, 0.5, 0.75, 0.95), guide = guide_colorbar(order = 1), labels = c("25%", "50%", "75%", "95%")) +
    scale_alpha_continuous(range = c(0, 0.5), guide = "none") +

    #coord_fixed() +
    theme_minimal(base_size = 14, base_family = "Helvetica") +
    labs(title = "Population Density at Generation 0") +
    theme(legend.position = "right", legend.box = "vertical", plot.title = element_text(hjust = 0.5))
  
  out_svg <- sprintf("hsm_density_t-%d.svg", gen)
  out_pdf <- sprintf("hsm_density_t-%d.pdf", gen)

  ggsave(out_svg, plot = p, width = 10, height = 10, units = "in", device = "svg")
  ggsave(out_pdf, plot = p, width = 10, height = 10, units = "in", device = cairo_pdf)
}


hsm_120 <- raster("/scratch/negishi/jeon96/monq/gnx/hsm_ssp585/120_monq_hsm.tif")
hsm_120_df <- as.data.frame(hsm_120, xy = TRUE, na.rm = FALSE)
colnames(hsm_120_df) <- c("x", "y", "suitability")

e <- ext(hsm_120)
x_ras <- c(e$xmin, e$xmax)
y_ras <- c(e$ymin, e$ymax)

gens <- c(130, 140, 149)

for (gen in gens) {
  loc_file <- sprintf("/scratch/negishi/jeon96/monq/gnx/GNX_mod-monq_gnx_params_ssp585/it--1/spp-MONQ/mod-monq_gnx_params_it--1_t-%d_spp-MONQ.csv", gen)
  location_120 <- read.csv(loc_file, header = TRUE)
  location_120 <- location_120 %>%
    mutate(
      x_norm = (x - 0) / (ncol_gnx - 1),
      y_norm = (y - 0) / (nrow_gnx - 1),
      x_projected = x_norm * diff(x_ras) + x_ras[1],
      y_projected = (1 - y_norm) * diff(y_ras) + y_ras[1]
    )
  
  p <- ggplot() +
    # Basemap
    geom_sf(data = admin_proj, fill=gray(0.95), color = NA) +
    geom_sf(data = admin_proj, fill=NA, color = "black", size = 0.2) +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"]), expand = TRUE) +
  
    # Dummy for legend
    geom_point(data = data.frame(x = -10*e$xmax, y = -10*e$ymax, level = c(0.25, 0.5, 0.75, 0.95)), aes(x = x, y = y, fill = factor(level)), shape = 22, size = 5, alpha = 1, show.legend = TRUE) +
    scale_fill_manual(values = c("0.25" = "#FFFFFF00", "0.5" = "#D8D2A8", "0.75" = "#B8AD70", "0.95" = "#7B6F2A"), name = "Relative density", labels = c("25%", "50%", "75%", "95%")) +

    # Reset the fill scale
    ggnewscale::new_scale_fill() +

    # Background raster layer
    geom_raster(data = dplyr::filter(hsm_120_df, !is.na(suitability) & suitability > 0), aes(x = x, y = y, fill = suitability), alpha = 0.6) +
    scale_fill_scico(name = "Suitability", palette = "vik", direction = -1, na.value = "transparent", limits = c(0, 1)) +

    # Reset the fill scale
    ggnewscale::new_scale_fill() +

    # 2D density layer over simulation points
    stat_density_2d(data = location_120, aes(x = x_projected, y = y_projected, fill = after_stat(level), alpha = after_stat(level) * 0.6), geom = "polygon", color = NA, linewidth = 0.15, show.legend = c(fill = TRUE, alpha = FALSE)) +
    scale_fill_gradientn(colours = c("#FFFFFF00", "#D8D2A8", "#B8AD70", "#7B6F2A"), name = "Relative density", breaks = c(0.25, 0.5, 0.75, 0.95), guide = guide_colorbar(order = 1), labels = c("25%", "50%", "75%", "95%")) +
    scale_alpha_continuous(range = c(0, 0.5), guide = "none") +

    #coord_fixed() +
    theme_minimal(base_size = 14, base_family = "Helvetica") +
    labs(title = "Population Density at Generation 0") +
    theme(legend.position = "right", legend.box = "vertical", plot.title = element_text(hjust = 0.5))
  
  out_svg <- sprintf("hsm_density_t-%d.svg", gen)
  out_pdf <- sprintf("hsm_density_t-%d.pdf", gen)

  ggsave(out_svg, plot = p, width = 10, height = 10, units = "in", device = "svg")
  ggsave(out_pdf, plot = p, width = 10, height = 10, units = "in", device = cairo_pdf)
}


df <- read.csv("/scratch/negishi/jeon96/monq/gnx/GNX_mod-monq_gnx_params_ssp585/it--1/spp-MONQ/mod-monq_gnx_params_it--1_spp-MONQ_OTHER_STATS_modif.csv")  # Replace with your actual file path; removed 0-30 generations which are part of non-neutral burn-in

df_long <- df %>%
  pivot_longer(cols = c(Nt, mean_fit),  # Add more if needed
               names_to = "type", values_to = "value") #%>%
#  filter(!(t > 15 & t <= 45))

df_long$type <- factor(df_long$type,
                      levels = c("Nt", "mean_fit"),
                      labels = c("Population Size", "Mean Fitness"))

# Generate the line plot
#p <- ggplot(df_long, aes(x = t, y = value, color = type)) +
p <- ggplot(df, aes(x = t, y = Nt)) +
  geom_line(size = 0.5) +
  #facet_wrap(~ type, scales = "free_y") +  # separate y-axis per panel
  theme_minimal(base_size = 14, base_family = "Helvetica") +
  theme(
    axis.text = element_text(size = 14),
    #legend.title = element_blank(),
    legend.text = element_text(size = 12),
    #plot.title = element_text(hjust = 0.5)
  ) +
#  labs(x = "Future Generation")
  labs(x = "Future Generation", y = "Population Size") +
  ylim(500000, 1400000)
ggsave("Pop_stats_gnx.svg", plot = p, width = 15, height = 10, units = "in", device = "svg")
ggsave("Pop_stats_gnx.pdf", plot = p, width = 15, height = 10, units = "in", device = cairo_pdf)

p <- ggplot(df, aes(x = t, y = Nt)) +
  geom_line(size = 0.5) +
  #facet_wrap(~ type, scales = "free_y") +  # separate y-axis per panel
  theme_minimal(base_size = 14, base_family = "Helvetica") +
  theme(
    axis.text = element_text(size = 14),
    #legend.title = element_blank(),
    legend.text = element_text(size = 12),
    #plot.title = element_text(hjust = 0.5)
  ) +
#  labs(x = "Future Generation")
  labs(x = "Future Generation", y = "Population Size") +
  xlim(0, 30) + ylim(900000, 1400000)
ggsave("Pop_stats_gnx_gen30.svg", plot = p, width = 6, height = 9, units = "in", device = "svg")
ggsave("Pop_stats_gnx_gen30.pdf", plot = p, width = 6, height = 9, units = "in", device = cairo_pdf)

p <- ggplot(df, aes(x = t, y = Nt)) +
  geom_line(size = 0.5) +
  #facet_wrap(~ type, scales = "free_y") +  # separate y-axis per panel
  theme_minimal(base_size = 14, base_family = "Helvetica") +
  theme(
    axis.text = element_text(size = 14),
    #legend.title = element_blank(),
    legend.text = element_text(size = 12),
    #plot.title = element_text(hjust = 0.5)
  ) +
#  labs(x = "Future Generation")
  labs(x = "Future Generation", y = "Population Size") +
  xlim(31, 60) + ylim(900000, 1400000)
ggsave("Pop_stats_gnx_gen60.svg", plot = p, width = 6, height = 9, units = "in", device = "svg")
ggsave("Pop_stats_gnx_gen60.pdf", plot = p, width = 6, height = 9, units = "in", device = cairo_pdf)

p <- ggplot(df, aes(x = t, y = Nt)) +
  geom_line(size = 0.5) +
  #facet_wrap(~ type, scales = "free_y") +  # separate y-axis per panel
  theme_minimal(base_size = 14, base_family = "Helvetica") +
  theme(
    axis.text = element_text(size = 14),
    #legend.title = element_blank(),
    legend.text = element_text(size = 12),
    #plot.title = element_text(hjust = 0.5)
  ) +
#  labs(x = "Future Generation")
  labs(x = "Future Generation", y = "Population Size") +
  xlim(61, 90) + ylim(900000, 1400000)
ggsave("Pop_stats_gnx_gen90.svg", plot = p, width = 6, height = 9, units = "in", device = "svg")
ggsave("Pop_stats_gnx_gen90.pdf", plot = p, width = 6, height = 9, units = "in", device = cairo_pdf)

p <- ggplot(df, aes(x = t, y = Nt)) +
  geom_line(size = 0.5) +
  #facet_wrap(~ type, scales = "free_y") +  # separate y-axis per panel
  theme_minimal(base_size = 14, base_family = "Helvetica") +
  theme(
    axis.text = element_text(size = 14),
    #legend.title = element_blank(),
    legend.text = element_text(size = 12),
    #plot.title = element_text(hjust = 0.5)
  ) +
#  labs(x = "Future Generation")
  labs(x = "Future Generation", y = "Population Size") +
  xlim(91, 120) + ylim(900000, 1400000)
ggsave("Pop_stats_gnx_gen120.svg", plot = p, width = 6, height = 9, units = "in", device = "svg")
ggsave("Pop_stats_gnx_gen120.pdf", plot = p, width = 6, height = 9, units = "in", device = cairo_pdf)


# SSP126
setwd("/scratch/negishi/jeon96/monq/gnx/GNX_mod-monq_gnx_params_ssp126/")
ncol_gnx <- 258
nrow_gnx <- 291

hsm_0 <- raster("/scratch/negishi/jeon96/monq/gnx/monq_hsm.tif")
hsm_0_df <- as.data.frame(hsm_0, xy = TRUE, na.rm = FALSE)
colnames(hsm_0_df) <- c("x", "y", "suitability")

e <- ext(hsm_0)
x_ras <- c(e$xmin, e$xmax)
y_ras <- c(e$ymin, e$ymax)

#hsm_0_df2 <- hsm_0_df[!is.na(hsm_0_df$suitability) & hsm_0_df$suitability > 0, ]
#ext <- extent(hsm_0_df2)
#x_ras <- c(ext@xmin, ext@xmax)
#y_ras <- c(ext@ymin, ext@ymax)

admin <- ne_countries(scale = "medium", returnclass = "sf")
admin_proj <- st_transform(admin, crs = st_crs(hsm_0))
range <- st_read("/scratch/negishi/jeon96/monq/gea-go/monq_450km_90y_onLand.shp") 
range_proj <- st_transform(range, crs = st_crs(hsm_0))
bbox <- st_bbox(range_proj)

gens <- c(0, 30, 40, 50, 60)

for (gen in gens) {
  loc_file <- sprintf("/scratch/negishi/jeon96/monq/gnx/GNX_mod-monq_gnx_params_ssp126/it--1/spp-MONQ/mod-monq_gnx_params_it--1_t-%d_spp-MONQ.csv", gen)
  location_0 <- read.csv(loc_file, header = TRUE)
  location_0 <- location_0 %>%
    mutate(
      x_norm = (x - 0) / (ncol_gnx - 1),
      y_norm = (y - 0) / (nrow_gnx - 1),
      x_projected = x_norm * diff(x_ras) + x_ras[1],
      y_projected = (1 - y_norm) * diff(y_ras) + y_ras[1]
    )
  
  p <- ggplot() +
    # Basemap
    geom_sf(data = admin_proj, fill=gray(0.95), color = NA) +
    geom_sf(data = admin_proj, fill=NA, color = "black", size = 0.2) +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"]), expand = TRUE) +
  
    # Dummy for legend
    geom_point(data = data.frame(x = -10*e$xmax, y = -10*e$ymax, level = c(0.25, 0.5, 0.75, 0.95)), aes(x = x, y = y, fill = factor(level)), shape = 22, size = 5, alpha = 1, show.legend = TRUE) +
    scale_fill_manual(values = c("0.25" = "#FFFFFF00", "0.5" = "#D8D2A8", "0.75" = "#B8AD70", "0.95" = "#7B6F2A"), name = "Relative density", labels = c("25%", "50%", "75%", "95%")) +

    # Reset the fill scale
    ggnewscale::new_scale_fill() +

    # Background raster layer
    geom_raster(data = dplyr::filter(hsm_0_df, !is.na(suitability) & suitability > 0), aes(x = x, y = y, fill = suitability), alpha = 0.6) +
    scale_fill_scico(name = "Suitability", palette = "vik", direction = -1, na.value = "transparent", limits = c(0, 1)) +

    # Reset the fill scale
    ggnewscale::new_scale_fill() +

    # 2D density layer over simulation points
    stat_density_2d(data = location_0, aes(x = x_projected, y = y_projected, fill = after_stat(level), alpha = after_stat(level) * 0.6), geom = "polygon", color = NA, linewidth = 0.15, show.legend = c(fill = TRUE, alpha = FALSE)) +
    scale_fill_gradientn(colours = c("#FFFFFF00", "#D8D2A8", "#B8AD70", "#7B6F2A"), name = "Relative density", breaks = c(0.25, 0.5, 0.75, 0.95), guide = guide_colorbar(order = 1), labels = c("25%", "50%", "75%", "95%")) +
    scale_alpha_continuous(range = c(0, 0.5), guide = "none") +

    #coord_fixed() +
    theme_minimal(base_size = 14, base_family = "Helvetica") +
    labs(title = "Population Density at Generation 0") +
    theme(legend.position = "right", legend.box = "vertical", plot.title = element_text(hjust = 0.5))
  
  out_svg <- sprintf("hsm_density_t-%d.svg", gen)
  out_pdf <- sprintf("hsm_density_t-%d.pdf", gen)

  ggsave(out_svg, plot = p, width = 10, height = 10, units = "in", device = "svg")
  ggsave(out_pdf, plot = p, width = 10, height = 10, units = "in", device = cairo_pdf)
}


hsm_60 <- raster("/scratch/negishi/jeon96/monq/gnx/hsm_ssp126/60_monq_hsm.tif")
hsm_60_df <- as.data.frame(hsm_60, xy = TRUE, na.rm = FALSE)
colnames(hsm_60_df) <- c("x", "y", "suitability")

e <- ext(hsm_60)
x_ras <- c(e$xmin, e$xmax)
y_ras <- c(e$ymin, e$ymax)

#hsm_30_df2 <- hsm_30_df[!is.na(hsm_30_df$suitability) & hsm_30_df$suitability > 0, ]
#ext <- extent(hsm_30_df2)
#x_ras <- c(ext@xmin, ext@xmax)
#y_ras <- c(ext@ymin, ext@ymax)

gens <- c(70, 80, 90)

for (gen in gens) {
  loc_file <- sprintf("/scratch/negishi/jeon96/monq/gnx/GNX_mod-monq_gnx_params_ssp126/it--1/spp-MONQ/mod-monq_gnx_params_it--1_t-%d_spp-MONQ.csv", gen)
  location_60 <- read.csv(loc_file, header = TRUE)
  location_60 <- location_60 %>%
    mutate(
      x_norm = (x - 0) / (ncol_gnx - 1),
      y_norm = (y - 0) / (nrow_gnx - 1),
      x_projected = x_norm * diff(x_ras) + x_ras[1],
      y_projected = (1 - y_norm) * diff(y_ras) + y_ras[1]
    )
  
  p <- ggplot() +
    # Basemap
    geom_sf(data = admin_proj, fill=gray(0.95), color = NA) +
    geom_sf(data = admin_proj, fill=NA, color = "black", size = 0.2) +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"]), expand = TRUE) +
  
    # Dummy for legend
    geom_point(data = data.frame(x = -10*e$xmax, y = -10*e$ymax, level = c(0.25, 0.5, 0.75, 0.95)), aes(x = x, y = y, fill = factor(level)), shape = 22, size = 5, alpha = 1, show.legend = TRUE) +
    scale_fill_manual(values = c("0.25" = "#FFFFFF00", "0.5" = "#D8D2A8", "0.75" = "#B8AD70", "0.95" = "#7B6F2A"), name = "Relative density", labels = c("25%", "50%", "75%", "95%")) +

    # Reset the fill scale
    ggnewscale::new_scale_fill() +

    # Background raster layer
    geom_raster(data = dplyr::filter(hsm_60_df, !is.na(suitability) & suitability > 0), aes(x = x, y = y, fill = suitability), alpha = 0.6) +
    scale_fill_scico(name = "Suitability", palette = "vik", direction = -1, na.value = "transparent", limits = c(0, 1)) +

    # Reset the fill scale
    ggnewscale::new_scale_fill() +

    # 2D density layer over simulation points
    stat_density_2d(data = location_60, aes(x = x_projected, y = y_projected, fill = after_stat(level), alpha = after_stat(level) * 0.6), geom = "polygon", color = NA, linewidth = 0.15, show.legend = c(fill = TRUE, alpha = FALSE)) +
    scale_fill_gradientn(colours = c("#FFFFFF00", "#D8D2A8", "#B8AD70", "#7B6F2A"), name = "Relative density", breaks = c(0.25, 0.5, 0.75, 0.95), guide = guide_colorbar(order = 1), labels = c("25%", "50%", "75%", "95%")) +
    scale_alpha_continuous(range = c(0, 0.5), guide = "none") +

    #coord_fixed() +
    theme_minimal(base_size = 14, base_family = "Helvetica") +
    labs(title = "Population Density at Generation 0") +
    theme(legend.position = "right", legend.box = "vertical", plot.title = element_text(hjust = 0.5))
  
  out_svg <- sprintf("hsm_density_t-%d.svg", gen)
  out_pdf <- sprintf("hsm_density_t-%d.pdf", gen)

  ggsave(out_svg, plot = p, width = 10, height = 10, units = "in", device = "svg")
  ggsave(out_pdf, plot = p, width = 10, height = 10, units = "in", device = cairo_pdf)
}


hsm_90 <- raster("/scratch/negishi/jeon96/monq/gnx/hsm_ssp126/90_monq_hsm.tif")
hsm_90_df <- as.data.frame(hsm_90, xy = TRUE, na.rm = FALSE)
colnames(hsm_90_df) <- c("x", "y", "suitability")

e <- ext(hsm_90)
x_ras <- c(e$xmin, e$xmax)
y_ras <- c(e$ymin, e$ymax)

gens <- c(100, 110, 120)

for (gen in gens) {
  loc_file <- sprintf("/scratch/negishi/jeon96/monq/gnx/GNX_mod-monq_gnx_params_ssp126/it--1/spp-MONQ/mod-monq_gnx_params_it--1_t-%d_spp-MONQ.csv", gen)
  location_90 <- read.csv(loc_file, header = TRUE)
  location_90 <- location_90 %>%
    mutate(
      x_norm = (x - 0) / (ncol_gnx - 1),
      y_norm = (y - 0) / (nrow_gnx - 1),
      x_projected = x_norm * diff(x_ras) + x_ras[1],
      y_projected = (1 - y_norm) * diff(y_ras) + y_ras[1]
    )
  
  p <- ggplot() +
    # Basemap
    geom_sf(data = admin_proj, fill=gray(0.95), color = NA) +
    geom_sf(data = admin_proj, fill=NA, color = "black", size = 0.2) +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"]), expand = TRUE) +
  
    # Dummy for legend
    geom_point(data = data.frame(x = -10*e$xmax, y = -10*e$ymax, level = c(0.25, 0.5, 0.75, 0.95)), aes(x = x, y = y, fill = factor(level)), shape = 22, size = 5, alpha = 1, show.legend = TRUE) +
    scale_fill_manual(values = c("0.25" = "#FFFFFF00", "0.5" = "#D8D2A8", "0.75" = "#B8AD70", "0.95" = "#7B6F2A"), name = "Relative density", labels = c("25%", "50%", "75%", "95%")) +

    # Reset the fill scale
    ggnewscale::new_scale_fill() +

    # Background raster layer
    geom_raster(data = dplyr::filter(hsm_90_df, !is.na(suitability) & suitability > 0), aes(x = x, y = y, fill = suitability), alpha = 0.6) +
    scale_fill_scico(name = "Suitability", palette = "vik", direction = -1, na.value = "transparent", limits = c(0, 1)) +

    # Reset the fill scale
    ggnewscale::new_scale_fill() +

    # 2D density layer over simulation points
    stat_density_2d(data = location_90, aes(x = x_projected, y = y_projected, fill = after_stat(level), alpha = after_stat(level) * 0.6), geom = "polygon", color = NA, linewidth = 0.15, show.legend = c(fill = TRUE, alpha = FALSE)) +
    scale_fill_gradientn(colours = c("#FFFFFF00", "#D8D2A8", "#B8AD70", "#7B6F2A"), name = "Relative density", breaks = c(0.25, 0.5, 0.75, 0.95), guide = guide_colorbar(order = 1), labels = c("25%", "50%", "75%", "95%")) +
    scale_alpha_continuous(range = c(0, 0.5), guide = "none") +

    #coord_fixed() +
    theme_minimal(base_size = 14, base_family = "Helvetica") +
    labs(title = "Population Density at Generation 0") +
    theme(legend.position = "right", legend.box = "vertical", plot.title = element_text(hjust = 0.5))
  
  out_svg <- sprintf("hsm_density_t-%d.svg", gen)
  out_pdf <- sprintf("hsm_density_t-%d.pdf", gen)

  ggsave(out_svg, plot = p, width = 10, height = 10, units = "in", device = "svg")
  ggsave(out_pdf, plot = p, width = 10, height = 10, units = "in", device = cairo_pdf)
}


hsm_120 <- raster("/scratch/negishi/jeon96/monq/gnx/hsm_ssp126/120_monq_hsm.tif")
hsm_120_df <- as.data.frame(hsm_120, xy = TRUE, na.rm = FALSE)
colnames(hsm_120_df) <- c("x", "y", "suitability")

e <- ext(hsm_120)
x_ras <- c(e$xmin, e$xmax)
y_ras <- c(e$ymin, e$ymax)

gens <- c(130, 140, 149)

for (gen in gens) {
  loc_file <- sprintf("/scratch/negishi/jeon96/monq/gnx/GNX_mod-monq_gnx_params_ssp126/it--1/spp-MONQ/mod-monq_gnx_params_it--1_t-%d_spp-MONQ.csv", gen)
  location_120 <- read.csv(loc_file, header = TRUE)
  location_120 <- location_120 %>%
    mutate(
      x_norm = (x - 0) / (ncol_gnx - 1),
      y_norm = (y - 0) / (nrow_gnx - 1),
      x_projected = x_norm * diff(x_ras) + x_ras[1],
      y_projected = (1 - y_norm) * diff(y_ras) + y_ras[1]
    )
  
  p <- ggplot() +
    # Basemap
    geom_sf(data = admin_proj, fill=gray(0.95), color = NA) +
    geom_sf(data = admin_proj, fill=NA, color = "black", size = 0.2) +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"]), expand = TRUE) +
  
    # Dummy for legend
    geom_point(data = data.frame(x = -10*e$xmax, y = -10*e$ymax, level = c(0.25, 0.5, 0.75, 0.95)), aes(x = x, y = y, fill = factor(level)), shape = 22, size = 5, alpha = 1, show.legend = TRUE) +
    scale_fill_manual(values = c("0.25" = "#FFFFFF00", "0.5" = "#D8D2A8", "0.75" = "#B8AD70", "0.95" = "#7B6F2A"), name = "Relative density", labels = c("25%", "50%", "75%", "95%")) +

    # Reset the fill scale
    ggnewscale::new_scale_fill() +

    # Background raster layer
    geom_raster(data = dplyr::filter(hsm_120_df, !is.na(suitability) & suitability > 0), aes(x = x, y = y, fill = suitability), alpha = 0.6) +
    scale_fill_scico(name = "Suitability", palette = "vik", direction = -1, na.value = "transparent", limits = c(0, 1)) +

    # Reset the fill scale
    ggnewscale::new_scale_fill() +

    # 2D density layer over simulation points
    stat_density_2d(data = location_120, aes(x = x_projected, y = y_projected, fill = after_stat(level), alpha = after_stat(level) * 0.6), geom = "polygon", color = NA, linewidth = 0.15, show.legend = c(fill = TRUE, alpha = FALSE)) +
    scale_fill_gradientn(colours = c("#FFFFFF00", "#D8D2A8", "#B8AD70", "#7B6F2A"), name = "Relative density", breaks = c(0.25, 0.5, 0.75, 0.95), guide = guide_colorbar(order = 1), labels = c("25%", "50%", "75%", "95%")) +
    scale_alpha_continuous(range = c(0, 0.5), guide = "none") +

    #coord_fixed() +
    theme_minimal(base_size = 14, base_family = "Helvetica") +
    labs(title = "Population Density at Generation 0") +
    theme(legend.position = "right", legend.box = "vertical", plot.title = element_text(hjust = 0.5))
  
  out_svg <- sprintf("hsm_density_t-%d.svg", gen)
  out_pdf <- sprintf("hsm_density_t-%d.pdf", gen)

  ggsave(out_svg, plot = p, width = 10, height = 10, units = "in", device = "svg")
  ggsave(out_pdf, plot = p, width = 10, height = 10, units = "in", device = cairo_pdf)
}



# SLIM
# Load the CSV file
setwd("/scratch/negishi/jeon96/monq/slim/output_K1000_SSP585/")
df <- read.csv("/scratch/negishi/jeon96/monq/slim/output_K1000_SSP585/stats_ssp585_merged.txt")  # Replace with your actual file path

# Convert to long format (specify your columns explicitly or with a common prefix)
df_long1 <- df %>%
  pivot_longer(cols = c(Num_vstrDel, Num_strDel, Num_modDel, Num_wkDel),  # Add more if needed
               names_to = "type", values_to = "value")

df_long1$type <- factor(df_long1$type,
                       levels = c("Num_vstrDel", "Num_strDel", "Num_modDel", "Num_wkDel"),
                       labels = c("Very strong", "Strong", "Moderate", "Weak"))

# Generate the line plot
p <- ggplot(df_long1, aes(x = Tick, y = value, color = type)) +
  geom_line(size = 0.5) +
  scale_color_manual(values = c("Very strong" = "#d73027", "Strong" = "#fc8d59", "Moderate" = "#91bfdb", "Weak" = "#4575b4"), name = NULL) +
  theme_minimal(base_size = 14, base_family = "Helvetica") +
  theme(
    axis.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(x = "Future Generation", y = "Number of Deleterious Mutations")
ggsave("Num_del_slim.svg", plot = p, width = 10, height = 6, units = "in", device = "svg")
ggsave("Num_del_slim.pdf", plot = p, width = 10, height = 6, units = "in", device = cairo_pdf)


# Convert to long format (specify your columns explicitly or with a common prefix)
df_long2 <- df %>%
  pivot_longer(cols = c(Num_vstrBen, Num_strBen, Num_modBen, Num_wkBen),  # Add more if needed
               names_to = "type", values_to = "value")

df_long2$type <- factor(df_long2$type,
                       levels = c("Num_vstrBen", "Num_strBen", "Num_modBen", "Num_wkBen"),
                       labels = c("Very strong", "Strong", "Moderate", "Weak"))

# Generate the line plot
p <- ggplot(df_long2, aes(x = Tick, y = value, color = type)) +
  geom_line(size = 0.5) +
  scale_color_manual(values = c("Very strong" = "#d73027", "Strong" = "#fc8d59", "Moderate" = "#91bfdb", "Weak" = "#4575b4"), name = NULL) +
  theme_minimal(base_size = 14, base_family = "Helvetica") +
  theme(
    axis.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(x = "Future Generation", y = "Number of Beneficial Mutations")
ggsave("Num_ben_slim.svg", plot = p, width = 10, height = 6, units = "in", device = "svg")
ggsave("Num_ben_slim.pdf", plot = p, width = 10, height = 6, units = "in", device = cairo_pdf)


# Convert to long format (specify your columns explicitly or with a common prefix)
df_long3 <- df %>%
  pivot_longer(cols = c(PopSize, MeanFitness, InbreedingLoad),  # Add more if needed
               names_to = "type", values_to = "value") %>%
  filter(Tick > 10)

df_long3$type <- factor(df_long3$type,
                       levels = c("PopSize", "MeanFitness", "InbreedingLoad"),
                       labels = c("Population Size", "Mean Fitness", "Inbreeding Load"))

# Generate the line plot
p <- ggplot(df_long3, aes(x = Tick, y = value, color = type)) +
  geom_line(size = 0.5) +
  facet_wrap(~ type, scales = "free_y") +  # separate y-axis per panel
  theme_minimal(base_size = 14, base_family = "Helvetica") +
  theme(
    axis.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(x = "Future Generation")
ggsave("Pop_stats_slim.svg", plot = p, width = 15, height = 6, units = "in", device = "svg")
ggsave("Pop_stats_slim.pdf", plot = p, width = 15, height = 6, units = "in", device = cairo_pdf)


# Convert to long format (specify your columns explicitly or with a common prefix)
df_long4 <- df %>%
  pivot_longer(cols = c(Heterozygosity,Froh100Kb),  # Add more if needed
               names_to = "type", values_to = "value") %>%
  filter(Tick > 10)

df_long4$type <- factor(df_long4$type,
                       levels = c("Heterozygosity", "Froh100Kb"),
                       labels = c("Heterozygosity", "Froh100Kb"))

# Generate the line plot
p <- ggplot(df_long4, aes(x = Tick, y = value, color = type)) +
  geom_line(size = 0.5) +
  facet_wrap(~ type, scales = "free_y") +  # separate y-axis per panel
  theme_minimal(base_size = 14, base_family = "Helvetica") +
  theme(
    axis.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(x = "Future Generation")


# Convert to long format (specify your columns explicitly or with a common prefix)
df_long5 <- df %>%
  pivot_longer(cols = c(LoadP_high,LoadP_mod,LoadP_low,LoadP_tot),  # Add more if needed
               names_to = "type", values_to = "value") %>%
  filter(Tick > 10)

df_long5$type <- factor(df_long5$type,
                       levels = c("LoadP_high", "LoadP_mod", "LoadP_low", "LoadP_tot"),
                       labels = c("LoadP_high", "LoadP_mod", "LoadP_low", "LoadP_tot"))

# Generate the line plot
p <- ggplot(df_long5, aes(x = Tick, y = value, color = type)) +
  geom_line(size = 0.5) +
  facet_wrap(~ type, scales = "free_y") +  # separate y-axis per panel
  theme_minimal(base_size = 14, base_family = "Helvetica") +
  theme(
    axis.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(x = "Future Generation")


# Convert to long format (specify your columns explicitly or with a common prefix)
df_long6 <- df %>%
  pivot_longer(cols = c(LoadR_high,LoadR_mod,LoadR_low,LoadR_tot),  # Add more if needed
               names_to = "type", values_to = "value") %>%
  filter(Tick > 10)

df_long6$type <- factor(df_long6$type,
                       levels = c("LoadR_high", "LoadR_mod", "LoadR_low", "LoadR_tot"),
                       labels = c("LoadR_high", "LoadR_mod", "LoadR_low", "LoadR_tot"))

# Generate the line plot
p <- ggplot(df_long6, aes(x = Tick, y = value, color = type)) +
  geom_line(size = 0.5) +
  facet_wrap(~ type, scales = "free_y") +  # separate y-axis per panel
  theme_minimal(base_size = 14, base_family = "Helvetica") +
  theme(
    axis.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(x = "Future Generation")



library(tidyverse)

# paths (adjust if needed)
f126 <- "stats_ssp126_merged.txt"
f585 <- "stats_ssp585_merged.txt"

read_stats <- function(path, scenario){
  readr::read_csv(path, show_col_types = FALSE) %>%
    mutate(Scenario = scenario)
}

df <- bind_rows(
  read_stats(f126, "SSP126"),
  read_stats(f585, "SSP585")
)

# Define years since burn-in separately for each scenario/pop if needed:
# (Usually Tick starts at the same burn-in tick in both, but this is safer.)
df <- df %>%
  group_by(Scenario, Population) %>%
  mutate(Year = Tick - min(Tick)) %>%
  ungroup()

#df_post <- df %>% filter(Year > 10) # If you want to exclude the very first few ticks:


df_f1 <- df %>%
  dplyr::select(Year, Scenario, Population, PopSize, MeanFitness) %>%
  tidyr::pivot_longer(cols = c(PopSize, MeanFitness),
               names_to = "Metric", values_to = "Value")

fig1 <- ggplot(df_f1, aes(x = Year, y = Value, color = Scenario)) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 1) +
  theme_minimal(base_size = 14, base_family = "Helvetica") +
  theme(legend.title = element_blank()) +
  labs(x = "Years since burn-in", y = NULL)

fig1


df_f2 <- df %>%
  dplyr::select(Year, Scenario, Population, Heterozygosity, Froh100Kb) %>%
  tidyr::pivot_longer(cols = c(Heterozygosity, Froh100Kb),
               names_to = "Metric", values_to = "Value")

fig2 <- ggplot(df_f2, aes(x = Year, y = Value, color = Scenario)) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 1) +
  theme_minimal(base_size = 14, base_family = "Helvetica") +
  theme(legend.title = element_blank()) +
  labs(x = "Years since burn-in", y = NULL)

fig2
ggsave("het_froh_slim.svg", plot = fig2, width = 8, height = 6, units = "in", device = "svg")
ggsave("het_froh_slim.pdf", plot = fig2, width = 8, height = 6, units = "in", device = cairo_pdf)


df_f3 <- df %>%
  dplyr::select(Year, Scenario, Population, LoadP_high, LoadP_mod, LoadP_low, LoadP_tot) %>%
  tidyr::pivot_longer(cols = starts_with("LoadP_"),
               names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric,
                         levels = c("LoadP_high","LoadP_mod","LoadP_low","LoadP_tot")))

fig3 <- ggplot(df_f3, aes(x = Year, y = Value, color = Scenario)) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
  #coord_cartesian(ylim = c(0, 1)) +
  #scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  theme_minimal(base_size = 14, base_family = "Helvetica") +
  theme(legend.title = element_blank()) +
  labs(x = "Years since burn-in", y = "Potential load (LoadP)")

fig3
ggsave("loadp_slim.svg", plot = fig3, width = 12, height = 6, units = "in", device = "svg")
ggsave("loadp_slim.pdf", plot = fig3, width = 12, height = 6, units = "in", device = cairo_pdf)


df_f4 <- df %>%
  dplyr::select(Year, Scenario, Population, LoadR_high, LoadR_mod, LoadR_low, LoadR_tot) %>%
  tidyr::pivot_longer(cols = starts_with("LoadR_"),
               names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric,
                         levels = c("LoadR_high","LoadR_mod","LoadR_low","LoadR_tot")))

fig4 <- ggplot(df_f4, aes(x = Year, y = Value, color = Scenario)) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
  #coord_cartesian(ylim = c(0, 1)) +
  #scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  theme_minimal(base_size = 14, base_family = "Helvetica") +
  theme(legend.title = element_blank()) +
  labs(x = "Years since burn-in", y = "Realized load (LoadR)")

fig4
ggsave("loadr_slim.svg", plot = fig4, width = 12, height = 6, units = "in", device = "svg")
ggsave("loadr_slim.pdf", plot = fig4, width = 12, height = 6, units = "in", device = cairo_pdf)


# Plotting locations
library(terra)
library(sf)
library(tidyverse)
library(rnaturalearth)
library(ggplot2)
library(dplyr)
library(scico)

hsm_0 <- rast("/scratch/negishi/jeon96/monq/slim/hsm/monq_hsm.tif")

ncol_slim <- 258
nrow_slim <- 291

e <- ext(hsm_0)
x_ras <- c(e$xmin, e$xmax)
y_ras <- c(e$ymin, e$ymax)

hsm_0_df <- as.data.frame(hsm_0, xy = TRUE, na.rm = FALSE)
colnames(hsm_0_df) <- c("x", "y", "suitability")

admin <- ne_countries(scale = "medium", returnclass = "sf")
admin_proj <- st_transform(admin, crs = st_crs(hsm_0))
range <- st_read("/scratch/negishi/jeon96/monq/gea-go/monq_450km_90y_onLand.shp") 
range_proj <- st_transform(range, crs = st_crs(hsm_0))
bbox <- st_bbox(range_proj)

gens <- c("01", "11", "21")

for (gen in gens) {
  loc_file <- sprintf("locations1120%s.txt", gen)
  loc <- read.table(loc_file, header = TRUE)
  loc2 <- loc %>%
  mutate(
    x_norm = (x - 0) / (ncol_slim - 1),
    y_norm = (y - 0) / (nrow_slim - 1),
    x_proj = x_norm * diff(x_ras) + x_ras[1],
    y_proj = y_norm * diff(y_ras) + y_ras[1]  
    )
  pts_proj <- st_as_sf(loc2, coords = c("x_proj", "y_proj"), crs = crs(hsm_0))
  pts_ll <- st_transform(pts_proj, 4326)  

  p <- ggplot() +
    geom_sf(data = admin_proj, fill = gray(0.95), color = NA) +
    geom_sf(data = admin_proj, fill = NA, color = "black", size = 0.2) +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
             ylim = c(bbox["ymin"], bbox["ymax"]),
             expand = TRUE) +

    geom_raster(data = dplyr::filter(hsm_0_df, !is.na(suitability) & suitability > 0),
                aes(x = x, y = y, fill = suitability),
                alpha = 0.6) +
    scale_fill_scico(
      name = "Suitability",
      palette = "vik",
      direction = -1,
      na.value = "transparent",
      limits = c(0, 1)
    ) +

    geom_point(data = loc2,
               aes(x = x_proj, y = y_proj),
               shape = 21,
               fill = "#FDE725",
               color = "black",
               stroke = 0.4,
               size = 2.2,
               alpha = 0.9) +

    theme_minimal(base_size = 14, base_family = "Helvetica") +
    labs(title = paste0("Locations of Gen ", gen))

  print(p)

  out_svg <- sprintf("slim_location_ssp585_gen%s.svg", gen)
  out_pdf <- sprintf("slim_location_ssp585_gen%s.pdf", gen)
  ggsave(out_svg, plot = p, width = 10, height = 10, units = "in", device = "svg")
  ggsave(out_pdf, plot = p, width = 10, height = 10, units = "in", device = cairo_pdf)
}
 

# SSP585
hsm_30 <- rast("/scratch/negishi/jeon96/monq/slim/hsm/monq_hsm_30_ssp585.tif")

e <- ext(hsm_30)
x_ras <- c(e$xmin, e$xmax)
y_ras <- c(e$ymin, e$ymax)

hsm_30_df <- as.data.frame(hsm_30, xy = TRUE, na.rm = FALSE)
colnames(hsm_30_df) <- c("x", "y", "suitability")

gens <- c("31", "41", "51")

for (gen in gens) {
  loc_file <- sprintf("locations1120%s.txt", gen)
  loc <- read.table(loc_file, header = TRUE)
  loc2 <- loc %>%
  mutate(
    x_norm = (x - 0) / (ncol_slim - 1),
    y_norm = (y - 0) / (nrow_slim - 1),
    x_proj = x_norm * diff(x_ras) + x_ras[1],
    y_proj = y_norm * diff(y_ras) + y_ras[1]  
    )
  pts_proj <- st_as_sf(loc2, coords = c("x_proj", "y_proj"), crs = crs(hsm_30))
  pts_ll <- st_transform(pts_proj, 4326)  

  p <- ggplot() +
    geom_sf(data = admin_proj, fill = gray(0.95), color = NA) +
    geom_sf(data = admin_proj, fill = NA, color = "black", size = 0.2) +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
             ylim = c(bbox["ymin"], bbox["ymax"]),
             expand = TRUE) +

    geom_raster(data = dplyr::filter(hsm_30_df, !is.na(suitability) & suitability > 0),
                aes(x = x, y = y, fill = suitability),
                alpha = 0.6) +
    scale_fill_scico(
      name = "Suitability",
      palette = "vik",
      direction = -1,
      na.value = "transparent",
      limits = c(0, 1)
    ) +

    geom_point(data = loc2,
               aes(x = x_proj, y = y_proj),
               shape = 21,
               fill = "#FDE725",
               color = "black",
               stroke = 0.4,
               size = 2.2,
               alpha = 0.9) +

    theme_minimal(base_size = 14, base_family = "Helvetica") +
    labs(title = paste0("Locations of Gen ", gen))

  print(p)

  out_svg <- sprintf("slim_location_ssp585_gen%s.svg", gen)
  out_pdf <- sprintf("slim_location_ssp585_gen%s.pdf", gen)
  ggsave(out_svg, plot = p, width = 10, height = 10, units = "in", device = "svg")
  ggsave(out_pdf, plot = p, width = 10, height = 10, units = "in", device = cairo_pdf)
}


hsm_60 <- rast("/scratch/negishi/jeon96/monq/slim/hsm/monq_hsm_60_ssp585.tif")

e <- ext(hsm_60)
x_ras <- c(e$xmin, e$xmax)
y_ras <- c(e$ymin, e$ymax)

hsm_60_df <- as.data.frame(hsm_60, xy = TRUE, na.rm = FALSE)
colnames(hsm_60_df) <- c("x", "y", "suitability")

gens <- c("61", "71", "81")

for (gen in gens) {
  loc_file <- sprintf("locations1120%s.txt", gen)
  loc <- read.table(loc_file, header = TRUE)
  loc2 <- loc %>%
  mutate(
    x_norm = (x - 0) / (ncol_slim - 1),
    y_norm = (y - 0) / (nrow_slim - 1),
    x_proj = x_norm * diff(x_ras) + x_ras[1],
    y_proj = y_norm * diff(y_ras) + y_ras[1]  
    )
  pts_proj <- st_as_sf(loc2, coords = c("x_proj", "y_proj"), crs = crs(hsm_60))
  pts_ll <- st_transform(pts_proj, 4326)  

  p <- ggplot() +
    geom_sf(data = admin_proj, fill = gray(0.95), color = NA) +
    geom_sf(data = admin_proj, fill = NA, color = "black", size = 0.2) +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
             ylim = c(bbox["ymin"], bbox["ymax"]),
             expand = TRUE) +

    geom_raster(data = dplyr::filter(hsm_60_df, !is.na(suitability) & suitability > 0),
                aes(x = x, y = y, fill = suitability),
                alpha = 0.6) +
    scale_fill_scico(
      name = "Suitability",
      palette = "vik",
      direction = -1,
      na.value = "transparent",
      limits = c(0, 1)
    ) +

    geom_point(data = loc2,
               aes(x = x_proj, y = y_proj),
               shape = 21,
               fill = "#FDE725",
               color = "black",
               stroke = 0.4,
               size = 2.2,
               alpha = 0.9) +

    theme_minimal(base_size = 14, base_family = "Helvetica") +
    labs(title = paste0("Locations of Gen ", gen))

  print(p)

  out_svg <- sprintf("slim_location_ssp585_gen%s.svg", gen)
  out_pdf <- sprintf("slim_location_ssp585_gen%s.pdf", gen)
  ggsave(out_svg, plot = p, width = 10, height = 10, units = "in", device = "svg")
  ggsave(out_pdf, plot = p, width = 10, height = 10, units = "in", device = cairo_pdf)
}


hsm_90 <- rast("/scratch/negishi/jeon96/monq/slim/hsm/monq_hsm_90_ssp585.tif")

e <- ext(hsm_90)
x_ras <- c(e$xmin, e$xmax)
y_ras <- c(e$ymin, e$ymax)

hsm_90_df <- as.data.frame(hsm_90, xy = TRUE, na.rm = FALSE)
colnames(hsm_90_df) <- c("x", "y", "suitability")

gens <- c(91, 101, 111, 121)

for (gen in gens) {
  gen3 <- sprintf("%03d", gen)  
  loc_file <- sprintf("locations112%s.txt", gen3)
  loc <- read.table(loc_file, header = TRUE)
  loc2 <- loc %>%
  mutate(
    x_norm = (x - 0) / (ncol_slim - 1),
    y_norm = (y - 0) / (nrow_slim - 1),
    x_proj = x_norm * diff(x_ras) + x_ras[1],
    y_proj = y_norm * diff(y_ras) + y_ras[1]  
    )
  pts_proj <- st_as_sf(loc2, coords = c("x_proj", "y_proj"), crs = crs(hsm_90))
  pts_ll <- st_transform(pts_proj, 4326)  

  p <- ggplot() +
    geom_sf(data = admin_proj, fill = gray(0.95), color = NA) +
    geom_sf(data = admin_proj, fill = NA, color = "black", size = 0.2) +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
             ylim = c(bbox["ymin"], bbox["ymax"]),
             expand = TRUE) +

    geom_raster(data = dplyr::filter(hsm_90_df, !is.na(suitability) & suitability > 0),
                aes(x = x, y = y, fill = suitability),
                alpha = 0.6) +
    scale_fill_scico(
      name = "Suitability",
      palette = "vik",
      direction = -1,
      na.value = "transparent",
      limits = c(0, 1)
    ) +

    geom_point(data = loc2,
               aes(x = x_proj, y = y_proj),
               shape = 21,
               fill = "#FDE725",
               color = "black",
               stroke = 0.4,
               size = 2.2,
               alpha = 0.9) +

    theme_minimal(base_size = 14, base_family = "Helvetica") +
    labs(title = paste0("Locations of Gen ", gen))

  print(p)

  out_svg <- sprintf("slim_location_ssp585_gen%s.svg", gen)
  out_pdf <- sprintf("slim_location_ssp585_gen%s.pdf", gen)
  ggsave(out_svg, plot = p, width = 10, height = 10, units = "in", device = "svg")
  ggsave(out_pdf, plot = p, width = 10, height = 10, units = "in", device = cairo_pdf)
}


# SSP126
setwd("/scratch/negishi/jeon96/monq/slim/output_K1000_SSP126/")

hsm_0 <- rast("/scratch/negishi/jeon96/monq/slim/hsm/monq_hsm.tif")

ncol_slim <- 258
nrow_slim <- 291

e <- ext(hsm_0)
x_ras <- c(e$xmin, e$xmax)
y_ras <- c(e$ymin, e$ymax)

hsm_0_df <- as.data.frame(hsm_0, xy = TRUE, na.rm = FALSE)
colnames(hsm_0_df) <- c("x", "y", "suitability")

admin <- ne_countries(scale = "medium", returnclass = "sf")
admin_proj <- st_transform(admin, crs = st_crs(hsm_0))
range <- st_read("/scratch/negishi/jeon96/monq/gea-go/monq_450km_90y_onLand.shp") 
range_proj <- st_transform(range, crs = st_crs(hsm_0))
bbox <- st_bbox(range_proj)

gens <- c("01", "11", "21")

for (gen in gens) {
  loc_file <- sprintf("locations1120%s.txt", gen)
  loc <- read.table(loc_file, header = TRUE)
  loc2 <- loc %>%
  mutate(
    x_norm = (x - 0) / (ncol_slim - 1),
    y_norm = (y - 0) / (nrow_slim - 1),
    x_proj = x_norm * diff(x_ras) + x_ras[1],
    y_proj = y_norm * diff(y_ras) + y_ras[1]  
    )
  pts_proj <- st_as_sf(loc2, coords = c("x_proj", "y_proj"), crs = crs(hsm_0))
  pts_ll <- st_transform(pts_proj, 4326)  

  p <- ggplot() +
    geom_sf(data = admin_proj, fill = gray(0.95), color = NA) +
    geom_sf(data = admin_proj, fill = NA, color = "black", size = 0.2) +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
             ylim = c(bbox["ymin"], bbox["ymax"]),
             expand = TRUE) +

    geom_raster(data = dplyr::filter(hsm_0_df, !is.na(suitability) & suitability > 0),
                aes(x = x, y = y, fill = suitability),
                alpha = 0.6) +
    scale_fill_scico(
      name = "Suitability",
      palette = "vik",
      direction = -1,
      na.value = "transparent",
      limits = c(0, 1)
    ) +

    geom_point(data = loc2,
               aes(x = x_proj, y = y_proj),
               shape = 21,
               fill = "#FDE725",
               color = "black",
               stroke = 0.4,
               size = 2.2,
               alpha = 0.9) +

    theme_minimal(base_size = 14, base_family = "Helvetica") +
    labs(title = paste0("Locations of Gen ", gen))

  print(p)

  out_svg <- sprintf("slim_location_ssp126_gen%s.svg", gen)
  out_pdf <- sprintf("slim_location_ssp126_gen%s.pdf", gen)
  ggsave(out_svg, plot = p, width = 10, height = 10, units = "in", device = "svg")
  ggsave(out_pdf, plot = p, width = 10, height = 10, units = "in", device = cairo_pdf)
}


hsm_30 <- rast("/scratch/negishi/jeon96/monq/slim/hsm/monq_hsm_30_ssp126.tif")

e <- ext(hsm_30)
x_ras <- c(e$xmin, e$xmax)
y_ras <- c(e$ymin, e$ymax)

hsm_30_df <- as.data.frame(hsm_30, xy = TRUE, na.rm = FALSE)
colnames(hsm_30_df) <- c("x", "y", "suitability")

gens <- c("31", "41", "51")

for (gen in gens) {
  loc_file <- sprintf("locations1120%s.txt", gen)
  loc <- read.table(loc_file, header = TRUE)
  loc2 <- loc %>%
  mutate(
    x_norm = (x - 0) / (ncol_slim - 1),
    y_norm = (y - 0) / (nrow_slim - 1),
    x_proj = x_norm * diff(x_ras) + x_ras[1],
    y_proj = y_norm * diff(y_ras) + y_ras[1]  
    )
  pts_proj <- st_as_sf(loc2, coords = c("x_proj", "y_proj"), crs = crs(hsm_30))
  pts_ll <- st_transform(pts_proj, 4326)  

  p <- ggplot() +
    geom_sf(data = admin_proj, fill = gray(0.95), color = NA) +
    geom_sf(data = admin_proj, fill = NA, color = "black", size = 0.2) +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
             ylim = c(bbox["ymin"], bbox["ymax"]),
             expand = TRUE) +

    geom_raster(data = dplyr::filter(hsm_30_df, !is.na(suitability) & suitability > 0),
                aes(x = x, y = y, fill = suitability),
                alpha = 0.6) +
    scale_fill_scico(
      name = "Suitability",
      palette = "vik",
      direction = -1,
      na.value = "transparent",
      limits = c(0, 1)
    ) +

    geom_point(data = loc2,
               aes(x = x_proj, y = y_proj),
               shape = 21,
               fill = "#FDE725",
               color = "black",
               stroke = 0.4,
               size = 2.2,
               alpha = 0.9) +

    theme_minimal(base_size = 14, base_family = "Helvetica") +
    labs(title = paste0("Locations of Gen ", gen))

  print(p)

  out_svg <- sprintf("slim_location_ssp126_gen%s.svg", gen)
  out_pdf <- sprintf("slim_location_ssp126_gen%s.pdf", gen)
  ggsave(out_svg, plot = p, width = 10, height = 10, units = "in", device = "svg")
  ggsave(out_pdf, plot = p, width = 10, height = 10, units = "in", device = cairo_pdf)
}


hsm_60 <- rast("/scratch/negishi/jeon96/monq/slim/hsm/monq_hsm_60_ssp126.tif")

e <- ext(hsm_60)
x_ras <- c(e$xmin, e$xmax)
y_ras <- c(e$ymin, e$ymax)

hsm_60_df <- as.data.frame(hsm_60, xy = TRUE, na.rm = FALSE)
colnames(hsm_60_df) <- c("x", "y", "suitability")

gens <- c("61", "71", "81")

for (gen in gens) {
  loc_file <- sprintf("locations1120%s.txt", gen)
  loc <- read.table(loc_file, header = TRUE)
  loc2 <- loc %>%
  mutate(
    x_norm = (x - 0) / (ncol_slim - 1),
    y_norm = (y - 0) / (nrow_slim - 1),
    x_proj = x_norm * diff(x_ras) + x_ras[1],
    y_proj = y_norm * diff(y_ras) + y_ras[1]  
    )
  pts_proj <- st_as_sf(loc2, coords = c("x_proj", "y_proj"), crs = crs(hsm_60))
  pts_ll <- st_transform(pts_proj, 4326)  

  p <- ggplot() +
    geom_sf(data = admin_proj, fill = gray(0.95), color = NA) +
    geom_sf(data = admin_proj, fill = NA, color = "black", size = 0.2) +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
             ylim = c(bbox["ymin"], bbox["ymax"]),
             expand = TRUE) +

    geom_raster(data = dplyr::filter(hsm_60_df, !is.na(suitability) & suitability > 0),
                aes(x = x, y = y, fill = suitability),
                alpha = 0.6) +
    scale_fill_scico(
      name = "Suitability",
      palette = "vik",
      direction = -1,
      na.value = "transparent",
      limits = c(0, 1)
    ) +

    geom_point(data = loc2,
               aes(x = x_proj, y = y_proj),
               shape = 21,
               fill = "#FDE725",
               color = "black",
               stroke = 0.4,
               size = 2.2,
               alpha = 0.9) +

    theme_minimal(base_size = 14, base_family = "Helvetica") +
    labs(title = paste0("Locations of Gen ", gen))

  print(p)

  out_svg <- sprintf("slim_location_ssp126_gen%s.svg", gen)
  out_pdf <- sprintf("slim_location_ssp126_gen%s.pdf", gen)
  ggsave(out_svg, plot = p, width = 10, height = 10, units = "in", device = "svg")
  ggsave(out_pdf, plot = p, width = 10, height = 10, units = "in", device = cairo_pdf)
}


hsm_90 <- rast("/scratch/negishi/jeon96/monq/slim/hsm/monq_hsm_90_ssp126.tif")

e <- ext(hsm_90)
x_ras <- c(e$xmin, e$xmax)
y_ras <- c(e$ymin, e$ymax)

hsm_90_df <- as.data.frame(hsm_90, xy = TRUE, na.rm = FALSE)
colnames(hsm_90_df) <- c("x", "y", "suitability")

gens <- c(91, 101, 111, 121)

for (gen in gens) {
  gen3 <- sprintf("%03d", gen)  
  loc_file <- sprintf("locations112%s.txt", gen3)
  loc <- read.table(loc_file, header = TRUE)
  loc2 <- loc %>%
  mutate(
    x_norm = (x - 0) / (ncol_slim - 1),
    y_norm = (y - 0) / (nrow_slim - 1),
    x_proj = x_norm * diff(x_ras) + x_ras[1],
    y_proj = y_norm * diff(y_ras) + y_ras[1]  
    )
  pts_proj <- st_as_sf(loc2, coords = c("x_proj", "y_proj"), crs = crs(hsm_90))
  pts_ll <- st_transform(pts_proj, 4326)  

  p <- ggplot() +
    geom_sf(data = admin_proj, fill = gray(0.95), color = NA) +
    geom_sf(data = admin_proj, fill = NA, color = "black", size = 0.2) +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
             ylim = c(bbox["ymin"], bbox["ymax"]),
             expand = TRUE) +

    geom_raster(data = dplyr::filter(hsm_90_df, !is.na(suitability) & suitability > 0),
                aes(x = x, y = y, fill = suitability),
                alpha = 0.6) +
    scale_fill_scico(
      name = "Suitability",
      palette = "vik",
      direction = -1,
      na.value = "transparent",
      limits = c(0, 1)
    ) +

    geom_point(data = loc2,
               aes(x = x_proj, y = y_proj),
               shape = 21,
               fill = "#FDE725",
               color = "black",
               stroke = 0.4,
               size = 2.2,
               alpha = 0.9) +

    theme_minimal(base_size = 14, base_family = "Helvetica") +
    labs(title = paste0("Locations of Gen ", gen))

  print(p)

  out_svg <- sprintf("slim_location_ssp126_gen%s.svg", gen)
  out_pdf <- sprintf("slim_location_ssp126_gen%s.pdf", gen)
  ggsave(out_svg, plot = p, width = 10, height = 10, units = "in", device = "svg")
  ggsave(out_pdf, plot = p, width = 10, height = 10, units = "in", device = cairo_pdf)
}


# Genetic rescue scenarios
# Load the CSV file
setwd("/scratch/negishi/jeon96/monq/slim/")
library(tidyverse)

# paths (adjust if needed)
f0000none <- "slim_rescue_0_0_0_none.stats"
f5031p2 <- "slim_rescue_50_3_1_p2.stats"
f10031p2 <- "slim_rescue_100_3_1_p2.stats"
f5031p3 <- "slim_rescue_50_3_1_p3.stats"
f10031p3 <- "slim_rescue_100_3_1_p3.stats"

read_stats <- function(path, scenario){
  readr::read_csv(path, show_col_types = FALSE) %>%
    mutate(scenario = scenario)
}

df <- bind_rows(
  read_stats(f0000none, "0_0_0_0_none"),
  read_stats(f5031p2, "50_3_1_NM"),
  read_stats(f10031p2, "100_3_1_NM"),
  read_stats(f5031p3, "50_3_1_SEMX"),
  read_stats(f10031p3, "100_3_1_SEMX")
)

# Define years since burn-in separately for each scenario/pop if needed:
# (Usually Tick starts at the same burn-in tick in both, but this is safer.)
df <- df %>%
  group_by(scenario) %>%
  mutate(year = tick - min(tick)) %>%
  ungroup()

#df_post <- df %>% filter(Year > 10) # If you want to exclude the very first few ticks:

df_f1 <- df %>%
  dplyr::select(year, scenario, pop_size, mean_fitness, mean_het, mean_het_A, mean_het_G) %>%
  tidyr::pivot_longer(cols = c(pop_size, mean_fitness, mean_het, mean_het_A, mean_het_G),
               names_to = "Metric", values_to = "Value")

df_f1 <- df_f1 %>% mutate(scenario = factor(scenario, levels = c("0_0_0_0_none", "50_3_1_NM", "100_3_1_NM", "50_3_1_SEMX", "100_3_1_SEMX"), labels = c("No rescue", "50 birds, 3 times, from NM", "100 birds, 3 times, from NM", "50 birds, 3 times, from SEMX", "100 birds, 3 times, from SEMX")))

fig1 <- ggplot(df_f1, aes(x = year, y = Value, color = scenario)) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 1, labeller = labeller(Metric = c(pop_size = "Population size", mean_fitness = "Mean fitness", mean_het = "Mean heterozygosity", mean_het_A = "Mean heterozygosity (empirical genome-wide)", mean_het_G = "Mean heterozygosity (empirical outlier SV)"))) +
  theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(legend.title = element_blank()) +
  labs(x = "Years since burn-in", y = NULL)

fig1
ggsave("het_rescue.svg", plot = fig1, width = 9, height = 6, units = "in", device = "svg")
ggsave("het_rescue.pdf", plot = fig1, width = 9, height = 6, units = "in", device = cairo_pdf)

df_f2 <- df %>%
  dplyr::select(year, scenario, LoadP_high, LoadP_mod, LoadP_low, LoadP_tot) %>%
  tidyr::pivot_longer(cols = c(LoadP_high, LoadP_mod, LoadP_low, LoadP_tot), names_to = "Metric", values_to = "Value") %>%
  mutate(LoadType = "Potential load", Class = case_when(str_detect(Metric, "_high") ~ "High", str_detect(Metric, "_mod")  ~ "Moderate", str_detect(Metric, "_low")  ~ "Low", str_detect(Metric, "_tot")  ~ "Total"), Site = "All sites", Class = factor(Class, levels = c("High", "Moderate", "Low", "Total")))

df_f2 <- df_f2 %>% mutate(scenario = factor(scenario, levels = c("0_0_0_0_none", "50_3_1_NM", "100_3_1_NM", "50_3_1_SEMX", "100_3_1_SEMX"), labels = c("No rescue", "50 birds, 3 times, from NM", "100 birds, 3 times, from NM", "50 birds, 3 times, from SEMX", "100 birds, 3 times, from SEMX")), Class = factor(Class, levels = c("High", "Moderate", "Low", "Total")), LoadType = factor(LoadType, levels = c("Potential load")))
df_f2 <- df_f2 %>% mutate(Metric = factor(Metric, levels = c("LoadP_high", "LoadP_mod", "LoadP_low", "LoadP_tot")))

fig2 <- ggplot(df_f2, aes(x = year, y = Value, color = scenario)) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 1, labeller = labeller(Metric = c(LoadP_high = "Potential Load (High, All sites)", LoadP_mod = "Potential Load (Moderate, All sites)", LoadP_low = "Potential Load (Low, All sites)", LoadP_tot = "Potential Load (Total, All sites)"))) +
  theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(legend.title = element_blank()) +
  labs(x = "Years since burn-in", y = NULL)

fig2
ggsave("loadp_all_recue.svg", plot = fig2, width = 9, height = 6, units = "in", device = "svg")
ggsave("loadp_all_recue.pdf", plot = fig2, width = 9, height = 6, units = "in", device = cairo_pdf)

df_f3 <- df %>%
  dplyr::select(year, scenario, LoadP_high_A, LoadP_mod_A, LoadP_low_A, LoadP_tot_A) %>%
  tidyr::pivot_longer(cols = c(LoadP_high_A, LoadP_mod_A, LoadP_low_A, LoadP_tot_A), names_to = "Metric", values_to = "Value") %>%
  mutate(LoadType = "Potential load", Class = case_when(str_detect(Metric, "_high") ~ "High", str_detect(Metric, "_mod")  ~ "Moderate", str_detect(Metric, "_low")  ~ "Low", str_detect(Metric, "_tot")  ~ "Total"), Site = "A sites", Class = factor(Class, levels = c("High", "Moderate", "Low", "Total")))

df_f3 <- df_f3 %>% mutate(scenario = factor(scenario, levels = c("0_0_0_0_none", "50_3_1_NM", "100_3_1_NM", "50_3_1_SEMX", "100_3_1_SEMX"), labels = c("No rescue", "50 birds, 3 times, from NM", "100 birds, 3 times, from NM", "50 birds, 3 times, from SEMX", "100 birds, 3 times, from SEMX")), Class = factor(Class, levels = c("High", "Moderate", "Low", "Total")), LoadType = factor(LoadType, levels = c("Potential load")))
df_f3 <- df_f3 %>% mutate(Metric = factor(Metric, levels = c("LoadP_high_A", "LoadP_mod_A", "LoadP_low_A", "LoadP_tot_A")))

fig3 <- ggplot(df_f3, aes(x = year, y = Value, color = scenario)) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 1, labeller = labeller(Metric = c(LoadP_high_A = "Potential Load (High, Empirical genome-wide)", LoadP_mod_A = "Potential Load (Moderate, Empirical genome-wide)", LoadP_low_A = "Potential Load (Low, Empirical genome-wide)", LoadP_tot_A = "Potential Load (Total, Empirical genome-wide)"))) +
  theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(legend.title = element_blank()) +
  labs(x = "Years since burn-in", y = NULL)

fig3
ggsave("loadp_a_recue.svg", plot = fig3, width = 9, height = 6, units = "in", device = "svg")
ggsave("loadp_a_recue.pdf", plot = fig3, width = 9, height = 6, units = "in", device = cairo_pdf)

df_f4 <- df %>%
  dplyr::select(year, scenario, LoadP_high_G, LoadP_mod_G, LoadP_low_G, LoadP_tot_G) %>%
  tidyr::pivot_longer(cols = c(LoadP_high_G, LoadP_mod_G, LoadP_low_G, LoadP_tot_G), names_to = "Metric", values_to = "Value") %>%
  mutate(LoadType = "Potential load", Class = case_when(str_detect(Metric, "_high") ~ "High", str_detect(Metric, "_mod")  ~ "Moderate", str_detect(Metric, "_low")  ~ "Low", str_detect(Metric, "_tot")  ~ "Total"), Site = "G sites", Class = factor(Class, levels = c("High", "Moderate", "Low", "Total")))

df_f4 <- df_f4 %>% mutate(scenario = factor(scenario, levels = c("0_0_0_0_none", "50_3_1_NM", "100_3_1_NM", "50_3_1_SEMX", "100_3_1_SEMX"), labels = c("No rescue", "50 birds, 3 times, from NM", "100 birds, 3 times, from NM", "50 birds, 3 times, from SEMX", "100 birds, 3 times, from SEMX")), Class = factor(Class, levels = c("High", "Moderate", "Low", "Total")), LoadType = factor(LoadType, levels = c("Potential load")))
df_f4 <- df_f4 %>% mutate(Metric = factor(Metric, levels = c("LoadP_high_G", "LoadP_mod_G", "LoadP_low_G", "LoadP_tot_G")))

fig4 <- ggplot(df_f4, aes(x = year, y = Value, color = scenario)) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 1, labeller = labeller(Metric = c(LoadP_high_G = "Potential Load (High, Empirical outlier SV)", LoadP_mod_G = "Potential Load (Moderate, Empirical outlier SV)", LoadP_low_G = "Potential Load (Low, Empirical outlier SV)", LoadP_tot_G = "Potential Load (Total, Empirical outlier SV)"))) +
  theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(legend.title = element_blank()) +
  labs(x = "Years since burn-in", y = NULL)

fig4
ggsave("loadp_g_recue.svg", plot = fig4, width = 9, height = 6, units = "in", device = "svg")
ggsave("loadp_g_recue.pdf", plot = fig4, width = 9, height = 6, units = "in", device = cairo_pdf)

df_f5 <- df %>%
  dplyr::select(year, scenario, LoadR_high, LoadR_mod, LoadR_low, LoadR_tot) %>%
  tidyr::pivot_longer(cols = c(LoadR_high, LoadR_mod, LoadR_low, LoadR_tot), names_to = "Metric", values_to = "Value") %>%
  mutate(LoadType = "Realized load", Class = case_when(str_detect(Metric, "_high") ~ "High", str_detect(Metric, "_mod")  ~ "Moderate", str_detect(Metric, "_low")  ~ "Low", str_detect(Metric, "_tot")  ~ "Total"), Site = "All sites", Class = factor(Class, levels = c("High", "Moderate", "Low", "Total")))

df_f5 <- df_f5 %>% mutate(scenario = factor(scenario, levels = c("0_0_0_0_none", "50_3_1_NM", "100_3_1_NM", "50_3_1_SEMX", "100_3_1_SEMX"), labels = c("No rescue", "50 birds, 3 times, from NM", "100 birds, 3 times, from NM", "50 birds, 3 times, from SEMX", "100 birds, 3 times, from SEMX")), Class = factor(Class, levels = c("High", "Moderate", "Low", "Total")), LoadType = factor(LoadType, levels = c("Realized load")))
df_f5 <- df_f5 %>% mutate(Metric = factor(Metric, levels = c("LoadR_high", "LoadR_mod", "LoadR_low", "LoadR_tot")))

fig5 <- ggplot(df_f5, aes(x = year, y = Value, color = scenario)) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 1, labeller = labeller(Metric = c(LoadR_high = "Realized Load (High, All sites)", LoadR_mod = "Realized Load (Moderate, All sites)", LoadR_low = "Realized Load (Low, All sites)", LoadR_tot = "Realized Load (Total, All sites)"))) +
  theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(legend.title = element_blank()) +
  labs(x = "Years since burn-in", y = NULL)

fig5
ggsave("loadr_all_recue.svg", plot = fig5, width = 9, height = 6, units = "in", device = "svg")
ggsave("loadr_all_recue.pdf", plot = fig5, width = 9, height = 6, units = "in", device = cairo_pdf)

df_f6 <- df %>%
  dplyr::select(year, scenario, LoadR_high_A, LoadR_mod_A, LoadR_low_A, LoadR_tot_A) %>%
  tidyr::pivot_longer(cols = c(LoadR_high_A, LoadR_mod_A, LoadR_low_A, LoadR_tot_A), names_to = "Metric", values_to = "Value") %>%
  mutate(LoadType = "Realized load", Class = case_when(str_detect(Metric, "_high") ~ "High", str_detect(Metric, "_mod")  ~ "Moderate", str_detect(Metric, "_low")  ~ "Low", str_detect(Metric, "_tot")  ~ "Total"), Site = "A sites", Class = factor(Class, levels = c("High", "Moderate", "Low", "Total")))

df_f6 <- df_f6 %>% mutate(scenario = factor(scenario, levels = c("0_0_0_0_none", "50_3_1_NM", "100_3_1_NM", "50_3_1_SEMX", "100_3_1_SEMX"), labels = c("No rescue", "50 birds, 3 times, from NM", "100 birds, 3 times, from NM", "50 birds, 3 times, from SEMX", "100 birds, 3 times, from SEMX")), Class = factor(Class, levels = c("High", "Moderate", "Low", "Total")), LoadType = factor(LoadType, levels = c("Realized load")))
df_f6 <- df_f6 %>% mutate(Metric = factor(Metric, levels = c("LoadR_high_A", "LoadR_mod_A", "LoadR_low_A", "LoadR_tot_A")))

fig6 <- ggplot(df_f6, aes(x = year, y = Value, color = scenario)) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 1, labeller = labeller(Metric = c(LoadR_high_A = "Realized Load (High, Empirical genome-wide)", LoadR_mod_A = "Realized Load (Moderate, Empirical genome-wide)", LoadR_low_A = "Realized Load (Low, Empirical genome-wide)", LoadR_tot_A = "Realized Load (Total, Empirical genome-wide)"))) +
  theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(legend.title = element_blank()) +
  labs(x = "Years since burn-in", y = NULL)

fig6
ggsave("loadr_a_recue.svg", plot = fig6, width = 9, height = 6, units = "in", device = "svg")
ggsave("loadr_a_recue.pdf", plot = fig6, width = 9, height = 6, units = "in", device = cairo_pdf)

df_f7 <- df %>%
  dplyr::select(year, scenario, LoadR_high_G, LoadR_mod_G, LoadR_low_G, LoadR_tot_G) %>%
  tidyr::pivot_longer(cols = c(LoadR_high_G, LoadR_mod_G, LoadR_low_G, LoadR_tot_G), names_to = "Metric", values_to = "Value") %>%
  mutate(LoadType = "Realized load", Class = case_when(str_detect(Metric, "_high") ~ "High", str_detect(Metric, "_mod")  ~ "Moderate", str_detect(Metric, "_low")  ~ "Low", str_detect(Metric, "_tot")  ~ "Total"), Site = "G sites", Class = factor(Class, levels = c("High", "Moderate", "Low", "Total")))

df_f7 <- df_f7 %>% mutate(scenario = factor(scenario, levels = c("0_0_0_0_none", "50_3_1_NM", "100_3_1_NM", "50_3_1_SEMX", "100_3_1_SEMX"), labels = c("No rescue", "50 birds, 3 times, from NM", "100 birds, 3 times, from NM", "50 birds, 3 times, from SEMX", "100 birds, 3 times, from SEMX")), Class = factor(Class, levels = c("High", "Moderate", "Low", "Total")), LoadType = factor(LoadType, levels = c("Realized load")))
df_f7 <- df_f7 %>% mutate(Metric = factor(Metric, levels = c("LoadR_high_G", "LoadR_mod_G", "LoadR_low_G", "LoadR_tot_G")))

fig7 <- ggplot(df_f7, aes(x = year, y = Value, color = scenario)) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 1, labeller = labeller(Metric = c(LoadR_high_G = "Realized Load (High, Empirical outlier SV)", LoadR_mod_G = "Realized Load (Moderate, Empirical outlier SV)", LoadR_low_G = "Realized Load (Low, Empirical outlier SV)", LoadR_tot_G = "Realized Load (Total, Empirical outlier SV)"))) +
  theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(legend.title = element_blank()) +
  labs(x = "Years since burn-in", y = NULL)

fig7
ggsave("loadr_g_recue.svg", plot = fig7, width = 9, height = 6, units = "in", device = "svg")
ggsave("loadr_g_recue.pdf", plot = fig7, width = 9, height = 6, units = "in", device = cairo_pdf)

df_f8 <- df %>%
  dplyr::select(year, scenario, ROH_100Kb, ROH_500Kb, migrant_prop_genomewide, migrant_prop_empirical) %>%
  tidyr::pivot_longer(cols = c(ROH_100Kb, ROH_500Kb, migrant_prop_genomewide, migrant_prop_empirical),
               names_to = "Metric", values_to = "Value")

df_f8 <- df_f8 %>% mutate(scenario = factor(scenario, levels = c("0_0_0_0_none", "50_3_1_NM", "100_3_1_NM", "50_3_1_SEMX", "100_3_1_SEMX"), labels = c("No rescue", "50 birds, 3 times, from NM", "100 birds, 3 times, from NM", "50 birds, 3 times, from SEMX", "100 birds, 3 times, from SEMX")))
df_f8 <- df_f8 %>% mutate(Metric = factor(Metric, levels = c("ROH_100Kb", "ROH_500Kb", "migrant_prop_genomewide", "migrant_prop_empirical")))

fig8 <- ggplot(df_f8, aes(x = year, y = Value, color = scenario)) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 1, labeller = labeller(Metric = c(ROH_100Kb = "ROH (> 100 kb)", ROH_500Kb = "ROH (> 500 kb)", migrant_prop_genomewide = "migrants' proportion (All sites)", migrant_prop_empirical = "migrants' proportion (Empirical genome-wide)"))) +
  theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(legend.title = element_blank()) +
  labs(x = "Years since burn-in", y = NULL)

fig8
ggsave("roh_migprop_rescue.svg", plot = fig8, width = 9, height = 6, units = "in", device = "svg")
ggsave("roh_migprop_rescue.pdf", plot = fig8, width = 9, height = 6, units = "in", device = cairo_pdf)

df_f9 <- df %>%
  dplyr::select(year, scenario, LoadP_high_G, LoadR_high_G, ROH_1Mb, migrant_prop_adaptive) %>%
  tidyr::pivot_longer(cols = c(LoadP_high_G, LoadR_high_G, ROH_1Mb, migrant_prop_adaptive),
               names_to = "Metric", values_to = "Value")

df_f9 <- df_f9 %>% mutate(scenario = factor(scenario, levels = c("0_0_0_0_none", "50_3_1_NM", "100_3_1_NM", "50_3_1_SEMX", "100_3_1_SEMX"), labels = c("No rescue", "50 birds, 3 times, from NM", "100 birds, 3 times, from NM", "50 birds, 3 times, from SEMX", "100 birds, 3 times, from SEMX")))
df_f9 <- df_f9 %>% mutate(Metric = factor(Metric, levels = c("ROH_1Mb", "LoadP_high_G", "LoadR_high_G", "migrant_prop_adaptive")))

fig9 <- ggplot(df_f9, aes(x = year, y = Value, color = scenario)) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 1, labeller = labeller(Metric = c(LoadP_high_G = "Potential Load (High, Empirical outlier SV)", LoadR_high_G = "Realized Load (High, Empirical outlier SV)", ROH_1Mb = "ROH (> 1 Mb)", migrant_prop_adaptive = "migrants' proportion (Empirical outlier SV)"))) +
  theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(legend.title = element_blank()) +
  labs(x = "Years since burn-in", y = NULL)

fig9
ggsave("main_rescue.svg", plot = fig9, width = 9, height = 6, units = "in", device = "svg")
ggsave("main_rescue.pdf", plot = fig9, width = 9, height = 6, units = "in", device = cairo_pdf)
