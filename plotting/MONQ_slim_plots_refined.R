library(ggplot2)
library(tidyr)
library(dplyr)
library(raster)
library(viridis)
library(ggspatial)
library(terra)
library(sf)
library(scales)
library(ggnewscale)
library(scico)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggsci)   # for Nature palette

# theme for figures
manuscript_theme <- theme_minimal(base_size = 16, base_family = "Helvetica") +
  theme(
    axis.text        = element_text(size = 14, color = "black"),
    axis.title       = element_text(size = 15),
    legend.title     = element_blank(),
    legend.text      = element_text(size = 13),
    strip.text       = element_text(size = 13, face = "bold"),
    strip.background = element_rect(fill = "grey92", color = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey85"),
    plot.title       = element_text(hjust = 0.5, size = 16, face = "bold")
  )

# Nature color palette
scenario_colors <- c(
  "No Rescue"                      = "#1B9E77",  # distinct green
  "50 Birds, 3 Times, from NM"     = "#9ECAE1",  # light blue
  "100 Birds, 3 Times, from NM"    = "#08519C",  # dark blue
  "50 Birds, 3 Times, from SE MX"  = "#FDBE85",  # light orange
  "100 Birds, 3 Times, from SE MX" = "#D94801"   # dark orange
)

scenario_levels <- c("0_0_0_0_none","50_3_1_NM","100_3_1_NM","50_3_1_SEMX","100_3_1_SEMX")
scenario_labels <- names(scenario_colors)

##### SLIM scenarios 
setwd("/Users/natal/Documents/Purdue/MONQ/sim_figures/")
df <- read.csv("/Users/natal/Documents/Purdue/MONQ/sim_figures/stats_ssp585_merged.txt")

# deleterious mutations
df_long1 <- df %>%
  pivot_longer(cols = c(Num_vstrDel, Num_strDel, Num_modDel, Num_wkDel),
               names_to = "type", values_to = "value")

df_long1$type <- factor(df_long1$type,
                        levels = c("Num_vstrDel","Num_strDel","Num_modDel","Num_wkDel"),
                        labels = c("Very Strong","Strong","Moderate","Weak"))

p <- ggplot(df_long1, aes(x = Tick, y = value, color = type)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(
    values = c("Very Strong" = "#d73027","Strong" = "#fc8d59",
               "Moderate" = "#91bfdb","Weak" = "#4575b4"),
    name = NULL) +
  manuscript_theme +
  labs(x = "Future Generation", y = "Number of Deleterious Mutations")

p
ggsave("Num_del_slim_NMA.svg", plot = p, width = 10, height = 6, units = "in", device = "svg")
ggsave("Num_del_slim_NMA.pdf", plot = p, width = 10, height = 6, units = "in", device = cairo_pdf)


# Beneficial mutations
df_long2 <- df %>%
  pivot_longer(cols = c(Num_vstrBen, Num_strBen, Num_modBen, Num_wkBen),
               names_to = "type", values_to = "value")

df_long2$type <- factor(df_long2$type,
                        levels = c("Num_vstrBen","Num_strBen","Num_modBen","Num_wkBen"),
                        labels = c("Very Strong","Strong","Moderate","Weak"))

p <- ggplot(df_long2, aes(x = Tick, y = value, color = type)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(
    values = c("Very Strong" = "#d73027","Strong" = "#fc8d59",
               "Moderate" = "#91bfdb","Weak" = "#4575b4"),
    name = NULL) +
  manuscript_theme +
  labs(x = "Future Generation", y = "Number of Beneficial Mutations")

p
ggsave("Num_ben_slim_NMA.svg", plot = p, width = 10, height = 6, units = "in", device = "svg")
ggsave("Num_ben_slim_NMA.pdf", plot = p, width = 10, height = 6, units = "in", device = cairo_pdf)


# Population stats (PopSize, MeanFitness, InbreedingLoad)
df_long3 <- df %>%
  pivot_longer(cols = c(PopSize, MeanFitness, InbreedingLoad),
               names_to = "type", values_to = "value") %>%
  filter(Tick > 10)

df_long3$type <- factor(df_long3$type,
                        levels = c("PopSize","MeanFitness","InbreedingLoad"),
                        labels = c("Population Size","Mean Fitness","Inbreeding Load"))

p <- ggplot(df_long3, aes(x = Tick, y = value, color = type)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ type, scales = "free_y") +
  manuscript_theme +
  theme(legend.position = "none") +
  labs(x = "Future Generation", y = NULL)

p
ggsave("Pop_stats_slim_NMA.svg", plot = p, width = 15, height = 6, units = "in", device = "svg")
ggsave("Pop_stats_slim_NMA.pdf", plot = p, width = 15, height = 6, units = "in", device = cairo_pdf)


# Heterozygosity & FROH
df_long4 <- df %>%
  pivot_longer(cols = c(Heterozygosity, Froh100Kb),
               names_to = "type", values_to = "value") %>%
  filter(Tick > 10)

df_long4$type <- factor(df_long4$type,
                        levels = c("Heterozygosity","Froh100Kb"),
                        labels = c("Heterozygosity","Froh100Kb"))

p <- ggplot(df_long4, aes(x = Tick, y = value, color = type)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ type, scales = "free_y") +
  manuscript_theme +
  theme(legend.position = "none") +
  labs(x = "Future Generation", y = NULL)

p


# Potential Load
df_long5 <- df %>%
  pivot_longer(cols = c(LoadP_high, LoadP_mod, LoadP_low, LoadP_tot),
               names_to = "type", values_to = "value") %>%
  filter(Tick > 10)

df_long5$type <- factor(df_long5$type,
                        levels = c("LoadP_high","LoadP_mod","LoadP_low","LoadP_tot"),
                        labels = c("High","Moderate","Low","Total"))

p <- ggplot(df_long5, aes(x = Tick, y = value, color = type)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ type, scales = "free_y") +
  manuscript_theme +
  theme(legend.position = "none") +
  labs(x = "Future Generation", y = "Potential Load")

p

# Realized Load
df_long6 <- df %>%
  pivot_longer(cols = c(LoadR_high, LoadR_mod, LoadR_low, LoadR_tot),
               names_to = "type", values_to = "value") %>%
  filter(Tick > 10)

df_long6$type <- factor(df_long6$type,
                        levels = c("LoadR_high","LoadR_mod","LoadR_low","LoadR_tot"),
                        labels = c("High","Moderate","Low","Total"))

p <- ggplot(df_long6, aes(x = Tick, y = value, color = type)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ type, scales = "free_y") +
  manuscript_theme +
  theme(legend.position = "none") +
  labs(x = "Future Generation", y = "Realized Load")

p

#### SSP126 vs SSP585 
library(tidyverse)

f126 <- "stats_ssp126_merged.txt"
f585 <- "stats_ssp585_merged.txt"

read_stats <- function(path, scenario){
  readr::read_csv(path, show_col_types = FALSE) %>%
    mutate(Scenario = scenario)
}

df <- bind_rows(
  read_stats(f126, "SSP126"),
  read_stats(f585, "SSP585")
) %>%
  group_by(Scenario, Population) %>%
  mutate(Year = Tick - min(Tick)) %>%
  ungroup()

# Define manual colors
scenario_cols <- c(
  "SSP126" = "#08519C",
  "SSP585" = "#D94801"
)

df <- bind_rows(
  read_stats(f126, "SSP126"),
  read_stats(f585, "SSP585")
) %>%
  group_by(Scenario, Population) %>%
  mutate(Year = Tick - min(Tick)) %>%
  ungroup()


# pop size and mean fitness
df_f1 <- df %>%
  dplyr::select(Year, Scenario, Population, PopSize, MeanFitness) %>%
  tidyr::pivot_longer(cols = c(PopSize, MeanFitness),
                      names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric,
                         levels = c("PopSize","MeanFitness"),
                         labels = c("Population Size","Mean Fitness")))

fig1 <- ggplot(df_f1, aes(x = Year, y = Value, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 1) +
  scale_color_manual(values = scenario_cols) +
  manuscript_theme +
  labs(x = "Years Since Burn-in", y = NULL)

fig1


# Het & FROH
df_f2 <- df %>%
  dplyr::select(Year, Scenario, Population, Heterozygosity, Froh100Kb) %>%
  tidyr::pivot_longer(cols = c(Heterozygosity, Froh100Kb),
                      names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric,
                         levels = c("Heterozygosity","Froh100Kb"),
                         labels = c("Heterozygosity","Froh100Kb")))

fig2 <- ggplot(df_f2, aes(x = Year, y = Value, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 1) +
  scale_color_manual(values = scenario_cols) +
  manuscript_theme +
  labs(x = "Years Since Burn-in", y = NULL)

fig2
ggsave("Fig_S15_NMA.svg", plot = fig2, width = 9, height = 6, units = "in", device = "svg")
ggsave("Fig_S15_NMA.pdf", plot = fig2, width = 9, height = 6, units = "in", device = cairo_pdf)


# Potential load
df_f3 <- df %>%
  dplyr::select(Year, Scenario, Population, LoadP_high, LoadP_mod, LoadP_low, LoadP_tot) %>%
  tidyr::pivot_longer(cols = starts_with("LoadP_"),
                      names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric,
                         levels = c("LoadP_high","LoadP_mod","LoadP_low","LoadP_tot")))

loadp_labels <- c(
  LoadP_high = "Potential Load (High)",
  LoadP_mod  = "Potential Load (Moderate)",
  LoadP_low  = "Potential Load (Low)",
  LoadP_tot  = "Potential Load (Total)"
)

fig3 <- ggplot(df_f3, aes(x = Year, y = Value, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 2,
             labeller = labeller(Metric = loadp_labels)) +
  scale_color_manual(values = scenario_cols) +
  manuscript_theme +
  labs(x = "Years Since Burn-in", y = "Potential Load")

fig3
ggsave("Fig_S16A_NMA.svg", plot = fig3, width = 9, height = 6, units = "in", device = "svg")
ggsave("Fig_S16A_NMA.pdf", plot = fig3, width = 9, height = 6, units = "in", device = cairo_pdf)


# Realized load
df_f4 <- df %>%
  dplyr::select(Year, Scenario, Population, LoadR_high, LoadR_mod, LoadR_low, LoadR_tot) %>%
  tidyr::pivot_longer(cols = starts_with("LoadR_"),
                      names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric,
                         levels = c("LoadR_high","LoadR_mod","LoadR_low","LoadR_tot")))

loadr_labels <- c(
  LoadR_high = "Realized Load (High)",
  LoadR_mod  = "Realized Load (Moderate)",
  LoadR_low  = "Realized Load (Low)",
  LoadR_tot  = "Realized Load (Total)"
)

fig4 <- ggplot(df_f4, aes(x = Year, y = Value, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 2,
             labeller = labeller(Metric = loadr_labels)) +
  scale_color_manual(values = scenario_cols) +
  manuscript_theme +
  labs(x = "Years Since Burn-in", y = "Realized Load")

fig4
ggsave("Fig_S16B_NMA.svg", plot = fig4, width = 9, height = 6, units = "in", device = "svg")
ggsave("Fig_S16B_NMA.pdf", plot = fig4, width = 9, height = 6, units = "in", device = cairo_pdf)


##### genetic rescue scenarios
setwd("/Users/natal/Documents/Purdue/MONQ/sim_figures/")
library(tidyverse)

f0000none <- "slim_rescue_0_0_0_none.stats"
f5031p2   <- "slim_rescue_50_3_1_p2.stats"
f10031p2  <- "slim_rescue_100_3_1_p2.stats"
f5031p3   <- "slim_rescue_50_3_1_p3.stats"
f10031p3  <- "slim_rescue_100_3_1_p3.stats"

read_stats <- function(path, scenario){
  readr::read_csv(path, show_col_types = FALSE) %>%
    mutate(scenario = scenario)
}

df <- bind_rows(
  read_stats(f0000none, "0_0_0_0_none"),
  read_stats(f5031p2,   "50_3_1_NM"),
  read_stats(f10031p2,  "100_3_1_NM"),
  read_stats(f5031p3,   "50_3_1_SEMX"),
  read_stats(f10031p3,  "100_3_1_SEMX")
) %>%
  group_by(scenario) %>%
  mutate(year = tick - min(tick)) %>%
  ungroup()

relabel_scenario <- function(x){
  factor(x,
         levels = scenario_levels,
         labels = scenario_labels)
}

# Pop size, fitness, heterozygosity
df_f1 <- df %>%
  dplyr::select(year, scenario, pop_size, mean_fitness, mean_het, mean_het_A, mean_het_G) %>%
  tidyr::pivot_longer(cols = c(pop_size, mean_fitness, mean_het, mean_het_A, mean_het_G),
                      names_to = "Metric", values_to = "Value") %>%
  mutate(scenario = relabel_scenario(scenario))

het_labels <- c(
  pop_size      = "Population Size",
  mean_fitness  = "Mean Fitness",
  mean_het      = "Mean Heterozygosity\nAll Variant Sites",
  mean_het_A    = "Mean Heterozygosity\nGenome-Wide Variant Sites",
  mean_het_G    = "Mean Heterozygosity\nOutlier SV Sites"
)

fig1 <- ggplot(df_f1, aes(x = year, y = Value, color = scenario)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = scenario_colors) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 1,
             labeller = labeller(Metric = het_labels)) +
  manuscript_theme +
  labs(x = "Years Since Burn-in", y = NULL)

fig1
ggsave("Fig_S17_NMA.svg", plot = fig1, width = 15, height = 14, units = "in", device = "svg")
ggsave("Fig_S17_NMA.pdf", plot = fig1, width = 15, height = 14, units = "in", device = cairo_pdf)

# Potential load (all variant sites)
# Two-row strip: line 1 = metric, line 2 = site type
df_f2 <- df %>%
  dplyr::select(year, scenario, LoadP_high, LoadP_mod, LoadP_low, LoadP_tot) %>%
  tidyr::pivot_longer(cols = c(LoadP_high, LoadP_mod, LoadP_low, LoadP_tot),
                      names_to = "Metric", values_to = "Value") %>%
  mutate(scenario = relabel_scenario(scenario),
         Metric   = factor(Metric, levels = c("LoadP_high","LoadP_mod","LoadP_low","LoadP_tot")))

loadp_all_labels <- c(
  LoadP_high = "Potential Load (High)\nAll Variant Sites",
  LoadP_mod  = "Potential Load (Moderate)\nAll Variant Sites",
  LoadP_low  = "Potential Load (Low)\nAll Variant Sites",
  LoadP_tot  = "Potential Load (Total)\nAll Variant Sites"
)

fig2 <- ggplot(df_f2, aes(x = year, y = Value, color = scenario)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = scenario_colors) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 2,
             labeller = labeller(Metric = loadp_all_labels)) +
  manuscript_theme +
  labs(x = "Years Since Burn-in", y = NULL)

fig2
ggsave("Fig_S19_NMA.svg", plot = fig2, width = 15, height = 8, units = "in", device = "svg")
ggsave("Fig_S19_NMA.pdf", plot = fig2, width = 15, height = 8, units = "in", device = cairo_pdf)


# Potential Load (genome-wide variant sites)
df_f3 <- df %>%
  dplyr::select(year, scenario, LoadP_high_A, LoadP_mod_A, LoadP_low_A, LoadP_tot_A) %>%
  tidyr::pivot_longer(cols = c(LoadP_high_A, LoadP_mod_A, LoadP_low_A, LoadP_tot_A),
                      names_to = "Metric", values_to = "Value") %>%
  mutate(scenario = relabel_scenario(scenario),
         Metric   = factor(Metric, levels = c("LoadP_high_A","LoadP_mod_A","LoadP_low_A","LoadP_tot_A")))

loadp_a_labels <- c(
  LoadP_high_A = "Potential Load (High)\nGenome-Wide Variant Sites",
  LoadP_mod_A  = "Potential Load (Moderate)\nGenome-Wide Variant Sites",
  LoadP_low_A  = "Potential Load (Low)\nGenome-Wide Variant Sites",
  LoadP_tot_A  = "Potential Load (Total)\nGenome-Wide Variant Sites"
)

fig3 <- ggplot(df_f3, aes(x = year, y = Value, color = scenario)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = scenario_colors) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 2,
             labeller = labeller(Metric = loadp_a_labels)) +
  manuscript_theme +
  labs(x = "Years Since Burn-in", y = NULL)

fig3
ggsave("Fig_S20_NMA.svg", plot = fig3, width = 15, height = 8, units = "in", device = "svg")
ggsave("Fig_S20_NMA.pdf", plot = fig3, width = 15, height = 8, units = "in", device = cairo_pdf)


# Potential Load (outlier SV sites)
df_f4 <- df %>%
  dplyr::select(year, scenario, LoadP_high_G, LoadP_mod_G, LoadP_low_G, LoadP_tot_G) %>%
  tidyr::pivot_longer(cols = c(LoadP_high_G, LoadP_mod_G, LoadP_low_G, LoadP_tot_G),
                      names_to = "Metric", values_to = "Value") %>%
  mutate(scenario = relabel_scenario(scenario),
         Metric   = factor(Metric, levels = c("LoadP_high_G","LoadP_mod_G","LoadP_low_G","LoadP_tot_G")))

loadp_g_labels <- c(
  LoadP_high_G = "Potential Load (High)\nOutlier SV Sites",
  LoadP_mod_G  = "Potential Load (Moderate)\nOutlier SV Sites",
  LoadP_low_G  = "Potential Load (Low)\nOutlier SV Sites",
  LoadP_tot_G  = "Potential Load (Total)\nOutlier SV Sites"
)

fig4 <- ggplot(df_f4, aes(x = year, y = Value, color = scenario)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = scenario_colors) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 2,
             labeller = labeller(Metric = loadp_g_labels)) +
  manuscript_theme +
  labs(x = "Years Since Burn-in", y = NULL)

fig4
ggsave("loadp_g_rescue_NMA.svg", plot = fig4, width = 15, height = 8, units = "in", device = "svg")
ggsave("loadp_g_rescue_NMA.pdf", plot = fig4, width = 15, height = 8, units = "in", device = cairo_pdf)


# Realized Load (all variant sites)
df_f5 <- df %>%
  dplyr::select(year, scenario, LoadR_high, LoadR_mod, LoadR_low, LoadR_tot) %>%
  tidyr::pivot_longer(cols = c(LoadR_high, LoadR_mod, LoadR_low, LoadR_tot),
                      names_to = "Metric", values_to = "Value") %>%
  mutate(scenario = relabel_scenario(scenario),
         Metric   = factor(Metric, levels = c("LoadR_high","LoadR_mod","LoadR_low","LoadR_tot")))

loadr_all_labels <- c(
  LoadR_high = "Realized Load (High)\nAll Variant Sites",
  LoadR_mod  = "Realized Load (Moderate)\nAll Variant Sites",
  LoadR_low  = "Realized Load (Low)\nAll Variant Sites",
  LoadR_tot  = "Realized Load (Total)\nAll Variant Sites"
)

fig5 <- ggplot(df_f5, aes(x = year, y = Value, color = scenario)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = scenario_colors) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 2,
             labeller = labeller(Metric = loadr_all_labels)) +
  manuscript_theme +
  labs(x = "Years Since Burn-in", y = NULL)

fig5
ggsave("Fig_S21_NMA.svg", plot = fig5, width = 15, height = 8, units = "in", device = "svg")
ggsave("Fig_S21_NMA.pdf", plot = fig5, width = 15, height = 8, units = "in", device = cairo_pdf)


# Realized Load (genome-wide variant sites)
df_f6 <- df %>%
  dplyr::select(year, scenario, LoadR_high_A, LoadR_mod_A, LoadR_low_A, LoadR_tot_A) %>%
  tidyr::pivot_longer(cols = c(LoadR_high_A, LoadR_mod_A, LoadR_low_A, LoadR_tot_A),
                      names_to = "Metric", values_to = "Value") %>%
  mutate(scenario = relabel_scenario(scenario),
         Metric   = factor(Metric, levels = c("LoadR_high_A","LoadR_mod_A","LoadR_low_A","LoadR_tot_A")))

loadr_a_labels <- c(
  LoadR_high_A = "Realized Load (High)\nGenome-Wide Variant Sites",
  LoadR_mod_A  = "Realized Load (Moderate)\nGenome-Wide Variant Sites",
  LoadR_low_A  = "Realized Load (Low)\nGenome-Wide Variant Sites",
  LoadR_tot_A  = "Realized Load (Total)\nGenome-Wide Variant Sites"
)

fig6 <- ggplot(df_f6, aes(x = year, y = Value, color = scenario)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = scenario_colors) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 2,
             labeller = labeller(Metric = loadr_a_labels)) +
  manuscript_theme +
  labs(x = "Years Since Burn-in", y = NULL)

fig6
ggsave("Fig_S22_NMA.svg", plot = fig6, width = 15, height = 8, units = "in", device = "svg")
ggsave("Fig_S22_NMA.pdf", plot = fig6, width = 15, height = 8, units = "in", device = cairo_pdf)


# Realized Load (outlier SV sites)
df_f7 <- df %>%
  dplyr::select(year, scenario, LoadR_high_G, LoadR_mod_G, LoadR_low_G, LoadR_tot_G) %>%
  tidyr::pivot_longer(cols = c(LoadR_high_G, LoadR_mod_G, LoadR_low_G, LoadR_tot_G),
                      names_to = "Metric", values_to = "Value") %>%
  mutate(scenario = relabel_scenario(scenario),
         Metric   = factor(Metric, levels = c("LoadR_high_G","LoadR_mod_G","LoadR_low_G","LoadR_tot_G")))

loadr_g_labels <- c(
  LoadR_high_G = "Realized Load (High)\nOutlier SV Sites",
  LoadR_mod_G  = "Realized Load (Moderate)\nOutlier SV Sites",
  LoadR_low_G  = "Realized Load (Low)\nOutlier SV Sites",
  LoadR_tot_G  = "Realized Load (Total)\nOutlier SV Sites"
)

fig7 <- ggplot(df_f7, aes(x = year, y = Value, color = scenario)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = scenario_colors) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 2,
             labeller = labeller(Metric = loadr_g_labels)) +
  manuscript_theme +
  labs(x = "Years Since Burn-in", y = NULL)

fig7
ggsave("loadr_g_rescue_NMA.svg", plot = fig7, width = 15, height = 8, units = "in", device = "svg")
ggsave("loadr_g_rescue_NMA.pdf", plot = fig7, width = 15, height = 8, units = "in", device = cairo_pdf)


# ROH & Donor Ancestry
df_f8 <- df %>%
  dplyr::select(year, scenario, ROH_100Kb, ROH_500Kb, migrant_prop_genomewide, migrant_prop_empirical) %>%
  tidyr::pivot_longer(cols = c(ROH_100Kb, ROH_500Kb, migrant_prop_genomewide, migrant_prop_empirical),
                      names_to = "Metric", values_to = "Value") %>%
  mutate(scenario = relabel_scenario(scenario),
         Metric   = factor(Metric, levels = c("ROH_100Kb","ROH_500Kb",
                                              "migrant_prop_genomewide","migrant_prop_empirical")))

roh_labels <- c(
  ROH_100Kb               = "ROH (>100 kb)",
  ROH_500Kb               = "ROH (>500 kb)",
  migrant_prop_genomewide = "Donor Ancestry\nAll Variant Sites",
  migrant_prop_empirical  = "Donor Ancestry\nGenome-Wide Variant Sites"
)

fig8 <- ggplot(df_f8, aes(x = year, y = Value, color = scenario)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = scenario_colors) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 2,
             labeller = labeller(Metric = roh_labels)) +
  manuscript_theme +
  labs(x = "Years Since Burn-in", y = NULL)

fig8
ggsave("Fig_S18_NMA.svg", plot = fig8, width = 15, height = 8, units = "in", device = "svg")
ggsave("Fig_S18_NMA.pdf", plot = fig8, width = 15, height = 8, units = "in", device = cairo_pdf)


# rescue summary figure
df_f9 <- df %>%
  dplyr::select(year, scenario, LoadP_high_G, LoadR_high_G, ROH_1Mb, migrant_prop_adaptive) %>%
  tidyr::pivot_longer(cols = c(LoadP_high_G, LoadR_high_G, ROH_1Mb, migrant_prop_adaptive),
                      names_to = "Metric", values_to = "Value") %>%
  mutate(
    scenario = relabel_scenario(scenario),
    Metric   = factor(Metric,
                      levels = c("ROH_1Mb",
                                 "migrant_prop_adaptive",
                                 "LoadP_high_G",
                                 "LoadR_high_G"))
  )

main_labels <- c(
  ROH_1Mb               = "ROH (>1 Mb)",
  migrant_prop_adaptive = "Donor Ancestry\nOutlier SV Sites",
  LoadP_high_G          = "Potential Load (High)\nOutlier SV Sites",
  LoadR_high_G          = "Realized Load (High)\nOutlier SV Sites"
)

fig9 <- ggplot(df_f9, aes(x = year, y = Value, color = scenario)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = scenario_colors) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 2,
             labeller = labeller(Metric = main_labels)) +
  manuscript_theme +
  labs(x = "Years Since Burn-in", y = NULL)

fig9
ggsave("Fig6_NMA.svg", plot = fig9, width = 15, height = 8, units = "in", device = "svg")
ggsave("Fig6_NMA.pdf", plot = fig9, width = 15, height = 8, units = "in", device = cairo_pdf)
