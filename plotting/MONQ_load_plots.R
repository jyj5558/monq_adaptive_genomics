library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggpubr)   # stat_compare_means
library(forcats)
library(scales)

xlsx_path <- "C:/Users/jyj55/OneDrive - purdue.edu/DeWoody_Lab/Dissertation/Montezuma quail/popgen/load/load_boxplot_data.xlsx"
sheets <- c(SNP = "SNP", SV = "SV", OUTLIER_SV = "outlier")

read_one <- function(sheet_nm, variant_type_label){
  df <- read_excel(xlsx_path, sheet = sheet_nm) %>%
    janitor::clean_names() %>%   # optional but helpful; remove if you want to keep original
    mutate(variant_type = variant_type_label)
  df
}

df_raw <- bind_rows(
  read_one(sheets["SNP"], "SNP"),
  read_one(sheets["SV"], "SV"),
  read_one(sheets["OUTLIER_SV"], "Outlier SV")
)

count_cols <- c("high_impact", "moderate_impact", "low_impact", "modifier")
prop_cols  <- c("high_prop",  "moderate_prop",  "low_prop",  "modifier_prop")

# Quick check (will warn if not found)
missing_counts <- setdiff(count_cols, names(df_raw))
missing_props  <- setdiff(prop_cols,  names(df_raw))
if(length(missing_counts) > 0) message("Missing count columns: ", paste(missing_counts, collapse=", "))
if(length(missing_props)  > 0) message("Missing prop columns: ",  paste(missing_props,  collapse=", "))

df_raw <- df_raw %>%
  mutate(
    high_prop = if_else(total_variants > 0, high_impact / total_variants, NA_real_),
    moderate_prop = if_else(total_variants > 0, moderate_impact / total_variants, NA_real_),
    low_prop = if_else(total_variants > 0, low_impact / total_variants, NA_real_),
    modifier_prop = if_else(total_variants > 0, modifier / total_variants, NA_real_),
    deleterious_prop = proportion_deleterious
)

id_col <- "sample"
pop_col <- "pop"


# Long format for Panel A (composition) using prop columns
df_prop_long <- df_raw %>%
  select(all_of(c(id_col, pop_col, "variant_type")), all_of(prop_cols)) %>%
  pivot_longer(
    cols = all_of(prop_cols),
    names_to = "impact_class",
    values_to = "prop"
  ) %>%
  mutate(
    impact_class = str_replace(impact_class, "_prop$", ""),
    impact_class = str_to_title(impact_class),
    impact_class = factor(impact_class, levels = c("High","Moderate","Low","Modifier")),
    variant_type = factor(variant_type, levels = c("SNP","SV","Outlier SV"))
  )

variant_fill <- c(
  "SNP" = "#9ecae1",
  "SV" = "#a1d99b",
  "Outlier SV" = "#fdae6b"
)

variant_point <- c(
  "SNP" = "#2171b5",
  "SV" = "#238b45",
  "Outlier SV" = "#d94801"
)

pop_fill <- c(
  "MX" = "#bcbddc",
  "WTX" = "#9ecae1",
  "CTX" = "#fdae6b"
)

pop_point <- c(
  "MX" = "#6a51a3",
  "WTX" = "#2171b5",
  "CTX" = "#d94801"
)

# For Panel A, we want the mean composition per variant_type (optionally per population too)
panelA_dat <- df_prop_long %>%
  group_by(variant_type, impact_class) %>%
  summarise(prop = mean(prop, na.rm = TRUE), .groups = "drop") %>%
  group_by(variant_type) %>%
  mutate(prop = prop / sum(prop, na.rm = TRUE)) %>%  # ensure sums to 1
  ungroup()

# Data for Panel B (individual deleterious burden across variant types)
panelB_dat <- df_raw %>%
  mutate(variant_type = factor(variant_type, levels = c("SNP", "SV", "Outlier SV"))) %>% 
  select(all_of(c(id_col, pop_col, "variant_type", "deleterious_prop")))

# Data for Panel C (population-wise deleterious burden within each variant type)
panelC_dat <- panelB_dat %>%
  mutate(
    pop = factor(pop, levels = c("MX", "WTX", "CTX"))
  )

# Plot for Panel A
pA <- ggplot(panelA_dat, aes(x = variant_type, y = prop, fill = impact_class)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.2) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Mean proportion of variants", fill = "Impact") +
  theme_minimal(base_size = 12)

pA_nolegend <- ggplot(panelA_dat, aes(x = variant_type, y = prop, fill = impact_class)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.2) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Mean proportion of variants", fill = "Impact") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

# Plot for Panel B
y_upper <- max(panelB_dat[["deleterious_prop"]], na.rm = TRUE) * 1.35
pB <- ggplot(panelB_dat, aes(x = variant_type, y = deleterious_prop)) +
  geom_boxplot(aes(fill = variant_type), color = "black", outlier.shape = NA, width = 0.6) +
  scale_fill_manual(values = variant_fill) +
  ggnewscale::new_scale_fill() +
  geom_point(aes(fill = variant_type), shape = 21, color = "black", stroke = 0.5, position = position_jitter(width = 0.15, height = 0), alpha = 0.35, size = 1.5) +
  scale_fill_manual(values = variant_point) +
  labs(x = NULL, y = "Proportion of Deleterious Variants") +
  theme_minimal(base_size = 12) +
  scale_y_continuous(limits = c(0, y_upper), expand = expansion(mult = c(0.02, 0.12))) +
  # Overall test + pairwise tests
  stat_compare_means(method = "kruskal.test", label = "p.format", label.y = y_upper * 0.98, size = 3.3, na.rm = TRUE) +
  stat_compare_means(
    comparisons = list(c("SNP","SV"), c("SNP","Outlier SV"), c("SV","Outlier SV")),
    method = "wilcox.test", label = "p.signif", na.rm = TRUE) +
  coord_cartesian(ylim = c(0, min(1.1, max(panelB_dat$deleterious_prop, na.rm = TRUE) * 1.15)))

pB_nolegend <- ggplot(panelB_dat, aes(x = variant_type, y = deleterious_prop)) +
  geom_boxplot(aes(fill = variant_type), color = "black", outlier.shape = NA, width = 0.6) +
  scale_fill_manual(values = variant_fill) +
  ggnewscale::new_scale_fill() +
  geom_point(aes(fill = variant_type), shape = 21, color = "black", stroke = 0.5, position = position_jitter(width = 0.15, height = 0), alpha = 0.35, size = 1.5) +
  scale_fill_manual(values = variant_point) +
  labs(x = NULL, y = "Proportion of Deleterious Variants") +
  theme_minimal(base_size = 12) +
  scale_y_continuous(limits = c(0, y_upper), expand = expansion(mult = c(0.02, 0.12))) +
  # Overall test + pairwise tests
  stat_compare_means(method = "kruskal.test", label = "p.format", label.y = y_upper * 0.9, size = 3.3,  na.rm = TRUE) +
  stat_compare_means(
    comparisons = list(c("SNP","SV"), c("SNP","Outlier SV"), c("SV","Outlier SV")),
    method = "wilcox.test", label = "p.signif", na.rm = TRUE) +
  coord_cartesian(ylim = c(0, min(1.1, max(panelB_dat$deleterious_prop, na.rm = TRUE) * 1.15))) +
  theme(legend.position = "none")

# Plot for Panel C
pC <- ggplot(panelC_dat, aes(x = pop, y = deleterious_prop)) +
  geom_boxplot(aes(fill = pop), color = "black", outlier.shape = NA, width = 0.6) +
  scale_fill_manual(values = pop_fill) +
  ggnewscale::new_scale_fill() +
  geom_point(aes(fill = pop), shape = 21, color = "black", stroke = 0.5, position = position_jitter(width = 0.15, height = 0), alpha = 0.35, size = 1.5) +
  facet_wrap(~ variant_type, nrow = 1, scale = "fixed") +
  scale_fill_manual(values = pop_point) +
  labs(x = NULL, y = "Proportion of Deleterious Variants") +
  theme_minimal(base_size = 12) +
  scale_y_continuous(limits = c(0, y_upper), expand = expansion(mult = c(0.02, 0.12))) +
  stat_compare_means(method = "kruskal.test", label = "p.format", label.y = y_upper * 0.99, size = 3.3, na.rm = TRUE) +
  stat_compare_means(comparisons = list(c("MX","WTX"), c("MX","CTX"), c("WTX","CTX")), method = "wilcox.test", label = "p.signif", na.rm = TRUE)

pC_nolegend <- ggplot(panelC_dat, aes(x = pop, y = deleterious_prop)) +
  geom_boxplot(aes(fill = pop), color = "black", outlier.shape = NA, width = 0.6) +
  scale_fill_manual(values = pop_fill) +
  ggnewscale::new_scale_fill() +
  geom_point(aes(fill = pop), shape = 21, color = "black", stroke = 0.5, position = position_jitter(width = 0.15, height = 0), alpha = 0.35, size = 1.5) +
  facet_wrap(~ variant_type, nrow = 1, scale = "fixed") +
  scale_fill_manual(values = pop_point) +
  labs(x = NULL, y = "Proportion of Deleterious Variants") +
  theme_minimal(base_size = 12) +
  scale_y_continuous(limits = c(0, y_upper), expand = expansion(mult = c(0.02, 0.12))) +
  stat_compare_means(method = "kruskal.test", label = "p.format", label.y = y_upper * 0.99, size = 3.3, na.rm = TRUE) +
  stat_compare_means(comparisons = list(c("MX","WTX"), c("MX","CTX"), c("WTX","CTX")), method = "wilcox.test", label = "p.signif", na.rm = TRUE) +
  theme(legend.position = "none")

ggsave("stackedbar_impact_composition.svg", pA, width = 4, height = 3)
ggsave("stackedbar_impact_composition_nolegend.svg", pA_nolegend, width = 4, height = 3)
ggsave("deleterious_prop_by_variant.svg", pB, width = 4, height = 3)
ggsave("deleterious_prop_by_variant_nolegend.svg", pB_nolegend, width = 4, height = 3)
ggsave("deleterious_prop_by_pop.svg", pC, width = 8, height = 3)
ggsave("deleterious_prop_by_pop_nolegend.svg", pC_nolegend, width = 8, height = 3)

