library(tidyverse)
source("functions/pivot_filter.R")
source("functions/wrangling.R")
source("functions/plotting.R")
Metadata <- read_csv(file = "data/Metadata.csv", show_col_types = FALSE)

# random effect ####
## alpha S#####
alpha_S_comparisons <- .wrangle_comparisons(
  scale = "alpha",
  metric = "S",
  study_level = TRUE
)

### alpha S plot comparison at study level####
studylevel <- alpha_S_comparisons$Comparisons |>
  summarise(
    mean = mean(Ratio),
    lower_ci = mean - qt(0.975, df = n() - 1) * (sd(Ratio) / sqrt(n())),
    upper_ci = mean + qt(0.975, df = n() - 1) * (sd(Ratio) / sqrt(n())),
    .by = c(Dataset_id, Comparison)
  ) |>
  left_join(Metadata, by = "Dataset_id")

studylevel <- studylevel |>
  mutate(Comparison = recode(Comparison,
    "Agriculture/Natural vegetation" = "Agriculture",
    "Forestry/Natural vegetation" = "Forestry",
    "Urban/Natural vegetation" = "Urban"
  ) |>
    fct_relevel("Agriculture", "Forestry", "Urban"))

lat_mean <- conditional_effects(lat_alpha_S_model4,
  effects = "Latitude"
)$Latitude

lat_mean_lu <- conditional_effects(lat_alpha_S_model5,
  effects = "Latitude:Comparison"
)$"Latitude:Comparison"

lat_mean_lu <- lat_mean_lu |>
  mutate(Comparison = recode(Comparison,
    "Agriculture" = "Agriculture",
    "Forestry" = "Forestry",
    "Urban" = "Urban"
  ) |>
    fct_relevel("Agriculture", "Forestry", "Urban"))

alpha_latitude_S <- studylevel |>
  mutate(Latitude = abs(Latitude), Dataset_id = as.character(Dataset_id)) |>
  .plot_comparison_latitude(lab_y = "α", subtitle = "S")

studylevel <- studylevel |>
  filter(Taxa %in% c("Algae", "Fish", "Macroinvertebrate")) |>
  mutate(
    Taxa = factor(
      Taxa,
      levels = c("Algae", "Fish", "Macroinvertebrate"),
      labels = c("Algae", "Fish", "Macroinvertebrates")
    )
  )

# set.seed(111)
# taxa_alpha_S1 <- brm(mean ~ Taxa + (1 | Dataset_id),
#   data = studylevel,
#   family = "gaussian",
#   iter = 2000,
#   warmup = 1000,
#   cores = 4,
#   chains = 4
# )
#
# save(taxa_alpha_S1, file = "taxa_alpha_S1.Rdata")
#
# set.seed(111)
# taxa_alpha_S2 <- brm(mean ~ Taxa * Comparison + (1 | Dataset_id),
#   data = studylevel,
#   family = "gaussian",
#   iter = 2000,
#   warmup = 1000,
#   cores = 4,
#   chains = 4
# )

save(taxa_alpha_S2, file = "taxa_alpha_S2.Rdata")

taxa_mean <- conditional_effects(taxa_alpha_S1, effects = "Taxa")$Taxa

taxa_mean_lu <- conditional_effects(taxa_alpha_S2, effects = "Taxa:Comparison")$`Taxa:Comparison`

unique(taxa_mean_lu$Comparison)

taxa_mean_lu <- taxa_mean_lu |>
  mutate(Comparison = recode(Comparison,
    "Agriculture/Natural vegetation" = "Agriculture",
    "Forestry/Natural vegetation" = "Forestry",
    "Urban/Natural vegetation" = "Urban"
  ) |>
    fct_relevel("Agriculture", "Forestry", "Urban"))

alpha_taxa_S <- .plot_comparison_taxa(
  studylevel = studylevel |> mutate(Dataset_id = as.character(Dataset_id)),
  taxamean = taxamean,
  lab_y = "α",
  subtitle = "S"
)

## alpha spie ####
alpha_Spie_comparisons <- .wrangle_comparisons(
  scale = "alpha",
  metric = "Spie",
  study_level = TRUE
)

### alpha Spie plot comparison at study level####
studylevel <- alpha_Spie_comparisons$Comparisons |>
  summarise(
    mean = mean(Ratio),
    lower_ci = mean - qt(0.975, df = n() - 1) * (sd(Ratio) / sqrt(n())),
    upper_ci = mean + qt(0.975, df = n() - 1) * (sd(Ratio) / sqrt(n())),
    .by = c(Dataset_id, Comparison)
  ) |>
  left_join(Metadata, by = "Dataset_id")

studylevel <- studylevel |>
  mutate(Comparison = recode(Comparison,
    "Agriculture/Natural vegetation" = "Agriculture",
    "Forestry/Natural vegetation" = "Forestry",
    "Urban/Natural vegetation" = "Urban"
  ) |>
    fct_relevel("Agriculture", "Forestry", "Urban"))

lat_mean <- conditional_effects(lat_alpha_Spie_model4,
  effects = "Latitude"
)$Latitude

lat_mean_lu <- conditional_effects(lat_alpha_Spie_model5,
  effects = "Latitude:Comparison"
)$"Latitude:Comparison"

lat_mean_lu <- lat_mean_lu |>
  mutate(Comparison = recode(Comparison,
    "Agriculture" = "Agriculture",
    "Forestry" = "Forestry",
    "Urban" = "Urban"
  ) |>
    fct_relevel("Agriculture", "Forestry", "Urban"))

alpha_latitude_Spie <- studylevel |>
  mutate(Latitude = abs(Latitude), Dataset_id = as.character(Dataset_id)) |>
  .plot_comparison_latitude(subtitle = expression(S[PIE]))

# Taxa
studylevel <- studylevel |>
  filter(Taxa %in% c("Algae", "Fish", "Macroinvertebrate")) |>
  mutate(
    Taxa = factor(
      Taxa,
      levels = c("Algae", "Fish", "Macroinvertebrate"),
      labels = c("Algae", "Fish", "Macroinvertebrates")
    )
  )

# set.seed(111)
# taxa_alpha_Spie1 <- brm(mean ~ Taxa + (1 | Dataset_id),
#   data = studylevel,
#   family = "gaussian",
#   iter = 2000,
#   warmup = 1000,
#   cores = 4,
#   chains = 4
# )
#
# save(taxa_alpha_Spie1, file = "taxa_alpha_Spie1.Rdata")
#
# set.seed(111)
# taxa_alpha_Spie2 <- brm(mean ~ Taxa * Comparison + (1 | Dataset_id),
#   data = studylevel,
#   family = "gaussian",
#   iter = 2000,
#   warmup = 1000,
#   cores = 4,
#   chains = 4
# )
#
# save(taxa_alpha_Spie2, file = "taxa_alpha_Spie2.Rdata")

taxa_mean <- conditional_effects(taxa_alpha_Spie1, effects = "Taxa")$Taxa

taxa_mean_lu <- conditional_effects(taxa_alpha_Spie2, effects = "Taxa:Comparison")$`Taxa:Comparison`

unique(taxa_mean_lu$Comparison)

taxa_mean_lu <- taxa_mean_lu |>
  mutate(Comparison = recode(Comparison,
    "Agriculture/Natural vegetation" = "Agriculture",
    "Forestry/Natural vegetation" = "Forestry",
    "Urban/Natural vegetation" = "Urban"
  ) |>
    fct_relevel("Agriculture", "Forestry", "Urban"))

alpha_taxa_Spie <- .plot_comparison_taxa(
  studylevel = studylevel |> mutate(Dataset_id = as.character(Dataset_id)),
  taxamean = taxamean,
  subtitle = expression(S[PIE])
)


## gamma S ####
gamma_S_comparison <- .wrangle_comparisons(
  scale = "gamma",
  metric = "S",
  study_level = TRUE
)

### gamma S plot comparison at study level####
# Latitude
studylevel <- gamma_S_comparison$Comparisons |>
  summarise(
    mean = mean(Ratio),
    lower_ci = mean -
      stats::qt(0.975, df = n() - 1) * (stats::sd(Ratio) / sqrt(n())),
    upper_ci = mean +
      stats::qt(0.975, df = n() - 1) * (stats::sd(Ratio) / sqrt(n())),
    .by = c(Dataset_id, Comparison)
  ) |>
  left_join(Metadata, by = join_by(Dataset_id))

studylevel <- studylevel |>
  mutate(Comparison = recode(Comparison,
    "Agriculture/Natural vegetation" = "Agriculture",
    "Forestry/Natural vegetation" = "Forestry",
    "Urban/Natural vegetation" = "Urban"
  ) |>
    fct_relevel("Agriculture", "Forestry", "Urban"))

lat_mean <- conditional_effects(lat_gamma_S_model4,
  effects = "Latitude"
)$Latitude

lat_mean_lu <- conditional_effects(lat_gamma_S_model5,
  effects = "Latitude:Comparison"
)$"Latitude:Comparison"

lat_mean_lu <- lat_mean_lu |>
  mutate(Comparison = recode(Comparison,
    "Agriculture/Natural vegetation" = "Agriculture",
    "Forestry/Natural vegetation" = "Forestry",
    "Urban/Natural vegetation" = "Urban"
  ) |>
    fct_relevel("Agriculture", "Forestry", "Urban"))

gamma_latitude_S <- studylevel |>
  mutate(Latitude = abs(Latitude), Dataset_id = as.character(Dataset_id)) |>
  .plot_comparison_latitude(lab_y = "γ")

# Taxa
studylevel <- studylevel |>
  mutate(
    Taxa = factor(
      x = Taxa,
      levels = c("Algae", "Fish", "Macroinvertebrate"),
      labels = c("Algae", "Fish", "Macroinvertebrates")
    )
  )

# set.seed(111)
# taxa_gamma_S1 <- brm(mean ~ Taxa + (1 | Dataset_id),
#   data = studylevel,
#   family = "gaussian",
#   iter = 2000,
#   warmup = 1000,
#   cores = 4,
#   chains = 4
# )
#
# save(taxa_gamma_S1, file = "taxa_gamma_S1.Rdata")
#
# taxa_gamma_S2 <- brm(mean ~ Taxa * Comparison + (1 | Dataset_id),
#   data = studylevel,
#   family = "gaussian",
#   iter = 2000,
#   warmup = 1000,
#   cores = 4,
#   chains = 4
# )
#
# save(taxa_gamma_S2, file = "taxa_gamma_S2.Rdata")

taxa_mean <- conditional_effects(taxa_gamma_S1, effects = "Taxa")$Taxa

taxa_mean_lu <- conditional_effects(taxa_gamma_S2, effects = "Taxa:Comparison")$`Taxa:Comparison`

unique(taxa_mean_lu$Comparison)

taxa_mean_lu <- taxa_mean_lu |>
  mutate(Comparison = recode(Comparison,
    "Agriculture/Natural vegetation" = "Agriculture",
    "Forestry/Natural vegetation" = "Forestry",
    "Urban/Natural vegetation" = "Urban"
  ) |>
    fct_relevel("Agriculture", "Forestry", "Urban"))

gamma_taxa_S <- .plot_comparison_taxa(
  studylevel = studylevel |> mutate(Dataset_id = as.character(Dataset_id)),
  taxamean = taxamean,
  lab_y = "γ"
)

## gamma Spie ####
gamma_Spie_comparison <- .wrangle_comparisons(
  scale = "gamma",
  metric = "Spie",
  study_level = TRUE
)

### gamma Spie plot comparison at study level####
# Latitude
studylevel <- gamma_Spie_comparison$Comparisons |>
  summarise(
    mean = mean(Ratio),
    lower_ci = mean -
      stats::qt(0.975, df = n() - 1) * (stats::sd(Ratio) / sqrt(n())),
    upper_ci = mean +
      stats::qt(0.975, df = n() - 1) * (stats::sd(Ratio) / sqrt(n())),
    .by = c(Dataset_id, Comparison)
  ) |>
  left_join(Metadata, by = join_by(Dataset_id))

studylevel <- studylevel |>
  mutate(Comparison = recode(Comparison,
    "Agriculture/Natural vegetation" = "Agriculture",
    "Forestry/Natural vegetation" = "Forestry",
    "Urban/Natural vegetation" = "Urban"
  ) |>
    fct_relevel("Agriculture", "Forestry", "Urban"))

lat_mean <- conditional_effects(lat_gamma_Spie_model4,
  effects = "Latitude"
)$Latitude

lat_mean_lu <- conditional_effects(lat_gamma_Spie_model5,
  effects = "Latitude:Comparison"
)$"Latitude:Comparison"

lat_mean_lu <- lat_mean_lu |>
  mutate(Comparison = recode(Comparison,
    "Agriculture/Natural vegetation" = "Agriculture",
    "Forestry/Natural vegetation" = "Forestry",
    "Urban/Natural vegetation" = "Urban"
  ) |>
    fct_relevel("Agriculture", "Forestry", "Urban"))

gamma_latitude_Spie <- studylevel |>
  mutate(Latitude = abs(Latitude), Dataset_id = as.character(Dataset_id)) |>
  .plot_comparison_latitude()

# Taxa
studylevel <- studylevel |>
  mutate(
    Taxa = factor(
      Taxa,
      levels = c("Algae", "Fish", "Macroinvertebrate"),
      labels = c("Algae", "Fish", "Macroinvertebrates")
    )
  )

# set.seed(111)
# taxa_gamma_Spie1 <- brm(mean ~ Taxa + (1 | Dataset_id),
#   data = studylevel,
#   family = "gaussian",
#   iter = 2000,
#   warmup = 1000,
#   cores = 4,
#   chains = 4
# )
#
# save(taxa_gamma_Spie1, file = "taxa_gamma_Spie1.Rdata")
#
# set.seed(111)
# taxa_gamma_Spie2 <- brm(mean ~ Taxa * Comparison + (1 | Dataset_id),
#   data = studylevel,
#   family = "gaussian",
#   iter = 2000,
#   warmup = 1000,
#   cores = 4,
#   chains = 4
# )
#
# save(taxa_gamma_Spie2, file = "taxa_gamma_Spie2.Rdata")

taxa_mean <- conditional_effects(taxa_gamma_Spie1, effects = "Taxa")$Taxa

taxa_mean_lu <- conditional_effects(taxa_gamma_Spie2, effects = "Taxa:Comparison")$`Taxa:Comparison`

unique(taxa_mean_lu$Comparison)

taxa_mean_lu <- taxa_mean_lu |>
  mutate(Comparison = recode(Comparison,
    "Agriculture/Natural vegetation" = "Agriculture",
    "Forestry/Natural vegetation" = "Forestry",
    "Urban/Natural vegetation" = "Urban"
  ) |>
    fct_relevel("Agriculture", "Forestry", "Urban"))

gamma_taxa_Spie <- .plot_comparison_taxa(
  studylevel = studylevel |> mutate(Dataset_id = as.character(Dataset_id)),
  taxamean = taxamean
)


## beta S ####
beta_S_comparison <- .wrangle_comparisons(
  scale = "beta",
  metric = "S",
  study_level = TRUE
)

### beta S plot comparison at study level####
# Latitude
studylevel <- beta_S_comparison$Comparisons |>
  summarise(
    mean = mean(Ratio),
    lower_ci = mean -
      stats::qt(0.975, df = n() - 1) * (stats::sd(Ratio) / sqrt(n())),
    upper_ci = mean +
      stats::qt(0.975, df = n() - 1) * (stats::sd(Ratio) / sqrt(n())),
    .by = c(Dataset_id, Comparison)
  ) |>
  left_join(Metadata, by = join_by(Dataset_id))

studylevel <- studylevel |>
  mutate(Comparison = recode(Comparison,
    "Agriculture/Natural vegetation" = "Agriculture",
    "Forestry/Natural vegetation" = "Forestry",
    "Urban/Natural vegetation" = "Urban"
  ) |>
    fct_relevel("Agriculture", "Forestry", "Urban"))

lat_mean <- conditional_effects(lat_beta_S_model4,
  effects = "Latitude"
)$Latitude

lat_mean_lu <- conditional_effects(lat_beta_S_model5,
  effects = "Latitude:Comparison"
)$"Latitude:Comparison"

lat_mean_lu <- lat_mean_lu |>
  mutate(Comparison = recode(Comparison,
    "Agriculture" = "Agriculture",
    "Forestry" = "Forestry",
    "Urban" = "Urban"
  ) |>
    fct_relevel("Agriculture", "Forestry", "Urban"))

beta_latitude_S <- studylevel |>
  mutate(Latitude = abs(Latitude), Dataset_id = as.character(Dataset_id)) |>
  .plot_comparison_latitude(lab_y = "β")

# Taxa
studylevel <- studylevel |>
  mutate(
    Taxa = factor(
      Taxa,
      levels = c("Algae", "Fish", "Macroinvertebrate"),
      labels = c("Algae", "Fish", "Macroinvertebrates")
    )
  )

# set.seed(111)
# taxa_beta_S1 <- brm(mean ~ Taxa + (1 | Dataset_id),
#   data = studylevel,
#   family = "gaussian",
#   iter = 2000,
#   warmup = 1000,
#   cores = 4,
#   chains = 4
# )
#
# save(taxa_beta_S1, file = "taxa_beta_S1.Rdata")
#
# set.seed(111)
# taxa_beta_S2 <- brm(mean ~ Taxa * Comparison + (1 | Dataset_id),
#   data = studylevel,
#   family = "gaussian",
#   iter = 2000,
#   warmup = 1000,
#   cores = 4,
#   chains = 4
# )
#
# save(taxa_beta_S2, file = "taxa_beta_S2.Rdata")

taxa_mean <- conditional_effects(taxa_beta_S1, effects = "Taxa")$Taxa

taxa_mean_lu <- conditional_effects(taxa_beta_S2, effects = "Taxa:Comparison")$`Taxa:Comparison`

unique(taxa_mean_lu$Comparison)

taxa_mean_lu <- taxa_mean_lu |>
  mutate(Comparison = recode(Comparison,
    "Agriculture/Natural vegetation" = "Agriculture",
    "Forestry/Natural vegetation" = "Forestry",
    "Urban/Natural vegetation" = "Urban"
  ) |>
    fct_relevel("Agriculture", "Forestry", "Urban"))

beta_taxa_S <- .plot_comparison_taxa(
  studylevel = studylevel |> mutate(Dataset_id = as.character(Dataset_id)),
  taxamean = taxamean,
  lab_y = "β"
)


## beta Spie #####
beta_Spie_comparison <- .wrangle_comparisons(
  scale = "beta",
  metric = "Spie",
  study_level = TRUE
)

### beta Spie plot comparison at study level####
# Latitude
studylevel <- beta_Spie_comparison$Comparisons |>
  summarise(
    mean = mean(Ratio),
    lower_ci = mean -
      stats::qt(0.975, df = n() - 1) * (stats::sd(Ratio) / sqrt(n())),
    upper_ci = mean +
      stats::qt(0.975, df = n() - 1) * (stats::sd(Ratio) / sqrt(n())),
    .by = c(Dataset_id, Comparison)
  ) |>
  left_join(Metadata, by = join_by(Dataset_id))

studylevel <- studylevel |>
  mutate(Comparison = recode(Comparison,
    "Agriculture/Natural vegetation" = "Agriculture",
    "Forestry/Natural vegetation" = "Forestry",
    "Urban/Natural vegetation" = "Urban"
  ) |>
    fct_relevel("Agriculture", "Forestry", "Urban"))

lat_mean <- conditional_effects(lat_beta_Spie_model4,
  effects = "Latitude"
)$Latitude

lat_mean_lu <- conditional_effects(lat_beta_Spie_model5,
  effects = "Latitude:Comparison"
)$"Latitude:Comparison"

lat_mean_lu <- lat_mean_lu |>
  mutate(Comparison = recode(Comparison,
    "Agriculture/Natural vegetation" = "Agriculture",
    "Forestry/Natural vegetation" = "Forestry",
    "Urban/Natural vegetation" = "Urban"
  ) |>
    fct_relevel("Agriculture", "Forestry", "Urban"))

beta_latitude_Spie <- studylevel |>
  mutate(Latitude = abs(Latitude), Dataset_id = as.character(Dataset_id)) |>
  .plot_comparison_latitude()

# Taxa
studylevel <- studylevel |>
  mutate(
    Taxa = factor(
      Taxa,
      levels = c("Algae", "Fish", "Macroinvertebrate"),
      labels = c("Algae", "Fish", "Macroinvertebrates")
    )
  )

# set.seed(111)
# taxa_beta_Spie1 <- brm(mean ~ Taxa + (1 | Dataset_id),
#   data = studylevel,
#   family = "gaussian",
#   iter = 2000,
#   warmup = 1000,
#   cores = 4,
#   chains = 4
# )
#
# save(taxa_beta_Spie1, file = "taxa_beta_Spie1.Rdata")
#
# set.seed(111)
# taxa_beta_Spie2 <- brm(mean ~ Taxa * Comparison + (1 | Dataset_id),
#   data = studylevel,
#   family = "gaussian",
#   iter = 2000,
#   warmup = 1000,
#   cores = 4,
#   chains = 4
# )
#
# save(taxa_beta_Spie2, file = "taxa_beta_Spie2.Rdata")

taxa_mean <- conditional_effects(taxa_beta_Spie1, effects = "Taxa")$Taxa

taxa_mean_lu <- conditional_effects(taxa_beta_Spie2, effects = "Taxa:Comparison")$`Taxa:Comparison`

unique(taxa_mean_lu$Comparison)

taxa_mean_lu <- taxa_mean_lu |>
  mutate(Comparison = recode(Comparison,
    "Agriculture/Natural vegetation" = "Agriculture",
    "Forestry/Natural vegetation" = "Forestry",
    "Urban/Natural vegetation" = "Urban"
  ) |>
    fct_relevel("Agriculture", "Forestry", "Urban"))

beta_taxa_Spie <- .plot_comparison_taxa(
  studylevel = studylevel |> mutate(Dataset_id = as.character(Dataset_id)),
  taxamean = taxamean
)

## combine plots #########
library(cowplot)
Fig2 <- cowplot::plot_grid(
  alpha_latitude_S,
  alpha_latitude_Spie,
  gamma_latitude_S,
  gamma_latitude_Spie,
  beta_latitude_S,
  beta_latitude_Spie,
  ncol = 2,
  align = "v",
  rel_heights = c(1, 1, 1),
  labels = c("A", "B", "C", "D", "E", "F"),
  label_size = 10,
  label_fontface = "plain"
) +
  cowplot::draw_label(
    label = "Absolute latitude",
    y = 0.001,
    hjust = 0.3,
    vjust = 1,
    size = 16
  ) +
  cowplot::draw_label(
    "Effect size",
    x = 0.001,
    y = 0.5,
    angle = 90,
    hjust = 0.4,
    vjust = 0,
    size = 16
  ) +
  theme(plot.margin = margin(t = 0, r = 0, b = 15, l = 19, unit = "pt"))

legend <- get_legend(
  beta_latitude_Spie +
    guides(color = guide_legend(nrow = 1)) +
    theme(
      legend.position = "top",
      legend.text = element_text(
        margin = margin(l = 8, r = 8),
        size = 14
      ),
      legend.key.width = unit(0.65, "cm"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    labs(color = " ")
)

Fig2 <- plot_grid(legend, Fig2,
  ncol = 1, rel_heights = c(0.08, 1)
)

Fig2 <- ggdraw(Fig2) +
  theme(plot.background = element_rect(fill = "white", color = NA))

cowplot::save_plot(
  "figures_revision/Fig2.png",
  Fig2,
  base_height = 11,
  base_width = 8.8,
  units = "in",
  dpi = 600
)

cowplot::save_plot(
  "figures_revision/Fig2.pdf",
  Fig2,
  base_height = 11,
  base_width = 8.8,
  units = "in",
  dpi = 600,
  device = cairo_pdf
)

Fig3 <- cowplot::plot_grid(
  alpha_taxa_S,
  alpha_taxa_Spie,
  gamma_taxa_S,
  gamma_taxa_Spie,
  beta_taxa_S,
  beta_taxa_Spie,
  ncol = 2,
  align = "v",
  rel_heights = c(1, 1, 1),
  labels = c("A", "B", "C", "D", "E", "F"),
  label_size = 10,
  label_fontface = "plain"
) +
  cowplot::draw_label(
    label = "Effect size",
    y = 0.001,
    hjust = 0.3,
    vjust = 1,
    size = 16
  ) +
  theme(plot.margin = margin(t = 0, r = 0, b = 15, l = 0, unit = "pt"))

legend <- get_legend(
  beta_taxa_Spie +
    guides(color = guide_legend(nrow = 1)) +
    theme(
      legend.position = "top",
      legend.text = element_text(
        margin = margin(l = 8, r = 8),
        size = 14
      ),
      legend.key.width = unit(0.65, "cm"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    labs(color = " ")
)

Fig3 <- plot_grid(legend, Fig3,
  ncol = 1, rel_heights = c(0.08, 1)
)

Fig3 <- ggdraw(Fig3) +
  theme(plot.background = element_rect(fill = "white", color = NA))

cowplot::save_plot(
  "figures_revision/Fig3.png",
  Fig3,
  base_height = 8,
  base_width = 9,
  units = "in",
  dpi = 600
)

cowplot::save_plot(
  "figures_revision/Fig3.pdf",
  Fig3,
  base_height = 8,
  base_width = 9,
  units = "in",
  dpi = 600,
  device = cairo_pdf
)


## sum model results for lat #####
# alpha_S
lat_alpha_S_sum1 <- as.data.frame(fixef(lat_alpha_S_model4, probs = c(
  0.10, 0.90,
  0.05, 0.95,
  0.025, 0.975
)))

lat_alpha_S_sum2 <- as.data.frame(fixef(lat_alpha_S_model5, probs = c(
  0.10, 0.90,
  0.05, 0.95,
  0.025, 0.975
)))


lat_alpha_S_sum3 <- as.data.frame(summary(lat_alpha_S_model4)$fixed)

lat_alpha_S_sum4 <- as.data.frame(summary(lat_alpha_S_model5)$fixed)

# alpha_Spie
lat_alpha_Spie_sum1 <- as.data.frame(fixef(lat_alpha_Spie_model4, probs = c(
  0.10, 0.90,
  0.05, 0.95,
  0.025, 0.975
)))

lat_alpha_Spie_sum2 <- as.data.frame(fixef(lat_alpha_Spie_model5, probs = c(
  0.10, 0.90,
  0.05, 0.95,
  0.025, 0.975
)))


lat_alpha_Spie_sum3 <- as.data.frame(summary(lat_alpha_Spie_model4)$fixed)

lat_alpha_Spie_sum4 <- as.data.frame(summary(lat_alpha_Spie_model5)$fixed)

# gamma_S
lat_gamma_S_sum1 <- as.data.frame(fixef(lat_gamma_S_model4, probs = c(
  0.10, 0.90,
  0.05, 0.95,
  0.025, 0.975
)))

lat_gamma_S_sum2 <- as.data.frame(fixef(lat_gamma_S_model5, probs = c(
  0.10, 0.90,
  0.05, 0.95,
  0.025, 0.975
)))


lat_gamma_S_sum3 <- as.data.frame(summary(lat_gamma_S_model4)$fixed)

lat_gamma_S_sum4 <- as.data.frame(summary(lat_gamma_S_model5)$fixed)

# gamma_Spie
lat_gamma_Spie_sum1 <- as.data.frame(fixef(lat_gamma_Spie_model4, probs = c(
  0.10, 0.90,
  0.05, 0.95,
  0.025, 0.975
)))

lat_gamma_Spie_sum2 <- as.data.frame(fixef(lat_gamma_Spie_model5, probs = c(
  0.10, 0.90,
  0.05, 0.95,
  0.025, 0.975
)))

lat_gamma_Spie_sum3 <- as.data.frame(summary(lat_gamma_Spie_model4)$fixed)

lat_gamma_Spie_sum4 <- as.data.frame(summary(lat_gamma_Spie_model5)$fixed)

# beta_S
lat_beta_S_sum1 <- as.data.frame(fixef(lat_beta_S_model4, probs = c(
  0.10, 0.90,
  0.05, 0.95,
  0.025, 0.975
)))

lat_beta_S_sum2 <- as.data.frame(fixef(lat_beta_S_model5, probs = c(
  0.10, 0.90,
  0.05, 0.95,
  0.025, 0.975
)))


lat_beta_S_sum3 <- as.data.frame(summary(lat_beta_S_model4)$fixed)

lat_beta_S_sum4 <- as.data.frame(summary(lat_beta_S_model5)$fixed)

# beta_Spie
lat_beta_Spie_sum1 <- as.data.frame(fixef(lat_beta_Spie_model4, probs = c(
  0.10, 0.90,
  0.05, 0.95,
  0.025, 0.975
)))

lat_beta_Spie_sum2 <- as.data.frame(fixef(lat_beta_Spie_model5, probs = c(
  0.10, 0.90,
  0.05, 0.95,
  0.025, 0.975
)))


lat_beta_Spie_sum3 <- as.data.frame(summary(lat_beta_Spie_model4)$fixed)

lat_beta_Spie_sum4 <- as.data.frame(summary(lat_beta_Spie_model5)$fixed)


sum_lat_model12 <- bind_rows(
  lat_alpha_S_sum1,
  lat_alpha_S_sum2,
  lat_alpha_Spie_sum1,
  lat_alpha_Spie_sum2,
  lat_gamma_S_sum1,
  lat_gamma_S_sum2,
  lat_gamma_Spie_sum1,
  lat_gamma_Spie_sum2,
  lat_beta_S_sum1,
  lat_beta_S_sum2,
  lat_beta_Spie_sum1,
  lat_beta_Spie_sum2,
)

sum_lat_model34 <- bind_rows(
  lat_alpha_S_sum3,
  lat_alpha_S_sum4,
  lat_alpha_Spie_sum3,
  lat_alpha_Spie_sum4,
  lat_gamma_S_sum3,
  lat_gamma_S_sum4,
  lat_gamma_Spie_sum3,
  lat_gamma_Spie_sum4,
  lat_beta_S_sum3,
  lat_beta_S_sum4,
  lat_beta_Spie_sum3,
  lat_beta_Spie_sum4
)

## sum model results for taxa #####
# alpha_S
taxa_alpha_S_sum1 <- as.data.frame(fixef(taxa_alpha_S1, probs = c(
  0.10, 0.90,
  0.05, 0.95,
  0.025, 0.975
)))

taxa_alpha_S_sum2 <- as.data.frame(fixef(taxa_alpha_S2, probs = c(
  0.10, 0.90,
  0.05, 0.95,
  0.025, 0.975
)))


taxa_alpha_S_sum3 <- as.data.frame(summary(taxa_alpha_S1)$fixed)

taxa_alpha_S_sum4 <- as.data.frame(summary(taxa_alpha_S2)$fixed)

# alpha_Spie
taxa_alpha_Spie_sum1 <- as.data.frame(fixef(taxa_alpha_Spie1, probs = c(
  0.10, 0.90,
  0.05, 0.95,
  0.025, 0.975
)))

taxa_alpha_Spie_sum2 <- as.data.frame(fixef(taxa_alpha_Spie2, probs = c(
  0.10, 0.90,
  0.05, 0.95,
  0.025, 0.975
)))


taxa_alpha_Spie_sum3 <- as.data.frame(summary(taxa_alpha_Spie1)$fixed)

taxa_alpha_Spie_sum4 <- as.data.frame(summary(taxa_alpha_Spie2)$fixed)

# gamma_S
taxa_gamma_S_sum1 <- as.data.frame(fixef(taxa_gamma_S1, probs = c(
  0.10, 0.90,
  0.05, 0.95,
  0.025, 0.975
)))

taxa_gamma_S_sum2 <- as.data.frame(fixef(taxa_gamma_S2, probs = c(
  0.10, 0.90,
  0.05, 0.95,
  0.025, 0.975
)))


taxa_gamma_S_sum3 <- as.data.frame(summary(taxa_gamma_S1)$fixed)

taxa_gamma_S_sum4 <- as.data.frame(summary(taxa_gamma_S2)$fixed)

# gamma_Spie
taxa_gamma_Spie_sum1 <- as.data.frame(fixef(taxa_gamma_Spie1, probs = c(
  0.10, 0.90,
  0.05, 0.95,
  0.025, 0.975
)))

taxa_gamma_Spie_sum2 <- as.data.frame(fixef(taxa_gamma_Spie2, probs = c(
  0.10, 0.90,
  0.05, 0.95,
  0.025, 0.975
)))

taxa_gamma_Spie_sum3 <- as.data.frame(summary(taxa_gamma_Spie1)$fixed)

taxa_gamma_Spie_sum4 <- as.data.frame(summary(taxa_gamma_Spie2)$fixed)

# beta_S
taxa_beta_S_sum1 <- as.data.frame(fixef(taxa_beta_S1, probs = c(
  0.10, 0.90,
  0.05, 0.95,
  0.025, 0.975
)))

taxa_beta_S_sum2 <- as.data.frame(fixef(taxa_beta_S2, probs = c(
  0.10, 0.90,
  0.05, 0.95,
  0.025, 0.975
)))


taxa_beta_S_sum3 <- as.data.frame(summary(taxa_beta_S1)$fixed)

taxa_beta_S_sum4 <- as.data.frame(summary(taxa_beta_S2)$fixed)

# beta_Spie
taxa_beta_Spie_sum1 <- as.data.frame(fixef(taxa_beta_Spie1, probs = c(
  0.10, 0.90,
  0.05, 0.95,
  0.025, 0.975
)))

taxa_beta_Spie_sum2 <- as.data.frame(fixef(taxa_beta_Spie2, probs = c(
  0.10, 0.90,
  0.05, 0.95,
  0.025, 0.975
)))


taxa_beta_Spie_sum3 <- as.data.frame(summary(taxa_beta_Spie1)$fixed)

taxa_beta_Spie_sum4 <- as.data.frame(summary(taxa_beta_Spie2)$fixed)


sum_taxa_model12 <- bind_rows(
  taxa_alpha_S_sum1,
  taxa_alpha_S_sum2,
  taxa_alpha_Spie_sum1,
  taxa_alpha_Spie_sum2,
  taxa_gamma_S_sum1,
  taxa_gamma_S_sum2,
  taxa_gamma_Spie_sum1,
  taxa_gamma_Spie_sum2,
  taxa_beta_S_sum1,
  taxa_beta_S_sum2,
  taxa_beta_Spie_sum1,
  taxa_beta_Spie_sum2,
)

sum_taxa_model34 <- bind_rows(
  taxa_alpha_S_sum3,
  taxa_alpha_S_sum4,
  taxa_alpha_Spie_sum3,
  taxa_alpha_Spie_sum4,
  taxa_gamma_S_sum3,
  taxa_gamma_S_sum4,
  taxa_gamma_Spie_sum3,
  taxa_gamma_Spie_sum4,
  taxa_beta_S_sum3,
  taxa_beta_S_sum4,
  taxa_beta_Spie_sum3,
  taxa_beta_Spie_sum4
)

sum_taxa_model34 <- rownames_to_column(sum_taxa_model34, var = "var") |>
  as_tibble()

readr::write_excel_csv(sum_taxa_model34, "sum_taxa_model34.csv")

sum_taxa_model12 <- rownames_to_column(sum_taxa_model12, var = "var") |>
  as_tibble()

readr::write_excel_csv(sum_taxa_model12,
  "sum_taxa_model12.csv",
  col_names = TRUE
)
