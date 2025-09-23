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
  left_join(Metadata, by = 'Dataset_id')

alpha_latitude_S <- studylevel |>
  mutate(Latitude = abs(Latitude), Dataset_id = as.character(Dataset_id)) |>
  .plot_comparison_latitude(lab_y = "α", subtitle = "S")

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
taxa_mean <- studylevel |>
  summarise(
    taxa_mean = mean(mean),
    lower_ci = taxa_mean -
      stats::qt(0.975, df = n() - 1) * (stats::sd(mean) / sqrt(n())),
    upper_ci = taxa_mean +
      stats::qt(0.975, df = n() - 1) * (stats::sd(mean) / sqrt(n())),
    .by = Taxa
  )

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
  left_join(Metadata, by = 'Dataset_id')

alpha_latitude_Spie <- studylevel |>
  mutate(Latitude = abs(Latitude), Dataset_id = as.character(Dataset_id)) |>
  .plot_comparison_latitude(subtitle = "Spie")

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
taxa_mean <- studylevel |>
  summarise(
    taxa_mean = mean(mean),
    lower_ci = taxa_mean -
      stats::qt(0.975, df = n() - 1) * (stats::sd(mean) / sqrt(n())),
    upper_ci = taxa_mean +
      stats::qt(0.975, df = n() - 1) * (stats::sd(mean) / sqrt(n())),
    .by = Taxa
  )

alpha_taxa_Spie <- .plot_comparison_taxa(
  studylevel = studylevel |> mutate(Dataset_id = as.character(Dataset_id)),
  taxamean = taxamean,
  subtitle = "Spie"
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

taxa_mean <- studylevel |>
  summarise(
    taxa_mean = mean(mean),
    lower_ci = taxa_mean -
      stats::qt(0.975, df = n() - 1) * (stats::sd(mean) / sqrt(n())),
    upper_ci = taxa_mean +
      stats::qt(0.975, df = n() - 1) * (stats::sd(mean) / sqrt(n())),
    .by = Taxa
  )

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
taxa_mean <- studylevel |>
  summarise(
    taxa_mean = mean(mean),
    lower_ci = taxa_mean -
      stats::qt(0.975, df = n() - 1) * (stats::sd(mean) / sqrt(n())),
    upper_ci = taxa_mean +
      stats::qt(0.975, df = n() - 1) * (stats::sd(mean) / sqrt(n())),
    .by = Taxa
  )

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
taxa_mean <- studylevel |>
  summarise(
    taxa_mean = mean(mean),
    lower_ci = taxa_mean -
      stats::qt(0.975, df = n() - 1) * (stats::sd(mean) / sqrt(n())),
    upper_ci = taxa_mean +
      stats::qt(0.975, df = n() - 1) * (stats::sd(mean) / sqrt(n())),
    .by = Taxa
  )

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
taxa_mean <- studylevel |>
  summarise(
    taxa_mean = mean(mean),
    lower_ci = taxa_mean -
      stats::qt(0.975, df = n() - 1) * (stats::sd(mean) / sqrt(n())),
    upper_ci = taxa_mean +
      stats::qt(0.975, df = n() - 1) * (stats::sd(mean) / sqrt(n())),
    .by = Taxa
  )

beta_taxa_Spie <- .plot_comparison_taxa(
  studylevel = studylevel |> mutate(Dataset_id = as.character(Dataset_id)),
  taxamean = taxamean
)

## combine plots #########
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
  labels = c('A', 'B', 'C', 'D', 'E', 'F'),
  label_size = 10,
  label_fontface = "plain"
) +
  cowplot::draw_label(
    label = 'Absolute latitude',
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

cowplot::save_plot(
  "figures/Fig2.png",
  Fig2,
  base_height = 11,
  base_width = 8.8,
  units = "in",
  dpi = 600
)

cowplot::save_plot(
  "figures/Fig2.pdf",
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
  labels = c('A', 'B', 'C', 'D', 'E', 'F'),
  label_size = 10,
  label_fontface = "plain"
) +
  cowplot::draw_label(
    label = 'Effect size',
    y = 0.001,
    hjust = 0.3,
    vjust = 1,
    size = 16
  ) +
  theme(plot.margin = margin(t = 0, r = 0, b = 15, l = 0, unit = "pt"))

cowplot::save_plot(
  "figures/Fig3.png",
  Fig3,
  base_height = 7,
  base_width = 8.8,
  units = "in",
  dpi = 600
)

cowplot::save_plot(
  "figures/Fig3.pdf",
  Fig3,
  base_height = 7,
  base_width = 8.8,
  units = "in",
  dpi = 600,
  device = cairo_pdf
)
