library(tidyverse)
source("functions/pivot_filter.R")
source("functions/wrangling.R")
source("functions/plotting.R")

# Plot ########
## S1 ####
#### alpha S ####
alpha_S_response <- .wrangle_responses(scale = "alpha", metric = "S") |>
  .plot_responses()

#### gamma S ####
gamma_S_response <- .wrangle_responses(scale = "gamma", metric = "S") |>
  .plot_responses()

### combine alpha gamma S plots####
plot_S1 <- cowplot::plot_grid(
  alpha_S_response,
  gamma_S_response,
  rel_widths = c(1, 0.7),
  labels = c('A', 'B')
) +
  cowplot::draw_label(y = 0.015, x = 0.60, label = 'Log Ratio')

ggsave(
  filename = "figures/S1_alpha_gamma_S_responses.png",
  width = 300,
  height = 200,
  units = 'mm'
)

## S2 ####
#### alpha N ####
alpha_N_response <- .wrangle_responses(scale = "alpha", metric = "N") |>
  .plot_responses()

#### gamma N ####
gamma_N_response <- .wrangle_responses(scale = "gamma", metric = "N") |>
  .plot_responses()

### combine alpha gamma N plots####
cowplot::plot_grid(
  alpha_N_response,
  gamma_N_response,
  rel_widths = c(1, 0.7),
  labels = c('A', 'B')
) +
  cowplot::draw_label(y = 0.015, x = 0.60, label = 'Log Ratio')

ggsave(
  filename = "figures/S2_alpha_gamma_N_response.png",
  width = 300,
  height = 200,
  units = 'mm'
)


# S3 ####
library(tidybayes)
library(ggridges)

## alpha N #####
alpha_N_comparison <- .wrangle_comparisons(scale = "alpha", metric = "N") |>
  .plot_comparisons()

## gamma N ####
gamma_N_comparison <- .wrangle_comparisons(scale = "gamma", metric = "N") |>
  .plot_comparisons()

## combine alpha gamma N comparisons ####
cowplot::plot_grid(
  alpha_N_comparison,
  gamma_N_comparison,
  rel_widths = c(1, 0.7),
  labels = c('A', 'B')
) +
  cowplot::draw_label(y = 0.015, x = 0.60, label = 'Log Ratio')
ggsave(
  filename = "figures/S3_alpha_gamma_N_comparisons.png",
  width = 300,
  height = 200,
  units = 'mm'
)

# Fig1 ####
## alpha S, Spie ####
alpha_S_comparison <- .wrangle_comparisons(scale = "alpha", metric = "S")
alpha_Spie_comparison <- .wrangle_comparisons(scale = "alpha", metric = "Spie")

df_alpha <- .wrangle_posterior(
  df0 = alpha_S_comparison$Comparisons,
  df2 = alpha_Spie_comparison$Comparisons,
  suffix = "A"
)
alpha_02_posterior <- .plot_posterior(
  model = df_alpha,
  suffix = "A",
  lab_y = "α",
  show_legend = TRUE
)
## beta S, Spie ####
beta_S_comparison <- .wrangle_comparisons(scale = "beta", metric = "S")
beta_Spie_comparison <- .wrangle_comparisons(scale = "beta", metric = "Spie")

df_beta <- .wrangle_posterior(
  df0 = beta_S_comparison$Comparisons,
  df2 = beta_Spie_comparison$Comparisons,
  suffix = "B"
)

beta_02_posterior <- .plot_posterior(
  model = df_beta,
  suffix = "B",
  lab_x = "Effect size",
  lab_y = "β"
)
## gamma S, Spie ####
gamma_S_comparison <- .wrangle_comparisons(scale = "gamma", metric = "S")
gamma_Spie_comparison <- .wrangle_comparisons(scale = "gamma", metric = "Spie")

df_gamma <- .wrangle_posterior(
  gamma_S_comparison$Comparisons,
  gamma_Spie_comparison$Comparisons,
  suffix = "G"
)
gamma_02_posterior <- .plot_posterior(
  model = df_gamma,
  suffix = "G",
  lab_y = "γ"
)

## combine the three panels vertically ####
Fig1 <- cowplot::plot_grid(
  alpha_02_posterior,
  gamma_02_posterior,
  beta_02_posterior,
  ncol = 1,
  align = "v",
  rel_heights = c(1.19, 1, 1),
  labels = c("A", "B", "C"),
  label_size = 12,
  label_fontface = "plain"
)

cowplot::save_plot(
  filename = "figures/Fig1.png",
  plot = Fig1,
  base_height = 11,
  base_width = 8.8,
  units = "in",
  dpi = 600
)

cowplot::save_plot(
  filename = "figures/Fig1.pdf",
  plot = Fig1,
  base_height = 11,
  base_width = 8.8,
  units = "in",
  dpi = 600,
  device = cairo_pdf
)
