library(tidyverse)

# Gamma-Beta S ----
ComparisonsG0 <- read.table(
  file = "results/temp/ComparisonsG0.txt",
  quote = "\"",
  comment.char = ""
) |>
  as_tibble() |>
  select(Comparison, RatioG0) |>
  mutate(Comparison = as.factor(Comparison))

ComparisonsB0 <- read.table(
  file = "results/temp/ComparisonsB0.txt",
  quote = "\"",
  comment.char = ""
) |>
  as_tibble() |>
  select(Comparison, RatioB0) |>
  mutate(Comparison = as.factor(Comparison))

alpha0 <- bind_cols(ComparisonsB0, ComparisonsG0) |>
  select(-Comparison...3) |>
  mutate(RatioA0 = RatioG0 - RatioB0)

# Gamma-Beta Spie ----
ComparisonsG2 <- read.table(
  file = "results/temp/ComparisonsG2.txt",
  quote = "\"",
  comment.char = ""
) |>
  as_tibble() |>
  select(Comparison, RatioG2) |>
  mutate(Comparison = as.factor(Comparison))

ComparisonsB2 <- read.table(
  file = "results/temp/ComparisonsB2.txt",
  quote = "\"",
  comment.char = ""
) |>
  as_tibble() |>
  select(Comparison, RatioB2) |>
  mutate(Comparison = as.factor(Comparison))

alpha2 <- bind_cols(ComparisonsB2, ComparisonsG2) |>
  select(-Comparison...3) |>
  mutate(RatioA2 = RatioG2 - RatioB2)

alpha <- bind_cols(alpha0, alpha2) |>
  select(Comparison = Comparison...1, RatioG0, RatioA2) |>
  pivot_longer(
    cols = starts_with("Ratio"),
    names_to = "ratio_type",
    names_transform = as.factor,
    values_to = "Ratio"
  ) |>
  group_by(Comparison, ratio_type) |>
  mutate(
    lower_90 = quantile(Ratio, probs = 0.05),
    upper_90 = quantile(Ratio, probs = 0.95),
    Comparison = fct_recode(
      .f = Comparison,
      "Forestry/Natural vegetation" = "FN",
      "Agriculture/Natural vegetation" = "AN",
      "Urban/Natural vegetation" = "UN"
    )
  )

write_csv(x = alpha, file = "alpha_calculated_with_beta_gamma.csv")

## plot alpha q=0 and q=2 together ####

alpha <- alpha |>
  mutate(
    Comparison = fct_recode(
      .f = Comparison,
      Forestry = "Forestry/Natural vegetation",
      Agriculture = "Agriculture/Natural vegetation",
      Urban = "Urban/Natural vegetation"
    )
  )

alpha_02_posterior <- ggplot() +
  ggridges::geom_density_ridges_gradient(
    data = alpha,
    aes(
      x = Ratio,
      y = Comparison,
      fill = ratio_type
    ),
    scale = 0.9,
    alpha = 0.1,
    linetype = 0
  ) +
  geom_point(
    data = alpha,
    aes(x = Ratio, y = Comparison, color = ratio_type),
    stat = ggstance:::StatSummaryh,
    fun.x = median,
    shape = 21,
    fill = "white",
    stroke = 2,
    size = 4.5
  ) +
  xlim(-1.0, 0.5) +
  geom_text(
    data = alpha |>
      filter(Ratio < 0) |>
      group_by(Comparison, ratio_type) |>
      summarise(Count = n()) |>
      mutate(d = Count / 1000) |>
      mutate(percentages = sprintf('%.2f', d)) |>
      ungroup() |>
      distinct(Comparison, ratio_type, percentages, .keep_all = T),
    aes(
      x = -0.85,
      y = as.numeric(Comparison) + if_else(ratio_type == "RatioA0", 0.45, 0.15),
      color = ratio_type,
      label = percentages
    ),
    size = 4.5,
    parse = T,
    show.legend = FALSE
  ) +
  geom_vline(xintercept = 0, linewidth = 0.5, linetype = 2) +
  labs(y = '', x = '') +
  scale_fill_manual(
    values = c(
      "RatioA0" = scales::alpha("#F2790F", 0.7),
      "RatioA2" = scales::alpha("#CAAE10", 0.7)
    ),
    name = " ",
    labels = c("RatioA0" = "S", "RatioA2" = "Spie")
  ) +
  scale_color_manual(
    values = c(
      "RatioA0" = scales::alpha("#F2790F", 0.9),
      "RatioA2" = scales::alpha("#CAAE10", 0.9)
    ),
    name = " ",
    labels = c("RatioA0" = "S", "RatioA2" = "Spie")
  ) +
  scale_y_discrete(labels = scales::wrap_format(9)) +
  theme_minimal() +
  guides(color = guide_legend(byrow = FALSE)) +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.key.spacing.x = unit(1, "cm")
  ) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.ticks.x = element_line(size = 0.5, color = "gray"),
    panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 13),
    legend.text = element_text(size = 16)
  )
ggsave(filename = "figures/figS5.png")
