ComparisonsB0 <- ComparisonsB0 |>
  rename(RatioB0 = Ratio) |>
  select(5:6, )

ComparisonsG0 <- ComparisonsG0 |>
  rename(RatioG0 = Ratio) |>
  select(5:6, )

alpha0 <- cbind(ComparisonsB0, ComparisonsG0)

alpha0 <- alpha0 |> select(-3, )

alpha0 <- alpha0 |>
  mutate(RatioA0 = RatioG0 - RatioB0)

ComparisonsB2 <- ComparisonsB2 |>
  rename(RatioB2 = Ratio) |>
  select(5:6, )

ComparisonsG2 <- ComparisonsG2 |>
  rename(RatioG2 = Ratio) |>
  select(5:6, )

alpha2 <- cbind(ComparisonsB2, ComparisonsG2)

alpha2 <- alpha2 |> select(-3, )

alpha2 <- alpha2 |>
  mutate(RatioA2 = RatioG2 - RatioB2)

alpha <- cbind(alpha0, alpha2)

alpha <- alpha |>
  select(c(1, 4, 8))

alpha <- alpha |>
  pivot_longer(
    cols = starts_with("Ratio"),
    names_to = "ratio_type",
    values_to = "Ratio"
  )

alpha$Comparison <- factor(
  alpha$Comparison,
  levels = c(
    "Forestry/Natural vegetation",
    "Agriculture/Natural vegetation",
    "Urban/Agriculture",
    "Urban/Natural vegetation"
  )
)

alpha <- alpha |>
  group_by(Comparison, ratio_type) |>
  summarise(
    lower_90 = quantile(Ratio, probs = 0.05),
    upper_90 = quantile(Ratio, probs = 0.95)
  ) |>
  left_join(alpha, by = c("Comparison", "ratio_type"))

write_csv(alpha, "alpha_calculated_with_beta_gamma.csv")

## plot alpha q=0 and q=2 together ####

Comparisons02 <- alpha_calculated_with_beta_gamma

Comparisons02 |>
  group_by(Comparison) |>
  summarise(n = n())

Comparisons02$Comparison <- factor(
  Comparisons02$Comparison,
  levels = c(
    "Forestry/Natural vegetation",
    "Agriculture/Natural vegetation",
    "Urban/Natural vegetation"
  ),
  labels = c('Forestry', 'Agriculture', 'Urban')
)

alpha_02_posterior <- ggplot() +
  geom_density_ridges_gradient(
    data = Comparisons02,
    aes(
      x = Ratio,
      # y = Comparison,
      y = Comparison,
      fill = ratio_type
      # y = fct_reorder(Comparison, Ratio, .fun = mean),
      # fill = factor(after_stat(x) > 0)
    ),
    scale = 0.9,
    alpha = 0.1,
    linetype = 0
  ) +
  geom_point(
    data = Comparisons02,
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
    data = Comparisons02 |>
      filter(Ratio < 0) |>
      group_by(Comparison, ratio_type) |>
      summarise(Count = n()) |>
      mutate(d = Count / 1000) |>
      mutate(percentages = sprintf('%.2f', d)) |>
      ungroup() |>
      distinct(Comparison, ratio_type, percentages, .keep_all = T),
    aes(
      x = -0.85,
      # y = Comparison,
      y = as.numeric(Comparison) + ifelse(ratio_type == "RatioA0", 0.45, 0.15),
      color = ratio_type,
      label = percentages
    ),
    size = 4.5,
    parse = T,
    show.legend = FALSE
  ) +
  geom_vline(xintercept = 0, size = 0.5, linetype = 2) +
  labs(y = '', x = '') +
  # guides(fill=guide_legend(title=" "))+
  scale_fill_manual(
    values = c(
      "RatioA0" = scales::alpha("#F2790F", 0.7),
      "RatioA2" = scales::alpha("#CAAE10", 0.7)
    ),
    name = " ",
    # labels = c("RatioA0" = " ",
    #            "RatioA2" = " ")
    labels = c("RatioA0" = "S", "RatioA2" = "Spie")
  ) +
  scale_color_manual(
    values = c(
      "RatioA0" = scales::alpha("#F2790F", 0.9),
      "RatioA2" = scales::alpha("#CAAE10", 0.9)
    ),
    name = " ",
    # labels = c("RatioA0" = " ",
    #            "RatioA2" = " ")
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
  # theme(legend.position = "none")+
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
