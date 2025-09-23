library(tidyverse)
library(tidybayes)

# column for stream order #ORD_STRA
#for alpha S: S_alpha_alt_3
#for alpha Spie: Spie_alpha_alt_2

mob_data_alpha_so <- read.table(
  file = "results/temp/mob_data_alpha_so.txt",
  quote = "\"",
  comment.char = ""
) |>
  select(Dataset_id, ORD_STRA) |>
  distinct()

Spie_alpha_alt_2 <- readRDS("results/models/Spie_alpha_alt_2.Rds")
S_alpha_alt_3 <- readRDS("results/models/S_alpha_alt_3.Rds")

# alpha S #####
#draw study level estimates
Forest <- S_alpha_alt_3 |>
  spread_draws(
    b_Land_useForest,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  ) |>
  mutate(Forest = b_Land_useForest + r_Dataset_id) |>
  filter(Land_use == "Land_useForest")

Agriculture <- S_alpha_alt_3 |>
  spread_draws(
    b_Land_useAgriculture,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  ) |>
  mutate(Agriculture = b_Land_useAgriculture + r_Dataset_id) |>
  filter(Land_use == "Land_useAgriculture")

Urban <- S_alpha_alt_3 |>
  spread_draws(
    b_Land_useUrban,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  ) |>
  mutate(Urban = b_Land_useUrban + r_Dataset_id) |>
  filter(Land_use == "Land_useUrban")

Forestry <- S_alpha_alt_3 |>
  spread_draws(
    b_Land_useForestry,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  ) |>
  mutate(Forestry = b_Land_useForestry + r_Dataset_id) |>
  filter(Land_use == "Land_useForestry")

Forest$AN <- Agriculture$Agriculture - Forest$Forest
Forest$UN <- Urban$Urban - Forest$Forest
Forest$FN <- Forestry$Forestry - Forest$Forest

ComparisonsA0 <- Forest |>
  pivot_longer(
    cols = c("AN", "UN", "FN"),
    names_to = "Comparison",
    values_to = "Ratio"
  ) |>
  left_join(
    y = mob_data_alpha_so,
    by = join_by(Dataset_id),
    relationship = "many-to-many"
  ) |>
  mutate(
    Dataset_id = as.factor(Dataset_id),
    Land_use = as.factor(Land_use),
    ORD_STRA = as.factor(ORD_STRA),
    Comparison = factor(
      Comparison,
      levels = c("AN", "UN", "FN"),
      labels = c(
        'Agriculture/Forest',
        'Urban/Forest',
        'Forestry/Forest'
      )
    )
  )

#stream order
alpha_stream_order_S <- ggplot() +
  ggridges::geom_density_ridges_gradient(
    data = ComparisonsA0,
    aes(
      x = Ratio,
      y = ORD_STRA
    ),
    scale = 0.9,
    alpha = 0.5,
    linetype = 0
  ) +
  xlim(-1, 1) +
  geom_vline(
    data = ComparisonsA0,
    aes(xintercept = mean(Ratio)),
    linewidth = 0.5,
    alpha = 0.5,
    lty = 2
  ) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  geom_point(
    data = ComparisonsA0,
    aes(x = Ratio, y = ORD_STRA),
    stat = ggstance:::StatSummaryh,
    fun.x = median,
    size = 1.25,
    shape = 18,
    colour = 'black'
  ) +
  labs(y = 'S', x = 'Log Ratio') +
  geom_text(
    data = ComparisonsA0 |>
      filter(Ratio < 0) |>
      group_by(ORD_STRA) |>
      summarise(Count = n()) |>
      mutate(percentages = Count / 1000) |>
      ungroup() |>
      distinct(ORD_STRA, percentages, .keep_all = T),
    aes(
      x = -0.85,
      y = ORD_STRA,
      label = paste(percentages)
    ),
    size = 3,
    nudge_y = 0.5,
    parse = T
  ) +
  theme(legend.position = "none") +
  scale_fill_manual(
    values = c("black", "gray"),
    labels = c("Below 0", "Above 0")
  ) +
  scale_y_discrete(labels = scales::wrap_format(9)) +
  theme(
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 7)
  )

# alpha Spie#####

Forest <- Spie_alpha_alt_2 |>
  spread_draws(
    b_Land_useForest,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  ) |>
  mutate(Forest = b_Land_useForest + r_Dataset_id) |>
  filter(Land_use == "Land_useForest")

Agriculture <- Spie_alpha_alt_2 |>
  spread_draws(
    b_Land_useAgriculture,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  ) |>
  mutate(Agriculture = b_Land_useAgriculture + r_Dataset_id) |>
  filter(Land_use == "Land_useAgriculture")

Urban <- Spie_alpha_alt_2 |>
  spread_draws(
    b_Land_useUrban,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  ) |>
  mutate(Urban = b_Land_useUrban + r_Dataset_id) |>
  filter(Land_use == "Land_useUrban")

Forestry <- Spie_alpha_alt_2 |>
  spread_draws(
    b_Land_useForestry,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  ) |>
  mutate(Forestry = b_Land_useForestry + r_Dataset_id) |>
  filter(Land_use == "Land_useForestry")

Forest$AN <- Agriculture$Agriculture -
  Forest$Forest
Forest$UN <- Urban$Urban -
  Forest$Forest
Forest$FN <- Forestry$Forestry -
  Forest$Forest

ComparisonsA2 <- Forest |>
  pivot_longer(
    cols = c("AN", "UN", "FN"),
    names_to = "Comparison",
    values_to = "Ratio"
  ) |>
  left_join(
    y = mob_data_alpha_so,
    by = join_by(Dataset_id),
    relationship = "many-to-many"
  ) |>
  mutate(
    Dataset_id = as.factor(Dataset_id),
    Land_use = as.factor(Land_use),
    ORD_STRA = as.factor(ORD_STRA),
    Comparison = factor(
      Comparison,
      levels = c("AN", "UN", "FN"),
      labels = c('Agriculture/Forest', 'Urban/Forest', 'Forestry/Forest')
    )
  )

### plot S and Spie in different stream order
Comparisons02 <- bind_cols(
  ComparisonsA0 |>
    rename(RatioA0 = Ratio) |>
    select(Dataset_id, Land_use, ORD_STRA, Comparison, RatioA0),
  ComparisonsA2 |>
    rename(RatioA2 = Ratio) |>
    select(Dataset_id, Land_use, ORD_STRA, Comparison, RatioA2)
) |>
  select(ORD_STRA = ORD_STRA...3, RatioA0, RatioA2) |>
  pivot_longer(
    cols = starts_with("ratio"),
    names_to = "ratio_type",
    values_to = "Ratio"
  )

n_total <- Comparisons02 |>
  summarise(n_total = n(), .by = c(ORD_STRA, ratio_type))

n_count <- Comparisons02 |>
  filter(Ratio < 0) |>
  summarise(n_count = n(), .by = c(ORD_STRA, ratio_type))

n_total_count <- left_join(
  x = n_total,
  y = n_count,
  by = join_by(ORD_STRA, ratio_type)
)

#plot
alpha_stream_order <- ggplot() +
  ggridges::geom_density_ridges_gradient(
    data = Comparisons02,
    aes(
      x = Ratio,
      y = ORD_STRA,
      fill = ratio_type
    ),
    scale = 0.9,
    alpha = 0.1,
    linetype = 0
  ) +
  geom_point(
    data = Comparisons02,
    aes(x = Ratio, y = ORD_STRA, color = ratio_type),
    stat = ggstance:::StatSummaryh,
    fun.x = median,
    shape = 21,
    fill = "white",
    stroke = 2,
    size = 4.5
  ) +
  xlim(-1.0, 1.0) +
  geom_text(
    data = n_total_count |>
      mutate(d = n_count / n_total) |>
      mutate(percentages = sprintf('%.2f', d)) |>
      distinct(ORD_STRA, ratio_type, percentages, .keep_all = T),
    aes(
      x = -0.85,
      y = as.numeric(ORD_STRA) + ifelse(ratio_type == "RatioA0", 0.55, 0.25),
      color = ratio_type,
      label = percentages
    ),
    size = 4.5,
    parse = T,
    show.legend = FALSE
  ) +
  geom_vline(xintercept = 0, linewidth = 0.5, linetype = 2) +
  labs(y = 'Basin size', x = 'Effect size') +
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

ggsave("figures/alpha_stream_order.pdf", alpha_stream_order)
