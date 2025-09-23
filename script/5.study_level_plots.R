library(tidyverse)
library(brms)
library(tidybayes)
library(ggridges)

Metadata <- read.csv(file = "data/Metadata.csv")

# random effect ####
## alpha S#####
alpha_S_model <- alpha_S_model2
rm(alpha_S_model2)
resid <- residuals(alpha_S_model, type = 'pearson', method = 'predict') |>
  dplyr::as_tibble() |>
  bind_cols(alpha_S_model$data)

fitted <- fitted(alpha_S_model, re_formula = NA)
predict <- predict(alpha_S_model)

resid$fitted <- fitted[, 'Estimate']
resid$predict <- predict[, 'Estimate']

plot(resid$Estimate ~ resid$fitted, ylab = 'Pearson residual')

plot(
  resid$Estimate ~ resid$Land_use,
  ylab = 'Pearson residual',
  xlab = "Land_use"
)

Natural_vegetation <- alpha_S_model |>
  spread_draws(
    b_Land_useNatural_vegetation,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  )

Natural_vegetation <- Natural_vegetation |>
  mutate(Natural_vegetation = b_Land_useNatural_vegetation + r_Dataset_id) |>
  filter(Land_use == "Land_useNatural_vegetation")

Agriculture <- alpha_S_model |>
  spread_draws(
    b_Land_useAgriculture,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  )

Agriculture <- Agriculture |>
  mutate(Agriculture = b_Land_useAgriculture + r_Dataset_id) |>
  filter(Land_use == "Land_useAgriculture")

Urban <- alpha_S_model |>
  spread_draws(
    b_Land_useUrban,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  )
Urban <- Urban |>
  mutate(Urban = b_Land_useUrban + r_Dataset_id) |>
  filter(Land_use == "Land_useUrban")

Forestry <- alpha_S_model |>
  spread_draws(
    b_Land_useForestry,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  )
Forestry <- Forestry |>
  mutate(Forestry = b_Land_useForestry + r_Dataset_id) |>
  filter(Land_use == "Land_useForestry")

Natural_vegetation$AN <- (Agriculture$Agriculture) -
  (Natural_vegetation$Natural_vegetation)
Natural_vegetation$UN <- (Urban$Urban) -
  (Natural_vegetation$Natural_vegetation)
Natural_vegetation$FN <- (Forestry$Forestry) -
  (Natural_vegetation$Natural_vegetation)


Comparisons <- Natural_vegetation |>
  pivot_longer(
    cols = c("AN", "UN", "FN"),
    names_to = "Comparison",
    values_to = "Ratio"
  )

Comparisons$Comparison <- factor(
  Comparisons$Comparison,
  levels = c("AN", "UN", "FN"),
  labels = c('Agriculture', 'Urban', 'Forestry')
)


Comparisons <- left_join(Comparisons, Metadata, by = 'Dataset_id')

###  fit model to latitudinal patterns####
install.packages(
  "cmdstanr",
  repos = c('https://stan-dev.r-universe.dev', getOption("repos")),
  lib = "~/share/groups/synthesis/Minghua"
)
library(cmdstanr)

# lat_alpha_S_model <- lm(
#   formula = Ratio ~ Latitude,
#   data = Comparisons
# )

lat_alpha_S_model <- brm(
  Ratio ~ Latitude + (1 | Dataset_id),
  data = Latitude_alpha_S_comparisons,
  family = "gaussian",
  iter = 2000,
  warmup = 1000,
  cores = 4,
  chains = 4
)

save(lat_alpha_S_model, file = "lat_alpha_S_model.R")

pp_check(lat_alpha_S_model) +
  scale_x_continuous(trans = 'log')
# mu <- Comparisons |>
#   group_by(Comparison) |>
#   summarise(mu=mean(Ratio))
#
# taxamu <- Comparisons |>
#   group_by(Taxa) |>
#   summarise(mu=mean(Ratio))
#
# continent <- alpha_metrics |>
#   select(c('Dataset_id','Continent'))
#
# continent <- continent[!duplicated(continent$Dataset_id),]
#
# Comparisons <- left_join(Comparisons,continent,by='Dataset_id')
#
# contmu <- Comparisons |>
#   group_by(Continent) |>
#   summarise(mu=mean(Ratio))

# Comparisons$Continent<- factor(Comparisons$Continent,
#                                 levels = c("Europe","Central_America",
#                                            "Africa","North_America",
#                                            "South_America","Oceania",
#                                            "Asia"),
#                                 labels = c("Europe","Central America",
#                                            "Africa","North America",
#                                            "South America","Oceania",
#                                            "Asia"
#                                 ))
### alpha S plot comparison at study level####
studylevel <- Comparisons |>
  group_by(Dataset_id, Comparison) |>
  summarise(
    mean = mean(Ratio),
    lower_ci = mean - qt(0.975, df = n() - 1) * (sd(Ratio) / sqrt(n())),
    upper_ci = mean + qt(0.975, df = n() - 1) * (sd(Ratio) / sqrt(n()))
  )

studylevel <- left_join(studylevel, Metadata, by = 'Dataset_id')
#calculate the absolute value of Latitude
studylevel$Latitude <- abs(as.numeric(studylevel$Latitude))
studylevel$Dataset_id <- as.character(studylevel$Dataset_id)

gradient_colors <- colorRampPalette(c(
  "#274659",
  "#CAAE10",
  "#F2790F",
  "#B93102"
))(length(unique(studylevel$Dataset_id)))
# absolute latitude

alpha_altitude_S <- ggplot() +
  geom_point(
    data = studylevel,
    aes(x = Latitude, y = mean, group = Dataset_id, colour = Dataset_id),
    size = 1,
    position = position_dodge(width = 0.5),
    alpha = 0.3
  ) +
  geom_linerange(
    data = studylevel,
    aes(
      x = Latitude,
      y = mean,
      ymin = lower_ci,
      ymax = upper_ci,
      group = Dataset_id,
      colour = Dataset_id
    ),
    position = position_dodge(width = 0.7),
    alpha = 0.5
  ) +
  geom_smooth(
    data = studylevel,
    aes(x = Latitude, y = mean),
    color = "#B93102",
    method = "glm"
  ) +
  geom_hline(yintercept = 0, linewidth = 0.5, lty = 2) +
  scale_color_manual(values = gradient_colors) +
  labs(x = '', y = 'α', subtitle = "S") +
  ylim(-1.0, 0.5) +
  theme_minimal() +
  theme(
    legend.position = 'none',
    panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 16)
    # axis.title.x = element_text(size=14)
  )

latitude_alpha_S <- lm(mean ~ Latitude, data = studylevel)
summary(latitude_alpha_S)

# Taxa
studylevel <- studylevel |>
  filter(Taxa %in% c("Macroinvertebrate", "Fish", "Algae"))

taxa_mean <- studylevel |>
  group_by(Taxa) |>
  summarise(
    taxa_mean = mean(mean),
    lower_ci = taxa_mean - qt(0.975, df = n() - 1) * (sd(mean) / sqrt(n())),
    upper_ci = taxa_mean + qt(0.975, df = n() - 1) * (sd(mean) / sqrt(n()))
  )

gradient_colors <- colorRampPalette(c(
  "#274659",
  "#CAAE10",
  "#F2790F",
  "#B93102"
))(length(unique(studylevel$Dataset_id)))

studylevel$Taxa <- factor(
  studylevel$Taxa,
  levels = c("Algae", "Fish", "Macroinvertebrate"),
  labels = c("Algae", "Fish", "Macroinvertebrates")
)

taxa_mean$Taxa <- factor(
  taxa_mean$Taxa,
  levels = c("Algae", "Fish", "Macroinvertebrate"),
  labels = c("Algae", "Fish", "Macroinvertebrates")
)

alpha_taxa_S <- ggplot() +
  geom_point(
    data = studylevel,
    aes(x = mean, y = Taxa, group = Dataset_id, colour = Dataset_id),
    size = 1,
    position = position_dodge(width = 0.7),
    alpha = 0.3
  ) +
  geom_linerange(
    data = studylevel,
    aes(
      x = mean,
      xmin = lower_ci,
      xmax = upper_ci,
      y = Taxa,
      group = Dataset_id,
      colour = Dataset_id
    ),
    position = position_dodge(width = 0.7),
    alpha = 0.5
  ) +
  geom_point(
    data = taxa_mean,
    aes(x = taxa_mean, y = Taxa),
    shape = 21,
    color = "black",
    fill = "white",
    stroke = 2,
    size = 4
  ) +
  geom_linerange(
    data = taxa_mean,
    aes(x = taxa_mean, xmin = lower_ci, xmax = upper_ci, y = Taxa),
    lwd = 2.5,
    position = position_dodge(width = 1)
  ) +
  geom_vline(xintercept = 0, linewidth = 0.5, lty = 2) +
  scale_color_manual(values = gradient_colors) +
  labs(y = 'α', x = '', subtitle = "S") +
  xlim(-0.75, 0.5) +
  theme_minimal() +
  theme(
    legend.position = 'none',
    panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 14)
  )

ggplot() +
  geom_point(
    data = studylevel,
    aes(
      x = mean,
      y = reorder(Dataset_id, Latitude),
      # group_by(Comparison),
      colour = Comparison
    ),
    size = 1,
    position = position_dodge(width = 0.5),
    alpha = 0.3
  ) +
  geom_linerange(
    data = studylevel,
    aes(
      x = mean,
      xmin = lower_ci,
      xmax = upper_ci,
      y = reorder(Dataset_id, Latitude),
      # group_by(Comparison),
      colour = Comparison
    ),
    position = position_dodge(width = 0.7),
    alpha = 0.5
  ) +
  geom_vline(xintercept = 0, size = 0.5, lty = 2) +
  theme(
    # legend.position="none",
    axis.text.y = element_blank()
  ) +
  labs(y = 'Study', x = 'Log Ratio')


#plot comparison#
alpha_S_continent <- ggplot() +
  geom_density_ridges_gradient(
    data = Comparisons,
    aes(
      x = Ratio,
      y = fct_reorder(Continent, Ratio, .fun = median),
      fill = factor(after_stat(x) > 0)
    ),
    scale = 0.9,
    alpha = 0.5,
    linetype = 0
  ) +
  xlim(-1, 1) +
  xlim(-1, 1) +
  geom_vline(
    data = Comparisons,
    aes(xintercept = median(Ratio)),
    size = 0.5,
    alpha = 0.5,
    lty = 2
  ) +
  geom_vline(xintercept = 0, size = 0.5) +
  geom_point(
    data = Comparisons,
    aes(x = Ratio, y = Continent),
    stat = ggstance:::StatSummaryh,
    fun.x = median,
    shape = 21,
    color = "black",
    fill = "white",
    stroke = 2,
    size = 4
  ) +
  labs(y = '', x = '') +
  geom_text(
    data = Comparisons |>
      group_by(Continent) |>
      summarise(Total = n(), Count = sum(Ratio < 0), d = (Count / Total)) |>
      mutate(percentages = sprintf('%.2f', d)) |>
      ungroup() |>
      distinct(Continent, percentages, .keep_all = T),
    aes(x = -0.85, y = Continent, label = paste(percentages)),
    size = 4,
    nudge_y = 0.5,
    parse = T
  ) +
  theme(legend.position = "none") +
  scale_fill_manual(
    values = c(
      "FALSE" = scales::alpha("#F2790F", 0.7),
      "TRUE" = scales::alpha("#CAAE10", 0.7)
    )
  ) +
  scale_y_discrete(labels = scales::wrap_format(9)) +
  theme_minimal() +
  theme(
    legend.position = 'none',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
    # axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 16),
    # axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 16)
  )
alpha_S_continent
ggsave("alpha_S_continent.png")

alpha_S_taxa <- ggplot() +
  geom_density_ridges_gradient(
    data = Comparisons,
    aes(
      x = Ratio,
      y = fct_reorder(Taxa, Ratio, .fun = median),
      fill = factor(after_stat(x) > 0)
    ),
    scale = 0.9,
    alpha = 0.5,
    linetype = 0
  ) +
  xlim(-1, 1) +
  xlim(-1, 1) +
  geom_vline(
    data = Comparisons,
    aes(xintercept = median(Ratio)),
    size = 0.5,
    alpha = 0.5,
    lty = 2
  ) +
  geom_vline(xintercept = 0, size = 0.5) +
  geom_point(
    data = Comparisons,
    aes(x = Ratio, y = Taxa),
    stat = ggstance:::StatSummaryh,
    fun.x = median,
    shape = 21,
    color = "black",
    fill = "white",
    stroke = 2,
    size = 4
  ) +
  labs(y = '', x = '') +
  geom_text(
    data = Comparisons |>
      group_by(Taxa) |>
      summarise(Total = n(), Count = sum(Ratio < 0), d = (Count / Total)) |>
      mutate(percentages = sprintf('%.2f', d)) |>
      ungroup() |>
      distinct(Taxa, percentages, .keep_all = T),
    aes(x = -0.85, y = Taxa, label = paste(percentages)),
    size = 4,
    nudge_y = 0.5,
    parse = T
  ) +
  theme(legend.position = "none") +
  scale_fill_manual(
    values = c(
      "FALSE" = scales::alpha("#F2790F", 0.7),
      "TRUE" = scales::alpha("#CAAE10", 0.7)
    )
  ) +
  scale_y_discrete(labels = scales::wrap_format(9)) +
  theme_minimal() +
  theme(
    legend.position = 'none',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
    # axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 16),
    # axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 16)
  )
alpha_S_taxa

cowplot::plot_grid(alpha_S_continent, alpha_S_taxa, labels = c('A', 'B')) +
  cowplot::draw_label(y = 0.02, label = 'Log Rotio', size = 18)
ggsave("random_continent_taxa.png", width = 300, height = 200, units = 'mm')

alpha_S_land <- ggplot() +
  geom_density_ridges_gradient(
    data = Comparisons,
    aes(
      x = Ratio,
      y = fct_reorder(Comparison, Ratio, .fun = median),
      fill = factor(after_stat(x) > 0)
    ),
    scale = 0.9,
    alpha = 0.5,
    linetype = 0
  ) +
  xlim(-1, 1) +
  xlim(-1, 1) +
  geom_vline(
    data = Comparisons,
    aes(xintercept = median(Ratio)),
    size = 0.5,
    alpha = 0.5,
    lty = 2
  ) +
  geom_vline(xintercept = 0, size = 0.5) +
  geom_point(
    data = Comparisons,
    aes(x = Ratio, y = Comparison),
    stat = ggstance:::StatSummaryh,
    fun.x = median,
    shape = 21,
    color = "black",
    fill = "white",
    stroke = 2,
    size = 4
  ) +
  labs(y = '', x = '') +
  geom_text(
    data = Comparisons |>
      group_by(Comparison) |>
      summarise(Total = n(), Count = sum(Ratio < 0), d = (Count / Total)) |>
      mutate(percentages = sprintf('%.2f', d)) |>
      ungroup() |>
      distinct(Comparison, percentages, .keep_all = T),
    aes(x = -0.85, y = Comparison, label = paste(percentages)),
    size = 4,
    nudge_y = 0.5,
    parse = T
  ) +
  theme(legend.position = "none") +
  scale_fill_manual(
    values = c(
      "FALSE" = scales::alpha("#F2790F", 0.7),
      "TRUE" = scales::alpha("#CAAE10", 0.7)
    )
  ) +
  scale_y_discrete(labels = scales::wrap_format(9)) +
  theme_minimal() +
  theme(
    legend.position = 'none',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
    # axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 16),
    # axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 16)
  )
alpha_S_land
ggsave("alpha_S_land.png")

## alpha N#####
alpha_N_model <- alpha_N_model2
rm(alpha_N_model2)
resid <- residuals(alpha_N_model, type = 'pearson', method = 'predict') |>
  dplyr::as_tibble() |>
  bind_cols(alpha_N_model$data)

fitted <- fitted(alpha_N_model, re_formula = NA)
predict <- predict(alpha_N_model)

resid$fitted <- fitted[, 'Estimate']
resid$predict <- predict[, 'Estimate']

plot(resid$Estimate ~ resid$fitted, ylab = 'Pearson residual')

plot(
  resid$Estimate ~ resid$Land_use,
  ylab = 'Pearson residual',
  xlab = "Land_use"
)

Natural_vegetation <- alpha_N_model |>
  spread_draws(
    b_Land_useNatural_vegetation,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  ) |>
  print(n = 100)

Natural_vegetation <- Natural_vegetation |>
  mutate(Natural_vegetation = b_Land_useNatural_vegetation + r_Dataset_id) |>
  filter(Land_use == "Land_useNatural_vegetation")

Agriculture <- alpha_N_model |>
  spread_draws(
    b_Land_useAgriculture,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  ) |>
  print(n = 100)
Agriculture <- Agriculture |>
  mutate(Agriculture = b_Land_useAgriculture + r_Dataset_id) |>
  filter(Land_use == "Land_useAgriculture")

Urban <- alpha_N_model |>
  spread_draws(
    b_Land_useUrban,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  )
Urban <- Urban |>
  mutate(Urban = b_Land_useUrban + r_Dataset_id) |>
  filter(Land_use == "Land_useUrban")

Forestry <- alpha_N_model |>
  spread_draws(
    b_Land_useForestry,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  )
Forestry <- Forestry |>
  mutate(Forestry = b_Land_useForestry + r_Dataset_id) |>
  filter(Land_use == "Land_useForestry")

Natural_vegetation$AN <- (Agriculture$Agriculture) -
  (Natural_vegetation$Natural_vegetation)
Natural_vegetation$UN <- (Urban$Urban) -
  (Natural_vegetation$Natural_vegetation)
Natural_vegetation$FN <- (Forestry$Forestry) -
  (Natural_vegetation$Natural_vegetation)
# Natural_vegetation$MN<-(Mining$Mining)-
#   (Natural_vegetation$Natural_vegetation)
# Natural_vegetation$UA<-(Urban$Urban)-
#   (Agriculture$Agriculture)

Comparisons <- Natural_vegetation |>
  pivot_longer(
    cols = c("AN", "UN", "FN"),
    names_to = "Comparison",
    values_to = "Ratio"
  )

Comparisons$Comparison <- factor(
  Comparisons$Comparison,
  levels = c("AN", "UN", "FN"),
  labels = c(
    'Agriculture/Natural vegetation',
    'Urban/Natural vegetation',
    'Forestry/Natural vegetation'
  )
)

Comparisons <- left_join(Comparisons, Metadata, by = 'Dataset_id')

# mu <- Comparisons |>
#   group_by(Comparison) |>
#   summarise(mu=mean(Ratio))
#
# taxamu <- Comparisons |>
#   group_by(Taxa) |>
#   summarise(mu=mean(Ratio))
#
# continent <- alpha_metrics |>
#   select(c('Dataset_id','Continent'))
#
# continent <- continent[!duplicated(continent$Dataset_id),]
#
# Comparisons <- left_join(Comparisons,continent,by='Dataset_id')
#
# contmu <- Comparisons |>
#   group_by(Continent) |>
#   summarise(mu=mean(Ratio))

Comparisons$Continent <- factor(
  Comparisons$Continent,
  levels = c(
    "Europe",
    "Central_America",
    "Africa",
    "North_America",
    "South_America",
    "Oceania",
    "Asia"
  ),
  labels = c(
    "Europe",
    "Central America",
    "Africa",
    "North America",
    "South America",
    "Oceania",
    "Asia"
  )
)

### alpha N plot comparison at study level####
studylevel <- Comparisons |>
  group_by(Dataset_id, Comparison) |>
  summarise(
    mean = mean(Ratio),
    lower_ci = mean - qt(0.975, df = n() - 1) * (sd(Ratio) / sqrt(n())),
    upper_ci = mean + qt(0.975, df = n() - 1) * (sd(Ratio) / sqrt(n()))
  )

studylevel <- left_join(studylevel, Metadata, by = 'Dataset_id')
#calculate the absolute value of Latitude
studylevel$Latitude <- abs(as.numeric(studylevel$Latitude))
studylevel$Dataset_id <- as.character(studylevel$Dataset_id)

gradient_colors <- colorRampPalette(c(
  "#274659",
  "#CAAE10",
  "#F2790F",
  "#B93102"
))(length(unique(studylevel$Dataset_id)))
# absolute latitude

alpha_altitude_N <- ggplot() +
  geom_point(
    data = studylevel,
    aes(x = Latitude, y = mean, group = Dataset_id, colour = Dataset_id),
    size = 1,
    position = position_dodge(width = 0.5),
    alpha = 0.3
  ) +
  geom_linerange(
    data = studylevel,
    aes(
      x = Latitude,
      y = mean,
      ymin = lower_ci,
      ymax = upper_ci,
      group = Dataset_id,
      colour = Dataset_id
    ),
    position = position_dodge(width = 0.7),
    alpha = 0.5
  ) +
  geom_smooth(
    data = studylevel,
    aes(x = Latitude, y = mean),
    color = "#B93102",
    method = "glm"
  ) +
  scale_color_manual(values = gradient_colors) +
  labs(x = 'Absolute latitude', y = 'α-N') +
  ylim(-0.5, 0.5) +
  theme_minimal() +
  theme(
    legend.position = 'none',
    panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14)
  )

# # Taxa
# studylevel <- studylevel |>
#   filter(Taxa %in% c("Macroinvertebrate","Fish","Algae"))
#
# taxa_mean <- studylevel |>
#   group_by(Taxa) |>
#   summarise(
#     taxa_mean = mean(mean),
#     lower_ci = taxa_mean - qt(0.975, df = n() - 1) * (sd(mean) / sqrt(n())),
#     upper_ci = taxa_mean + qt(0.975, df = n() - 1) * (sd(mean) / sqrt(n()))
#   )
#
# studylevel$Taxa<- factor(studylevel$Taxa,
#                          levels = c("Macrophyte",
#                                     "Amphibian",
#                                     "Algae",
#                                     "Zooplankton",
#                                     "Fish",
#                                     "Macroinvertebrate"
#                          ))
# taxa_mean$Taxa<- factor(taxa_mean$Taxa,
#                         levels = c("Macrophyte",
#                                    "Amphibian",
#                                    "Algae",
#                                    "Zooplankton",
#                                    "Fish",
#                                    "Macroinvertebrate"))
#
# alpha_taxa_N <- ggplot() +
#   geom_point(data=studylevel, aes(x =mean,y = Taxa,
#                                   group = Dataset_id,
#                                   colour = Dataset_id),
#              size = 1,
#              position = position_dodge(width = 0.7),
#              alpha=0.3)+
#   geom_linerange(data = studylevel,
#                  aes(x =mean, xmin = lower_ci, xmax = upper_ci,
#                      y= Taxa,
#                      group = Dataset_id,
#                      colour = Dataset_id),
#                  position = position_dodge(width = 0.7),
#                  alpha = 0.5)+
#   geom_point(data=taxa_mean, aes(x =taxa_mean, y = Taxa),
#              shape = 21,
#              color = "black",
#              fill = "white",
#              stroke = 2,
#              size = 4)+
#   geom_linerange(data = taxa_mean,
#                  aes(x = taxa_mean, xmin =lower_ci,  xmax =upper_ci, y= Taxa),
#                  lwd = 2.5, position = position_dodge(width = 1))+
#   geom_vline(xintercept = 0,  size = 0.5,lty = 2)+
#   scale_color_manual(values = gradient_colors)+
#   labs(y = '',
#        x = 'α-N')+
#   xlim(-0.5,0.5)+
#   theme_minimal() +
#   theme(legend.position = 'none',
#         panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
#         panel.grid.minor.y = element_blank(),
#         axis.text.x = element_text(size = 12),
#         axis.text.y = element_text(size = 12),
#         axis.title.y = element_text(size=14),
#         axis.title.x = element_text(size=14))
#

#plot comparison#
alpha_N_continent <- ggplot() +
  geom_density_ridges_gradient(
    data = Comparisons,
    aes(
      x = Ratio,
      y = fct_reorder(Continent, Ratio, .fun = median),
      fill = factor(after_stat(x) > 0)
    ),
    scale = 0.9,
    alpha = 0.5,
    linetype = 0
  ) +
  xlim(-1, 1) +
  xlim(-1, 1) +
  geom_vline(
    data = Comparisons,
    aes(xintercept = median(Ratio)),
    size = 0.5,
    alpha = 0.5,
    lty = 2
  ) +
  geom_vline(xintercept = 0, size = 0.5) +
  geom_point(
    data = Comparisons,
    aes(x = Ratio, y = Continent),
    stat = ggstance:::StatSummaryh,
    fun.x = median,
    shape = 21,
    color = "black",
    fill = "white",
    stroke = 2,
    size = 4
  ) +
  labs(y = '', x = '') +
  geom_text(
    data = Comparisons |>
      group_by(Continent) |>
      summarise(Total = n(), Count = sum(Ratio < 0), d = (Count / Total)) |>
      mutate(percentages = sprintf('%.2f', d)) |>
      ungroup() |>
      distinct(Continent, percentages, .keep_all = T),
    aes(x = -0.85, y = Continent, label = paste(percentages)),
    size = 4,
    nudge_y = 0.5,
    parse = T
  ) +
  theme(legend.position = "none") +
  scale_fill_manual(
    values = c(
      "FALSE" = scales::alpha("#F2790F", 0.7),
      "TRUE" = scales::alpha("#CAAE10", 0.7)
    )
  ) +
  scale_y_discrete(labels = scales::wrap_format(9)) +
  theme_minimal() +
  theme(
    legend.position = 'none',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
    # axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 16),
    # axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 16)
  )
alpha_N_continent
ggsave("alpha_N_continent.png")

alpha_N_taxa <- ggplot() +
  geom_density_ridges_gradient(
    data = Comparisons,
    aes(
      x = Ratio,
      y = fct_reorder(Taxa, Ratio, .fun = median),
      fill = factor(after_stat(x) > 0)
    ),
    scale = 0.9,
    alpha = 0.5,
    linetype = 0
  ) +
  xlim(-1, 1) +
  xlim(-1, 1) +
  geom_vline(
    data = Comparisons,
    aes(xintercept = median(Ratio)),
    size = 0.5,
    alpha = 0.5,
    lty = 2
  ) +
  geom_vline(xintercept = 0, size = 0.5) +
  geom_point(
    data = Comparisons,
    aes(x = Ratio, y = Taxa),
    stat = ggstance:::StatSummaryh,
    fun.x = median,
    shape = 21,
    color = "black",
    fill = "white",
    stroke = 2,
    size = 4
  ) +
  labs(y = '', x = '') +
  geom_text(
    data = Comparisons |>
      group_by(Taxa) |>
      summarise(Total = n(), Count = sum(Ratio < 0), d = (Count / Total)) |>
      mutate(percentages = sprintf('%.2f', d)) |>
      ungroup() |>
      distinct(Taxa, percentages, .keep_all = T),
    aes(x = -0.85, y = Taxa, label = paste(percentages)),
    size = 4,
    nudge_y = 0.5,
    parse = T
  ) +
  theme(legend.position = "none") +
  scale_fill_manual(
    values = c(
      "FALSE" = scales::alpha("#F2790F", 0.7),
      "TRUE" = scales::alpha("#CAAE10", 0.7)
    )
  ) +
  scale_y_discrete(labels = scales::wrap_format(9)) +
  theme_minimal() +
  theme(
    legend.position = 'none',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
    # axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 16),
    # axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 16)
  )
alpha_N_taxa

cowplot::plot_grid(alpha_N_continent, alpha_N_taxa, labels = c('A', 'B')) +
  cowplot::draw_label(y = 0.02, label = 'Log Rotio', size = 18)
ggsave("random_continent_taxa_N.png", width = 300, height = 200, units = 'mm')

alpha_N_land <- ggplot() +
  geom_density_ridges_gradient(
    data = Comparisons,
    aes(
      x = Ratio,
      y = fct_reorder(Comparison, Ratio, .fun = median),
      fill = factor(after_stat(x) > 0)
    ),
    scale = 0.9,
    alpha = 0.5,
    linetype = 0
  ) +
  xlim(-1, 1) +
  xlim(-1, 1) +
  geom_vline(
    data = Comparisons,
    aes(xintercept = median(Ratio)),
    size = 0.5,
    alpha = 0.5,
    lty = 2
  ) +
  geom_vline(xintercept = 0, size = 0.5) +
  geom_point(
    data = Comparisons,
    aes(x = Ratio, y = Comparison),
    stat = ggstance:::StatSummaryh,
    fun.x = median,
    shape = 21,
    color = "black",
    fill = "white",
    stroke = 2,
    size = 4
  ) +
  labs(y = '', x = '') +
  geom_text(
    data = Comparisons |>
      group_by(Comparison) |>
      summarise(Total = n(), Count = sum(Ratio < 0), d = (Count / Total)) |>
      mutate(percentages = sprintf('%.2f', d)) |>
      ungroup() |>
      distinct(Comparison, percentages, .keep_all = T),
    aes(x = -0.85, y = Comparison, label = paste(percentages)),
    size = 4,
    nudge_y = 0.5,
    parse = T
  ) +
  theme(legend.position = "none") +
  scale_fill_manual(
    values = c(
      "FALSE" = scales::alpha("#F2790F", 0.7),
      "TRUE" = scales::alpha("#CAAE10", 0.7)
    )
  ) +
  scale_y_discrete(labels = scales::wrap_format(9)) +
  theme_minimal() +
  theme(
    legend.position = 'none',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
    # axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 16),
    # axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 16)
  )
alpha_N_land
ggsave("alpha_N_land.png")

## alpha spie ####
alpha_Spie_model <- alpha_Spie_model2
rm(alpha_Spie_model2)

Natural_vegetation <- alpha_Spie_model |>
  spread_draws(
    b_Land_useNatural_vegetation,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  )

Natural_vegetation <- Natural_vegetation |>
  mutate(Natural_vegetation = b_Land_useNatural_vegetation + r_Dataset_id) |>
  filter(Land_use == "Land_useNatural_vegetation")

Agriculture <- alpha_Spie_model |>
  spread_draws(
    b_Land_useAgriculture,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  ) |>
  print(n = 100)
Agriculture <- Agriculture |>
  mutate(Agriculture = b_Land_useAgriculture + r_Dataset_id) |>
  filter(Land_use == "Land_useAgriculture")

Urban <- alpha_Spie_model |>
  spread_draws(
    b_Land_useUrban,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  )
Urban <- Urban |>
  mutate(Urban = b_Land_useUrban + r_Dataset_id) |>
  filter(Land_use == "Land_useUrban")

Forestry <- alpha_Spie_model |>
  spread_draws(
    b_Land_useForestry,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  )
Forestry <- Forestry |>
  mutate(Forestry = b_Land_useForestry + r_Dataset_id) |>
  filter(Land_use == "Land_useForestry")

Natural_vegetation$AN <- (Agriculture$Agriculture) -
  (Natural_vegetation$Natural_vegetation)
Natural_vegetation$UN <- (Urban$Urban) -
  (Natural_vegetation$Natural_vegetation)
Natural_vegetation$FN <- (Forestry$Forestry) -
  (Natural_vegetation$Natural_vegetation)


Comparisons <- Natural_vegetation |>
  pivot_longer(
    cols = c("AN", "UN", "FN"),
    names_to = "Comparison",
    values_to = "Ratio"
  )

Comparisons$Comparison <- factor(
  Comparisons$Comparison,
  levels = c("AN", "UN", "FN"),
  labels = c(
    'Agriculture/Natural vegetation',
    'Urban/Natural vegetation',
    'Forestry/Natural vegetation'
  )
)

Comparisons <- left_join(Comparisons, Metadata, by = 'Dataset_id')


### alpha Spie plot comparison at study level####
studylevel <- Comparisons |>
  group_by(Dataset_id, Comparison) |>
  summarise(
    mean = mean(Ratio),
    lower_ci = mean - qt(0.975, df = n() - 1) * (sd(Ratio) / sqrt(n())),
    upper_ci = mean + qt(0.975, df = n() - 1) * (sd(Ratio) / sqrt(n()))
  )

studylevel <- left_join(studylevel, Metadata, by = 'Dataset_id')
#calculate the absolute value of Latitude
studylevel$Latitude <- abs(as.numeric(studylevel$Latitude))
studylevel$Dataset_id <- as.character(studylevel$Dataset_id)

gradient_colors <- colorRampPalette(c(
  "#274659",
  "#CAAE10",
  "#F2790F",
  "#B93102"
))(length(unique(studylevel$Dataset_id)))
# absolute latitude

alpha_altitude_Spie <- ggplot() +
  geom_point(
    data = studylevel,
    aes(x = Latitude, y = mean, group = Dataset_id, colour = Dataset_id),
    size = 1,
    position = position_dodge(width = 0.5),
    alpha = 0.3
  ) +
  geom_linerange(
    data = studylevel,
    aes(
      x = Latitude,
      y = mean,
      ymin = lower_ci,
      ymax = upper_ci,
      group = Dataset_id,
      colour = Dataset_id
    ),
    position = position_dodge(width = 0.7),
    alpha = 0.5
  ) +
  geom_smooth(
    data = studylevel,
    aes(x = Latitude, y = mean),
    color = "#B93102",
    method = "glm"
  ) +
  geom_hline(yintercept = 0, linewidth = 0.5, lty = 2) +
  scale_color_manual(values = gradient_colors) +
  labs(x = '', y = '', subtitle = "Spie") +
  ylim(-1.0, 0.5) +
  theme_minimal() +
  theme(
    legend.position = 'none',
    panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14)
    # axis.title.x = element_text(size=14)
  )
latitude_alpha_Spie <- lm(mean ~ Latitude, data = studylevel)
summary(latitude_alpha_Spie)

# Taxa
studylevel <- studylevel |>
  filter(Taxa %in% c("Macroinvertebrate", "Fish", "Algae"))

taxa_mean <- studylevel |>
  group_by(Taxa) |>
  summarise(
    taxa_mean = mean(mean),
    lower_ci = taxa_mean - qt(0.975, df = n() - 1) * (sd(mean) / sqrt(n())),
    upper_ci = taxa_mean + qt(0.975, df = n() - 1) * (sd(mean) / sqrt(n()))
  )

studylevel$Taxa <- factor(
  studylevel$Taxa,
  levels = c("Algae", "Fish", "Macroinvertebrate"),
  labels = c("Algae", "Fish", "Macroinvertebrates")
)

taxa_mean$Taxa <- factor(
  taxa_mean$Taxa,
  levels = c("Algae", "Fish", "Macroinvertebrate"),
  labels = c("Algae", "Fish", "Macroinvertebrates")
)


alpha_taxa_Spie <- ggplot() +
  geom_point(
    data = studylevel,
    aes(x = mean, y = Taxa, group = Dataset_id, colour = Dataset_id),
    size = 1,
    position = position_dodge(width = 0.7),
    alpha = 0.3
  ) +
  geom_linerange(
    data = studylevel,
    aes(
      x = mean,
      xmin = lower_ci,
      xmax = upper_ci,
      y = Taxa,
      group = Dataset_id,
      colour = Dataset_id
    ),
    position = position_dodge(width = 0.7),
    alpha = 0.5
  ) +
  geom_point(
    data = taxa_mean,
    aes(x = taxa_mean, y = Taxa),
    shape = 21,
    color = "black",
    fill = "white",
    stroke = 2,
    size = 4
  ) +
  geom_linerange(
    data = taxa_mean,
    aes(x = taxa_mean, xmin = lower_ci, xmax = upper_ci, y = Taxa),
    lwd = 2.5,
    position = position_dodge(width = 1)
  ) +
  geom_vline(xintercept = 0, size = 0.5, lty = 2) +
  scale_color_manual(values = gradient_colors) +
  labs(y = '', x = '', subtitle = "Spie") +
  xlim(-0.75, 0.5) +
  theme_minimal() +
  theme(
    legend.position = 'none',
    panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14)
  )

## gamma S ####
gamma_S_model <- readRDS("results/models/gamma_S_model.rds")

Natural_vegetation <- gamma_S_model |>
  tidybayes::spread_draws(
    b_Land_useNatural_vegetation,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  )

Natural_vegetation <- Natural_vegetation |>
  mutate(Natural_vegetation = b_Land_useNatural_vegetation + r_Dataset_id) |>
  filter(Land_use == "Land_useNatural_vegetation")

Agriculture <- gamma_S_model |>
  tidybayes::spread_draws(
    b_Land_useAgriculture,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  ) |>
  print(n = 100)
Agriculture <- Agriculture |>
  mutate(Agriculture = b_Land_useAgriculture + r_Dataset_id) |>
  filter(Land_use == "Land_useAgriculture")

Urban <- gamma_S_model |>
  tidybayes::spread_draws(
    b_Land_useUrban,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  )
Urban <- Urban |>
  mutate(Urban = b_Land_useUrban + r_Dataset_id) |>
  filter(Land_use == "Land_useUrban")

Forestry <- gamma_S_model |>
  tidybayes::spread_draws(
    b_Land_useForestry,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  )
Forestry <- Forestry |>
  mutate(Forestry = b_Land_useForestry + r_Dataset_id) |>
  filter(Land_use == "Land_useForestry")

Natural_vegetation$AN <- (Agriculture$b_Land_useAgriculture) -
  (Natural_vegetation$b_Land_useNatural_vegetation)
Natural_vegetation$UN <- (Urban$b_Land_useUrban) -
  (Natural_vegetation$b_Land_useNatural_vegetation)
Natural_vegetation$FN <- (Forestry$b_Land_useForestry) -
  (Natural_vegetation$b_Land_useNatural_vegetation)

Comparisons <- Natural_vegetation |>
  pivot_longer(
    cols = c("AN", "UN", "FN"),
    names_to = "Comparison",
    values_to = "Ratio"
  )

Comparisons$Comparison <- factor(
  Comparisons$Comparison,
  levels = c("AN", "UN", "FN"),
  labels = c(
    'Agriculture/Natural vegetation',
    'Urban/Natural vegetation',
    'Forestry/Natural vegetation'
  )
)

Comparisons <- left_join(Comparisons, Metadata, by = 'Dataset_id')

### gamma S plot comparison at study level####
studylevel <- Comparisons |>
  group_by(Dataset_id, Comparison) |>
  summarise(
    mean = mean(Ratio),
    lower_ci = mean - qt(0.975, df = n() - 1) * (sd(Ratio) / sqrt(n())),
    upper_ci = mean + qt(0.975, df = n() - 1) * (sd(Ratio) / sqrt(n()))
  )

studylevel <- left_join(studylevel, Metadata, by = 'Dataset_id')
#calculate the absolute value of Latitude
studylevel$Latitude <- abs(as.numeric(studylevel$Latitude))
studylevel$Dataset_id <- as.character(studylevel$Dataset_id)

gradient_colors <- colorRampPalette(c(
  "#274659",
  "#CAAE10",
  "#F2790F",
  "#B93102"
))(length(unique(studylevel$Dataset_id)))
# absolute latitude

gamma_altitude_S <- ggplot() +
  geom_point(
    data = studylevel,
    aes(x = Latitude, y = mean, group = Dataset_id, colour = Dataset_id),
    size = 1,
    position = position_dodge(width = 0.5),
    alpha = 0.3
  ) +
  geom_linerange(
    data = studylevel,
    aes(
      x = Latitude,
      y = mean,
      ymin = lower_ci,
      ymax = upper_ci,
      group = Dataset_id,
      colour = Dataset_id
    ),
    position = position_dodge(width = 0.7),
    alpha = 0.5
  ) +
  geom_smooth(
    data = studylevel,
    aes(x = Latitude, y = mean),
    color = "#B93102",
    method = "glm"
  ) +
  geom_hline(yintercept = 0, linewidth = 0.5, lty = 2) +
  scale_color_manual(values = gradient_colors) +
  labs(x = '', y = 'γ') +
  ylim(-1.0, 0.5) +
  theme_minimal() +
  theme(
    legend.position = 'none',
    panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
    panel.grid.minor.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 16)
    # axis.title.x = element_text(size=14)
  )

latitude_gamma_S <- lm(mean ~ Latitude, data = studylevel)
summary(latitude_gamma_S)

# Taxa
taxa_mean <- studylevel |>
  group_by(Taxa) |>
  summarise(
    taxa_mean = mean(mean),
    lower_ci = taxa_mean - qt(0.975, df = n() - 1) * (sd(mean) / sqrt(n())),
    upper_ci = taxa_mean + qt(0.975, df = n() - 1) * (sd(mean) / sqrt(n()))
  )

studylevel$Taxa <- factor(
  studylevel$Taxa,
  levels = c("Algae", "Fish", "Macroinvertebrate"),
  labels = c("Algae", "Fish", "Macroinvertebrates")
)

taxa_mean$Taxa <- factor(
  taxa_mean$Taxa,
  levels = c("Algae", "Fish", "Macroinvertebrate"),
  labels = c("Algae", "Fish", "Macroinvertebrates")
)


gamma_taxa_S <- ggplot() +
  geom_point(
    data = studylevel,
    aes(x = mean, y = Taxa, group = Dataset_id, colour = Dataset_id),
    size = 1,
    position = position_dodge(width = 0.7),
    alpha = 0.3
  ) +
  geom_linerange(
    data = studylevel,
    aes(
      x = mean,
      xmin = lower_ci,
      xmax = upper_ci,
      y = Taxa,
      group = Dataset_id,
      colour = Dataset_id
    ),
    position = position_dodge(width = 0.7),
    alpha = 0.5
  ) +
  geom_point(
    data = taxa_mean,
    aes(x = taxa_mean, y = Taxa),
    shape = 21,
    color = "black",
    fill = "white",
    stroke = 2,
    size = 4
  ) +
  geom_linerange(
    data = taxa_mean,
    aes(x = taxa_mean, xmin = lower_ci, xmax = upper_ci, y = Taxa),
    lwd = 2.5,
    position = position_dodge(width = 1)
  ) +
  geom_vline(xintercept = 0, size = 0.5, lty = 2) +
  scale_color_manual(values = gradient_colors) +
  labs(y = 'γ', x = '') +
  xlim(-0.75, 0.5) +
  theme_minimal() +
  theme(
    legend.position = 'none',
    panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 14)
  )

## gamma Spie ####
gamma_Spie_model <- readRDS(file = "results/models/gamma_Spie_model.rds")

Natural_vegetation <- gamma_Spie_model |>
  spread_draws(
    b_Land_useNatural_vegetation,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  )

Natural_vegetation <- Natural_vegetation |>
  mutate(Natural_vegetation = b_Land_useNatural_vegetation + r_Dataset_id) |>
  filter(Land_use == "Land_useNatural_vegetation")

Agriculture <- gamma_Spie_model |>
  spread_draws(
    b_Land_useAgriculture,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  ) |>
  print(n = 100)
Agriculture <- Agriculture |>
  mutate(Agriculture = b_Land_useAgriculture + r_Dataset_id) |>
  filter(Land_use == "Land_useAgriculture")

Urban <- gamma_Spie_model |>
  spread_draws(
    b_Land_useUrban,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  )
Urban <- Urban |>
  mutate(Urban = b_Land_useUrban + r_Dataset_id) |>
  filter(Land_use == "Land_useUrban")

Forestry <- gamma_Spie_model |>
  spread_draws(
    b_Land_useForestry,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  )
Forestry <- Forestry |>
  mutate(Forestry = b_Land_useForestry + r_Dataset_id) |>
  filter(Land_use == "Land_useForestry")

Natural_vegetation$AN <- (Agriculture$Agriculture) -
  (Natural_vegetation$Natural_vegetation)
Natural_vegetation$UN <- (Urban$Urban) -
  (Natural_vegetation$Natural_vegetation)
Natural_vegetation$FN <- (Forestry$Forestry) -
  (Natural_vegetation$Natural_vegetation)

Comparisons <- Natural_vegetation |>
  pivot_longer(
    cols = c("AN", "UN", "FN"),
    names_to = "Comparison",
    values_to = "Ratio"
  )

Comparisons$Comparison <- factor(
  Comparisons$Comparison,
  levels = c("AN", "UN", "FN"),
  labels = c(
    'Agriculture/Natural vegetation',
    'Urban/Natural vegetation',
    'Forestry/Natural vegetation'
  )
)

Comparisons <- left_join(Comparisons, Metadata, by = 'Dataset_id')

### gamma Spie plot comparison at study level####
studylevel <- Comparisons |>
  group_by(Dataset_id, Comparison) |>
  summarise(
    mean = mean(Ratio),
    lower_ci = mean - qt(0.975, df = n() - 1) * (sd(Ratio) / sqrt(n())),
    upper_ci = mean + qt(0.975, df = n() - 1) * (sd(Ratio) / sqrt(n()))
  )

studylevel <- left_join(studylevel, Metadata, by = 'Dataset_id')
#calculate the absolute value of Latitude
studylevel$Latitude <- abs(as.numeric(studylevel$Latitude))
studylevel$Dataset_id <- as.character(studylevel$Dataset_id)

gradient_colors <- colorRampPalette(c(
  "#274659",
  "#CAAE10",
  "#F2790F",
  "#B93102"
))(length(unique(studylevel$Dataset_id)))
# absolute latitude

gamma_altitude_Spie <- ggplot() +
  geom_point(
    data = studylevel,
    aes(x = Latitude, y = mean, group = Dataset_id, colour = Dataset_id),
    size = 1,
    position = position_dodge(width = 0.5),
    alpha = 0.3
  ) +
  geom_linerange(
    data = studylevel,
    aes(
      x = Latitude,
      y = mean,
      ymin = lower_ci,
      ymax = upper_ci,
      group = Dataset_id,
      colour = Dataset_id
    ),
    position = position_dodge(width = 0.7),
    alpha = 0.5
  ) +
  geom_smooth(
    data = studylevel,
    aes(x = Latitude, y = mean),
    color = "#B93102",
    method = "glm"
  ) +
  geom_hline(yintercept = 0, linewidth = 0.5, lty = 2) +
  scale_color_manual(values = gradient_colors) +
  labs(x = '', y = '') +
  ylim(-1.0, 0.5) +
  theme_minimal() +
  theme(
    legend.position = 'none',
    panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
    panel.grid.minor.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14)
    # axis.title.x = element_text(size=14)
  )

latitude_gamma_Spie <- lm(mean ~ Latitude, data = studylevel)
summary(latitude_gamma_Spie)

# Taxa
taxa_mean <- studylevel |>
  group_by(Taxa) |>
  summarise(
    taxa_mean = mean(mean),
    lower_ci = taxa_mean - qt(0.975, df = n() - 1) * (sd(mean) / sqrt(n())),
    upper_ci = taxa_mean + qt(0.975, df = n() - 1) * (sd(mean) / sqrt(n()))
  )

studylevel$Taxa <- factor(
  studylevel$Taxa,
  levels = c("Algae", "Fish", "Macroinvertebrate"),
  labels = c("Algae", "Fish", "Macroinvertebrates")
)

taxa_mean$Taxa <- factor(
  taxa_mean$Taxa,
  levels = c("Algae", "Fish", "Macroinvertebrate"),
  labels = c("Algae", "Fish", "Macroinvertebrates")
)


gamma_taxa_Spie <- ggplot() +
  geom_point(
    data = studylevel,
    aes(x = mean, y = Taxa, group = Dataset_id, colour = Dataset_id),
    size = 1,
    position = position_dodge(width = 0.7),
    alpha = 0.3
  ) +
  geom_linerange(
    data = studylevel,
    aes(
      x = mean,
      xmin = lower_ci,
      xmax = upper_ci,
      y = Taxa,
      group = Dataset_id,
      colour = Dataset_id
    ),
    position = position_dodge(width = 0.7),
    alpha = 0.5
  ) +
  geom_point(
    data = taxa_mean,
    aes(x = taxa_mean, y = Taxa),
    shape = 21,
    color = "black",
    fill = "white",
    stroke = 2,
    size = 4
  ) +
  geom_linerange(
    data = taxa_mean,
    aes(x = taxa_mean, xmin = lower_ci, xmax = upper_ci, y = Taxa),
    lwd = 2.5,
    position = position_dodge(width = 1)
  ) +
  geom_vline(xintercept = 0, size = 0.5, lty = 2) +
  scale_color_manual(values = gradient_colors) +
  labs(y = '', x = '') +
  xlim(-0.75, 0.5) +
  theme_minimal() +
  theme(
    legend.position = 'none',
    panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14)
  )

## gamma N ####
gamma_N_model <- gamma_N_model2
rm(gamma_N_model2)
resid <- residuals(gamma_N_model, type = 'pearson', method = 'predict') |>
  dplyr::as_tibble() |>
  bind_cols(gamma_N_model$data)

fitted <- fitted(gamma_N_model, re_formula = NA)
predict <- predict(gamma_N_model)

resid$fitted <- fitted[, 'Estimate']
resid$predict <- predict[, 'Estimate']

plot(resid$Estimate ~ resid$fitted, ylab = 'Pearson residual')

plot(
  resid$Estimate ~ resid$Land_use,
  ylab = 'Pearson residual',
  xlab = "Land_use"
)

Natural_vegetation <- gamma_N_model |>
  spread_draws(
    b_Land_useNatural_vegetation,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  )

Natural_vegetation <- Natural_vegetation |>
  mutate(Natural_vegetation = b_Land_useNatural_vegetation + r_Dataset_id) |>
  filter(Land_use == "Land_useNatural_vegetation")

Agriculture <- gamma_N_model |>
  spread_draws(
    b_Land_useAgriculture,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  ) |>
  print(n = 100)
Agriculture <- Agriculture |>
  mutate(Agriculture = b_Land_useAgriculture + r_Dataset_id) |>
  filter(Land_use == "Land_useAgriculture")

Urban <- gamma_N_model |>
  spread_draws(
    b_Land_useUrban,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  )
Urban <- Urban |>
  mutate(Urban = b_Land_useUrban + r_Dataset_id) |>
  filter(Land_use == "Land_useUrban")

Forestry <- gamma_N_model |>
  spread_draws(
    b_Land_useForestry,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  )
Forestry <- Forestry |>
  mutate(Forestry = b_Land_useForestry + r_Dataset_id) |>
  filter(Land_use == "Land_useForestry")

Natural_vegetation$AN <- (Agriculture$Agriculture) -
  (Natural_vegetation$Natural_vegetation)
Natural_vegetation$UN <- (Urban$Urban) -
  (Natural_vegetation$Natural_vegetation)
Natural_vegetation$FN <- (Forestry$Forestry) -
  (Natural_vegetation$Natural_vegetation)

Comparisons <- Natural_vegetation |>
  pivot_longer(
    cols = c("AN", "UN", "FN"),
    names_to = "Comparison",
    values_to = "Ratio"
  )

Comparisons$Comparison <- factor(
  Comparisons$Comparison,
  levels = c("AN", "UN", "FN"),
  labels = c(
    'Agriculture/Natural vegetation',
    'Urban/Natural vegetation',
    'Forestry/Natural vegetation'
  )
)

Comparisons <- left_join(Comparisons, Metadata, by = 'Dataset_id')

### gamma N plot comparison at study level####
studylevel <- Comparisons |>
  group_by(Dataset_id, Comparison) |>
  summarise(
    mean = mean(Ratio),
    lower_ci = mean - qt(0.975, df = n() - 1) * (sd(Ratio) / sqrt(n())),
    upper_ci = mean + qt(0.975, df = n() - 1) * (sd(Ratio) / sqrt(n()))
  )

studylevel <- left_join(studylevel, Metadata, by = 'Dataset_id')
#calculate the absolute value of Latitude
studylevel$Latitude <- abs(as.numeric(studylevel$Latitude))
studylevel$Dataset_id <- as.character(studylevel$Dataset_id)

gradient_colors <- colorRampPalette(c(
  "#274659",
  "#CAAE10",
  "#F2790F",
  "#B93102"
))(length(unique(studylevel$Dataset_id)))
# absolute latitude

gamma_altitude_N <- ggplot() +
  geom_point(
    data = studylevel,
    aes(x = Latitude, y = mean, group = Dataset_id, colour = Dataset_id),
    size = 1,
    position = position_dodge(width = 0.5),
    alpha = 0.3
  ) +
  geom_linerange(
    data = studylevel,
    aes(
      x = Latitude,
      y = mean,
      ymin = lower_ci,
      ymax = upper_ci,
      group = Dataset_id,
      colour = Dataset_id
    ),
    position = position_dodge(width = 0.7),
    alpha = 0.5
  ) +
  geom_smooth(
    data = studylevel,
    aes(x = Latitude, y = mean),
    color = "#B93102",
    method = "glm"
  ) +
  geom_hline(yintercept = 0, linewidth = 0.5, lty = 2) +
  scale_color_manual(values = gradient_colors) +
  labs(x = 'Absolute latitude', y = 'γ-N') +
  ylim(-1.0, 0.5) +
  theme_minimal() +
  theme(
    legend.position = 'none',
    panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14)
  )

# Taxa
studylevel <- studylevel |>
  filter(Taxa %in% c("Macroinvertebrate", "Fish", "Algae"))

taxa_mean <- studylevel |>
  group_by(Taxa) |>
  summarise(
    taxa_mean = mean(mean),
    lower_ci = taxa_mean - qt(0.975, df = n() - 1) * (sd(mean) / sqrt(n())),
    upper_ci = taxa_mean + qt(0.975, df = n() - 1) * (sd(mean) / sqrt(n()))
  )

studylevel$Taxa <- factor(
  studylevel$Taxa,
  levels = c("Algae", "Fish", "Macroinvertebrate"),
  labels = c("Algae", "Fish", "Macroinvertebrates")
)

taxa_mean$Taxa <- factor(
  taxa_mean$Taxa,
  levels = c("Algae", "Fish", "Macroinvertebrate"),
  labels = c("Algae", "Fish", "Macroinvertebrates")
)


gamma_taxa_N <- ggplot() +
  geom_point(
    data = studylevel,
    aes(x = mean, y = Taxa, group = Dataset_id, colour = Dataset_id),
    size = 1,
    position = position_dodge(width = 0.7),
    alpha = 0.3
  ) +
  geom_linerange(
    data = studylevel,
    aes(
      x = mean,
      xmin = lower_ci,
      xmax = upper_ci,
      y = Taxa,
      group = Dataset_id,
      colour = Dataset_id
    ),
    position = position_dodge(width = 0.7),
    alpha = 0.5
  ) +
  geom_point(
    data = taxa_mean,
    aes(x = taxa_mean, y = Taxa),
    shape = 21,
    color = "black",
    fill = "white",
    stroke = 2,
    size = 4
  ) +
  geom_linerange(
    data = taxa_mean,
    aes(x = taxa_mean, xmin = lower_ci, xmax = upper_ci, y = Taxa),
    lwd = 2.5,
    position = position_dodge(width = 1)
  ) +
  geom_vline(xintercept = 0, size = 0.5, lty = 2) +
  scale_color_manual(values = gradient_colors) +
  labs(y = '', x = 'γ-N') +
  xlim(-0.75, 0.5) +
  theme_minimal() +
  theme(
    legend.position = 'none',
    panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14)
  )

## beta S ####
beta_S_model <- beta_S_model2
rm(beta_S_model2)

Natural_vegetation <- beta_S_model |>
  spread_draws(
    b_Land_useNatural_vegetation,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  )

Natural_vegetation <- Natural_vegetation |>
  mutate(Natural_vegetation = b_Land_useNatural_vegetation + r_Dataset_id) |>
  filter(Land_use == "Land_useNatural_vegetation")

Agriculture <- beta_S_model |>
  spread_draws(
    b_Land_useAgriculture,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  ) |>
  print(n = 100)
Agriculture <- Agriculture |>
  mutate(Agriculture = b_Land_useAgriculture + r_Dataset_id) |>
  filter(Land_use == "Land_useAgriculture")

Urban <- beta_S_model |>
  spread_draws(
    b_Land_useUrban,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  )
Urban <- Urban |>
  mutate(Urban = b_Land_useUrban + r_Dataset_id) |>
  filter(Land_use == "Land_useUrban")

Forestry <- beta_S_model |>
  spread_draws(
    b_Land_useForestry,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  )
Forestry <- Forestry |>
  mutate(Forestry = b_Land_useForestry + r_Dataset_id) |>
  filter(Land_use == "Land_useForestry")

Natural_vegetation$AN <- (Agriculture$Agriculture) -
  (Natural_vegetation$Natural_vegetation)
Natural_vegetation$UN <- (Urban$Urban) -
  (Natural_vegetation$Natural_vegetation)
Natural_vegetation$FN <- (Forestry$Forestry) -
  (Natural_vegetation$Natural_vegetation)

Comparisons <- Natural_vegetation |>
  pivot_longer(
    cols = c("AN", "UN", "FN"),
    names_to = "Comparison",
    values_to = "Ratio"
  )

Comparisons$Comparison <- factor(
  Comparisons$Comparison,
  levels = c("AN", "UN", "FN"),
  labels = c('Agriculture', 'Urban', 'Forestry')
)

Comparisons <- left_join(Comparisons, Metadata, by = 'Dataset_id')

### beta S plot comparison at study level####
studylevel <- Comparisons |>
  group_by(Dataset_id, Comparison) |>
  summarise(
    mean = mean(Ratio),
    lower_ci = mean - qt(0.975, df = n() - 1) * (sd(Ratio) / sqrt(n())),
    upper_ci = mean + qt(0.975, df = n() - 1) * (sd(Ratio) / sqrt(n()))
  )

studylevel <- left_join(studylevel, Metadata, by = 'Dataset_id')
#calculate the absolute value of Latitude
studylevel$Latitude <- abs(as.numeric(studylevel$Latitude))
studylevel$Dataset_id <- as.character(studylevel$Dataset_id)

gradient_colors <- colorRampPalette(c(
  "#274659",
  "#CAAE10",
  "#F2790F",
  "#B93102"
))(length(unique(studylevel$Dataset_id)))
# absolute latitude

beta_altitude_S <- ggplot() +
  geom_point(
    data = studylevel,
    aes(x = Latitude, y = mean, group = Dataset_id, colour = Dataset_id),
    size = 1,
    position = position_dodge(width = 0.5),
    alpha = 0.3
  ) +
  geom_linerange(
    data = studylevel,
    aes(
      x = Latitude,
      y = mean,
      ymin = lower_ci,
      ymax = upper_ci,
      group = Dataset_id,
      colour = Dataset_id
    ),
    position = position_dodge(width = 0.7),
    alpha = 0.5
  ) +
  geom_smooth(
    data = studylevel,
    aes(x = Latitude, y = mean),
    color = "#B93102",
    method = "glm"
  ) +
  geom_hline(yintercept = 0, linewidth = 0.5, lty = 2) +
  scale_color_manual(values = gradient_colors) +
  labs(x = '', y = 'β') +
  ylim(-1.0, 0.5) +
  theme_minimal() +
  theme(
    legend.position = 'none',
    panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
    panel.grid.minor.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 16)
    # axis.title.x = element_text(size=14)
  )

latitude_beta_S <- lm(mean ~ Latitude, data = studylevel)
summary(latitude_beta_S)

# Taxa
taxa_mean <- studylevel |>
  group_by(Taxa) |>
  summarise(
    taxa_mean = mean(mean),
    lower_ci = taxa_mean - qt(0.975, df = n() - 1) * (sd(mean) / sqrt(n())),
    upper_ci = taxa_mean + qt(0.975, df = n() - 1) * (sd(mean) / sqrt(n()))
  )

studylevel$Taxa <- factor(
  studylevel$Taxa,
  levels = c("Algae", "Fish", "Macroinvertebrate"),
  labels = c("Algae", "Fish", "Macroinvertebrates")
)

taxa_mean$Taxa <- factor(
  taxa_mean$Taxa,
  levels = c("Algae", "Fish", "Macroinvertebrate"),
  labels = c("Algae", "Fish", "Macroinvertebrates")
)


beta_taxa_S <- ggplot() +
  geom_point(
    data = studylevel,
    aes(x = mean, y = Taxa, group = Dataset_id, colour = Dataset_id),
    size = 1,
    position = position_dodge(width = 0.7),
    alpha = 0.3
  ) +
  geom_linerange(
    data = studylevel,
    aes(
      x = mean,
      xmin = lower_ci,
      xmax = upper_ci,
      y = Taxa,
      group = Dataset_id,
      colour = Dataset_id
    ),
    position = position_dodge(width = 0.7),
    alpha = 0.5
  ) +
  geom_point(
    data = taxa_mean,
    aes(x = taxa_mean, y = Taxa),
    shape = 21,
    color = "black",
    fill = "white",
    stroke = 2,
    size = 4
  ) +
  geom_linerange(
    data = taxa_mean,
    aes(x = taxa_mean, xmin = lower_ci, xmax = upper_ci, y = Taxa),
    lwd = 2.5,
    position = position_dodge(width = 1)
  ) +
  geom_vline(xintercept = 0, size = 0.5, lty = 2) +
  scale_color_manual(values = gradient_colors) +
  labs(y = 'β', x = '') +
  xlim(-0.75, 0.5) +
  theme_minimal() +
  theme(
    legend.position = 'none',
    panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 14)
  )
## beta Spie #####
beta_Spie_model <- beta_Spie_model2
rm(beta_Spie_model2)

Natural_vegetation <- beta_Spie_model |>
  spread_draws(
    b_Land_useNatural_vegetation,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  )

Natural_vegetation <- Natural_vegetation |>
  mutate(Natural_vegetation = b_Land_useNatural_vegetation + r_Dataset_id) |>
  filter(Land_use == "Land_useNatural_vegetation")

Agriculture <- beta_Spie_model |>
  spread_draws(
    b_Land_useAgriculture,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  ) |>
  print(n = 100)
Agriculture <- Agriculture |>
  mutate(Agriculture = b_Land_useAgriculture + r_Dataset_id) |>
  filter(Land_use == "Land_useAgriculture")

Urban <- beta_Spie_model |>
  spread_draws(
    b_Land_useUrban,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  )
Urban <- Urban |>
  mutate(Urban = b_Land_useUrban + r_Dataset_id) |>
  filter(Land_use == "Land_useUrban")

Forestry <- beta_Spie_model |>
  spread_draws(
    b_Land_useForestry,
    r_Dataset_id[Dataset_id, Land_use],
    ndraws = 1000,
    seed = 111
  )
Forestry <- Forestry |>
  mutate(Forestry = b_Land_useForestry + r_Dataset_id) |>
  filter(Land_use == "Land_useForestry")

Natural_vegetation$AN <- (Agriculture$Agriculture) -
  (Natural_vegetation$Natural_vegetation)
Natural_vegetation$UN <- (Urban$Urban) -
  (Natural_vegetation$Natural_vegetation)
Natural_vegetation$FN <- (Forestry$Forestry) -
  (Natural_vegetation$Natural_vegetation)

Comparisons <- Natural_vegetation |>
  pivot_longer(
    cols = c("AN", "UN", "FN"),
    names_to = "Comparison",
    values_to = "Ratio"
  )

Comparisons$Comparison <- factor(
  Comparisons$Comparison,
  levels = c("AN", "UN", "FN"),
  labels = c(
    'Agriculture/Natural vegetation',
    'Urban/Natural vegetation',
    'Forestry/Natural vegetation'
  )
)

Comparisons <- left_join(Comparisons, Metadata, by = 'Dataset_id')

### beta Spie plot comparison at study level####
studylevel <- Comparisons |>
  group_by(Dataset_id, Comparison) |>
  summarise(
    mean = mean(Ratio),
    lower_ci = mean - qt(0.975, df = n() - 1) * (sd(Ratio) / sqrt(n())),
    upper_ci = mean + qt(0.975, df = n() - 1) * (sd(Ratio) / sqrt(n()))
  )

studylevel <- left_join(studylevel, Metadata, by = 'Dataset_id')
#calculate the absolute value of Latitude
studylevel$Latitude <- abs(as.numeric(studylevel$Latitude))
studylevel$Dataset_id <- as.character(studylevel$Dataset_id)

gradient_colors <- colorRampPalette(c(
  "#274659",
  "#CAAE10",
  "#F2790F",
  "#B93102"
))(length(unique(studylevel$Dataset_id)))
# absolute latitude

beta_altitude_Spie <- ggplot() +
  geom_point(
    data = studylevel,
    aes(x = Latitude, y = mean, group = Dataset_id, colour = Dataset_id),
    size = 1,
    position = position_dodge(width = 0.5),
    alpha = 0.3
  ) +
  geom_linerange(
    data = studylevel,
    aes(
      x = Latitude,
      y = mean,
      ymin = lower_ci,
      ymax = upper_ci,
      group = Dataset_id,
      colour = Dataset_id
    ),
    position = position_dodge(width = 0.7),
    alpha = 0.5
  ) +
  geom_smooth(
    data = studylevel,
    aes(x = Latitude, y = mean),
    color = "#B93102",
    method = "glm"
  ) +
  geom_hline(yintercept = 0, linewidth = 0.5, lty = 2) +
  scale_color_manual(values = gradient_colors) +
  labs(x = '', y = '') +
  ylim(-1.0, 0.5) +
  theme_minimal() +
  theme(
    legend.position = 'none',
    panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
    panel.grid.minor.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14)
    # axis.title.x = element_text(size=14)
  )

latitude_beta_Spie <- lm(mean ~ Latitude, data = studylevel)
summary(latitude_beta_Spie)

# Taxa
taxa_mean <- studylevel |>
  group_by(Taxa) |>
  summarise(
    taxa_mean = mean(mean),
    lower_ci = taxa_mean - qt(0.975, df = n() - 1) * (sd(mean) / sqrt(n())),
    upper_ci = taxa_mean + qt(0.975, df = n() - 1) * (sd(mean) / sqrt(n()))
  )

studylevel$Taxa <- factor(
  studylevel$Taxa,
  levels = c("Algae", "Fish", "Macroinvertebrate"),
  labels = c("Algae", "Fish", "Macroinvertebrates")
)

taxa_mean$Taxa <- factor(
  taxa_mean$Taxa,
  levels = c("Algae", "Fish", "Macroinvertebrate"),
  labels = c("Algae", "Fish", "Macroinvertebrates")
)


beta_taxa_Spie <- ggplot() +
  geom_point(
    data = studylevel,
    aes(x = mean, y = Taxa, group = Dataset_id, colour = Dataset_id),
    size = 1,
    position = position_dodge(width = 0.7),
    alpha = 0.3
  ) +
  geom_linerange(
    data = studylevel,
    aes(
      x = mean,
      xmin = lower_ci,
      xmax = upper_ci,
      y = Taxa,
      group = Dataset_id,
      colour = Dataset_id
    ),
    position = position_dodge(width = 0.7),
    alpha = 0.5
  ) +
  geom_point(
    data = taxa_mean,
    aes(x = taxa_mean, y = Taxa),
    shape = 21,
    color = "black",
    fill = "white",
    stroke = 2,
    size = 4
  ) +
  geom_linerange(
    data = taxa_mean,
    aes(x = taxa_mean, xmin = lower_ci, xmax = upper_ci, y = Taxa),
    lwd = 2.5,
    position = position_dodge(width = 1)
  ) +
  geom_vline(xintercept = 0, size = 0.5, lty = 2) +
  scale_color_manual(values = gradient_colors) +
  labs(y = '', x = '') +
  xlim(-0.75, 0.5) +
  theme_minimal() +
  theme(
    legend.position = 'none',
    panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14)
  )

## combine plots #########
Fig2 <- cowplot::plot_grid(
  alpha_altitude_S,
  alpha_altitude_Spie,
  gamma_altitude_S,
  gamma_altitude_Spie,
  beta_altitude_S,
  beta_altitude_Spie,
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
  "reference_Fig2.png",
  Fig2,
  base_height = 11,
  base_width = 8.8,
  units = "in",
  dpi = 600
)

cowplot::save_plot(
  "reference_Fig2.pdf",
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
  "reference_Fig3.png",
  Fig3,
  base_height = 7,
  base_width = 8.8,
  units = "in",
  dpi = 600
)

cowplot::save_plot(
  "reference_Fig3.pdf",
  Fig3,
  base_height = 7,
  base_width = 8.8,
  units = "in",
  dpi = 600,
  device = cairo_pdf
)
