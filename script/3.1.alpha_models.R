library(tidyverse)
library(brms)
library(tidybayes)
library(ggridges)
set.seed(42)

# models ####
alpha_metrics <- read.csv("results/temp/alpha_metrics.txt", sep="")

#remove mining from our data selection
alpha_metrics <- alpha_metrics |> 
  filter(!Land_use %in% "Mining")

#choose the datasets that have natural vegetation as 'reference'
lu_com_alpha <- alpha_metrics |> 
  summarise(n=n(), .by = c(Land_use, Dataset_id)) |>
  pivot_wider(names_from = Land_use, values_from = n)

lu_com_alpha |> 
  filter(!is.na(Natural_vegetation)) |> 
  summarise(n = n_distinct(Dataset_id))

alpha_dataset_id <- lu_com_alpha |> 
  filter(!is.na(Natural_vegetation)) |> 
  distinct(Dataset_id, .keep_all = TRUE) |> 
  pull(Dataset_id)

alpha_metrics <- alpha_metrics |>
  filter(Dataset_id %in% alpha_dataset_id)

alpha_metrics |> 
  summarise(n = n_distinct(Dataset_id))

alpha_metrics |> 
  summarise(n = n_distinct(Site_id))

#count countries numbers
Metadata <- read.csv(file = "data/Metadata.csv")

alpha_metrics_meta <-left_join(
  x = alpha_metrics,
  y = Metadata,
  by = join_by("Dataset_id", "Taxa"))

data_clean <- dplyr::distinct(alpha_metrics_meta, Latitude.y, Longitude.y, .keep_all = TRUE)

sf::sf_use_s2(FALSE)

country_id <- bdc::bdc_country_from_coordinates(
  data = data_clean,
  lat = "Latitude.y",
  lon = "Longitude.y"
)
country_id|> 
  summarise(n = n_distinct(country))

alpha_metrics$Land_use <- factor(alpha_metrics$Land_use,
                                levels = c("Natural_vegetation","Forestry","Agriculture", "Urban"))


alpha_S_model <- brm(
  S  ~ 0 + Land_use +
    (0 + Land_use | Dataset_id/Block/Site_id/Plot_id), 
  family = 'poisson',
  data = alpha_metrics,
  # control = list(adapt_delta = 0.999, max_treedepth = 12),
  # iter = 4000,
  # warmup = 1000,
  cores = 12,
  chains = 4
)
saveRDS(object = alpha_S_model, file="results/models/alpha_S_model.rds")

alpha_N_model <- brm(
  N  ~ 0 + Land_use +
    (0 + Land_use | Dataset_id/Block/Site_id/Plot_id),
  family = 'lognormal',
  data = alpha_metrics,
  # control = list(adapt_delta = 0.999, max_treedepth = 12),
  # iter = 4000,
  # warmup = 1000,
  cores = 12,
  chains = 4
)
saveRDS(object = alpha_N_model, file="results/models/alpha_N_model.rds")


alpha_Spie_model <- brm(
  Spie  ~ 0 + Land_use +
    (0 + Land_use | Dataset_id/Block/Site_id/Plot_id),
  family = 'lognormal',
  data = alpha_metrics,
  control = list(adapt_delta = 0.999, max_treedepth = 12),
  iter = 4000,
  warmup = 1000,
  cores = 12,
  chains = 4
)
saveRDS(object = alpha_Spie_model,file="results/models/alpha_Spie_model.rds")
pp_check(alpha_Spie_model) +  scale_x_continuous(trans = 'log')

