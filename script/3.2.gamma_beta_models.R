library(tidyverse)
library(brms)
library(tidybayes)
library(ggridges)
set.seed(42)

# Models ####
gamma_bet_metrics <- read.delim(file = 'results/temp/gamma_bet_metrics.txt',
                                sep = " ") |> 
  as_tibble()

gamma_bet_metrics <- gamma_bet_metrics |> 
  filter(!Land_use %in% "Mining")


lu_com_gamma <- gamma_bet_metrics |>
  group_by(Land_use,Dataset_id) |>
  summarise(n=n()) |>
  pivot_wider(names_from = Land_use,values_from = n)

lu_com_gamma |> 
  filter(!is.na(Natural_vegetation)) |> 
  summarise(n = n_distinct(Dataset_id))

gamma_dataset_id <- lu_com_gamma |> 
  filter(!is.na(Natural_vegetation)) |> 
  distinct(Dataset_id, .keep_all = TRUE) |> 
  pull(Dataset_id)

gamma_bet_metrics <- gamma_bet_metrics |> 
  filter(Dataset_id %in% gamma_dataset_id)

gamma_bet_metrics |> 
  distinct(Dataset_id)

gamma_bet_metrics$Land_use <- factor(gamma_bet_metrics$Land_use,
                                     levels = c("Natural_vegetation","Forestry",
                                                "Agriculture", "Urban"))
## gamma S ####
gamma_bet_metrics |> 
  summarise(mu = mean(gamma_S), .by = Land_use)

get_prior(gamma_S  ~ 0 + Land_use + (0 + Land_use|Dataset_id/Block),
          data = gamma_bet_metrics)
gamma_S_model2 <- brm(
  gamma_S  ~ 0 + Land_use + (0 + Land_use|Dataset_id/Block),
  # prior = c(prior('normal(0,1)', class = 'sd'),
  # prior('normal(55,10)', class = 'b')),
  family = 'lognormal',
  data = gamma_bet_metrics,
  control = list(adapt_delta = 0.999, max_treedepth = 12),
  iter = 4000,
  warmup = 1000,
  cores = 12,
  chains = 4
  # ,
  # sample_prior = 'only'
)
saveRDS(gamma_S_model2,file="results/models/gamma_S_model.rds")

## gamma N ####
gamma_N_model2 <- brm(
  gamma_N  ~ 0 + Land_use +
    (0 + Land_use | Dataset_id/Block),
  family = 'lognormal',
  data = gamma_bet_metrics,
  control = list(adapt_delta = 0.999, max_treedepth = 12),
  iter = 4000,
  warmup = 1000,
  cores = 12,
  chains = 4
)
saveRDS(gamma_N_model2,file="results/models/gamma_N_model.rds")

## gamma Spie ####
gamma_Spie_model2 <- brm(
  gamma_Spie  ~ 0 + Land_use +
    (0 + Land_use | Dataset_id/Block),
  family = 'lognormal',
  data = gamma_bet_metrics,
  control = list(adapt_delta = 0.999, max_treedepth = 12),
  iter = 4000,
  warmup = 1000,
  cores = 12,
  chains = 4
)
saveRDS(gamma_Spie_model2,file = "results/models/gamma_Spie_model.rds")
## beta S ####
beta_S_model2 <- brm(
  beta_S  ~ 0 + Land_use +
    (0 + Land_use | Dataset_id/Block),
  family = 'lognormal',
  data = gamma_bet_metrics,
  control = list(adapt_delta = 0.999, max_treedepth = 12),
  iter = 4000,
  warmup = 1000,
  cores = 12,
  chains = 4
)
saveRDS(beta_S_model2,file = "results/models/beta_S_model.rds")

## beta Spie ####
beta_Spie_model2 <- brm(
  beta_Spie  ~ 0 + Land_use +
    (0 + Land_use | Dataset_id/Block),
  family = 'lognormal',
  data = gamma_bet_metrics,
  control = list(adapt_delta = 0.999, max_treedepth = 12),
  iter = 4000,
  warmup = 1000,
  cores = 12,
  chains = 4
)
saveRDS(beta_Spie_model2,file = "results/models/beta_Spie_model.rds")


