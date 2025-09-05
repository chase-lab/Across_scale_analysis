#' load model results and prepares two objects for plotting
#' @param scale  a string "alpha" or "gamma"
#' @param metric a string "N", "S" or "Spie"
#' @returns a list of two tibbles called `fitted` and `study` and the `metric` argument.
#' @examples
#' alpha_N_data <- wrangle_responses(scale = "alpha", metric = "N")
#'

.wrangle_responses <- function(scale, metric) {
  stopifnot(
    "scale must be alpha, beta or gamma" = is.element(
      scale,
      c("alpha", "beta", "gamma")
    )
  )
  stopifnot(
    "metric must be N, S or Spie" = is.element(metric, c("N", "S", "Spie"))
  )
  model <- readRDS(
    file = paste0("results/models/", scale, "_", metric, "_model.rds")
  )

  fitted <- fitted(model, re_formula = NA) |>
    as_tibble() |>
    bind_cols(model$data) |>
    mutate(
      Dataset_id = as.character(Dataset_id),
      Filter = paste(Dataset_id, Land_use, sep = "_")
    )
  if (scale == "gamma") {
    if (metric == "N") {
      fitted <- fitted |> rename(N = gamma_N)
    } else {
      fitted <- fitted |> rename(S = gamma_S)
    }
  }
  study <- coef(model)$Dataset_id |>
    pivot_filter(fitted$Filter) |>
    mutate(Land_use = Land_use |> as.factor())
  return(list(fitted = fitted, study = study, metric = metric))
}

#' load model results and prepares two objects for plotting
#' @param scale  a string "alpha" or "gamma"
#' @param metric a string "N" or "S" "Spie"
#' @returns a list of two tibbles called `fitted` and `study` and the `metric` argument.
#' @examples
#' alpha_N_data <- wrangle_responses(scale = "alpha", metric = "N")
#'
.wrangle_comparisons <- function(scale, metric, study_level = FALSE) {
  stopifnot(
    "scale must be alpha, beta or gamma" = is.element(
      scale,
      c("alpha", "beta", "gamma")
    )
  )
  stopifnot(
    "metric must be N, S or Spie" = is.element(metric, c("N", "S", "Spie"))
  )
  model <- readRDS(
    file = paste0("results/models/", scale, "_", metric, "_model.rds")
  )

  if (isTRUE(study_level)) {
    # spread_draws <- model |>
    #   tidybayes::spread_draws(
    #     b_Land_useNatural_vegetation,
    #     b_Land_useAgriculture,
    #     b_Land_useUrban,
    #     b_Land_useForestry,
    #     r_Dataset_id[Dataset_id, Land_use],
    #     ndraws = 1000,
    #     seed = 111
    #   ) |>
    #   ungroup() |>
    #   mutate(across(
    #     .cols = starts_with("b_Land_use"),
    #     .fns = function(x) x + r_Dataset_id
    #   ))
    Natural_vegetation <- model |>
      tidybayes::spread_draws(
        b_Land_useNatural_vegetation,
        r_Dataset_id[Dataset_id, Land_use],
        ndraws = 1000,
        seed = 111
      )

    Natural_vegetation <- Natural_vegetation |>
      mutate(
        Natural_vegetation = b_Land_useNatural_vegetation + r_Dataset_id
      ) |>
      filter(Land_use == "Land_useNatural_vegetation")

    Agriculture <- model |>
      tidybayes::spread_draws(
        b_Land_useAgriculture,
        r_Dataset_id[Dataset_id, Land_use],
        ndraws = 1000,
        seed = 111
      )
    Agriculture <- Agriculture |>
      mutate(Agriculture = b_Land_useAgriculture + r_Dataset_id) |>
      filter(Land_use == "Land_useAgriculture")

    Urban <- model |>
      tidybayes::spread_draws(
        b_Land_useUrban,
        r_Dataset_id[Dataset_id, Land_use],
        ndraws = 1000,
        seed = 111
      )
    Urban <- Urban |>
      mutate(Urban = b_Land_useUrban + r_Dataset_id) |>
      filter(Land_use == "Land_useUrban")

    Forestry <- model |>
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
      ) |>
      ungroup()

    Comparisons$Comparison <- factor(
      Comparisons$Comparison,
      levels = c("AN", "UN", "FN"),
      labels = c(
        'Agriculture/Natural vegetation',
        'Urban/Natural vegetation',
        'Forestry/Natural vegetation'
      )
    )
  } else {
    spread_draws <- model |>
      tidybayes::spread_draws(
        b_Land_useNatural_vegetation,
        b_Land_useAgriculture,
        b_Land_useUrban,
        b_Land_useForestry,
        ndraws = 1000,
        seed = 111
      )

    Comparisons <- spread_draws |>
      mutate(
        AN = b_Land_useAgriculture - b_Land_useNatural_vegetation,
        UN = b_Land_useUrban - b_Land_useNatural_vegetation,
        FN = b_Land_useForestry - b_Land_useNatural_vegetation
      ) |>
      select(any_of("Dataset_id"), AN, UN, FN) |>
      pivot_longer(
        cols = all_of(c("AN", "UN", "FN")),
        names_to = "Comparison",
        values_to = "Ratio"
      ) |>
      mutate(
        Comparison = factor(
          Comparison,
          levels = c("AN", "UN", "FN"),
          labels = c(
            "Agriculture/Natural vegetation",
            "Urban/Natural vegetation",
            "Forestry/Natural vegetation"
          )
        )
      )
  }

  return(list(Comparisons = Comparisons, metric = metric))
}

.wrangle_comparisons_backup <- function(scale, metric, study_level = FALSE) {
  stopifnot(
    "scale must be alpha, beta or gamma" = is.element(
      scale,
      c("alpha", "beta", "gamma")
    )
  )
  stopifnot(
    "metric must be N, S or Spie" = is.element(metric, c("N", "S", "Spie"))
  )
  model <- readRDS(
    file = paste0("results/models/", scale, "_", metric, "_model.rds")
  )

  # residuals <- stats::residuals(
  #   model,
  #   type = "pearson",
  #   method = "predict"
  # ) |>
  #   as_tibble() |>
  #   bind_cols(model$data)
  #
  # residuals$fitted <- stats::fitted(model, re_formula = NA)[, "Estimate"]
  # residuals$predict <- stats::predict(model)[, "Estimate"]

  if (isTRUE(study_level)) {
    spread_draws <- model |>
      tidybayes::spread_draws(
        b_Land_useNatural_vegetation,
        b_Land_useAgriculture,
        b_Land_useUrban,
        b_Land_useForestry,
        r_Dataset_id[Dataset_id, Land_use],
        ndraws = 1000,
        seed = 111
      ) |>
      ungroup() |>
      mutate(across(
        .cols = starts_with("b_Land_use"),
        .fns = function(x) x + r_Dataset_id
      ))
  } else {
    spread_draws <- model |>
      tidybayes::spread_draws(
        b_Land_useNatural_vegetation,
        b_Land_useAgriculture,
        b_Land_useUrban,
        b_Land_useForestry,
        ndraws = 1000,
        seed = 111
      )
  }

  Comparisons <- spread_draws |>
    mutate(
      AN = b_Land_useAgriculture - b_Land_useNatural_vegetation,
      UN = b_Land_useUrban - b_Land_useNatural_vegetation,
      FN = b_Land_useForestry - b_Land_useNatural_vegetation
    ) |>
    select(any_of("Dataset_id"), AN, UN, FN) |>
    pivot_longer(
      cols = all_of(c("AN", "UN", "FN")),
      names_to = "Comparison",
      values_to = "Ratio"
    ) |>
    mutate(
      Comparison = factor(
        Comparison,
        levels = c("AN", "UN", "FN"),
        labels = c(
          "Agriculture/Natural vegetation",
          "Urban/Natural vegetation",
          "Forestry/Natural vegetation"
        )
      )
    )

  return(list(Comparisons = Comparisons, metric = metric))
}


#' Prepare a tidy comparison of posteriors data frame
#'
#' This helper takes the two raw data frames that contain a column named
#' `Ratio` (one for the *0*‑scenario and one for the *2*‑scenario) and
#' reshapes them into a single long‑format tibble ready for plotting.
#'
#' @param df0 A data frame (or tibble) containing the *0*‑scenario results.
#'   It must have at least the columns `Comparison` and `Ratio`.
#' @param df2 A data frame (or tibble) containing the *2*‑scenario results.
#'   It must have at least the columns `Comparison` and `Ratio`.
#' @param scale A character string that identifies the metric (e.g. `"alpha"`,
#'   `"beta"` or `"gamma"`).  The scale is appended to the new column names
#'   (`Ratio<scale>0` and `Ratio<scale>2`) so they stay distinct after the
#'   bind‑cols operation.
#'
#' @return A tibble with three columns:
#'   \describe{
#'     \item{Comparison}{Factor with ordered levels
#'       `"Forestry/Natural vegetation"`,
#'       `"Agriculture/Natural vegetation"`,
#'       `"Urban/Natural vegetation"` and the corresponding short
#'       labels `"Forestry"`, `"Agriculture"`, `"Urban"` .}
#'     \item{ratio_type}{Character indicating which ratio the row refers to
#'       (`"Ratio<suffix>0"` or `"Ratio<suffix>2"`).}
#'     \item{Ratio}{Numeric effect‑size value.}
#'   }
#'   The result is suitable for direct use with `ggplot2` (e.g. via
#'   `geom_density_ridges_gradient`).
#'   No side‑effects are produced.
#' @examples
#'   df_long <- prep_comparisons(ComparisonsA0, ComparisonsA2, "alpha")

.wrangle_posterior <- function(df0, df2, suffix) {
  # rename the Ratio column so we can keep both versions
  bind_rows(df0, df2, .id = "ratio_type") |>
    mutate(
      ratio_type = stringi::stri_replace_all_fixed(
        str = ratio_type,
        pattern = c("1", "2"),
        replacement = c(
          paste0("ratio", suffix, "0"),
          paste0("ratio", suffix, "2")
        ),
        vectorize_all = FALSE
      ),
      Comparison = factor(
        Comparison,
        levels = c(
          "Forestry/Natural vegetation",
          "Agriculture/Natural vegetation",
          "Urban/Natural vegetation"
        ),
        labels = c("Forestry", "Agriculture", "Urban")
      )
    )
}
