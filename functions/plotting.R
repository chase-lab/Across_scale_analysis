#' Plot responses
#' @param model a list of model data given by .wrangle_responses.
#' @param metric Optional string "N" or "S". Taken from `metric` by default.

.plot_responses <- function(model, metric = NULL) {
  if (is.null(metric)) {
    metric <- model$metric
  }
  gradient_colors <- grDevices::colorRampPalette(c(
    "#274659",
    "#CAAE10",
    "#F2790F",
    "#B93102"
  ))(length(unique(model$fitted$Dataset_id)))

  ggplot() +
    # raw data
    geom_point(
      data = model$fitted,
      aes(
        x = Land_use,
        y = get(metric),
        group = Dataset_id,
        color = Dataset_id
      ),
      size = 0.5,
      position = position_dodge(width = 0.5),
      alpha = 0.3
    ) +
    geom_point(
      data = model$study,
      aes(
        x = Land_use,
        y = exp(Estimate),
        color = Dataset_id,
        group = Dataset_id
      ),
      position = position_dodge(width = 0.5),
      alpha = 0.5,
      size = 2
    ) +
    geom_linerange(
      data = model$study,
      aes(
        x = Land_use,
        ymin = exp(lower),
        ymax = exp(upper),
        color = Dataset_id,
        group = Dataset_id
      ),
      position = position_dodge(width = 0.5),
      alpha = 0.5
    ) +
    geom_linerange(
      data = model$fitted |>
        distinct(Land_use, Estimate, .keep_all = TRUE),
      aes(x = Land_use, ymin = Q2.5, ymax = Q97.5),
      lwd = 2.5,
      position = position_dodge(width = 1)
    ) +
    geom_point(
      data = model$fitted |>
        distinct(Land_use, Estimate),
      aes(x = Land_use, y = Estimate),
      shape = 21,
      colour = "black",
      fill = "white",
      stroke = 2,
      size = 4
    ) +
    scale_x_discrete(labels = ~ sub("_", " ", .x, fixed = TRUE)) +
    scale_fill_manual(values = gradient_colors) +
    labs(
      x = '',
      y = case_match(metric, "N" ~ "Number of individuals", "S" ~ "q=0")
    ) +
    scale_y_continuous(trans = 'log2') +
    scale_colour_manual(values = gradient_colors) +
    theme_minimal() +
    theme(
      legend.position = 'none',
      panel.border = element_rect(color = "gray", fill = NA, linewidth = 0.5),
      panel.grid.minor.y = element_blank(),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 11),
      axis.title.y = element_text(size = 14)
    )
}

#' Plot comparisons
#' @param model a list of model data given by .wrangle_comparisons.
#' @param metric Optional string "N" or "S". Taken from `metric` by default.
#'
.plot_comparisons <- function(model, metric = NULL) {
  if (is.null(metric)) {
    metric <- model$metric
  }
  ggplot() +
    geom_density_ridges_gradient(
      data = model$Comparisons,
      aes(
        x = Ratio,
        y = fct_reorder(Comparison, Ratio, .fun = mean),
        fill = factor(after_stat(x) > 0)
      ),
      scale = 0.9,
      alpha = 0.1,
      linetype = 0
    ) +
    xlim(-0.8, 1) +
    geom_vline(
      data = model$Comparisons,
      aes(xintercept = median(Ratio)),
      size = 0.5,
      alpha = 0.5,
      lty = 2
    ) +
    geom_vline(xintercept = 0, size = 0.5) +
    geom_point(
      data = model$Comparisons,
      aes(x = Ratio, y = Comparison),
      stat = ggstance:::StatSummaryh,
      fun.x = median,
      shape = 21,
      color = "black",
      fill = "white",
      stroke = 2,
      size = 4
    ) +
    labs(
      y = case_match(metric, "N" ~ "Number of individuals", "S" ~ "q=0"),
      x = 'Log Ratio'
    ) +
    geom_text(
      data = model$Comparisons |>
        filter(Ratio < 0) |>
        group_by(Comparison) |>
        summarise(Count = n()) |>
        mutate(percentages = Count / 1000) |>
        ungroup() |>
        distinct(Comparison, percentages, .keep_all = T),
      aes(x = -0.55, y = Comparison, label = paste(percentages)),
      size = 3,
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
      panel.border = element_rect(color = "gray", fill = NA, linewidth = 0.5),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 16),
      axis.title.y = element_blank(),
      axis.text.y = element_blank()
    )
}


#' @param model a data.frame
#' @param suffix "description"A", "B" or "C"
.plot_posterior <- function(
  model,
  suffix,
  lab_x = "",
  lab_y = "",
  show_legend = FALSE
) {
  # ----- colour / fill palettes (same hues for all, just different names) -----
  fill_pal <- c(
    scales::alpha("#F2790F", 0.7),
    scales::alpha("#CAAE10", 0.7)
  )
  colour_pal <- c(
    scales::alpha("#F2790F", 0.9),
    scales::alpha("#CAAE10", 0.9)
  )
  names(fill_pal) <- names(colour_pal) <- c(
    paste0("ratio", suffix, "0"),
    paste0("ratio", suffix, "2")
  )

  ggplot() +
    # density ridges -------------------------------------------------
    ggridges::geom_density_ridges_gradient(
      data = model,
      aes(
        x = Ratio,
        y = Comparison,
        fill = ratio_type
      ),
      scale = 0.9,
      alpha = 0.1,
      linetype = 0
    ) +

    # median points --------------------------------------------------
    geom_point(
      data = model,
      aes(
        x = Ratio,
        y = Comparison,
        colour = ratio_type
      ),
      stat = ggstance:::StatSummaryh,
      fun.x = median,
      shape = 21,
      fill = "white",
      stroke = 2,
      size = 4.5
    ) +

    # text annotations (percentage of negative ratios) -------------
    geom_text(
      data = model |>
        dplyr::filter(Ratio < 0) |>
        summarise(Count = n(), .by = c(Comparison, ratio_type)) |>
        mutate(
          d = Count / 1000,
          percentages = sprintf("%.2f", d),
          # shift the label a little so the two ratio types do not overlap
          y_adj = as.numeric(Comparison) +
            if_else(grepl("0$", ratio_type), 0.45, 0.15)
        ),
      aes(
        x = -0.85,
        y = y_adj,
        colour = ratio_type,
        label = percentages
      ),
      size = 4.5,
      parse = TRUE,
      show.legend = FALSE
    ) +

    # vertical line at zero -----------------------------------------
    geom_vline(xintercept = 0, size = 0.5, linetype = 2) +

    # limits & labels ------------------------------------------------
    xlim(-1.0, 0.5) +
    labs(x = lab_x, y = lab_y) +

    # manual palettes ------------------------------------------------
    scale_fill_manual(
      values = fill_pal,
      name = " ",
      labels = c("ratioA0" = "S", "ratioA2" = "Spie")
    ) +
    scale_colour_manual(
      values = colour_pal,
      name = " ",
      labels = c("ratioA0" = "S", "ratioA2" = "Spie")
    ) +

    # tidy axis ------------------------------------------------------
    scale_y_discrete(labels = scales::wrap_format(9)) +

    # theme -----------------------------------------------------------
    theme_minimal() +
    theme(
      legend.position = if (show_legend) "top" else "none",
      legend.direction = "horizontal",
      legend.key.spacing.x = unit(1, "cm"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.ticks.x = element_line(size = 0.5, colour = "gray"),
      panel.border = element_rect(colour = "gray", fill = NA, linewidth = 0.5),
      axis.title = element_text(size = 16),
      axis.title.y = element_text(
        angle = 0,
        vjust = 0.5,
        margin = margin(r = 15)
      ),
      axis.text = element_text(size = 13),
      legend.text = element_text(size = 16)
    )
}

#' @param studylevel a data.frame
#' @returns a `ggplot` object
.plot_comparison_latitude <- function(
  studylevel,
  lab_x = "",
  lab_y = "",
  subtitle = NULL
) {
  gradient_colors <- grDevices::colorRampPalette(c(
    "#274659",
    "#CAAE10",
    "#F2790F",
    "#B93102"
  ))(length(unique(studylevel$Dataset_id)))

  ggplot() +
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
    labs(x = lab_x, y = lab_y, subtitle = subtitle) +
    ylim(-1.0, 0.5) +
    theme_minimal() +
    theme(
      legend.position = 'none',
      panel.border = element_rect(color = "gray", fill = NA, linewidth = 0.5),
      panel.grid.minor.y = element_blank(),
      axis.title.x = element_blank(),
      plot.subtitle = element_text(hjust = 0.5, size = 14),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 16)
    )
}

#' @param studylevel a data.frame
#' @returns a `ggplot` object
.plot_comparison_taxa <- function(
  studylevel,
  taxamean,
  lab_x = "",
  lab_y = "",
  subtitle = NULL
) {
  gradient_colors <- grDevices::colorRampPalette(c(
    "#274659",
    "#CAAE10",
    "#F2790F",
    "#B93102"
  ))(length(unique(studylevel$Dataset_id)))

  ggplot() +
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
    labs(x = lab_x, y = lab_y, subtitle = subtitle) +
    xlim(-0.75, 0.5) +
    theme_minimal() +
    theme(
      legend.position = 'none',
      panel.border = element_rect(color = "gray", fill = NA, linewidth = 0.5),
      panel.grid.minor.y = element_blank(),
      plot.subtitle = element_text(hjust = 0.5, size = 14),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 16),
      axis.title.x = element_text(size = 14)
    )
}
