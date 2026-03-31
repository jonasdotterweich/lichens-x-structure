# =============================================================================
# 04_predictor_engineering.R
# Functions for feature engineering and collinearity checking
# =============================================================================
# Provides: create_deadwood_composites(), select_core_predictors(),
#           calculate_correlation_matrix(), calculate_vif(),
#           identify_collinear_pairs()
# =============================================================================

library(dplyr)
library(tidyr)
library(tibble)

source(here::here("utils.R"))


# -----------------------------------------------------------------------------
# PUBLIC FUNCTIONS
# -----------------------------------------------------------------------------

#' Create composite deadwood variables
#'
#' Adds three derived columns to the data:
#' * `deadwood_decay2to5` – sum of decay stages 2-5
#' * `deadwood_decay4to5` – sum of decay stages 4-5 (highly decomposed)
#' * `pct_decay_advanced` – percentage of total deadwood in stages 4-5
#'
#' @param data A data frame or tibble containing individual decay-stage columns
#'   (pipe-friendly). Expected columns: `deadwood_decay2`, `deadwood_decay3`,
#'   `deadwood_decay4`, `deadwood_decay5`, `deadwood_total`.
#' @return The input tibble with three new columns appended.
#' @examples
#' structure_fe <- structure_clean |> create_deadwood_composites()
create_deadwood_composites <- function(data) {
  stopifnot(is.data.frame(data))
  
  required <- c("deadwood_decay2", "deadwood_decay3",
                "deadwood_decay4", "deadwood_decay5", "deadwood_total")
  missing_cols <- setdiff(required, colnames(data))
  if (length(missing_cols) > 0) {
    stop("create_deadwood_composites() requires these columns: ",
         paste(missing_cols, collapse = ", "))
  }
  
  data <- data |>
    dplyr::mutate(
      deadwood_decay2to5 = deadwood_decay2 + deadwood_decay3 +
        deadwood_decay4 + deadwood_decay5,
      deadwood_decay4to5 = deadwood_decay4 + deadwood_decay5,
      pct_decay_advanced = dplyr::if_else(
        deadwood_total > 0,
        (deadwood_decay4to5 / deadwood_total) * 100,
        0
      )
    )
  
  lichen_message("Created composite deadwood variables: ",
                 "deadwood_decay2to5, deadwood_decay4to5, pct_decay_advanced")
  tibble::as_tibble(data)
}


#' Select the core predictor set for modeling
#'
#' Retains the project ID and a curated set of predictors chosen on the basis
#' of ecological relevance and low collinearity (as established during
#' exploratory analysis). Scaled versions of continuous predictors are also
#' retained when present.
#'
#' @param data A data frame or tibble (pipe-friendly).
#' @param id_col Character. Name of the plot ID column. Default "project_id".
#' @param extra_cols Character vector of additional columns to keep (e.g.
#'   categorical encoded columns). Default `NULL`.
#' @return A tibble with the core predictor columns.
#' @examples
#' structure_core <- structure_fe |> select_core_predictors()
select_core_predictors <- function(data,
                                   id_col     = "project_id",
                                   extra_cols = NULL) {
  stopifnot(is.data.frame(data))
  
  core <- c(
    id_col,
    # Old-growth proxies
    "dbh_max", "n_living_trees_80cm", "n_dead_trees_50cm",
    # Deadwood
    "deadwood_decay4to5", "deadwood_total",
    # Canopy
    "canopy_cover",
    # Tree composition
    "ba_spruce", "ba_beech", "ba_late_successional",
    # Stand structure
    "tree_height_median", "elevation",
    # Management
    "logging_intensity", "logged"
  )
  
  all_keep <- unique(c(core, extra_cols))
  keep     <- intersect(all_keep, colnames(data))
  
  # Also keep any _scaled variants if they exist
  scaled_keep <- grep("_scaled$", colnames(data), value = TRUE)
  keep <- unique(c(keep, scaled_keep))
  
  data <- data |> dplyr::select(dplyr::all_of(keep))
  
  lichen_message("Core predictor set selected: ", ncol(data) - 1,
                 " predictors (+ ID)")
  tibble::as_tibble(data)
}


#' Scale continuous predictors (mean = 0, SD = 1)
#'
#' Adds `_scaled` suffix columns for every numeric column that does **not**
#' already end in `_scaled`, is not the ID column, and is not a binary
#' 0/1 column.
#'
#' @param data A data frame or tibble (pipe-friendly).
#' @param id_col Character. Name of the plot ID column to exclude from scaling.
#'   Default "project_id".
#' @return The input tibble with additional `*_scaled` columns appended.
#' @examples
#' structure_scaled <- structure_core |> scale_predictors()
scale_predictors <- function(data, id_col = "project_id") {
  stopifnot(is.data.frame(data))
  
  num_cols <- data |>
    dplyr::select(-dplyr::any_of(id_col)) |>
    dplyr::select(dplyr::where(is.numeric)) |>
    colnames()
  
  # Exclude already-scaled and binary columns
  to_scale <- num_cols[
    !grepl("_scaled$", num_cols) &
      sapply(num_cols, function(col) {
        vals <- unique(stats::na.omit(data[[col]]))
        length(vals) > 2 || !all(vals %in% c(0, 1))
      })
  ]
  
  data <- data |>
    dplyr::mutate(dplyr::across(
      dplyr::all_of(to_scale),
      list(scaled = ~as.numeric(scale(.))),
      .names = "{.col}_scaled"
    ))
  
  lichen_message("Scaled ", length(to_scale), " continuous predictor(s)")
  tibble::as_tibble(data)
}


#' Calculate the correlation matrix for numeric predictors
#'
#' @param data A data frame or tibble (pipe-friendly). Non-numeric columns are
#'   silently dropped before computation.
#' @param use Character. Passed to `cor()`. Default "complete.obs".
#' @param method Character. Correlation method. Default "pearson".
#' @return A tibble in long format: var1, var2, correlation.
#' @examples
#' cor_long <- structure_clean |>
#'   dplyr::select(where(is.numeric)) |>
#'   calculate_correlation_matrix()
calculate_correlation_matrix <- function(data,
                                         use    = "complete.obs",
                                         method = "pearson") {
  stopifnot(is.data.frame(data))
  
  num_data <- data |> dplyr::select(dplyr::where(is.numeric)) |> tidyr::drop_na()
  
  if (ncol(num_data) < 2) {
    stop("Need at least 2 numeric columns to compute correlations.")
  }
  
  cor_mat <- stats::cor(num_data, use = use, method = method)
  
  # Convert to long tibble (upper triangle only to avoid duplicates)
  cor_long <- cor_mat |>
    as.data.frame() |>
    tibble::rownames_to_column("var1") |>
    tidyr::pivot_longer(-var1, names_to = "var2", values_to = "correlation") |>
    dplyr::filter(var1 < var2) |>
    dplyr::mutate(correlation = round(correlation, 4)) |>
    dplyr::arrange(dplyr::desc(abs(correlation)))
  
  lichen_message("Correlation matrix computed for ", ncol(num_data),
                 " variables (", nrow(cor_long), " pairs)")
  cor_long
}


#' Calculate Variance Inflation Factors (VIF) for a set of predictors
#'
#' Uses `car::vif()` on an OLS regression of each predictor against all others.
#'
#' @param data A data frame or tibble (pipe-friendly).
#' @param predictors Character vector of column names to include in the VIF
#'   calculation. If `NULL` (default), all numeric columns are used.
#' @param threshold Numeric. VIF threshold above which a variable is flagged.
#'   Default 10.
#' @return A tibble with columns: variable, VIF, status.
#' @examples
#' vif_result <- structure_clean |>
#'   calculate_vif(predictors = c("dbh_max", "canopy_cover", "elevation"),
#'                 threshold = 10)
calculate_vif <- function(data, predictors = NULL, threshold = 10) {
  if (!requireNamespace("car", quietly = TRUE)) {
    stop("Package 'car' is required for VIF calculation. ",
         "Install it with: install.packages('car')")
  }
  stopifnot(is.data.frame(data))
  
  if (is.null(predictors)) {
    predictors <- data |> dplyr::select(dplyr::where(is.numeric)) |> colnames()
  }
  
  available <- intersect(predictors, colnames(data))
  if (length(available) < 2) {
    stop("VIF requires at least 2 numeric predictor columns.")
  }
  
  vif_data <- data |> dplyr::select(dplyr::all_of(available)) |> tidyr::drop_na()
  
  # Arbitrary response (VIF is response-independent)
  response_col <- available[1]
  predictors_rhs <- setdiff(available, response_col)
  
  formula_vif <- stats::as.formula(
    paste(response_col, "~", paste(predictors_rhs, collapse = " + "))
  )
  
  lm_vif <- stats::lm(formula_vif, data = vif_data)
  vif_vals <- car::vif(lm_vif)
  
  result <- tibble::tibble(
    variable = names(vif_vals),
    VIF      = round(unname(vif_vals), 2),
    status   = dplyr::case_when(
      VIF < 5           ~ "OK",
      VIF < threshold   ~ "moderate",
      TRUE              ~ "high collinearity"
    )
  ) |>
    dplyr::arrange(dplyr::desc(VIF))
  
  n_high <- sum(result$VIF >= threshold)
  if (n_high > 0) {
    lichen_warning(n_high, " variable(s) with VIF >= ", threshold)
  } else {
    lichen_message("All VIF values below threshold (", threshold, ")")
  }
  
  result
}


#' Identify highly correlated predictor pairs
#'
#' Returns a filtered subset of the correlation matrix where |r| exceeds
#' the threshold.
#'
#' @param data A data frame or tibble (pipe-friendly). Non-numeric columns are
#'   dropped.
#' @param threshold Numeric. |r| threshold. Default 0.7.
#' @param use Character. Passed to `cor()`. Default "complete.obs".
#' @return A tibble with columns: var1, var2, correlation, abs_r.
#' @examples
#' high_cor <- structure_clean |> identify_collinear_pairs(threshold = 0.7)
identify_collinear_pairs <- function(data, threshold = 0.7, use = "complete.obs") {
  stopifnot(is.data.frame(data))
  
  cor_long <- calculate_correlation_matrix(data, use = use)
  
  high_cor <- cor_long |>
    dplyr::mutate(abs_r = abs(correlation)) |>
    dplyr::filter(abs_r >= threshold) |>
    dplyr::arrange(dplyr::desc(abs_r))
  
  if (nrow(high_cor) == 0) {
    lichen_message("No predictor pairs with |r| >= ", threshold)
  } else {
    lichen_warning(nrow(high_cor), " highly correlated pair(s) (|r| >= ",
                   threshold, ")")
    print(high_cor, n = 20)
  }
  
  high_cor
}