# =============================================================================
# 03_data_cleaning.R
# Functions for missing-value handling and outlier detection
# =============================================================================
# Provides: assess_missing_values(), remove_high_missing_plots(),
#           impute_missing_values(), detect_and_handle_outliers()
# =============================================================================

library(dplyr)
library(tidyr)
library(tibble)


source(here::here("utils.R"))


# -----------------------------------------------------------------------------
# PUBLIC FUNCTIONS
# -----------------------------------------------------------------------------

#' Assess missing values in a data frame
#'
#' Computes per-column missing-value counts, percentages, and a plain-text
#' status flag.
#'
#' @param data A data frame or tibble (first argument – pipe-friendly).
#' @return A tibble with columns: variable, n_missing, pct_missing, status.
#' @examples
#' missing_report <- structure_clean |> assess_missing_values()
#' print(missing_report)
assess_missing_values <- function(data) {
  stopifnot(is.data.frame(data))
  
  lichen_banner("Missing Value Assessment")
  
  report <- data |>
    dplyr::summarise(dplyr::across(dplyr::everything(), ~sum(is.na(.)))) |>
    tidyr::pivot_longer(dplyr::everything(),
                        names_to  = "variable",
                        values_to = "n_missing") |>
    dplyr::mutate(
      pct_missing = round(n_missing / nrow(data) * 100, 1),
      status = dplyr::case_when(
        n_missing == 0      ~ "complete",
        pct_missing < 5     ~ "<5% missing",
        pct_missing < 20    ~ "5-20% missing",
        TRUE                ~ ">20% missing"
      )
    ) |>
    dplyr::arrange(dplyr::desc(n_missing))
  
  n_complete <- sum(report$n_missing == 0)
  n_any_na   <- sum(report$n_missing > 0)
  
  lichen_message("Variables: ", nrow(report),
                 " | Complete: ", n_complete,
                 " | With missing: ", n_any_na)
  
  report
}


#' Remove plots that exceed a threshold of missing values in key predictors
#'
#' @param data A data frame or tibble with a plot-ID column (pipe-friendly).
#' @param key_predictors Character vector of column names to check for
#'   missingness. Set to \code{NULL} (default) to skip this step and return
#'   the data unchanged.  Pass the predictor columns relevant to your dataset,
#'   e.g. \code{c("dbh_max", "canopy_cover", "deadwood_total")}.
#' @param max_missing Integer. Maximum number of missing key predictors
#'   allowed before a plot is removed. Default 2.
#' @param id_col Character. Name of the plot ID column.
#'   Default from \code{get_project_config()$data$id_col}.
#' @return The data with high-missing plots removed (tibble).
#' @examples
#' structure_clean <- structure_raw |>
#'   remove_high_missing_plots(
#'     key_predictors = c("dbh_max", "canopy_cover", "deadwood_total"),
#'     max_missing    = 2
#'   )
remove_high_missing_plots <- function(
    data,
    key_predictors = NULL,
    max_missing    = 2L,
    id_col         = get_project_config()$data$id_col) {
  stopifnot(is.data.frame(data))
  
  if (is.null(key_predictors) || length(key_predictors) == 0) {
    lichen_message("key_predictors not specified – skipping plot removal by missing count")
    return(tibble::as_tibble(data))
  }
  
  # Compute once – reused inside the mutate
  present_preds <- intersect(key_predictors, colnames(data))
  
  plots_to_remove <- data |>
    dplyr::mutate(
      .n_na = rowSums(is.na(data[, present_preds, drop = FALSE]))
    ) |>
    dplyr::filter(.n_na > max_missing) |>
    dplyr::pull(dplyr::all_of(id_col))
  
  if (length(plots_to_remove) > 0) {
    lichen_warning(length(plots_to_remove),
                   " plot(s) removed (>", max_missing,
                   " missing key predictors): ",
                   paste(plots_to_remove, collapse = ", "))
    data <- data |>
      dplyr::filter(!(.data[[id_col]] %in% plots_to_remove))
  } else {
    lichen_message("No plots removed – all within missing threshold")
  }
  
  lichen_message("Plots remaining: ", nrow(data))
  tibble::as_tibble(data)
}


#' Impute missing values in numeric columns
#'
#' @param data A data frame or tibble (pipe-friendly).
#' @param strategy Character. Imputation strategy: "median" (default) or "mean".
#' @param cols Tidy-select expression. Columns to impute. Default is all
#'   numeric columns (`where(is.numeric)`).
#' @return The data with NAs replaced by the chosen statistic (tibble).
#' @examples
#' structure_imputed <- structure_clean |> impute_missing_values(strategy = "median")
impute_missing_values <- function(data,
                                  strategy = c("median", "mean"),
                                  cols = dplyr::where(is.numeric)) {
  stopifnot(is.data.frame(data))
  strategy <- match.arg(strategy)
  
  stat_fn <- if (strategy == "median") {
    function(x) median(x, na.rm = TRUE)
  } else {
    function(x) mean(x, na.rm = TRUE)
  }
  
  n_before <- sum(is.na(dplyr::select(data, {{ cols }})))
  
  data <- data |>
    dplyr::mutate(dplyr::across({{ cols }},
                                ~ifelse(is.na(.), stat_fn(.), .)))
  
  n_after <- sum(is.na(dplyr::select(data, {{ cols }})))
  lichen_message("Imputed ", n_before - n_after,
                 " missing value(s) using ", strategy)
  tibble::as_tibble(data)
}

#' Conditionally impute missing values when missingness is below a threshold
#'
#' Implements an explicit imputation policy for scientific transparency:
#'   - If any selected column has pct_missing > threshold_pct (and has NAs),
#'     the function STOPS and prints an offender table to force a manual decision.
#'   - If 0 < pct_missing <= threshold_pct, missing values are imputed using the
#'     specified strategy (median/mean).
#'   - If there are no missing values, the data is returned unchanged.
#'   - In all cases, an audit CSV report is written so it is clear whether the
#'     model was fit on original vs imputed data.
#'
#' IMPORTANT: You should pass only predictor columns in `cols` (do not impute
#' response variables unless you have a strong justification).
#'
#' @param data A data frame or tibble (pipe-friendly).
#' @param cols Tidy-select expression. Columns eligible for imputation.
#'   Default: all numeric columns (`where(is.numeric)`).
#' @param threshold_pct Numeric. Maximum allowed percent missing per column
#'   before stopping for manual decision. Default 5.
#' @param strategy Character. Imputation strategy: "median" (default) or "mean".
#' @param id_col Character. Plot ID column name (used only in messages).
#'   Default from \code{get_project_config()$data$id_col}.
#' @param report_path Character. Path to write the audit CSV report.
#'   Default: \code{here::here("outputs", "Sumava", "reports", "qc_imputation_report.csv")}.
#' @return The data (tibble), possibly imputed. An attribute \code{imputation_audit}
#'   is attached when the function runs.
#' @examples
#' # Impute only if <= 5% missing; stop otherwise
#' df2 <- modeling_data |>
#'   impute_missing_values_if_needed(
#'     cols = dplyr::where(is.numeric) &
#'            !dplyr::ends_with("_presence") &
#'            !dplyr::matches("richness$"),
#'     threshold_pct = 5,
#'     strategy = "median"
#'   )
impute_missing_values_if_needed <- function(
    data,
    cols          = dplyr::where(is.numeric),
    threshold_pct = 5,
    strategy      = c("median", "mean"),
    id_col        = get_project_config()$data$id_col,
    report_path   = here::here("outputs", "Sumava", "reports", "qc_imputation_report.csv")
) {
  stopifnot(is.data.frame(data))
  strategy <- match.arg(strategy)
  
  # Determine which columns are selected
  selected_df <- dplyr::select(data, {{ cols }})
  selected_cols <- colnames(selected_df)
  
  if (length(selected_cols) == 0) {
    lichen_message("No columns selected for imputation policy – skipping")
    return(tibble::as_tibble(data))
  }
  
  # Missingness report for selected columns only
  miss <- data |>
    dplyr::select(dplyr::all_of(selected_cols)) |>
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ sum(is.na(.)))) |>
    tidyr::pivot_longer(dplyr::everything(),
                        names_to  = "variable",
                        values_to = "n_missing") |>
    dplyr::mutate(
      pct_missing   = round(n_missing / nrow(data) * 100, 3),
      threshold_pct = threshold_pct,
      strategy      = strategy,
      action = dplyr::case_when(
        n_missing == 0 ~ "no_missing",
        pct_missing <= threshold_pct ~ paste0("impute_", strategy),
        TRUE ~ "STOP_exceeds_threshold"
      )
    ) |>
    dplyr::arrange(dplyr::desc(pct_missing))
  
  # Always write audit report (best-effort)
  dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(miss, report_path)
  
  offenders <- miss |>
    dplyr::filter(n_missing > 0, pct_missing > threshold_pct)
  
  if (nrow(offenders) > 0) {
    lichen_warning("Missingness exceeds ", threshold_pct,
                   "% for ", nrow(offenders), " column(s). Manual decision required.")
    print(offenders, n = nrow(offenders))
    
    stop(
      "Imputation policy stop: column(s) exceed missingness threshold (",
      threshold_pct, "%). See audit report: ", report_path, "\n",
      "Decide manually: fix upstream data, drop predictors, drop plots, or adjust threshold."
    )
  }
  
  cols_to_impute <- miss |>
    dplyr::filter(n_missing > 0, pct_missing <= threshold_pct) |>
    dplyr::pull(variable)
  
  if (length(cols_to_impute) == 0) {
    lichen_message("No imputation needed (0 missing in selected columns). ",
                   "Audit saved: ", report_path)
    
    out <- tibble::as_tibble(data)
    attr(out, "imputation_audit") <- list(
      threshold_pct   = threshold_pct,
      strategy        = strategy,
      imputed_columns = character(),
      report_path     = report_path
    )
    return(out)
  }
  
  lichen_message("Imputing ", length(cols_to_impute),
                 " column(s) with <= ", threshold_pct, "% missing using ", strategy,
                 ". Audit saved: ", report_path)
  
  out <- impute_missing_values(
    data,
    strategy = strategy,
    cols     = dplyr::all_of(cols_to_impute)
  )
  
  attr(out, "imputation_audit") <- list(
    threshold_pct   = threshold_pct,
    strategy        = strategy,
    imputed_columns = cols_to_impute,
    report_path     = report_path
  )
  
  tibble::as_tibble(out)
}

#' Detect and optionally cap extreme outliers
#'
#' Identifies observations with |Z-score| > threshold. By default the function
#' flags them; with `action = "cap"` it replaces them with the 1st/99th
#' percentile values; with `action = "remove"` it drops the rows.
#'
#' @param data A data frame or tibble (pipe-friendly).
#' @param cols Tidy-select expression. Numeric columns to check.
#'   Default: `where(is.numeric)`.
#' @param z_threshold Numeric. |Z-score| threshold. Default 3.
#' @param action Character. One of "flag" (default), "cap", or "remove".
#' @param id_col Character. Name of the plot ID column used in reports.
#'   Default "project_id".
#' @return A tibble, optionally modified. If `action = "flag"`, a `.outlier`
#'   logical column is added; otherwise the data is modified and returned.
#'   The function always prints a report of affected plots.
#' @examples
#' structure_clean <- structure_imputed |>
#'   detect_and_handle_outliers(z_threshold = 3, action = "flag")
detect_and_handle_outliers <- function(
    data,
    cols         = dplyr::where(is.numeric),
    z_threshold  = 3,
    action       = c("flag", "cap", "remove"),
    id_col       = get_project_config()$data$id_col) {
  stopifnot(is.data.frame(data))
  action <- match.arg(action)
  
  # Compute Z-scores for selected columns
  numeric_cols <- data |> dplyr::select({{ cols }}) |> colnames()
  
  outlier_summary <- data |>
    dplyr::select(dplyr::any_of(id_col), dplyr::all_of(numeric_cols)) |>
    tidyr::pivot_longer(-dplyr::any_of(id_col),
                        names_to  = "variable",
                        values_to = "value") |>
    dplyr::group_by(variable) |>
    dplyr::mutate(
      z_score = (value - mean(value, na.rm = TRUE)) / sd(value, na.rm = TRUE)
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(abs(z_score) > z_threshold)
  
  if (nrow(outlier_summary) == 0) {
    lichen_message("No extreme outliers detected (|Z| > ", z_threshold, ")")
    return(tibble::as_tibble(data))
  }
  
  lichen_warning(nrow(outlier_summary), " extreme outlier(s) detected (|Z| > ",
                 z_threshold, ")")
  print(outlier_summary, n = 30)
  
  if (action == "flag") {
    outlier_ids <- unique(outlier_summary[[id_col]])
    data <- data |>
      dplyr::mutate(.outlier = if (id_col %in% colnames(data)) {
        .data[[id_col]] %in% outlier_ids
      } else {
        rep(FALSE, nrow(data))
      })
    lichen_message("Flagged outlier rows with `.outlier` column")
    
  } else if (action == "cap") {
    for (col in unique(outlier_summary$variable)) {
      lo <- quantile(data[[col]], 0.01, na.rm = TRUE)
      hi <- quantile(data[[col]], 0.99, na.rm = TRUE)
      data[[col]] <- pmin(pmax(data[[col]], lo), hi)
    }
    lichen_message("Capped outliers to 1st/99th percentile")
    
  } else if (action == "remove") {
    if (id_col %in% colnames(data)) {
      remove_ids <- unique(outlier_summary[[id_col]])
      data <- data |>
        dplyr::filter(!(.data[[id_col]] %in% remove_ids))
      lichen_message("Removed ", length(remove_ids), " outlier plot(s)")
    } else {
      lichen_warning("id_col '", id_col, "' not found – cannot remove rows")
    }
  }
  
  tibble::as_tibble(data)
}