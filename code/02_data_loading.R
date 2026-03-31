# =============================================================================
# 02_data_loading.R
# Input data specification and validation
# =============================================================================
#
# PURPOSE
# -------
# The framework does NOT load or pre-process raw data files.  Data preparation
# (loading, translating, and initial cleaning) is handled in distinct scripts
# (e.g. 01a_lichen_prep.R, 01b_plot_prep.R).  Modify those scripts for
# every new dataset without touching this framework.
#
# This file provides two utilities:
#   1. validate_input_data()  – confirms the data frame has the right structure
#   2. describe_data()        – prints a column-by-column overview for inspection
#
# ─────────────────────── REQUIRED INPUT DATA STRUCTURE ───────────────────────
#
# Pass a SINGLE "modelling data frame" (one row per plot) that contains:
#
#   IDENTIFIER
#     <id_col>         character/integer – unique plot ID
#                      (set id_col in utils.R → get_project_config()$data$id_col)
#
#   COORDINATES  (needed for spatial diagnostics in 07_model_diagnostics.R)
#     X                numeric – easting
#     Y                numeric – northing
#
#   PREDICTORS   (numeric; scale them first with scale_predictors() in 03)
#     ...any columns you pass to predictors = c(...) in fit_glm_model()
#
#   CATEGORICAL PREDICTORS (optional, processed by 04_categorical_encoding.R)
#     management, exposure, dominant_species, etc.
#     Column names and category labels should match utils.R → categorical$*
#
#   RESPONSE VARIABLES
#     binary cols      integer 0/1  – for binomial GLMs / GLMMs
#     count  cols      integer >= 0 – for NB / Poisson GLMs / GLMMs
#     Column names should match utils.R → lichen_groups$binary / $count
#
# ──────────────────────────────────────────────────────────────────────────────

library(dplyr)
library(tibble)
library(purrr)

source(here::here("utils.R"))


# -----------------------------------------------------------------------------
# PUBLIC FUNCTIONS
# -----------------------------------------------------------------------------

#' Validate the structure of a user-prepared modelling data frame
#'
#' Call this after your own data-preparation scripts to confirm that the
#' resulting data frame has the structure expected by the modelling framework
#' (R/03 onwards).  The function prints a formatted report and returns results
#' invisibly so it can be used inside a pipeline.
#'
#' @param data A data frame or tibble prepared by the user (pipe-friendly).
#' @param required_cols Character vector of column names that MUST be present.
#'   Default: empty (no columns required beyond the id column).
#' @param id_col Character. Name of the plot identifier column. Default from
#'   \code{get_project_config()$data$id_col}.
#' @param check_coords Logical. Warn if X or Y coordinate columns are missing.
#'   Default TRUE (coordinates are needed for spatial diagnostics).
#' @return A named list (returned invisibly): \code{valid} (logical),
#'   \code{missing_cols} (character), \code{n_rows}, \code{n_cols},
#'   \code{messages} (character vector of issues found).
#' @examples
#' # After your own prep scripts produce modeling_data:
#' v <- validate_input_data(
#'   modeling_data,
#'   required_cols = c("X", "Y", "parmelia_agg_presence",
#'                     "calicioids_richness", "dbh_max_scaled")
#' )
#' stopifnot(v$valid)
validate_input_data <- function(data,
                                required_cols = character(),
                                id_col        = get_project_config()$data$id_col,
                                check_coords  = TRUE) {
  stopifnot(is.data.frame(data))
  
  lichen_banner("Validating Input Data Structure")
  
  messages     <- character()
  valid        <- TRUE
  missing_cols <- character()
  
  # 1. ID column
  if (!is.null(id_col) && nchar(id_col) > 0 && !id_col %in% colnames(data)) {
    msg <- paste0("ID column '", id_col, "' not found. ",
                  "Set data$id_col in utils.R to match your data.")
    lichen_warning(msg)
    messages <- c(messages, msg)
    valid    <- FALSE
  }
  
  # 2. Required columns
  if (length(required_cols) > 0) {
    missing_cols <- setdiff(required_cols, colnames(data))
    if (length(missing_cols) > 0) {
      msg <- paste("Missing required column(s):",
                   paste(missing_cols, collapse = ", "))
      lichen_warning(msg)
      messages <- c(messages, msg)
      valid    <- FALSE
    }
  }
  
  # 3. Coordinate columns for spatial diagnostics
  if (check_coords) {
    missing_xy <- setdiff(c("X", "Y"), colnames(data))
    if (length(missing_xy) > 0) {
      msg <- paste0("Coordinate column(s) missing: ",
                    paste(missing_xy, collapse = ", "),
                    " – needed by calculate_morans_i() / run_dharma_tests()")
      lichen_warning(msg)
      messages <- c(messages, msg)
    }
  }
  
  # 4. All-NA columns
  all_na <- names(which(sapply(data, function(x) all(is.na(x)))))
  if (length(all_na) > 0) {
    msg <- paste("All-NA column(s) detected:", paste(all_na, collapse = ", "))
    lichen_warning(msg)
    messages <- c(messages, msg)
  }
  
  cat(sprintf("  Rows: %d  |  Columns: %d\n", nrow(data), ncol(data)))
  
  if (valid && length(messages) == 0) {
    lichen_message("Input data is valid and ready for the modelling framework")
  } else if (!valid) {
    cat("\n  Please fix the following before proceeding:\n")
    for (m in messages) cat("  \u2022 ", m, "\n", sep = "")
  }
  
  invisible(list(
    valid        = valid,
    missing_cols = missing_cols,
    n_rows       = nrow(data),
    n_cols       = ncol(data),
    messages     = messages
  ))
}


#' Print a concise structural summary of a data frame
#'
#' Helps you verify that your prepared data has the expected columns, types,
#' and missingness before passing it to the modelling framework.  The summary
#' is always printed; the tibble is also returned invisibly for downstream use.
#'
#' @param data A data frame or tibble (pipe-friendly).
#' @return A tibble (one row per column) with: column, type, n_missing,
#'   pct_missing, min_val, max_val, n_unique.  Returned invisibly.
#' @examples
#' describe_data(modeling_data)
describe_data <- function(data) {
  stopifnot(is.data.frame(data))
  
  lichen_banner("Data Structure Summary")
  cat(sprintf("  %d rows \u00d7 %d columns\n\n", nrow(data), ncol(data)))
  
  result <- purrr::map_dfr(colnames(data), function(col) {
    x      <- data[[col]]
    n_miss <- sum(is.na(x))
    tibble::tibble(
      column      = col,
      type        = class(x)[1],
      n_missing   = n_miss,
      pct_missing = round(n_miss / length(x) * 100, 1),
      min_val     = if (is.numeric(x)) as.character(round(min(x, na.rm = TRUE), 3)) else NA_character_,
      max_val     = if (is.numeric(x)) as.character(round(max(x, na.rm = TRUE), 3)) else NA_character_,
      n_unique    = dplyr::n_distinct(x, na.rm = TRUE)
    )
  })
  
  print(result, n = nrow(result))
  invisible(result)
}