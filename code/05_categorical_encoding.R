# =============================================================================
# 05_categorical_encoding.R
# Generic functions for encoding categorical predictor variables
# =============================================================================
# Provides: encode_ordered_factor(), encode_management(), encode_exposure(),
#           encode_dominant_species(), encode_all_categoricals()
#
# Default category labels are read from utils.R → get_project_config().
# To adapt for a new dataset:
#   1. Update categorical$management_levels and categorical$main_tree_species
#      in utils.R.
#   2. OR pass your own values directly via the function arguments.
# =============================================================================

library(dplyr)
library(stringr)
library(tibble)

source(here::here("utils.R"))


# -----------------------------------------------------------------------------
# GENERIC HELPER
# -----------------------------------------------------------------------------

#' Encode any column as an ordered factor with derived integer and binary columns
#'
#' A generic building block used by `encode_management()` but applicable to any
#' ordinal categorical variable (e.g. disturbance intensity, stand age class).
#'
#' @param data A data frame or tibble (pipe-friendly).
#' @param col Character. Name of the column to encode.
#' @param ordered_levels Character vector of category labels ordered from
#'   lowest to highest.
#' @param binary_threshold_level Character. Level at or below which the binary
#'   column is 0 (e.g. the "none" category). Default: first element of
#'   \code{ordered_levels}.
#' @param suffix Character. Suffix appended to new column names (e.g. "logging"
#'   produces \code{logging_factor}, \code{logged}, \code{logging_intensity}).
#'   Default derived from \code{col}.
#' @return The input tibble with three new columns:
#'   \code{<suffix>_factor} (ordered factor),
#'   \code{<binary_col>} (0/1 integer),
#'   \code{<suffix>_intensity} (integer 0 … n−1).
#' @examples
#' df |> encode_ordered_factor(
#'   col            = "disturbance",
#'   ordered_levels = c("none", "low", "moderate", "high"),
#'   suffix         = "disturbance"
#' )
encode_ordered_factor <- function(data,
                                  col,
                                  ordered_levels,
                                  binary_threshold_level = ordered_levels[1],
                                  suffix = gsub("[^a-z0-9_]", "_",
                                                tolower(col))) {
  stopifnot(is.data.frame(data), is.character(col),
            is.character(ordered_levels), length(ordered_levels) > 1)
  
  if (!col %in% colnames(data)) {
    stop("Column '", col, "' not found in data.")
  }
  
  factor_col    <- paste0(suffix, "_factor")
  intensity_col <- paste0(suffix, "_intensity")
  binary_col    <- paste0(suffix, "_binary")
  
  data <- data |>
    dplyr::mutate(
      !!factor_col    := factor(!!rlang::sym(col),
                                levels  = ordered_levels,
                                ordered = TRUE),
      !!binary_col    := dplyr::if_else(
        !!rlang::sym(col) == binary_threshold_level, 0L, 1L
      ),
      !!intensity_col := as.integer(!!rlang::sym(factor_col)) - 1L
    )
  
  cat_counts <- table(data[[factor_col]], useNA = "ifany")
  lichen_message("Encoded '", col, "' as ordered factor. Categories:")
  print(cat_counts)
  
  tibble::as_tibble(data)
}


# -----------------------------------------------------------------------------
# DOMAIN-SPECIFIC WRAPPERS
# -----------------------------------------------------------------------------

#' Encode the past-management (logging history) column
#'
#' Wraps \code{encode_ordered_factor()} for the management variable and creates
#' three named columns:
#' * \code{management_factor} – ordered factor (levels in \code{ordered_levels})
#' * \code{logged}            – integer 0/1 (0 = unlogged / first level)
#' * \code{logging_intensity} – integer 0 … n−1
#'
#' Default level labels are read from
#' \code{get_project_config()$categorical$management_levels}.  Pass your own
#' vector to \code{ordered_levels} when your data uses different labels.
#'
#' @param data A data frame or tibble (pipe-friendly).
#' @param col Character. Name of the management column. Default "management".
#' @param ordered_levels Character vector of management category labels ordered
#'   from least to most intensive.  Default from project config.
#' @return The input tibble with three new columns.
#' @examples
#' # Using project config defaults
#' structure_enc <- structure_clean |> encode_management()
#'
#' # Custom labels for a different dataset
#' structure_enc <- structure_clean |>
#'   encode_management(ordered_levels = c("unmanaged", "selective", "clearcut"))
encode_management <- function(
    data,
    col            = "management",
    ordered_levels = get_project_config()$categorical$management_levels) {
  
  stopifnot(is.data.frame(data))
  
  if (is.null(ordered_levels) || length(ordered_levels) < 2) {
    stop("encode_management() requires 'ordered_levels' with at least 2 elements.\n",
         "Update categorical$management_levels in utils.R or pass ordered_levels directly.")
  }
  
  data <- encode_ordered_factor(
    data,
    col                    = col,
    ordered_levels         = ordered_levels,
    binary_threshold_level = ordered_levels[1],
    suffix                 = "management"
  ) |>
    # Rename generic binary/intensity columns to domain-specific names
    dplyr::rename(
      logged            = management_binary,
      logging_intensity = management_intensity
    )
  
  tibble::as_tibble(data)
}


#' Encode an aspect / exposure column into directional groups
#'
#' Simplifies cardinal and intercardinal direction labels into groups and
#' creates:
#' * \code{exposure_simple} – simplified character label
#' * \code{exposure_factor} – unordered factor
#'
#' By default, groups are derived by matching common direction strings (N, S,
#' E, W, NE, SE, NW, SW, flat/evaluated).  Pass \code{grouping_map} to
#' override with a custom named vector for your own labels.
#'
#' @param data A data frame or tibble (pipe-friendly).
#' @param col Character. Name of the exposure column. Default "exposure".
#' @param grouping_map Named character vector where **names** are the output
#'   group labels and **values** are regex patterns matched against the column.
#'   Rows not matching any pattern receive the label "other".
#'   Default: cardinal direction patterns (N/S/E/W/NE_SE/NW_SW/flat).
#' @return The input tibble with two new columns appended.
#' @examples
#' # Default (cardinal direction grouping)
#' structure_enc <- structure_enc |> encode_exposure()
#'
#' # Custom grouping for a dataset with three slope classes
#' structure_enc <- structure_enc |>
#'   encode_exposure(
#'     grouping_map = c(north_facing = "north|NW|NE",
#'                      south_facing = "south|SW|SE",
#'                      flat         = "flat|level")
#'   )
encode_exposure <- function(data,
                            col          = "exposure",
                            grouping_map = NULL) {
  stopifnot(is.data.frame(data))
  
  if (!col %in% colnames(data)) {
    stop("Column '", col, "' not found in data.")
  }
  
  # Default grouping uses cardinal direction patterns
  if (is.null(grouping_map)) {
    grouping_map <- c(
      flat   = "not evaluated|flat",
      N      = "^N ",
      S      = "^S ",
      E      = "^E ",
      W      = "^W ",
      NE_SE  = "NE|SE",
      NW_SW  = "NW|SW"
    )
  }
  
  # Build the case_when call dynamically from the named vector
  col_sym <- rlang::sym(col)
  case_exprs <- purrr::imap(grouping_map, function(pattern, label) {
    rlang::expr(stringr::str_detect(!!col_sym, !!pattern) ~ !!label)
  })
  
  data <- data |>
    dplyr::mutate(
      exposure_simple = dplyr::case_when(
        !!!case_exprs,
        TRUE ~ "other"
      ),
      exposure_factor = factor(exposure_simple)
    )
  
  cat_counts <- table(data$exposure_simple, useNA = "ifany")
  lichen_message("Exposure encoded. Categories:")
  print(cat_counts)
  
  tibble::as_tibble(data)
}


#' Encode the dominant tree species column
#'
#' Groups rare species codes into "other" and creates:
#' * \code{dominant_grouped} – character with main species or "other"
#' * \code{dominant_factor}  – factor with main species + "other" levels
#'
#' Default main species are read from
#' \code{get_project_config()$categorical$main_tree_species}.
#'
#' @param data A data frame or tibble (pipe-friendly).
#' @param col Character. Name of the dominant species column.
#'   Default "dominant_species".
#' @param main_species Character vector of species codes to keep as individual
#'   levels.  Default from project config
#'   (\code{get_project_config()$categorical$main_tree_species}).
#' @return The input tibble with two new columns appended.
#' @examples
#' # Project-config defaults (SM, BK, JD)
#' structure_enc <- structure_enc |> encode_dominant_species()
#'
#' # Custom species codes
#' structure_enc <- structure_enc |>
#'   encode_dominant_species(main_species = c("Picea", "Fagus", "Abies"))
encode_dominant_species <- function(
    data,
    col          = "dominant_species",
    main_species = get_project_config()$categorical$main_tree_species) {
  
  stopifnot(is.data.frame(data))
  
  if (!col %in% colnames(data)) {
    stop("Column '", col, "' not found in data.")
  }
  
  if (is.null(main_species) || length(main_species) == 0) {
    stop("encode_dominant_species() requires 'main_species'.\n",
         "Update categorical$main_tree_species in utils.R or pass main_species directly.")
  }
  
  factor_levels <- c(main_species, "other")
  
  data <- data |>
    dplyr::mutate(
      dominant_grouped = dplyr::if_else(
        !!rlang::sym(col) %in% main_species,
        !!rlang::sym(col),
        "other"
      ),
      dominant_factor = factor(dominant_grouped, levels = factor_levels)
    )
  
  cat_counts <- table(data$dominant_grouped, useNA = "ifany")
  lichen_message("Dominant species encoded. Categories:")
  print(cat_counts)
  
  tibble::as_tibble(data)
}


#' Apply all categorical encoding steps in sequence
#'
#' Convenience wrapper that calls \code{encode_management()},
#' \code{encode_exposure()}, and \code{encode_dominant_species()} only when the
#' respective column is present in the data.
#'
#' @param data A data frame or tibble (pipe-friendly).
#' @param management_col Character. Column name for management. Default "management".
#' @param exposure_col Character. Column name for exposure. Default "exposure".
#' @param dominant_col Character. Column name for dominant species.
#'   Default "dominant_species".
#' @param ... Additional arguments forwarded to the individual encoding
#'   functions (e.g. \code{ordered_levels}, \code{grouping_map},
#'   \code{main_species}).
#' @return The input tibble with all categorical encoding columns appended.
#' @examples
#' structure_enc <- structure_clean |> encode_all_categoricals()
encode_all_categoricals <- function(data,
                                    management_col = "management",
                                    exposure_col   = "exposure",
                                    dominant_col   = "dominant_species",
                                    ...) {
  stopifnot(is.data.frame(data))
  
  lichen_banner("Encoding Categorical Variables")
  
  if (management_col %in% colnames(data)) {
    data <- encode_management(data, col = management_col, ...)
  } else {
    lichen_warning("Column '", management_col, "' not found – skipping management encoding")
  }
  
  if (exposure_col %in% colnames(data)) {
    data <- encode_exposure(data, col = exposure_col, ...)
  } else {
    lichen_warning("Column '", exposure_col, "' not found – skipping exposure encoding")
  }
  
  if (dominant_col %in% colnames(data)) {
    data <- encode_dominant_species(data, col = dominant_col, ...)
  } else {
    lichen_warning("Column '", dominant_col, "' not found – skipping dominant species encoding")
  }
  
  lichen_message("All categorical encoding complete")
  tibble::as_tibble(data)
}