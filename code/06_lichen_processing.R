# =============================================================================
# 06_lichen_processing.R
# Functions for preparing lichen response variables
# =============================================================================
# Provides: check_species_prevalence(), create_species_groups(),
#           validate_response_variables()
# =============================================================================

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

source(here::here("utils.R"))


# -----------------------------------------------------------------------------
# PUBLIC FUNCTIONS
# -----------------------------------------------------------------------------

#' Check prevalence and modeling suitability of lichen species
#'
#' Computes per-species occurrence statistics across all plots and categorises
#' each species by its suitability for binomial GLM modeling.
#'
#' @param data A data frame or tibble of raw lichen records (pipe-friendly).
#'   Expected columns: `Plot` (plot identifier) and `species`.
#' @param n_plots Integer. Total number of surveyed plots (denominator for
#'   prevalence). Default 120.
#' @param min_prevalence Numeric. Lower bound (%) for "good" modeling range.
#'   Default 20.
#' @param max_prevalence Numeric. Upper bound (%) for "good" modeling range.
#'   Default 80.
#' @return A tibble with one row per species: species, n_plots, prevalence_pct,
#'   mean_cover, median_cover, max_cover, total_observations, model_suitability.
#' @examples
#' prevalence <- lichens_raw |> check_species_prevalence()
#' prevalence |> dplyr::filter(model_suitability == "good_balanced")
check_species_prevalence <- function(data,
                                     n_plots        = get_project_config()$data$n_plots,
                                     min_prevalence = get_project_config()$modeling$prevalence_min,
                                     max_prevalence = get_project_config()$modeling$prevalence_max) {
  stopifnot(is.data.frame(data))
  
  required <- c("Plot", "species")
  missing_cols <- setdiff(required, colnames(data))
  if (length(missing_cols) > 0) {
    stop("check_species_prevalence() requires columns: ",
         paste(missing_cols, collapse = ", "))
  }
  
  cover_col <- intersect(c("cover", "Cover", "abundance"), colnames(data))[1]
  has_cover <- !is.na(cover_col)
  
  .cover_stat <- function(fn) {
    if (has_cover) fn(.data[[cover_col]], na.rm = TRUE) else NA_real_
  }
  
  summary <- data |>
    dplyr::group_by(species) |>
    dplyr::summarise(
      n_plots            = dplyr::n_distinct(Plot),
      total_observations = dplyr::n(),
      mean_cover         = .cover_stat(mean),
      median_cover       = .cover_stat(median),
      max_cover          = .cover_stat(max),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      prevalence_pct    = round(n_plots / !!n_plots * 100, 1),
      model_suitability = dplyr::case_when(
        prevalence_pct >= min_prevalence & prevalence_pct <= max_prevalence
        ~ "good_balanced",
        prevalence_pct > max_prevalence
        ~ "ubiquitous_low_variance",
        prevalence_pct >= 10 & prevalence_pct < min_prevalence
        ~ "rare_consider_grouping",
        TRUE
        ~ "very_rare_not_suitable"
      )
    ) |>
    dplyr::arrange(dplyr::desc(n_plots))
  
  n_good <- sum(summary$model_suitability == "good_balanced")
  lichen_message("Species assessed: ", nrow(summary),
                 " | Suitable for modeling (", min_prevalence, "-", max_prevalence,
                 "%): ", n_good)
  summary
}


#' Create lichen species groups / response variables for modelling
#'
#' Aggregates raw lichen records into plot-level presence/absence and species-
#' richness response variables.  The grouping is fully user-defined: pass a
#' named list where each element describes one response variable.
#'
#' This replaces the previous hard-coded Šumava species lists, making the
#' function applicable to any lichen dataset.
#'
#' @param data A data frame or tibble of raw lichen records (pipe-friendly).
#'   Required columns: \code{Plot} and \code{species}.
#' @param groups A named list.  Each element must be a list with two items:
#'   \itemize{
#'     \item \code{type}    – \code{"binary"} (presence/absence) or
#'                            \code{"count"} (species richness).
#'     \item \code{species} – character vector of species names to include.
#'   }
#'   The name of the list element becomes the output column name.
#'   See the example below.
#' @param plot_col Character. Name of the plot identifier column in \code{data}.
#'   Default "Plot".
#' @param species_col Character. Name of the species name column in \code{data}.
#'   Default "species".
#' @return A tibble with one row per plot and one column per entry in
#'   \code{groups}, plus the plot identifier column.
#' @examples
#' # Define your own species groups
#' my_groups <- list(
#'   parmelia_agg_presence = list(
#'     type    = "binary",
#'     species = c("Parmelia sulcata", "Parmelia saxatilis")
#'   ),
#'   calicioids_richness = list(
#'     type    = "count",
#'     species = c("Calicium viride", "Chaenotheca brunneola",
#'                 "Chaenotheca ferruginea")
#'   )
#' )
#' lichen_groups <- lichens_raw |> create_species_groups(groups = my_groups)
create_species_groups <- function(data,
                                  groups,
                                  plot_col    = "Plot",
                                  species_col = "species") {
  
  stopifnot(is.data.frame(data))
  
  if (missing(groups) || !is.list(groups) || length(groups) == 0) {
    stop(
      "create_species_groups() requires a non-empty named 'groups' list.\n",
      "Each element must have $type ('binary' or 'count') and $species.\n",
      "Example:\n",
      "  groups <- list(\n",
      "    my_presence = list(type = 'binary', species = c('Sp. A', 'Sp. B')),\n",
      "    my_richness = list(type = 'count',  species = c('Sp. C', 'Sp. D'))\n",
      "  )"
    )
  }
  
  required <- c(plot_col, species_col)
  missing_cols <- setdiff(required, colnames(data))
  if (length(missing_cols) > 0) {
    stop("create_species_groups() requires columns: ",
         paste(missing_cols, collapse = ", "))
  }
  
  # Rename to standard names for internal processing
  if (plot_col != "Plot" || species_col != "species") {
    data <- data |>
      dplyr::rename(Plot = dplyr::all_of(plot_col),
                    species = dplyr::all_of(species_col))
  }
  
  all_plots <- tibble::tibble(Plot = unique(data$Plot))
  
  for (grp_name in names(groups)) {
    grp <- groups[[grp_name]]
    
    if (!all(c("type", "species") %in% names(grp))) {
      stop("Group '", grp_name,
           "' must have both $type and $species elements.")
    }
    if (!grp$type %in% c("binary", "count")) {
      stop("Group '", grp_name,
           "': $type must be 'binary' or 'count', not '", grp$type, "'.")
    }
    
    if (grp$type == "binary") {
      present_plots <- data |>
        dplyr::filter(species %in% grp$species) |>
        dplyr::distinct(Plot) |>
        dplyr::mutate(.present = 1L)
      
      col_data <- all_plots |>
        dplyr::left_join(present_plots, by = "Plot") |>
        dplyr::mutate(!!grp_name := tidyr::replace_na(.present, 0L)) |>
        dplyr::select(Plot, dplyr::all_of(grp_name))
      
    } else {  # count / richness
      col_data <- data |>
        dplyr::filter(species %in% grp$species) |>
        dplyr::group_by(Plot) |>
        dplyr::summarise(!!grp_name := dplyr::n_distinct(species),
                         .groups = "drop")
      col_data <- all_plots |>
        dplyr::left_join(col_data, by = "Plot") |>
        dplyr::mutate(!!grp_name := tidyr::replace_na(!!rlang::sym(grp_name), 0L))
    }
    
    all_plots <- all_plots |>
      dplyr::left_join(col_data, by = "Plot")
  }
  
  lichen_message("Species groups created: ", length(groups),
                 " response variable(s) for ", nrow(all_plots), " plots")
  tibble::as_tibble(all_plots)
}


#' Validate response variables in the modeling dataset
#'
#' Checks that all expected response columns are present, reports
#' prevalence/range, and flags potential issues (e.g. zero-inflation, extreme
#' imbalance).
#'
#' @param data A data frame or tibble that includes response variable columns
#'   (pipe-friendly).
#' @param binary_responses Character vector of binary (0/1) response column
#'   names. Defaults to the six binomial response variables.
#' @param count_responses Character vector of count response column names.
#'   Default `"calicioids_richness"`.
#' @return A tibble with one row per response variable: variable, type,
#'   n_present (or mean for counts), prevalence_pct (binary only), min, max,
#'   status.
#' @examples
#' response_check <- modeling_data |> validate_response_variables()
validate_response_variables <- function(
    data,
    binary_responses = get_project_config()$lichen_groups$binary,
    count_responses  = get_project_config()$lichen_groups$count) {
  
  stopifnot(is.data.frame(data))
  
  lichen_banner("Validating Response Variables")
  
  all_responses <- c(binary_responses, count_responses)
  missing_vars  <- setdiff(all_responses, colnames(data))
  
  if (length(missing_vars) > 0) {
    lichen_warning("Missing response columns: ",
                   paste(missing_vars, collapse = ", "))
  }
  
  present_vars <- intersect(all_responses, colnames(data))
  
  report <- purrr::map_dfr(present_vars, function(var) {
    x    <- data[[var]]
    type <- if (var %in% binary_responses) "binary" else "count"
    
    if (type == "binary") {
      n_present      <- sum(x == 1, na.rm = TRUE)
      prevalence_pct <- round(n_present / length(x) * 100, 1)
      status <- dplyr::case_when(
        prevalence_pct < 10 | prevalence_pct > 90 ~ "imbalanced",
        prevalence_pct < 20 | prevalence_pct > 80 ~ "moderate_imbalance",
        TRUE                                        ~ "ok"
      )
      tibble::tibble(variable = var, type = type,
                     n_present = n_present, prevalence_pct = prevalence_pct,
                     min = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE),
                     status = status)
    } else {
      pct_zero <- mean(x == 0, na.rm = TRUE) * 100
      status   <- if (pct_zero > 50) "zero_inflated" else "ok"
      tibble::tibble(variable = var, type = type,
                     n_present = sum(x > 0, na.rm = TRUE),
                     prevalence_pct = NA_real_,
                     min = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE),
                     status = status)
    }
  })
  
  print(report, n = 20)
  lichen_message("Response variable validation complete")
  report
}