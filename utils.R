# =============================================================================
# utils.R
# Helper functions and project configuration
# =============================================================================
# Provides: path helpers (here::here()), configuration loading,
#           message formatting, and generic file I/O utilities.
#
# ─────────────────────── HOW TO USE FOR A NEW DATASET ────────────────────────
#
# utils.R is the single place to change when you move to a new dataset.
# Update the sections below and all other scripts will pick up the changes:
#
#   1. data$id_col           → column name used as the plot identifier
#   2. data$n_plots          → total number of surveyed plots (prevalence denominator)
#   3. output$root           → where model outputs are written
#   4. categorical$*         → management levels / species codes for your system
#   5. lichen_groups$*       → response variable column names in your modelling data
#   6. glmm$random_effects   → add random effects here to use fit_glmm_model()
#
# All R/0x_*.R functions read their defaults from get_project_config(), so
# changing this file is usually enough – no need to edit the other scripts.
# =============================================================================

library(here)
library(tibble)
library(dplyr)
library(readr)

# -----------------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------------

#' Return the project configuration as a named list
#'
#' Edit the values in this function to adapt the framework to a new dataset.
#' All downstream functions read their defaults from this config, so changing
#' values here propagates through the entire workflow.
#'
#' @return A named list with data settings, output paths, modelling thresholds,
#'   categorical encoding rules, response variable names, and GLMM settings.
#' @examples
#' cfg <- get_project_config()
#' cfg$data$id_col
#' cfg$categorical$management_levels
get_project_config <- function() {
  list(
    
    # ── DATA ────────────────────────────────────────────────────────────────
    # Paths are only used by the helper functions in 01_data_loading.R.
    # Your own prep scripts (01a, 01b) can ignore these.
    data = list(
      id_col       = "plot_id",               # plot identifier column name
      n_plots      = 120L,                        # total number of surveyed plots
      structural   = here::here("data", "Sumava", "Structural_data.xlsx"),
      lichen       = here::here("data","Sumava", "Lichen_data.xlsx"),
      coordinates  = here::here("data","Sumava","Biodiversity_120 plot.xlsx")
    ),
    
    # ── OUTPUT DIRECTORIES ──────────────────────────────────────────────────
    output = list(
      root     = here::here("outputs","Sumava", "model_outputs"),
      standard = here::here("outputs", "Sumava","model_outputs", "standard"),
      reduced  = here::here("outputs", "Sumava", "model_outputs", "reduced")
    ),
    
    # ── MODELLING THRESHOLDS ────────────────────────────────────────────────
    modeling = list(
      prevalence_min    = 20,    # minimum species prevalence (%) for modeling
      prevalence_max    = 80,    # maximum species prevalence (%) for modeling
      missing_threshold = 2L,    # max missing key predictors before removing a plot
      outlier_z         = 3,     # |Z-score| threshold to flag extreme outliers
      vif_threshold     = 10,    # VIF threshold for collinearity
      cor_threshold     = 0.7,   # |r| threshold for high correlation
      n_sim_dharma      = 1000L  # DHARMa simulation iterations
    ),
    
    # ── CATEGORICAL ENCODING ────────────────────────────────────────────────
    # These values drive encode_management(), encode_dominant_species(), etc.
    # Update them to match the category labels in YOUR cleaned data.
    categorical = list(
      # Management levels ordered from least → most intensive logging.
      # The first element is treated as "unlogged" (logged = 0).
      management_levels = c(
        "no logging, no stumps",
        "logging with most biomass left in place",
        "logging with minority of biomass left in place",
        "logging without leaving biomass"
      ),
      # Dominant tree species codes to keep as individual factor levels.
      # All other codes are collapsed to "other".
      main_tree_species = c("SM", "BK", "JD")
    ),
    
    # ── RESPONSE VARIABLES ──────────────────────────────────────────────────
    # Column names in the modelling data frame.
    # Update to match your own response variables.
    lichen_groups = list(
      binary = c(
        "parmelia_agg_presence",
        "ochrolechia_presence",
        "core_ogf_presence",
        "mycoblastus_presence",
        "xylographa_presence",
        "elite_rare_presence"
      ),
      count  = "calicioids_richness"
    ),
    
    # ── GLMM SETTINGS ───────────────────────────────────────────────────────
    # Used by fit_glmm_model() in 06_model_fitting.R.
    # Add random effects here when you want to account for clustering/grouping.
    #
    # glmmTMB syntax examples:
    #   "(1|region)"           – random intercept per geographic region
    #   "(1|plot_cluster)"     – random intercept per plot cluster
    #   "(1|observer)"         – random intercept per observer/field crew
    #   "(dbh_max|region)"     – random slope + intercept per region
    #
    # Leave as character() for standard GLMs (no random effects).
    glmm = list(
      random_effects = character()
    )
  )
}


# -----------------------------------------------------------------------------
# PATH HELPERS
# -----------------------------------------------------------------------------

#' Build a project-relative output path using here::here()
#'
#' @param ... Path components passed to here::here().
#' @return A character string with the full path.
#' @examples
#' project_path("Lichens", "model_outputs", "results.csv")
project_path <- function(...) {
  here::here(...)
}


#' Ensure an output directory exists, creating it recursively if needed
#'
#' @param path Character. Directory path to create.
#' @return The path invisibly.
#' @examples
#' ensure_dir(project_path("Lichens", "model_outputs", "reduced"))
ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    lichen_message("Created directory: ", path)
  }
  invisible(path)
}


# -----------------------------------------------------------------------------
# MESSAGING / FORMATTING
# -----------------------------------------------------------------------------

#' Print a formatted project message with a check-mark prefix
#'
#' @param ... Character strings passed to cat().
#' @return NULL invisibly.
#' @examples
#' lichen_message("Data loaded: 120 plots")
lichen_message <- function(...) {
  cat("\u2713 ", ..., "\n", sep = "")
  invisible(NULL)
}


#' Print a formatted warning with a flag prefix
#'
#' @param ... Character strings passed to cat().
#' @return NULL invisibly.
#' @examples
#' lichen_warning("3 plots removed due to missing data")
lichen_warning <- function(...) {
  cat("\u26A0\uFE0F  WARNING: ", ..., "\n", sep = "")
  invisible(NULL)
}


#' Print a section header banner
#'
#' @param title Character. Section title.
#' @return NULL invisibly.
#' @examples
#' lichen_banner("STEP 1: Data Loading")
lichen_banner <- function(title) {
  line <- paste(rep("=", 60), collapse = "")
  cat("\n", line, "\n", title, "\n", line, "\n\n", sep = "")
  invisible(NULL)
}


# -----------------------------------------------------------------------------
# FILE I/O
# -----------------------------------------------------------------------------

#' Save a data frame to CSV with a standard message
#'
#' @param data A data frame or tibble.
#' @param path Character. Full file path (use project_path()).
#' @param ... Additional arguments passed to readr::write_csv().
#' @return The data invisibly.
#' @examples
#' save_csv(results, project_path("Lichens", "model_outputs", "results.csv"))
save_csv <- function(data, path, ...) {
  ensure_dir(dirname(path))
  readr::write_csv(data, path, ...)
  lichen_message("Saved CSV: ", basename(path))
  invisible(data)
}


#' Save an R object to an .RData file with a standard message
#'
#' @param ... Objects to save (passed by name to save()).
#' @param path Character. Full .RData file path.
#' @return NULL invisibly.
#' @examples
#' save_rdata(models_glm, diagnostics_summary, path = project_path("Lichens", "models.RData"))
save_rdata <- function(..., path) {
  ensure_dir(dirname(path))
  save(..., file = path)
  lichen_message("Saved RData: ", basename(path))
  invisible(NULL)
}


#' Save a ggplot to PNG with sensible defaults
#'
#' @param plot A ggplot object.
#' @param path Character. Full PNG file path.
#' @param width Numeric. Width in inches (default 10).
#' @param height Numeric. Height in inches (default 7).
#' @param dpi Numeric. Resolution (default 300).
#' @return NULL invisibly.
#' @examples
#' save_plot(my_plot, project_path("Lichens", "model_outputs", "figure.png"))
save_plot <- function(plot, path, width = 10, height = 7, dpi = 300) {
  ensure_dir(dirname(path))
  ggplot2::ggsave(filename = path, plot = plot,
                  width = width, height = height, dpi = dpi)
  lichen_message("Saved plot: ", basename(path))
  invisible(NULL)
}


# -----------------------------------------------------------------------------
# MISC HELPERS
# -----------------------------------------------------------------------------

#' Significance star string for a p-value
#'
#' @param p Numeric p-value.
#' @return Character: "***", "**", "*", "." or "ns".
#' @examples
#' sig_stars(0.001)
#' sig_stars(0.06)
sig_stars <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    p < 0.1   ~ ".",
    TRUE       ~ "ns"
  )
}