# =============================================================================
# 08_model_diagnostics.R
# Functions for spatial autocorrelation and DHARMa diagnostic testing
# =============================================================================
# Provides: calculate_morans_i(), run_dharma_tests(),
#           test_spatial_pattern_comparison(), compile_diagnostic_report()
# =============================================================================

library(dplyr)
library(tibble)
library(ape)

source(here::here("utils.R"))


# -----------------------------------------------------------------------------
# INTERNAL HELPERS
# -----------------------------------------------------------------------------

#' Build an inverse-distance weight matrix from coordinate columns
#'
#' @param coords A data frame or matrix with columns X and Y.
#' @return A square numeric matrix with zeros on the diagonal.
.build_inv_dist_weights <- function(coords) {
  coord_mat  <- as.matrix(coords[, c("X", "Y")])
  dist_mat   <- as.matrix(stats::dist(coord_mat))
  inv_w      <- 1 / dist_mat
  diag(inv_w) <- 0
  inv_w
}


# -----------------------------------------------------------------------------
# PUBLIC FUNCTIONS
# -----------------------------------------------------------------------------

#' Calculate Moran's I for one or more response variables
#'
#' Uses `ape::Moran.I()` with an inverse-distance weight matrix.
#'
#' @param data A data frame or tibble containing the response variable(s) and
#'   coordinate columns X and Y (pipe-friendly).
#' @param variables Character vector of column names to test.
#' @param coords_cols Character vector of length 2. Names of the easting and
#'   northing columns. Default `c("X", "Y")`.
#' @return A tibble with one row per variable: variable, morans_i, expected_i,
#'   sd, p_value, sig.
#' @examples
#' raw_spatial <- data |>
#'   calculate_morans_i(variables = c("parmelia_agg_presence",
#'                                    "calicioids_richness"))
calculate_morans_i <- function(data,
                               variables,
                               coords_cols = c("X", "Y")) {
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package 'ape' is required. Install with: install.packages('ape')")
  }
  stopifnot(is.data.frame(data), is.character(variables))
  
  missing_vars  <- setdiff(variables, colnames(data))
  missing_coord <- setdiff(coords_cols, colnames(data))
  if (length(missing_vars) > 0) {
    stop("Variable(s) not found in data: ", paste(missing_vars, collapse = ", "))
  }
  if (length(missing_coord) > 0) {
    stop("Coordinate column(s) not found: ", paste(missing_coord, collapse = ", "))
  }
  
  lichen_banner("Moran's I – Spatial Autocorrelation")
  
  weights <- .build_inv_dist_weights(data[, coords_cols])
  
  results <- purrr::map_dfr(variables, function(var) {
    moran <- ape::Moran.I(data[[var]], weights)
    tibble::tibble(
      variable  = var,
      morans_i  = round(moran$observed, 4),
      expected_i = round(moran$expected, 4),
      sd        = round(moran$sd, 4),
      p_value   = round(moran$p.value, 4),
      sig       = sig_stars(moran$p.value)
    )
  })
  
  cat(sprintf("%-30s | I = %7.4f | p = %.4f  %s\n",
              results$variable, results$morans_i,
              results$p_value, results$sig))
  
  n_sig <- sum(results$p_value < 0.05)
  lichen_message("Significant spatial clustering (p < 0.05): ",
                 n_sig, "/", nrow(results), " variables")
  results
}


#' Run the full DHARMa diagnostic battery for a fitted model
#'
#' Simulates scaled residuals via `DHARMa::simulateResiduals()` then runs four
#' tests: dispersion, zero-inflation, outliers, and spatial autocorrelation.
#' Optionally saves the DHARMa plot to a PNG file.
#'
#' @param model A fitted `glm` or `glmmTMB` object.
#' @param coords_df A data frame or tibble with columns X and Y, one row per
#'   observation (same order as the modeling data).
#' @param model_name Character. Label used in output messages and the saved
#'   PNG filename. Default: `attr(model, "model_name")`.
#' @param output_dir Character. Directory for the diagnostic PNG. Set to
#'   `NULL` (default) to skip saving.
#' @param n_sim Integer. Number of DHARMa simulations. Default 1000.
#' @return A named list: `diagnostics` (tibble), `sim_resid` (DHARMa object).
#' @examples
#' diag <- run_dharma_tests(m_parmelia,
#'                          coords_df  = data[, c("X", "Y")],
#'                          model_name = "parmelia_agg_presence",
#'                          output_dir = project_path("Lichens", "model_outputs"))
run_dharma_tests <- function(model,
                             coords_df,
                             model_name = NULL,
                             output_dir = NULL,
                             n_sim      = get_project_config()$modeling$n_sim_dharma) {
  if (!requireNamespace("DHARMa", quietly = TRUE)) {
    stop("Package 'DHARMa' is required. Install with: install.packages('DHARMa')")
  }
  if (is.null(model_name)) {
    model_name <- attr(model, "model_name") %||% "model"
  }
  
  lichen_banner(paste("DHARMa Diagnostics:", model_name))
  
  sim_resid <- DHARMa::simulateResiduals(fittedModel = model,
                                         n = n_sim, plot = FALSE)
  
  # Optionally save plot
  if (!is.null(output_dir)) {
    ensure_dir(output_dir)
    plot_path <- file.path(output_dir,
                           paste0("DHARMa_", gsub("\\s+", "_", model_name), ".png"))
    grDevices::png(filename = plot_path, width = 10, height = 6,
                   units = "in", res = 300)
    graphics::plot(sim_resid,
                   main = paste("DHARMa residuals:", model_name))
    grDevices::dev.off()
    lichen_message("DHARMa plot saved: ", basename(plot_path))
  }
  
  # Four core tests
  disp_test    <- DHARMa::testDispersion(sim_resid, plot = FALSE)
  zi_test      <- DHARMa::testZeroInflation(sim_resid, plot = FALSE)
  outlier_test <- DHARMa::testOutliers(sim_resid, plot = FALSE)
  spatial_test <- DHARMa::testSpatialAutocorrelation(
    sim_resid, x = coords_df$X, y = coords_df$Y, plot = FALSE)
  
  cat(sprintf("  Dispersion     p = %.4f %s\n",
              disp_test$p.value, sig_stars(disp_test$p.value)))
  cat(sprintf("  Zero-inflation p = %.4f %s\n",
              zi_test$p.value, sig_stars(zi_test$p.value)))
  cat(sprintf("  Outliers       p = %.4f %s\n",
              outlier_test$p.value, sig_stars(outlier_test$p.value)))
  cat(sprintf("  Spatial autocorr p = %.4f %s\n",
              spatial_test$p.value, sig_stars(spatial_test$p.value)))
  
  diagnostics <- tibble::tibble(
    model              = model_name,
    dispersion_p       = round(disp_test$p.value,    4),
    zero_inflation_p   = round(zi_test$p.value,       4),
    outliers_p         = round(outlier_test$p.value,  4),
    spatial_p          = round(spatial_test$p.value,  4),
    dispersion_ok      = disp_test$p.value    > 0.05,
    zero_inflation_ok  = zi_test$p.value       > 0.05,
    outliers_ok        = outlier_test$p.value  > 0.05,
    spatial_ok         = spatial_test$p.value  > 0.05
  )
  
  all_ok <- all(unlist(diagnostics[, 6:9]))
  if (all_ok) {
    lichen_message("All DHARMa tests PASSED")
  } else {
    lichen_warning("Some DHARMa tests FAILED – see diagnostics tibble")
  }
  
  list(diagnostics = diagnostics, sim_resid = sim_resid)
}


#' Compare spatial autocorrelation before and after model fitting
#'
#' Checks whether the model has explained any spatial clustering that was
#' present in the raw response variable.
#'
#' @param raw_spatial A tibble returned by `calculate_morans_i()` for the
#'   **raw** response variables (before modeling).
#' @param dharma_result A list returned by `run_dharma_tests()` for the same
#'   model (contains `diagnostics`).
#' @param model_name Character. Response variable / model name to look up in
#'   `raw_spatial$variable`.
#' @return A tibble with columns: model, raw_morans_i, raw_p, residual_p,
#'   interpretation.
#' @examples
#' spatial_cmp <- test_spatial_pattern_comparison(
#'   raw_spatial   = raw_spatial_results,
#'   dharma_result = diag_parmelia,
#'   model_name    = "parmelia_agg_presence"
#' )
test_spatial_pattern_comparison <- function(raw_spatial,
                                            dharma_result,
                                            model_name) {
  stopifnot(is.data.frame(raw_spatial))
  
  raw_row <- raw_spatial |>
    dplyr::filter(variable == model_name)
  
  if (nrow(raw_row) == 0) {
    lichen_warning("Model '", model_name,
                   "' not found in raw_spatial – skipping comparison")
    return(NULL)
  }
  
  raw_p      <- raw_row$p_value
  residual_p <- dharma_result$diagnostics$spatial_p
  
  interpretation <- dplyr::case_when(
    raw_p < 0.05 & residual_p >= 0.05 ~ "fully_explained",
    raw_p < 0.05 & residual_p  < 0.05 ~ "partially_explained",
    raw_p >= 0.05                      ~ "no_spatial_issue"
  )
  
  icon <- dplyr::case_when(
    interpretation == "fully_explained"    ~ "\u2705 Fully explained",
    interpretation == "partially_explained" ~ "\u26A0\uFE0F Partially explained",
    TRUE                                    ~ "\u2705 No spatial issues"
  )
  
  cat(sprintf("  %-25s: Raw I = %.4f (p=%.4f) → Residual p = %.4f  %s\n",
              model_name, raw_row$morans_i, raw_p, residual_p, icon))
  
  tibble::tibble(
    model         = model_name,
    raw_morans_i  = raw_row$morans_i,
    raw_p         = raw_p,
    residual_p    = residual_p,
    interpretation = interpretation
  )
}


#' Compile a full diagnostic report for all fitted models
#'
#' Loops over a named list of fitted models, running pre-DHARMa checks and
#' DHARMa tests for each, and assembles a single summary tibble.
#'
#' @param models Named list of fitted `glm` or `glmmTMB` objects. Names are
#'   used as model identifiers.
#' @param data A data frame that includes X and Y coordinate columns.
#' @param raw_spatial A tibble from `calculate_morans_i()` on the raw data.
#'   Pass `NULL` to skip spatial-pattern comparison.
#' @param output_dir Character. Directory for DHARMa plots. Default `NULL`
#'   (plots not saved).
#' @param n_sim Integer. DHARMa simulations. Default 1000.
#' @return A named list:
#'   * `dharma_summary` – tibble with one row per model
#'   * `spatial_comparison` – tibble from `test_spatial_pattern_comparison()`
#'   * `pre_dharma` – named list of pre-DHARMa check results
#' @examples
#' report <- compile_diagnostic_report(
#'   models      = models_glm,
#'   data        = modeling_data,
#'   raw_spatial = raw_spatial_results,
#'   output_dir  = project_path("Lichens", "model_outputs")
#' )
compile_diagnostic_report <- function(models,
                                      data,
                                      raw_spatial  = NULL,
                                      output_dir   = NULL,
                                      n_sim        = get_project_config()$modeling$n_sim_dharma) {
  stopifnot(is.list(models), length(models) > 0,
            is.data.frame(data))
  
  lichen_banner("Compiling Full Diagnostic Report")
  
  coords_df <- data[, c("X", "Y")]
  
  dharma_summary    <- tibble::tibble()
  spatial_comparison <- tibble::tibble()
  pre_dharma        <- list()
  
  for (nm in names(models)) {
    model <- models[[nm]]
    
    # Pre-DHARMa
    pre_dharma[[nm]] <- run_pre_dharma_checks(model, model_name = nm)
    
    if (pre_dharma[[nm]]$status == "INVALID") {
      lichen_warning("Skipping DHARMa for '", nm, "' (pre-DHARMa INVALID)")
      next
    }
    
    # DHARMa
    diag_result <- run_dharma_tests(model,
                                    coords_df  = coords_df,
                                    model_name = nm,
                                    output_dir = output_dir,
                                    n_sim      = n_sim)
    
    dharma_summary <- dplyr::bind_rows(dharma_summary, diag_result$diagnostics)
    
    # Spatial comparison
    if (!is.null(raw_spatial)) {
      cmp <- test_spatial_pattern_comparison(raw_spatial, diag_result, nm)
      if (!is.null(cmp)) {
        spatial_comparison <- dplyr::bind_rows(spatial_comparison, cmp)
      }
    }
  }
  
  n_valid   <- sum(sapply(pre_dharma, function(x) x$status == "VALID"))
  n_passed  <- if (nrow(dharma_summary) > 0) {
    sum(apply(dharma_summary[, c("dispersion_ok", "zero_inflation_ok",
                                 "outliers_ok", "spatial_ok")], 1, all))
  } else 0L
  
  lichen_message("Pre-DHARMa VALID: ", n_valid, "/", length(models))
  lichen_message("DHARMa all-pass:  ", n_passed, "/", nrow(dharma_summary))
  
  list(
    dharma_summary     = dharma_summary,
    spatial_comparison = spatial_comparison,
    pre_dharma         = pre_dharma
  )
}


# Internal null-coalescing operator (re-declared for standalone sourcing)
`%||%` <- function(a, b) if (!is.null(a)) a else b