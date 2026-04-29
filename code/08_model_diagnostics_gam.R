# =============================================================================
# 08_model_diagnostics_gam.R
# Diagnostics and validation utilities for Generalized Additive Models (GAMs)
# =============================================================================
# Provides:
#   - save_gamcheck_plot()
#   - save_concurvity()
#   - calculate_residual_morans_i()
#   - run_gam_diagnostics()
#
# Purpose
# -------
# This script implements the GAM-equivalent diagnostic workflow to the GLM/GLMM
# (glmmTMB + DHARMa) checks in the original framework. GAMs (mgcv) have a
# different set of standard validation tools, so this module focuses on:
#
#   1) Basis dimension adequacy & residual patterns:
#        mgcv::gam.check()
#      - Produces diagnostic plots (QQ / residuals vs fitted) and the "k-index"
#        table to assess whether basis dimensions (k) are sufficient.
#
#   2) Nonlinear multicollinearity (GAM analogue of VIF):
#        mgcv::concurvity(full = TRUE)
#      - Quantifies how well each smooth/term can be approximated by other terms.
#        High concurvity can reduce interpretability and destabilize partial effects.
#
#   3) Spatial autocorrelation diagnostics:
#        Moran's I on model residuals (typically Pearson residuals)
#      - Helps assess whether there is remaining spatial structure not captured by
#        predictors (and whether adding a spatial smooth s(X, Y) may be warranted).
#
# Output files
# ------------
# run_gam_diagnostics() writes standardized outputs to a user-specified reports
# directory (e.g. outputs/Sumava/reports_gam/), including:
#
#   - gamcheck_<model_name>.png
#       Diagnostic plots produced by mgcv::gam.check()
#
#   - concurvity_<model_name>.csv
#       A tidy table of concurvity values (worst/observed/estimate) per term
#
#   - fit_<model_name>.csv
#       One-row overview from summarize_gam_fit() (AIC, deviance explained, etc.)
#
#   - morans_resid_<model_name>.csv   (optional / depends on calculate_morans_i())
#       Moran's I results for residual spatial autocorrelation
#
# Notes / limitations
# -------------------
# - This module uses mgcv-native diagnostics rather than DHARMa simulations.
#   DHARMa may work for some mgcv GAMs, but mgcv::gam.check() + concurvity are
#   the canonical checks and integrate directly with smoothness assumptions.
# - Spatial residual checks here assume coordinate columns (default: X, Y) exist
#   in the modeling data. If coordinates differ, pass coords_cols explicitly.
# - Moran's I implementation is delegated to calculate_morans_i() if it exists
#   in the current session (typically defined/sourced via utils.R in this repo).
# =============================================================================

# =============================================================================
# 08_model_diagnostics_gam.R
# Diagnostics for GAM models (mgcv)
# =============================================================================

library(dplyr)
library(readr)
library(tibble)
library(mgcv)

source(here::here("utils.R"))
source(here::here("code", "07_model_fitting_gam.R"))


# -----------------------------------------------------------------------------
# Spatial autocorrelation (Moran's I) helpers
# -----------------------------------------------------------------------------

#' Build an inverse-distance weight matrix from X/Y coordinates
#' - Diagonal is set to 0 (no self-weight)
#' - Uses 1 / distance
#' - Handles duplicate coordinates by adding a tiny jitter if needed
.build_inv_dist_weights <- function(coords_df, eps = 1e-12) {
  coords <- as.matrix(coords_df)
  if (ncol(coords) != 2) stop("coords_df must have exactly two columns (X, Y).")
  
  # Compute pairwise Euclidean distances
  d <- as.matrix(stats::dist(coords))
  
  # Avoid division by zero (duplicates): jitter very slightly if needed
  if (any(d == 0 & row(d) != col(d))) {
    coords2 <- coords
    coords2 <- coords2 + matrix(stats::rnorm(length(coords2), sd = 1e-6), ncol = 2)
    d <- as.matrix(stats::dist(coords2))
  }
  
  w <- 1 / pmax(d, eps)
  diag(w) <- 0
  w
}

#' Calculate Moran's I for one or more variables using inverse-distance weights
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
      variable   = var,
      morans_i   = round(moran$observed, 4),
      expected_i = round(moran$expected, 4),
      sd         = round(moran$sd, 4),
      p_value    = round(moran$p.value, 4),
      sig        = sig_stars(moran$p.value)
    )
  })
  
  cat(sprintf("%-30s | I = %7.4f | p = %.4f  %s\n",
              results$variable, results$morans_i,
              results$p_value, results$sig))
  
  n_sig <- sum(results$p_value < 0.05, na.rm = TRUE)
  lichen_message("Significant spatial clustering (p < 0.05): ",
                 n_sig, "/", nrow(results), " variables")
  results
}


#' Save gam.check diagnostic plots to a PNG
save_gamcheck_plot <- function(model, file, width = 1200, height = 900, res = 150) {
  grDevices::png(filename = file, width = width, height = height, res = res)
  on.exit(grDevices::dev.off(), add = TRUE)
  
  # mgcv::gam.check prints and plots
  mgcv::gam.check(model)
}

#' Tidy concurvity output into a data.frame and save
save_concurvity <- function(model, file_csv) {
  conc <- mgcv::concurvity(model, full = TRUE)
  
  # conc is a list/matrix-like structure; we’ll store all 3 rows (worst/observed/estimate)
  conc_df <- as.data.frame(t(conc))
  conc_df$term <- rownames(conc_df)
  rownames(conc_df) <- NULL
  
  # Reorder columns: term first
  conc_df <- conc_df |>
    dplyr::select(term, dplyr::everything())
  
  readr::write_csv(conc_df, file_csv)
  conc_df
}

#' Residual Moran's I (simple distance-based neighbors)
#' NOTE: this uses the same "calculate_morans_i" concept but on residuals.
#' If you already have calculate_morans_i() in utils.R, we can reuse it instead.
calculate_residual_morans_i <- function(model, data, coords_cols = c("X", "Y")) {
  if (!all(coords_cols %in% names(data))) {
    stop("Coordinate column(s) not found: ", paste(coords_cols, collapse = ", "))
  }
  
  tmp <- data
  tmp$resid_pearson <- residuals(model, type = "pearson")
  
  calculate_morans_i(
    tmp,
    variables = c("resid_pearson"),
    coords_cols = coords_cols
  )
}


#' Run standard GAM diagnostics and write outputs
run_gam_diagnostics <- function(model,
                                model_name,
                                data,
                                out_dir_reports,
                                coords_cols = c("X", "Y"),
                                run_morans_i = TRUE) {
  
  dir.create(out_dir_reports, recursive = TRUE, showWarnings = FALSE)
  
  # 1) Save gam.check plot
  file_gamcheck <- file.path(out_dir_reports, paste0("gamcheck_", model_name, ".png"))
  save_gamcheck_plot(model, file_gamcheck)
  
  # 2) Save concurvity table
  file_conc <- file.path(out_dir_reports, paste0("concurvity_", model_name, ".csv"))
  conc_df <- save_concurvity(model, file_conc)
  
  # 3) Save fit overview row
  fit_row <- summarize_gam_fit(model, model_name)
  # caller (runner) can bind_rows these; we also write per-model snapshot
  readr::write_csv(fit_row, file.path(out_dir_reports, paste0("fit_", model_name, ".csv")))
  
  # 4) Residual Moran's I (optional)
  moran_out <- NULL
  if (isTRUE(run_morans_i)) {
    moran_out <- calculate_residual_morans_i(model, data, coords_cols = coords_cols)
    # moran_out could be a tibble/data.frame; write if it is
    try({
      readr::write_csv(as.data.frame(moran_out),
                       file.path(out_dir_reports, paste0("morans_resid_", model_name, ".csv")))
    }, silent = TRUE)
  }
  
  invisible(list(
    gamcheck_png = file_gamcheck,
    concurvity_csv = file_conc,
    fit = fit_row,
    concurvity = conc_df,
    morans = moran_out
  ))
}