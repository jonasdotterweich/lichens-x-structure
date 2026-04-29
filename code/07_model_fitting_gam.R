# =============================================================================
# 07_model_fitting_gam.R
# Functions for fitting Generalized Additive Models (GAMs)
# =============================================================================
# Provides:
#   - fit_gam_model()
#   - extract_gam_parametric_summary()
#   - extract_gam_smooth_summary()
#   - summarize_gam_fit()
#
# This script is the GAM analogue of the GLM/GLMM model-fitting module.
# It uses the mgcv package as the modeling engine and supports:
#
#   - Binary (presence/absence) responses:
#       family = binomial(link = "logit")
#
#   - Count (richness) responses with overdispersion:
#       family = mgcv::nb()   # negative binomial with theta estimated by mgcv
#
# Key GAM concepts used here:
#   - Predictors can be included as:
#       * parametric terms (linear effects and factors)  → linear_terms
#       * smooth terms (nonlinear effects)               → smooth_terms, fitted via s(x)
#   - Smoothing is estimated via REML by default (method = "REML"), with optional
#     shrinkage/term selection enabled (select = TRUE).
#   - Optional spatial control can be added via a 2D smooth s(X, Y) (spatial = TRUE).
#
# Outputs:
#   - The extract_* functions return tidy summaries for:
#       * parametric effects (estimates, SE, t/z, p-values)
#       * smooth effects (EDF, reference DF, F/Chi.sq, p-values)
#   - summarize_gam_fit() returns a compact one-row overview (n, AIC, deviance explained).
#
# Notes:
#   - Random effects are not handled in this module (GAMMs are a later extension).
#     The current scope is GAMs for robust, interpretable nonlinear relationships.
#   - Diagnostics (gam.check, concurvity, Moran's I on residuals) are implemented
#     in the GAM diagnostics script (08_model_diagnostics_gam.R).
# =============================================================================

library(dplyr)
library(tibble)
library(mgcv)

source(here::here("utils.R"))


# -----------------------------------------------------------------------------
# PUBLIC FUNCTIONS
# -----------------------------------------------------------------------------
#' Fit a GAM using mgcv::gam()
#'
#' @param data data.frame
#' @param response character response column name
#' @param linear_terms character vector of linear/factor terms
#' @param smooth_terms character vector of continuous predictors to smooth
#' @param family mgcv family, e.g. binomial() or mgcv::nb()
#' @param spatial logical; add s(X,Y)
#' @param coords_cols character length-2; default c("X","Y")
#' @param k integer basis dimension for 1D smooths
#' @param k_spatial integer basis dimension for spatial smooth (if spatial = TRUE)
#' @param method smoothing parameter method; default "REML"
#' @param select logical; shrinkage/term selection
fit_gam_model <- function(data,
                          response,
                          linear_terms = character(0),
                          smooth_terms = character(0),
                          family,
                          spatial = FALSE,
                          coords_cols = c("X", "Y"),
                          k = 5,
                          k_spatial = 30,
                          method = "REML",
                          select = TRUE) {
  
  stopifnot(is.data.frame(data), is.character(response), length(response) == 1)
  
  if (!response %in% names(data)) {
    stop("Response column '", response, "' not found.")
  }
  
  all_terms <- unique(c(linear_terms, smooth_terms))
  missing_terms <- setdiff(all_terms, names(data))
  if (length(missing_terms) > 0) {
    stop("Term(s) not found in data: ", paste(missing_terms, collapse = ", "))
  }
  
  # Build smooth part
  smooth_part <- character(0)
  if (length(smooth_terms) > 0) {
    smooth_part <- paste0("s(", smooth_terms, ", k = ", k, ")")
  }
  
  # Optional spatial smooth
  spatial_part <- character(0)
  if (isTRUE(spatial)) {
    if (!all(coords_cols %in% names(data))) {
      stop("Spatial coords missing: ", paste(setdiff(coords_cols, names(data)), collapse = ", "))
    }
    spatial_part <- paste0("s(", coords_cols[1], ", ", coords_cols[2], ", k = ", k_spatial, ")")
  }
  
  rhs_terms <- c(linear_terms, smooth_part, spatial_part)
  rhs_terms <- rhs_terms[nzchar(rhs_terms)]
  
  if (length(rhs_terms) == 0) stop("No predictors specified for GAM.")
  
  formula_str <- paste(response, "~", paste(rhs_terms, collapse = " + "))
  model_formula <- stats::as.formula(formula_str)
  
  model <- mgcv::gam(
    formula = model_formula,
    data = data,
    family = family,
    method = method,
    select = select
  )
  
  attr(model, "model_name") <- response
  lichen_message("Fitted GAM: ", response)
  model
}

#' Extract parametric coefficient summary from a GAM
extract_gam_parametric_summary <- function(model, model_name = NULL) {
  if (is.null(model_name)) model_name <- attr(model, "model_name") %||% "model"
  sm <- summary(model)
  
  if (is.null(sm$p.table) || nrow(sm$p.table) == 0) {
    return(tibble::tibble(
      model = character(0), term = character(0),
      estimate = numeric(0), std_error = numeric(0),
      statistic = numeric(0), p_value = numeric(0), sig = character(0)
    ))
  }
  
  df <- as.data.frame(sm$p.table)
  df$term <- rownames(df); rownames(df) <- NULL
  
  # mgcv uses either t or z columns depending on family
  stat_col <- intersect(c("t value", "z value"), names(df))[1]
  p_col    <- intersect(c("Pr(>|t|)", "Pr(>|z|)"), names(df))[1]
  
  if (is.na(stat_col) || is.na(p_col)) {
    stop(
      "Unexpected mgcv parametric table columns. Found: ",
      paste(names(df), collapse = ", ")
    )
  }
  
  out <- tibble::as_tibble(df) |>
    dplyr::transmute(
      model     = model_name,
      term      = term,
      estimate  = .data[["Estimate"]],
      std_error = .data[["Std. Error"]],
      statistic = .data[[stat_col]],
      p_value   = .data[[p_col]],
      sig       = sig_stars(p_value)
    )
  out
}

#' Extract smooth term summary from a GAM
extract_gam_smooth_summary <- function(model, model_name = NULL) {
  if (is.null(model_name)) model_name <- attr(model, "model_name") %||% "model"
  sm <- summary(model)
  
  if (is.null(sm$s.table) || nrow(sm$s.table) == 0) {
    return(tibble::tibble(
      model = character(0), smooth = character(0),
      edf = numeric(0), ref_df = numeric(0),
      statistic = numeric(0), p_value = numeric(0), sig = character(0)
    ))
  }
  
  df <- as.data.frame(sm$s.table)
  df$smooth <- rownames(df); rownames(df) <- NULL
  
  # Column names can differ slightly; handle common cases
  stat_col <- intersect(c("F", "Chi.sq"), names(df))[1]
  p_col <- intersect(c("p-value", "p.value"), names(df))[1]
  
  out <- tibble::as_tibble(df) |>
    dplyr::transmute(
      model = model_name,
      smooth = smooth,
      edf = edf,
      ref_df = `Ref.df`,
      statistic = .data[[stat_col]],
      p_value = .data[[p_col]],
      sig = sig_stars(p_value)
    )
  
  out
}

# Null-coalescing operator
`%||%` <- function(a, b) if (!is.null(a)) a else b



summarize_gam_fit <- function(model, model_name = NULL) {
  if (is.null(model_name)) model_name <- attr(model, "model_name") %||% "model"
  
  sm <- summary(model)
  
  aic <- tryCatch(AIC(model), error = function(e) NA_real_)
  dev <- tryCatch(sm$dev.expl, error = function(e) NA_real_)
  r2  <- tryCatch(sm$r.sq, error = function(e) NA_real_)
  
  n <- tryCatch(as.integer(model$df.null + 1), error = function(e) NA_integer_)
  if (is.na(n)) n <- tryCatch(nrow(model$model), error = function(e) NA_integer_)
  
  fam <- tryCatch(model$family$family, error = function(e) NA_character_)
  lnk <- tryCatch(model$family$link, error = function(e) NA_character_)
  
  tibble::tibble(
    model = model_name,
    n = n,
    AIC = aic,
    dev_explained = dev,
    r_sq = r2,
    family = fam,
    link = lnk
  )
}