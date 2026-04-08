# =============================================================================
# 07_model_fitting.R
# Functions for fitting Generalized Linear Mixed Models (GLMMs)
# =============================================================================
# Provides: fit_glmm_model(), run_pre_dharma_checks(), extract_model_summary()
#
# All models are fitted with glmmTMB, which handles both:
#   - Binary (presence/absence)  → family = binomial(link = "logit")  [default]
#   - Count (richness)           → family = glmmTMB::nbinom2
#
# Random effects are set in utils.R → glmm$random_effects.
# See that section for syntax examples and what to change per dataset.
# =============================================================================

library(dplyr)
library(tibble)
library(glmmTMB)

source(here::here("utils.R"))


# -----------------------------------------------------------------------------
# PUBLIC FUNCTIONS
# -----------------------------------------------------------------------------

#' Fit a model (GLM-first, GLMM-optional) using glmmTMB
#'
#' This is the single model-fitting function for this framework.
#' It always uses glmmTMB as the engine.
#'
#' - If cfg$model$use_random_effects is FALSE (default): fits a fixed-effects-only
#'   model (i.e. GLM behaviour) with glmmTMB.
#' - If TRUE: fits a mixed model with the random effects specified in cfg$glmm$random_effects
#'   or passed explicitly.
#'
#' @param data A data frame or tibble (pipe-friendly).
#' @param response Character. Name of the response column.
#' @param predictors Character vector of fixed-effect predictor column names.
#' @param random_effects Character vector of random effect terms in
#'   lme4/glmmTMB formula syntax. Default: get_project_config()$glmm$random_effects.
#' @param family A family object. Default binomial(link="logit").
#'   For richness/count responses, use glmmTMB::nbinom2.
#' @return A fitted glmmTMB object with attribute model_name.
fit_model <- function(data,
                      response,
                      predictors,
                      random_effects = get_project_config()$glmm$random_effects,
                      family         = stats::binomial(link = "logit")) {
  
  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    stop("Package 'glmmTMB' is required. Install from CRAN: install.packages(\"glmmTMB\")")
  }
  stopifnot(is.data.frame(data), is.character(response),
            is.character(predictors), length(predictors) > 0)
  
  cfg <- get_project_config()
  use_re <- isTRUE(cfg$model$use_random_effects)
  
  if (!response %in% colnames(data)) {
    stop("Response column '", response, "' not found in data.")
  }
  missing_preds <- setdiff(predictors, colnames(data))
  if (length(missing_preds) > 0) {
    stop("Predictor column(s) not found: ", paste(missing_preds, collapse = ", "))
  }
  
  fixed_part <- paste(predictors, collapse = " + ")
  
  if (!use_re) {
    # GLM-first mode: no random effects in the formula
    lichen_message("Random effects disabled (GLM mode) → fitting fixed-effects-only model: ", response)
    
    formula_str   <- paste0(response, " ~ ", fixed_part)
    model_formula <- stats::as.formula(formula_str)
    
    model <- glmmTMB::glmmTMB(model_formula, data = data, family = family)
    attr(model, "model_name") <- response
    
    lichen_message("Fitted model: ", response, " (AIC = ", round(stats::AIC(model), 1), ")")
    return(model)
  }
  
  # GLMM mode
  if (length(random_effects) == 0) {
    stop(
      "cfg$model$use_random_effects is TRUE, but no random effects were provided.\n",
      "Set glmm$random_effects in utils.R or pass random_effects directly."
    )
  }
  
  re_part       <- paste(random_effects, collapse = " + ")
  formula_str   <- paste0(response, " ~ ", fixed_part, " + ", re_part)
  model_formula <- stats::as.formula(formula_str)
  
  model <- glmmTMB::glmmTMB(model_formula, data = data, family = family)
  attr(model, "model_name") <- response
  
  lichen_message("Fitted GLMM: ", response,
                 " (AIC = ", round(stats::AIC(model), 1), ")")
  lichen_message("Random effects: ", paste(random_effects, collapse = ", "))
  model
}



#' Run pre-DHARMa model sanity checks
#'
#' Evaluates four diagnostics on a fitted \code{glmmTMB} model before running
#' the full DHARMa simulation:
#' \enumerate{
#'   \item Maximum absolute coefficient
#'   \item Maximum standard error
#'   \item Convergence
#'   \item Dispersion (via \code{sigma()})
#' }
#'
#' @param model A fitted \code{glmmTMB} object (from \code{fit_glmm_model()}).
#' @param model_name Character. Human-readable name for messages.
#'   If \code{NULL} (default), uses \code{attr(model, "model_name")} or the
#'   response variable name.
#' @return A named list: \code{status} ("VALID", "BORDERLINE", or "INVALID"),
#'   \code{max_coef}, \code{max_se}, \code{flags} (character vector).
#' @examples
#' checks <- run_pre_dharma_checks(m_parmelia)
#' checks$status
run_pre_dharma_checks <- function(model, model_name = NULL) {
  family_type <- .detect_family(model)
  
  if (is.null(model_name)) {
    model_name <- attr(model, "model_name") %||%
      tryCatch(deparse(stats::formula(model)[[2]]), error = function(e) "model")
  }
  
  lichen_banner(paste("Pre-DHARMa Checks:", model_name))
  
  # Extract coefficient table
  coef_table <- .extract_coef_table(model, family_type)
  
  # 1. Coefficients
  max_coef      <- max(abs(coef_table[, "Estimate"]))
  max_coef_name <- rownames(coef_table)[which.max(abs(coef_table[, "Estimate"]))]
  coef_flag <- .check_threshold(max_coef, 10, 5,
                                label = paste0("Max coef: ", round(max_coef, 2),
                                               " (", max_coef_name, ")"))
  
  # 2. Standard errors
  max_se      <- max(coef_table[, "Std. Error"])
  max_se_name <- rownames(coef_table)[which.max(coef_table[, "Std. Error"])]
  se_flag <- .check_threshold(max_se, 10, 5,
                              label = paste0("Max SE: ", round(max_se, 2),
                                             " (", max_se_name, ")"))
  
  # 3. Convergence
  conv_flag <- .check_convergence(model, family_type)
  
  # 4. Pseudo-R² or dispersion
  r2_flag <- .check_fit_quality(model, family_type)
  
  all_flags <- c(coef_flag, se_flag, conv_flag, r2_flag)
  status <- dplyr::case_when(
    any(all_flags == "FAIL")    ~ "INVALID",
    any(all_flags == "WARNING") ~ "BORDERLINE",
    TRUE                         ~ "VALID"
  )
  
  cat("→ Overall:", status, "\n")
  
  list(status   = status,
       max_coef = max_coef,
       max_se   = max_se,
       flags    = all_flags)
}


#' Extract a tidy model summary tibble
#'
#' @param model A fitted \code{glmmTMB} object (from \code{fit_glmm_model()}).
#' @param model_name Character. Label for the \code{model} column in the output.
#'   Default: \code{attr(model, "model_name")}.
#' @return A tibble with columns: model, term, estimate, std_error, statistic,
#'   p_value, sig.
#' @examples
#' tidy_coefs <- extract_model_summary(m_parmelia)
extract_model_summary <- function(model, model_name = NULL) {
  if (is.null(model_name)) {
    model_name <- attr(model, "model_name") %||% "model"
  }
  
  family_type <- .detect_family(model)
  coef_table  <- .extract_coef_table(model, family_type)
  
  result <- tibble::as_tibble(coef_table, rownames = "term") |>
    dplyr::rename_with(
      ~dplyr::case_when(
        . == "Estimate"   ~ "estimate",
        . == "Std. Error" ~ "std_error",
        grepl("z|t value", .) ~ "statistic",
        grepl("Pr\\(", .) ~ "p_value",
        TRUE              ~ .
      )
    ) |>
    dplyr::mutate(
      model = model_name,
      sig   = sig_stars(p_value)
    ) |>
    dplyr::select(model, term, estimate, std_error, statistic, p_value, sig)
  
  result
}


# -----------------------------------------------------------------------------
# INTERNAL HELPERS
# -----------------------------------------------------------------------------

# Null-coalescing operator for base R
`%||%` <- function(a, b) if (!is.null(a)) a else b

# All models are glmmTMB; return "nb" for negative-binomial, "glmm" otherwise.
.detect_family <- function(model) {
  fam_name <- tryCatch(model$modelInfo$family$family, error = function(e) "unknown")
  if (grepl("nbinom", fam_name, ignore.case = TRUE)) "nb" else "glmm"
}

# glmmTMB stores fixed-effect coefficients under $coefficients$cond for both
# binomial and NB families.
.extract_coef_table <- function(model, family_type) {
  summary(model)$coefficients$cond
}

.check_threshold <- function(value, fail_at, warn_at, label) {
  flag <- dplyr::case_when(
    value > fail_at ~ "FAIL",
    value > warn_at ~ "WARNING",
    TRUE             ~ "PASS"
  )
  icon <- dplyr::case_when(
    flag == "FAIL"    ~ "\U1F6A8",
    flag == "WARNING" ~ "\u26A0\uFE0F",
    TRUE              ~ "\u2705"
  )
  cat("  ", label, icon, "\n", sep = " ")
  flag
}

# glmmTMB convergence is stored in model$fit$convergence (0 = converged).
.check_convergence <- function(model, family_type) {
  conv      <- tryCatch(model$fit$convergence, error = function(e) 1L)
  converged <- isTRUE(conv == 0L)
  icon <- if (converged) "\u2705" else "\U1F6A8"
  cat("  Converged:", converged, icon, "\n")
  if (converged) "PASS" else "FAIL"
}

# Use sigma() (residual SD / dispersion parameter) for all glmmTMB models.
.check_fit_quality <- function(model, family_type) {
  disp <- tryCatch(sigma(model), error = function(e) NA_real_)
  if (is.na(disp)) {
    cat("  Dispersion: not available\n")
    return("PASS")
  }
  flag <- if (disp < 0.1 || disp > 100) "WARNING" else "PASS"
  icon <- if (flag == "PASS") "\u2705" else "\u26A0\uFE0F"
  cat("  Dispersion:", round(disp, 3), icon, "\n")
  flag
}