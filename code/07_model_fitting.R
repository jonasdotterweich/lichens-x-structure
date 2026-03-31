# =============================================================================
# 07_model_fitting.R
# Functions for fitting GLM and NB models
# =============================================================================
# Provides: fit_glm_model(), fit_nb_model(), run_pre_dharma_checks(),
#           extract_model_summary()
# =============================================================================

library(dplyr)
library(tibble)
library(glmmTMB)

source(here::here("utils.R"))


# -----------------------------------------------------------------------------
# PUBLIC FUNCTIONS
# -----------------------------------------------------------------------------

#' Fit a binomial GLM (logit link) for presence/absence response
#'
#' @param data A data frame or tibble containing all response and predictor
#'   columns (first argument – pipe-friendly).
#' @param response Character. Name of the binary (0/1) response column.
#' @param predictors Character vector of predictor column names.
#' @param family A `family` object. Default `binomial(link = "logit")`.
#' @return A fitted `glm` object with an additional attribute `model_name` set
#'   to `response`.
#' @examples
#' m_parmelia <- fit_glm_model(modeling_data,
#'                             response   = "parmelia_agg_presence",
#'                             predictors = predictors_scaled)
fit_glm_model <- function(data,
                          response,
                          predictors,
                          family = stats::binomial(link = "logit")) {
  stopifnot(is.data.frame(data), is.character(response),
            is.character(predictors), length(predictors) > 0)
  
  if (!response %in% colnames(data)) {
    stop("Response column '", response, "' not found in data.")
  }
  missing_preds <- setdiff(predictors, colnames(data))
  if (length(missing_preds) > 0) {
    stop("Predictor column(s) not found: ", paste(missing_preds, collapse = ", "))
  }
  
  formula_str <- paste(response, "~", paste(predictors, collapse = " + "))
  model_formula <- stats::as.formula(formula_str)
  
  model <- stats::glm(model_formula, data = data, family = family)
  attr(model, "model_name") <- response
  
  lichen_message("Fitted binomial GLM: ", response,
                 " (AIC = ", round(stats::AIC(model), 1), ")")
  model
}


#' Fit a negative-binomial GLM for count response (calicioids richness)
#'
#' Uses `glmmTMB::glmmTMB()` with `family = nbinom2` (NB2 parameterisation).
#'
#' @param data A data frame or tibble (pipe-friendly).
#' @param response Character. Name of the count response column.
#'   Default `"calicioids_richness"`.
#' @param predictors Character vector of predictor column names.
#' @return A fitted `glmmTMB` object with attribute `model_name` set to
#'   `response`.
#' @examples
#' m_calicioids <- fit_nb_model(modeling_data,
#'                              response   = "calicioids_richness",
#'                              predictors = predictors_scaled)
fit_nb_model <- function(data,
                         response   = "calicioids_richness",
                         predictors) {
  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    stop("Package 'glmmTMB' is required. Install from CRAN with: install.packages(\"glmmTMB\")")
  }
  stopifnot(is.data.frame(data), is.character(response),
            is.character(predictors), length(predictors) > 0)
  
  if (!response %in% colnames(data)) {
    stop("Response column '", response, "' not found in data.")
  }
  missing_preds <- setdiff(predictors, colnames(data))
  if (length(missing_preds) > 0) {
    stop("Predictor column(s) not found: ", paste(missing_preds, collapse = ", "))
  }
  
  formula_str   <- paste(response, "~", paste(predictors, collapse = " + "))
  model_formula <- stats::as.formula(formula_str)
  
  model <- glmmTMB::glmmTMB(model_formula, data = data,
                            family = glmmTMB::nbinom2)
  attr(model, "model_name") <- response
  
  lichen_message("Fitted negative-binomial GLM: ", response,
                 " (AIC = ", round(stats::AIC(model), 1), ")")
  model
}


#' Run pre-DHARMa model sanity checks
#'
#' Evaluates four model diagnostics before running the full DHARMa simulation:
#' 1. Maximum absolute coefficient
#' 2. Maximum standard error
#' 3. Convergence
#' 4. Pseudo-R² (binomial) or dispersion parameter (NB)
#'
#' @param model A fitted `glm` or `glmmTMB` object.
#' @param model_name Character. Human-readable name for messages.
#'   If `NULL` (default), uses `attr(model, "model_name")` or the response
#'   variable name.
#' @return A named list: `status` ("VALID", "BORDERLINE", or "INVALID"),
#'   `max_coef`, `max_se`, `flags` (character vector of individual results).
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
#' @param model A fitted `glm` or `glmmTMB` object.
#' @param model_name Character. Label for the `model` column in the output.
#'   Default: `attr(model, "model_name")`.
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

#' Fit a GLMM (Generalized Linear Mixed Model) using glmmTMB
#'
#' Extends \code{fit_glm_model()} by adding random effects to the formula.
#' Uses \code{glmmTMB::glmmTMB()} for all families, including binomial and
#' negative-binomial.
#'
#' **Why use a GLMM instead of a GLM?**
#' Add random effects when your plots are grouped (e.g. by forest region,
#' survey cluster, or observer) and you want to account for non-independence
#' within groups.  The syntax change is minimal: specify the grouping
#' structure in \code{random_effects} and \code{glmmTMB} handles the rest.
#' DHARMa diagnostics work identically for GLMMs.
#'
#' @param data A data frame or tibble (pipe-friendly).
#' @param response Character. Name of the response column.
#' @param predictors Character vector of fixed-effect predictor column names.
#' @param random_effects Character vector of random effect terms in
#'   \code{lme4}/\code{glmmTMB} formula syntax.  Examples:
#'   \code{c("(1|region)")} – random intercept by region;
#'   \code{c("(1|region)", "(1|observer)")} – two crossed random intercepts;
#'   \code{c("(dbh_max|region)")} – random slope + intercept.
#'   Default: \code{get_project_config()$glmm$random_effects}.
#' @param family A family object.  Default \code{binomial(link = "logit")} for
#'   binary responses.  Use \code{glmmTMB::nbinom2} for count responses.
#' @return A fitted \code{glmmTMB} object with attribute \code{model_name}.
#' @examples
#' m_parmelia_glmm <- fit_glmm_model(
#'   modeling_data,
#'   response       = "parmelia_agg_presence",
#'   predictors     = predictors_scaled,
#'   random_effects = c("(1|region)")
#' )
#'
#' # Count response with NB family
#' m_calicioids_glmm <- fit_glmm_model(
#'   modeling_data,
#'   response       = "calicioids_richness",
#'   predictors     = predictors_scaled,
#'   random_effects = c("(1|region)"),
#'   family         = glmmTMB::nbinom2
#' )
fit_glmm_model <- function(data,
                           response,
                           predictors,
                           random_effects = get_project_config()$glmm$random_effects,
                           family         = stats::binomial(link = "logit")) {
  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    stop("Package 'glmmTMB' is required. Install from CRAN with: install.packages(\"glmmTMB\")")
  }
  stopifnot(is.data.frame(data), is.character(response),
            is.character(predictors), length(predictors) > 0)
  
  if (length(random_effects) == 0) {
    stop("fit_glmm_model() requires at least one random effect.\n",
         "Set glmm$random_effects in utils.R or pass random_effects directly.\n",
         "For a standard GLM with no random effects use fit_glm_model().")
  }
  
  if (!response %in% colnames(data)) {
    stop("Response column '", response, "' not found in data.")
  }
  missing_preds <- setdiff(predictors, colnames(data))
  if (length(missing_preds) > 0) {
    stop("Predictor column(s) not found: ", paste(missing_preds, collapse = ", "))
  }
  
  fixed_part    <- paste(predictors, collapse = " + ")
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




# Null-coalescing operator for base R
`%||%` <- function(a, b) if (!is.null(a)) a else b

.detect_family <- function(model) {
  if (inherits(model, "glmmTMB")) {
    fam_name <- tryCatch(model$modelInfo$family$family, error = function(e) "unknown")
    if (grepl("nbinom", fam_name, ignore.case = TRUE)) return("nb")
    return("glmm")
  }
  fam <- tryCatch(stats::family(model)$family, error = function(e) "unknown")
  if (fam == "binomial") "binomial" else fam
}

.extract_coef_table <- function(model, family_type) {
  if (family_type == "nb") {
    summary(model)$coefficients$cond
  } else {
    summary(model)$coefficients
  }
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

.check_convergence <- function(model, family_type) {
  if (family_type == "nb") {
    conv <- model$fit$convergence
    converged <- conv == 0
  } else {
    converged <- model$converged
    conv <- ifelse(converged, 0, 1)
  }
  icon <- if (converged) "\u2705" else "\U1F6A8"
  cat("  Converged:", converged, icon, "\n")
  if (converged) "PASS" else "FAIL"
}

.check_fit_quality <- function(model, family_type) {
  if (family_type == "binomial") {
    pseudo_r2 <- 1 - (model$deviance / model$null.deviance)
    flag <- dplyr::case_when(
      pseudo_r2 > 0.95 ~ "FAIL",
      pseudo_r2 > 0.85 ~ "WARNING",
      TRUE              ~ "PASS"
    )
    icon <- if (flag == "PASS") "\u2705" else if (flag == "WARNING") "\u26A0\uFE0F" else "\U1F6A8"
    cat("  Pseudo-R\u00b2:", round(pseudo_r2, 3), icon, "\n")
    flag
  } else if (family_type %in% c("nb", "glmm")) {
    disp <- tryCatch(sigma(model), error = function(e) NA_real_)
    if (is.na(disp)) {
      cat("  Dispersion: not available\n")
      return("PASS")
    }
    flag <- if (disp < 0.1 || disp > 100) "WARNING" else "PASS"
    icon <- if (flag == "PASS") "\u2705" else "\u26A0\uFE0F"
    cat("  Dispersion:", round(disp, 3), icon, "\n")
    flag
  } else {
    "PASS"
  }
}