# =============================================================================
# 09.5_model_tuning_gam.R
# Full GAM tuning across covariates, smooth/linear shape, and spatial options
# =============================================================================

library(dplyr)
library(readr)
library(tibble)
library(tidyr)
library(purrr)

source(here::here("utils.R"))
source(here::here("code", "02_data_loading.R"))
source(here::here("code", "03_data_cleaning.R"))
source(here::here("code", "06_lichen_processing.R"))
source(here::here("code", "07_model_fitting_gam.R"))
source(here::here("code", "08_model_diagnostics_gam.R"))

cfg <- get_project_config()


# ---- Output directories -------------------------------------------------------
OUT_DIR_MODELS_GAM  <- here::here("outputs", "Sumava", "models_gam")
OUT_DIR_REPORTS_GAM <- here::here("outputs", "Sumava", "reports_gam")

OUT_DIR_TUNING_MODELS  <- file.path(OUT_DIR_MODELS_GAM,  "tuning_gam")
OUT_DIR_TUNING_REPORTS <- file.path(OUT_DIR_REPORTS_GAM, "tuning_gam")

dir.create(OUT_DIR_TUNING_MODELS,  recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR_TUNING_REPORTS, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# SECTION 1 – Load and preprocess data
# =============================================================================

lichen_banner("Loading and Preprocessing Data")

PATH_LICHEN    <- here::here("outputs", "Sumava", "lichen_clean.csv")
PATH_STRUCTURE <- here::here("outputs", "Sumava", "structure_clean.csv")

lichen_clean    <- readr::read_csv(PATH_LICHEN,    show_col_types = FALSE)
structure_clean <- readr::read_csv(PATH_STRUCTURE, show_col_types = FALSE)

coords_check <- dplyr::inner_join(
  lichen_clean    |> dplyr::select(plot_id, X, Y),
  structure_clean |> dplyr::select(plot_id, X, Y),
  by     = "plot_id",
  suffix = c(".lichen", ".structure")
)

tol <- 1e-6
coords_check <- coords_check |>
  dplyr::mutate(
    dx           = X.lichen - X.structure,
    dy           = Y.lichen - Y.structure,
    coords_match = (abs(dx) <= tol) & (abs(dy) <= tol)
  )

n_mismatch <- sum(!coords_check$coords_match, na.rm = TRUE)
if (n_mismatch > 0) {
  print(coords_check |> dplyr::filter(!coords_match), n = 50)
  stop("Coordinate mismatch between lichen_clean and structure_clean for ",
       n_mismatch, " plot(s). Fix before continuing.")
} else {
  cat("✓ Coordinates match for all joined plots (tol =", tol, ")\n")
}

structure_clean <- structure_clean |> dplyr::select(-X, -Y)
modeling_data <- dplyr::inner_join(lichen_clean, structure_clean, by = "plot_id")

validate_input_data(
  modeling_data,
  required_cols = c("plot_id", "X", "Y"),
  id_col        = "plot_id",
  check_coords  = TRUE
)

validate_response_variables(modeling_data)

modeling_data <- impute_missing_values_if_needed(
  modeling_data,
  cols = dplyr::where(is.numeric) &
    !dplyr::ends_with("_presence") &
    !dplyr::matches("richness$"),
  threshold_pct = 5,
  strategy      = "median",
  report_path   = file.path(OUT_DIR_TUNING_REPORTS, "qc_imputation_report_tuning.csv")
)


# =============================================================================
# SECTION 2 – Responses and candidate predictor sets
# =============================================================================

responses_bin <- intersect(cfg$lichen_groups$binary, names(modeling_data))
responses_cnt <- intersect(cfg$lichen_groups$count,  names(modeling_data))

if (length(responses_bin) == 0 && length(responses_cnt) == 0) {
  stop("No response columns found in modeling_data. Check utils.R config.")
}

lichen_message("Binary responses:  ", paste(responses_bin, collapse = ", "))
lichen_message("Count  responses:  ", paste(responses_cnt, collapse = ", "))

.keep <- function(...) intersect(c(...), names(modeling_data))
`%||%` <- function(a, b) if (!is.null(a)) a else b

decay_block   <- c("decay2", "decay3", "decay4", "decay5")
decay_present <- intersect(decay_block, names(modeling_data))
smooth_base   <- intersect(c("elevation", "volume_snags", "canopy_cover"), names(modeling_data))

# fixed_linear_terms are always linear.
# nonlinear_candidates are iterated through all smooth/linear allocations.
candidate_predictor_sets <- list(
  structure_core = list(
    fixed_linear_terms   = .keep("dbh_max", "dbh_sd", "n_dead_50cm"),
    nonlinear_candidates = smooth_base
  ),
  structure_core_decay = list(
    fixed_linear_terms   = .keep("dbh_max", "dbh_sd", "n_dead_50cm", decay_present),
    nonlinear_candidates = smooth_base
  ),
  composition = list(
    fixed_linear_terms   = .keep("ba_spruce", "ba_beech", "n_dead_50cm"),
    nonlinear_candidates = smooth_base
  ),
  composition_decay = list(
    fixed_linear_terms   = .keep("ba_spruce", "ba_beech", "n_dead_50cm", decay_present),
    nonlinear_candidates = smooth_base
  ),
  full_nodecay = list(
    fixed_linear_terms   = .keep("dbh_max", "dbh_sd", "n_dead_50cm", "ba_spruce", "ba_beech"),
    nonlinear_candidates = smooth_base
  ),
  full_decay = list(
    fixed_linear_terms   = .keep("dbh_max", "dbh_sd", "n_dead_50cm", "ba_spruce", "ba_beech", decay_present),
    nonlinear_candidates = smooth_base
  ),
  deadwood_centric = list(
    fixed_linear_terms   = .keep("dbh_max", "n_dead_50cm", "ba_spruce"),
    nonlinear_candidates = .keep("volume_snags", "canopy_cover")
  ),
  deadwood_centric_decay = list(
    fixed_linear_terms   = .keep("dbh_max", "n_dead_50cm", "ba_spruce", decay_present),
    nonlinear_candidates = .keep("volume_snags", "canopy_cover")
  )
)

candidate_predictor_sets <- Filter(
  function(x) {
    length(unique(c(x$fixed_linear_terms, x$nonlinear_candidates))) > 0
  },
  candidate_predictor_sets
)

k_values         <- c(4L, 5L, 7L, 10L)
spatial_values   <- c(FALSE, TRUE)
k_spatial_values <- c(15L, 30L, 60L)
TOP_N            <- 3L
DELTA_AIC_TOL    <- 2


# =============================================================================
# SECTION 3 – Model-spec grid + diagnostics helpers
# =============================================================================

.all_subsets <- function(x) {
  if (length(x) == 0) return(list(character(0)))
  out <- vector("list", 2^length(x))
  idx <- 0L
  for (mask in 0:(2^length(x) - 1)) {
    idx <- idx + 1L
    bits <- as.logical(intToBits(mask))[seq_along(x)]
    out[[idx]] <- x[bits]
  }
  out
}

.term_string <- function(x) {
  x <- sort(unique(x))
  if (length(x) == 0) "" else paste(x, collapse = "; ")
}

.make_spec_id <- function(set_name, smooth_terms, linear_terms, k, spatial, k_spatial) {
  paste0(
    "set=", set_name,
    "__s(", paste(sort(unique(smooth_terms)), collapse = "+"), ")",
    "__l(", paste(sort(unique(linear_terms)), collapse = "+"), ")",
    "__k=", k,
    "__sp=", ifelse(isTRUE(spatial), "1", "0"),
    "__ksp=", ifelse(isTRUE(spatial), as.character(k_spatial), "NA")
  )
}

.extract_diag_metrics <- function(model, data) {
  k_index_min <- NA_real_
  k_pvalue_min <- NA_real_
  k_problem_any <- NA
  concurvity_max_worst <- NA_real_
  concurvity_ok <- NA
  resid_morans_i <- NA_real_
  resid_morans_p <- NA_real_
  resid_spatial_ok <- NA

  kcheck <- tryCatch(mgcv::k.check(model), error = function(e) NULL)
  if (!is.null(kcheck) && nrow(kcheck) > 0) {
    k_index <- suppressWarnings(as.numeric(kcheck[, "k-index"]))
    p_vals  <- suppressWarnings(as.numeric(kcheck[, "p-value"]))
    k_index_min  <- suppressWarnings(min(k_index, na.rm = TRUE))
    k_pvalue_min <- suppressWarnings(min(p_vals, na.rm = TRUE))
    k_problem_any <- any(k_index < 1 & p_vals < 0.05, na.rm = TRUE)
  }

  conc <- tryCatch(mgcv::concurvity(model, full = TRUE), error = function(e) NULL)
  if (!is.null(conc) && "worst" %in% rownames(conc)) {
    worst_vals <- suppressWarnings(as.numeric(conc["worst", ]))
    concurvity_max_worst <- suppressWarnings(max(worst_vals, na.rm = TRUE))
    concurvity_ok <- is.finite(concurvity_max_worst) && concurvity_max_worst <= 0.9
  }

  moran <- tryCatch(
    calculate_residual_morans_i(model, data, coords_cols = c("X", "Y")),
    error = function(e) NULL
  )
  if (!is.null(moran) && nrow(moran) > 0) {
    resid_morans_i <- suppressWarnings(as.numeric(moran$morans_i[1]))
    resid_morans_p <- suppressWarnings(as.numeric(moran$p_value[1]))
    resid_spatial_ok <- is.na(resid_morans_p) || resid_morans_p >= 0.05
  }

  k_ok <- ifelse(is.na(k_problem_any), NA, !k_problem_any)
  diagnostics_ok <- isTRUE(k_ok) && isTRUE(concurvity_ok) &&
    (isTRUE(resid_spatial_ok) || is.na(resid_spatial_ok))

  tibble::tibble(
    k_index_min = k_index_min,
    k_pvalue_min = k_pvalue_min,
    k_ok = k_ok,
    concurvity_max_worst = concurvity_max_worst,
    concurvity_ok = concurvity_ok,
    resid_morans_i = resid_morans_i,
    resid_morans_p = resid_morans_p,
    resid_spatial_ok = resid_spatial_ok,
    diagnostics_ok = diagnostics_ok
  )
}

.build_model_specs <- function() {
  specs <- list()
  idx <- 0L

  for (set_name in names(candidate_predictor_sets)) {
    set_cfg <- candidate_predictor_sets[[set_name]]
    fixed_lin <- sort(unique(set_cfg$fixed_linear_terms))
    nonlin_pool <- sort(unique(set_cfg$nonlinear_candidates))
    smooth_subsets <- .all_subsets(nonlin_pool)

    for (smooth_terms in smooth_subsets) {
      smooth_terms <- sort(unique(smooth_terms))
      linear_terms <- sort(unique(c(fixed_lin, setdiff(nonlin_pool, smooth_terms))))

      if (length(smooth_terms) + length(linear_terms) == 0) next

      for (k_val in k_values) {
        for (spatial_val in spatial_values) {
          ksp_values <- if (isTRUE(spatial_val)) k_spatial_values else NA_integer_
          for (ksp in ksp_values) {
            idx <- idx + 1L
            spec_id <- .make_spec_id(
              set_name,
              smooth_terms = smooth_terms,
              linear_terms = linear_terms,
              k = k_val,
              spatial = spatial_val,
              k_spatial = ksp
            )

            specs[[idx]] <- list(
              spec_index = idx,
              spec_id = spec_id,
              candidate_name = set_name,
              smooth_terms = smooth_terms,
              linear_terms = linear_terms,
              k = as.integer(k_val),
              spatial = isTRUE(spatial_val),
              k_spatial = ifelse(isTRUE(spatial_val), as.integer(ksp), NA_integer_)
            )
          }
        }
      }
    }
  }

  specs
}

model_specs <- .build_model_specs()
lichen_message("Candidate predictor sets: ", length(candidate_predictor_sets))
lichen_message("Total model specs: ", length(model_specs))


# =============================================================================
# SECTION 4 – Tuning execution
# =============================================================================

.tune_response <- function(resp, family_fn) {
  tuning_rows <- vector("list", length(model_specs))
  model_store <- list()

  for (i in seq_along(model_specs)) {
    spec <- model_specs[[i]]
    model_key <- sprintf("spec_%05d", i)

    fit_result <- tryCatch({
      m <- fit_gam_model(
        data         = modeling_data,
        response     = resp,
        linear_terms = spec$linear_terms,
        smooth_terms = spec$smooth_terms,
        family       = family_fn(),
        spatial      = spec$spatial,
        k            = spec$k,
        k_spatial    = spec$k_spatial %||% 30L
      )
      sm  <- summarize_gam_fit(m, model_name = resp)
      diag <- .extract_diag_metrics(m, modeling_data)
      list(model = m, summary = sm, diag = diag, error = NULL)
    }, error = function(e) {
      lichen_warning("FAILED ", resp, " / ", spec$spec_id, ": ", conditionMessage(e))
      list(model = NULL, summary = NULL, diag = NULL, error = conditionMessage(e))
    })

    if (!is.null(fit_result$model)) {
      model_store[[model_key]] <- fit_result$model
    }

    sm <- fit_result$summary
    diag <- fit_result$diag %||% tibble::tibble(
      k_index_min = NA_real_, k_pvalue_min = NA_real_, k_ok = NA,
      concurvity_max_worst = NA_real_, concurvity_ok = NA,
      resid_morans_i = NA_real_, resid_morans_p = NA_real_,
      resid_spatial_ok = NA, diagnostics_ok = FALSE
    )

    tuning_rows[[i]] <- tibble::tibble(
      response          = resp,
      spec_index        = spec$spec_index,
      spec_id           = spec$spec_id,
      candidate_name    = spec$candidate_name,
      k                 = spec$k,
      spatial           = spec$spatial,
      k_spatial         = spec$k_spatial,
      smooth_terms_used = .term_string(spec$smooth_terms),
      linear_terms_used = .term_string(spec$linear_terms),
      AIC               = if (!is.null(sm)) sm$AIC           else NA_real_,
      dev_explained     = if (!is.null(sm)) sm$dev_explained else NA_real_,
      r_sq              = if (!is.null(sm)) sm$r_sq          else NA_real_,
      n                 = if (!is.null(sm)) sm$n             else NA_integer_,
      family            = if (!is.null(sm)) sm$family        else NA_character_,
      link              = if (!is.null(sm)) sm$link          else NA_character_,
      fit_error         = fit_result$error %||% NA_character_
    ) |>
      dplyr::bind_cols(diag)
  }

  resp_tbl <- dplyr::bind_rows(tuning_rows)
  valid_tbl <- resp_tbl |> dplyr::filter(!is.na(AIC))
  min_aic <- if (nrow(valid_tbl) > 0) min(valid_tbl$AIC, na.rm = TRUE) else NA_real_

  resp_tbl <- resp_tbl |>
    dplyr::mutate(
      delta_AIC = AIC - min_aic,
      diagnostics_rank_penalty = dplyr::if_else(diagnostics_ok %in% TRUE, 0L, 1L, missing = 1L)
    ) |>
    dplyr::arrange(diagnostics_rank_penalty, AIC)

  list(table = resp_tbl, models = model_store)
}

lichen_banner("GAM Tuning: Full model specification grid")

.run_group <- function(responses, family_fn, family_label) {
  group_tables <- list()
  for (resp in responses) {
    lichen_banner(paste("Tuning:", resp, "|", family_label))

    result <- .tune_response(resp, family_fn)
    resp_tbl <- result$table
    models <- result$models

    readr::write_csv(resp_tbl, file.path(OUT_DIR_TUNING_REPORTS, paste0("tuning_", resp, ".csv")))
    lichen_message("Saved tuning table: tuning_", resp, ".csv")

    top_rows <- resp_tbl |>
      dplyr::filter(!is.na(AIC)) |>
      dplyr::arrange(diagnostics_rank_penalty, AIC) |>
      dplyr::slice_head(n = TOP_N)

    for (i in seq_len(nrow(top_rows))) {
      mdl_key <- sprintf("spec_%05d", top_rows$spec_index[i])
      mdl_obj <- models[[mdl_key]]
      if (!is.null(mdl_obj)) {
        rds_file <- file.path(
          OUT_DIR_TUNING_MODELS,
          sprintf("tuning_%s_rank%02d_spec%05d.rds", resp, i, top_rows$spec_index[i])
        )
        saveRDS(mdl_obj, rds_file)
      }
    }

    group_tables[[resp]] <- resp_tbl
  }
  group_tables
}

all_tuning_tables <- c(
  .run_group(responses_bin, function() stats::binomial(link = "logit"), "binomial"),
  .run_group(responses_cnt, function() mgcv::nb(link = "log"), "nb")
)


# =============================================================================
# SECTION 5 – Cross-response comparison and recommendations
# =============================================================================

lichen_banner("Cross-Response Full-Spec Summary")

all_tbl <- dplyr::bind_rows(all_tuning_tables) |>
  dplyr::group_by(response) |>
  dplyr::mutate(
    rank_by_aic = rank(AIC, ties.method = "first", na.last = "keep"),
    rank_diag_then_aic = rank(diagnostics_rank_penalty * 1e9 + AIC, ties.method = "first", na.last = "keep"),
    within_delta2 = !is.na(delta_AIC) & delta_AIC <= DELTA_AIC_TOL
  ) |>
  dplyr::ungroup()

readr::write_csv(all_tbl, file.path(OUT_DIR_TUNING_REPORTS, "all_model_specs_long.csv"))

delta_wide <- all_tbl |>
  dplyr::select(response, spec_id, delta_AIC) |>
  tidyr::pivot_wider(names_from = spec_id, values_from = delta_AIC)
readr::write_csv(delta_wide, file.path(OUT_DIR_TUNING_REPORTS, "response_by_spec_deltaAIC.csv"))

best_spec_per_response <- all_tbl |>
  dplyr::filter(!is.na(AIC)) |>
  dplyr::group_by(response) |>
  dplyr::arrange(diagnostics_rank_penalty, AIC, .by_group = TRUE) |>
  dplyr::slice_head(n = 1) |>
  dplyr::ungroup()
readr::write_csv(best_spec_per_response, file.path(OUT_DIR_TUNING_REPORTS, "best_spec_per_response.csv"))

consensus_spec_ranking <- all_tbl |>
  dplyr::filter(!is.na(AIC)) |>
  dplyr::group_by(spec_id, candidate_name, smooth_terms_used, linear_terms_used, k, spatial, k_spatial) |>
  dplyr::summarise(
    n_responses_evaluated = dplyr::n_distinct(response),
    n_best = sum(rank_by_aic == 1, na.rm = TRUE),
    n_best_diag_first = sum(rank_diag_then_aic == 1, na.rm = TRUE),
    n_within_delta2 = sum(within_delta2, na.rm = TRUE),
    mean_delta_AIC = mean(delta_AIC, na.rm = TRUE),
    median_delta_AIC = stats::median(delta_AIC, na.rm = TRUE),
    max_delta_AIC = max(delta_AIC, na.rm = TRUE),
    n_diag_ok = sum(diagnostics_ok, na.rm = TRUE),
    responses_best = paste(sort(unique(response[rank_diag_then_aic == 1])), collapse = "; "),
    responses_within_delta2 = paste(sort(unique(response[within_delta2])), collapse = "; "),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(n_best_diag_first), dplyr::desc(n_within_delta2), mean_delta_AIC)
readr::write_csv(consensus_spec_ranking, file.path(OUT_DIR_TUNING_REPORTS, "consensus_spec_ranking.csv"))

recommended_spec_groups <- best_spec_per_response |>
  dplyr::group_by(spec_id, candidate_name, smooth_terms_used, linear_terms_used, k, spatial, k_spatial) |>
  dplyr::summarise(
    n_responses = dplyr::n(),
    responses = paste(sort(unique(response)), collapse = "; "),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(n_responses))
readr::write_csv(recommended_spec_groups, file.path(OUT_DIR_TUNING_REPORTS, "recommended_spec_groups.csv"))

selected_spec_manifest <- best_spec_per_response |>
  dplyr::transmute(
    response,
    family,
    link,
    spec_id,
    candidate_name,
    smooth_terms_used,
    linear_terms_used,
    k,
    spatial,
    k_spatial,
    diagnostics_ok,
    delta_AIC,
    selected_at = as.character(Sys.time())
  ) |>
  dplyr::arrange(response)
readr::write_csv(selected_spec_manifest, file.path(OUT_DIR_TUNING_REPORTS, "selected_spec_manifest.csv"))


# =============================================================================
# SECTION 6 – Console summary
# =============================================================================

lichen_banner("Tuning Complete – Summary")

cat(sprintf("Responses tuned          : %d binary + %d count\n",
            length(responses_bin), length(responses_cnt)))
cat(sprintf("Candidate sets           : %d\n", length(candidate_predictor_sets)))
cat(sprintf("Total model specs        : %d\n", length(model_specs)))
cat(sprintf("k values                 : %s\n", paste(k_values, collapse = ", ")))
cat(sprintf("Spatial values           : %s\n", paste(spatial_values, collapse = ", ")))
cat(sprintf("k_spatial values         : %s\n", paste(k_spatial_values, collapse = ", ")))

cat("\nBest spec per response (diag-first ranking):\n")
print(
  best_spec_per_response |>
    dplyr::select(response, spec_id, AIC, delta_AIC, diagnostics_ok, k, spatial, k_spatial),
  n = nrow(best_spec_per_response)
)

cat("\nTop consensus specs:\n")
print(
  consensus_spec_ranking |>
    dplyr::select(spec_id, n_best_diag_first, n_within_delta2, mean_delta_AIC, responses_best) |>
    dplyr::slice_head(n = 10),
  n = 10
)

cat("\n✅ 09.5_model_tuning_gam.R complete.\n")
cat("   Tuning reports : ", OUT_DIR_TUNING_REPORTS, "\n", sep = "")
cat("   Tuning models  : ", OUT_DIR_TUNING_MODELS,  "\n", sep = "")
