# =============================================================================
# 09.5_model_tuning_gam.R
# GAM predictor-set and basis-dimension tuning
# =============================================================================
#
# Stage 1 tuning: for every response variable (binary + count), fits a grid of
# candidate non-spatial GAMs that vary:
#
#   1. Candidate predictor blocks (>= 6 sets, ecologically motivated).
#      Decay dummy variables (decay2–decay5) are toggled as a block.
#
#   2. k – 1D smooth basis dimension – across c(4, 5, 7, 10).
#
# Per-response outputs
#   - Tuning table (CSV):  outputs/Sumava/reports_gam/tuning_gam/tuning_<resp>.csv
#   - Top-3 models (RDS):  outputs/Sumava/models_gam/tuning_gam/tuning_<resp>_rank<N>_*.rds
#
# Cross-response output
#   - Covariate-signature frequency table (CSV):
#       outputs/Sumava/reports_gam/tuning_gam/cross_response_signature_summary.csv
#
# Usage
#   source(here::here("code", "09.5_model_tuning_gam.R"))
# =============================================================================

library(dplyr)
library(readr)
library(tibble)

source(here::here("utils.R"))
source(here::here("code", "07_model_fitting_gam.R"))

cfg <- get_project_config()


# ---- Output directories (created only if needed) ----------------------------

OUT_DIR_MODELS_GAM  <- here::here("outputs", "Sumava", "models_gam")
OUT_DIR_REPORTS_GAM <- here::here("outputs", "Sumava", "reports_gam")

OUT_DIR_TUNING_MODELS  <- file.path(OUT_DIR_MODELS_GAM,  "tuning_gam")
OUT_DIR_TUNING_REPORTS <- file.path(OUT_DIR_REPORTS_GAM, "tuning_gam")

dir.create(OUT_DIR_TUNING_MODELS,  recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR_TUNING_REPORTS, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# SECTION 1 – Load, merge, and preprocess data
# (Steps identical to 09_run_models_gam.R to ensure modeling_data is the same)
# =============================================================================

lichen_banner("Loading and Preprocessing Data")

PATH_LICHEN    <- here::here("outputs", "Sumava", "lichen_clean.csv")
PATH_STRUCTURE <- here::here("outputs", "Sumava", "structure_clean.csv")

lichen_clean    <- readr::read_csv(PATH_LICHEN,    show_col_types = FALSE)
structure_clean <- readr::read_csv(PATH_STRUCTURE, show_col_types = FALSE)

# Verify coordinate consistency before merging
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
  cat("\u2713 Coordinates match for all joined plots (tol =", tol, ")\n")
}

# Drop duplicate coord columns from structure before join
structure_clean <- structure_clean |> dplyr::select(-X, -Y)

modeling_data <- dplyr::inner_join(lichen_clean, structure_clean, by = "plot_id")

# Validate merged data
validate_input_data(
  modeling_data,
  required_cols = c("plot_id", "X", "Y"),
  id_col        = "plot_id",
  check_coords  = TRUE
)

validate_response_variables(modeling_data)

# Impute missing predictor values (median, max 5 % missing; identical policy to 09_run)
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
# SECTION 2 – Define response groups and candidate predictor sets
# =============================================================================

# Response groups – filter to columns that actually exist in modeling_data
responses_bin <- intersect(cfg$lichen_groups$binary, names(modeling_data))
responses_cnt <- intersect(cfg$lichen_groups$count,  names(modeling_data))

if (length(responses_bin) == 0 && length(responses_cnt) == 0) {
  stop("No response columns found in modeling_data. Check utils.R config.")
}

lichen_message("Binary responses:  ", paste(responses_bin, collapse = ", "))
lichen_message("Count  responses:  ", paste(responses_cnt, collapse = ", "))

# ---- Decay block (toggled on/off as a unit) ---------------------------------
decay_block   <- c("decay2", "decay3", "decay4", "decay5")
decay_present <- intersect(decay_block, names(modeling_data))

# ---- Smooth terms base (same as 09_run_models_gam.R) -----------------------
smooth_base <- intersect(c("elevation", "volume_snags", "canopy_cover"), names(modeling_data))

# ---- Helper: intersect safely with modeling_data columns --------------------
.keep <- function(...) intersect(c(...), names(modeling_data))

# ---- Candidate predictor sets (8 ecologically motivated blocks) -------------
#
# Each element is a named list with smooth_terms and linear_terms.
# Decay dummies are included/excluded as a complete block.
#
candidate_predictor_sets <- list(

  # 1. Forest structure core (diameter, dead-wood density) — no decay, no tree comp
  structure_core = list(
    smooth_terms = smooth_base,
    linear_terms = .keep("dbh_max", "dbh_sd", "n_dead_50cm")
  ),

  # 2. Forest structure core + decay block
  structure_core_decay = list(
    smooth_terms = smooth_base,
    linear_terms = .keep("dbh_max", "dbh_sd", "n_dead_50cm", decay_present)
  ),

  # 3. Tree species composition + dead wood — no decay
  composition = list(
    smooth_terms = smooth_base,
    linear_terms = .keep("ba_spruce", "ba_beech", "n_dead_50cm")
  ),

  # 4. Tree species composition + dead wood + decay block
  composition_decay = list(
    smooth_terms = smooth_base,
    linear_terms = .keep("ba_spruce", "ba_beech", "n_dead_50cm", decay_present)
  ),

  # 5. Full linear set (structure + composition) — no decay
  full_nodecay = list(
    smooth_terms = smooth_base,
    linear_terms = .keep("dbh_max", "dbh_sd", "n_dead_50cm", "ba_spruce", "ba_beech")
  ),

  # 6. Full linear set + decay block (mirrors 09_run_models_gam.R default)
  full_decay = list(
    smooth_terms = smooth_base,
    linear_terms = .keep("dbh_max", "dbh_sd", "n_dead_50cm", "ba_spruce", "ba_beech",
                         decay_present)
  ),

  # 7. Deadwood-centric: reduced smooth set, no tree-comp, no decay
  deadwood_centric = list(
    smooth_terms = .keep("volume_snags", "canopy_cover"),
    linear_terms = .keep("dbh_max", "n_dead_50cm", "ba_spruce")
  ),

  # 8. Deadwood-centric + decay block
  deadwood_centric_decay = list(
    smooth_terms = .keep("volume_snags", "canopy_cover"),
    linear_terms = .keep("dbh_max", "n_dead_50cm", "ba_spruce", decay_present)
  )
)

# Remove any candidate where both term vectors are empty (data-driven safety)
candidate_predictor_sets <- Filter(
  function(x) length(x$smooth_terms) + length(x$linear_terms) > 0,
  candidate_predictor_sets
)

lichen_message("Candidate predictor sets: ", length(candidate_predictor_sets))

# ---- k values for 1D smooths ------------------------------------------------
k_values <- c(4L, 5L, 7L, 10L)

# ---- Top-N models to save per response --------------------------------------
TOP_N <- 3L


# =============================================================================
# SECTION 3 – Tuning helper functions
# =============================================================================

#' Build a canonical covariate-signature string for a model candidate
.make_signature <- function(smooth_terms, linear_terms) {
  paste0(
    "s(", paste(sort(smooth_terms), collapse = "+"), ")",
    "__l(", paste(sort(linear_terms), collapse = "+"), ")"
  )
}


#' Fit all candidates × k for a single response; return tuning table + model list
.tune_response <- function(resp, family_fn) {

  tuning_rows  <- vector("list", length(candidate_predictor_sets) * length(k_values))
  model_store  <- list()   # key: "candname_k<val>"
  row_idx      <- 0L

  for (cand_name in names(candidate_predictor_sets)) {
    cand <- candidate_predictor_sets[[cand_name]]

    for (k_val in k_values) {
      row_idx    <- row_idx + 1L
      model_key  <- paste0(cand_name, "_k", k_val)

      fit_result <- tryCatch({
        m   <- fit_gam_model(
          data         = modeling_data,
          response     = resp,
          linear_terms = cand$linear_terms,
          smooth_terms = cand$smooth_terms,
          family       = family_fn(),
          spatial      = FALSE,
          k            = k_val
        )
        sm  <- summarize_gam_fit(m, model_name = resp)
        list(model = m, summary = sm, error = NULL)
      }, error = function(e) {
        lichen_warning("FAILED  ", resp, " / ", cand_name, " / k=", k_val,
                       ": ", conditionMessage(e))
        list(model = NULL, summary = NULL, error = conditionMessage(e))
      })

      if (!is.null(fit_result$model)) {
        model_store[[model_key]] <- fit_result$model
      }

      sm <- fit_result$summary
      tuning_rows[[row_idx]] <- tibble::tibble(
        response          = resp,
        candidate_id      = row_idx,
        candidate_name    = cand_name,
        k                 = k_val,
        smooth_terms_used = paste(cand$smooth_terms, collapse = "; "),
        linear_terms_used = paste(cand$linear_terms,  collapse = "; "),
        spatial           = FALSE,
        AIC           = if (!is.null(sm)) sm$AIC           else NA_real_,
        dev_explained = if (!is.null(sm)) sm$dev_explained else NA_real_,
        r_sq          = if (!is.null(sm)) sm$r_sq          else NA_real_,
        n             = if (!is.null(sm)) sm$n             else NA_integer_,
        family        = if (!is.null(sm)) sm$family        else NA_character_,
        link          = if (!is.null(sm)) sm$link          else NA_character_,
        fit_error     = if (!is.null(fit_result$error)) fit_result$error else NA_character_
      )
    }
  }

  resp_tbl <- dplyr::bind_rows(tuning_rows)

  # Compute delta_AIC (0 = best-fitting candidate for this response)
  min_aic  <- min(resp_tbl$AIC, na.rm = TRUE)
  resp_tbl <- resp_tbl |>
    dplyr::mutate(delta_AIC = AIC - min_aic) |>
    dplyr::arrange(AIC)

  list(table = resp_tbl, models = model_store)
}


# =============================================================================
# SECTION 4 – Run tuning grid
# =============================================================================

lichen_banner("GAM Tuning: Candidate Predictor Sets \u00d7 k")

all_tuning_tables <- list()   # accumulate per-response tables for cross-summary

.run_group <- function(responses, family_fn, family_label) {
  for (resp in responses) {
    lichen_banner(paste("Tuning:", resp, "|", family_label))

    result <- .tune_response(resp, family_fn)
    resp_tbl <- result$table
    models   <- result$models

    # ---- Save per-response tuning table ----
    readr::write_csv(
      resp_tbl,
      file.path(OUT_DIR_TUNING_REPORTS, paste0("tuning_", resp, ".csv"))
    )
    lichen_message("Saved tuning table: tuning_", resp, ".csv")

    # ---- Save top-N model objects as RDS ----
    top_rows <- resp_tbl |>
      dplyr::filter(!is.na(AIC)) |>
      dplyr::slice_head(n = TOP_N)

    for (i in seq_len(nrow(top_rows))) {
      mdl_key <- paste0(top_rows$candidate_name[i], "_k", top_rows$k[i])
      mdl_obj <- models[[mdl_key]]
      if (!is.null(mdl_obj)) {
        rds_file <- file.path(
          OUT_DIR_TUNING_MODELS,
          sprintf("tuning_%s_rank%02d_%s_k%d.rds",
                  resp, i, top_rows$candidate_name[i], top_rows$k[i])
        )
        saveRDS(mdl_obj, rds_file)
        lichen_message("Saved model RDS: ", basename(rds_file))
      }
    }

    all_tuning_tables[[resp]] <<- resp_tbl
  }
}

.run_group(responses_bin, function() stats::binomial(link = "logit"), "binomial")
.run_group(responses_cnt, function() mgcv::nb(link = "log"),          "nb")


# =============================================================================
# SECTION 5 – Cross-response covariate-signature summary
# =============================================================================

lichen_banner("Cross-Response Covariate Signature Summary")

all_tbl <- dplyr::bind_rows(all_tuning_tables)

# Top-3 per response
top3_tbl <- all_tbl |>
  dplyr::filter(!is.na(AIC)) |>
  dplyr::group_by(response) |>
  dplyr::slice_min(order_by = AIC, n = TOP_N, with_ties = FALSE) |>
  dplyr::ungroup()

# Attach canonical signature
top3_tbl <- top3_tbl |>
  dplyr::rowwise() |>
  dplyr::mutate(
    signature = .make_signature(
      strsplit(smooth_terms_used, "; ")[[1]],
      strsplit(linear_terms_used, "; ")[[1]]
    )
  ) |>
  dplyr::ungroup()

# Summarise signature frequency and delta_AIC across responses
sig_summary <- top3_tbl |>
  dplyr::group_by(signature, smooth_terms_used, linear_terms_used) |>
  dplyr::summarise(
    n_responses_in_top3 = dplyr::n(),
    mean_delta_AIC      = mean(delta_AIC, na.rm = TRUE),
    max_delta_AIC       = max(delta_AIC, na.rm = TRUE),
    responses           = paste(sort(unique(response)), collapse = "; "),
    .groups             = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(n_responses_in_top3), mean_delta_AIC)

readr::write_csv(
  sig_summary,
  file.path(OUT_DIR_TUNING_REPORTS, "cross_response_signature_summary.csv")
)
lichen_message("Saved cross-response signature summary.")


# =============================================================================
# SECTION 6 – Console summary
# =============================================================================

lichen_banner("Tuning Complete – Summary")

cat(sprintf("Responses tuned  : %d binary + %d count\n",
            length(responses_bin), length(responses_cnt)))
cat(sprintf("Candidate sets   : %d\n", length(candidate_predictor_sets)))
cat(sprintf("k values         : %s\n", paste(k_values, collapse = ", ")))
cat(sprintf("Grid size/resp   : %d\n",
            length(candidate_predictor_sets) * length(k_values)))

cat("\nTop-3 models per response (sorted by AIC):\n")
cat(strrep("-", 72), "\n")
for (resp in names(all_tuning_tables)) {
  tbl <- all_tuning_tables[[resp]] |>
    dplyr::filter(!is.na(AIC)) |>
    dplyr::slice_head(n = TOP_N)
  cat(sprintf("\n  %s\n", resp))
  for (i in seq_len(nrow(tbl))) {
    cat(sprintf("    [%d] %-26s  k=%2d  AIC=%9.2f  \u0394AIC=%6.2f  dev=%.3f\n",
                i,
                tbl$candidate_name[i],
                tbl$k[i],
                tbl$AIC[i],
                tbl$delta_AIC[i],
                tbl$dev_explained[i]))
  }
}

cat("\nCovariate signature frequency across top-3 sets:\n")
cat(strrep("-", 72), "\n")
print(
  sig_summary |>
    dplyr::select(n_responses_in_top3, mean_delta_AIC, max_delta_AIC, signature),
  n = 20
)

cat("\n\u2705 09.5_model_tuning_gam.R complete.\n")
cat("   Tuning tables : ", OUT_DIR_TUNING_REPORTS, "\n", sep = "")
cat("   Top models    : ", OUT_DIR_TUNING_MODELS,  "\n", sep = "")
