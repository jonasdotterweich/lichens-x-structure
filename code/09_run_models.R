# =============================================================================
# 09_run_models.R
# Runner: merge canonical data, preprocess, fit models, run diagnostics, save outputs
# =============================================================================

library(dplyr)
library(readr)

source(here::here("utils.R"))
source(here::here("code", "02_data_loading.R"))
source(here::here("code", "03_data_cleaning.R"))
source(here::here("code", "04_predictor_engineering.R"))
source(here::here("code", "05_categorical_encoding.R"))
source(here::here("code", "06_lichen_processing.R"))
source(here::here("code", "07_model_fitting.R"))
source(here::here("code", "08_model_diagnostics.R"))

cfg <- get_project_config()

# ---- Paths ----
PATH_LICHEN    <- here::here("outputs", "Sumava", "lichen_clean.csv")
PATH_STRUCTURE <- here::here("outputs", "Sumava", "structure_clean.csv")
#PATH_MODELING_DATA <- here::here("outputs", "Sumava", "modeling_data_merged.csv")

OUT_DIR_MODELS  <- here::here("outputs", "Sumava", "models")
OUT_DIR_REPORTS <- here::here("outputs", "Sumava", "reports")
OUT_DIR_MODELING_DATA <- here::here("outputs", "Sumava")
dir.create(OUT_DIR_MODELS,  recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR_REPORTS, recursive = TRUE, showWarnings = FALSE)

# ---- Load + merge ----
lichen_clean    <- readr::read_csv(PATH_LICHEN, show_col_types = FALSE)
structure_clean <- readr::read_csv(PATH_STRUCTURE, show_col_types = FALSE)


# Pre-check: are coordinates identical in both canonical datasets?
coords_check <- dplyr::inner_join(
  lichen_clean  |> dplyr::select(plot_id, X, Y),
  structure_clean |> dplyr::select(plot_id, X, Y),
  by = "plot_id",
  suffix = c(".lichen", ".structure")
)

# Differences (tolerance protects against tiny floating-point/rounding issues)
tol <- 1e-6
coords_check <- coords_check |>
  dplyr::mutate(
    dx = X.lichen - X.structure,
    dy = Y.lichen - Y.structure,
    coords_match = (abs(dx) <= tol) & (abs(dy) <= tol)
  )

n_mismatch <- sum(!coords_check$coords_match, na.rm = TRUE)

if (n_mismatch > 0) {
  print(coords_check |> dplyr::filter(!coords_match), n = 50)
  stop("Coordinate mismatch between lichen_clean and structure_clean for ",
       n_mismatch, " plot(s). Fix before dropping X/Y from either table.")
} else {
  cat("✓ Coordinates match for all joined plots (tol =", tol, ")\n")
}


##to not end up with double the coordinate columns while merging:

structure_clean <- structure_clean|>
  dplyr::select(-X, -Y)

modeling_data <- dplyr::inner_join(lichen_clean, structure_clean, by = "plot_id")



# ---- Validate merged data ----
validate_input_data(
  modeling_data,
  required_cols = c("plot_id", "X", "Y"),
  id_col = "plot_id",
  check_coords = TRUE
)

# Validate response columns exist / prevalence OK
resp_report <- validate_response_variables(modeling_data)
readr::write_csv(resp_report,
                 file.path(OUT_DIR_REPORTS, "qc_response_variables.csv"))

# ---- Categorical encoding (optional but useful) ----
# (only runs if those columns exist)
#modeling_data <- encode_all_categoricals(
 # modeling_data,
 # management_col = "past_management",
 # exposure_col   = "exposure",
 # dominant_col   = "dominant_tree_species"
#)

# ---- Missingness policy (STOP if any predictor >5% missing) ----
# Impute only predictors (NOT responses) IF missingness <=5%
modeling_data <- impute_missing_values_if_needed(
  modeling_data,
  cols = dplyr::where(is.numeric) &
    !dplyr::ends_with("_presence") &
    !dplyr::matches("richness$"),
  threshold_pct = 5,
  strategy = "median",
  report_path = file.path(OUT_DIR_REPORTS, "qc_imputation_report.csv")
)


# ---- Choose predictors ----
predictor_cols <- c( "elevation", "volume_snags", "canopy_cover","dbh_max", 
                     "dbh_sd", "n_dead_50cm",
                     "ba_spruce", "ba_beech", "decay2", "decay3", "decay4", 
                     "decay5"
)
 
  
 

# ---- Scaling ---

modeling_data <- modeling_data |>
  dplyr::mutate(dplyr::across(
    dplyr::all_of(predictor_cols),
    ~ as.numeric(scale(.)),
    .names = "{.col}_scaled"
  ))


# ---- Predictors to use in models (scaled column names) ----
predictors_bin   <- paste0(predictor_cols, "_scaled")
predictors_count <- predictors_bin



# Save the merged modeling dataset (useful for reproducibility)
readr::write_csv(modeling_data,
                 file.path(OUT_DIR_MODELING_DATA, "modeling_data_merged.csv"))




# ---- Fit models ----
responses_bin  <- cfg$lichen_groups$binary
responses_cnt  <- cfg$lichen_groups$count

models <- list()

# Binary models (binomial)
for (resp in responses_bin) {
  models[[resp]] <- fit_model(
    data       = modeling_data,
    response   = resp,
    predictors = predictors_bin,
    family     = stats::binomial(link = "logit")
  )
}

# Count models (negative binomial)
for (resp in responses_cnt) {
  models[[resp]] <- fit_model(
    data       = modeling_data,
    response   = resp,
    predictors = predictors_count,
    family     = glmmTMB::nbinom2
  )
}

# Save model objects
for (nm in names(models)) {
  saveRDS(models[[nm]],
          file.path(OUT_DIR_MODELS, paste0("model_", nm, ".rds")))
}

# Save coefficient summaries
coef_tbl <- dplyr::bind_rows(lapply(names(models), function(nm) {
  extract_model_summary(models[[nm]], model_name = nm)
}))

readr::write_csv(coef_tbl,
                 file.path(OUT_DIR_MODELS, "model_coefficients_all.csv"))

# ---- Spatial autocorrelation (raw responses) ----
raw_spatial <- calculate_morans_i(
  modeling_data,
  variables = c(responses_bin, responses_cnt),
  coords_cols = c("X", "Y")
)
readr::write_csv(raw_spatial,
                 file.path(OUT_DIR_REPORTS, "qc_morans_i_raw.csv"))

# ---- Diagnostics (DHARMa etc.) ----
diag_report <- compile_diagnostic_report(
  models      = models,
  data        = modeling_data,
  raw_spatial = raw_spatial,
  output_dir  = OUT_DIR_REPORTS
)

readr::write_csv(diag_report$dharma_summary,
                 file.path(OUT_DIR_REPORTS, "qc_dharma_summary.csv"))
readr::write_csv(diag_report$spatial_comparison,
                 file.path(OUT_DIR_REPORTS, "qc_spatial_comparison.csv"))

cat("\n✅ 09_run_models.R complete.\n")
cat("   Models saved to:  ", OUT_DIR_MODELS, "\n", sep = "")
cat("   Reports saved to: ", OUT_DIR_REPORTS, "\n", sep = "")