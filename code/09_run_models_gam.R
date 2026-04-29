# =============================================================================
# 09_run_models_gam.R
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
source(here::here("code", "07_model_fitting_gam.R"))
source(here::here("code", "08_model_diagnostics_gam.R"))

cfg <- get_project_config()

# ---- Paths ----
PATH_LICHEN    <- here::here("outputs", "Sumava", "lichen_clean.csv")
PATH_STRUCTURE <- here::here("outputs", "Sumava", "structure_clean.csv")
#PATH_MODELING_DATA <- here::here("outputs", "Sumava", "modeling_data_merged.csv")

OUT_DIR_MODELS_GAM  <- here::here("outputs", "Sumava", "models_gam")
OUT_DIR_REPORTS_GAM  <- here::here("outputs", "Sumava", "reports_gam")
FIG_DIR_GAM          <- here::here("figures", "gam")
OUT_DIR_MODELING_DATA_GAM  <- here::here("outputs", "Sumava", "modeling_data_gam")
dir.create(OUT_DIR_MODELS_GAM,  recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR_REPORTS_GAM, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR_GAM,  recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR_MODELING_DATA_GAM, recursive = TRUE, showWarnings = FALSE)


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
                 file.path(OUT_DIR_REPORTS_GAM, "qc_response_variables_gam.csv"))

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
  report_path = file.path(OUT_DIR_REPORTS_GAM, "qc_imputation_report_gam.csv")
)


# ---- Choose predictors ----
predictor_cols <- c( "elevation", "volume_snags", "canopy_cover","dbh_max", 
                     "dbh_sd", "n_dead_50cm",
                     "ba_spruce", "ba_beech", "decay2", "decay3", "decay4", 
                     "decay5"
)
 
# Decide which ones are smooth vs linear
smooth_terms <- c("elevation", "volume_snags", "canopy_cover")
linear_terms <- setdiff(predictor_cols, smooth_terms)  
 

## A reduced set of predictors for testing nonfitting models:

linear_terms_core <- c("dbh_max", "n_dead_50cm", "logged", "logging_intensity")
smooth_terms_core <- smooth_terms  # keep elevation, volume_snags, canopy_cover


# ---- Scaling ---
## Not utilized in GAMs  
###########

#modeling_data <- modeling_data |>
#  dplyr::mutate(dplyr::across(
#    dplyr::all_of(predictor_cols),
#    ~ as.numeric(scale(.)),
#   .names = "{.col}_scaled"
#  ))

# ---- Predictors to use in models (scaled column names) ----
#predictors_bin   <- paste0(predictor_cols, "_scaled")
#predictors_count <- predictors_bin



# Save the merged modeling dataset (useful for reproducibility)
#readr::write_csv(modeling_data,
#                 file.path(OUT_DIR_MODELING_DATA, "modeling_data_merged.csv"))

#######


# ---- Fit models ----
responses_bin  <- cfg$lichen_groups$binary
responses_cnt  <- cfg$lichen_groups$count

# For some responses, we want to fit a spatial model (with s(X, Y)) in addition to the non-spatial version, to check for spatial autocorrelation in residuals and robustness of predictors. List those responses here:
responses_spatial_control <- c("ochrolechia_presence", "core_ogf_presence", #mycoblastus_presence",
                              "calicioids_richness", "parmelia_agg_richness")

models <- list()

# Binary models (binomial)
for (resp in responses_bin) {
  
  lin <- linear_terms
  sm  <- smooth_terms
  
 #as of now mycoblastus_presence is counted out
 # if (resp == "mycoblastus_presence") {    # this is to use a reduced predictor on lichen groups were the fitting of the model failed
 #   lin <- linear_terms_core
 #   sm  <- smooth_terms_core
  }
  
  models[[resp]] <- fit_gam_model(
    data = modeling_data,
    response = resp,
    linear_terms = lin,
    smooth_terms = sm,
    family = stats::binomial(link = "logit"),
    spatial = FALSE
  )
  
  run_gam_diagnostics(
    model = models[[resp]],
    model_name = resp,
    data = modeling_data,
    out_dir_reports = OUT_DIR_REPORTS_GAM,
    coords_cols = c("X", "Y"),
    run_morans_i = TRUE
  )
  
  
  if (resp %in% responses_spatial_control) {
    nm_sp <- paste0(resp, "_sxy")
    
    ksp <- 30                   
    if (resp == "core_ogf_presence") ksp <- 15
    #if (resp == "mycoblastus_presence") ksp <- 80
    
    models[[nm_sp]] <- fit_gam_model(
      data = modeling_data,
      response = resp,
      linear_terms = lin,
      smooth_terms = sm,
      family = stats::binomial(link = "logit"),
      spatial = TRUE,
      coords_cols = c("X", "Y"),
      k_spatial = ksp
    )
    
    run_gam_diagnostics(
      model = models[[nm_sp]],
      model_name = nm_sp,
      data = modeling_data,
      out_dir_reports = OUT_DIR_REPORTS_GAM,
      coords_cols = c("X", "Y"),
      run_morans_i = TRUE
    )
  }
  
  



# Count models (negative binomial)
for (resp in responses_cnt) {
  models[[resp]] <- fit_gam_model(
    data = modeling_data,
    response = resp,
    linear_terms = linear_terms,
    smooth_terms = smooth_terms,
    family = mgcv::nb(link = "log"),
    spatial = FALSE
  )
  
  run_gam_diagnostics(
    model = models[[resp]],
    model_name = resp,
    data = modeling_data,
    out_dir_reports = OUT_DIR_REPORTS_GAM,
    coords_cols = c("X", "Y"),
    run_morans_i = TRUE
  )
  
}

# Save model objects
for (nm in names(models)) {
  saveRDS(models[[nm]],
          file.path(OUT_DIR_MODELS_GAM, paste0("model_gam_", nm, ".rds")))
}

# ---- Save coefficient summaries ---- 
param_tbl <- dplyr::bind_rows(lapply(names(models), function(nm) {
  extract_gam_parametric_summary(models[[nm]], model_name = nm)
}))

smooth_tbl <- dplyr::bind_rows(lapply(names(models), function(nm) {
  extract_gam_smooth_summary(models[[nm]], model_name = nm)
}))

fit_tbl <- dplyr::bind_rows(lapply(names(models), function(nm) {
  summarize_gam_fit(models[[nm]], model_name = nm)
}))

readr::write_csv(param_tbl,  file.path(OUT_DIR_MODELS_GAM, "model_parametric_terms_all.csv"))
readr::write_csv(smooth_tbl, file.path(OUT_DIR_MODELS_GAM, "model_smooth_terms_all.csv"))
readr::write_csv(fit_tbl,    file.path(OUT_DIR_MODELS_GAM, "model_fit_overview_all.csv"))

# ---- Spatial autocorrelation (raw responses) ----
raw_spatial <- calculate_morans_i(
  modeling_data,
  variables = c(responses_bin, responses_cnt),
  coords_cols = c("X", "Y")
)
readr::write_csv(raw_spatial,
                 file.path(OUT_DIR_REPORTS_GAM, "qc_morans_i_raw_GAM.csv"))



cat("\n✅ 09_run_models_gam.R complete.\n")
cat("   Models saved to:  ", OUT_DIR_MODELS_GAM, "\n", sep = "")
cat("   Reports saved to: ", OUT_DIR_REPORTS_GAM, "\n", sep = "")