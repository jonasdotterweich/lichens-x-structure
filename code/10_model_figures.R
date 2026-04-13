# =============================================================================
# 10_model_figures.R
# Create coefficient "importance" barplots + response curves for all fitted models
# =============================================================================

library(dplyr)
library(readr)
library(ggplot2)

# Needed for marginal effect curves
if (!requireNamespace("ggeffects", quietly = TRUE)) install.packages("ggeffects")
library(ggeffects)

source(here::here("utils.R"))
source(here::here("code", "07_model_fitting.R"))  # for extract_model_summary()

# ---- Paths ----
MODELS_DIR  <- here::here("outputs", "Sumava", "models")
FIG_DIR     <- here::here("figures")  # as requested: working dir
REPORTS_DIR <- here::here("outputs", "Sumava", "reports")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(REPORTS_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Load models ----
model_files <- list.files(MODELS_DIR, pattern = "^model_.*\\.rds$", full.names = TRUE)
stopifnot(length(model_files) > 0)

models <- setNames(
  lapply(model_files, readRDS),
  gsub("^model_|\\.rds$", "", basename(model_files))
)

# ---- Extract coefficient tables (one big tibble) ----
coef_all <- dplyr::bind_rows(lapply(names(models), function(nm) {
  extract_model_summary(models[[nm]], model_name = nm)
}))

readr::write_csv(coef_all, file.path(REPORTS_DIR, "model_coefficients_all_tidy.csv"))

# Helper: drop intercept, compute CI, compute "importance" = |estimate|
prep_coef <- function(df) {
  df |>
    dplyr::filter(term != "(Intercept)") |>
    dplyr::mutate(
      conf_low  = estimate - 1.96 * std_error,
      conf_high = estimate + 1.96 * std_error,
      importance = abs(estimate),
      direction = ifelse(estimate >= 0, "positive", "negative")
    ) |>
    dplyr::arrange(dplyr::desc(importance))
}

# ---- 1) Bar plot per response: standardized coefficient importance ----
for (resp in unique(coef_all$model)) {
  df <- coef_all |>
    dplyr::filter(model == resp) |>
    prep_coef()
  
  if (nrow(df) == 0) next
  
  p <- ggplot(df, aes(x = reorder(term, importance), y = estimate, fill = direction)) +
    geom_col(width = 0.7) +
    geom_errorbar(aes(ymin = conf_low, ymax = conf_high), width = 0.2, linewidth = 0.4) +
    coord_flip() +
    scale_fill_manual(values = c(positive = "#2C7FB8", negative = "#D95F0E")) +
    labs(
      title = paste("Predictor effects (scaled):", resp),
      subtitle = "Bars = coefficient estimate; error bars = ~95% CI. Predictors are scaled so magnitudes are comparable.",
      x = NULL,
      y = "Coefficient estimate"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
  
  out <- file.path(FIG_DIR, paste0("coef_bar_", resp, ".png"))
  ggsave(out, p, width = 9, height = 5, dpi = 300)
}

# ---- 2) Response curves (marginal effects) ----
# Pick key predictors you want curves for (must match your model terms)
# Use the scaled names if that’s what your models used.
key_terms <- c(
  "dbh_max_scaled",
  "canopy_cover_scaled",
  "volume_snags_scaled",
  "elevation_scaled"
)

for (resp in names(models)) {
  m <- models[[resp]]
  
  for (term in key_terms) {
    # Skip terms not in this model
    if (!term %in% all.vars(stats::formula(m))) next
    
    # ggpredict term specification:
    # use a nice range of the scaled variable (-2..2 SD is typical)
    eff <- tryCatch(
      ggeffects::ggpredict(m, terms = c(paste0(term, " [-2:2]"))),
      error = function(e) NULL
    )
    if (is.null(eff)) next
    
    p <- plot(eff) +
      ggtitle(paste("Response curve:", resp)) +
      labs(subtitle = paste("Marginal effect of", term, "(others held constant)")) +
      theme_bw(base_size = 11)
    
    out <- file.path(FIG_DIR, paste0("curve_", resp, "__", term, ".png"))
    ggsave(out, p, width = 7.5, height = 5, dpi = 300)
  }
}

cat("\n✅ Figures written to:", FIG_DIR, "\n")
cat("✅ Tidy coefficient table written to:", REPORTS_DIR, "\n")