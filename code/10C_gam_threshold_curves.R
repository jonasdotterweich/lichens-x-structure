# =============================================================================
# 10C_gam_threshold_curves.R
# Plot probability curves with a cutoff (e.g., p=0.2) and mark threshold crossings
# =============================================================================

library(dplyr)
library(ggplot2)
library(readr)

if (!requireNamespace("ggeffects", quietly = TRUE)) install.packages("ggeffects")
library(ggeffects)
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---- Settings ----
indicator <- "elite_rare_presence" # "calicioids_richness"   "parmelia_agg_presence" "ochrolechia_presence"  "core_ogf_presence"     "xylographa_presence"   "elite_rare_presence"
cutoff_p  <- 0.25

# covariates to draw curves for (must be in the model)
covariates <- c("canopy_cover", "volume_snags", "elevation", "ba_beech", "ba_spruce", "n_dead_50cm")

# If you want the other covariates held at something other than median, change here
hold_at <- list()  # e.g. list(elevation = 900)

# ---- Paths ----
MODELS_DIR <- here::here("outputs", "Sumava", "models_gam")
FIG_DIR    <- here::here("figures", "gam_thresholds", indicator)
MANIFEST_PATH <- here::here("outputs", "Sumava", "reports_gam", "tuning_gam", "selected_spec_manifest.csv")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Load model ----
model_path <- file.path(MODELS_DIR, paste0("model_gam_", indicator, ".rds"))
stopifnot(file.exists(model_path))
m <- readRDS(model_path)

.split_terms <- function(x) {
  x <- x %||% ""
  if (length(x) == 0 || is.na(x) || !nzchar(x)) return(character(0))
  trimws(strsplit(x, ";", fixed = TRUE)[[1]])
}

if (file.exists(MANIFEST_PATH)) {
  manifest <- readr::read_csv(MANIFEST_PATH, show_col_types = FALSE)
  if (indicator %in% manifest$response) {
    row <- manifest |>
      dplyr::filter(response == indicator) |>
      dplyr::slice(1)
    covariates <- unique(c(
      .split_terms(row$smooth_terms_used),
      .split_terms(row$linear_terms_used)
    ))
    covariates <- covariates[!(covariates %in% c("X", "Y"))]
  }
}

# Helper: find approximate crossings of y(x) with cutoff
find_crossings <- function(x, y, cutoff) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 2) return(numeric(0))
  
  # ensure x ordered
  o <- order(x)
  x <- x[o]; y <- y[o]
  
  s <- y - cutoff
  idx <- which(s[-1] * s[-length(s)] <= 0) # sign change or exact touch
  
  if (length(idx) == 0) return(numeric(0))
  
  # linear interpolation for crossing location
  x0 <- x[idx]; x1 <- x[idx + 1]
  y0 <- y[idx]; y1 <- y[idx + 1]
  xc <- x0 + (cutoff - y0) * (x1 - x0) / (y1 - y0)
  
  # drop NaN from flat segments
  xc[is.finite(xc)]
}

# ---- Make curves ----
for (v in covariates) {
  
  # Skip covariates not in model frame (works for smooth or linear)
  if (!v %in% names(m$model)) next
  
  # ggeffects curve over the observed range
  eff <- tryCatch(
    ggeffects::ggpredict(m, terms = paste0(v, " [all]")),
    error = function(e) NULL
  )
  if (is.null(eff)) next
  
  eff_df <- as.data.frame(eff)
  
  # Apply "hold_at" overrides by re-predicting with explicit condition if desired.
  # (Optional advanced: for now we keep ggeffects default which holds others at typical values.)
  
  # Find crossings using mean prediction
  xs <- find_crossings(eff_df$x, eff_df$predicted, cutoff_p)
  
  # Build plot
  p <- ggplot(eff_df, aes(x = x, y = predicted)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.20) +
    geom_line(linewidth = 0.9) +
    geom_hline(yintercept = cutoff_p, linetype = "dashed", linewidth = 0.3) +
    labs(
      title = paste0(indicator, ": threshold curve for ", v),
      subtitle = paste0("Cutoff line at p = ", cutoff_p),
      x = v,
      y = "Predicted probability"
    ) +
    theme_bw(base_size = 11)
  
  # Add vertical lines at crossings (if any)
  if (length(xs) > 0) {
    p <- p +
      geom_vline(xintercept = xs, linetype = "dotted", linewidth = 0.3, color = "grey20")
    
    # annotate with numbers (first crossing only, to keep it clean)
    p <- p +
      annotate("text",
               x = xs[1], y = cutoff_p,
               label = paste0("x≈", signif(xs[1], 3)),
               vjust = -0.8, hjust = -0.05, size = 3)
  } else {
    p <- p + labs(caption = "No crossing of cutoff within observed range.")
  }
  
  out <- file.path(FIG_DIR, paste0("threshold_", indicator, "__", v, "_p", cutoff_p, ".png"))
  ggsave(out, p, width = 7.5, height = 5, dpi = 300)
  
  # Write crossings table per covariate
  out_csv <- file.path(FIG_DIR, paste0("threshold_", indicator, "__", v, "_p", cutoff_p, ".csv"))
  readr::write_csv(
    tibble(indicator = indicator, covariate = v, cutoff = cutoff_p, crossing_x = xs),
    out_csv
  )
}

cat("\n✅ Threshold curves written to:", FIG_DIR, "\n")
