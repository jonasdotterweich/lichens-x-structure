# =============================================================================
# 10B_gam_figures_sensitivity_heatmap.R
# Meeting figure: sensitivity heatmap (low->high effect size) across indicators
# =============================================================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)

MODELS_DIR <- here::here("outputs", "Sumava", "models_gam")
FIG_DIR    <- here::here("figures", "gam_meeting")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

model_files <- list.files(MODELS_DIR, pattern = "^model_gam_.*\\.rds$", full.names = TRUE)
stopifnot(length(model_files) > 0)

models <- setNames(
  lapply(model_files, readRDS),
  gsub("^model_gam_|\\.rds$", "", basename(model_files))
)

## for excluding spatial control models 

models <- models[!grepl("_sxy$", names(models))]

vars_to_test <- c("canopy_cover", "elevation", "volume_snags", "ba_beech", "ba_spruce", "n_dead_50cm")

typical_row <- function(m) {
  mf <- m$model
  out <- as.list(rep(NA, ncol(mf)))
  names(out) <- names(mf)
  for (v in names(mf)) {
    if (is.numeric(mf[[v]])) out[[v]] <- median(mf[[v]], na.rm = TRUE)
    if (is.factor(mf[[v]]))  out[[v]] <- levels(mf[[v]])[1]
  }
  as.data.frame(out)
}

sens_one <- function(m, model_name) {
  mf <- m$model
  nd <- typical_row(m)
  
  res <- lapply(vars_to_test, function(v) {
    if (!v %in% names(mf)) return(NULL)
    lo <- as.numeric(quantile(mf[[v]], 0.2, na.rm = TRUE))
    hi <- as.numeric(quantile(mf[[v]], 0.8, na.rm = TRUE))
    
    nd_lo <- nd; nd_lo[[v]] <- lo
    nd_hi <- nd; nd_hi[[v]] <- hi
    
    p_lo <- as.numeric(predict(m, newdata = nd_lo, type = "response"))
    p_hi <- as.numeric(predict(m, newdata = nd_hi, type = "response"))
    
    tibble(model = model_name, var = v, effect_20_80 = p_hi - p_lo, p20 = p_lo, p80 = p_hi)
  })
  
  bind_rows(res)
}

tbl <- bind_rows(lapply(names(models), function(nm) sens_one(models[[nm]], nm)))

# Heatmap
p <- ggplot(tbl, aes(x = var, y = model, fill = effect_20_80)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#D95F0E", mid = "white", high = "#2C7FB8", midpoint = 0) +
  labs(
    title = "Indicator sensitivity to structure covariates",
    subtitle = "Fill = predicted change moving covariate from 20th → 80th percentile (others held at median).",
    x = NULL, y = NULL, fill = "Δpred"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(FIG_DIR, "sensitivity_heatmap.png"), p, width = 9, height = 5.5, dpi = 300)
readr::write_csv(tbl, file.path(FIG_DIR, "sensitivity_heatmap_table.csv"))

cat("✅ Wrote:", file.path(FIG_DIR, "sensitivity_heatmap.png"), "\n")