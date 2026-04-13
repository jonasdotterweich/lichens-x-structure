# ==============================================================================
# 01b_plot_prep.R
# ==============================================================================
# PURPOSE  : Load, translate, clean, and export structural predictor data for
#            the Sumava NP case study.
# OUTPUT   : outputs/Sumava/structure_clean.csv
#              Canonical schema: plot_id | X | Y | <predictors...>
# AUTHOR   : Jonas Dotterweich
# STUDY    : Sumava NP, Czech Republic (120 plots)
# NOTE     : All study-specific decisions (Czech translations, variable
#            rename map, decay stage groupings) live inside this file.
#            They are NOT in a shared config for the same reason as 01a:
#            every new dataset will have different column names and categories.
# DEPENDS  : 01a_lichen_prep.R must have been run first
#            (lichen_clean.csv is used for the ID cross-check in Step 13).
# ==============================================================================


# ==============================================================================
# 0. SETUP
# ==============================================================================

library(tidyverse)
library(readxl)
library(corrplot)
library(naniar)
library(car)

source(here::here("utils.R"))

# Pull paths and settings from the shared config
cfg <- get_project_config()

PATH_STRUCTURE_RAW <- cfg$data$structural
PATH_COORDS_RAW    <- cfg$data$coordinates
PATH_LICHEN_CLEAN  <- here::here("outputs", "Sumava", "lichen_clean.csv")
PATH_OUT_CLEAN     <- here::here("outputs", "Sumava", "structure_clean.csv")
PATH_OUT_REPORTS   <- here::here("outputs", "Sumava", "reports")

N_PLOTS <- cfg$data$n_plots

dir.create(PATH_OUT_REPORTS, recursive = TRUE, showWarnings = FALSE)

# Collinearity thresholds (these ARE study-independent statistical conventions,
# but kept here as named constants for easy review / override per study)

VIF_MODERATE <- cfg$modeling$vif_threshold / 2   # → 5
VIF_HIGH     <- cfg$modeling$vif_threshold        # → 10
COR_HIGH     <- cfg$modeling$cor_threshold        # → 0.7


# ==============================================================================
# 1. LOAD RAW STRUCTURAL DATA
# ==============================================================================

structure_raw <- read_xlsx(PATH_STRUCTURE_RAW)

# Remove whitespace / hidden characters from column names (common in xlsx exports)
colnames(structure_raw) <- trimws(gsub("[\r\n]", "", colnames(structure_raw)))

cat("Raw structure data:", nrow(structure_raw), "plots x",
    ncol(structure_raw), "columns\n")


# ==============================================================================
# 2. TRANSLATE CZECH LABELS
# ------------------------------------------------------------------------------
# All translation maps are study-specific and defined here.
# The pattern (named character vector + recode) is reusable; the vocabulary is not.
# ==============================================================================

# --- 2a. Aspect / exposure ---------------------------------------------------
EXPOSURE_MAP <- c(
  "JV - jihovychodni"             = "SE - southeast",
  "J - jizni"                     = "S - south",
  "V - vychodni"                  = "E - east",
  "S - severni"                   = "N - north",
  "SV - severovychodni"           = "NE - northeast",
  "nehodnoceno (rovina do 5 st.)" = "not evaluated (flat land up to 5 deg.)",
  "SZ - severozapadni"            = "NW - northwest",
  "Z - zapadni"                   = "W - west",
  "JZ - jihozapadni"              = "SW - southwest"
)

# --- 2b. Coverage classes (shared across all E1 / E2 coverage columns) -------
COVERAGE_MAP <- c(
  "nevyskytuje se"                 = "absent",
  "bez vyskytu"                    = "absent",
  "ojedinely vyskyt"               = "isolated occurrence",
  "velmi ridky vyskyt (do 0.2%)"   = "very sparse occurrence (up to 0.2%)",
  "ridky vyskyt (0.2-1%)"          = "sparse occurrence (0.2-1%)",
  "ridky vyslyt (0.2-1%)"          = "sparse occurrence (0.2-1%)",
  "ridky vyskyt (od 0.2 do 1%)"    = "sparse occurrence (from 0.2 to 1%)",
  "malocetny vyskyt (1-5%)"        = "infrequent occurrence (1-5%)",
  "malocetny vyskyt (od 1-5%)"     = "infrequent occurrence (from 1-5%)",
  "hojny vyskyt (6-25%)"           = "abundant occurrence (6-25%)",
  "velmi hojny vyskyt (25-50%)"    = "very common occurrence (25-50%)",
  "velmi hojny vyskyt (26-50%)"    = "very common occurrence (26-50%)",
  "velkoplosny vyskyt (51-75%)"    = "large-scale occurrence (51-75%)",
  "dominantni vyskyt (76-100%)"    = "dominant occurrence (76-100%)"
)

# --- 2c. Management history --------------------------------------------------
MANAGEMENT_MAP <- c(
  "bez tezby, zadne parezy"             = "no logging, no stumps",
  "tezba s ponechanim vetsiny hmoty"    = "logging with most biomass left in place",
  "tezba bez ponechani hmoty"           = "logging without leaving biomass",
  "tezba s ponechanim mensiny hmoty"    = "logging with minority of biomass left in place"
)

# --- 2d. Plot surface appearance ---------------------------------------------
APPEARANCE_MAP <- c(
  "zivy les"                        = "living forest",
  "zadna z predchozich moznosti"    = "none of the previous options",
  "niva potoka/reky"                = "stream/river floodplain",
  "velkoplosne vyvraty a polomy"    = "large-scale uprooting and windthrow",
  "kurovec"                         = "bark beetle",
  "porostni mezera"                 = "canopy gap",
  "redina po pastve"                = "pasture woodland",
  "holina po tezbe"                 = "clear-cut area"
)

# Apply all translations
coverage_cols <- c(
  "coverage E1 layer - woody regeneration [%]",
  "coverage E1 layer - shrubs and subshrubs [%]",
  "coverage E1 layer - shrubs [%]",
  "coverage E1 layer - grasses [%]",
  "coverage E1 layer - herbs [%]",
  "coverage E1 layer - total vegetation [%]",
  "coverage - dwarf pine [%]",
  "coverage E2 layer - woody plants (1. 3-5m) [%]"
)

structure <- structure_raw %>%
  mutate(
    exposure = dplyr::recode(exposure, !!!EXPOSURE_MAP),
    `Past management (logging history)` =
      dplyr::recode(`Past management (logging history)`, !!!MANAGEMENT_MAP),
    `Plot surface appearance`  =
      dplyr::recode(`Plot surface appearance`, !!!APPEARANCE_MAP),
    across(
      any_of(coverage_cols),
      ~ dplyr::recode(., !!!COVERAGE_MAP)
    )
  )

cat("Translations applied: exposure, coverage classes, management, appearance\n")


# ==============================================================================
# 3. SELECT AND RENAME TO CANONICAL NAMES
# ------------------------------------------------------------------------------
# This rename map is the single point where raw column names (verbose, in Czech
# or mixed language) are mapped to the canonical snake_case names required by
# the modeling pipeline. Edit this map when column names change between studies.
# ==============================================================================

structure_clean <- structure %>%
  select(
    plot_id = Project_ID,
    
    # --- Deadwood ---
    deadwood_total     = `volume of decaying wood (logs+stumps+snags) [m3/ha]`,
    decay1             = `volume of decaying wood - decay stage 1 [m3/ha]`,
    decay2             = `volume of decaying wood - decay stage 2 [m3/ha]`,
    decay3             = `volume of decaying wood - decay stage 3 [m3/ha]`,
    decay4             = `volume of decaying wood - decay stage 4 [m3/ha]`,
    decay5             = `volume of decaying wood - decay stage 5 [m3/ha]`,
    volume_snags       = `volume of standing deadwood [m3/ha]`,
    volume_lying_logs  = `volume of lying logs [m3/ha]`,
    stump_volume       = `stump volume [m3/ha]`,
    
    # --- Tree size (old-growth proxies) ---
    dbh_mean            = `mean DBH [mm] (live+dead)`,
    dbh_median          = `median DBH [mm] (live+dead)`,
    dbh_max             = `max_DBH [mm] (live+dead)`,
    dbh_sd              = `sd_DBH (live+dead)`,
    n_dead_50cm         = `number of dead trees with DBH>50cm`,
    n_living_80cm       = `number of living with DBH>80cm`,
    
    # --- Canopy ---
    canopy_cover        = `canopy cover of main tree layer [%]`,
    total_cover         = `total canopy coverage [%]`,
    canopy_e3           = `canopy cover E3 layer estimate [%]`,
    
    # --- Stand characteristics ---
    tree_height_median  = `median height of main canopy (live+dead trees) [m]`,
    regeneration        = `regeneration density [pcs/ha]`,
    
    # --- Tree species composition ---
    ba_spruce           = `% basal area - spruce (live)`,
    ba_beech            = `% basal area - beech (live)`,
    ba_fir              = `% basal area - fir (live)`,
    ba_late_succ        = `% basal area - late successional spp (maple,sycamore,oak,ash,elm,linden) (live)`,
    ba_early_succ       = `% basal area - early successional spp (birch,aspen,poplar,pine,alder,rowan) (live)`,
    dominant_tree_species = `dominant tree species (by basal area)`,
    
    # --- Site descriptors ---
    past_management     = `Past management (logging history)`,
    exposure            = exposure,
    elevation           = `Elevation (m above sea level)`,
    habitat_type        = `Habitat type (Natura 2000)`,
    forest_type         = `Forest type`
  )

cat("\n\u2713 Selected", ncol(structure_clean) - 1, "variables\n")


# ==============================================================================
# 4. MISSING VALUE ASSESSMENT
# ==============================================================================

cat("\n=== MISSING VALUE REPORT ===\n")

missing_summary <- structure_clean %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
  mutate(
    pct_missing = round(n_missing / N_PLOTS * 100, 1),
    status = case_when(
      n_missing == 0      ~ "\u2705 Complete",
      pct_missing <  5    ~ "\u26a0\ufe0f  <5% missing",
      pct_missing < 20    ~ "\U0001f534 5-20% missing",
      TRUE                ~ "\u274c >20% missing"
    )
  ) %>%
  arrange(desc(n_missing))

print(missing_summary, n = 40)

gg_miss_var(structure_clean, show_pct = TRUE) +
  labs(title    = "Missing data by variable",
       subtitle = paste0("Forest structural predictors (n = ", N_PLOTS, " plots)"))

write_csv(missing_summary,
          file.path(PATH_OUT_REPORTS, "qc_missing_values.csv"))
cat("  \u2713 Saved: qc_missing_values.csv\n")


# ==============================================================================
# 5. DEADWOOD COMPOSITE VARIABLES
# ------------------------------------------------------------------------------
# Which decay stages to group is an ecological decision made per study, informed
# by the correlation analysis in Step 6 (and originally tested in the old
# 01b exploratory section). The groupings below reflect the Bohubin study
# conclusion. A new study should re-run Step 6 before deciding.
# ==============================================================================

structure_clean <- structure_clean %>%
  mutate(
    deadwood_decay2to5   = decay2 + decay3 + decay4 + decay5,
    deadwood_decay4to5   = decay4 + decay5,
    pct_decay_advanced   = if_else(deadwood_total > 0,
                                   (deadwood_decay4to5 / deadwood_total) * 100,
                                   0)
  )

cat("\n\u2713 Composite deadwood variables created:",
    "deadwood_decay2to5, deadwood_decay4to5, pct_decay_advanced\n")


# ==============================================================================
# 6. DATA-DRIVEN DECAY STAGE CORRELATION CHECK
# ------------------------------------------------------------------------------
# Retained from original script — guides the grouping decision above.
# Reads lichen_clean.csv produced by 01a to perform the cross-variable test.
# ==============================================================================

cat("\n=== DECAY STAGE CORRELATION ANALYSIS ===\n")

lichen_clean_for_check <- tryCatch(
  read_csv(PATH_LICHEN_CLEAN, show_col_types = FALSE),
  error = function(e) {
    warning("Could not load lichen_clean.csv — skipping decay correlation check.\n",
            "Run 01a_lichen_prep.R first.")
    NULL
  }
)

if (!is.null(lichen_clean_for_check)) {
  
  decay_check_data <- lichen_clean_for_check %>%
    inner_join(structure_clean, by = "plot_id")
  
  decay_correlations <- decay_check_data %>%
    summarise(
      r_decay1 = cor(decay1, calicioids_richness, use = "complete.obs"),
      r_decay2 = cor(decay2, calicioids_richness, use = "complete.obs"),
      r_decay3 = cor(decay3, calicioids_richness, use = "complete.obs"),
      r_decay4 = cor(decay4, calicioids_richness, use = "complete.obs"),
      r_decay5 = cor(decay5, calicioids_richness, use = "complete.obs")
    ) %>%
    pivot_longer(everything(),
                 names_to  = "decay_stage",
                 values_to = "r_calicioids") %>%
    mutate(
      abs_r    = abs(r_calicioids),
      strength = case_when(
        abs_r > 0.5 ~ "\U0001f534 Strong",
        abs_r > 0.3 ~ "\u26a0\ufe0f  Moderate",
        abs_r > 0.1 ~ "\u2705 Weak",
        TRUE        ~ "\u26aa Negligible"
      )
    )
  
  cat("Pearson r with calicioids_richness:\n")
  print(decay_correlations, n = Inf)
  
  write_csv(decay_correlations,
            file.path(PATH_OUT_REPORTS, "qc_decay_correlations.csv"))
  cat("  \u2713 Saved: qc_decay_correlations.csv\n")
}


# ==============================================================================
# 7. REMOVE PLOTS WITH EXCESSIVE MISSING KEY PREDICTORS
# ==============================================================================

KEY_PREDICTORS <- c("dbh_max", "deadwood_total", "canopy_cover",
                    "ba_spruce", "ba_beech", "tree_height_median")
NA_KEY_MAX     <- 2   # Maximum allowed missing KEY predictors per plot

plots_to_drop <- structure_clean %>%
  mutate(n_na_key = rowSums(is.na(select(., all_of(KEY_PREDICTORS))))) %>%
  filter(n_na_key > NA_KEY_MAX) %>%
  pull(plot_id)

if (length(plots_to_drop) > 0) {
  cat("\n\u26a0\ufe0f  Dropping", length(plots_to_drop),
      "plot(s) with >", NA_KEY_MAX, "missing key predictors:\n")
  print(plots_to_drop)
  structure_clean <- filter(structure_clean, !plot_id %in% plots_to_drop)
}

# Impute remaining NAs with column median
structure_clean <- structure_clean %>%
  mutate(across(where(is.numeric),
                ~ if_else(is.na(.), median(., na.rm = TRUE), .)))

cat("\u2713 After NA handling:", nrow(structure_clean), "plots remaining\n")


# ==============================================================================
# 8. ENCODE CATEGORICAL VARIABLES
# ==============================================================================

# Management — ordered factor + numeric proxy
management_levels <- c(
  "no logging, no stumps",
  "logging with most biomass left in place",
  "logging with minority of biomass left in place",
  "logging without leaving biomass"
)

structure_clean <- structure_clean %>%
  mutate(
    management_factor = factor(past_management,
                               levels  = management_levels,
                               ordered = TRUE),
    logged            = if_else(past_management == "no logging, no stumps", 0L, 1L),
    logging_intensity = as.integer(management_factor) - 1L
  )

cat("\nManagement categories:\n")
print(table(structure_clean$past_management))

# Exposure — simplify to 8-direction + flat
structure_clean <- structure_clean %>%
  mutate(
    exposure_simple = case_when(
      str_detect(exposure, "not evaluated") ~ "flat",
      str_detect(exposure, "^N ")           ~ "N",
      str_detect(exposure, "^S ")           ~ "S",
      str_detect(exposure, "^E ")           ~ "E",
      str_detect(exposure, "^W ")           ~ "W",
      str_detect(exposure, "NE|SE")         ~ "NE_SE",
      str_detect(exposure, "NW|SW")         ~ "NW_SW",
      TRUE                                  ~ "other"
    ),
    exposure_factor = factor(exposure_simple)
  )

cat("\nExposure categories:\n")
print(table(structure_clean$exposure_simple))

# Dominant tree species — keep SM (spruce), BK (beech), JD (fir), other
structure_clean <- structure_clean %>%
  mutate(
    dominant_grouped = case_when(
      dominant_tree_species %in% c("SM", "BK", "JD") ~ dominant_tree_species,
      TRUE ~ "other"
    ),
    dominant_factor = factor(dominant_grouped,
                             levels = c("SM", "BK", "JD", "other"))
  )

cat("\nDominant species:\n")
print(table(structure_clean$dominant_grouped))

cat("\n\u2713 Categorical variables encoded\n")


# ==============================================================================
# 9. COLLINEARITY CHECK
# ==============================================================================

numeric_for_vif <- structure_clean %>%
  select(
    deadwood_total, deadwood_decay2to5, deadwood_decay4to5,
    dbh_mean, dbh_median, dbh_max, dbh_sd,
    canopy_cover, total_cover,
    tree_height_median,
    ba_spruce, ba_beech, ba_late_succ, ba_early_succ,
    n_dead_50cm, n_living_80cm,
    logging_intensity,
    elevation
  ) %>%
  drop_na()

# Correlation matrix
cor_matrix <- cor(numeric_for_vif, use = "complete.obs")

png(file.path(PATH_OUT_REPORTS, "qc_correlation_matrix.png"),
    width = 1000, height = 1000)
corrplot(cor_matrix, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, tl.cex = 0.8,
         title = "Correlation Matrix — Structural Predictors",
         mar   = c(0, 0, 2, 0))
dev.off()
cat("\n\u2713 Saved: qc_correlation_matrix.png\n")

high_cor_pairs <- which(abs(cor_matrix) > COR_HIGH & cor_matrix != 1,
                        arr.ind = TRUE)
if (nrow(high_cor_pairs) > 0) {
  cat("\u26a0\ufe0f  Highly correlated pairs (|r| >", COR_HIGH, "):\n")
  data.frame(
    var1        = rownames(cor_matrix)[high_cor_pairs[, 1]],
    var2        = colnames(cor_matrix)[high_cor_pairs[, 2]],
    correlation = round(cor_matrix[high_cor_pairs], 3)
  ) %>%
    filter(var1 < var2) %>%
    arrange(desc(abs(correlation))) %>%
    print(row.names = FALSE)
}


# ==============================================================================
# 10. VIF CALCULATION
# ==============================================================================

cat("\n=== VIF ANALYSIS ===\n")

vif_model  <- lm(deadwood_total ~ ., data = numeric_for_vif)
vif_values <- car::vif(vif_model)

vif_table <- data.frame(
  variable = names(vif_values),
  VIF      = round(vif_values, 2),
  status   = case_when(
    vif_values >= VIF_HIGH ~ "\U0001f534 High collinearity",
    vif_values >= VIF_MODERATE ~ "\u26a0\ufe0f  Moderate",
    TRUE ~ "\u2705 OK"
  ),
  stringsAsFactors = FALSE,
  row.names        = NULL
) %>% arrange(desc(VIF))

print(vif_table, row.names = FALSE)
cat(sprintf("\nMean VIF: %.2f  |  Max VIF: %.2f\n",
            mean(vif_values), max(vif_values)))

if (max(vif_values) >= VIF_HIGH) {
  warning(sum(vif_values >= VIF_HIGH),
          " predictor(s) with VIF >= ", VIF_HIGH,
          " — resolve before modeling.")
}

write_csv(vif_table, file.path(PATH_OUT_REPORTS, "qc_vif_report.csv"))
cat("  \u2713 Saved: qc_vif_report.csv\n")


# ==============================================================================
# 11. FINAL PREDICTOR SELECTION
# ------------------------------------------------------------------------------
# Variables retained here reflect VIF + ecological relevance decisions made for
# the Bohubin study. The decision is documented inline.
# tree_height_median and regeneration are dropped (high VIF / low relevance).
# ==============================================================================

structure_final <- structure_clean %>%
  select(
    plot_id,
    
    # --- Old-growth proxies ---
    dbh_max,           # Keep max; drop mean/median (collinear)
    dbh_sd,
    n_living_80cm,
    n_dead_50cm,
    
    # --- Deadwood ---
    volume_snags,
    decay2, decay3, decay4, decay5,   # Individual stages retained for 02 predictor sets
    
    # --- Canopy ---
    canopy_cover,
    
    # --- Composition ---
    ba_spruce, ba_beech, ba_late_succ,
    
    # --- Site ---
    elevation,
    
    # --- Management (encoded) ---
    logging_intensity,
    logged,
    
    # --- Categorical (for contextual models / descriptive tables) ---
    past_management,
    management_factor,
    exposure_simple,
    exposure_factor,
    dominant_tree_species,
    dominant_factor,
    
    # --- Contextual (retained in output but not in core predictor set) ---
    habitat_type,
    forest_type
  )

cat("\n\u2713 Final predictor set:", ncol(structure_final) - 1, "variables\n")


# ==============================================================================
# 12. DISTRIBUTION CHECKS
# ==============================================================================

# Zero-inflation
zero_pct <- structure_final %>%
  summarise(
    pct_zero_snags     = mean(volume_snags  == 0, na.rm = TRUE) * 100,
    pct_zero_decay4    = mean(decay4         == 0, na.rm = TRUE) * 100,
    pct_zero_decay5    = mean(decay5         == 0, na.rm = TRUE) * 100,
    pct_zero_n_dead    = mean(n_dead_50cm    == 0, na.rm = TRUE) * 100
  )
cat("\n=== ZERO-INFLATION CHECK ===\n")
print(zero_pct)

if (zero_pct$pct_zero_decay5 > 50) {
  cat("\u26a0\ufe0f  decay5 is >50% zero (", round(zero_pct$pct_zero_decay5, 1),
      "%) — consider zero-inflated model\n", sep = "")
}

# Histograms
cont_vars <- c("dbh_max", "dbh_sd", "n_dead_50cm", "n_living_80cm",
               "volume_snags", "canopy_cover",
               "ba_spruce", "ba_beech", "ba_late_succ", "elevation",
               "decay4", "decay5")

png(file.path(PATH_OUT_REPORTS, "qc_predictor_distributions.png"),
    width = 1400, height = 1000)
par(mfrow = c(3, 4))
for (v in cont_vars) {
  hist(structure_final[[v]], main = v, xlab = "", col = "lightblue", breaks = 20)
}
par(mfrow = c(1, 1))
dev.off()
cat("  \u2713 Saved: qc_predictor_distributions.png\n")


# ==============================================================================
# 13. ATTACH COORDINATES
# ==============================================================================

coords_raw <- read_excel(PATH_COORDS_RAW)
colnames(coords_raw) <- trimws(colnames(coords_raw))

coords_clean <- coords_raw %>%
  rename(plot_id = č., coordinates = souradnice) %>%
  select(plot_id, X, Y)

structure_with_coords <- structure_final %>%
  left_join(coords_clean, by = "plot_id")

missing_coord_plots <- filter(structure_with_coords, is.na(X))$plot_id
if (length(missing_coord_plots) > 0) {
  warning("These plots have no coordinates: ",
          paste(missing_coord_plots, collapse = ", "))
} else {
  cat("\n\u2713 All plots matched to coordinates\n")
}

cat("Coordinate ranges (S-JTSK):\n")
cat("  X:", range(structure_with_coords$X, na.rm = TRUE), "\n")
cat("  Y:", range(structure_with_coords$Y, na.rm = TRUE), "\n")


# ==============================================================================
# 14. ID CROSS-CHECK AGAINST lichen_clean.csv
# ------------------------------------------------------------------------------
# Ensures that both canonical outputs cover the same set of plots before saving.
# ==============================================================================

if (!is.null(lichen_clean_for_check)) {
  
  lichen_ids   <- lichen_clean_for_check$plot_id
  structure_ids <- structure_with_coords$plot_id
  
  only_in_lichen    <- setdiff(lichen_ids, structure_ids)
  only_in_structure <- setdiff(structure_ids, lichen_ids)
  
  if (length(only_in_lichen) > 0) {
    warning("Plots in lichen_clean but NOT in structure_clean:\n  ",
            paste(only_in_lichen, collapse = ", "),
            "\n  These plots will be lost at the merge step in 02.")
  }
  if (length(only_in_structure) > 0) {
    cat("\u26a0\ufe0f  Plots in structure_clean but NOT in lichen_clean (no lichen data):\n  ",
        paste(only_in_structure, collapse = ", "),
        "\n  These will be excluded from the modelling dataset.\n")
  }
  if (length(only_in_lichen) == 0 && length(only_in_structure) == 0) {
    cat("\u2713 Plot ID sets match perfectly between lichen_clean and structure_clean\n")
  }
}


# ==============================================================================
# 15. BUILD CANONICAL structure_clean.csv
# ------------------------------------------------------------------------------
# Schema (per contract):
#   plot_id | X | Y | <canonical predictors...>
# Columns NOT included here:
#   - scaled columns (added at runtime in 02_glm_setup.R, never saved)
#   - raw Czech column names (already renamed in Step 3)
#   - collinear variables dropped in Step 11 (dbh_mean, dbh_median, total_cover, etc.)
# ==============================================================================

structure_clean_final <- structure_with_coords %>%
  select(
    plot_id,
    X,
    Y,
    # Predictors in canonical names
    elevation,
    volume_snags,
    canopy_cover,
    decay2, decay3, decay4, decay5,
    dbh_max, dbh_sd,
    n_dead_50cm,
    n_living_80cm,
    ba_spruce, ba_beech, ba_late_succ,
    # Management
    logging_intensity,
    logged,
    # Categorical — English labels
    exposure            = exposure_simple,
    habitat_type,
    forest_type,
    dominant_tree_species,
    past_management
  )

# Validate: plot_id is unique
if (anyDuplicated(structure_clean_final$plot_id)) {
  stop("Duplicate plot_id values in structure_clean_final")
}

cat("\n\u2550\u2550\u2550 structure_clean SUMMARY \u2550\u2550\u2550\n")
cat("Dimensions:", nrow(structure_clean_final), "plots x",
    ncol(structure_clean_final), "columns\n")
cat("Columns   :", paste(names(structure_clean_final), collapse = ", "), "\n")

# Summary statistics for continuous predictors
cat("\n=== CONTINUOUS PREDICTOR SUMMARY ===\n")
structure_clean_final %>%
  select(where(is.numeric), -X, -Y, -plot_id) %>%
  summary() %>%
  print()

# Save canonical output
write_csv(structure_clean_final, PATH_OUT_CLEAN)
cat(sprintf("\n  \u2713 Saved: %s  [%d rows x %d cols]\n",
            PATH_OUT_CLEAN,
            nrow(structure_clean_final),
            ncol(structure_clean_final)))

cat("\n\u2705 01b_plot_prep.R complete.\n")
cat("   Both canonical outputs are ready:\n")
cat("     \u2022", PATH_LICHEN_CLEAN, "\n")
cat("     \u2022", PATH_OUT_CLEAN, "\n")
cat("   Next: run 02_data_loading.R\n")    


