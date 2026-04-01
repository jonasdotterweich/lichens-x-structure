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
    exposure = recode(exposure, !!!EXPOSURE_MAP),
    `Past management (logging history)` =
      recode(`Past management (logging history)`, !!!MANAGEMENT_MAP),
    `Plot surface appearance`  =
      recode(`Plot surface appearance`, !!!APPEARANCE_MAP),
    across(
      any_of(coverage_cols),
      ~ recode(., !!!COVERAGE_MAP)
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
cat("   Next: run 02_glm_setup.R\n")                              "bez vyskytu" = "absent",
"ojedinely vyskyt" = "isolated occurrence",
"dominantni vyskyt (76-100%)" = "dominant occurrence (76-100%)",
"ridky vyslyt (0.2-1%)" = "sparse occurrence (0.2-1%)"
))



structure <- structure %>%
  mutate(`coverage E1 layer - shrubs and subshrubs [%]` = recode(`coverage E1 layer - shrubs and subshrubs [%]`,
                                                                 "nevyskytuje se" = "absent",
                                                                 "ojedinely vyskyt" = "isolated occurrence",
                                                                 "ridky vyskyt (od 0.2 do 1%)" = "sparse occurrence (from 0.2 to 1%)",
                                                                 "velmi ridky vyskyt (do 0.2%)" = "very sparse occurrence (up to 0.2%)",
                                                                 "malocetny vyskyt (od 1-5%)" = "infrequent occurrence (from 1-5%)",
                                                                 "hojny vyskyt (6-25%)" = "abundant occurrence (6-25%)"
  ))




structure <- structure %>%
  mutate(`coverage E1 layer - shrubs [%]` = recode(`coverage E1 layer - shrubs [%]`,
                                                   "hojny vyskyt (6-25%)" = "abundant occurrence (6-25%)",
                                                   "nevyskytuje se" = "absent",
                                                   "ridky vyskyt (od 0.2 do 1%)" = "sparse occurrence (from 0.2 to 1%)",
                                                   "velmi hojny vyskyt (25-50%)" = "very common occurrence (25-50%)",
                                                   "velmi ridky vyskyt (do 0.2%)" = "very sparse occurrence (up to 0.2%)",
                                                   "dominantni vyskyt (76-100%)" = "dominant occurrence (76-100%)",
                                                   "malocetny vyskyt (od 1-5%)" = "infrequent occurrence (from 1-5%)",
                                                   "velkoplosny vyskyt (51-75%)" = "large-scale occurrence (51-75%)",
                                                   "ojedinely vyskyt" = "isolated occurrence",
                                                   "velmi hojny vyskyt (26-50%)" = "very common occurrence (26-50%)"
  ))



structure <- structure %>%
  mutate(`coverage E1 layer - grasses [%]` = recode(`coverage E1 layer - grasses [%]`,
                                                    "malocetny vyskyt (od 1-5%)" = "infrequent occurrence (from 1-5%)",
                                                    "velmi ridky vyskyt (do 0.2%)" = "very sparse occurrence (up to 0.2%)",
                                                    "velmi hojny vyskyt (25-50%)" = "very common occurrence (25-50%)",
                                                    "nevyskytuje se" = "absent",
                                                    "dominantni vyskyt (76-100%)" = "dominant occurrence (76-100%)",
                                                    "ridky vyskyt (od 0.2 do 1%)" = "sparse occurrence (from 0.2 to 1%)",
                                                    "velkoplosny vyskyt (51-75%)" = "large-scale occurrence (51-75%)",
                                                    "hojny vyskyt (6-25%)" = "abundant occurrence (6-25%)",
                                                    "ojedinely vyskyt" = "isolated occurrence",
                                                    "velmi hojny vyskyt (26-50%)" = "very common occurrence (26-50%)"
  ))


structure <- structure %>%
  mutate(`coverage E1 layer - herbs [%]` = recode(`coverage E1 layer - herbs [%]`,
                                                  "velmi ridky vyskyt (do 0.2%)" = "very sparse occurrence (up to 0.2%)",
                                                  "malocetny vyskyt (od 1-5%)" = "infrequent occurrence (from 1-5%)",
                                                  "velkoplosny vyskyt (51-75%)" = "large-scale occurrence (51-75%)",
                                                  "ridky vyskyt (od 0.2 do 1%)" = "sparse occurrence (from 0.2 to 1%)",
                                                  "ojedinely vyskyt" = "isolated occurrence",
                                                  "hojny vyskyt (6-25%)" = "abundant occurrence (6-25%)",
                                                  "nevyskytuje se" = "absent",
                                                  "velmi hojny vyskyt (25-50%)" = "very common occurrence (25-50%)",
                                                  "dominantni vyskyt (76-100%)" = "dominant occurrence (76-100%)",
                                                  "velmi hojny vyskyt (26-50%)" = "very common occurrence (26-50%)"
  ))



structure <- structure %>%
  mutate(`coverage E1 layer - total vegetation [%]` = recode(`coverage E1 layer - total vegetation [%]`,
                                                             "velmi ridky vyskyt (do 0.2%)" = "very sparse occurrence (up to 0.2%)",
                                                             "malocetny vyskyt (od 1-5%)" = "infrequent occurrence (from 1-5%)",
                                                             "velkoplosny vyskyt (51-75%)" = "large-scale occurrence (51-75%)",
                                                             "ridky vyskyt (od 0.2 do 1%)" = "sparse occurrence (from 0.2 to 1%)",
                                                             "ojedinely vyskyt" = "isolated occurrence",
                                                             "hojny vyskyt (6-25%)" = "abundant occurrence (6-25%)",
                                                             "nevyskytuje se" = "absent",
                                                             "velmi hojny vyskyt (25-50%)" = "very common occurrence (25-50%)",
                                                             "dominantni vyskyt (76-100%)" = "dominant occurrence (76-100%)",
                                                             "velmi hojny vyskyt (26-50%)" = "very common occurrence (26-50%)"
  ))


structure <- structure %>%
  mutate(`coverage - dwarf pine [%]` = recode(`coverage - dwarf pine [%]`,
                                              "nevyskytuje se" = "absent",
                                              "malocetny vyskyt (1-5%)" = "infrequent occurrence (1-5%)"
  ))


structure <- structure %>%
  mutate(`coverage E2 layer - woody plants (1. 3-5m) [%]` = recode(`coverage E2 layer - woody plants (1. 3-5m) [%]`,
                                                                   "hojny vyskyt (6-25%)" = "abundant occurrence (6-25%)",
                                                                   "velmi hojny vyskyt (25-50%)" = "very common occurrence (25-50%)",
                                                                   "bez vyskytu" = "absent",
                                                                   "malocetny vyskyt (1-5%)" = "infrequent occurrence (1-5%)",
                                                                   "ojedinely vyskyt" = "isolated occurrence",
                                                                   "dominantni vyskyt (76-100%)" = "dominant occurrence (76-100%)",
                                                                   "velkoplosny vyskyt (51-75%)" = "large-scale occurrence (51-75%)",
                                                                   "ridky vyskyt (0.2-1%)" = "sparse occurrence (0.2-1%)",
                                                                   "ridky vyslyt (0.2-1%)" = "sparse occurrence (0.2-1%)"
  ))


structure <- structure %>%
  mutate(`Past management (logging history)` = recode(`Past management (logging history)`,
                                                      "bez tezby, zadne parezy" = "no logging, no stumps",
                                                      "tezba s ponechanim vetsiny hmoty" = "logging with most biomass left in place",
                                                      "tezba bez ponechani hmoty" = "logging without leaving biomass",
                                                      "tezba s ponechanim mensiny hmoty" = "logging with minority of biomass left in place"
  ))




structure <- structure %>%  
  mutate(`Plot surface appearance` = recode(`Plot surface appearance`,
                                            "zivy les" = "living forest",
                                            "zadna z predchozich moznosti" = "none of the previous options",
                                            "niva potoka/reky" = "stream/river floodplain",
                                            "velkoplosne vyvraty a polomy" = "large-scale uprooting and windthrow",
                                            "kurovec" = "bark beetle",
                                            "porostni mezera" = "canopy gap",
                                            "redina po pastve" = "pasture woodland",
                                            "holina po tezbe" = "clear-cut area"
  ))







# For structure_df
str(structure)
head(structure, 20)
colnames(structure)



write.csv(structure, "Lichens/Structural_data_english.csv", row.names = FALSE)









# =============================================================================
# LICHEN MODELING - DATA CLEANING WORKFLOW
# =============================================================================
# Goal: Clean structural predictors and merge with lichen response variables
# Author: Jonas Dotterweich
# Date: 2026-01-08
# =============================================================================

library(tidyverse)
library(usdm)        # for VIF calculation
library(corrplot)   # for correlation visualization
library(naniar)     # for missing data visualization

# =============================================================================
# STEP 1: LOAD DATA
# =============================================================================

#structure <- read_csv("structure.csv")
#lichen <- read_csv("lichen_species_groups.csv")

lichen <- lichen_groups

# Clean column names (remove hidden characters)
colnames(structure) <- trimws(colnames(structure))
colnames(structure) <- gsub("[\r\n]", "", colnames(structure))

cat("Structure data:", nrow(structure), "plots ×", ncol(structure), "variables\n")
cat("Lichen data:", nrow(lichen), "plots ×", ncol(lichen), "variables\n")

# Verify ID matching
missing_plots <- lichen$Plot[!lichen$Plot %in% structure$Project_ID]
if(length(missing_plots) > 0) {
  cat("\n⚠️  WARNING:", length(missing_plots), "lichen plots not found in structure data:\n")
  print(missing_plots)
}

# =============================================================================
# STEP 2: SELECT & RENAME PRIORITY PREDICTORS
# =============================================================================

structure_clean <- structure %>%
  select(
    project_id = Project_ID,
    
    # DEADWOOD (critical for calicioids - lignicolous species)
    deadwood_total = `volume of decaying wood (logs+stumps+snags) [m3/ha]`,
    deadwood_decay1 = `volume of decaying wood - decay stage 1 [m3/ha]`,
    deadwood_decay2 = `volume of decaying wood - decay stage 2 [m3/ha]`,
    deadwood_decay3 = `volume of decaying wood - decay stage 3 [m3/ha]`,
    deadwood_decay4 = `volume of decaying wood - decay stage 4 [m3/ha]`,
    deadwood_decay5 = `volume of decaying wood - decay stage 5 [m3/ha]`,
    volume_snags = `volume of standing deadwood [m3/ha]`,
    volume_lying_logs = `volume of lying logs [m3/ha]`,
    stump_volume = `stump volume [m3/ha]`,
    
    # TREE SIZE (old-growth proxy)
    dbh_mean = `mean DBH [mm] (live+dead)`,
    dbh_median = `median DBH [mm] (live+dead)`,
    dbh_max = `max_DBH [mm] (live+dead)`,
    dbh_sd = `sd_DBH (live+dead)`,
    n_dead_trees_50cm = `number of dead trees with DBH>50cm`,
    n_living_trees_80cm = `number of living with DBH>80cm`,
    
    # CANOPY STRUCTURE (microclimate)
    canopy_cover = `canopy cover of main tree layer [%]`,
    total_cover = `total canopy coverage [%]`,
    canopy_e3_estimate = `canopy cover E3 layer estimate [%]`,
    
    # STAND CHARACTERISTICS
    tree_height_median = `median height of main canopy (live+dead trees) [m]`,
    regeneration_density = `regeneration density [pcs/ha]`,
    
    # COMPOSITION (tree species - affects bark substrate availability for lichens)
    ba_spruce = `% basal area - spruce (live)`,
    ba_beech = `% basal area - beech (live)`,
    ba_fir = `% basal area - fir (live)`,
    ba_late_successional = `% basal area - late successional spp (maple,sycamore,oak,ash,elm,linden) (live)`,
    ba_early_successional = `% basal area - early successional spp (birch,aspen,poplar,pine,alder,rowan) (live)`,
    dominant_species = `dominant tree species (by basal area)`,
    
    # UNDERSTORY LAYERS (potential microclimate indicators)
    e1_woody_regen = `coverage E1 layer - woody regeneration [%]`,
    e1_herbs = `coverage E1 layer - herbs [%]`,
    e1_total_veg = `coverage E1 layer - total vegetation [%]`,
    
    # MANAGEMENT & SITE
    management = `Past management (logging history)`,
    exposure = exposure,
    elevation = `Elevation (m above sea level)`
  )

cat("\n✅ Selected", ncol(structure_clean)-1, "priority predictors\n")

# =============================================================================
# STEP 3: MISSING VALUE ASSESSMENT
# =============================================================================

cat("\n=== MISSING VALUE REPORT ===\n")

missing_summary <- structure_clean %>%
  summarise(across(everything(), ~sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
  mutate(
    pct_missing = round(n_missing / nrow(structure_clean) * 100, 1),
    status = case_when(
      n_missing == 0 ~ "✅ Complete",
      pct_missing < 5 ~ "⚠️  <5% missing",
      pct_missing < 20 ~ "🔴 5-20% missing",
      TRUE ~ "❌ >20% missing"
    )
  ) %>%
  arrange(desc(n_missing))

print(missing_summary, n = 35)

# Visualize missing data pattern
gg_miss_var(structure_clean, show_pct = TRUE) +
  labs(title = "Missing Data by Variable",
       subtitle = "Forest structural predictors (n=120 plots)")

# =============================================================================
# STEP 4: CREATE DEADWOOD COMPOSITE VARIABLES
# =============================================================================

# Sum decay stages 2-5 (exclude decay1 = fresh deadwood)
structure_clean <- structure_clean %>%
  mutate(
    deadwood_decay2to5 = deadwood_decay2 + deadwood_decay3 + 
      deadwood_decay4 + deadwood_decay5,
    
    # Also create decay 4-5 (highly decomposed = best for specialist lichens)
    deadwood_decay4to5 = deadwood_decay4 + deadwood_decay5,
    
    # Calculate % of deadwood in advanced decay
    pct_decay_advanced = ifelse(deadwood_total > 0, 
                                (deadwood_decay4to5 / deadwood_total) * 100, 
                                0)
  )

cat("\n✅ Created composite deadwood variables:\n")
cat("   • deadwood_decay2to5 (sum of stages 2-5)\n")
cat("   • deadwood_decay4to5 (sum of stages 4-5)\n")
cat("   • pct_decay_advanced (% of deadwood in stages 4-5)\n")

# =============================================================================
# STEP 5: HANDLE MISSING VALUES
# =============================================================================

# Strategy: Remove plots with >3 missing values in key predictors, impute median for others

key_predictors <- c("dbh_max", "deadwood_total", "canopy_cover", 
                    "ba_spruce", "ba_beech", "tree_height_median")

plots_with_many_na <- structure_clean %>%
  mutate(n_na_key = rowSums(is.na(select(., all_of(key_predictors))))) %>%
  filter(n_na_key > 2) %>%
  pull(project_id)

if(length(plots_with_many_na) > 0) {
  cat("\n⚠️  Removing", length(plots_with_many_na), "plots with >2 missing KEY predictors:\n")
  print(plots_with_many_na)
  structure_clean <- structure_clean %>%
    filter(! project_id %in% plots_with_many_na)
}

# Impute remaining missing values with median (for numeric variables)
structure_clean <- structure_clean %>%
  mutate(across(where(is.numeric), ~ifelse(is.na(.), median(., na.rm = TRUE), .)))

cat("\n✅ After cleaning:", nrow(structure_clean), "plots remaining\n")

# =============================================================================
# STEP 6: ENCODE CATEGORICAL VARIABLES
# =============================================================================

# MANAGEMENT:   Convert to ordered factor
cat("\n=== MANAGEMENT CATEGORIES ===\n")
print(table(structure_clean$management))

structure_clean <- structure_clean %>%
  mutate(
    management_factor = factor(management, 
                               levels = c(
                                 "no logging, no stumps",
                                 "logging with most biomass left in place",
                                 "logging with partial biomass removal",
                                 "logging with complete biomass removal"
                               ),
                               ordered = TRUE
    ),
    
    # Binary:   logged vs not logged
    logged = ifelse(management == "no logging, no stumps", 0, 1),
    
    # Logging intensity (0-3)
    logging_intensity = as.numeric(management_factor) - 1
  )

# EXPOSURE:  Simplify to main directions
structure_clean <- structure_clean %>%
  mutate(
    exposure_simple = case_when(
      str_detect(exposure, "not evaluated") ~ "flat",
      str_detect(exposure, "^N ") ~ "N",
      str_detect(exposure, "^S ") ~ "S",
      str_detect(exposure, "^E ") ~ "E",
      str_detect(exposure, "^W ") ~ "W",
      str_detect(exposure, "NE|SE") ~ "NE_SE",
      str_detect(exposure, "NW|SW") ~ "NW_SW",
      TRUE ~ "other"
    ),
    exposure_factor = factor(exposure_simple)
  )

cat("\n=== EXPOSURE CATEGORIES ===\n")
print(table(structure_clean$exposure_simple))

# DOMINANT SPECIES:  Group rare species
structure_clean <- structure_clean %>%
  mutate(
    dominant_grouped = case_when(
      dominant_species %in% c("SM", "BK", "JD") ~ dominant_species,
      TRUE ~ "other"
    ),
    dominant_factor = factor(dominant_grouped, levels = c("SM", "BK", "JD", "other"))
  )

cat("\n=== DOMINANT SPECIES ===\n")
print(table(structure_clean$dominant_grouped))

cat("\n✅ Categorical variables encoded\n")

# =============================================================================
# STEP 7: COLLINEARITY CHECK
# =============================================================================

# Select only numeric predictors for correlation analysis
numeric_predictors <- structure_clean %>%
  select(
    deadwood_total, deadwood_decay2to5, deadwood_decay4to5,
    dbh_mean, dbh_median, dbh_max, dbh_sd,
    canopy_cover, total_cover,
    tree_height_median,
    ba_spruce, ba_beech, ba_late_successional, ba_early_successional,
    n_dead_trees_50cm, n_living_trees_80cm,
    logging_intensity,
    elevation
  ) %>%
  drop_na()

# Correlation matrix
cor_matrix <- cor(numeric_predictors, use = "complete.obs")

# Visualize
png("correlation_matrix.png", width = 1000, height = 1000)
corrplot(cor_matrix, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45, tl.cex = 0.8,
         title = "Correlation Matrix - Structural Predictors",
         mar = c(0,0,2,0))
dev.off()

cat("\n✅ Correlation matrix saved to correlation_matrix.png\n")

# Find highly correlated pairs (|r| > 0.7)
high_cor <- which(abs(cor_matrix) > 0.7 & cor_matrix != 1, arr.ind = TRUE)
if(nrow(high_cor) > 0) {
  cat("\n⚠️  HIGHLY CORRELATED PAIRS (|r| > 0.7):\n")
  high_cor_pairs <- data.frame(
    var1 = rownames(cor_matrix)[high_cor[,1]],
    var2 = colnames(cor_matrix)[high_cor[,2]],
    correlation = round(cor_matrix[high_cor], 3)
  ) %>%
    filter(var1 < var2) %>%  # Remove duplicates
    arrange(desc(abs(correlation)))
  
  print(high_cor_pairs, row.names = FALSE)
}

# =============================================================================
# STEP 8: VARIANCE INFLATION FACTOR (VIF) - USING CAR PACKAGE
# =============================================================================

cat("\n=== CALCULATING VIF ===\n")

# Install car if needed (uncomment if not installed):
# install.packages("car")

library(car)

# Fit a linear model with all numeric predictors to calculate VIF
# Note: Response doesn't matter for VIF, we just need the predictor structure
vif_model <- lm(deadwood_total ~ 
                  deadwood_decay2to5 + deadwood_decay4to5 +
                  dbh_mean + dbh_median + dbh_max + dbh_sd +
                  canopy_cover + total_cover +
                  tree_height_median +
                  ba_spruce + ba_beech + ba_late_successional + ba_early_successional +
                  n_dead_trees_50cm + n_living_trees_80cm +
                  logging_intensity +
                  elevation,
                data = numeric_predictors)

# Calculate VIF
vif_values <- vif(vif_model)

# Create formatted table
vif_table <- data.frame(
  variable = names(vif_values),
  VIF = round(vif_values, 2),
  status = case_when(
    vif_values < 5 ~ "✅ OK",
    vif_values < 10 ~ "⚠️  Moderate",
    TRUE ~ "🔴 High collinearity"
  )
) %>%
  arrange(desc(VIF))

print(vif_table, row.names = FALSE)

cat("\n✅ VIF calculated using car package\n")
cat("\nInterpretation:\n")
cat("  VIF < 5:  ✅ No collinearity concern\n")
cat("  VIF 5-10: ⚠️  Moderate collinearity (consider removing)\n")
cat("  VIF > 10: 🔴 High collinearity (MUST remove one of the pair)\n")

# Identify problematic variables
problematic_vars <- vif_table %>%
  filter(VIF > 10)

if(nrow(problematic_vars) > 0) {
  cat("\n🔴 VARIABLES WITH VIF > 10 (REMOVE THESE):\n")
  print(problematic_vars, row.names = FALSE)
}


# =============================================================================
# STEP 9: VARIABLE SELECTION - FINAL PREDICTOR SET
# =============================================================================

# Based on collinearity + ecological relevance, select core predictors

structure_final <- structure_clean %>%
  select(
    project_id,
    
    # OLD-GROWTH PROXIES (keep max, drop mean/median due to collinearity)
    dbh_max,
    n_living_trees_80cm,  # Presence of very large living trees
    n_dead_trees_50cm,    # Large dead trees (snags)
    
    # DEADWOOD (key for calicioids)
    deadwood_decay4to5,  # Highly decomposed = specialist habitat
    deadwood_total,      # Overall deadwood availability
    
    # CANOPY (microclimate)
    canopy_cover,
    
    # COMPOSITION (tree species composition affects bark substrate)
    ba_spruce,
    ba_beech,
    ba_late_successional,
    
    # STAND STRUCTURE
    tree_height_median,
    elevation,
    
    # MANAGEMENT
    logging_intensity,
    logged,
    
    # CATEGORICAL
    management_factor,
    exposure_factor,
    dominant_factor
  )

cat("\n✅ Final predictor set:", ncol(structure_final)-1, "variables\n")
cat("   • 12 continuous predictors\n")
cat("   • 3 categorical predictors\n")

# =============================================================================
# STEP 10: SCALE CONTINUOUS PREDICTORS
# =============================================================================

structure_scaled <- structure_final %>%
  mutate(across(
    c(dbh_max, n_living_trees_80cm, n_dead_trees_50cm,
      deadwood_decay4to5, deadwood_total,
      canopy_cover, ba_spruce, ba_beech, ba_late_successional,
      tree_height_median, elevation),
    list(scaled = ~as.numeric(scale(.))),
    .names = "{.col}_scaled"
  ))

# Check for extreme outliers (>3 SD)
outlier_check <- structure_scaled %>%
  select(project_id, ends_with("_scaled")) %>%
  pivot_longer(-project_id, names_to = "variable", values_to = "z_score") %>%
  filter(abs(z_score) > 3)

if(nrow(outlier_check) > 0) {
  cat("\n⚠️  EXTREME OUTLIERS DETECTED (|Z| > 3):\n")
  print(outlier_check)
  cat("\nConsider investigating these plots or using robust scaling (IQR)\n")
} else {
  cat("\n✅ No extreme outliers detected\n")
}

cat("\n✅ Continuous predictors scaled (mean=0, SD=1)\n")

# =============================================================================
# STEP 11: MERGE WITH LICHEN DATA
# =============================================================================

modeling_data <- lichen %>%
  rename(project_id = Plot) %>%
  inner_join(structure_scaled, by = "project_id")

cat("\n=== FINAL MODELING DATASET ===\n")
cat("Dimensions:", nrow(modeling_data), "plots ×", ncol(modeling_data), "variables\n")
cat("\nResponse variables:\n")
cat("  • calicioids_richness (count:  0-8 species)\n")
cat("  • parmelia_agg_presence (binary)\n")
cat("  • ochrolechia_presence (binary)\n")
cat("  • core_ogf_presence (binary)\n")
cat("  • mycoblastus_presence (binary)\n")
cat("  • xylographa_presence (binary)\n")
cat("  • elite_rare_presence (binary)\n")
cat("\nPredictors:   12 continuous (scaled) + 3 categorical\n")

# =============================================================================
# STEP 12: SAVE CLEANED DATA
# =============================================================================

write_csv(modeling_data, "modeling_data_final.csv")
write_csv(missing_summary, "data_cleaning_report_missing.csv")

# Save VIF table (handle both usdm and manual formats)
if(exists("vif_table")) {
  write_csv(as.data.frame(vif_table), "data_cleaning_report_vif.csv")
}

cat("\n✅ FILES SAVED:\n")
cat("   • modeling_data_final.csv (", nrow(modeling_data), " plots)\n", sep = "")
cat("   • data_cleaning_report_missing.csv\n")
cat("   • data_cleaning_report_vif.csv\n")
cat("   • correlation_matrix.png\n")

# =============================================================================
# STEP 13: SUMMARY STATISTICS
# =============================================================================

cat("\n=== SUMMARY STATISTICS - CONTINUOUS PREDICTORS ===\n")
structure_final %>%
  select(dbh_max: elevation) %>%
  summary() %>%
  print()

cat("\n=== SUMMARY STATISTICS - CATEGORICAL PREDICTORS ===\n")
cat("\nManagement:\n")
print(table(structure_final$management_factor))

cat("\nExposure:\n")
print(table(structure_final$exposure_factor))

cat("\nDominant species:\n")
print(table(structure_final$dominant_factor))

# =============================================================================
# STEP 14: DISTRIBUTION CHECKS (for GLM assumptions)
# =============================================================================

# Check for zero-inflation in deadwood variables
cat("\n=== ZERO-INFLATION CHECK ===\n")
zero_counts <- structure_final %>%
  summarise(
    pct_zero_deadwood_total = mean(deadwood_total == 0, na.rm = TRUE) * 100,
    pct_zero_deadwood_decay45 = mean(deadwood_decay4to5 == 0, na.rm = TRUE) * 100,
    pct_zero_n_dead_trees = mean(n_dead_trees_50cm == 0, na.rm = TRUE) * 100
  )
print(zero_counts)

if(zero_counts$pct_zero_deadwood_decay45 > 50) {
  cat("⚠️  High zero-inflation in deadwood_decay4to5 (", 
      round(zero_counts$pct_zero_deadwood_decay45, 1), 
      "%) - consider zero-inflated models or hurdle models\n", sep = "")
}

# Histograms of key predictors
png("predictor_distributions.png", width = 1400, height = 1000)
par(mfrow = c(3, 4))
hist(structure_final$dbh_max, main = "Max DBH", xlab = "mm", col = "lightblue", breaks = 20)
hist(structure_final$deadwood_total, main = "Total Deadwood", xlab = "m³/ha", col = "lightblue", breaks = 20)
hist(structure_final$deadwood_decay4to5, main = "Deadwood Decay 4-5", xlab = "m³/ha", col = "lightblue", breaks = 20)
hist(structure_final$canopy_cover, main = "Canopy Cover", xlab = "%", col = "lightblue", breaks = 20)
hist(structure_final$ba_spruce, main = "Spruce BA %", xlab = "%", col = "lightblue", breaks = 20)
hist(structure_final$ba_beech, main = "Beech BA %", xlab = "%", col = "lightblue", breaks = 20)
hist(structure_final$ba_late_successional, main = "Late Successional BA %", xlab = "%", col = "lightblue", breaks = 20)
hist(structure_final$tree_height_median, main = "Tree Height (median)", xlab = "m", col = "lightblue", breaks = 20)
hist(structure_final$elevation, main = "Elevation", xlab = "m a.s.l.", col = "lightblue", breaks = 20)
hist(structure_final$n_living_trees_80cm, main = "N Large Living Trees", xlab = "count", col = "lightblue", breaks = 20)
hist(structure_final$n_dead_trees_50cm, main = "N Large Dead Trees", xlab = "count", col = "lightblue", breaks = 20)
hist(structure_final$logging_intensity, main = "Logging Intensity", xlab = "0-3 scale", col = "lightblue", breaks = 4)
par(mfrow = c(1, 1))
dev.off()

cat("\n✅ Predictor distributions saved to predictor_distributions.png\n")

cat("\n=================================================\n")
cat("✅ DATA CLEANING COMPLETE!\n")
cat("=================================================\n")
cat("\nNext steps:\n")
cat("1. Review VIF and correlation reports\n")
cat("2. Decide on final predictor set based on collinearity\n")
cat("3. Begin GLM modeling with modeling_data_final.csv\n")
cat("4. Consider transformations for highly skewed predictors\n")
cat("\nSuggested first models:\n\n")
cat("# CALICIOIDS RICHNESS (negative binomial):\n")
cat("library(MASS)\n")
cat("model_calicioids <- glm. nb(calicioids_richness ~ \n")
cat("                    dbh_max_scaled + deadwood_decay4to5_scaled + \n")
cat("                    canopy_cover_scaled + ba_spruce_scaled + \n")
cat("                    logging_intensity, \n")
cat("                    data = modeling_data)\n")
cat("summary(model_calicioids)\n\n")
cat("# PARMELIA AGG PRESENCE (binomial):\n")
cat("model_parmelia <- glm(parmelia_agg_presence ~ \n")
cat("                  dbh_max_scaled + ba_beech_scaled + \n")
cat("                  elevation_scaled + logged, \n")
cat("                  family = binomial, \n")
cat("                  data = modeling_data)\n")
cat("summary(model_parmelia)\n")











# =============================================================================
# DATA-DRIVEN DECAY STAGE GROUPING ANALYSIS
# =============================================================================



# Merge for analysis
test_data <- lichen %>%
  rename(project_id = Plot) %>%
  inner_join(structure_clean, by = "project_id")

# =============================================================================
# TEST 1: CORRELATIONS BETWEEN DECAY STAGES AND LICHEN RESPONSES
# =============================================================================

cat("=== CORRELATION ANALYSIS:  DECAY STAGES vs LICHEN RICHNESS ===\n\n")

# Calculate correlations for each decay stage
decay_correlations <- test_data %>%
  summarise(
    decay1_vs_calicioids = cor(deadwood_decay1, calicioids_richness, use = "complete.obs"),
    decay2_vs_calicioids = cor(deadwood_decay2, calicioids_richness, use = "complete.obs"),
    decay3_vs_calicioids = cor(deadwood_decay3, calicioids_richness, use = "complete.obs"),
    decay4_vs_calicioids = cor(deadwood_decay4, calicioids_richness, use = "complete.obs"),
    decay5_vs_calicioids = cor(deadwood_decay5, calicioids_richness, use = "complete. obs")
  ) %>%
  pivot_longer(everything(), names_to = "decay_stage", values_to = "correlation") %>%
  mutate(decay_stage = str_extract(decay_stage, "decay[0-9]"))

print(decay_correlations)

cat("\n🔍 INTERPRETATION:\n")
cat("Which decay stages have STRONGEST correlation with calicioids?\n")
cat("→ Those are the ones we should GROUP TOGETHER\n\n")

# =============================================================================
# TEST 2: COMPARE DIFFERENT GROUPING STRATEGIES
# =============================================================================

test_data <- test_data %>%
  mutate(
    # Option A:  Decay 4-5 only (my suggestion - but is it justified?)
    decay_4to5 = deadwood_decay4 + deadwood_decay5,
    
    # Option B: Decay 3-5 (your suggestion!)
    decay_3to5 = deadwood_decay3 + deadwood_decay4 + deadwood_decay5,
    
    # Option C:  Decay 2-5 (exclude only fresh)
    decay_2to5 = deadwood_decay2 + deadwood_decay3 + deadwood_decay4 + deadwood_decay5,
    
    # Option D:  Keep all separate (no grouping)
    # (use individual stages)
  )

# Compare correlations
grouping_comparison <- test_data %>%
  summarise(
    total = cor(deadwood_total, calicioids_richness, use = "complete.obs"),
    decay_2to5 = cor(decay_2to5, calicioids_richness, use = "complete.obs"),
    decay_3to5 = cor(decay_3to5, calicioids_richness, use = "complete.obs"),
    decay_4to5 = cor(decay_4to5, calicioids_richness, use = "complete.obs"),
    decay5_only = cor(deadwood_decay5, calicioids_richness, use = "complete.obs")
  ) %>%
  pivot_longer(everything(), names_to = "grouping", values_to = "correlation") %>%
  arrange(desc(abs(correlation)))

cat("=== WHICH GROUPING HAS STRONGEST CORRELATION? ===\n")
print(grouping_comparison)

cat("\n🎯 THE DATA WILL TELL US which grouping is best!\n\n")

# =============================================================================
# TEST 3: VISUALIZE RELATIONSHIPS
# =============================================================================

library(ggplot2)

# Plot each decay stage vs calicioid richness
p1 <- ggplot(test_data, aes(x = deadwood_decay1, y = calicioids_richness)) +
  geom_point() + geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Decay Stage 1", x = "Volume [m³/ha]", y = "Calicioid Richness")

p2 <- ggplot(test_data, aes(x = deadwood_decay2, y = calicioids_richness)) +
  geom_point() + geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Decay Stage 2", x = "Volume [m³/ha]", y = "Calicioid Richness")

p3 <- ggplot(test_data, aes(x = deadwood_decay3, y = calicioids_richness)) +
  geom_point() + geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Decay Stage 3", x = "Volume [m³/ha]", y = "Calicioid Richness")

p4 <- ggplot(test_data, aes(x = deadwood_decay4, y = calicioids_richness)) +
  geom_point() + geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Decay Stage 4", x = "Volume [m³/ha]", y = "Calicioid Richness")

p5 <- ggplot(test_data, aes(x = deadwood_decay5, y = calicioids_richness)) +
  geom_point() + geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Decay Stage 5", x = "Volume [m³/ha]", y = "Calicioid Richness")

library(patchwork)
(p1 + p2 + p3) / (p4 + p5)

ggsave("decay_stage_correlations.png", width = 12, height = 8)

cat("✅ Saved plot:  decay_stage_correlations.png\n")
cat("👁️  LOOK AT THE PLOTS:  Which stages show clear positive relationships?\n\n")

# =============================================================================
# TEST 4: MODEL COMPARISON (AIC-based selection)
# =============================================================================

library(MASS)

cat("=== MODEL COMPARISON:  WHICH GROUPING FITS BEST? ===\n\n")

# Model A: Total deadwood
model_total <- glm.nb(calicioids_richness ~ deadwood_total + dbh_max + 
                        canopy_cover + logging_intensity, 
                      data = test_data)

# Model B: Decay 2-5
model_2to5 <- glm.nb(calicioids_richness ~ decay_2to5 + dbh_max + 
                       canopy_cover + logging_intensity, 
                     data = test_data)

# Model C: Decay 3-5 (YOUR suggestion)
model_3to5 <- glm.nb(calicioids_richness ~ decay_3to5 + dbh_max + 
                       canopy_cover + logging_intensity, 
                     data = test_data)

# Model D:  Decay 4-5 (MY suggestion)
model_4to5 <- glm.nb(calicioids_richness ~ decay_4to5 + dbh_max + 
                       canopy_cover + logging_intensity, 
                     data = test_data)

# Model E: Individual stages (no grouping)
model_separate <- glm.nb(calicioids_richness ~ deadwood_decay3 + deadwood_decay4 + 
                           deadwood_decay5 + dbh_max + canopy_cover + logging_intensity, 
                         data = test_data)

# Compare AIC
aic_comparison <- data.frame(
  Model = c("Total deadwood", "Decay 2-5", "Decay 3-5", "Decay 4-5", "Separate 3,4,5"),
  AIC = c(AIC(model_total), AIC(model_2to5), AIC(model_3to5), 
          AIC(model_4to5), AIC(model_separate)),
  Delta_AIC = NA
) %>%
  mutate(Delta_AIC = AIC - min(AIC)) %>%
  arrange(AIC)

cat("AIC COMPARISON (lower = better fit):\n")
print(aic_comparison)

cat("\n📊 INTERPRETATION:\n")
cat("  • Delta_AIC < 2:  Models are equivalent\n")
cat("  • Delta_AIC 2-7:   Moderate support for lower AIC model\n")
cat("  • Delta_AIC > 10: Strong support for lower AIC model\n\n")

cat("🔬 THE MODEL WITH LOWEST AIC IS THE BEST GROUPING!\n")

# =============================================================================

cat("\n=== ARGUMENTS FOR DECAY 3-5 vs 4-5 ===\n\n")

cat("FOR DECAY 3-5:\n")
cat("  ✅ Includes 'mid-decay' wood that may still support lichens\n")
cat("  ✅ Larger sample (more wood volume → more statistical power)\n")
cat("  ✅ May be conservative cutoff for 'old' wood in Czech forests\n")
cat("  ✅ Stage 3 = bark mostly gone = substrate exposed\n")
cat("  ✅ If correlation analysis shows decay3 also correlates, include it!\n\n")

cat("FOR DECAY 4-5:\n")
cat("  ✅ More restrictive = stronger old-growth signal?\n")
cat("  ✅ If decay3 has NO correlation, don't include it (dilutes signal)\n")
cat("  ⚠️  But this is an ASSUMPTION until we test it!\n\n")

cat("🧪 EMPIRICAL TEST:\n")
cat("Run the correlation analysis above.\n")
cat("If cor(decay3, calicioids) > 0.3 AND significant → INCLUDE IT\n")
cat("If cor(decay3, calicioids) ≈ 0 → EXCLUDE IT\n")

cat("\n=== MY REVISED RECOMMENDATION ===\n\n")

cat("1. RUN THE CORRELATION ANALYSIS FIRST\n")
cat("   → See which decay stages actually predict calicioid richness\n\n")

cat("2. TEST MULTIPLE GROUPINGS WITH AIC COMPARISON\n")
cat("   → Let model fit decide, not assumptions\n\n")

cat("3. CHECK THE LITERATURE:\n")
cat("   → Does Malíček et al. (2019) mention decay stages?\n")
cat("   → Are there Czech forestry standards for what each stage means?\n\n")

cat("4. START CONSERVATIVE:\n")
cat("   → If uncertain, use decay_3to5 (broader definition)\n")
cat("   → You can always test decay_4to5 as sensitivity analysis\n\n")

cat("5. REPORT TRANSPARENTLY:\n")
cat("   → 'We grouped decay 3-5 based on [DATA-DRIVEN REASON]'\n")
cat("   → NOT 'We assumed stage 3 is too young'\n")





# Quick test script to paste into R: 

# Calculate correlations
cor_test <- lichen_groups %>%
  left_join(structure, by = c("Plot" = "Project_ID")) %>%
  summarise(
    r_decay1 = cor(calicioids_richness, `volume of decaying wood - decay stage 1 [m3/ha]`, use = "complete.obs"),
    r_decay2 = cor(calicioids_richness, `volume of decaying wood - decay stage 2 [m3/ha]`, use = "complete.obs"),
    r_decay3 = cor(calicioids_richness, `volume of decaying wood - decay stage 3 [m3/ha]`, use = "complete.obs"),
    r_decay4 = cor(calicioids_richness, `volume of decaying wood - decay stage 4 [m3/ha]`, use = "complete.obs"),
    r_decay5 = cor(calicioids_richness, `volume of decaying wood - decay stage 5 [m3/ha]`, use = "complete.obs")
  )

print(cor_test)







### testing all Lichen Groups

# Test all lichen groups
all_group_correlations <- test_data %>%
  summarise(
    # CALICIOIDS (we already know:  weak)
    calicioids_decay45 = cor(calicioids_richness, decay_4to5, use = "complete.obs"),
    
    # PARMELIA (bark specialist - mature trees)
    parmelia_decay45 = cor(parmelia_agg_presence, decay_4to5, use = "complete.obs"),
    parmelia_dbh = cor(parmelia_agg_presence, dbh_max, use = "complete.obs"),
    
    # OCHROLECHIA (also bark/mature forest)
    ochrolechia_decay45 = cor(ochrolechia_presence, decay_4to5, use = "complete.obs"),
    ochrolechia_canopy = cor(ochrolechia_presence, canopy_cover, use = "complete.obs"),
    
    # CORE OGF (multi-species assemblage)
    core_decay45 = cor(core_ogf_presence, decay_4to5, use = "complete.obs"),
    
    # MYCOBLASTUS (bark specialist)
    myco_decay45 = cor(mycoblastus_presence, decay_4to5, use = "complete.obs"),
    
    # XYLOGRAPHA (lignicolous - DEADWOOD specialist!)
    xylo_decay45 = cor(xylographa_presence, decay_4to5, use = "complete.obs"),
    xylo_decay5 = cor(xylographa_presence, deadwood_decay5, use = "complete.obs"),
    
    # ELITE RARE
    elite_decay45 = cor(elite_rare_presence, decay_4to5, use = "complete.obs")
  )

print(all_group_correlations)

print(all_group_correlations, width = Inf)

all_group_correlations %>%
  pivot_longer(everything(), names_to = "relationship", values_to = "correlation") %>%
  mutate(
    abs_r = abs(correlation),
    strength = case_when(
      abs_r > 0.5 ~ "🔴 Strong",
      abs_r > 0.3 ~ "⚠️  Moderate", 
      abs_r > 0.1 ~ "✅ Weak",
      TRUE ~ "⚪ Negligible"
    )
  ) %>%
  arrange(desc(abs_r)) %>%
  print(n = Inf)










#=============================================================================
# CREATE SUBSET WITH SELECTED PREDICTORS (HIGH + MODERATE + 2 LOW)
# =============================================================================

# Define variables to keep based on your criteria
variables_to_keep <- c(
  # ID for joining
  "Project_ID",
  
  # HIGH USABILITY (16 variables)
  "Elevation (m above sea level)",
  "exposure",
  "volume of standing deadwood [m3/ha]",
  "volume of lying logs [m3/ha]",
  "volume of decaying wood - decay stage 4 [m3/ha]",
  "volume of decaying wood - decay stage 5 [m3/ha]",
  "median height of main canopy (live+dead trees) [m]",
  "canopy cover of main tree layer [%]",
  "number of dead trees with DBH>50cm",
  "number of living with DBH>80cm",
  "max_DBH [mm] (live+dead)",
  "dominant tree species (by basal area)",
  "% basal area - spruce (live)",
  "% basal area - beech (live)",
  "% basal area - fir (live)",
  "% basal area - late successional spp (maple,sycamore,oak,ash,elm,linden) (live)",
  "Past management (logging history)",
  
  # MODERATE USABILITY (11 variables)
  "volume of standing deadwood DBH>7cm in inner circle r<=7m [m3]",
  "volume of standing deadwood DBH>30cm in outer ring r>7m [m3]",
  "stump volume [m3/ha]",
  "total canopy coverage [%]",
  "regeneration density [pcs/ha]",
  "volume of decaying wood - decay stage 2 [m3/ha]",
  "volume of decaying wood - decay stage 3 [m3/ha]",
  "sd_DBH (live+dead)",
  "canopy cover E3 layer estimate [%]",
  "% basal area - early successional spp (birch,aspen,poplar,pine,alder,rowan) (live)",
  
  # LOW USABILITY (only 2 kept for context)
  "Habitat type (Natura 2000)",
  "Forest type"
)

# Check which variables exist in structure
missing_vars <- setdiff(variables_to_keep, colnames(structure))
if(length(missing_vars) > 0) {
  cat("⚠️  WARNING: These variables not found in structure:\n")
  print(missing_vars)
}

# Create subset dataframe
structure_selected <- structure %>%
  dplyr::select(all_of(variables_to_keep))

cat("\n✅ SELECTED PREDICTOR SUBSET CREATED\n")
cat("Dimensions:", nrow(structure_selected), "plots ×", ncol(structure_selected), "variables\n\n")

cat("Variables included:\n")
cat("  • HIGH usability: 16 core predictors\n")
cat("  • MODERATE usability: 11 variables to test\n")
cat("  • LOW usability: 2 contextual variables (habitat, forest type)\n")
cat("  • Total: 29 predictors + 1 ID = 30 columns\n\n")

# Summary by category
cat("═══ VARIABLE BREAKDOWN ═══\n\n")

cat("🟢 HIGH USABILITY (n=16):\n")
high_vars <- c("Elevation (m above sea level)", "exposure", 
               "volume of standing deadwood [m3/ha]", "volume of lying logs [m3/ha]",
               "volume of decaying wood - decay stage 4 [m3/ha]", 
               "volume of decaying wood - decay stage 5 [m3/ha]",
               "median height of main canopy (live+dead trees) [m]",
               "canopy cover of main tree layer [%]",
               "number of dead trees with DBH>50cm", "number of living with DBH>80cm",
               "max_DBH [mm] (live+dead)", "dominant tree species (by basal area)",
               "% basal area - spruce (live)", "% basal area - beech (live)",
               "% basal area - fir (live)", 
               "% basal area - late successional spp (maple,sycamore,oak,ash,elm,linden) (live)",
               "Past management (logging history)")
for(v in high_vars) cat("  •", v, "\n")

cat("\n🟡 MODERATE USABILITY (n=11):\n")
mod_vars <- c("volume of standing deadwood DBH>7cm in inner circle r<=7m [m3]",
              "volume of standing deadwood DBH>30cm in outer ring r>7m [m3]",
              "stump volume [m3/ha]", "total canopy coverage [%]",
              "regeneration density [pcs/ha]",
              "volume of decaying wood - decay stage 2 [m3/ha]",
              "volume of decaying wood - decay stage 3 [m3/ha]",
              "sd_DBH (live+dead)", "canopy cover E3 layer estimate [%]",
              "% basal area - early successional spp (birch,aspen,poplar,pine,alder,rowan) (live)")
for(v in mod_vars) cat("  •", v, "\n")

cat("\n📋 CONTEXTUAL (n=2):\n")
cat("  • Habitat type (Natura 2000)\n")
cat("  • Forest type\n")



# Show column names for verification
cat("\n═══ COLUMN NAMES IN NEW DATAFRAME ═══\n")


print(colnames(structure_selected))

## subset not saved yet as coords will be added


#### adding the plot coordinates to the selected structure

coords <- read_excel("Lichens/Biodiversity_120 plot.xlsx")

#Check structure
cat("Coordinate file dimensions:", nrow(coords), "rows ×", ncol(coords), "columns\n")
cat("Column names:\n")
print(colnames(coords))

# Clean column names (remove special characters)
colnames(coords) <- trimws(colnames(coords))

# Rename ID column for joining
coords_clean <- coords %>%
  rename(Project_ID = č.,
         coordinates = souradnice) %>%  # Or whatever the exact column name is
  dplyr::select(
    Project_ID,
    X,
    Y,
    coordinates  # lat/long string
  )

cat("\n✅ Coordinate data cleaned:\n")
print(head(coords_clean))

# Join to structure_selected
structure_selected_with_coords <- structure_selected %>%
  left_join(coords_clean, by = "Project_ID")

# Verify the join
cat("\n═══ JOIN VERIFICATION ═══\n")
cat("Rows in structure_selected:", nrow(structure_selected), "\n")
cat("Rows after join:", nrow(structure_selected_with_coords), "\n")
cat("Rows with missing X:", sum(is.na(structure_selected_with_coords$X)), "\n")
cat("Rows with missing Y:", sum(is.na(structure_selected_with_coords$Y)), "\n")

# Check for plots that didn't match
if(sum(is.na(structure_selected_with_coords$X)) > 0) {
  cat("\n⚠️  WARNING: These plots have no coordinates:\n")
  missing_coords <- structure_selected_with_coords %>%
    filter(is.na(X)) %>%
    pull(Project_ID)
  print(missing_coords)
} else {
  cat("\n✅ SUCCESS: All", nrow(structure_selected_with_coords), "plots matched!\n")
}



# Check coordinate ranges (S-JTSK system)
cat("\n═══ COORDINATE RANGES (S-JTSK) ═══\n")
cat("X range:", range(structure_selected_with_coords$X, na.rm = TRUE), "\n")
cat("Y range:", range(structure_selected_with_coords$Y, na.rm = TRUE), "\n")
cat("(Negative values are normal for S-JTSK coordinate system)\n")





# Update structure_selected
structure_selected <- structure_selected_with_coords

# Save updated dataset
write.csv(structure_selected, "Lichens/structure_selected_with_coordinates.csv", row.names = FALSE)

cat("\n✅ Coordinates added successfully!\n")
cat("✅ Saved to:  structure_selected_with_coordinates. csv\n")
cat("\nNew dimensions:", nrow(structure_selected), "plots ×", ncol(structure_selected), "variables\n")
cat("(Added 3 columns: X, Y, souradnice)\n")



####------------------------------------
# from here on, we continue in the script: PRE-MODELING DIAGNOSTICS