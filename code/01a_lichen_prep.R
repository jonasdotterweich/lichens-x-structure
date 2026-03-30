# ==============================================================================
# 01a_lichen_prep.R
# ==============================================================================
# PURPOSE  : Load, translate, explore, group, and export lichen data for
#            the Sumava NP case study.
# OUTPUT   : outputs/Sumava_CZ/lichen_clean.csv
#              Canonical schema: plot_id | X | Y | <group>_presence (0/1 int)
#                                              | <group>_richness (int)
# AUTHOR   : Jonas Dotterweich
# STUDY    : Sumava NP, Czech Republic (120 plots)
# NOTE     : All study-specific decisions (translations, species groupings,
#            indicator lists) live inside this file by design. They are NOT
#            moved to a shared config because they are dataset-specific and
#            will differ for every new case study.
# ==============================================================================


# ==============================================================================
# 0. SETUP
# ==============================================================================

library(tidyverse)
library(readxl)
library(ggplot2)
library(pheatmap)


#individual working directory for loading the data

# Paths (relative to project root — edit only these two lines for a new study)
PATH_LICHEN_RAW  <- "Lichens/Licen_data.xlsx"
PATH_COORDS_RAW  <- "Lichens/Biodiversity_120 plot.xlsx"
PATH_OUT_CLEAN   <- "studies/sumava_CZ/outputs/lichen_clean.csv"
PATH_OUT_REPORTS <- "studies/sumava_CZ/outputs/"   # QC / diagnostic CSVs

N_PLOTS <- 120   # Total survey plots in this study

dir.create(PATH_OUT_REPORTS, recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# 1. LOAD RAW DATA
# ==============================================================================

lichens_raw <- read_xlsx(PATH_LICHEN_RAW)

cat("Raw lichen data loaded:", nrow(lichens_raw), "records,",
    n_distinct(lichens_raw$species), "species,",
    n_distinct(lichens_raw$Plot), "plots\n")


# ==============================================================================
# 2. TRANSLATE CZECH LABELS
# ------------------------------------------------------------------------------
# This translation map is study-specific (Czech NP Bohubin data).
# A German or Slovak dataset will have different source strings.
# Keep it here, not in a shared config.
# ==============================================================================

OBJECT_ID_MAP <- c(
  "3 kláda"                      = "3 log",
  "8 kameny"                     = "8 stones",
  "4 pahýl"                      = "4 stump",
  "5 kládas"                     = "5 logs",
  "7 kameny"                     = "7 stones",
  "6 kláda"                      = "6 log",
  "7 spadlá větev"               = "7 fallen branch",
  "5 kláda"                      = "5 log",
  "6 nízký pahýl"                = "6 low stump",
  "3 pahýl"                      = "3 stump",
  "4 kláda"                      = "4 log",
  "5 pařez"                      = "5 stump",
  "6 pahýl"                      = "6 stump",
  "4 Picea abies mladý"          = "4 Picea abies young",
  "Pinus sylvestris padlá větev" = "Pinus sylvestris fallen branch",
  "6 vývrat"                     = "6 uprooted tree",
  "7 pařez"                      = "7 stump",
  "6 pařez"                      = "6 stump",
  "5 pahýl"                      = "5 stump",
  "5 nízký pahýl"                = "5 low stump",
  "10 nízký pahýl"               = "10 low stump",
  "9 pahýl"                      = "9 stump",
  "8 kláda"                      = "8 log",
  "7 vývrat"                     = "7 uprooted tree",
  "7 nízký pahýl"                = "7 low stump",
  "4 nízký pahýl"                = "4 low stump",
  "1 pahýl"                      = "1 stump",
  "7 pahýl"                      = "7 stump",
  "7 kláda"                      = "7 log",
  "2 kláda"                      = "2 log",
  "4 pařez"                      = "4 stump",
  "8 pařez"                      = "8 stump"
)

lichens <- lichens_raw %>%
  mutate(`object ID` = recode(`object ID`, !!!OBJECT_ID_MAP))

cat("Translation applied to 'object ID'\n")


# ==============================================================================
# 3. SPECIES-LEVEL PREVALENCE SUMMARY
# ------------------------------------------------------------------------------
# Used to make grouping decisions. Printed as a QC report.
# ==============================================================================

species_summary <- lichens %>%
  group_by(species) %>%
  summarise(
    n_plots            = n_distinct(Plot),
    prevalence_pct     = (n_plots / N_PLOTS) * 100,
    total_observations = n(),
    mean_cover         = mean(cover, na.rm = TRUE),
    median_cover       = median(cover, na.rm = TRUE),
    max_cover          = max(cover, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_plots)) %>%
  mutate(
    model_suitability = case_when(
      prevalence_pct >= 20 & prevalence_pct <= 80 ~ "Good - Balanced",
      prevalence_pct >  80                         ~ "Ubiquitous - Low variance",
      prevalence_pct >= 10 & prevalence_pct <  20  ~ "Rare - Consider grouping",
      prevalence_pct <  10                         ~ "Very rare - Not suitable"
    )
  )

cat("\n=== SPECIES PREVALENCE SUMMARY (top 20) ===\n")
print(species_summary, n = 20)

cat("\nModeling suitability breakdown:\n")
print(species_summary %>% count(model_suitability) %>% arrange(desc(n)))

# Save QC report
write_csv(species_summary,
          file.path(PATH_OUT_REPORTS, "qc_species_prevalence.csv"))
cat("  \u2713 Saved: qc_species_prevalence.csv\n")


# ==============================================================================
# 4. PLOT-LEVEL DIVERSITY SUMMARY (QC)
# ==============================================================================

plot_summary <- lichens %>%
  group_by(Plot) %>%
  summarise(
    n_species        = n_distinct(species),
    total_cover      = sum(cover, na.rm = TRUE),
    mean_cover       = mean(cover, na.rm = TRUE),
    n_observations   = n(),
    dominant_species = species[which.max(cover)],
    dominant_cover   = max(cover, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_species))

cat("\n=== PLOT-LEVEL DIVERSITY (top 20) ===\n")
print(plot_summary, n = 20)


# ==============================================================================
# 5. FULL PRESENCE/ABSENCE MATRIX (intermediate object — not saved to disk)
# ==============================================================================

presence_absence_matrix <- lichens %>%
  distinct(Plot, species) %>%
  mutate(presence = 1L) %>%
  pivot_wider(
    names_from  = species,
    values_from = presence,
    values_fill = 0L
  )

cat("\nFull PA matrix dimensions:", dim(presence_absence_matrix), "\n")


# ==============================================================================
# 6. CANDIDATE SPECIES (prevalence 20–80%)
# ==============================================================================

candidate_species <- species_summary %>%
  filter(prevalence_pct >= 20 & prevalence_pct <= 80) %>%
  pull(species)

cat("\n=== CANDIDATE SPECIES FOR DIRECT BINOMIAL MODELING ===\n")
cat("Count:", length(candidate_species), "\n")
print(candidate_species)


# ==============================================================================
# 7. EXPLORATORY FUNCTIONS
# ------------------------------------------------------------------------------
# These functions support interactive exploration during prep.
# They do NOT produce any outputs that feed the modeling pipeline.
# ==============================================================================

# --- 7a. Single-species viability report -------------------------------------

check_species_viability <- function(species_in_question,
                                    species_summary_df = species_summary,
                                    lichen_df          = lichens,
                                    min_prevalence     = 20,
                                    max_prevalence     = 80) {
  #' Report modeling viability for a single lichen species.
  #'
  #' @param species_in_question character(1) — Binomial species name.
  #' @param species_summary_df  data.frame   — Output of species_summary block.
  #' @param lichen_df           data.frame   — Raw lichen records (Plot, species, cover).
  #' @param min_prevalence      numeric(1)   — Lower prevalence threshold (%).
  #' @param max_prevalence      numeric(1)   — Upper prevalence threshold (%).
  #' @return invisible list(species, summary, plot_data), or invisible NULL if not found.
  
  if (!species_in_question %in% species_summary_df$species) {
    cat("\u274c Species NOT FOUND:", species_in_question, "\n")
    possible <- agrep(species_in_question, species_summary_df$species,
                      max.distance = 0.3, value = TRUE)
    if (length(possible) > 0) {
      cat("  Did you mean:\n")
      print(possible)
    }
    return(invisible(NULL))
  }
  
  sp_data        <- filter(species_summary_df, species == species_in_question)
  plot_occurrence <- lichen_df %>%
    filter(species == species_in_question) %>%
    group_by(Plot) %>%
    summarise(
      n_observations = n(),
      mean_cover     = mean(cover, na.rm = TRUE),
      max_cover      = max(cover, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_cover))
  
  n_present <- sp_data$n_plots
  n_absent  <- N_PLOTS - n_present
  ratio     <- round(n_present / n_absent, 2)
  
  cat("\n\u2550\u2550\u2550 SPECIES VIABILITY REPORT \u2550\u2550\u2550\n")
  cat("Species      :", species_in_question, "\n")
  cat("Prevalence   :", n_present, "/", N_PLOTS, "plots (",
      round(sp_data$prevalence_pct, 1), "%)\n")
  cat("Suitability  :", sp_data$model_suitability, "\n")
  cat("P/A balance  :", n_present, "present /", n_absent, "absent | ratio", ratio, "\n")
  cat("Top 5 plots  :\n")
  print(head(plot_occurrence, 5), row.names = FALSE)
  cat("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n")
  
  # Plots
  p1 <- plot_occurrence %>%
    ggplot(aes(x = reorder(Plot, mean_cover), y = mean_cover)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    coord_flip() +
    labs(title = paste("Cover:", species_in_question),
         x = "Plot", y = "Mean cover") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 6))
  
  p2 <- lichen_df %>%
    filter(species == species_in_question) %>%
    ggplot(aes(x = cover)) +
    geom_histogram(binwidth = 0.5, fill = "forestgreen",
                   color = "white", alpha = 0.7) +
    labs(title = "Cover distribution", x = "Cover class", y = "Frequency") +
    theme_minimal()
  
  print(p1)
  print(p2)
  
  invisible(list(species   = species_in_question,
                 summary   = sp_data,
                 plot_data = plot_occurrence))
}


# --- 7b. Batch viability check -----------------------------------------------

check_multiple_species_detailed_nopause <- function(species_list,
                                                    species_summary_df = species_summary,
                                                    lichen_df          = lichens,
                                                    min_prevalence     = 20,
                                                    max_prevalence     = 80) {
  #' Run check_species_viability() for a vector of species names.
  #'
  #' @param species_list     character — Vector of species names to assess.
  #' @param species_summary_df data.frame — Pre-computed species summary.
  #' @param lichen_df        data.frame — Raw lichen records.
  #' @param min_prevalence   numeric(1) — Lower prevalence threshold (%).
  #' @param max_prevalence   numeric(1) — Upper prevalence threshold (%).
  #' @return data.frame with one row per species: name, n_plots, prevalence_pct,
  #'   model_suitability.
  
  cat("\u2554 BATCH VIABILITY CHECK | n =", length(species_list), "\u255d\n")
  
  results_summary <- tibble()
  
  for (i in seq_along(species_list)) {
    sp <- species_list[i]
    cat("\n\u2500\u2500\u2500 Species", i, "of", length(species_list), "\u2500\u2500\u2500\n")
    result <- check_species_viability(
      species_in_question = sp,
      species_summary_df  = species_summary_df,
      lichen_df           = lichen_df,
      min_prevalence      = min_prevalence,
      max_prevalence      = max_prevalence
    )
    if (!is.null(result)) {
      results_summary <- bind_rows(
        results_summary,
        select(result$summary, species, n_plots, prevalence_pct, model_suitability)
      )
    } else {
      results_summary <- bind_rows(
        results_summary,
        tibble(species           = sp,
               n_plots           = NA_integer_,
               prevalence_pct    = NA_real_,
               model_suitability = "NOT FOUND")
      )
    }
  }
  
  cat("\n\u2554 BATCH SUMMARY \u255d\n")
  print(results_summary, n = Inf)
  cat("\nSuitable (20-80%):",
      sum(results_summary$prevalence_pct >= 20 &
            results_summary$prevalence_pct <= 80, na.rm = TRUE), "\n")
  cat("Not found        :", sum(is.na(results_summary$prevalence_pct)), "\n\n")
  
  write_csv(results_summary,
            file.path(PATH_OUT_REPORTS, "qc_batch_species_check.csv"))
  cat("  \u2713 Saved: qc_batch_species_check.csv\n")
  
  results_summary
}


# --- 7c. Co-occurrence analysis ----------------------------------------------

analyze_cooccurrence <- function(species_list, lichen_df = lichens) {
  #' Compute Jaccard similarity, co-occurrence counts, and shared plots for a
  #' set of species.
  #'
  #' @param species_list character — Vector of species names.
  #' @param lichen_df    data.frame — Raw lichen records (Plot, species, cover).
  #' @return list(pa_matrix, jaccard_matrix, cooccur_counts, shared_plots).
  
  pa_matrix <- lichen_df %>%
    filter(species %in% species_list) %>%
    distinct(Plot, species) %>%
    mutate(presence = 1L) %>%
    pivot_wider(names_from = species, values_from = presence, values_fill = 0L)
  
  species_cols  <- select(pa_matrix, -Plot)
  found_species <- intersect(species_list, names(species_cols))
  
  # Jaccard matrix
  jaccard_matrix <- matrix(NA_real_,
                           nrow = length(found_species),
                           ncol = length(found_species),
                           dimnames = list(found_species, found_species))
  
  for (i in seq_along(found_species)) {
    for (j in seq_along(found_species)) {
      s1 <- found_species[i]
      s2 <- found_species[j]
      a  <- sum(species_cols[[s1]] == 1L & species_cols[[s2]] == 1L)
      b  <- sum(species_cols[[s1]] == 1L & species_cols[[s2]] == 0L)
      cc <- sum(species_cols[[s1]] == 0L & species_cols[[s2]] == 1L)
      jaccard_matrix[i, j] <- a / (a + b + cc)
    }
  }
  
  # Pairwise counts
  cooccur_counts <- expand_grid(species1 = found_species,
                                species2 = found_species) %>%
    filter(species1 != species2) %>%
    rowwise() %>%
    mutate(
      both_present = sum(species_cols[[species1]] == 1L &
                           species_cols[[species2]] == 1L),
      sp1_only     = sum(species_cols[[species1]] == 1L &
                           species_cols[[species2]] == 0L),
      sp2_only     = sum(species_cols[[species1]] == 0L &
                           species_cols[[species2]] == 1L),
      neither      = sum(species_cols[[species1]] == 0L &
                           species_cols[[species2]] == 0L),
      jaccard      = both_present / (both_present + sp1_only + sp2_only)
    ) %>%
    ungroup() %>%
    arrange(desc(jaccard))
  
  # Shared plot lists
  shared_plots <- list()
  if (length(found_species) >= 2) {
    for (i in 1:(length(found_species) - 1)) {
      for (j in (i + 1):length(found_species)) {
        s1 <- found_species[i]
        s2 <- found_species[j]
        shared <- pa_matrix %>%
          filter(!!sym(s1) == 1L & !!sym(s2) == 1L) %>%
          pull(Plot)
        if (length(shared) > 0) {
          shared_plots[[paste(s1, "\u00d7", s2)]] <- shared
        }
      }
    }
  }
  
  list(pa_matrix      = pa_matrix,
       jaccard_matrix = jaccard_matrix,
       cooccur_counts = cooccur_counts,
       shared_plots   = shared_plots)
}


# --- 7d. Co-occurrence visualisation -----------------------------------------

plot_cooccurrence <- function(cooccur_analysis) {
  #' Plot a Jaccard heatmap and a tile-plot for pairs with jaccard > 0.3.
  #'
  #' @param cooccur_analysis list — Output of analyze_cooccurrence().
  #' @return invisible NULL. Prints plots as side effect.
  
  jmat <- cooccur_analysis$jaccard_matrix
  valid <- rowSums(is.na(jmat)) < ncol(jmat)
  jmat  <- jmat[valid, valid, drop = FALSE]
  
  if (nrow(jmat) < 2) {
    cat("\u26a0\ufe0f  Not enough species for co-occurrence heatmap\n")
    return(invisible(NULL))
  }
  
  pheatmap(
    jmat,
    cluster_rows    = TRUE, cluster_cols = TRUE,
    display_numbers = TRUE, number_format = "%.2f",
    color           = colorRampPalette(c("white", "yellow", "orange", "red"))(50),
    main            = "Species Co-occurrence (Jaccard Similarity)",
    fontsize = 10, fontsize_number = 8, angle_col = 45
  )
  
  high_pairs <- cooccur_analysis$cooccur_counts %>%
    filter(!is.na(jaccard), jaccard > 0.3) %>%
    select(species1, species2, jaccard, both_present)
  
  if (nrow(high_pairs) > 0) {
    p <- ggplot(high_pairs,
                aes(x = reorder(species1, jaccard),
                    y = reorder(species2, jaccard),
                    fill = jaccard, label = both_present)) +
      geom_tile(color = "white") +
      geom_text(color = "black", size = 3) +
      scale_fill_gradient2(low = "white", mid = "yellow", high = "red",
                           midpoint = 0.5, limits = c(0, 1),
                           name = "Jaccard") +
      labs(title = "Species pairs with Jaccard > 0.3",
           subtitle = "Numbers = shared plot count",
           x = "Species 1", y = "Species 2") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            axis.text.y = element_text(size = 8))
    print(p)
  }
  
  invisible(NULL)
}


# --- 7e. Co-occurrence report ------------------------------------------------

print_cooccurrence_report <- function(cooccur_analysis) {
  #' Print a structured co-occurrence report to the console.
  #'
  #' @param cooccur_analysis list — Output of analyze_cooccurrence().
  #' @return invisible NULL.
  
  cat("\n\u2554 CO-OCCURRENCE REPORT \u255d\n\n")
  
  cat("Top co-occurring pairs (Jaccard: 0=never, 1=always):\n")
  cooccur_analysis$cooccur_counts %>%
    filter(!is.na(jaccard)) %>%
    arrange(desc(jaccard)) %>%
    head(10) %>%
    select(species1, species2, jaccard, both_present) %>%
    print(n = 10)
  
  cat("\nShared plot details:\n")
  if (length(cooccur_analysis$shared_plots) > 0) {
    for (pair in names(cooccur_analysis$shared_plots)) {
      plots <- cooccur_analysis$shared_plots[[pair]]
      cat(" •", pair, ":", length(plots), "plots |",
          paste(plots, collapse = ", "), "\n")
    }
  } else {
    cat("  No species pairs share plots\n")
  }
  
  pa_matrix    <- cooccur_analysis$pa_matrix
  species_cols <- select(pa_matrix, -Plot)
  
  cat("\nIsolation stats (% of occurrences where species is alone):\n")
  isolation <- tibble(
    species          = names(species_cols),
    n_plots_present  = map_int(species_cols, sum),
    n_plots_alone    = map_int(names(species_cols), function(sp) {
      sum(species_cols[[sp]] == 1L &
            rowSums(select(species_cols, -all_of(sp))) == 0L)
    })
  ) %>%
    mutate(pct_alone = round((n_plots_alone / n_plots_present) * 100, 1)) %>%
    arrange(desc(pct_alone))
  print(isolation, n = Inf)
  
  invisible(NULL)
}


# --- 7f. Integrated batch + co-occurrence ------------------------------------

check_multiple_species_with_cooccurrence <- function(species_list,
                                                     species_summary_df = species_summary,
                                                     lichen_df          = lichens,
                                                     min_prevalence     = 20,
                                                     max_prevalence     = 80) {
  #' Run batch viability check followed by co-occurrence analysis.
  #'
  #' @param species_list character — Vector of species to assess.
  #' @return list(species_summary, cooccurrence).
  
  cat("\n\u2554 COMPREHENSIVE SPECIES ANALYSIS \u255d\n")
  
  summary_tbl <- check_multiple_species_detailed_nopause(
    species_list       = species_list,
    species_summary_df = species_summary_df,
    lichen_df          = lichen_df,
    min_prevalence     = min_prevalence,
    max_prevalence     = max_prevalence
  )
  
  found <- summary_tbl %>%
    filter(!is.na(prevalence_pct)) %>%
    pull(species)
  
  if (length(found) < 2) {
    cat("\u26a0\ufe0f  Less than 2 species found — skipping co-occurrence analysis\n")
    return(summary_tbl)
  }
  
  cat("\n\u2500\u2500\u2500 CO-OCCURRENCE ANALYSIS \u2500\u2500\u2500\n")
  cooccur <- analyze_cooccurrence(found, lichen_df)
  print_cooccurrence_report(cooccur)
  plot_cooccurrence(cooccur)
  
  write_csv(cooccur$cooccur_counts,
            file.path(PATH_OUT_REPORTS, "qc_cooccurrence_matrix.csv"))
  cat("  \u2713 Saved: qc_cooccurrence_matrix.csv\n")
  
  list(species_summary = summary_tbl, cooccurrence = cooccur)
}


# ==============================================================================
# 8. RUN EXPLORATORY CHECKS
# ------------------------------------------------------------------------------
# The OGF indicator list from Malíček et al. (2019) is defined here because it
# is a literature reference specific to the Czech / Bohubin study context.
# ==============================================================================

ogf_indicators <- c(
  "Xylographa vitiligo",
  "Chaenotheca sphaerocephala",
  "Parmelia saxatilis",
  "Lecanora subintricata",
  "Mycoblastus sanguinarius",
  "Mycoblastus affinis",
  "Chaenothecopsis pusilla",
  "Calicium viride",
  "Ochrolechia alboflavescens",
  "Chaenothecopsis viridireagens",
  "Bryoria capillaris/nadvornikiana",
  "Micarea globulosella"
)

ogf_results        <- check_multiple_species_with_cooccurrence(ogf_indicators)
ogf_species_info   <- ogf_results$species_summary
ogf_cooccur_data   <- ogf_results$cooccurrence


# ==============================================================================
# 9. DEFINE SPECIES GROUPINGS
# ------------------------------------------------------------------------------
# These groupings are an ecological decision made AFTER examining the data above.
# They are defined here and not in a shared config because:
#   (a) They depend on the species present in THIS dataset.
#   (b) A future study may have completely different species or no viable groups.
# ==============================================================================

SPECIES_GROUPS <- list(
  
  # Group 1: Parmelia aggregate — justified by co-occurrence and Malíček #3
  parmelia_agg = c("Parmelia saxatilis", "Parmelia sulcata"),
  
  # Group 2: Mycoblastus — core OGF specialists
  mycoblastus = c("Mycoblastus sanguinarius", "Mycoblastus affinis"),
  
  # Group 3: Ochrolechia — OGF associates
  ochrolechia = c(
    "Ochrolechia alboflavescens",
    "Ochrolechia microstictoides",
    "Ochrolechia androgyna",
    "Ochrolechia mahluensis"
  ),
  
  # Group 4: Calicioids — lignicolous OGF indicators (count used as richness response)
  calicioids = c(
    "Calicium viride", "Calicium glaucellum", "Calicium trabinellum",
    "Chaenotheca sphaerocephala", "Chaenotheca ferruginea",
    "Chaenotheca chrysocephala", "Chaenotheca xyloxena",
    "Chaenotheca trichialis", "Chaenotheca brunneola",
    "Chaenotheca furfuracea", "Chaenotheca stemonea",
    "Chaenothecopsis pusilla", "Mycocalicium subtile"
  ),
  
  # Group 5: Xylographa — wood substrate specialists
  xylographa = c("Xylographa vitiligo", "Xylographa paralella"),
  
  # Group 6: Elite rare indicators — Malíček top-ranked, low prevalence
  elite_rare = c(
    "Xylographa vitiligo", "Chaenotheca sphaerocephala",
    "Calicium viride", "Chaenothecopsis pusilla",
    "Micarea globulosella"
  ),
  
  # Group 7: Core OGF complex — Mycoblastus + Ochrolechia + Lecanora
  core_ogf = c(
    "Mycoblastus sanguinarius", "Mycoblastus affinis",
    "Ochrolechia alboflavescens", "Lecanora subintricata"
  )
)


# ==============================================================================
# 10. CREATE GROUP PRESENCE / RICHNESS VARIABLES
# ==============================================================================

create_group_variables <- function(lichen_df, species_groups_list) {
  #' Build a plot-level data frame of presence (0/1) and richness (int) for
  #' each species group.
  #'
  #' @param lichen_df          data.frame — Raw lichen records (Plot, species, cover).
  #' @param species_groups_list list      — Named list; each element is a character
  #'   vector of species belonging to that group.
  #' @return data.frame — One row per plot. Columns: Plot, then for each group:
  #'   <group>_presence (integer 0/1), <group>_richness (integer >= 0).
  
  all_plots  <- tibble(Plot = unique(lichen_df$Plot))
  group_data <- all_plots
  
  for (group_name in names(species_groups_list)) {
    spp <- species_groups_list[[group_name]]
    
    pres <- lichen_df %>%
      filter(species %in% spp) %>%
      distinct(Plot) %>%
      mutate(!!paste0(group_name, "_presence") := 1L)
    
    rich <- lichen_df %>%
      filter(species %in% spp) %>%
      group_by(Plot) %>%
      summarise(!!paste0(group_name, "_richness") := n_distinct(species),
                .groups = "drop")
    
    group_data <- group_data %>%
      left_join(pres, by = "Plot") %>%
      left_join(rich, by = "Plot") %>%
      mutate(
        across(ends_with("_presence"), ~ replace_na(., 0L)),
        across(ends_with("_richness"), ~ replace_na(., 0L))
      )
  }
  
  group_data
}

lichen_groups <- create_group_variables(lichens, SPECIES_GROUPS)


# ==============================================================================
# 11. GROUP PREVALENCE REPORT
# ==============================================================================

cat("\n\u2554 SPECIES GROUP PREVALENCE SUMMARY \u255d\n\n")

group_summary <- tibble(
  group       = names(SPECIES_GROUPS),
  n_plots     = map_int(names(SPECIES_GROUPS),
                        ~ sum(lichen_groups[[paste0(.x, "_presence")]])),
  prevalence  = (n_plots / N_PLOTS) * 100,
  max_richness = map_int(names(SPECIES_GROUPS),
                         ~ max(lichen_groups[[paste0(.x, "_richness")]]))
) %>%
  arrange(desc(prevalence))

for (i in seq_len(nrow(group_summary))) {
  g    <- group_summary$group[i]
  n    <- group_summary$n_plots[i]
  prev <- group_summary$prevalence[i]
  cat(sprintf("  %-22s  %3d plots (%5.1f%%)  max richness: %d\n",
              toupper(g), n, prev, group_summary$max_richness[i]))
  if (prev >= 20 & prev <= 80) cat("    \u2705 Suitable for binomial GLM\n")
  else if (prev >= 10)         cat("    \u26a0\ufe0f  Marginal — consider richness model\n")
  else if (prev < 10)          cat("    \u274c Too rare — use MaxEnt or exclude\n")
  else                         cat("    \u26a0\ufe0f  Ubiquitous — limited variation\n")
  cat("\n")
}

# Save group summary
write_csv(group_summary,
          file.path(PATH_OUT_REPORTS, "qc_group_prevalence.csv"))
cat("  \u2713 Saved: qc_group_prevalence.csv\n")

# Prevalence barplot
ggplot(group_summary,
       aes(x = reorder(group, prevalence), y = prevalence,
           fill = cut(prevalence, breaks = c(0, 10, 20, 80, 100)))) +
  geom_col(alpha = 0.85) +
  geom_hline(yintercept = c(20, 80), linetype = "dashed", color = "red") +
  geom_text(aes(label = n_plots), hjust = -0.3, size = 4) +
  scale_fill_manual(
    values = c("(0,10]" = "red", "(10,20]" = "orange",
               "(20,80]" = "darkgreen", "(80,100]" = "gold"),
    labels = c("Too rare", "Marginal", "Suitable", "Ubiquitous"),
    name   = "Suitability"
  ) +
  coord_flip() +
  labs(title    = "Species Group Prevalence",
       subtitle = paste0("Red dashed lines = 20-80% threshold | n = ", N_PLOTS, " plots"),
       x = "Group", y = "Prevalence (%)") +
  theme_minimal()


# ==============================================================================
# 12. LOAD COORDINATES AND ATTACH TO LICHEN GROUPS
# ==============================================================================

coords_raw <- read_excel(PATH_COORDS_RAW)
colnames(coords_raw) <- trimws(colnames(coords_raw))

cat("\nCoordinate file:", nrow(coords_raw), "rows,",
    ncol(coords_raw), "cols\n")
cat("Columns:", paste(colnames(coords_raw), collapse = ", "), "\n")

# Rename the Czech column headers — study-specific
coords_clean <- coords_raw %>%
  rename(
    plot_id     = č.,
    coordinates = souradnice
  ) %>%
  select(plot_id, X, Y)

cat("Coordinate ranges:\n")
cat("  X:", range(coords_clean$X, na.rm = TRUE), "\n")
cat("  Y:", range(coords_clean$Y, na.rm = TRUE), "\n")


# ==============================================================================
# 13. BUILD CANONICAL lichen_clean.csv
# ------------------------------------------------------------------------------
# Schema (per contract defined in Step 1):
#   plot_id  | X  | Y  | <group>_presence (int 0/1) | calicioids_richness (int)
#
# IMPORTANT: Only the response columns that will be modeled are kept.
#   - All _richness columns except calicioids_richness are dropped because the
#     count distribution for other groups is not suitable for count models
#     in this study. Future studies should review this decision.
#   - mycoblastus_presence is retained but flagged: low prevalence means it
#     may be excluded during modeling (decision recorded in study_config.R).
# ==============================================================================

# Map raw 'Plot' to canonical 'plot_id'
lichen_with_id <- lichen_groups %>%
  rename(plot_id = Plot)

# Verify all plots have coordinates
missing_coords <- setdiff(lichen_with_id$plot_id, coords_clean$plot_id)
if (length(missing_coords) > 0) {
  warning("These plots have no coordinates and will be DROPPED: ",
          paste(missing_coords, collapse = ", "))
}

# Join coordinates
lichen_clean <- lichen_with_id %>%
  inner_join(coords_clean, by = "plot_id") %>%
  select(
    plot_id,
    X,
    Y,
    # Count response
    calicioids_richness,
    # Binary responses — canonical _presence columns
    parmelia_agg_presence,
    ochrolechia_presence,
    core_ogf_presence,
    mycoblastus_presence,    # Low prevalence — flag in study_config.R
    xylographa_presence,
    elite_rare_presence
  ) %>%
  # Enforce integer type for all response columns (contract requirement)
  mutate(
    across(ends_with("_presence"), as.integer),
    calicioids_richness = as.integer(calicioids_richness)
  )

# Validate: no NAs in any response column
na_check <- lichen_clean %>%
  summarise(across(c(calicioids_richness, ends_with("_presence")),
                   ~ sum(is.na(.)))) %>%
  pivot_longer(everything())
if (any(na_check$value > 0)) {
  stop("NA values in response columns — check data:\n",
       paste(filter(na_check, value > 0)$name, collapse = ", "))
}

# Validate: plot_id is unique
if (anyDuplicated(lichen_clean$plot_id)) {
  stop("Duplicate plot_id values in lichen_clean — check create_group_variables()")
}

cat("\n\u2550\u2550\u2550 lichen_clean SUMMARY \u2550\u2550\u2550\n")
cat("Dimensions:", nrow(lichen_clean), "plots x", ncol(lichen_clean), "columns\n")
cat("Columns   :", paste(names(lichen_clean), collapse = ", "), "\n")
cat("NA check  : all zero\n")

# Save canonical output
write_csv(lichen_clean, PATH_OUT_CLEAN)
cat(sprintf("  \u2713 Saved: %s  [%d rows x %d cols]\n",
            PATH_OUT_CLEAN, nrow(lichen_clean), ncol(lichen_clean)))

cat("\n\u2705 01a_lichen_prep.R complete.\n")
cat("   Next: run 01b_plot_prep.R to produce structure_clean.csv\n")
cat("   Then: run 02_glm_setup.R to merge and begin modeling\n")