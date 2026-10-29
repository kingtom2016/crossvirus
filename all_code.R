#!/usr/bin/env Rscript
# IMG/VR cross-habitat analysis for Supplementary Figures S1–S2.
# Single file: no cache, no genome re-annotation, no intermediate RDS files.
# Usage: Rscript run_analysis.R IMGVR_metadata.tsv viral_core_features.tsv results

.libPaths(c("/mnt/h/agent_analysis/test/IMGVR_crosshabitat/Rlib", .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(rstatix) 
  library(jintao)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript run_analysis.R <IMGVR_metadata.tsv> <core_features.tsv> <output_dir>")
}
meta_file <- args[1L]
feature_file <- args[2L]
out <- args[3L]
for (f in c(meta_file, feature_file)) {
  if (!file.exists(f)) stop("Input file not found: ", f)
}
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# These thresholds reproduce the curated4_ge4 analysis settings (config.R).
MIN_GENOMES <- 4L
MIN_COMPLETENESS <- 50
MAX_TEST_PER_GROUP <- 20000L  # original Wilcoxon sample cap; see README
SEED <- 20261019L

# 1. Read IMG/VR metadata. Habitat breadth uses ALL eligible genome records,
#    regardless of completeness or presence of annotated core-genome features.
meta_cols <- c("UVIG", "vOTU", "Ecosystem classification",
               "Length", "Topology", "Estimated completeness")
d <- fread(meta_file, select = meta_cols, showProgress = TRUE)
setnames(d, old = meta_cols,
         new = c("genome", "votu", "ecosystem", "size", "topology", "compl"))
d[, `:=`(size = as.numeric(size), compl = as.numeric(compl))]
d <- d[!is.na(ecosystem) & !is.na(genome) & nzchar(genome) & !is.na(votu) &
       nzchar(votu) & votu != "vOTU_00217951"]

# Habitat mapping: EXACTLY the four keys of config.R scheme_defs$curated4_ge4$map.
# Every other level-1;level-2 combination is out of the analysis, as in apply_scheme()
# (none of the four keys is removed by the row-level filters, so the mapping is equivalent).
ns <- nchar(gsub("[^;]", "", d$ecosystem))
d[, `:=`(
  e1 = sub(";.*$", "", ecosystem),
  e2 = ifelse(ns >= 1, sub("^[^;]*;([^;]*).*$", "\\1", ecosystem), NA_character_),
  e3 = ifelse(ns >= 2, sub("^[^;]*;[^;]*;([^;]*).*$", "\\1", ecosystem), NA_character_),
  e4 = ifelse(ns >= 3, sub("^[^;]*;[^;]*;[^;]*;([^;]*).*$", "\\1", ecosystem), NA_character_)
)]
# Level-2 recode, verbatim from 00_build_slim.R: the habitat labels the four-class scheme
# is built on are the RECODED ones (e.g. "Mammals: Human" -> Human, "Mammals" -> Animals),
# and Soil is re-tagged only for records whose level-3 label really is Soil.
e2v <- d$e2
e2_map <- data.table::fcase(
  e2v %in% c("Human", "Mammals: Human"), "Human",
  # Water = the water column only; Fish/Algae/Plankton are host-associated, not aquatic habitats
  e2v == "Aquatic", "Aquatic",
  grepl("Arthropoda", e2v, fixed = TRUE), "Arthropoda",
  e2v %in% c("Mammals", "Birds", "Reptilia", "Amphibia", "Mollusca", "Porifera",
             "Annelida", "Tunicates", "Cnidaria", "Cephalochordata", "Nematoda",
             "Invertebrates", "Fish"), "Animals",
  e2v %in% c("Protists", "Protozoa", "Amoebozoa"), "Protists",
  e2v %in% c("Lab enrichment", "Lab synthesis", "Laboratory developed", "Lab culture"), "Laboratory",
  e2v %in% c("Wastewater", "WWTP", "Sewage treatment plant"), "Wastewater treatment",
  e2v %in% c("Food production", "Animal feed production", "Feedstock"), "Food production",
  e2v %in% c("", "Unclassified") | is.na(e2v), "Unclassified",
  default = NA_character_
)
d[, e2 := data.table::fifelse(is.na(e2_map), e2v, e2_map)]
soil_idx <- !is.na(d$e1) & d$e1 == "Environmental" & !is.na(d$e3) & d$e3 == "Soil"
d[soil_idx, e2 := "Soil"]
# Sediment (level 3 or 4), Unclassified and Aquaculture are not water-column samples.
drop_water <- !is.na(d$e1) & d$e1 == "Environmental" & !is.na(d$e2) & d$e2 == "Aquatic" &
              ((!is.na(d$e3) & d$e3 %in% c("Sediment", "Unclassified", "Aquaculture")) |
               (!is.na(d$e4) & d$e4 == "Sediment"))
cat("water records dropped (sediment / unclassified / aquaculture):", sum(drop_water), "\n")
d[drop_water, e2 := "ExcludedWater"]
habitat_map <- c(
  "Host-associated;Animals" = "Animals",
  "Host-associated;Human"   = "Human",
  "Environmental;Aquatic"   = "Water",
  "Environmental;Soil"      = "Soil"
)
d[, habitat_raw := paste0(e1, ";", e2)]
d <- d[habitat_raw %in% names(habitat_map)]
d[, habitat := unname(habitat_map[habitat_raw])]
d <- d[!is.na(habitat), .(genome, votu, habitat, size, topology, compl)]
if (anyDuplicated(d$genome)) {
  stop("Duplicate genome IDs after habitat selection: check source metadata; do not double-count.")
}
breadth <- d[, .(n_hab = uniqueN(habitat), n_genomes = .N), by = votu]
breadth <- breadth[n_genomes >= MIN_GENOMES]
d <- d[votu %in% breadth$votu]
cat("vOTUs with at least", MIN_GENOMES, "genomes:", nrow(breadth), "\n")

# 2. Join the precomputed gene features (from the SAME IMG/VR genome records).
#    The merged feature TSV should include the previous Prodigal-GV annotations.
feature_cols <- c("genome_id", "GC_percent", "ORF_count",
                  "ORF_density_per_kb", "coding_density_pct")
core <- fread(feature_file, select = feature_cols, showProgress = TRUE)
core[, genome := sub("\\|.*$", "", genome_id)]
n_dup <- sum(duplicated(core$genome))
if (n_dup) warning(n_dup, " duplicated core-feature IDs: retaining first occurrence (original rule).")
core[, genome_id := NULL]
core <- unique(core, by = "genome")
d <- merge(d, core, by = "genome", all.x = TRUE, sort = FALSE)
rm(core)

# 3. One annotated >=50%-complete representative per vOTU for S2 panels A–E.
#    Tie-breaker: longer assembled genome. This also determines the cohort for S1.
candidates <- d[is.finite(compl) & compl >= MIN_COMPLETENESS &
                is.finite(GC_percent)]
setorder(candidates, votu, -compl, -size)
rep <- candidates[, .SD[1L], by = votu]
rm(candidates)
if (!nrow(rep)) stop("No eligible annotated representative genomes. Check feature input and IDs.")

# Panel F: per-vOTU provirus fraction over EVERY genome of the population (no completeness
# filter on individual genomes), i.e. the same quantity as provirus_pct in the pipeline's
# breadth4_votu_host_info.csv.
prov <- d[, .(n_all = .N, n_provirus = sum(topology == "Provirus", na.rm = TRUE)), by = votu]
cohort <- merge(
  rep[, .(votu, representative_genome = genome,
          rep_completeness = compl,
          genome_length = size / (compl / 100),
          GC_percent, ORF_count = ORF_count / (compl / 100),
          ORF_density_per_kb, coding_density_pct)],
  breadth, by = "votu", all.x = TRUE
)
cohort <- merge(cohort, prov, by = "votu", all.x = TRUE)
cohort[, pct_provirus := 100 * n_provirus / n_all]
setorder(cohort, n_hab, votu)
stopifnot(all(cohort$n_hab %in% 1:4), all(cohort$n_all >= 1L))
cat("S1/S2 cohort:", nrow(cohort), "vOTUs\n")
print(cohort[, .(vOTUs = .N), by = n_hab][order(n_hab)])

# Cohort check:
panel_file <- "/mnt/h/agent_analysis/test/IMGVR_crosshabitat/results/panel_counts_curated4_ge4.csv"
if (file.exists(panel_file)) {
  pc <- fread(panel_file)
  got <- cohort[, .(n = .N), by = n_hab][order(n_hab)]
  if (sum(got$n) != sum(pc$n_votu) || any(got$n != pc$n_votu[order(pc$bin)])) {
    warning("Cohort differs from results/panel_counts_curated4_ge4.csv: this run has ",
            paste(sprintf("%d=%d", got$n_hab, got$n), collapse = " "),
            " but the pipeline reports ",
            paste(sprintf("%d=%d", pc$bin, pc$n_votu), collapse = " "),
            ". Check metadata and core-feature table versions.")
  }
}
fwrite(cohort, file.path(out, "votu_analysis_table.csv"))

metrics <- list(
  list(key = "genome_length",      y = "Completeness-corrected viral genome length (bp)"),
  list(key = "GC_percent",         y = "GC (%)"),
  list(key = "ORF_count",          y = "ORFs per genome"),
  list(key = "ORF_density_per_kb", y = "ORFs per kb"),
  list(key = "coding_density_pct", y = "Coding density (%)"),
  list(key = "pct_provirus",       y = "Provirus genomes (%)")
)
set.seed(SEED)
test_tables <- list()
letter_tables <- list()
for (m in metrics) {
  x <- cohort[is.finite(get(m$key)),
              .(value = get(m$key), group = as.character(n_hab))]
  x <- x[, .SD[if (.N > MAX_TEST_PER_GROUP) sample.int(.N, MAX_TEST_PER_GROUP)
                else seq_len(.N)], by = group]
  levels_present <- sort(unique(x$group))
  if (length(levels_present) < 2L) stop("At least two habitat-breadth groups are needed.")
  st <- as.data.table(jintao::statistic_letter_singlelevel(
          x, "value", "group", use_adjust_p = TRUE, wilcox = TRUE, output_statis = TRUE))
  lt <- as.data.table(jintao::statistic_letter_singlelevel(
          x, "value", "group", use_adjust_p = TRUE, wilcox = TRUE))
  lcol <- grep("^letter$", names(lt), ignore.case = TRUE)[1]
  if (is.na(lcol)) lcol <- grep("letter", names(lt), ignore.case = TRUE)[1]
  gcol <- if ("group" %in% names(lt)) "group" else names(lt)[1]
  test_tables[[m$key]] <- st[, .(metric = m$key,
                                 group_1 = group1_ori1, group_2 = group1_ori2,
                                 p_value = p, p_BH = p.adj)]
  letter_tables[[m$key]] <- data.table(metric = m$key,
                                       n_hab = as.integer(lt[[gcol]]),
                                       letter = lt[[lcol]])
  cat(sprintf("  %-20s letters: %s\n", m$key,
      paste(sprintf("%s=%s", letter_tables[[m$key]]$n_hab, letter_tables[[m$key]]$letter),
            collapse = " ")))
}
tests <- rbindlist(test_tables)
letters <- rbindlist(letter_tables)
fwrite(tests, file.path(out, "wilcoxon_BH.csv"))
fwrite(letters, file.path(out, "significance_letters.csv"))

# 5. Exploratory adjusted associations 
logistic <- glm(I(n_hab > 1L) ~ log10(n_genomes),
                data = cohort, family = binomial())
capture_coefficients <- function(fit, model_name) {
  s <- summary(fit)$coefficients
  data.table(model = model_name, term = rownames(s),
             estimate = s[, 1L], std_error = s[, 2L],
             statistic = s[, 3L], p_value = s[, 4L],
             n = nobs(fit))
}
models <- list(capture_coefficients(logistic, "multi_habitat_logistic"))
for (m in metrics) {
  md <- data.frame(y = cohort[[m$key]],
                   breadth = factor(cohort$n_hab, levels = 1:4),
                   log_genomes = log10(cohort$n_genomes),
                   completeness = cohort$rep_completeness)
  fit <- lm(y ~ breadth + log_genomes + completeness, data = md)
  models[[length(models) + 1L]] <- capture_coefficients(fit, paste0(m$key, "_adjusted_lm"))
}
fwrite(rbindlist(models), file.path(out, "adjusted_models.csv"))

# 6. Figure S1
counts <- cohort[, .N, by = n_hab][order(n_hab)]
p1 <- ggplot(counts, aes(x = n_hab, y = N)) +
  geom_col(width = 0.65, fill = "white", colour = "black", linewidth = 0.5) +
  scale_x_continuous(breaks = 1:4) +
  scale_y_log10(breaks = 10^(0:ceiling(log10(max(counts$N))))) +
  labs(x = "Number of observed habitat categories per vOTU",
       y = "Number of vOTUs (log10 scale)") +
  theme_classic(base_size = 13)
ggsave(file.path(out, "Figure_S1.png"), p1, width = 5.5, height = 5.2, dpi = 300)
ggsave(file.path(out, "Figure_S1.pdf"), p1, width = 5.5, height = 5.2)

# 7. Figure S2
cohort[, group := factor(n_hab, levels = 1:4)]
small_group <- cohort[, .N, by = n_hab][N <= 60L, n_hab]
make_panel <- function(m) {
  lt <- letters[metric == m$key]
  if (m$key == "pct_provirus") {
    st <- cohort[, .(
      mean_p = mean(pct_provirus),
      se = sd(pct_provirus) / sqrt(.N),
      weighted_p = 100 * sum(n_provirus) / sum(n_all)
    ), by = n_hab]
    st[, `:=`(group = factor(n_hab, levels = 1:4),
              letter = lt$letter[match(n_hab, lt$n_hab)])]
    st[, label_y := mean_p + se + max(2, 0.05 * max(mean_p + se))]
    return(ggplot(st, aes(group, mean_p)) +
      geom_col(fill = "grey75", colour = "black", width = 0.65) +
      geom_errorbar(aes(ymin = pmax(0, mean_p - se), ymax = mean_p + se), width = 0.17) +
      geom_point(aes(y = weighted_p), shape = 23, fill = "white", size = 2.8) +
      geom_text(aes(y = label_y, label = letter), size = 4.5) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
      labs(x = "Habitat categories per vOTU", y = "Provirus (%) — mean ± SE") +
      theme_classic(base_size = 12))
  }
  position <- cohort[is.finite(get(m$key)),
                     .(y = max(get(m$key))), by = n_hab]
  yy <- cohort[[m$key]]
  span <- diff(range(yy[is.finite(yy)]))
  offset <- if (span > 0) 0.06 * span else 0.1 * max(1, abs(yy[1L]))
  position[, `:=`(
    group = factor(n_hab, levels = 1:4),
    y = y + offset,
    letter = lt$letter[match(n_hab, lt$n_hab)]
  )]
  ggplot(cohort, aes(x = group, y = .data[[m$key]])) +
    geom_violin(fill = "grey75", trim = TRUE, scale = "width", na.rm = TRUE) +
    geom_boxplot(width = 0.13, fill = "white", outlier.shape = NA, na.rm = TRUE) +
    geom_jitter(data = cohort[n_hab %in% small_group], width = 0.06,
                size = 0.65, alpha = 0.6, na.rm = TRUE) +
    geom_text(data = position, aes(x = group, y = y, label = letter),
              inherit.aes = FALSE, size = 4.5) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.14))) +
    labs(x = "Habitat categories per vOTU", y = m$y) +
    theme_classic(base_size = 12)
}
panels <- lapply(metrics, make_panel)
p2 <- patchwork::wrap_plots(panels, ncol = 3) +
  patchwork::plot_annotation(tag_levels = "A")
ggsave(file.path(out, "Figure_S2.png"), p2, width = 13, height = 8, dpi = 300)
ggsave(file.path(out, "Figure_S2.pdf"), p2, width = 13, height = 8)

cat("Done. Results written to:", normalizePath(out), "\n")
