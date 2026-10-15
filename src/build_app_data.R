# Builds the CSVs the Shiny app reads (src/ShinyAFTrajectories/output/data_<pop>.csv)
# from the raw annotated timecourse files in data/.
#
# Run from the project root:
#   Rscript src/build_app_data.R
#
# Masking follows Good et al. 2017 (LTEE-metagenomic-master/parse_file.py
# parse_annotated_timecourse and timecourse_utils.py mask_timepoints), so the
# trajectories match their figures. A frequency is left blank (NA) when:
#   1. the mutation's read depth at that timepoint is < MIN_DEPTH
#   2. the population's average depth at that timepoint (the "Depth" row,
#      first row of the file) is < MIN_DEPTH
#   3. Deletion P-value < 0.01 and the timepoint is at or after Deletion index,
#      i.e. the region was deleted from the genome and later reads are noise
#
# Rows marked FAIL are dropped, and only generations <= 60,000 are kept
# (larger numbers encode clone samples).
#
# CloneAllele column: Good et al. merge every indel variant seen at a site into
# one Allele string (e.g. "indel;+G;+GG;-1;+GGGGGGG"), most of which are read
# errors. For indel rows, this looks up breseq INS/DEL calls from sequenced
# clones of the same population (data/LTEE-clone-curated.tsv) within
# CLONE_WINDOW bp. A clone is a single genome, so its call is the real mutation.
# Matching calls are written as plain text, e.g. "(G)₅→₆", joined with " / ".
# The cell is left blank when no clone call matches.

suppressMessages({
  library(dplyr)
  library(readr)
})

MIN_DEPTH          <- 5
DELETION_P_CUTOFF  <- 0.01
MAX_GEN            <- 60000
CLONE_WINDOW       <- 10
DATA_DIR           <- "data"
OUT_DIR            <- "src/ShinyAFTrajectories/output"

populations <- c("m1", "m2", "m3", "m4", "m5", "m6",
                 "p1", "p2", "p3", "p4", "p5", "p6")

population_labels <- c(m1 = "Ara-1", m2 = "Ara-2", m3 = "Ara-3",
                       m4 = "Ara-4", m5 = "Ara-5", m6 = "Ara-6",
                       p1 = "Ara+1", p2 = "Ara+2", p3 = "Ara+3",
                       p4 = "Ara+4", p5 = "Ara+5", p6 = "Ara+6")

# breseq's HTML mutation text -> plain text, e.g.
#   "(G)<sub>5&rarr;6</sub>"   -> "(G)₅→₆"
#   "&Delta;1,199&nbsp;bp"     -> "Δ1,199 bp"
breseq_html_to_text <- function(html) {
  subscript <- function(s) chartr("0123456789→", "₀₁₂₃₄₅₆₇₈₉→", s)
  out <- gsub("&rarr;", "→", html, fixed = TRUE)
  out <- gsub("&Delta;", "Δ", out, fixed = TRUE)
  out <- gsub("&nbsp;", " ", out, fixed = TRUE)
  out <- gsub("&#8211;", "–", out, fixed = TRUE)
  for (i in seq_along(out)) {
    m <- regmatches(out[i], gregexpr("<sub>[^<]*</sub>", out[i]))[[1]]
    for (s in m) out[i] <- sub(s, subscript(gsub("</?sub>", "", s)), out[i], fixed = TRUE)
  }
  gsub("<[^>]+>", "", out)
}

clone_calls <- read_tsv(file.path(DATA_DIR, "LTEE-clone-curated.tsv"),
                        show_col_types = FALSE, progress = FALSE) %>%
  filter(treatment == "LTEE", type %in% c("INS", "DEL")) %>%
  distinct(population, start_position, html_mutation) %>%
  mutate(text = breseq_html_to_text(html_mutation)) %>%
  arrange(population, start_position, text)

clone_allele <- function(pop, position, allele) {
  calls <- clone_calls %>% filter(population == population_labels[[pop]])
  vapply(seq_along(position), function(i) {
    if (!grepl("^indel;", allele[i])) return(NA_character_)
    hits <- unique(calls$text[abs(calls$start_position - position[i]) <= CLONE_WINDOW])
    if (length(hits) == 0) NA_character_ else paste(hits, collapse = " / ")
  }, character(1))
}

build_population <- function(pop) {
  raw <- read_csv(file.path(DATA_DIR, paste0(pop, "_annotated_timecourse.txt")),
                  show_col_types = FALSE, progress = FALSE)

  # Population timepoints, in file order. "Deletion index" counts these.
  gens <- sub("^AC:", "", grep("^AC:", names(raw), value = TRUE))
  gens <- gens[as.numeric(gens) < 1000000]

  avg_depth <- unlist(raw[1, paste0("DP:", gens)])
  muts      <- raw %>% filter(`Passed?` != "FAIL")

  ac <- as.matrix(muts[, paste0("AC:", gens)])
  dp <- as.matrix(muts[, paste0("DP:", gens)])

  keep <- dp >= MIN_DEPTH                       # rule 1
  keep[, avg_depth < MIN_DEPTH] <- FALSE        # rule 2
  deleted <- which(muts$`Deletion P-value` < DELETION_P_CUTOFF)
  for (i in deleted) {                          # rule 3
    # Deletion index is 0-based in Good et al.'s code.
    first_masked <- muts$`Deletion index`[i] + 1
    if (first_masked <= length(gens)) keep[i, first_masked:length(gens)] <- FALSE
  }

  af <- ac / dp
  af[!keep] <- NA
  colnames(af) <- paste0("AF:", gens)
  af <- af[, as.numeric(gens) <= MAX_GEN, drop = FALSE]

  meta <- muts %>%
    select(Gene, Allele, Annotation, Position) %>%
    mutate(CloneAllele = clone_allele(pop, Position, Allele))

  out <- bind_cols(meta, as_tibble(round(af, 4)))
  write_csv(out, file.path(OUT_DIR, paste0("data_", pop, ".csv")), na = "")

  n_indel <- sum(grepl("^indel;", meta$Allele))
  cat(sprintf("%s: %d mutations, %d trimmed for deletions, %.1f%% of points blank, %d/%d indels matched to clones\n",
              pop, nrow(out), length(deleted), 100 * mean(is.na(af)),
              sum(!is.na(meta$CloneAllele)), n_indel))
}

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
for (pop in populations) build_population(pop)
