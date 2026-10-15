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

suppressMessages({
  library(dplyr)
  library(readr)
})

MIN_DEPTH          <- 5
DELETION_P_CUTOFF  <- 0.01
MAX_GEN            <- 60000
DATA_DIR           <- "data"
OUT_DIR            <- "src/ShinyAFTrajectories/output"

populations <- c("m1", "m2", "m3", "m4", "m5", "m6",
                 "p1", "p2", "p3", "p4", "p5", "p6")

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

  out <- bind_cols(muts %>% select(Gene, Allele, Annotation, Position),
                   as_tibble(round(af, 4)))
  write_csv(out, file.path(OUT_DIR, paste0("data_", pop, ".csv")), na = "")

  cat(sprintf("%s: %d mutations, %d trimmed for deletions, %.1f%% of points blank\n",
              pop, nrow(out), length(deleted), 100 * mean(is.na(af))))
}

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
for (pop in populations) build_population(pop)
