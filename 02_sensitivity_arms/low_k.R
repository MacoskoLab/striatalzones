#!/usr/bin/env Rscript
# ==============================================================================
# low_k.R
#
# The fixed-neighbour profile peaked at k = 5 for MSN, Astrocyte and OPC
# receivers — the LOWEST value tested — so the peak (if there is one) sits below
# the range we ran. This extends downward: k = 1, 3, and 5.
#
#   k = 1  the single nearest partner: the most contact-like definition possible
#   k = 3  between that and the existing grid
#   k = 5  recomputed as a bridge, to measure the Monte-Carlo offset against the
#          existing k=5 numbers (the permutation seed differs between runs)
#
# All six senders, the same fixed candidate set, 200 shuffles — so the counts
# are directly comparable to the existing k = 5/10/20/40 results.
#
# Only the NULL is needed: run_targeted_null_for_donor() computes the real
# effect and its null together at each grid cell.
#
# Run inside RStudio.  Output: pipeline_v3/08_low_k/
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dbscan); library(qs)
  library(dplyr); library(tibble); library(readr); library(tidyr); library(stringr)
})

CODE_DIR <- "/broad/macosko/akraft/CAP/rerun/processed_pucks/aggregated/NoReg_Final_by_donor/Final_Cell_Arrangements_Oct9/Fig6"
PIPE_OUT <- file.path(CODE_DIR, "LR_sensitivity_analysis", "pipeline_v3")
OUT_DIR  <- file.path(PIPE_OUT, "08_low_k")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

KSTR <- Sys.getenv("LR_LOW_K", "1,3,5")
CFG <- list(ks = as.numeric(strsplit(KSTR, ",")[[1]]), trim = 0.1, n_shuffles = 200, seed = 42)
KTAG <- gsub(",", "-", KSTR)

logf <- file.path(OUT_DIR, "low_k.log")
say <- function(...) { m <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
                       cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }

say("=== LOW-k SWEEP ===  k = ", paste(CFG$ks, collapse = ", "),
    " | trim ", CFG$trim, " | ", CFG$n_shuffles, " shuffles | all senders")

tgt <- readr::read_csv(file.path(PIPE_OUT, "03_null", "target_interactions_for_null.csv"),
                       show_col_types = FALSE)
say("candidate interactions: ", nrow(tgt))

grid <- tibble(Method = "knn", Scale = as.numeric(CFG$ks), Trim = CFG$trim)

quiet <- list(log_msg = function(..., level = "INFO") invisible(NULL))
cellchat_db <- load_cellchatdb_human(PIPE_OUT, quiet)
db <- prepare_interaction_db(cellchat_db, c("Secreted Signaling", "Cell-Cell Contact"),
                             OUT_DIR, quiet)

donors <- sort(setdiff(names(donor_objs), "s20"))
say("donors: ", length(donors))

for (di in seq_along(donors)) {
  d_name <- donors[di]
  fout <- file.path(OUT_DIR, sprintf("%s_lowk_%s.csv", d_name, KTAG))
  if (file.exists(fout)) { say("skip ", d_name, " (done)"); next }
  say("donor ", d_name, " (", di, "/", length(donors), ")")

  obj <- if (is.character(donor_objs[[d_name]])) qs::qread(donor_objs[[d_name]]) else donor_objs[[d_name]]
  d_data <- prep_donor_data(obj, db$required_genes,
                            c("MSN","Astrocyte","Oligo","OPC","Microglia","Vascular"),
                            1:6, quiet, d_name)
  rm(obj); invisible(gc(verbose = FALSE))
  if (is.null(d_data)) { say("  skipped"); next }

  donor_db <- filter_interactions_per_donor(d_data, db$interaction, db$entity_to_genes,
                                            0.05, quiet, d_name)
  tgt_d <- dplyr::semi_join(tgt, donor_db, by = c("Ligand" = "ligand", "Receptor" = "receptor"))
  if (nrow(tgt_d) == 0) { say("  no candidates pass detection"); next }

  res <- run_targeted_null_for_donor(
    d_name = d_name, d_data = d_data, target_df = tgt_d, complex_db = db$complex,
    scale_grid = grid, radii_um = numeric(0), knn_k = CFG$ks, trim = NULL,
    n_shuffles = CFG$n_shuffles, shuffle_seed = CFG$seed, donor_index = di,
    zone_a = 1, zone_b = 4, min_senders_per_zone = 5, logger = quiet)

  readr::write_csv(res, fout)
  say("  wrote ", basename(fout), " rows=", nrow(res))
  rm(d_data, res); invisible(gc(verbose = FALSE))
}

fs <- list.files(OUT_DIR, pattern = "_lowk.*\\.csv$", full.names = TRUE)
all <- bind_rows(lapply(fs, readr::read_csv, show_col_types = FALSE)) %>%
  distinct(Donor, Method, Scale, Sender, Receiver, Interaction, .keep_all = TRUE)
readr::write_csv(all, file.path(OUT_DIR, "lowk_all_donors.csv"))

summ <- all %>%
  group_by(Method, Scale, Sender, Receiver, Interaction) %>%
  summarise(Donors = n(),
            Robust_95 = mean(abs(Real_Delta) > Null_95_Abs, na.rm = TRUE),
            .groups = "drop")
readr::write_csv(summ, file.path(OUT_DIR, "lowk_robustness.csv"))

say("ROBUST COUNTS (Robust_95 >= 0.5):")
tot <- summ %>% filter(Robust_95 >= 0.5) %>% count(Scale)
for (i in seq_len(nrow(tot))) say("   ALL      k=", tot$Scale[i], " : ", tot$n[i])
bys <- summ %>% filter(Robust_95 >= 0.5) %>% count(Sender, Scale)
for (i in seq_len(nrow(bys))) say("   ", bys$Sender[i], "  k=", bys$Scale[i], " : ", bys$n[i])
byr <- summ %>% filter(Robust_95 >= 0.5) %>% count(Receiver, Scale)
for (i in seq_len(nrow(byr))) say("   ->", byr$Receiver[i], "  k=", byr$Scale[i], " : ", byr$n[i])
say("DONE")
