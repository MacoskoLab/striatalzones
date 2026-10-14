#!/usr/bin/env Rscript
# ==============================================================================
# matched_radii.R
#
# Density-matched radius arm.  Each sender's kNN saturation plateaus at a
# characteristic k* (MSN 2, Astrocyte 2, Oligo 3, OPC 4, Vascular 4, Microglia 5).
# The physical distance to that k*-th partner, rounded to the nearest 100 um, is
# that sender's "density-matched" radius:
#
#     MSN 200 | Oligo 200 | Astrocyte 300 | OPC 500 | Vascular 500 | Microglia 700
#
# This runs the radius arm at the UNION of those radii {200, 300, 500, 700} for
# ALL senders, on the SAME fixed candidate set and permutation null as the main
# run, so counts are directly comparable.  200 is recomputed here as a bridge /
# consistency check against the main pipeline's 200 um numbers.
#
# Because all senders are kept in the target set, the per-sender shuffle indexing
# matches the main run (unlike extend_radii.R, which restricted senders and so
# drew an independent null).  Each sender is later read off at its own matched
# radius; the full sender x radius matrix is a by-product.
#
# Reuses run_targeted_null_for_donor() from the main pipeline unchanged.
# Run inside RStudio (needs Seurat + dbscan).  Output: pipeline_v3/10_matched_radii/
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dbscan); library(qs)
  library(dplyr); library(tibble); library(readr); library(tidyr); library(stringr)
})

CODE_DIR <- "/broad/macosko/akraft/CAP/rerun/processed_pucks/aggregated/NoReg_Final_by_donor/Final_Cell_Arrangements_Oct9/Fig6"
PIPE_OUT <- file.path(CODE_DIR, "LR_sensitivity_analysis", "pipeline_v3")
OUT_DIR  <- file.path(PIPE_OUT, "16_radius600")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

MR <- list(
  radii      = c(600),
  trim       = 0.1,
  n_shuffles = 200,
  seed       = 42
)

logf <- file.path(OUT_DIR, "radius600.log")
say <- function(...) { m <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
                       cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }

say("=== SINGLE RADIUS 600um (completes Vascular K=5) ===")
say("radii: ", paste(MR$radii, collapse = ", "), " um | trim ", MR$trim,
    " | ", MR$n_shuffles, " shuffles | ALL senders")

# ---- the same fixed candidate set the original null used ---------------------
tgt <- readr::read_csv(file.path(PIPE_OUT, "03_null", "target_interactions_for_null.csv"),
                       show_col_types = FALSE)
say("candidate interactions: ", nrow(tgt))

grid <- tibble(Method = "radius", Scale = as.numeric(MR$radii), Trim = MR$trim)

# ---- database (cached from the main run) ------------------------------------
quiet <- list(log_msg = function(..., level = "INFO") invisible(NULL))
cellchat_db <- load_cellchatdb_human(PIPE_OUT, quiet)
db <- prepare_interaction_db(cellchat_db, c("Secreted Signaling", "Cell-Cell Contact"),
                             OUT_DIR, quiet)

donors <- sort(setdiff(names(donor_objs), "s20"))
say("donors: ", length(donors))

for (di in seq_along(donors)) {
  d_name <- donors[di]
  fout <- file.path(OUT_DIR, sprintf("%s_r600.csv", d_name))
  if (file.exists(fout)) { say("skip ", d_name, " (done)"); next }
  say("donor ", d_name, " (", di, "/", length(donors), ")")

  obj <- if (is.character(donor_objs[[d_name]])) qs::qread(donor_objs[[d_name]]) else donor_objs[[d_name]]
  d_data <- prep_donor_data(obj, db$required_genes,
                            c("MSN","Astrocyte","Oligo","OPC","Microglia","Vascular"),
                            1:6, quiet, d_name)
  rm(obj); invisible(gc(verbose = FALSE))
  if (is.null(d_data)) { say("  skipped (prep returned NULL)"); next }

  # same per-donor detection filter the pipeline applies (audit B4)
  donor_db <- filter_interactions_per_donor(d_data, db$interaction, db$entity_to_genes,
                                            0.05, quiet, d_name)
  tgt_d <- dplyr::semi_join(tgt, donor_db, by = c("Ligand" = "ligand", "Receptor" = "receptor"))
  if (nrow(tgt_d) == 0) { say("  no candidates pass detection"); next }

  res <- run_targeted_null_for_donor(
    d_name = d_name, d_data = d_data, target_df = tgt_d, complex_db = db$complex,
    scale_grid = grid, radii_um = MR$radii, knn_k = NULL, trim = NULL,
    n_shuffles = MR$n_shuffles, shuffle_seed = MR$seed, donor_index = di,
    zone_a = 1, zone_b = 4, min_senders_per_zone = 5, logger = quiet)

  readr::write_csv(res, fout)
  say("  wrote ", basename(fout), " rows=", nrow(res))
  rm(d_data, res); invisible(gc(verbose = FALSE))
}

# ---- combine and summarise --------------------------------------------------
fs <- list.files(OUT_DIR, pattern = "_matched\\.csv$", full.names = TRUE)
all <- bind_rows(lapply(fs, readr::read_csv, show_col_types = FALSE))
readr::write_csv(all, file.path(OUT_DIR, "r600_all_donors.csv"))

summ <- all %>%
  group_by(Method, Scale, Sender, Receiver, Interaction) %>%
  summarise(Donors = n(),
            Robust_95 = mean(abs(Real_Delta) > Null_95_Abs, na.rm = TRUE),
            .groups = "drop")
readr::write_csv(summ, file.path(OUT_DIR, "r600_robustness.csv"))

cnt <- summ %>% filter(Robust_95 >= 0.5) %>% count(Sender, Scale)
say("ROBUST COUNTS (Robust_95 >= 0.5), sender x radius:")
for (i in seq_len(nrow(cnt))) say("   ", cnt$Sender[i], "  ", cnt$Scale[i], " um : ", cnt$n[i])
say("DONE")
