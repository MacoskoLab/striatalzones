#!/usr/bin/env Rscript
# ==============================================================================
# run_explode_compare.R — launcher for the explode-vs-full-complex comparison
# at radius = 200 um, with the custom DB (TGFB2_TGFBR1_TGFBR3 injected;
# SHH_PTCH1/PTCH2 confirmed stock) and Rob65/Shuff99 reporting downstream.
#
# Sources LR_pipeline_radius_only.R (functions only), overrides CONFIG to a
# single radius + a fresh output dir, wires the 18 filtered donor .qs as lazy
# inputs, and runs main(). Complex handling is chosen by the LR_EXPLODE env var
# (read inside main()): unset/false -> full complex; true -> exploded.
#
# Usage (headless, one arm per invocation):
#   LR_OUT=<dir> LR_EXPLODE=false Rscript run_explode_compare.R   # full complex
#   LR_OUT=<dir> LR_EXPLODE=true  Rscript run_explode_compare.R   # exploded
# ==============================================================================
suppressWarnings(suppressMessages({
  library(Seurat); library(qs); library(data.table)
  library(dplyr); library(tidyr); library(readr); library(tibble); library(Matrix)
}))

FIG6 <- "/broad/macosko/akraft/CAP/rerun/processed_pucks/aggregated/NoReg_Final_by_donor/Final_Cell_Arrangements_Oct9/Fig6"
LRDIR <- file.path(FIG6, "LR_sensitivity_analysis")
PIPE  <- file.path(LRDIR, "LR_pipeline_radius_only.R")

options(LR_PIPELINE_NO_RUN = TRUE)   # load functions without auto-running main()
source(PIPE)

# ---- overrides -------------------------------------------------------------
out <- Sys.getenv("LR_OUT")
if (!nzchar(out)) stop("LR_OUT env var must be set")
rad <- as.numeric(Sys.getenv("LR_RADIUS", "200"))   # single scoring radius (um)
CONFIG$output_dir     <- out
CONFIG$radii_um       <- rad
CONFIG$primary_radius <- rad
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# reuse the already-downloaded CellChatDB.human.rda (no outbound internet needed)
src_rda <- file.path(LRDIR, "pipeline_v3", "00_db", "CellChatDB.human.rda")
dst_db  <- file.path(out, "00_db")
dir.create(dst_db, recursive = TRUE, showWarnings = FALSE)
if (file.exists(src_rda) && !file.exists(file.path(dst_db, "CellChatDB.human.rda")))
  file.copy(src_rda, file.path(dst_db, "CellChatDB.human.rda"))

# ---- donors: named vector of .qs paths (s20 auto-excluded inside main) ------
fdir <- file.path(LRDIR, "filtered_donor_objects")
qs_files <- list.files(fdir, pattern = "_filtered\\.qs$", full.names = TRUE)
nm <- sub("_filtered\\.qs$", "", basename(qs_files))
donor_objs <- setNames(qs_files, nm)

cat(sprintf("Launcher: LR_OUT=%s | LR_EXPLODE=%s | donors=%d\n",
            out, Sys.getenv("LR_EXPLODE", "false"), length(donor_objs)))

main()
cat("=== LAUNCHER DONE ===\n")
