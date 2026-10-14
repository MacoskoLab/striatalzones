#!/usr/bin/env Rscript

# ==============================================================================
# build_filtered_donor_objs.R
#
# Loader + one-time filter/cache for the LR sensitivity analysis.
#
# Source data: per-donor Seurat objects (all cell types) saved as .qs in
#   /broad/macosko/akraft/CAP/rerun/processed_pucks/aggregated/
#     NoReg_Final_by_donor/Final_Cell_Arrangements_Oct9/FullQC/donor_objects/
#   named  s<NN>_with_master.qs  (s5, s12, s15, ... s38; includes s20).
#
# Filtering (to shrink objects and speed up the LR pipeline):
#   1) Remove non-spatial cells           -> drop cells with NA x_um2 / y_um2.
#   2) Keep only zones 1-6                 -> Final_Zone_Assignments_No_Smooth %in% 1:6.
#      (White matter lies OUTSIDE zones 1-6, so this step also removes all WM.)
#   ALL cell types are retained (the LR pipeline does its own type assignment).
#
# Caching: each filtered object is written to
#   .../Fig6/LR_sensitivity_analysis/filtered_donor_objects/s<NN>_filtered.qs
# so subsequent runs reload the small filtered objects instead of re-filtering.
#
# Usage (interactive or via Rscript):
#   source("build_filtered_donor_objs.R")
#   donor_objs <- load_donor_objs()            # builds cache on first call, reloads after
#   # ... then run the LR pipeline, which expects `donor_objs` in the session.
#
# To force a rebuild (e.g. after changing filter rules):
#   donor_objs <- load_donor_objs(rebuild = TRUE)
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
})

# ---- Paths & filter config ----
RAW_DIR  <- "/broad/macosko/akraft/CAP/rerun/processed_pucks/aggregated/NoReg_Final_by_donor/Final_Cell_Arrangements_Oct9/FullQC/donor_objects"
OUT_BASE <- "/broad/macosko/akraft/CAP/rerun/processed_pucks/aggregated/NoReg_Final_by_donor/Final_Cell_Arrangements_Oct9/Fig6/LR_sensitivity_analysis"
FILT_DIR <- file.path(OUT_BASE, "filtered_donor_objects")

# Domain/zone label column. The generalized public name "Spatial_Domain" is
# preferred when present; the internal as-run name is the fallback.
ZONE_COL_CHOICES <- c("Spatial_Domain", "Final_Zone_Assignments_No_Smooth")
COORD_X     <- "x_um2"
COORD_Y     <- "y_um2"
VALID_ZONES <- 1:6                                 # WM is outside 1..6 -> dropped here

# Number of donor objects we expect to find (s5,s12,s15,s16,s17,s18,s20,s21,
# s23,s26,s29,s31,s32,s33,s34,s35,s36,s37,s38 = 19; s20 is dropped later by the
# LR pipeline, leaving 18).
EXPECTED_N_DONORS <- 19

# ------------------------------------------------------------------------------
# Filter a single Seurat object to spatial cells in zones 1-6.
# Returns the subset object; errors loudly if required columns are missing.
# ------------------------------------------------------------------------------
filter_donor_obj <- function(obj, donor_id) {
  meta <- obj@meta.data
  ZONE_COL <- intersect(ZONE_COL_CHOICES, colnames(meta))[1]
  if (is.na(ZONE_COL)) ZONE_COL <- ZONE_COL_CHOICES[2]

  req_cols <- c(ZONE_COL, COORD_X, COORD_Y)
  missing  <- setdiff(req_cols, colnames(meta))
  if (length(missing) > 0) {
    zone_like <- grep("zone", colnames(meta), ignore.case = TRUE, value = TRUE)
    stop("Donor ", donor_id, ": missing metadata column(s): ",
         paste(missing, collapse = ", "),
         ". Zone-like columns present: ",
         if (length(zone_like)) paste(zone_like, collapse = ", ") else "(none)")
  }

  # Coerce zone to integer regardless of factor/character storage.
  zone_int <- suppressWarnings(as.integer(as.character(meta[[ZONE_COL]])))

  keep_mask <- !is.na(meta[[COORD_X]]) &
               !is.na(meta[[COORD_Y]]) &
               !is.na(zone_int) &
               zone_int %in% VALID_ZONES

  keep_cells <- rownames(meta)[keep_mask]

  n_before <- ncol(obj)
  n_after  <- length(keep_cells)
  if (n_after == 0) stop("Donor ", donor_id, ": 0 cells survive filtering.")

  obj <- subset(obj, cells = keep_cells)

  attr(obj, ".filter_stats") <- c(before = n_before, after = n_after)
  obj
}

# ------------------------------------------------------------------------------
# Build (if needed) and load the filtered donor list.
#   rebuild = TRUE  -> re-filter from raw even if a cached file exists.
# Returns a NAMED list `donor_objs` keyed by donor id (s5, s12, ...).
# Also writes a filter summary CSV to OUT_BASE.
# ------------------------------------------------------------------------------
# lazy = TRUE (default): build/refresh the per-donor filtered cache WITHOUT
#   retaining the objects, and return a named CHARACTER VECTOR of .qs paths.
#   The pipeline loads each donor on demand and frees it, so peak memory is one
#   donor (~2-25 GB) instead of all 19 at once (~145 GB measured).
# lazy = FALSE: legacy behaviour — returns the actual Seurat objects.
load_donor_objs <- function(rebuild = FALSE, lazy = TRUE) {

  if (!dir.exists(RAW_DIR)) stop("RAW_DIR not found: ", RAW_DIR)
  dir.create(FILT_DIR, recursive = TRUE, showWarnings = FALSE)

  # Match exactly s<digits>_with_master.qs  (excludes s5_with_master_with_disease_scores.qs,
  # subdirectories, and other .qs in the folder).
  raw_files <- list.files(RAW_DIR, pattern = "^s[0-9]+_with_master\\.qs$", full.names = TRUE)
  if (length(raw_files) == 0) stop("No s*_with_master.qs files found in ", RAW_DIR)

  donor_ids <- sub("_with_master\\.qs$", "", basename(raw_files))

  message(sprintf("[%s] Found %d donor object files in %s",
                  format(Sys.time(), "%H:%M:%S"), length(raw_files), RAW_DIR))
  if (length(raw_files) != EXPECTED_N_DONORS) {
    warning(sprintf("Expected %d donor files, found %d (%s).",
                    EXPECTED_N_DONORS, length(raw_files),
                    paste(donor_ids, collapse = ",")))
  }

  donor_objs <- vector("list", length(raw_files))
  names(donor_objs) <- donor_ids
  stat_rows <- list()

  for (i in seq_along(raw_files)) {
    id        <- donor_ids[i]
    filt_path <- file.path(FILT_DIR, paste0(id, "_filtered.qs"))

    # Try the cache first; a truncated/corrupt .qs (e.g. job killed mid-write)
    # returns NULL here and falls through to a clean rebuild instead of erroring.
    cached_obj <- NULL
    if (lazy && !rebuild && file.exists(filt_path)) {
      # trust the cache; do not read it just to count cells
      message(sprintf("[%s] %s: cached (lazy, not loaded)", format(Sys.time(), "%H:%M:%S"), id))
      donor_objs[[id]] <- filt_path
      stat_rows[[id]] <- data.frame(Donor = id, Cells_Before = NA_integer_,
                                    Cells_After = NA_integer_, Source = "cache")
      next
    }
    if (!rebuild && file.exists(filt_path)) {
      cached_obj <- tryCatch(qs::qread(filt_path), error = function(e) {
        message(sprintf("[%s] %s: cached file unreadable (%s) -> rebuilding",
                        format(Sys.time(), "%H:%M:%S"), id, conditionMessage(e)))
        NULL
      })
    }

    if (!is.null(cached_obj)) {
      message(sprintf("[%s] %s: cached filtered object OK", format(Sys.time(), "%H:%M:%S"), id))
      obj <- cached_obj
      stat_rows[[id]] <- data.frame(Donor = id, Cells_Before = NA_integer_,
                                    Cells_After = ncol(obj), Source = "cache")
    } else {
      message(sprintf("[%s] %s: reading raw and filtering ...", format(Sys.time(), "%H:%M:%S"), id))
      obj <- qs::qread(raw_files[i])
      obj <- filter_donor_obj(obj, id)
      st  <- attr(obj, ".filter_stats"); attr(obj, ".filter_stats") <- NULL
      qs::qsave(obj, filt_path)
      message(sprintf("[%s] %s: %d -> %d cells (saved %s)",
                      format(Sys.time(), "%H:%M:%S"), id, st[["before"]], st[["after"]],
                      basename(filt_path)))
      stat_rows[[id]] <- data.frame(Donor = id, Cells_Before = st[["before"]],
                                    Cells_After = st[["after"]], Source = "filtered")
    }

    donor_objs[[id]] <- if (lazy) filt_path else obj
    rm(obj)
    invisible(gc(verbose = FALSE))
  }
  if (lazy) donor_objs <- unlist(donor_objs)

  # Filter summary
  summary_df <- do.call(rbind, stat_rows)
  rownames(summary_df) <- NULL
  summary_path <- file.path(OUT_BASE, "filtered_donor_summary.csv")
  utils::write.csv(summary_df, summary_path, row.names = FALSE)
  message(sprintf("[%s] Wrote filter summary: %s", format(Sys.time(), "%H:%M:%S"), summary_path))

  message(sprintf("[%s] Prepared %d donors%s: %s",
                  format(Sys.time(), "%H:%M:%S"), length(donor_objs),
                  if (lazy) " (lazy: paths only)" else "",
                  paste(names(donor_objs), collapse = ", ")))
  donor_objs
}

# If run non-interactively as a standalone build step (Rscript build_filtered_donor_objs.R),
# build the cache and report. When sourced interactively, just defines the functions.
if (!interactive() && sys.nframe() == 0) {
  donor_objs <- load_donor_objs()
  message("donor_objs is ready (", length(donor_objs), " donors).")
}
