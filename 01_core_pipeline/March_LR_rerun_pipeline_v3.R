#!/usr/bin/env Rscript

# ==============================================================================
# March_LR_rerun_pipeline_v3.R   (v3.1 — post-audit)
#
# Spatial ligand-receptor analysis, Zone 1 vs Zone 4, with a RADIUS SENSITIVITY
# SWEEP and a DENSITY-CONTROLLED kNN ARM, to address reviewer concerns:
#
#   (1) "Is the MSN-glia pattern stable across windows?" -> score every
#       interaction at radii {50,100,200,400,1000} um and report Delta/p/FDR
#       per window.
#   (2) "A single radius is not biologically equivalent across cell types."
#       -> kNN arm holds neighbor COUNT constant across cell types, and records
#       the realized neighbor DISTANCE so the spatial extent of each k is
#       reportable. Radius arm records neighbor COUNTS. The two arms bracket the
#       confound: radius holds distance fixed / count varies; kNN holds count
#       fixed / distance varies.
#
# ------------------------------------------------------------------------------
# CHANGES FROM THE FIRST v3 DRAFT (all from the code audit; see notes inline):
#
#  [B1] trim = 0 (was 0.10).  CRITICAL, and it also applies to the v2 results.
#       R's mean(x, trim=0.10) does NOTHING for n<=9 neighbors and then trims
#       from n=10 up, so the estimator's DEFINITION changes with neighbor count.
#       Verified: with identical receptor distributions in both zones and only
#       the neighbor count differing (6 vs 12), trimmed Delta = +0.267 while
#       plain Delta = +0.003 — a purely density-driven false signal, which the
#       label-permutation null CANNOT detect (it fixes per-sender scores).
#       Since the reviewer's concern is precisely density asymmetry, this had to
#       go. trim=0 is also what makes the potential a single sparse matmul.
#  [B2] BH family size is now FIXED across windows (was: per-window m, which
#       silently shrank when t.test returned NaN, making FDRs incomparable
#       across the sweep — i.e. corrupting the actual deliverable). t.test is
#       wrapped in tryCatch so a constant-delta interaction returns NA instead
#       of aborting the whole summarise after hours of scoring.
#  [B3] Per-donor checkpoints are keyed on the scoring CONFIG (not just the
#       scale grid) and written ATOMICALLY; the combine step globs by donor
#       name. Without this, changing `trim` and re-running would silently reuse
#       the old scores, and a job killed mid-write would be read as complete.
#  [B4] The null now applies the same per-donor detection filter as scoring, so
#       Donors_Tested matches N_Donors and Robust_95 has the right denominator.
#  [B5] Permutations are drawn ONCE per (donor, sender type) and reused across
#       every grid cell -> cross-scale comparison is paired under common random
#       numbers (better statistics AND faster), and adding the kNN arm in Stage
#       2 no longer perturbs the Stage-1 radius nulls.
#  [B6] CellChatDB is cached: a failed download no longer destroys the local
#       copy and no longer hard-stops a resumed job.
#  [B7] knn_k is read from the LR_KNN_K env var (re-sourcing the file used to
#       reset CONFIG$knn_k to NULL, making the documented Stage-2 instruction a
#       no-op). kNN feasibility is checked PER SCALE.
#  [P1/P2/P3] Vectorized potentials (sparse matmul), senders restricted to
#       zones A u B, columnar emit. ~26h -> well under 1h for scoring.
#  [G3] Frac0_ZA/ZB recorded: at small radii a sender with no neighbors scores
#       exactly 0 and is averaged in, so Mean_ZA partly measures "is a partner
#       nearby" rather than signaling. This cannot be recovered post hoc.
#  [G5] Empirical p uses the (1+r)/(1+n) estimator; the across-donor summary
#       reports a Fisher-combined p instead of an invalid median-of-p.
#  [G6] kNN arm carries kn$dist and reports realized k-th-neighbor distance.
#
# Assumptions: `donor_objs` is a named list of Seurat objects in the session
# (see build_filtered_donor_objs.R). s20 is excluded; 18 donors must remain.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(Matrix)
  library(dbscan)
})

# Seurat is needed ONLY to read donor objects (get_assay_matrix / prep_donor_data).
# Loading it lazily lets the scoring / stats / null functions run — and be
# unit-tested — in a minimal R library that has no Seurat.
# A donor entry may be a Seurat object OR a path to a filtered .qs (lazy mode,
# the default). Resolving on demand keeps peak memory at ONE donor instead of
# all 19 (~145 GB measured) — see build_filtered_donor_objs.R.
.resolve_donor <- function(x) {
  if (is.character(x) && length(x) == 1) {
    if (!requireNamespace("qs", quietly = TRUE))
      stop("Donor entries are file paths but the 'qs' package is unavailable.")
    return(qs::qread(x))
  }
  x
}

.need_seurat <- function() {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat is required to read donor Seurat objects, but is not installed ",
         "in this R library (.libPaths(): ", paste(.libPaths(), collapse = ", "), ")")
  }
}

# ==============================================================================
# computeExpr_LR — single genes get their expression; complexes get the
# geometric mean of subunit expression. (Unchanged; replaces CellChat's version,
# which has a dimension-drop bug when only one single-gene entity is present.)
# ==============================================================================
computeExpr_LR <- function(geneLR, data.use, complex_input) {
  nLR <- length(geneLR)
  nCells <- ncol(data.use)

  result <- matrix(0, nrow = nLR, ncol = nCells)
  rownames(result) <- geneLR
  colnames(result) <- colnames(data.use)

  complex_names <- rownames(complex_input)
  subunit_cols <- grep("^subunit", colnames(complex_input), value = TRUE)

  for (i in seq_len(nLR)) {
    entity <- geneLR[i]

    if (entity %in% complex_names) {
      subs <- unlist(complex_input[entity, subunit_cols, drop = TRUE], use.names = FALSE)
      subs <- subs[!is.na(subs) & subs != ""]
      if (length(subs) == 0) next

      subs_present <- subs[subs %in% rownames(data.use)]
      if (length(subs_present) != length(subs)) next

      if (length(subs_present) == 1) {
        result[i, ] <- as.numeric(data.use[subs_present[1], ])
      } else {
        sub_mat <- as.matrix(data.use[subs_present, , drop = FALSE])
        n_subs <- nrow(sub_mat)
        result[i, ] <- apply(sub_mat, 2, prod)^(1 / n_subs)
      }

    } else if (entity %in% rownames(data.use)) {
      result[i, ] <- as.numeric(data.use[entity, ])
    }
  }

  result
}

# ==============================================================================
# CONFIG
# ==============================================================================
CONFIG <- list(
  output_dir = "/broad/macosko/akraft/CAP/rerun/processed_pucks/aggregated/NoReg_Final_by_donor/Final_Cell_Arrangements_Oct9/Fig6/LR_sensitivity_analysis/pipeline_v3",

  # ---- SCALE GRID ----
  radii_um        = c(50, 100, 200, 400, 1000),
  primary_radius  = 200,
  primary_method  = "radius",

  # [B7] Density-controlled kNN arm, read from the environment so that
  # re-sourcing this file does not reset it. Stage 1: leave unset.
  # Stage 2:  Sys.setenv(LR_KNN_K = "3,10,20,40")   (or export LR_KNN_K=... )
  knn_k = {
    v <- Sys.getenv("LR_KNN_K", "")
    if (nzchar(v)) as.numeric(strsplit(v, ",")[[1]]) else NULL
  },

  # [B1, revised] TRIM IS NOW A SWEPT DIMENSION, not a fixed choice.
  #
  # trim = 0.10 (as in v2 / prior publications) is NEUTRAL for the Zone1-vs-Zone4
  # contrast whenever the two zones have comparable receiver-neighbour counts —
  # verified: matched n, no biological difference -> Delta ~ 0 at every density.
  #
  # The narrow risk is the n = ceiling(1/trim) = 10 THRESHOLD, where R's
  # mean(x, trim=) switches from keeping 100% of values to keeping 80%. A pair
  # whose mean neighbour count is BELOW the threshold in one zone and ABOVE it
  # in the other picks up a large spurious Delta (measured +0.25 for 8 vs 12
  # neighbours with identical underlying distributions). Both-below or
  # both-well-above is fine. Frac0/N_nbr and the Straddle flag in
  # 01_neighbor_density_summary.csv identify affected pairs.
  #
  # For the RADIUS SWEEP specifically, a proportional trim is not comparable
  # across windows (neighbourhoods cross the threshold as radius grows), so we
  # score every scale at BOTH trim levels and report them side by side.
  trim_levels  = c(0.10, 0),
  primary_trim = 0.10,

  zone_a = 1,
  zone_b = 4,
  valid_zones = c(1, 2, 3, 4, 5, 6),

  cell_types = c("MSN", "Astrocyte", "OPC", "Oligo", "Microglia", "Vascular"),

  min_pct_expr = 0.05,
  valid_annotations = c("Secreted Signaling", "Cell-Cell Contact"),

  # Minimum senders per zone for a donor-level mean to be reported
  min_senders_per_zone = 5,

  # Permutation null
  n_shuffles = 200,
  shuffle_seed = 42,
  fdr_for_null = 0.05,

  log_file = "pipeline_progress.log",
  status_file = "pipeline_status.txt"
)

# Grid is Method x Scale x Trim. Trim is swept so the sensitivity analysis is
# not confounded by the estimator changing definition across windows.
build_scale_grid <- function(cfg) {
  base <- tibble(Method = "radius", Scale = as.numeric(cfg$radii_um))
  if (!is.null(cfg$knn_k) && length(cfg$knn_k) > 0) {
    base <- bind_rows(base, tibble(Method = "knn", Scale = as.numeric(cfg$knn_k)))
  }
  trims <- unique(as.numeric(cfg$trim_levels))
  # Radius arm sweeps every trim level (that is where the trim question lives).
  # kNN arm runs at the PRIMARY trim only: for k < 1/trim the trimmed and plain
  # means are identical anyway, and the trim comparison is already settled by the
  # radius arm. Halves the Stage-2 cost.
  g <- tidyr::crossing(base %>% filter(Method == "radius"), Trim = trims)
  if (any(base$Method == "knn")) {
    g <- bind_rows(g, base %>% filter(Method == "knn") %>%
                     mutate(Trim = as.numeric(cfg$primary_trim)))
  }
  g %>% arrange(match(Method, c("radius", "knn")), Scale, desc(Trim))
}

# n at which mean(x, trim=t) starts discarding values: floor(n*t) >= 1.
trim_threshold_n <- function(trim) if (trim <= 0) Inf else ceiling(1 / trim)

# ==============================================================================
# [B3] Config-keyed, atomic checkpointing
# ==============================================================================

# Dependency-free deterministic short hash (no `digest` requirement).
cfg_hash <- function(parts) {
  s <- paste(vapply(parts, function(v) paste(as.character(v), collapse = "|"),
                    character(1)), collapse = "||")
  iv <- utf8ToInt(s)
  h <- 0
  for (i in iv) h <- (h * 31 + i) %% 2147483647
  sprintf("%07x", h)
}

# Key on everything that changes the SCORES. Trim is NOT included: it is part of
# the grid, so adding a trim level is handled by the grid-coverage check.
make_cfg_key <- function(cfg) {
  paste0("z", cfg$zone_a, "v", cfg$zone_b, "_",
         cfg_hash(list(cfg$zone_a, cfg$zone_b, cfg$valid_zones,
                       cfg$cell_types, cfg$min_pct_expr, cfg$valid_annotations,
                       cfg$min_senders_per_zone)))
}

write_csv_atomic <- function(x, path) {
  tmp <- paste0(path, ".tmp", Sys.getpid())
  readr::write_csv(x, tmp)
  if (!file.rename(tmp, path)) {
    unlink(tmp)
    stop("atomic rename failed for ", path)
  }
  invisible(path)
}

# A checkpoint counts as complete only if it contains every (Method, Scale) in
# the current grid. Unreadable/partial files -> recompute.
checkpoint_covers_grid <- function(path, scale_grid) {
  if (!file.exists(path)) return(FALSE)
  existing <- tryCatch(
    suppressMessages(readr::read_csv(path, show_col_types = FALSE,
                                     col_select = c("Method", "Scale", "Trim"))),
    error = function(e) NULL
  )
  if (is.null(existing) || nrow(existing) == 0) return(FALSE)
  have <- existing %>%
    mutate(Scale = as.numeric(Scale), Trim = as.numeric(Trim)) %>%
    distinct(Method, Scale, Trim)
  want <- scale_grid %>%
    mutate(Scale = as.numeric(Scale), Trim = as.numeric(Trim)) %>%
    distinct(Method, Scale, Trim)
  nrow(dplyr::anti_join(want, have, by = c("Method", "Scale", "Trim"))) == 0
}

# ==============================================================================
# Logging
# ==============================================================================
make_logger <- function(output_dir, log_file, status_file) {
  log_path <- file.path(output_dir, log_file)
  status_path <- file.path(output_dir, status_file)

  if (!file.exists(log_path)) writeLines(character(), log_path)
  if (!file.exists(status_path)) writeLines(character(), status_path)

  log_msg <- function(..., level = "INFO") {
    msg <- paste0(..., collapse = "")
    line <- sprintf("[%s] [%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, msg)
    cat(line, "\n")
    cat(line, "\n", file = log_path, append = TRUE)
    writeLines(line, status_path)
    invisible(line)
  }

  list(log_path = log_path, status_path = status_path, log_msg = log_msg)
}

get_assay_matrix <- function(obj, assay = "RNA", layer = "data") {
  .need_seurat()
  if ("layer" %in% names(formals(Seurat::GetAssayData))) {
    return(Seurat::GetAssayData(obj, assay = assay, layer = layer))
  }
  Seurat::GetAssayData(obj, assay = assay, slot = layer)
}

# ==============================================================================
# [B6] CellChatDB — cache if present; a failed fetch never destroys the cache
# ==============================================================================
load_cellchatdb_human <- function(output_dir, logger) {
  url <- "https://github.com/jinworks/CellChat/raw/main/data/CellChatDB.human.rda"
  db_dir <- file.path(output_dir, "00_db")
  if (!dir.exists(db_dir)) dir.create(db_dir, recursive = TRUE)
  rda_path <- file.path(db_dir, "CellChatDB.human.rda")

  if (file.exists(rda_path)) {
    logger$log_msg("Using CACHED CellChatDB: ", rda_path,
                   " (md5=", as.character(tools::md5sum(rda_path)), ")")
  } else {
    logger$log_msg("Downloading CellChatDB.human: ", url)
    tmp <- paste0(rda_path, ".part")
    ok <- tryCatch(
      utils::download.file(url, destfile = tmp, mode = "wb", quiet = TRUE) == 0,
      error = function(e) FALSE
    )
    if (!ok || !file.exists(tmp)) {
      unlink(tmp)
      stop("Failed to download CellChatDB.human.rda. If this node has no outbound ",
           "internet, copy the .rda to: ", rda_path)
    }
    file.rename(tmp, rda_path)
    logger$log_msg("Downloaded DB (md5=", as.character(tools::md5sum(rda_path)), ")")
  }

  e <- new.env(parent = emptyenv())
  load(rda_path, envir = e)
  if (!exists("CellChatDB.human", envir = e)) stop("CellChatDB.human not found in .rda")
  db <- get("CellChatDB.human", envir = e)

  meta_out <- tibble(
    db_source = rda_path,
    used_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    md5 = as.character(tools::md5sum(rda_path))
  )
  # append, so provenance across restarts is preserved
  readr::write_csv(meta_out, file.path(db_dir, "CellChatDB_human_download_metadata.csv"),
                   append = file.exists(file.path(db_dir, "CellChatDB_human_download_metadata.csv")))
  readr::write_csv(as_tibble(db$interaction), file.path(db_dir, "CellChatDB_human_interaction_full.csv"))
  readr::write_csv(tibble::rownames_to_column(as.data.frame(db$complex), "complex_name"),
                   file.path(db_dir, "CellChatDB_human_complex_full.csv"))
  db
}

# ==============================================================================
# Filtered interaction DB + gene requirements
# ==============================================================================
prepare_interaction_db <- function(cellchat_db, valid_annotations, output_dir, logger) {
  interaction <- cellchat_db$interaction %>%
    dplyr::filter(.data$annotation %in% valid_annotations) %>%
    dplyr::select(interaction_name, ligand, receptor, annotation, pathway_name) %>%
    dplyr::mutate(
      pathway_name = ifelse(is.na(.data$pathway_name) | .data$pathway_name == "",
                            "Unknown", .data$pathway_name)
    )

  complex <- cellchat_db$complex
  entities <- unique(c(interaction$ligand, interaction$receptor))
  complex_names <- rownames(complex)
  subunit_cols <- grep("^subunit", colnames(complex), value = TRUE)

  get_subunits <- function(entity) {
    if (!(entity %in% complex_names)) return(entity)
    v <- unlist(complex[entity, subunit_cols, drop = TRUE], use.names = FALSE)
    v <- v[!is.na(v) & v != ""]
    unique(v)
  }

  entity_to_genes <- lapply(entities, get_subunits)
  names(entity_to_genes) <- entities
  required_genes <- unique(unlist(entity_to_genes, use.names = FALSE))

  out_dir <- file.path(output_dir, "00_db")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  readr::write_csv(interaction, file.path(out_dir, "CellChatDB_human_interactions_filtered.csv"))
  readr::write_csv(interaction %>% distinct(pathway_name) %>% arrange(pathway_name) %>%
                     mutate(pathway_id = row_number()),
                   file.path(out_dir, "pathway_name_to_id.csv"))
  readr::write_csv(interaction %>% distinct(interaction_name, ligand, receptor,
                                            annotation, pathway_name),
                   file.path(out_dir, "interaction_to_pathway_mapping.csv"))
  readr::write_csv(tibble(
    entity = names(entity_to_genes),
    is_complex = names(entity_to_genes) %in% complex_names,
    subunits = vapply(entity_to_genes, function(v) paste(v, collapse = ";"), character(1))
  ), file.path(out_dir, "entity_to_subunits_mapping.csv"))

  logger$log_msg("Filtered CellChat DB interactions: ", nrow(interaction),
                 " | Entities: ", length(entities),
                 " | Required genes: ", length(required_genes))

  list(interaction = interaction, complex = complex, entities = entities,
       entity_to_genes = entity_to_genes, required_genes = required_genes)
}

# ==============================================================================
# Public-release input adapter
# ==============================================================================
# Generalized metadata columns take precedence over the internal (as-run) names:
#   Spatial_Domain -> used as the domain/zone label (integer-coercible; the
#                     as-run internal name Final_Zone_Assignments_No_Smooth is
#                     still accepted for exact reproduction of the paper runs)
#   Cell_types     -> used directly as the analysis cell type (bypasses the
#                     internal Class_name / AK_celltype mapping below)
normalize_input_meta <- function(meta) {
  if ("Spatial_Domain" %in% colnames(meta))
    meta$Final_Zone_Assignments_No_Smooth <- meta$Spatial_Domain
  meta
}

# ==============================================================================
# Input audit
# ==============================================================================
audit_inputs <- function(donor_objects_list, output_dir, logger) {
  out_dir <- file.path(output_dir, "00_qc")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  rows <- lapply(names(donor_objects_list), function(d_name) {
    o <- .resolve_donor(donor_objects_list[[d_name]])
    meta <- normalize_input_meta(o@meta.data)
    on.exit({ rm(o); invisible(gc(verbose = FALSE)) }, add = TRUE)
    valid <- !is.na(meta$x_um2) & !is.na(meta$y_um2) &
      !is.na(meta$Final_Zone_Assignments_No_Smooth)
    tibble(
      Donor = d_name,
      Cells_Total = nrow(meta),
      Cells_WithCoordsAndZone = sum(valid),
      Cells_Zone1 = sum(valid & meta$Final_Zone_Assignments_No_Smooth == 1),
      Cells_Zone4 = sum(valid & meta$Final_Zone_Assignments_No_Smooth == 4)
    )
  })
  audit_df <- bind_rows(rows) %>% arrange(Donor)
  readr::write_csv(audit_df, file.path(out_dir, "donor_input_audit.csv"))
  logger$log_msg("Wrote input audit")
  audit_df
}

# ==============================================================================
# Donor prep
# ==============================================================================
assign_analysis_type <- function(meta) {
  # Generalized input: a Cell_types column is used directly as the analysis type.
  if ("Cell_types" %in% colnames(meta)) {
    meta$Analysis_Type <- as.character(meta$Cell_types)
    return(meta)
  }
  if (!"Class_name" %in% colnames(meta))
    stop("No cell-type metadata found: provide a 'Cell_types' column ",
         "(or the internal Class_name / AK_celltype columns).")
  meta$Analysis_Type <- NA_character_
  msn_subtypes <- c("STR D2 Matrix MSN", "STR D1 Matrix MSN",
                    "STRd D1 Striosome MSN", "STRd D2 Striosome MSN")
  is_msn <- meta$Class_name == "CN LGE GABA" & meta$AK_celltype %in% msn_subtypes
  meta$Analysis_Type[is_msn] <- "MSN"
  meta$Analysis_Type[meta$Class_name == "Astro-Epen"] <- "Astrocyte"
  meta$Analysis_Type[meta$Class_name == "Immune"] <- "Microglia"
  meta$Analysis_Type[meta$Class_name == "Vascular"] <- "Vascular"
  if ("AK_celltype" %in% colnames(meta)) {
    meta$Analysis_Type[meta$Class_name == "OPC-Oligo" & meta$AK_celltype == "OPC"] <- "OPC"
    meta$Analysis_Type[meta$Class_name == "OPC-Oligo" & meta$AK_celltype == "Oligo"] <- "Oligo"
  }
  meta
}

build_data_use <- function(obj, cells, required_genes, assay = "RNA", layer = "data") {
  feats <- rownames(obj)
  feat_map <- setNames(feats, toupper(feats))
  match_feat <- feat_map[toupper(required_genes)]
  present <- !is.na(match_feat)

  expr <- get_assay_matrix(obj, assay = assay, layer = layer)
  present_mat <- expr[match_feat[present], cells, drop = FALSE]
  rownames(present_mat) <- required_genes[present]

  missing_genes <- required_genes[!present]
  if (length(missing_genes) > 0) {
    zero_mat <- Matrix::Matrix(0, nrow = length(missing_genes),
                               ncol = length(cells), sparse = TRUE)
    rownames(zero_mat) <- missing_genes
    colnames(zero_mat) <- cells
    data.use <- rbind(present_mat, zero_mat)
  } else {
    data.use <- present_mat
  }
  data.use[required_genes, , drop = FALSE]
}

prep_donor_data <- function(obj, required_genes, cell_types, valid_zones, logger, donor_name) {
  meta <- normalize_input_meta(obj@meta.data)
  zone_int <- suppressWarnings(as.integer(as.character(meta$Final_Zone_Assignments_No_Smooth)))

  valid_cells <- rownames(meta)[
    !is.na(meta$x_um2) & !is.na(meta$y_um2) &
      !is.na(zone_int) & zone_int %in% valid_zones
  ]
  if (length(valid_cells) < 100) {
    logger$log_msg("Skipping donor ", donor_name, ": <100 valid cells")
    return(NULL)
  }

  meta_valid <- assign_analysis_type(meta[valid_cells, , drop = FALSE])
  target_cells <- rownames(meta_valid)[meta_valid$Analysis_Type %in% cell_types]
  if (length(target_cells) < 50) {
    logger$log_msg("Skipping donor ", donor_name, ": <50 target cells after type filtering")
    return(NULL)
  }

  meta_final <- meta_valid[target_cells, , drop = FALSE]
  meta_final$Final_Zone_Assignments_No_Smooth <-
    suppressWarnings(as.integer(as.character(meta_final$Final_Zone_Assignments_No_Smooth)))

  list(
    meta = meta_final,
    coords = as.matrix(meta_final[, c("x_um2", "y_um2"), drop = FALSE]),
    zones = meta_final$Final_Zone_Assignments_No_Smooth,
    type_indices = split(seq_len(nrow(meta_final)), meta_final$Analysis_Type),
    data.use = build_data_use(obj, target_cells, required_genes, assay = "RNA", layer = "data"),
    cell_ids = rownames(meta_final)
  )
}

# ==============================================================================
# Per-donor interaction filtering (detection across ALL zones x cell types)
# ==============================================================================
filter_interactions_per_donor <- function(d_data, interaction_db, entity_to_genes,
                                          min_pct_expr, logger, donor_name) {
  required_genes <- rownames(d_data$data.use)
  gene_max_pct <- setNames(rep(0, length(required_genes)), required_genes)

  det <- d_data$data.use > 0
  zones <- d_data$zones

  for (ctype in names(d_data$type_indices)) {
    idx <- d_data$type_indices[[ctype]]
    for (z in sort(unique(zones))) {
      idx_z <- idx[zones[idx] == z]
      if (length(idx_z) > 10) {
        gene_max_pct <- pmax(gene_max_pct, Matrix::rowMeans(det[, idx_z, drop = FALSE]))
      }
    }
  }

  valid_genes <- names(gene_max_pct)[gene_max_pct >= min_pct_expr]
  entities <- names(entity_to_genes)
  entity_valid <- vapply(entities, function(e) all(entity_to_genes[[e]] %in% valid_genes),
                         logical(1))
  names(entity_valid) <- entities

  # isTRUE-style guard: an entity absent from the map must not produce NA
  lig_ok <- entity_valid[interaction_db$ligand];  lig_ok[is.na(lig_ok)] <- FALSE
  rec_ok <- entity_valid[interaction_db$receptor]; rec_ok[is.na(rec_ok)] <- FALSE
  filtered <- interaction_db[lig_ok & rec_ok, , drop = FALSE]

  logger$log_msg("Donor ", donor_name, ": interactions passing detection filter: ",
                 nrow(filtered), " / ", nrow(interaction_db))
  filtered
}

# ==============================================================================
# Neighbor construction (shared by scoring and the null)
#
# A "neighbor set" for (sender type, receiver type, scale) is a list over sender
# cells, each element an integer vector of LOCAL receiver indices (positions in
# r_rows = d_data$type_indices[[receiver_type]]).
# ==============================================================================

make_cell_type_vec <- function(d_data) {
  v <- character(nrow(d_data$meta))
  for (ctype in names(d_data$type_indices)) v[d_data$type_indices[[ctype]]] <- ctype
  v
}

# ONE frNN per sender type at the LARGEST radius, self-excluded, keeping $dist.
frnn_for_sender <- function(d_data, s_rows, max_radius) {
  nn_res <- dbscan::frNN(d_data$coords, eps = max_radius,
                         query = d_data$coords[s_rows, , drop = FALSE])
  ids <- nn_res$id; dsts <- nn_res$dist
  for (i in seq_along(ids)) {
    keep <- ids[[i]] != s_rows[i]      # self-exclusion
    ids[[i]] <- ids[[i]][keep]
    dsts[[i]] <- dsts[[i]][keep]
  }
  list(id = ids, dist = dsts)
}

restrict_frnn_to_receiver <- function(frnn, cell_type_vec, rtype, global_to_local) {
  ids <- frnn$id; dsts <- frnn$dist
  n <- length(ids)
  local <- vector("list", n); dist_out <- vector("list", n)
  for (i in seq_len(n)) {
    g <- ids[[i]]
    if (length(g) == 0) { local[[i]] <- integer(0); dist_out[[i]] <- numeric(0); next }
    keep <- cell_type_vec[g] == rtype
    local[[i]] <- global_to_local[g[keep]]
    dist_out[[i]] <- dsts[[i]][keep]
  }
  list(local = local, dist = dist_out)
}

radius_subset <- function(restricted, radius) {
  mapply(function(loc, dd) loc[dd <= radius],
         restricted$local, restricted$dist, SIMPLIFY = FALSE)
}
radius_subset_dist <- function(restricted, radius) {
  mapply(function(dd) dd[dd <= radius], restricted$dist, SIMPLIFY = FALSE)
}

# [B7.2 / G6] kNN for one (sender type, receiver type). Returns ORDERED
# nearest-first local receiver indices AND their distances, self-excluded when
# sender_type == receiver_type. `k_feasible` reports the largest usable k so the
# grid loop can skip only the scales that are actually infeasible.
knn_for_pair <- function(d_data, s_rows, r_rows, sender_eq_receiver, max_k) {
  n_r <- length(r_rows)
  k_feasible <- min(max_k, n_r - as.integer(sender_eq_receiver))
  if (k_feasible < 1) return(NULL)

  need <- k_feasible + as.integer(sender_eq_receiver)
  need <- min(need, n_r)

  kn <- dbscan::kNN(x = d_data$coords[r_rows, , drop = FALSE],
                    query = d_data$coords[s_rows, , drop = FALSE],
                    k = need)
  id_mat <- matrix(kn$id, nrow = length(s_rows))
  d_mat  <- matrix(kn$dist, nrow = length(s_rows))

  if (!sender_eq_receiver) {
    ids <- lapply(seq_len(nrow(id_mat)), function(i) id_mat[i, seq_len(k_feasible)])
    dss <- lapply(seq_len(nrow(d_mat)),  function(i) d_mat[i, seq_len(k_feasible)])
  } else {
    self_local <- match(s_rows, r_rows)
    ids <- vector("list", nrow(id_mat)); dss <- vector("list", nrow(id_mat))
    for (i in seq_len(nrow(id_mat))) {
      keep <- id_mat[i, ] != self_local[i]
      ii <- id_mat[i, keep]; dd <- d_mat[i, keep]
      take <- seq_len(min(k_feasible, length(ii)))
      ids[[i]] <- ii[take]; dss[[i]] <- dd[take]
    }
  }
  list(id = ids, dist = dss, k_feasible = k_feasible)
}

# [P1] All receptors at once. Exactly equivalent to a plain (trim=0) mean of
# receptor expression over each sender's neighbor set.
#   nn_local : list over senders of local receiver indices
#   Rmat     : [n_receivers x n_receptors]
potentials_all <- function(nn_local, Rmat) {
  n <- length(nn_local)
  out <- matrix(0, nrow = n, ncol = ncol(Rmat),
                dimnames = list(NULL, colnames(Rmat)))
  if (n == 0 || ncol(Rmat) == 0) return(out)
  len <- lengths(nn_local)
  nz <- which(len > 0)
  if (length(nz) == 0) return(out)

  i <- rep.int(nz, len[nz])
  j <- unlist(nn_local[nz], use.names = FALSE)
  S <- Matrix::sparseMatrix(i = i, j = j, x = 1 / len[i],
                            dims = c(n, nrow(Rmat)))
  out[] <- as.matrix(S %*% Rmat)
  colnames(out) <- colnames(Rmat)
  out
}

# Reference (naive) implementation — kept for the correctness test only.
potential_from_neighbors <- function(nn_local, r_vec, trim = 0) {
  vapply(nn_local, function(idx) {
    if (length(idx) == 0) return(0)
    mean(r_vec[idx], trim = trim)
  }, numeric(1))
}

# Build the [n_receivers x n_receptors] receptor matrix (sparse: expression is
# zero-inflated, and this is the big transient).
build_Rmat <- function(entity_expr, receptors, r_rows) {
  m <- entity_expr[receptors, r_rows, drop = FALSE]
  R <- Matrix::t(Matrix::Matrix(m, sparse = TRUE))
  colnames(R) <- receptors
  R
}

# ==============================================================================
# Core scoring across the Method x Scale grid
# ==============================================================================
score_donor_zone_comparison <- function(
    d_name, d_data, interaction_db_donor, complex_db,
    scale_grid, radii_um, knn_k, trim, cell_types, zone_a, zone_b,
    min_senders_per_zone, logger) {

  if (nrow(interaction_db_donor) == 0) return(tibble())

  receptors <- unique(interaction_db_donor$receptor)
  all_entities <- unique(c(interaction_db_donor$ligand, receptors))
  entity_expr <- computeExpr_LR(all_entities, d_data$data.use, complex_db)
  rownames(entity_expr) <- all_entities
  active_receptors <- intersect(receptors, rownames(entity_expr))

  zones <- d_data$zones
  zA_cells <- which(zones == zone_a)
  zB_cells <- which(zones == zone_b)
  zoned_all <- c(zA_cells, zB_cells)

  n_cells <- nrow(d_data$meta)
  cell_type_vec <- make_cell_type_vec(d_data)

  r_g2l <- list()
  for (rtype in cell_types) {
    if (!(rtype %in% names(d_data$type_indices))) next
    rr <- d_data$type_indices[[rtype]]
    m <- integer(n_cells); m[rr] <- seq_along(rr)
    r_g2l[[rtype]] <- m
  }

  # Derive the work from the GRID actually passed in (supports additive re-runs
  # that compute only the missing cells).
  radius_scales <- sort(unique(scale_grid$Scale[scale_grid$Method == "radius"]))
  knn_scales    <- sort(unique(scale_grid$Scale[scale_grid$Method == "knn"]))
  need_radius <- length(radius_scales) > 0
  do_knn      <- length(knn_scales) > 0
  max_radius  <- if (need_radius) max(radius_scales) else NA_real_
  max_k       <- if (do_knn) max(knn_scales) else NA_integer_

  results <- list()

  for (sender_type in cell_types) {
    if (!(sender_type %in% names(d_data$type_indices))) next
    s_rows_all <- d_data$type_indices[[sender_type]]
    if (length(s_rows_all) < 10) next

    # [P2] Ligand-nonzero mask over ALL senders of this type, to preserve the
    # v2 row set, BEFORE restricting the query to zoned senders.
    lig_v <- unique(interaction_db_donor$ligand)
    lig_nz <- setNames(
      Matrix::rowSums(entity_expr[lig_v, s_rows_all, drop = FALSE]) > 0, lig_v)

    # [P2] Only zone-A / zone-B senders are ever scored -> ~3x less work.
    s_rows <- s_rows_all[s_rows_all %in% zoned_all]
    if (length(s_rows) < 2) next
    s_in_zA <- s_rows %in% zA_cells
    s_in_zB <- s_rows %in% zB_cells
    nA_s <- sum(s_in_zA); nB_s <- sum(s_in_zB)
    has_zA <- nA_s >= min_senders_per_zone
    has_zB <- nB_s >= min_senders_per_zone
    if (!has_zA && !has_zB) next
    A <- which(s_in_zA); B <- which(s_in_zB)

    frnn <- NULL
    if (need_radius) {
      logger$log_msg("Donor ", d_name, ": frNN(max=", max_radius, ") sender ", sender_type,
                     " (", length(s_rows), " zoned senders of ", length(s_rows_all), ")")
      frnn <- frnn_for_sender(d_data, s_rows, max_radius)
    }

    for (receiver_type in cell_types) {
      if (!(receiver_type %in% names(r_g2l))) next
      r_rows <- d_data$type_indices[[receiver_type]]
      if (length(r_rows) < 10) next

      logger$log_msg("Donor ", d_name, ": scoring ", sender_type, " -> ", receiver_type)

      restricted <- if (need_radius)
        restrict_frnn_to_receiver(frnn, cell_type_vec, receiver_type, r_g2l[[receiver_type]]) else NULL
      Rmat <- build_Rmat(entity_expr, active_receptors, r_rows)

      knn_pair <- NULL
      if (do_knn) {
        knn_pair <- knn_for_pair(d_data, s_rows, r_rows,
                                 sender_eq_receiver = (sender_type == receiver_type),
                                 max_k = max_k)
        if (is.null(knn_pair)) {
          logger$log_msg("Donor ", d_name, ": kNN infeasible for ", sender_type, " -> ",
                         receiver_type, " (", length(r_rows), " receivers)", level = "WARN")
        }
      }

      # Neighbour sets depend on (Method, Scale) only — build once, reuse for
      # every trim level.
      ms <- scale_grid %>% distinct(Method, Scale)

      for (mi in seq_len(nrow(ms))) {
        method <- ms$Method[mi]
        scale  <- ms$Scale[mi]

        d_kth <- NULL
        if (method == "radius") {
          nn_local <- radius_subset(restricted, scale)
        } else {
          if (is.null(knn_pair)) next
          if (scale > knn_pair$k_feasible) {
            logger$log_msg("Donor ", d_name, ": skip k=", scale, " for ", sender_type,
                           " -> ", receiver_type, " (k_feasible=", knn_pair$k_feasible, ")",
                           level = "WARN")
            next
          }
          nn_local <- lapply(knn_pair$id, function(v) v[seq_len(min(scale, length(v)))])
          # [G6] realized distance of the k-th (farthest kept) neighbor
          d_kth <- vapply(knn_pair$dist, function(v) {
            vv <- v[seq_len(min(scale, length(v)))]
            if (length(vv) == 0) NA_real_ else max(vv)
          }, numeric(1))
        }

        nbr_counts <- lengths(nn_local)
        n_nbr_ZA <- if (has_zA) mean(nbr_counts[A]) else NA_real_
        n_nbr_ZB <- if (has_zB) mean(nbr_counts[B]) else NA_real_
        # [G3] fraction of senders with NO partner — at small radii this, not
        # expression, drives the zone difference.
        frac0_ZA <- if (has_zA) mean(nbr_counts[A] == 0) else NA_real_
        frac0_ZB <- if (has_zB) mean(nbr_counts[B] == 0) else NA_real_
        # [G6] realized spatial extent of the kNN arm
        dk_ZA <- if (!is.null(d_kth) && has_zA) stats::median(d_kth[A], na.rm = TRUE) else NA_real_
        dk_ZB <- if (!is.null(d_kth) && has_zB) stats::median(d_kth[B], na.rm = TRUE) else NA_real_
        dk90_ZA <- if (!is.null(d_kth) && has_zA)
          as.numeric(stats::quantile(d_kth[A], 0.90, na.rm = TRUE, names = FALSE)) else NA_real_
        dk90_ZB <- if (!is.null(d_kth) && has_zB)
          as.numeric(stats::quantile(d_kth[B], 0.90, na.rm = TRUE, names = FALSE)) else NA_real_

        keep_k <- which(unname(lig_nz[interaction_db_donor$ligand]) &
                          interaction_db_donor$receptor %in% active_receptors)
        if (length(keep_k) == 0) next
        ivar <- interaction_db_donor[keep_k, , drop = FALSE]
        Lk <- entity_expr[ivar$ligand, s_rows, drop = FALSE]          # [K x nS]

        trims_here <- sort(unique(scale_grid$Trim[scale_grid$Method == method &
                                                   scale_grid$Scale == scale]),
                           decreasing = TRUE)

        for (trim_v in trims_here) {
          # [P1] trim = 0 collapses to one sparse matmul over all receptors.
          # trim > 0 needs the per-receptor path (a trimmed mean requires the
          # neighbourhood's order statistics, so it cannot be a matmul).
          potentials_mat <- if (trim_v == 0) {
            potentials_all(nn_local, Rmat)
          } else {
            pm <- matrix(0, nrow = length(nn_local), ncol = length(active_receptors),
                         dimnames = list(NULL, active_receptors))
            for (rc in unique(ivar$receptor)) {
              r_vec <- as.numeric(entity_expr[rc, r_rows])
              if (sum(r_vec) == 0) next
              pm[, rc] <- potential_from_neighbors(nn_local, r_vec, trim = trim_v)
            }
            pm
          }

          tPk <- t(potentials_mat[, ivar$receptor, drop = FALSE])     # [K x nS]

          mean_zA <- if (has_zA)
            rowSums(Lk[, A, drop = FALSE] * tPk[, A, drop = FALSE]) / length(A)
          else rep(NA_real_, nrow(ivar))
          mean_zB <- if (has_zB)
            rowSums(Lk[, B, drop = FALSE] * tPk[, B, drop = FALSE]) / length(B)
          else rep(NA_real_, nrow(ivar))

          # Does this pair straddle the trim threshold between zones? If so the
          # trimmed estimator is not comparable between the two zones.
          thr <- trim_threshold_n(trim_v)
          straddle <- is.finite(thr) && !is.na(n_nbr_ZA) && !is.na(n_nbr_ZB) &&
            ((n_nbr_ZA < thr) != (n_nbr_ZB < thr))

          results[[length(results) + 1L]] <- tibble(
            Donor = d_name, Method = method, Scale = scale, Trim = trim_v,
            Sender = sender_type, Receiver = receiver_type,
            Interaction = ivar$interaction_name,
            Ligand = ivar$ligand, Receptor = ivar$receptor,
            Pathway = ivar$pathway_name, Annotation = ivar$annotation,
            Mean_ZA = unname(mean_zA), Mean_ZB = unname(mean_zB),
            Delta_ZA_ZB = unname(mean_zA - mean_zB),
            N_senders_ZA = nA_s, N_senders_ZB = nB_s,
            N_nbr_ZA = n_nbr_ZA, N_nbr_ZB = n_nbr_ZB,
            Frac0_ZA = frac0_ZA, Frac0_ZB = frac0_ZB,
            Trim_Straddle = straddle,
            Dkth_med_ZA = dk_ZA, Dkth_med_ZB = dk_ZB,
            Dkth_p90_ZA = dk90_ZA, Dkth_p90_ZB = dk90_ZB
          )
        }
      }
    }
  }

  bind_rows(results)
}

# ==============================================================================
# [B2] Paired t-test per (Method, Scale), with a FIXED BH family so FDRs are
# comparable across the sweep (that comparison is the deliverable).
# ==============================================================================
compute_paired_stats <- function(raw_df, primary_method = "radius", primary_scale = 200,
                                 primary_trim = 0.10) {
  base <- raw_df %>%
    filter(!is.na(.data$Mean_ZA) & !is.na(.data$Mean_ZB)) %>%
    group_by(Method, Scale, Trim, Sender, Receiver, Interaction, Ligand, Receptor,
             Pathway, Annotation) %>%
    summarise(
      N_Donors = n(),
      Avg_ZA = mean(Mean_ZA), Avg_ZB = mean(Mean_ZB),
      Mean_Delta = mean(Delta_ZA_ZB), SD_Delta = sd(Delta_ZA_ZB),
      Mean_Nnbr_ZA = mean(N_nbr_ZA, na.rm = TRUE),
      Mean_Nnbr_ZB = mean(N_nbr_ZB, na.rm = TRUE),
      Mean_Frac0_ZA = mean(Frac0_ZA, na.rm = TRUE),
      Mean_Frac0_ZB = mean(Frac0_ZB, na.rm = TRUE),
      Any_Trim_Straddle = any(Trim_Straddle, na.rm = TRUE),
      Mean_Dkth_ZA = mean(Dkth_med_ZA, na.rm = TRUE),
      Mean_Dkth_ZB = mean(Dkth_med_ZB, na.rm = TRUE),
      P_Value = if (n() >= 3) {
        tryCatch(t.test(Mean_ZA, Mean_ZB, paired = TRUE)$p.value,
                 error = function(e) NA_real_)
      } else NA_real_,
      .groups = "drop"
    )

  m_fixed <- base %>%
    filter(Method == primary_method, Scale == primary_scale, Trim == primary_trim,
           !is.na(P_Value)) %>%
    nrow()

  base %>%
    group_by(Method, Scale, Trim) %>%
    mutate(
      N_Tests_Testable = sum(!is.na(P_Value)),
      N_Tests_In_Family = max(m_fixed, N_Tests_Testable[1]),
      FDR = p.adjust(P_Value, method = "BH", n = N_Tests_In_Family[1])
    ) %>%
    ungroup() %>%
    arrange(Method, Scale, desc(Trim), FDR)
}

# Trim concordance: for each (Method, Scale), how well does the primary
# (trimmed) estimator agree with the untrimmed one? This is the table the
# rebuttal quotes — "results are unchanged with and without trimming" —
# and it isolates any disagreement onto the threshold-straddling pairs.
build_trim_concordance <- function(stats_df, primary_trim, fdr_thresh = 0.05) {
  key <- c("Method", "Scale", "Sender", "Receiver", "Interaction", "Ligand", "Receptor")
  trims <- unique(stats_df$Trim)
  others <- setdiff(trims, primary_trim)
  if (length(others) == 0) return(tibble())

  bind_rows(lapply(others, function(ot) {
    a <- stats_df %>% filter(Trim == primary_trim) %>%
      dplyr::select(all_of(key), Mean_Delta, FDR, Any_Trim_Straddle)
    b <- stats_df %>% filter(Trim == ot) %>%
      dplyr::select(all_of(key), Mean_Delta, FDR)
    inner_join(a, b, by = key, suffix = c("_ref", "_alt")) %>%
      group_by(Method, Scale) %>%
      summarise(
        Trim_Primary = primary_trim, Trim_Alt = ot,
        N_Compared = n(),
        N_Sig_Primary = sum(FDR_ref <= fdr_thresh, na.rm = TRUE),
        N_Sig_Alt     = sum(FDR_alt <= fdr_thresh, na.rm = TRUE),
        N_Sig_Both    = sum(FDR_ref <= fdr_thresh & FDR_alt <= fdr_thresh, na.rm = TRUE),
        Jaccard_Sig = ifelse(
          (N_Sig_Primary + N_Sig_Alt - N_Sig_Both) > 0,
          N_Sig_Both / (N_Sig_Primary + N_Sig_Alt - N_Sig_Both), NA_real_),
        Sign_Agree_Frac = mean(sign(Mean_Delta_ref) == sign(Mean_Delta_alt), na.rm = TRUE),
        Spearman_Delta = suppressWarnings(
          stats::cor(Mean_Delta_ref, Mean_Delta_alt, method = "spearman",
                     use = "complete.obs")),
        # disagreements attributable to the trim threshold
        N_Primary_Only = sum(FDR_ref <= fdr_thresh & !(FDR_alt <= fdr_thresh), na.rm = TRUE),
        N_Primary_Only_Straddle = sum(FDR_ref <= fdr_thresh & !(FDR_alt <= fdr_thresh) &
                                        Any_Trim_Straddle, na.rm = TRUE),
        .groups = "drop"
      )
  })) %>% arrange(Method, Scale)
}

build_stability_table <- function(stats_df) {
  key <- c("Sender", "Receiver", "Interaction", "Ligand", "Receptor", "Pathway", "Annotation")
  tag <- function(m, s, tr) paste0(m, "_", s, "_trim", format(tr))
  delta_wide <- stats_df %>%
    mutate(col = paste0("Delta_", tag(Method, Scale, Trim))) %>%
    dplyr::select(all_of(key), col, Mean_Delta) %>%
    tidyr::pivot_wider(names_from = col, values_from = Mean_Delta,
                       values_fn = function(x) x[1])
  fdr_wide <- stats_df %>%
    mutate(col = paste0("FDR_", tag(Method, Scale, Trim))) %>%
    dplyr::select(all_of(key), col, FDR) %>%
    tidyr::pivot_wider(names_from = col, values_from = FDR,
                       values_fn = function(x) x[1])
  delta_wide %>% left_join(fdr_wide, by = key)
}

# ==============================================================================
# [B5] Targeted permutation null. Permutations drawn ONCE per (donor, sender
# type) and reused across every grid cell -> across-scale comparison is paired
# under common random numbers, and Stage 2 does not perturb Stage-1 nulls.
# ==============================================================================
run_targeted_null_for_donor <- function(
    d_name, d_data, target_df, complex_db,
    scale_grid, radii_um, knn_k, trim,
    n_shuffles, shuffle_seed, donor_index, zone_a, zone_b,
    min_senders_per_zone, logger) {

  zones <- d_data$zones
  zA <- which(zones == zone_a); zB <- which(zones == zone_b)
  zoned_all <- c(zA, zB)

  n_cells <- nrow(d_data$meta)
  cell_type_vec <- make_cell_type_vec(d_data)

  unique_sender_types <- sort(unique(target_df$Sender))
  r_g2l <- list()
  for (rtype in unique(target_df$Receiver)) {
    if (!(rtype %in% names(d_data$type_indices))) next
    rr <- d_data$type_indices[[rtype]]
    m <- integer(n_cells); m[rr] <- seq_along(rr)
    r_g2l[[rtype]] <- m
  }

  # Derive the work from the GRID actually passed in (supports additive re-runs
  # that compute only the missing cells).
  radius_scales <- sort(unique(scale_grid$Scale[scale_grid$Method == "radius"]))
  knn_scales    <- sort(unique(scale_grid$Scale[scale_grid$Method == "knn"]))
  need_radius <- length(radius_scales) > 0
  do_knn      <- length(knn_scales) > 0
  max_radius  <- if (need_radius) max(radius_scales) else NA_real_
  max_k       <- if (do_knn) max(knn_scales) else NA_integer_

  results <- list()

  for (s_type in unique_sender_types) {
    if (!(s_type %in% names(d_data$type_indices))) next
    s_rows_all <- d_data$type_indices[[s_type]]
    if (length(s_rows_all) < 10) next

    s_rows <- s_rows_all[s_rows_all %in% zoned_all]
    isA <- s_rows %in% zA; isB <- s_rows %in% zB
    nA <- sum(isA); nB <- sum(isB)
    if (nA < min_senders_per_zone || nB < min_senders_per_zone) next
    n_zoned <- length(s_rows)

    # [B5] grid-independent seed; one permutation set reused for all scales
    set.seed(shuffle_seed + donor_index * 1000L + match(s_type, unique_sender_types))
    W <- matrix(-1 / nB, nrow = n_shuffles, ncol = n_zoned)
    for (p in seq_len(n_shuffles)) {
      perm <- sample.int(n_zoned)
      W[p, perm[seq_len(nA)]] <- 1 / nA
    }

    # entity expression hoisted out of the receiver loop [P6]
    req_all <- unique(c(target_df$Ligand[target_df$Sender == s_type],
                        target_df$Receptor[target_df$Sender == s_type]))
    entity_expr <- computeExpr_LR(req_all, d_data$data.use, complex_db)
    rownames(entity_expr) <- req_all

    # [B4 FIX] Scoring drops an interaction whose ligand has ZERO expression in
    # THIS sender type. The null did not, so those donors were still counted ->
    # Donors_Tested could exceed N_Donors and deflate Robust_95.
    lig_v_s  <- unique(target_df$Ligand[target_df$Sender == s_type])
    lig_nz_s <- setNames(
      Matrix::rowSums(entity_expr[lig_v_s, s_rows_all, drop = FALSE]) > 0, lig_v_s)

    frnn <- NULL
    if (need_radius) {
      logger$log_msg("Donor ", d_name, ": NULL frNN sender ", s_type,
                     " (", n_zoned, " zoned senders, ", n_shuffles, " shuffles)")
      frnn <- frnn_for_sender(d_data, s_rows, max_radius)
    }

    for (r_type in unique(target_df$Receiver[target_df$Sender == s_type])) {
      if (!(r_type %in% names(d_data$type_indices))) next
      r_rows <- d_data$type_indices[[r_type]]
      if (length(r_rows) < 10) next

      pair_int <- target_df %>% filter(Sender == s_type, Receiver == r_type)
      keep_lig <- unname(lig_nz_s[pair_int$Ligand]); keep_lig[is.na(keep_lig)] <- FALSE
      pair_int <- pair_int[keep_lig, , drop = FALSE]
      if (nrow(pair_int) == 0) next

      restricted <- if (need_radius)
        restrict_frnn_to_receiver(frnn, cell_type_vec, r_type, r_g2l[[r_type]]) else NULL
      req_rec <- unique(pair_int$Receptor)
      Rmat <- build_Rmat(entity_expr, req_rec, r_rows)

      knn_pair <- NULL
      if (do_knn) {
        knn_pair <- knn_for_pair(d_data, s_rows, r_rows,
                                 sender_eq_receiver = (s_type == r_type), max_k = max_k)
      }

      ms <- scale_grid %>% distinct(Method, Scale)

      for (mi in seq_len(nrow(ms))) {
        method <- ms$Method[mi]
        scale  <- ms$Scale[mi]

        if (method == "radius") {
          nn_local <- radius_subset(restricted, scale)
        } else {
          if (is.null(knn_pair) || scale > knn_pair$k_feasible) next
          nn_local <- lapply(knn_pair$id, function(v) v[seq_len(min(scale, length(v)))])
        }

        Lk <- entity_expr[pair_int$Ligand, s_rows, drop = FALSE]
        trims_here <- sort(unique(scale_grid$Trim[scale_grid$Method == method &
                                                    scale_grid$Scale == scale]),
                           decreasing = TRUE)

        for (trim_v in trims_here) {
        potentials_mat <- if (trim_v == 0) {
          potentials_all(nn_local, Rmat)
        } else {
          pm <- matrix(0, nrow = length(nn_local), ncol = length(req_rec),
                       dimnames = list(NULL, req_rec))
          for (rc in req_rec) {
            r_vec <- as.numeric(entity_expr[rc, r_rows])
            if (sum(r_vec) == 0) next
            pm[, rc] <- potential_from_neighbors(nn_local, r_vec, trim = trim_v)
          }
          pm
        }

        # scores [n_zoned x K]
        score_matrix <- t(Lk) * potentials_mat[, pair_int$Receptor, drop = FALSE]

        real_delta <- colMeans(score_matrix[isA, , drop = FALSE]) -
                      colMeans(score_matrix[isB, , drop = FALSE])
        null_deltas <- as.matrix(W %*% score_matrix)     # [n_shuffles x K]

        abs_real <- abs(real_delta)
        for (k in seq_len(nrow(pair_int))) {
          abs_null <- abs(null_deltas[, k])
          results[[length(results) + 1L]] <- tibble(
            Donor = d_name, Method = method, Scale = scale, Trim = trim_v,
            Sender = s_type, Receiver = r_type,
            Interaction = pair_int$Interaction[k],
            Ligand = pair_int$Ligand[k], Receptor = pair_int$Receptor[k],
            Pathway = pair_int$Pathway[k], Annotation = pair_int$Annotation[k],
            Real_Delta = unname(real_delta[k]),
            Null_Mean = mean(abs_null), Null_SD = sd(null_deltas[, k]),
            Null_75_Abs = as.numeric(stats::quantile(abs_null, 0.75, names = FALSE)),
            Null_90_Abs = as.numeric(stats::quantile(abs_null, 0.90, names = FALSE)),
            Null_95_Abs = as.numeric(stats::quantile(abs_null, 0.95, names = FALSE)),
            Null_99_Abs = as.numeric(stats::quantile(abs_null, 0.99, names = FALSE)),
            # [G5] (1+r)/(1+n): never reports exactly 0
            Empirical_P = (1 + sum(abs_null >= abs_real[k])) / (n_shuffles + 1),
            N_Shuffles = n_shuffles
          )
        }
        }
      }
    }
  }

  bind_rows(results)
}

summarize_null_robustness <- function(null_df, stats_df) {
  robust <- null_df %>%
    group_by(Method, Scale, Trim, Sender, Receiver, Interaction, Ligand, Receptor,
             Pathway, Annotation) %>%
    summarise(
      Donors_Tested = n(),
      Pass_75 = sum(abs(Real_Delta) > Null_75_Abs, na.rm = TRUE),
      Pass_90 = sum(abs(Real_Delta) > Null_90_Abs, na.rm = TRUE),
      Pass_95 = sum(abs(Real_Delta) > Null_95_Abs, na.rm = TRUE),
      Pass_99 = sum(abs(Real_Delta) > Null_99_Abs, na.rm = TRUE),
      Robust_75 = Pass_75 / Donors_Tested,
      Robust_90 = Pass_90 / Donors_Tested,
      Robust_95 = Pass_95 / Donors_Tested,
      Robust_99 = Pass_99 / Donors_Tested,
      Mean_Real_Delta = mean(Real_Delta, na.rm = TRUE),
      # [G5] median-of-p is NOT a combined p-value; report it as a description
      # and give a proper Fisher combination alongside.
      Median_Donor_Empirical_P = median(Empirical_P, na.rm = TRUE),
      Fisher_Combined_P = {
        p <- Empirical_P[!is.na(Empirical_P)]
        if (length(p) < 2) NA_real_
        else stats::pchisq(-2 * sum(log(p)), df = 2 * length(p), lower.tail = FALSE)
      },
      .groups = "drop"
    )

  robust %>%
    left_join(
      stats_df %>%
        dplyr::select(Method, Scale, Trim, Sender, Receiver, Interaction, Ligand, Receptor,
                      Pathway, Annotation, N_Donors, Mean_Delta, SD_Delta, P_Value, FDR),
      by = c("Method", "Scale", "Trim", "Sender", "Receiver", "Interaction",
             "Ligand", "Receptor", "Pathway", "Annotation")
    ) %>%
    arrange(Method, Scale, desc(Trim), FDR, desc(Robust_95))
}

# ==============================================================================
# MAIN
# ==============================================================================
main <- function() {

  if (!exists("donor_objs")) stop("donor_objs not loaded")
  if (!dir.exists(CONFIG$output_dir)) dir.create(CONFIG$output_dir, recursive = TRUE)

  logger <- make_logger(CONFIG$output_dir, CONFIG$log_file, CONFIG$status_file)
  scale_grid <- build_scale_grid(CONFIG)
  CFG_KEY <- make_cfg_key(CONFIG)

  logger$log_msg("=== START PIPELINE v3.1 (post-audit) ===")
  logger$log_msg("Output dir: ", CONFIG$output_dir)
  logger$log_msg("Scale grid (", nrow(scale_grid), " cells): ",
                 paste(sprintf("%s:%s/trim%s", scale_grid$Method, scale_grid$Scale,
                               format(scale_grid$Trim)), collapse = ", "))
  logger$log_msg("Primary cell: ", CONFIG$primary_method, ":", CONFIG$primary_radius,
                 "/trim", format(CONFIG$primary_trim), " | config key=", CFG_KEY)
  for (tv in setdiff(unique(scale_grid$Trim), 0)) {
    logger$log_msg("trim=", format(tv), " switches on at n>=", trim_threshold_n(tv),
                   " neighbours; pairs straddling that between zones are flagged ",
                   "Trim_Straddle=TRUE.")
  }

  qc_dir <- file.path(CONFIG$output_dir, "00_qc")
  if (!dir.exists(qc_dir)) dir.create(qc_dir, recursive = TRUE)
  writeLines(capture.output(sessionInfo()), file.path(qc_dir, "sessionInfo.txt"))
  readr::write_csv(scale_grid, file.path(qc_dir, "scale_grid.csv"))
  readr::write_csv(tibble(
    Parameter = names(CONFIG),
    Value = vapply(CONFIG, function(x) paste(x, collapse = ", "), character(1))
  ), file.path(qc_dir, paste0("pipeline_config_", CFG_KEY, ".csv")))

  # ---- Donors ----
  donors_all <- names(donor_objs)
  if (is.null(donors_all) || length(donors_all) == 0)
    stop("donor_objs must be a *named* list (of Seurat objects) or a *named* character vector of .qs paths")
  donor_objs_use <- if ("s20" %in% donors_all) donor_objs[setdiff(donors_all, "s20")] else donor_objs
  if (length(donor_objs_use) != 18) {
    stop("Expected 18 donors after excluding s20. Found: ", length(donor_objs_use),
         " (", paste(names(donor_objs_use), collapse = ","), ")")
  }
  logger$log_msg("Donors: ", length(donor_objs_use), " (s20 excluded)")
  writeLines(names(donor_objs_use), file.path(qc_dir, "donors_used.txt"))

  audit_inputs(donor_objs_use, CONFIG$output_dir, logger)

  cellchat_db <- load_cellchatdb_human(CONFIG$output_dir, logger)
  db <- prepare_interaction_db(cellchat_db, CONFIG$valid_annotations, CONFIG$output_dir, logger)

  # =========================================================================
  # Scoring
  # =========================================================================
  raw_dir <- file.path(CONFIG$output_dir, "01_raw_scores")
  if (!dir.exists(raw_dir)) dir.create(raw_dir, recursive = TRUE)
  logger$log_msg("--- Scoring across scale grid ---")

  for (d_name in names(donor_objs_use)) {
    donor_out <- file.path(raw_dir, sprintf("%s_raw_scores_%s.csv", d_name, CFG_KEY))

    # ADDITIVE checkpointing: compute only the (Method, Scale, Trim) cells that
    # this donor's existing file does not already contain, then append. Adding an
    # arm (e.g. the Stage-2 kNN grid) therefore does NOT redo the radius sweep.
    prev_rows <- NULL
    grid_todo <- scale_grid
    if (file.exists(donor_out)) {
      prev_rows <- tryCatch(
        suppressMessages(readr::read_csv(donor_out, show_col_types = FALSE)),
        error = function(e) NULL)
      if (!is.null(prev_rows) && nrow(prev_rows) > 0 &&
          all(c("Method", "Scale", "Trim") %in% colnames(prev_rows))) {
        have <- prev_rows %>%
          mutate(Scale = as.numeric(Scale), Trim = as.numeric(Trim)) %>%
          distinct(Method, Scale, Trim)
        grid_todo <- dplyr::anti_join(
          scale_grid %>% mutate(Scale = as.numeric(Scale), Trim = as.numeric(Trim)),
          have, by = c("Method", "Scale", "Trim"))
      } else {
        prev_rows <- NULL
      }
    }
    if (nrow(grid_todo) == 0) {
      logger$log_msg("Skipping donor ", d_name, " (checkpoint covers grid, key ", CFG_KEY, ")")
      next
    }
    logger$log_msg("Starting donor ", d_name, " — computing ", nrow(grid_todo), " of ",
                   nrow(scale_grid), " grid cells: ",
                   paste(sprintf("%s:%s/t%s", grid_todo$Method, grid_todo$Scale,
                                 format(grid_todo$Trim)), collapse = ", "))

    obj_i <- .resolve_donor(donor_objs_use[[d_name]])
    d_data <- prep_donor_data(obj_i, db$required_genes,
                              CONFIG$cell_types, CONFIG$valid_zones, logger, d_name)
    rm(obj_i); invisible(gc(verbose = FALSE))
    if (is.null(d_data)) next

    donor_db <- filter_interactions_per_donor(d_data, db$interaction, db$entity_to_genes,
                                              CONFIG$min_pct_expr, logger, d_name)

    donor_raw <- score_donor_zone_comparison(
      d_name = d_name, d_data = d_data, interaction_db_donor = donor_db,
      complex_db = db$complex, scale_grid = grid_todo,
      radii_um = CONFIG$radii_um, knn_k = CONFIG$knn_k, trim = NULL,
      cell_types = CONFIG$cell_types, zone_a = CONFIG$zone_a, zone_b = CONFIG$zone_b,
      min_senders_per_zone = CONFIG$min_senders_per_zone, logger = logger)

    if (!is.null(prev_rows) && nrow(donor_raw) > 0) donor_raw <- bind_rows(prev_rows, donor_raw)
    else if (!is.null(prev_rows)) donor_raw <- prev_rows
    write_csv_atomic(donor_raw, donor_out)
    logger$log_msg("CHECKPOINT: ", basename(donor_out), " (rows=", nrow(donor_raw), ")")
    rm(d_data, donor_raw); invisible(gc(verbose = FALSE))
  }

  # [B3] combine BY NAME, not by glob
  logger$log_msg("--- Combining donor raw outputs ---")
  raw_files <- file.path(raw_dir, sprintf("%s_raw_scores_%s.csv",
                                          names(donor_objs_use), CFG_KEY))
  raw_files <- raw_files[file.exists(raw_files)]
  if (length(raw_files) == 0) stop("No donor raw score files found in ", raw_dir)
  raw_all <- suppressMessages(bind_rows(lapply(raw_files, readr::read_csv, show_col_types = FALSE)))

  missing_donors <- setdiff(names(donor_objs_use), unique(raw_all$Donor))
  if (length(missing_donors) > 0) {
    logger$log_msg("Donors with no scored rows: ", paste(missing_donors, collapse = ","),
                   level = "WARN")
  }

  raw_all_path <- file.path(CONFIG$output_dir, paste0("01_raw_scores_all_donors_", CFG_KEY, ".csv"))
  write_csv_atomic(raw_all, raw_all_path)
  logger$log_msg("CHECKPOINT: combined raw scores (rows=", nrow(raw_all), ")")

  # Neighbor-density / zero-partner diagnostics
  nbr_summary <- raw_all %>%
    distinct(Donor, Method, Scale, Trim, Sender, Receiver,
             N_nbr_ZA, N_nbr_ZB, Frac0_ZA, Frac0_ZB, Trim_Straddle,
             Dkth_med_ZA, Dkth_med_ZB) %>%
    group_by(Method, Scale, Trim, Sender, Receiver) %>%
    summarise(
      Mean_Nnbr_ZA = mean(N_nbr_ZA, na.rm = TRUE), Mean_Nnbr_ZB = mean(N_nbr_ZB, na.rm = TRUE),
      SD_Nnbr_ZA = sd(N_nbr_ZA, na.rm = TRUE), SD_Nnbr_ZB = sd(N_nbr_ZB, na.rm = TRUE),
      Mean_Frac0_ZA = mean(Frac0_ZA, na.rm = TRUE), Mean_Frac0_ZB = mean(Frac0_ZB, na.rm = TRUE),
      Mean_Dkth_ZA = mean(Dkth_med_ZA, na.rm = TRUE), Mean_Dkth_ZB = mean(Dkth_med_ZB, na.rm = TRUE),
      # [G7] for the kNN arm neighbor count is fixed by construction
      Nnbr_Fixed_By_Design = first(Method) == "knn",
      # Donors where this pair straddles the trim threshold between zones
      N_Donors_Trim_Straddle = sum(Trim_Straddle, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(Method, Scale, desc(Trim), Sender, Receiver)
  readr::write_csv(nbr_summary, file.path(CONFIG$output_dir, "01_neighbor_density_summary.csv"))
  logger$log_msg("Wrote neighbor-density + zero-partner summary")

  # Suggested k for the Stage-2 kNN arm
  prim <- nbr_summary %>% filter(Method == "radius", Scale == CONFIG$primary_radius,
                                 Trim == CONFIG$primary_trim)
  if (nrow(prim) > 0) {
    prim_pair <- prim %>%
      mutate(Mean_Nnbr = rowMeans(cbind(Mean_Nnbr_ZA, Mean_Nnbr_ZB), na.rm = TRUE))
    readr::write_csv(
      prim_pair %>% group_by(Receiver) %>%
        summarise(Median_Nnbr_at_primary = median(Mean_Nnbr, na.rm = TRUE),
                  Min_Nnbr = min(Mean_Nnbr, na.rm = TRUE),
                  Max_Nnbr = max(Mean_Nnbr, na.rm = TRUE), .groups = "drop"),
      file.path(CONFIG$output_dir, "01_suggested_k_by_receiver.csv"))
    sk <- max(3, round(median(prim_pair$Mean_Nnbr, na.rm = TRUE)))
    logger$log_msg("SUGGESTED k ~ ", sk, ". For Stage 2 run a VECTOR spanning the ",
                   "per-receiver spread, e.g.  Sys.setenv(LR_KNN_K='3,10,20,40')  then re-source.")
  }

  # =========================================================================
  # Stats
  # =========================================================================
  logger$log_msg("--- Paired t-tests per scale (fixed BH family) ---")
  stats_df <- compute_paired_stats(raw_all, CONFIG$primary_method, CONFIG$primary_radius,
                                   CONFIG$primary_trim)

  stats_dir <- file.path(CONFIG$output_dir, "02_stats")
  if (!dir.exists(stats_dir)) dir.create(stats_dir, recursive = TRUE)
  write_csv_atomic(stats_df, file.path(stats_dir, paste0("paired_t_stats_by_scale_", CFG_KEY, ".csv")))
  write_csv_atomic(build_stability_table(stats_df),
                   file.path(stats_dir, paste0("delta_fdr_stability_wide_", CFG_KEY, ".csv")))
  logger$log_msg("CHECKPOINT: stats (rows=", nrow(stats_df), ")")

  # Trim concordance — the "unchanged with and without trimming" table
  conc <- build_trim_concordance(stats_df, CONFIG$primary_trim, CONFIG$fdr_for_null)
  if (nrow(conc) > 0) {
    write_csv_atomic(conc, file.path(stats_dir, paste0("trim_concordance_", CFG_KEY, ".csv")))
    pc <- conc %>% filter(Method == CONFIG$primary_method, Scale == CONFIG$primary_radius)
    if (nrow(pc) > 0) {
      logger$log_msg("TRIM CONCORDANCE at primary window (", CONFIG$primary_radius, "um): ",
                     "sig trim", format(CONFIG$primary_trim), "=", pc$N_Sig_Primary[1],
                     ", sig trim0=", pc$N_Sig_Alt[1],
                     ", both=", pc$N_Sig_Both[1],
                     ", Jaccard=", round(pc$Jaccard_Sig[1], 3),
                     ", sign-agreement=", round(pc$Sign_Agree_Frac[1], 3),
                     ", Spearman(Delta)=", round(pc$Spearman_Delta[1], 3))
      if (!is.na(pc$N_Primary_Only[1]) && pc$N_Primary_Only[1] > 0) {
        logger$log_msg("  ", pc$N_Primary_Only_Straddle[1], " of ", pc$N_Primary_Only[1],
                       " trim-only hits are threshold-STRADDLING pairs (inspect these).")
      }
    }
  }

  # =========================================================================
  # Targeted null: discovery at the primary window, evaluated at every scale
  # =========================================================================
  sig_primary <- stats_df %>%
    filter(Method == CONFIG$primary_method, Scale == CONFIG$primary_radius,
           Trim == CONFIG$primary_trim,
           !is.na(FDR), FDR <= CONFIG$fdr_for_null, N_Donors >= 3)

  null_dir <- file.path(CONFIG$output_dir, "03_null")
  if (!dir.exists(null_dir)) dir.create(null_dir, recursive = TRUE)
  readr::write_csv(sig_primary, file.path(null_dir, "significant_at_primary_window.csv"))
  logger$log_msg("Discovery hits at primary window: ", nrow(sig_primary))

  if (nrow(sig_primary) == 0) {
    logger$log_msg("No significant interactions; skipping null.", level = "WARN")
    logger$log_msg("=== PIPELINE COMPLETE ===")
    return(invisible(NULL))
  }

  target_df <- sig_primary %>%
    dplyr::transmute(Sender, Receiver, Interaction, Ligand, Receptor, Pathway, Annotation) %>%
    distinct()
  readr::write_csv(target_df, file.path(null_dir, "target_interactions_for_null.csv"))
  TGT_KEY <- cfg_hash(list(target_df$Sender, target_df$Receiver, target_df$Interaction))

  null_by_donor_dir <- file.path(null_dir, "by_donor")
  if (!dir.exists(null_by_donor_dir)) dir.create(null_by_donor_dir, recursive = TRUE)
  logger$log_msg("--- Targeted null (", CONFIG$n_shuffles, " shuffles x ",
                 nrow(scale_grid), " scales) ---")

  donor_names_sorted <- sort(names(donor_objs_use))
  for (di in seq_along(donor_names_sorted)) {
    d_name <- donor_names_sorted[di]
    donor_null_out <- file.path(null_by_donor_dir,
                                sprintf("%s_null_%s_%s.csv", d_name, CFG_KEY, TGT_KEY))
    if (checkpoint_covers_grid(donor_null_out, scale_grid)) {
      logger$log_msg("Skipping donor ", d_name, " (null checkpoint covers grid)")
      next
    }
    logger$log_msg("Null for donor ", d_name, " (index=", di, ")")

    obj_i <- .resolve_donor(donor_objs_use[[d_name]])
    d_data <- prep_donor_data(obj_i, db$required_genes,
                              CONFIG$cell_types, CONFIG$valid_zones, logger, d_name)
    rm(obj_i); invisible(gc(verbose = FALSE))
    if (is.null(d_data)) next

    # [B4] apply the SAME detection filter the scoring step used, so
    # Donors_Tested cannot exceed N_Donors
    donor_db <- filter_interactions_per_donor(d_data, db$interaction, db$entity_to_genes,
                                              CONFIG$min_pct_expr, logger, d_name)
    target_d <- dplyr::semi_join(target_df, donor_db,
                                 by = c("Ligand" = "ligand", "Receptor" = "receptor"))
    if (nrow(target_d) == 0) {
      logger$log_msg("Donor ", d_name, ": no target interactions pass detection; skipping")
      next
    }

    donor_null <- run_targeted_null_for_donor(
      d_name = d_name, d_data = d_data, target_df = target_d, complex_db = db$complex,
      scale_grid = scale_grid, radii_um = CONFIG$radii_um, knn_k = CONFIG$knn_k,
      trim = NULL, n_shuffles = CONFIG$n_shuffles, shuffle_seed = CONFIG$shuffle_seed,
      donor_index = di, zone_a = CONFIG$zone_a, zone_b = CONFIG$zone_b,
      min_senders_per_zone = CONFIG$min_senders_per_zone, logger = logger)

    write_csv_atomic(donor_null, donor_null_out)
    logger$log_msg("CHECKPOINT: ", basename(donor_null_out), " (rows=", nrow(donor_null), ")")
    rm(d_data, donor_null); invisible(gc(verbose = FALSE))
  }

  null_files <- file.path(null_by_donor_dir,
                          sprintf("%s_null_%s_%s.csv", donor_names_sorted, CFG_KEY, TGT_KEY))
  null_files <- null_files[file.exists(null_files)]
  null_all <- suppressMessages(bind_rows(lapply(null_files, readr::read_csv, show_col_types = FALSE)))
  write_csv_atomic(null_all, file.path(null_dir, paste0("targeted_null_all_", CFG_KEY, "_", TGT_KEY, ".csv")))

  robust_df <- summarize_null_robustness(null_all, stats_df)
  write_csv_atomic(robust_df, file.path(null_dir,
                                        paste0("validated_interactions_all_scales_", CFG_KEY, ".csv")))
  if (nrow(robust_df) > 0 && any(robust_df$Donors_Tested > robust_df$N_Donors, na.rm = TRUE)) {
    logger$log_msg("Donors_Tested exceeds N_Donors for some rows (audit B4 regression).",
                   level = "WARN")
  }
  logger$log_msg("Wrote robustness summary (rows=", nrow(robust_df), ")")
  logger$log_msg("=== PIPELINE COMPLETE ===")
}

# Runs on source(), as before. The option is only used by the correctness test
# harness (test_pipeline_correctness.R) to load the functions without executing.
if (!isTRUE(getOption("LR_PIPELINE_NO_RUN", FALSE))) main()
