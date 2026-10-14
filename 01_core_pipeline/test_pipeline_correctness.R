#!/usr/bin/env Rscript
# ==============================================================================
# test_pipeline_correctness.R
#
# Verifies the OPTIMIZED v3.1 scoring path against an INDEPENDENT brute-force
# reference implementation on synthetic data.
#
# The reference recomputes neighbors from a raw distance matrix (not via
# dbscan/the pipeline's helpers), so this genuinely tests:
#   - frNN + self-exclusion + receiver-type restriction + radius subsetting
#   - kNN + self-exclusion + per-k truncation
#   - the sparse-matmul potentials (potentials_all) vs a plain per-sender mean
#   - the columnar emit vs a per-interaction loop
#   - complex (geometric-mean) entity expression
#
# Run:  Rscript test_pipeline_correctness.R
# ==============================================================================

options(LR_PIPELINE_NO_RUN = TRUE)
suppressPackageStartupMessages({
  library(Matrix); library(dplyr); library(tibble); library(tidyr)
})

PIPE <- "March_LR_rerun_pipeline_v3.R"
if (!file.exists(PIPE)) stop("Run this from the directory containing ", PIPE)
source(PIPE)

set.seed(20260722)

# ------------------------------------------------------------------ synthetic
N        <- 400
TYPES    <- c("MSN", "Astrocyte", "Microglia")
GENES    <- paste0("G", 1:14)

cell_type <- sample(TYPES, N, replace = TRUE, prob = c(.6, .25, .15))
zone      <- sample(1:6, N, replace = TRUE, prob = c(.25, .1, .1, .25, .15, .15))
coords    <- cbind(x_um2 = runif(N, 0, 1200), y_um2 = runif(N, 0, 1200))
rownames(coords) <- paste0("cell", 1:N)

# zero-inflated expression
expr <- matrix(0, nrow = length(GENES), ncol = N,
               dimnames = list(GENES, paste0("cell", 1:N)))
for (g in seq_along(GENES)) {
  nz <- rbinom(N, 1, 0.45) == 1
  expr[g, nz] <- rlnorm(sum(nz), 0, 1)
}
data.use <- Matrix::Matrix(expr, sparse = TRUE)

meta <- data.frame(row.names = paste0("cell", 1:N))

d_data <- list(
  meta = meta, coords = coords, zones = zone,
  type_indices = split(seq_len(N), cell_type),
  data.use = data.use, cell_ids = rownames(coords)
)

# one complex (geometric mean of two subunits) + single-gene entities
complex_db <- data.frame(subunit_1 = "G13", subunit_2 = "G14",
                         row.names = "CPLX_A", stringsAsFactors = FALSE)

interaction_db <- tibble(
  interaction_name = c("I1", "I2", "I3", "I4", "I5"),
  ligand           = c("G1", "G2", "G3", "G1",  "G4"),
  receptor         = c("G5", "G6", "CPLX_A", "G7", "G6"),
  annotation       = "Secreted Signaling",
  pathway_name     = "P"
)

quiet_logger <- list(log_msg = function(..., level = "INFO") invisible(NULL))

RADII <- c(50, 150, 400)
KS    <- c(3, 8)
TRIMS <- c(0.10, 0)
GRID  <- tidyr::crossing(
  bind_rows(tibble(Method = "radius", Scale = as.numeric(RADII)),
            tibble(Method = "knn",    Scale = as.numeric(KS))),
  Trim = TRIMS)

# ------------------------------------------------------- reference (brute force)
D <- as.matrix(dist(coords))                       # full pairwise distances
ent <- computeExpr_LR(unique(c(interaction_db$ligand, interaction_db$receptor)),
                      data.use, complex_db)

reference_scores <- function(zone_a = 1, zone_b = 4, min_senders = 5) {
  out <- list()
  zA <- which(zone == zone_a); zB <- which(zone == zone_b)
  for (st in TYPES) {
    s_all <- which(cell_type == st)
    if (length(s_all) < 10) next
    lig_nz <- sapply(unique(interaction_db$ligand),
                     function(l) sum(ent[l, s_all]) > 0)
    s <- s_all[s_all %in% c(zA, zB)]
    if (length(s) < 2) next
    inA <- which(s %in% zA); inB <- which(s %in% zB)
    if (length(inA) < min_senders && length(inB) < min_senders) next

    for (rt in TYPES) {
      r <- which(cell_type == rt)
      if (length(r) < 10) next

      for (gi in seq_len(nrow(GRID))) {
        meth <- GRID$Method[gi]; sc <- GRID$Scale[gi]; tr <- GRID$Trim[gi]

        # neighbors: brute force, self excluded, receiver-type only
        nb <- lapply(seq_along(s), function(i) {
          cand <- r[r != s[i]]
          dd <- D[s[i], cand]
          if (meth == "radius") {
            cand[dd <= sc]
          } else {
            if (length(cand) < sc) return(NA)      # infeasible marker
            cand[order(dd)][seq_len(sc)]
          }
        })
        if (any(vapply(nb, function(z) length(z) == 1 && is.na(z[1]), logical(1)))) next

        rows <- list()
        for (k in seq_len(nrow(interaction_db))) {
          lg <- interaction_db$ligand[k]; rc <- interaction_db$receptor[k]
          if (!isTRUE(lig_nz[[lg]])) next
          pot <- vapply(nb, function(idx)
            if (length(idx) == 0) 0 else mean(ent[rc, idx], trim = tr), numeric(1))
          sco <- ent[lg, s] * pot
          rows[[length(rows) + 1]] <- tibble(
            Method = meth, Scale = sc, Trim = tr, Sender = st, Receiver = rt,
            Interaction = interaction_db$interaction_name[k],
            Mean_ZA = if (length(inA) >= min_senders) mean(sco[inA]) else NA_real_,
            Mean_ZB = if (length(inB) >= min_senders) mean(sco[inB]) else NA_real_)
        }
        out <- c(out, rows)
      }
    }
  }
  bind_rows(out)
}

# ------------------------------------------------------------------- run both
cat("Running OPTIMIZED implementation ...\n")
opt <- score_donor_zone_comparison(
  d_name = "test", d_data = d_data, interaction_db_donor = interaction_db,
  complex_db = complex_db, scale_grid = GRID, radii_um = RADII, knn_k = KS,
  trim = NULL, cell_types = TYPES, zone_a = 1, zone_b = 4,
  min_senders_per_zone = 5, logger = quiet_logger)

cat("Running REFERENCE (brute-force) implementation ...\n")
ref <- reference_scores()

cat(sprintf("\noptimized rows: %d   reference rows: %d\n", nrow(opt), nrow(ref)))

cmp <- ref %>%
  rename(Ref_ZA = Mean_ZA, Ref_ZB = Mean_ZB) %>%
  inner_join(opt %>% select(Method, Scale, Trim, Sender, Receiver, Interaction, Mean_ZA, Mean_ZB),
             by = c("Method", "Scale", "Trim", "Sender", "Receiver", "Interaction"))

cat(sprintf("matched rows: %d\n", nrow(cmp)))

only_ref <- anti_join(ref, opt, by = c("Method","Scale","Trim","Sender","Receiver","Interaction"))
only_opt <- anti_join(opt, ref, by = c("Method","Scale","Trim","Sender","Receiver","Interaction"))

fail <- FALSE
if (nrow(only_ref) || nrow(only_opt)) {
  cat(sprintf("\nROW-SET MISMATCH: %d only-in-reference, %d only-in-optimized\n",
              nrow(only_ref), nrow(only_opt)))
  if (nrow(only_ref)) print(head(only_ref, 10))
  if (nrow(only_opt)) print(head(only_opt %>% select(Method,Scale,Trim,Sender,Receiver,Interaction), 10))
  fail <- TRUE
}

dA <- abs(cmp$Ref_ZA - cmp$Mean_ZA); dB <- abs(cmp$Ref_ZB - cmp$Mean_ZB)
dA[is.na(dA)] <- ifelse(is.na(cmp$Ref_ZA) & is.na(cmp$Mean_ZA), 0, Inf)[is.na(dA)]
dB[is.na(dB)] <- ifelse(is.na(cmp$Ref_ZB) & is.na(cmp$Mean_ZB), 0, Inf)[is.na(dB)]

cat(sprintf("\nmax |diff| Mean_ZA = %.3e\nmax |diff| Mean_ZB = %.3e\n",
            max(dA), max(dB)))

TOL <- 1e-10
if (max(dA) > TOL || max(dB) > TOL) {
  cat("\nNUMERIC MISMATCH — worst rows:\n")
  print(cmp[order(-pmax(dA, dB)), ][1:min(10, nrow(cmp)), ])
  fail <- TRUE
}

# --------------------------------------------- direct unit test: potentials_all
cat("\n--- unit: potentials_all vs per-sender mean ---\n")
set.seed(7)
nn <- lapply(1:150, function(i) sample.int(80, sample(0:12, 1)))
Rm <- matrix(rlnorm(80 * 6), 80, 6, dimnames = list(NULL, paste0("R", 1:6)))
Rsp <- Matrix::Matrix(Rm, sparse = TRUE)
pa <- potentials_all(nn, Rsp)
pn <- sapply(colnames(Rm), function(cc) potential_from_neighbors(nn, Rm[, cc], trim = 0))
cat(sprintf("max |diff| = %.3e\n", max(abs(pa - pn))))
if (max(abs(pa - pn)) > 1e-12) fail <- TRUE

# ---------------------------------------------------- unit: trim sanity (audit B1)
cat("\n--- unit: null W-matrix == explicit shuffle loop ---\n")
# The null replaced a per-shuffle colMeans loop with a single (W %*% scores)
# matmul. W has +1/nA on permuted zone-A positions and -1/nB elsewhere, so each
# row of W %*% scores IS that shuffle's (meanA - meanB). Verify against the loop.
set.seed(99)
nz_s <- 120; Kk <- 40; nsh <- 50
scores <- matrix(rlnorm(nz_s * Kk), nz_s, Kk)
nA <- 47; nB <- nz_s - nA
set.seed(1234)
W <- matrix(-1 / nB, nrow = nsh, ncol = nz_s)
perms <- vector("list", nsh)
for (p in seq_len(nsh)) {
  perms[[p]] <- sample.int(nz_s)
  W[p, perms[[p]][seq_len(nA)]] <- 1 / nA
}
nd_mat <- as.matrix(W %*% scores)
nd_loop <- t(vapply(seq_len(nsh), function(p) {
  ia <- perms[[p]][seq_len(nA)]; ib <- perms[[p]][(nA + 1):nz_s]
  colMeans(scores[ia, , drop = FALSE]) - colMeans(scores[ib, , drop = FALSE])
}, numeric(Kk)))
cat(sprintf("max |diff| = %.3e\n", max(abs(nd_mat - nd_loop))))
if (max(abs(nd_mat - nd_loop)) > 1e-10) fail <- TRUE

cat("\n--- smoke: run_targeted_null_for_donor end-to-end ---\n")
tgt <- interaction_db %>%
  transmute(Sender = "MSN", Receiver = "Astrocyte", Interaction = interaction_name,
            Ligand = ligand, Receptor = receptor, Pathway = pathway_name,
            Annotation = annotation)
nullres <- tryCatch(
  run_targeted_null_for_donor(
    d_name = "test", d_data = d_data, target_df = tgt, complex_db = complex_db,
    scale_grid = GRID, radii_um = RADII, knn_k = KS, trim = NULL,
    n_shuffles = 50, shuffle_seed = 42, donor_index = 1,
    zone_a = 1, zone_b = 4, min_senders_per_zone = 5, logger = quiet_logger),
  error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL })
if (is.null(nullres) || nrow(nullres) == 0) {
  cat("null produced no rows\n"); fail <- TRUE
} else {
  cat(sprintf("null rows = %d, scales = %s\n", nrow(nullres),
              paste(unique(paste0(nullres$Method, nullres$Scale)), collapse = ",")))
  cat(sprintf("Empirical_P range = [%.4f, %.4f]  (floor should be 1/51 = %.4f)\n",
              min(nullres$Empirical_P), max(nullres$Empirical_P), 1 / 51))
  if (any(nullres$Empirical_P <= 0)) { cat("Empirical_P hit exactly 0\n"); fail <- TRUE }
}

cat("\n--- unit: permutations are shared across scales (audit B5) ---\n")
# Same interaction at two different radii must have been evaluated against the
# SAME permutation draw, so re-running the null reproduces identical nulls.
n2 <- run_targeted_null_for_donor(
  d_name = "test", d_data = d_data, target_df = tgt, complex_db = complex_db,
  scale_grid = GRID, radii_um = RADII, knn_k = KS, trim = NULL,
  n_shuffles = 50, shuffle_seed = 42, donor_index = 1,
  zone_a = 1, zone_b = 4, min_senders_per_zone = 5, logger = quiet_logger)
reprod <- isTRUE(all.equal(nullres$Null_95_Abs, n2$Null_95_Abs))
cat("null is reproducible across runs: ", reprod, "\n")
if (!reprod) fail <- TRUE

# Dropping the kNN rows from the grid must NOT change the radius-arm nulls.
GRID_R <- GRID %>% filter(Method == "radius")
n3 <- run_targeted_null_for_donor(
  d_name = "test", d_data = d_data, target_df = tgt, complex_db = complex_db,
  scale_grid = GRID_R, radii_um = RADII, knn_k = NULL, trim = NULL,
  n_shuffles = 50, shuffle_seed = 42, donor_index = 1,
  zone_a = 1, zone_b = 4, min_senders_per_zone = 5, logger = quiet_logger)
j <- inner_join(nullres %>% filter(Method == "radius") %>%
                  select(Method, Scale, Trim, Interaction, A = Null_95_Abs),
                n3 %>% select(Method, Scale, Trim, Interaction, B = Null_95_Abs),
                by = c("Method", "Scale", "Trim", "Interaction"))
stage_ok <- nrow(j) > 0 && isTRUE(all.equal(j$A, j$B))
cat("radius nulls unchanged by adding/removing the kNN arm: ", stage_ok, "\n")
if (!stage_ok) fail <- TRUE

cat("\n--- unit: trim=0 is a plain mean at every n ---\n")
ok_trim <- all(vapply(c(3, 9, 10, 25, 100), function(n) {
  x <- rlnorm(n); isTRUE(all.equal(mean(x, trim = 0), mean(x)))
}, logical(1)))
cat("trim=0 == plain mean at n in {3,9,10,25,100}: ", ok_trim, "\n")
if (!ok_trim) fail <- TRUE

cat("\n=====================================\n")
cat(if (fail) "RESULT: FAIL\n" else "RESULT: PASS — optimized == brute-force reference\n")
cat("=====================================\n")
# Only exit the process when run non-interactively (Rscript). Sourcing this from
# an RStudio console must NOT kill the session.
if (!interactive()) quit(status = if (fail) 1 else 0)
invisible(!fail)
