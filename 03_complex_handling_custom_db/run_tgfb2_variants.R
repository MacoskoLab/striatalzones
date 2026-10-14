#!/usr/bin/env Rscript
# Targeted: score every TGFB2 receptor definition at r=200 (full-complex geomean
# for multi-subunit; single-gene otherwise), Zone1-vs-Zone4, same null as the
# main pipeline. Overrides prepare_interaction_db to a 6-row DB and audit_inputs
# to a no-op; forces ALL scored interactions through the null (fdr_for_null=1) so
# every variant gets a robustness value regardless of discovery FDR.
suppressWarnings(suppressMessages({
  library(Seurat); library(qs); library(data.table)
  library(dplyr); library(tidyr); library(readr); library(tibble); library(Matrix)
}))
FIG6 <- "/broad/macosko/akraft/CAP/rerun/processed_pucks/aggregated/NoReg_Final_by_donor/Final_Cell_Arrangements_Oct9/Fig6"
LRDIR <- file.path(FIG6, "LR_sensitivity_analysis")
options(LR_PIPELINE_NO_RUN = TRUE)
source(file.path(LRDIR, "LR_pipeline_radius_only.R"))

out <- file.path(LRDIR, "explode_compare", "tgfb2_variants")
CONFIG$output_dir <- out; CONFIG$radii_um <- c(200); CONFIG$primary_radius <- 200
CONFIG$fdr_for_null <- 1.0            # null every scored interaction
dir.create(file.path(out,"00_db"), recursive = TRUE, showWarnings = FALSE)
src_rda <- file.path(LRDIR,"pipeline_v3","00_db","CellChatDB.human.rda")
if (file.exists(src_rda)) file.copy(src_rda, file.path(out,"00_db","CellChatDB.human.rda"), overwrite = FALSE)

# no-op audit (skip the upfront full donor load)
audit_inputs <- function(donor_objects_list, output_dir, logger) invisible(NULL)

# tiny DB: all TGFB2 receptor definitions
prepare_interaction_db <- function(cellchat_db, valid_annotations, output_dir, logger) {
  complex <- cellchat_db$complex
  scn <- grep("^subunit", colnames(complex), value = TRUE)
  addc <- function(complex, nm, subs) {
    if (nm %in% rownames(complex)) return(complex)
    r <- complex[1,,drop=FALSE]; r[] <- ""
    for (j in seq_along(subs)) r[[scn[j]]] <- subs[j]
    rownames(r) <- nm; rbind(complex, r)
  }
  complex <- addc(complex, "TGFBR1_TGFBR3", c("TGFBR1","TGFBR3"))
  complex <- addc(complex, "TGFBR2_TGFBR3", c("TGFBR2","TGFBR3"))
  # canonical TGFBR1+TGFBR2 is stock complex "TGFbR1_R2"
  interaction <- data.frame(
    interaction_name = c("TGFB2_TGFBR1_TGFBR3","TGFB2_TGFBR2_TGFBR3","TGFB2_TGFBR1_TGFBR2",
                         "TGFB2_TGFBR3","TGFB2_TGFBR2","TGFB2_TGFBR1"),
    ligand   = "TGFB2",
    receptor = c("TGFBR1_TGFBR3","TGFBR2_TGFBR3","TGFbR1_R2","TGFBR3","TGFBR2","TGFBR1"),
    annotation = "Secreted Signaling", pathway_name = "TGFb", stringsAsFactors = FALSE)

  complex_names <- rownames(complex); subunit_cols <- scn
  get_subunits <- function(e) { if (!(e %in% complex_names)) return(e)
    v <- unlist(complex[e, subunit_cols, drop=TRUE], use.names=FALSE); v <- v[!is.na(v)&v!=""]; unique(v) }
  entities <- unique(c(interaction$ligand, interaction$receptor))
  entity_to_genes <- lapply(entities, get_subunits); names(entity_to_genes) <- entities
  required_genes <- unique(unlist(entity_to_genes, use.names=FALSE))
  logger$log_msg("TGFB2-variant DB: ", nrow(interaction), " interactions; genes: ",
                 paste(required_genes, collapse=","))
  for (e in entities) logger$log_msg("  entity ", e, " -> ", paste(entity_to_genes[[e]], collapse="+"))
  list(interaction = interaction, complex = complex, entities = entities,
       entity_to_genes = entity_to_genes, required_genes = required_genes)
}

fdir <- file.path(LRDIR, "filtered_donor_objects")
qsf <- list.files(fdir, pattern="_filtered\\.qs$", full.names=TRUE)
donor_objs <- setNames(qsf, sub("_filtered\\.qs$","", basename(qsf)))
cat("TGFB2-variant run: out=", out, " donors=", length(donor_objs), "\n", sep="")
main()
cat("=== TGFB2 VARIANTS DONE ===\n")
