suppressPackageStartupMessages({ library(qs); library(SeuratObject) })
FILT <- "/broad/macosko/akraft/CAP/rerun/processed_pucks/aggregated/NoReg_Final_by_donor/Final_Cell_Arrangements_Oct9/Fig6/LR_sensitivity_analysis/filtered_donor_objects"
OUT  <- "/broad/macosko/akraft/CAP/rerun/processed_pucks/aggregated/NoReg_Final_by_donor/Final_Cell_Arrangements_Oct9/Fig6/LR_sensitivity_analysis/pipeline_v3/09_coords"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
fs <- list.files(FILT, pattern="_filtered\\.qs$", full.names=TRUE)
cat("donors:", length(fs), "\n")
for (fp in fs) {
  id <- sub("_filtered\\.qs$","",basename(fp))
  s <- qs::qread(fp); m <- s@meta.data
  keep <- intersect(c("x_um2","y_um2","Spatial_Domain","Cell_types","Class_name","AK_celltype","Final_Zone_Assignments_No_Smooth"), colnames(m))
  df <- m[, keep, drop=FALSE]
  write.csv(df, file.path(OUT, paste0(id,"_coords.csv")), row.names=FALSE)
  cat(id, nrow(df), "cells\n")
}
cat("DONE\n")
