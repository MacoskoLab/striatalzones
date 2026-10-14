# Stage 2 — add the density-controlled kNN arm (k = 5,10,20,40 at the primary trim).
# Scoring is ADDITIVE: only the new kNN cells are computed; the radius sweep is reused.
# The null fully recomputes (the grid changed), which also applies the B4 fix.
Sys.setenv(LR_KNN_K = "5,10,20,40")
CODE_DIR <- "/broad/macosko/akraft/CAP/rerun/processed_pucks/aggregated/NoReg_Final_by_donor/Final_Cell_Arrangements_Oct9/Fig6"
source(file.path(CODE_DIR, "build_filtered_donor_objs.R"))
donor_objs <- load_donor_objs()
source(file.path(CODE_DIR, "March_LR_rerun_pipeline_v3.R"))
