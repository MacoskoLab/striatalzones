options(LR_PIPELINE_NO_RUN = TRUE)
CODE_DIR <- "/broad/macosko/akraft/CAP/rerun/processed_pucks/aggregated/NoReg_Final_by_donor/Final_Cell_Arrangements_Oct9/Fig6"
source(file.path(CODE_DIR, "build_filtered_donor_objs.R")); donor_objs <- load_donor_objs()
source(file.path(CODE_DIR, "March_LR_rerun_pipeline_v3.R"))
source(file.path(CODE_DIR, "radius600.R"))
