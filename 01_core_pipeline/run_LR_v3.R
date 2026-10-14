#!/usr/bin/env Rscript

# ==============================================================================
# run_LR_v3.R
#
# One-shot driver for the LR sensitivity analysis:
#   1) Build/load the filtered per-donor objects -> `donor_objs`
#   2) Source the v3 pipeline (radius sweep; kNN arm added in Stage 2)
#
# Run from the directory that contains both scripts, e.g. inside tmux on a
# bigmem node, or via sbatch:
#   Rscript run_LR_v3.R
#
# First run builds the filtered cache (slow); later runs reload it (fast).
# ==============================================================================

CODE_DIR <- "/broad/macosko/akraft/CAP/rerun/processed_pucks/aggregated/NoReg_Final_by_donor/Final_Cell_Arrangements_Oct9/Fig6"

# 1) Filtered donors -> donor_objs
source(file.path(CODE_DIR, "build_filtered_donor_objs.R"))
donor_objs <- load_donor_objs()          # add rebuild = TRUE to re-filter from raw

# 2) Run the pipeline (its trailing main() executes on source)
source(file.path(CODE_DIR, "March_LR_rerun_pipeline_v3.R"))
