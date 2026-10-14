# Spatial ligand–receptor (LR) sensitivity analysis — Zone 1 vs Zone 4

Code release accompanying Kraft, Lee et al., *Mesoscale molecular architecture of the striatum across cell types and adult aging* (Cell, 2026). Nomenclature note: "MSN" in script and variable names corresponds to "SPN" (spiny projection neuron) in the paper.

Code for the spatially constrained ligand–receptor analysis of human striatal
slide-tags data (18 donors), comparing MSN/glial signaling between striatal
Zone 1 and Zone 4, including the full radius/neighbor sensitivity battery.

## Method overview

Each candidate interaction (CellChat v2 database, `CellChatDB.human`; Jin et al.
2025) is scored per sender cell as **ligand expression × mean receptor expression
among neighbor cells within a spatial radius**. Multi-subunit receptors use the
geometric mean of subunit expression (full-complex scoring). Per donor, the
Zone 1 − Zone 4 difference (Δ) of this score is computed on a 10% trimmed mean.

Statistics, in two stages:

1. **Discovery** — for each interaction, a donor-paired t-test on the per-donor
   Δ values, Benjamini–Hochberg FDR across interactions, q < 0.05.
2. **Robustness** — for each discovered interaction, a within-donor permutation
   null (zone labels shuffled within donor, 200 shuffles) yields each donor's
   95th-percentile null Δ. `Robust_95` = fraction of donors whose observed |Δ|
   exceeds their own 95th-percentile null. Reported interactions require
   **Robust_95 ≥ 0.5** (≥ half of donors beat their own null).

Sensitivity design: the primary radius is 200 µm; the sweep covers
50–1,000 µm, a fixed-neighbor-count (kNN) arm holds neighbor count constant
across cell types, and a density-matched arm assigns each sender its own
radius from its kNN saturation distance.

## Directory guide

### `01_core_pipeline/` — the main pipeline
| script | role |
|---|---|
| `March_LR_rerun_pipeline_v3.R` | v3.1 pipeline: scoring, radius sweep {50,100,200,400,1000} µm, kNN arm, paired-t + BH discovery, permutation-null robustness. All downstream arms source their functions from here (or from the radius-only variant). |
| `build_filtered_donor_objs.R` | one-time loader/filter/cache: per-donor Seurat objects → filtered `.qs` inputs (cell types, zone labels, spatial coordinates). |
| `run_LR_v3.R` | one-shot driver (build/load donor objects, source pipeline, run). |
| `run_LR_v3_stage2.R` | stage 2: adds the fixed-neighbor-count (kNN) arm (k = 5,10,20,40) additively; recomputes the null on the combined grid. |
| `test_pipeline_correctness.R` | unit tests of the scoring/stat/null machinery on synthetic data. |

### `02_sensitivity_arms/` — radius/neighbor sensitivity follow-ups
| script | role |
|---|---|
| `low_k.R` + `run_lowk.R`, `run_lowk2.R` | extends the kNN arm downward (k = 1, 3, 5). |
| `matched_radii.R` + `run_matched.R` | density-matched radius arm: each sender scored at the physical distance of its own kNN-saturation k*, rounded to 100 µm. |
| `radius600.R` + `run_radius600.R` | single fixed 600 µm run (completes the per-sender radius grid). (Header comment is inherited from `matched_radii.R`; the script runs radius = 600 µm.) |
| `extract_coords.R` | exports per-cell spatial coordinates used by the neighbor-distance summaries. |

### `03_complex_handling_custom_db/` — receptor-complex handling & curated additions
| script | role |
|---|---|
| `LR_pipeline_radius_only.R` | the v3.1 pipeline stripped to the radius arm only (identical scoring/stat/null functions); sourced by the two launchers below. |
| `run_explode_compare.R` | full-complex (geometric mean) vs exploded (per-subunit) receptor scoring at 200 µm, stock + curated database. |
| `run_tgfb2_variants.R` | targeted scoring of every TGFB2 receptor definition (including the curated `TGFB2_TGFBR3` betaglycan single-receptor and `TGFB2_TGFBR1_TGFBR3`) through the same permutation null. |

### `04_figures/` — figure panels
R scripts for the published panels: the 4-panel sensitivity figure
(`fig_4panel_combo.R`), the composite G/H + sensitivity row
(`Fig_composite_GH_sensitivity.R`), pathway-family and ephrin/semaphorin panels,
custom Hedgehog/TGFβ and MSN-ligand panels, and Zone 1/Zone 4 network/chord
maps. These consume CSV outputs of the pipeline (paths need adjusting; see
below).

## Environment

Run on a Linux compute cluster via a conda environment (`seurat_env`):

- R 4.3.3
- Seurat 5.3.0, SeuratObject
- qs 0.27.3, data.table 1.17.8, dbscan 1.1.11
- dplyr, tidyr, readr, tibble, Matrix, stringr
- Figure scripts additionally: ggplot2, patchwork, scales, ragg, extrafont

Interaction database: `CellChatDB.human` from CellChat v2 (Jin et al. 2025),
subset to the "Secreted Signaling" and "Cell-Cell Contact" categories.

## Execution order

1. `01_core_pipeline/run_LR_v3.R` (then `run_LR_v3_stage2.R` for the kNN arm)
2. Sensitivity arms in `02_sensitivity_arms/` (each `run_*.R` is standalone;
   all source the core pipeline for its functions)
3. `03_complex_handling_custom_db/run_explode_compare.R` (once with
   `LR_EXPLODE=false`, once with `LR_EXPLODE=true`) and `run_tgfb2_variants.R`
4. Figure scripts in `04_figures/`

## Input schema

The pipeline runs on per-donor Seurat objects saved as `.qs`, one file per
donor, named `s<N>_with_master.qs`, in the directory set by `RAW_DIR` at the
top of `build_filtered_donor_objs.R`. Each object's `@meta.data` must carry:

| column | content |
|---|---|
| `x_um2`, `y_um2` | spatial coordinates in µm (cells with NA are dropped) |
| `Spatial_Domain` | domain/zone assignment; integer-coercible 1–6 (cells outside 1–6, incl. white matter, are dropped). The two compared domains are set in `CONFIG` (as run: 1 vs 4). |
| `Cell_types` | cell type per cell; used directly as the analysis type. The sender/receiver set is defined in `CONFIG$cell_types` (as run: MSN, Astrocyte, OPC, Oligo, Microglia, Vascular). |

For exact reproduction of the paper runs, the internal as-run column names are
also accepted and take effect when the generalized ones are absent:
`Final_Zone_Assignments_No_Smooth` (≙ `Spatial_Domain`) and `Class_name` +
`AK_celltype` (collapsed by the pipeline into the six analysis types).

The interaction database is **not** bundled: `CellChatDB.human.rda` is
downloaded automatically from the CellChat GitHub repository on first run and
cached under `00_db/` (an internet connection is needed once; subsequent runs
use the cache).

## Notes for reuse

- Scripts contain absolute paths to the original cluster filesystem
  (`/broad/macosko/...`); adjust the path constants at the top of each script
  (`FIG6`, `LRDIR`, `OUT_DIR`, and `RAW_DIR`/`OUT_BASE` in
  `build_filtered_donor_objs.R`) to your local layout. Paths were deliberately
  left as-run for provenance.
- To run on your own data, provide donor objects with `x_um2`/`y_um2`,
  `Spatial_Domain`, and `Cell_types` metadata — no other renaming needed.
- Scripts are as-run except for a small input adapter (added for this release)
  that accepts the generalized `Spatial_Domain`/`Cell_types` columns; the
  scoring, statistics, and null machinery are untouched.
- The permutation null shuffles zone labels **within donor**, so donor identity
  and per-donor density structure are preserved under the null.
