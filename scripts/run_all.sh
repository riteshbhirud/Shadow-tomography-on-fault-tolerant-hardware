#!/usr/bin/env bash
# Reproduce every Task A, end-to-end, dataset and Task B result. Long runs (~h) resume when rerun.
# Usage: bash scripts/run_all.sh [--fast]   (--fast: the test suite only, about 25 min)
# Output: results/*_julia/ and handover/data_julia/, so the stored results are not overwritten.
set -euo pipefail
cd "$(dirname "$0")/.."
export JULIA_PYTHONCALL_EXE="${JULIA_PYTHONCALL_EXE:-$PWD/.venv/bin/python}"
export JULIA_CONDAPKG_BACKEND="${JULIA_CONDAPKG_BACKEND:-Null}"
JL="julia --project=."
mkdir -p results
$JL -e 'using Pkg; Pkg.test()' 2>&1 | tee results/test_output_julia.txt
if [[ "${1:-}" == "--fast" ]]; then exit 0; fi

# noise-model calibration against Bravyi et al. 2308.07915
$JL scripts/calib_bravyi_reproduction.jl          # ~h  results/calibration_julia/bravyi_reproduction.json
$JL scripts/calib_bravyi_decoder_sensitivity.jl   #     results/calibration_julia/bravyi_decoder_sensitivity.json

# surface code
$JL scripts/surface_calibration.jl                # ~h  results/taskA_julia/surface_qr.csv, surface_state_independence.csv

# BB72 (stim distance searches need up to ~10 GB each, so run them one at a time)
$JL scripts/taskA_bb72_memory_sanity.jl           #     bb72_memory_sanity.json
$JL scripts/taskA_bb_distances.jl                 #     bb72_distances.json (X, Z, adopted Y, minimal Y)
$JL scripts/taskA_bb_gauss_order_screen.jl Y '{"target_cheeger": 1.5}' 8   # bb72_gauss_order_screen_Y.json
$JL scripts/taskA_bb72_calibration.jl             # ~h  bb72_qr.csv
$JL scripts/taskA_bb72_gap_decoder_comparison.jl  #     bb72_gap_decoder_comparison.json
$JL scripts/taskA_bb72_postselection.jl           # ~h  bb72_postselection.csv, bb72_post_*.npz
$JL scripts/taskA_bb72_multi.jl                   #     bb72_multi.json
$JL scripts/taskA_bb72_multi_distance.jl          # ~h  bb72_multi_distance.json
$JL scripts/taskA_bb72_Y_incomplete_flux.jl       #     bb72_Y_incomplete_flux.json
$JL scripts/taskA_bb72_Y_incomplete_flux.jl --spectators --rounds 3
$JL scripts/taskA_bb72_Y_incomplete_flux.jl --spectators --control --rounds 3
$JL scripts/taskA_bb72_small_graph_search.jl      #     bb72_small_graph_search.json (measured logical only)
$JL scripts/taskA_bb72_small_graph_search.jl --spectator-screen --workers 2   # ~h
$JL scripts/taskA_bb72_small_graph_q.jl           #     bb72_small_graph_q.csv
$JL scripts/taskA_bb72_small_graph_verify.jl      #     bb72_small_graph_verify.json
$JL scripts/taskA_bb72_calibration.jl --small-graphs   # ~h  bb72_qr_E6.csv

# surface-code Y (Gidney), BB72 local settings, Bell/GHZ datasets
$JL scripts/surface_y_gidney_reproduction.jl      #     surface_y_gidney_reproduction.json
$JL scripts/surface_y_verify.jl                   #     surface_y_verify.json (distances, noiseless)
$JL scripts/surface_y_calibration.jl              #     surface_y_qr.csv
$JL scripts/surface_local_calibration.jl          #     surface_local_q.csv, surface_local_independence.csv
$JL scripts/taskA_bb72_local_noiseless.jl         #     bb72_local_noiseless.json
$JL scripts/taskA_bb72_local_calibration.jl --ghz # ~h  bb72_local_q.csv (resumable)
$JL scripts/taskA_bb72_local_distance.jl --only 2 # ~h  bb72_local_distance.json (Bell, with spectators)
$JL scripts/taskA_bb72_local_distance.jl --only 3 --exclude 2 --out results/taskA_julia/bb72_local_distance_ghz.json
$JL scripts/local_setting_costs.jl                #     local_setting_costs.csv
$JL scripts/b2_prep_vs_measure.jl                 #     results/b2_julia/prep_vs_measure.csv
for s in code circuit sim summary; do $JL scripts/b3_lacross.jl $s; done   # results/b3_julia/ (La-cross)

# end-to-end and plots
$JL scripts/end_to_end_surface.jl                 #     results/e2e_julia/
julia --project=scripts/plots scripts/plots/plot_taskA.jl --in results/taskA_julia   # PNGs

# datasets (handover/data_julia/)
$JL scripts/calib_dataset_settings.jl             #     results/calibration_julia/dataset_settings_q.csv
$JL scripts/make_datasets.jl                      #     surface and BB72 shot records
$JL scripts/make_calibration.jl                   #     calibration.csv
$JL scripts/make_bases.jl                         #     bases/
$JL scripts/make_calibration_local.jl             #     calibration_local.csv
$JL scripts/make_local_datasets.jl                # ~h  nov20/
$JL scripts/nov20_postselection.jl                #     results/postselection_julia/

# Task B (needs Oscar: the small scripts/oscar environment), then the native ensembles built from it
julia --project=scripts/oscar scripts/taskb_run_all.jl   # ~h  results/taskB_julia/
$JL scripts/make_native_ensembles.jl              #     native_ensembles/
