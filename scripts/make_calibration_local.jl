# calibration_local.csv: the flip rates behind the noise-expected values of the Nov-20 datasets, joined from results/taskA.
# Run: julia --project=. scripts/make_calibration_local.jl [OUT]   (default handover/data_julia/calibration_local.csv)
const REPO = dirname(@__DIR__)
using FTShadows

out = isempty(ARGS) ? joinpath(REPO, "handover", "data_julia", "calibration_local.csv") : ARGS[1]
abspath(out) == abspath(joinpath(REPO, "handover", "calibration_local.csv")) && error("refusing to overwrite handover/calibration_local.csv")
mkpath(dirname(out))
res(f) = FTShadows.upstream(REPO, "results", "taskA", f)
rows = FTShadows.calibration_local_join(res("surface_local_q.csv"), res("bb72_local_q.csv"))
FTShadows.write_csv_table(out, FTShadows.pandas_table(rows, FTShadows.CALIBRATION_LOCAL_COLUMNS))
println(length(rows), " rows")
