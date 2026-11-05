# calibration.csv: q per setting type, code and noise model with 95% Wilson CIs (sections dataset_settings and scaling).
# Run: julia --project=. scripts/make_calibration.jl [OUT]   (default handover/data_julia/calibration.csv)
const REPO = dirname(@__DIR__)
using FTShadows

out = isempty(ARGS) ? joinpath(REPO, "handover", "data_julia", "calibration.csv") : ARGS[1]
abspath(out) == abspath(joinpath(REPO, "handover", "calibration.csv")) && error("refusing to overwrite handover/calibration.csv")
mkpath(dirname(out))
rows = FTShadows.calibration_rows(REPO; input = FTShadows.upstream)
FTShadows.write_csv_table(out, FTShadows.pandas_table(rows, FTShadows.CALIBRATION_COLUMNS))
println(length(rows), " rows -> ", out)
