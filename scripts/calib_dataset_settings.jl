# Flip rates of the settings in the v1.0/v1.1 datasets, per string and per subset product (flips within a setting are correlated).
# Run: julia --project=. scripts/calib_dataset_settings.jl [--workers 4]
# Output: results/calibration_julia/dataset_settings_q.csv
const REPO = dirname(@__DIR__)
get!(ENV, "JULIA_PYTHONCALL_EXE", joinpath(REPO, ".venv", "bin", "python"))
get!(ENV, "JULIA_CONDAPKG_BACKEND", "Null")
using Distributed, FTShadows

argval(flag, default) = (i = findfirst(==(flag), ARGS); i === nothing ? default : ARGS[i+1])
out = joinpath(REPO, "results", "calibration_julia", "dataset_settings_q.csv")
mkpath(dirname(out))
wp = addprocs(parse(Int, argval("--workers", "4")); exeflags = "--project=$(Base.active_project())")
@everywhere wp using FTShadows
try
    rows = reduce(vcat, pmap(FTShadows.dataset_settings_job, WorkerPool(wp), FTShadows.dataset_settings_jobs()))
    FTShadows.pandas_csv(out, rows)
    for r in rows
        println(join((r[c] for c in ("code", "setting_id", "n_logical", "target", "kind", "shots", "flips", "q",
                                     "q_if_independent")), "  "))
    end
finally
    rmprocs(wp)
end
