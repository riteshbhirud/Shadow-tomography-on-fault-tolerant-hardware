# Single-patch flip rates of the surface local settings, as rows of calibration_local.csv.
# Run: julia --project=. scripts/surface_local_calibration.jl [--d 3] [--workers 6] [--out FILE]
# Output: handover/data_julia/calibration_local_surface_d3.csv unless --out is given.
const REPO = dirname(@__DIR__)
get!(ENV, "JULIA_PYTHONCALL_EXE", joinpath(REPO, ".venv", "bin", "python"))
get!(ENV, "JULIA_CONDAPKG_BACKEND", "Null")
using Distributed, FTShadows

argval(flag, default) = (i = findfirst(==(flag), ARGS); i === nothing ? default : ARGS[i+1])
ds_ = parse.(Int, split(argval("--d", "3"), ','))
out = argval("--out", joinpath(REPO, "handover", "data_julia", "calibration_local_surface_d$(join(ds_, '_')).csv"))
mkpath(dirname(out))
wp = addprocs(parse(Int, argval("--workers", "6")); exeflags = "--project=$(Base.active_project())")
@everywhere wp using FTShadows
try
    cfgs = [(d, string(b), r) for d in ds_ for b in "XZY" for r in (0, 3, 6)]
    st = FTShadows.surface_calibration(cfgs; pmap = (f, xs) -> pmap(f, WorkerPool(wp), xs))
    rows = FTShadows.calibration_local_rows(st)
    cols = collect(keys(rows[1]))
    FTShadows.write_csv_table(out, FTShadows.CSVTable(cols, [Any[r[c] for c in cols] for r in rows], "\n"))
    println(length(rows), " rows written to ", out)
finally
    rmprocs(wp)
end
