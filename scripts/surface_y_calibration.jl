# Flip rate q(b) of the surface-code in-place Y measurement as a single-patch setting:
# d = 3, 5, 7 and b = 0 .. d÷2 + 2; at the
# default b also the -1 eigenstate and a Z eigenstate (<Y> = 0).
# Run: julia --project=. scripts/surface_y_calibration.jl [--workers 8]
# Output: results/taskA_julia/surface_y_qr.csv
const REPO = dirname(@__DIR__)
get!(ENV, "JULIA_PYTHONCALL_EXE", joinpath(REPO, ".venv", "bin", "python"))
get!(ENV, "JULIA_CONDAPKG_BACKEND", "Null")
using Distributed, FTShadows

argval(flag, default) = (i = findfirst(==(flag), ARGS); i === nothing ? default : ARGS[i+1])
out = joinpath(REPO, "results", "taskA_julia", "surface_y_qr.csv")
mkpath(dirname(out))
wp = addprocs(parse(Int, argval("--workers", "8")); exeflags = "--project=$(Base.active_project())")
@everywhere wp using FTShadows
try
    rows = FTShadows.surface_y_calibration(; pmap = (f, xs) -> pmap(f, WorkerPool(wp), xs))
    FTShadows.pandas_csv(out, rows)
finally
    rmprocs(wp)
end
