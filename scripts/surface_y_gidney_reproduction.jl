# Gidney's Y / X / Z memory numbers (arXiv:2302.07395) with his circuits and plain matching; X / Z give the decoder factor for Y.
# Run: julia --project=. scripts/surface_y_gidney_reproduction.jl [--workers 8]
# Output: results/taskA_julia/surface_y_gidney_reproduction.json
const REPO = dirname(@__DIR__)
get!(ENV, "JULIA_PYTHONCALL_EXE", joinpath(REPO, ".venv", "bin", "python"))
get!(ENV, "JULIA_CONDAPKG_BACKEND", "Null")
using Distributed, FTShadows

argval(flag, default) = (i = findfirst(==(flag), ARGS); i === nothing ? default : ARGS[i+1])
out = joinpath(REPO, "results", "taskA_julia", "surface_y_gidney_reproduction.json")
mkpath(dirname(out))
nw = parse(Int, argval("--workers", "8"))
stats = FTShadows.gidney_reference_stats(joinpath(REPO, "results", "taskA", "gidney_reference", "stats.csv"))
wp = addprocs(nw; exeflags = "--project=$(Base.active_project())")
@everywhere wp using FTShadows
try
    rows = FTShadows.gidney_reproduction(stats; pmap = (f, xs) -> pmap(f, WorkerPool(wp), xs), workers = nw)
    res = FTShadows.OrderedDict("rows" => rows, "noise" => "SI1000", "p" => 1e-3,
                                "decoder" => "pymatching (MatchingGapDecoder, no correlations)")
    write(out, FTShadows.pyjson(res; indent = 1))
finally
    rmprocs(wp)
end
