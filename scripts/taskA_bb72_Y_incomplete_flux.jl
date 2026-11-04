# Test whether the earlier 33-qubit Ybar graph was an incomplete-flux construction.
# Keeps 2 flux checks of the adopted
# 20-edge and the minimal 18-edge Ybar graphs and searches circuit distance at r = 6.
# --spectators [--rounds r]: the 33-qubit reconstruction with all readable spectator logicals observed
# (default r = 3: the r = 6 search grew past 36 GB); --control keeps every flux check.
# Run: julia --project=. scripts/taskA_bb72_Y_incomplete_flux.jl [--spectators] [--control] [--rounds R]
# Output: results/taskA_julia/bb72_Y_incomplete_flux*.json
const REPO = dirname(@__DIR__)
get!(ENV, "JULIA_PYTHONCALL_EXE", joinpath(REPO, ".venv", "bin", "python"))
get!(ENV, "JULIA_CONDAPKG_BACKEND", "Null")
using FTShadows

const OUT = joinpath(REPO, "results", "taskA_julia")
const R = 6
mkpath(OUT)
if "--spectators" in ARGS
    i = findfirst(==("--rounds"), ARGS)
    control = "--control" in ARGS
    row = FTShadows.incomplete_flux_spectators(i === nothing ? 3 : parse(Int, ARGS[i+1]); control)
    println(FTShadows.pyjson(row))
    FTShadows.write_json(joinpath(OUT, "bb72_Y_incomplete_flux_spectators" * (control ? "_control" : "") * ".json"), row)
else
    FTShadows.write_json(joinpath(OUT, "bb72_Y_incomplete_flux.json"), FTShadows.incomplete_flux_rows(R))
end
