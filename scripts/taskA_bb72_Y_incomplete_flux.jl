# Was the earlier 33-qubit Ybar graph an incomplete-flux construction? Circuit distance with only 2 flux checks kept.
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
    # the 33-qubit graph with all readable spectators observed; default r = 3 since the r = 6 search went past 36 GB.
    # --control keeps every flux check.
    i = findfirst(==("--rounds"), ARGS)
    control = "--control" in ARGS
    row = FTShadows.incomplete_flux_spectators(i === nothing ? 3 : parse(Int, ARGS[i+1]); control)
    println(FTShadows.pyjson(row))
    FTShadows.write_json(joinpath(OUT, "bb72_Y_incomplete_flux_spectators" * (control ? "_control" : "") * ".json"), row)
else
    FTShadows.write_json(joinpath(OUT, "bb72_Y_incomplete_flux.json"), FTShadows.incomplete_flux_rows(R))
end
