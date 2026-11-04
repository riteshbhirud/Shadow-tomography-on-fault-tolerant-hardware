# Deformed-code and circuit distances of BB72 X, Z, Y gauging, with all spectator logicals.
# "Y_minimal" is the rejected h = 0.6 graph, kept to
# reproduce its circuit distance 4.
# Run: julia --project=. scripts/taskA_bb_distances.jl    Output: results/taskA_julia/bb72_distances.json
const REPO = dirname(@__DIR__)
get!(ENV, "JULIA_PYTHONCALL_EXE", joinpath(REPO, ".venv", "bin", "python"))
get!(ENV, "JULIA_CONDAPKG_BACKEND", "Null")
using FTShadows

const R = 6
const OUT = joinpath(REPO, "results", "taskA_julia")
mkpath(OUT)
# Sequential and resumable (each search needs several GB): finished entries are kept in
# bb72_distances_partial.json until all four are done.
part = joinpath(OUT, "bb72_distances_partial.json")
res = FTShadows.read_json_or_empty(part)
for name in ["Y_minimal", "X", "Z", "Y"]
    haskey(res, name) && continue
    _, row = FTShadows.bb_distance_entry(name; rounds = R)
    res[name] = row
    FTShadows.write_json(part, res)
    println(name, " ", FTShadows.pyjson(FTShadows.jsonable(row)))
    flush(stdout)
end
FTShadows.write_json(joinpath(OUT, "bb72_distances.json"), res)
rm(part)
println(FTShadows.pyjson(FTShadows.jsonable(res); indent = 1))
