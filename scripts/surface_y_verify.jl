# Checks of the surface-code in-place Y measurement: Gidney's circuits, determinism, noiseless exactness, distance per (d, b).
# Run: julia --project=. scripts/surface_y_verify.jl
# Output: results/taskA_julia/surface_y_verify.json
const REPO = dirname(@__DIR__)
get!(ENV, "JULIA_PYTHONCALL_EXE", joinpath(REPO, ".venv", "bin", "python"))
get!(ENV, "JULIA_CONDAPKG_BACKEND", "Null")
using FTShadows
using FTShadows: OrderedDict, pyjson

out = joinpath(REPO, "results", "taskA_julia", "surface_y_verify.json")
mkpath(dirname(out))
t = time()
res = OrderedDict{String,Any}(
    "reproduction_byte_identical" => FTShadows.y_reproduction(joinpath(REPO, "results", "taskA", "gidney_reference", "circuits")),
    "determinism" => FTShadows.y_determinism(),
    "noiseless" => FTShadows.y_noiseless(),
    "readout_crosscheck" => FTShadows.y_readout_crosscheck())
println(pyjson(res; indent = 1))
res["distances"] = FTShadows.y_distances()
res["runtime_s"] = round(Int, time() - t)
write(out, pyjson(res; indent = 1))
