# Verification of the surface-code in-place Y measurement as an oracle setting:
# byte reproduction of Gidney's circuits, detector
# determinism, noiseless exactness on random states, a cross-check of the X / Z fragments against the
# builder's readout, and the circuit distance per (d, b).
# Run: julia --project=. scripts/surface_y_verify.jl    Output: results/taskA_julia/surface_y_verify.json
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
