# Noiseless exactness of every compiled BB72 local setting, 2 and 3 logical qubits.
# Random stabilizer states and joint
# eigenstates; every subset product must match the exact expectation.
# Run: julia --project=. scripts/taskA_bb72_local_noiseless.jl    Output: results/taskA_julia/bb72_local_noiseless.json
const REPO = dirname(@__DIR__)
get!(ENV, "JULIA_PYTHONCALL_EXE", joinpath(REPO, ".venv", "bin", "python"))
get!(ENV, "JULIA_CONDAPKG_BACKEND", "Null")
using FTShadows

const OUT = joinpath(REPO, "results", "taskA_julia")
mkpath(OUT)
t0 = time()
rows = []
for m in (2, 3), bases in FTShadows.product_strings("XYZ", m)
    push!(rows, FTShadows.local_noiseless_row(bases))
    println(FTShadows.pyjson(rows[end]), " ", round(Int, time() - t0))
    flush(stdout)
end
FTShadows.write_json(joinpath(OUT, "bb72_local_noiseless.json"),
    FTShadows.OrderedDict("rows" => rows, "all_ok" => all(r["mismatches"] == 0 for r in rows), "runtime_s" => round(Int, time() - t0)))
