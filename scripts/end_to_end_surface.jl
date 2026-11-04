# End-to-end check: surface-code Bell / GHZ
# states, random native settings, reference estimator vs exact values. Ensembles: "local" (uniform
# per-patch X/Z readouts) and "native" (1/2 local, 1/4 ZZ-chain, 1/4 XX-chain surgery).
# Run: julia --project=. scripts/end_to_end_surface.jl
# Output: results/e2e_julia/e2e_surface.csv, results/e2e_julia/e2e_surface.json
const REPO = dirname(@__DIR__)
get!(ENV, "JULIA_PYTHONCALL_EXE", joinpath(REPO, ".venv", "bin", "python"))
get!(ENV, "JULIA_CONDAPKG_BACKEND", "Null")
using FTShadows
using FTShadows: OrderedDict

out = joinpath(REPO, "results", "e2e_julia")
mkpath(out)
t = time()
rows = OrderedDict{String,Any}[]
for st in keys(FTShadows.E2E_STATES), en in ("local", "native"), p in (0.0, 1e-3)
    append!(rows, FTShadows.e2e_run(st, en, p, 20000, 1 + length(rows)))
end
FTShadows.pandas_csv(joinpath(out, "e2e_surface.csv"), rows)
summ = FTShadows.e2e_summary(rows)
res = OrderedDict("summary" => summ, "d" => FTShadows.E2E_D, "rounds" => FTShadows.E2E_R, "runtime_s" => time() - t)
write(joinpath(out, "e2e_surface.json"), FTShadows.pyjson(res; indent = 1))
foreach(println, summ)
foreach(println, filter(r -> r["lam"] > 0, rows))
