# Compiled cost of every local Pauli setting.
# Run: julia --project=. scripts/local_setting_costs.jl [OUT]   (default results/taskA_julia/local_setting_costs.csv)
const REPO = dirname(@__DIR__)
get!(ENV, "JULIA_PYTHONCALL_EXE", joinpath(REPO, ".venv", "bin", "python"))
get!(ENV, "JULIA_CONDAPKG_BACKEND", "Null")
using FTShadows

out = isempty(ARGS) ? joinpath(REPO, "results", "taskA_julia", "local_setting_costs.csv") : ARGS[1]
mkpath(dirname(out))
res(f) = FTShadows.upstream(REPO, "results", "taskA", f)
rows = FTShadows.local_setting_costs(res("surface_local_q.csv"), res("bb72_local_q.csv"))
FTShadows.write_setting_costs(out, rows)
println(length(rows), " rows written to ", out)
