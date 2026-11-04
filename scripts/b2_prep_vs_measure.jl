# Cost of preparing a logical state vs measuring it with our local Pauli settings.
# Unit: physical qubit x round.
# Run: julia --project=. scripts/b2_prep_vs_measure.jl [OUT]   (default results/b2_julia/prep_vs_measure.csv)
const REPO = dirname(@__DIR__)
get!(ENV, "JULIA_PYTHONCALL_EXE", joinpath(REPO, ".venv", "bin", "python"))
get!(ENV, "JULIA_CONDAPKG_BACKEND", "Null")
using FTShadows, Printf

out = isempty(ARGS) ? joinpath(REPO, "results", "b2_julia", "prep_vs_measure.csv") : ARGS[1]
mkpath(dirname(out))
rows = FTShadows.prep_vs_measure_rows(FTShadows.upstream(REPO, "results", "taskA", "local_setting_costs.csv"),
                                      joinpath(REPO, "scripts", "data", "b2_cost_sources.json"))
FTShadows.write_prep_vs_measure(out, rows)
# a short summary table on stdout (not in pandas formatting)
fmt(x) = x === nothing ? "NaN" : x isa Integer ? string(x) : @sprintf("%.6g", x)
for (i, r) in enumerate(rows)
    println(i - 1, "  ", rpad(r["task"], 96), rpad(r["architecture"], 30), lpad(r["n_logical"], 4), "  ",
            lpad(fmt(r["prep_volume"]), 12), lpad(fmt(r["measure_volume"]), 12), lpad(fmt(r["ratio_measure_over_total"]), 13))
end
println(length(rows), " rows written to ", out)
