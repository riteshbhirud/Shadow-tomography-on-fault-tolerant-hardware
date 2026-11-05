# Circuit distance of the two-Xbar simultaneous gauging circuit at r = 6, with spectators.
# Run: julia --project=. scripts/taskA_bb72_multi_distance.jl
# Output: results/taskA_julia/bb72_multi_distance.json
const REPO = dirname(@__DIR__)
get!(ENV, "JULIA_PYTHONCALL_EXE", joinpath(REPO, ".venv", "bin", "python"))
get!(ENV, "JULIA_CONDAPKG_BACKEND", "Null")
using FTShadows
using FTShadows: OrderedDict

const OUT = joinpath(REPO, "results", "taskA_julia")
const R = 6
t = time()
code, Ls = FTShadows.multi_pair()
spec = FTShadows.multi_spectators(code, Ls)
res = OrderedDict{String,Any}()
for (label, extra) in [("L_only", []), ("L_plus_spectators", spec)]
    exp = FTShadows.multi_experiment(code, Ls, R; extra_ops = extra, name = "bb72_multi_XX_distance")
    d = FTShadows.circuit_distance(exp)
    res[label] = OrderedDict{String,Any}("n_targets" => length(exp.targets), "circuit_distance" => d, "graph" => exp.graph_info,
                                         "seconds" => time() - t)
    println(label, " ", FTShadows.pyjson(res[label]))
    flush(stdout)
end
res["n_spectators"] = length(spec)
mkpath(OUT)
FTShadows.write_json(joinpath(OUT, "bb72_multi_distance.json"), res)
