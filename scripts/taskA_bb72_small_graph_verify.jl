# Full check of the 6-edge Xbar_0 / Zbar_0 gauging graphs proposed as defaults:
# deformed-code distance, circuit distance
# with spectators, Cheeger constant and weights, noiseless exactness. The graphs are the first
# distance-6 rows of results/taskA/bb72_small_graph_search_spectators.json. Resumable.
# Run: julia --project=. scripts/taskA_bb72_small_graph_verify.jl    Output: results/taskA_julia/bb72_small_graph_verify.json
const REPO = dirname(@__DIR__)
get!(ENV, "JULIA_PYTHONCALL_EXE", joinpath(REPO, ".venv", "bin", "python"))
get!(ENV, "JULIA_CONDAPKG_BACKEND", "Null")
using FTShadows, JSON
using FTShadows: OrderedDict, Circuits

input_file(f) = FTShadows.upstream(REPO, "results", "taskA", f)
const OUT = joinpath(REPO, "results", "taskA_julia")
const R = 6

function chosen()
    d = JSON.parsefile(input_file("bb72_small_graph_search_spectators.json"))
    out = OrderedDict{Char,Vector{Tuple{Int,Int}}}()
    for P in "XZ", blk in d[string(P)], row in blk["rows"]
        row["distance"] == 6 && !haskey(out, P) && (out[P] = [Tuple(Int.(e)) for e in row["edges"]])
    end
    return out
end

mkpath(OUT)
fout = joinpath(OUT, "bb72_small_graph_verify.json")
res = FTShadows.read_json_or_empty(fout)
for (P, edges) in chosen()
    haskey(res, string(P)) && haskey(res[string(P)], "circuit_distance_all_targets") && continue
    t = time()
    exp = Circuits.bb_gauging(P, R; extra_targets = true, edges)
    d_isd, d_stim, k = FTShadows.deformed_distances(exp)
    gi = exp.graph_info
    row = OrderedDict{String,Any}("edges" => [collect(e) for e in edges], "graph" => gi,
        "cheeger" => Float64(Circuits.cheeger(exp.graph)), "extra_qubits" => gi["V"] + gi["E"] + gi["flux"], "deformed_k" => k,
        "deformed_distance_isd_upper" => d_isd, "deformed_distance_stim" => d_stim,
        "noiseless" => FTShadows.noiseless_gauging_check(P, edges), "rounds" => R)
    println(P, " ", FTShadows.pyjson(row))
    row["circuit_distance_L_only"] = FTShadows.circuit_distance(Circuits.bb_gauging(P, R; edges))
    println(P, " L only ", row["circuit_distance_L_only"])
    row["circuit_distance_all_targets"] = FTShadows.circuit_distance(exp)
    row["n_targets"] = length(exp.targets)
    row["seconds"] = time() - t
    println(P, " ", FTShadows.pyjson(row))
    flush(stdout)
    res[string(P)] = row
    FTShadows.write_json(fout, res)
end
