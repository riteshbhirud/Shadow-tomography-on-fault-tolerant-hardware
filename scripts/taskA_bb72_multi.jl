# Simultaneous gauging of two commuting logicals {Xbar_a, Xbar_b} on BB72: noiseless exactness, deformed distance, flip rates.
# Run: julia --project=. scripts/taskA_bb72_multi.jl
# Output: results/taskA_julia/bb72_multi.json
const REPO = dirname(@__DIR__)
get!(ENV, "JULIA_PYTHONCALL_EXE", joinpath(REPO, ".venv", "bin", "python"))
get!(ENV, "JULIA_CONDAPKG_BACKEND", "Null")
using FTShadows
using FTShadows: OrderedDict, bitmean

const OUT = joinpath(REPO, "results", "taskA_julia")
const SHOTS = 30000
t = time()
code, Ls = FTShadows.multi_pair()
out = OrderedDict{String,Any}("supports" => [sort(collect(keys(L))) for L in Ls])
exp_for(signs, rounds) = FTShadows.multi_experiment(code, Ls, rounds; state = collect(zip(Ls, signs)), name = "bb72_multi_XX")

ok = true
for signs in [(0, 0), (0, 1), (1, 0), (1, 1)]
    res = run_oracle(exp_for(signs, 2), 0.0, 100; seed = 1)
    global ok &= all(res.outcomes .== permutedims(collect(signs)))
end
out["noiseless_exact_all_signs"] = ok
e = exp_for((0, 0), 3)
d, k = FTShadows.Circuits.code_distance_stim(e.n_data, FTShadows.epoch_checks(e))
out["graph"], out["deformed_distance"], out["deformed_k"] = e.graph_info, d, k
res = run_oracle(e, 1e-3, SHOTS; decoder = "relay", seed = 7, with_gap = false)  # q only; the gap is measured separately
err = Int.(res.outcomes)  # expected all zeros
q1, q2 = bitmean(err[:, 1]), bitmean(err[:, 2])
merge!(out, OrderedDict{String,Any}("shots" => size(err, 1), "cost" => res.cost, "q" => [q1, q2],
    "q_ci" => [collect(wilson(count(!iszero, err[:, j]), size(err, 1))) for j in 1:2],
    "q_any" => bitmean(err[:, 1] .| err[:, 2]), "q_both" => bitmean(err[:, 1] .& err[:, 2]), "q_both_if_independent" => q1 * q2,
    "flip_correlation" => any(!iszero, err) ? FTShadows.np_corrcoef(err[:, 1], err[:, 2]) : nothing,
    "runtime_s" => time() - t))
mkpath(OUT)
FTShadows.write_json(joinpath(OUT, "bb72_multi.json"), out)
println(FTShadows.pyjson(FTShadows.jsonable(out); indent = 1))
