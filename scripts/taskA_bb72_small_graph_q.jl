# Flip rate q of the 6-edge Xbar_0 / Zbar_0 gauging graphs that passed the spectator screen.
# Same r, p, noise, decoder and state as the
# 9-edge defaults in results/taskA/bb72_qr.csv. The screen and the defaults are read from results/taskA.
# Run: julia --project=. scripts/taskA_bb72_small_graph_q.jl [--workers 6] [--shots 30000]
# Output: results/taskA_julia/bb72_small_graph_q.csv
const REPO = dirname(@__DIR__)
get!(ENV, "JULIA_PYTHONCALL_EXE", joinpath(REPO, ".venv", "bin", "python"))
get!(ENV, "JULIA_CONDAPKG_BACKEND", "Null")
using Distributed, FTShadows, JSON

input_file(f) = FTShadows.upstream(REPO, "results", "taskA", f)
const OUT = joinpath(REPO, "results", "taskA_julia")
const P_PHYS, R, CHUNK = 1e-3, 3, 2500

argval(flag, default) = (i = findfirst(==(flag), ARGS); i === nothing ? default : ARGS[i+1])
shots = parse(Int, argval("--shots", "30000"))
cdir = joinpath(OUT, "bb72_small_graph_q_chunks")
mkpath(cdir)
cands = FTShadows.small_graph_q_candidates(JSON.parsefile(input_file("bb72_small_graph_search_spectators.json")))
jobs = [(i - 1, P, es, c, cdir, CHUNK) for (i, (P, es)) in enumerate(cands) for c in 0:shots÷CHUNK-1]  # 0-based ids name the chunk files
errs = zeros(Int, length(cands))
wp = addprocs(parse(Int, argval("--workers", "6")); exeflags = "--project=$(Base.active_project())")
try
    @everywhere wp using FTShadows
    for (i, e) in pmap(FTShadows.small_graph_q_chunk, WorkerPool(wp), jobs)
        errs[i+1] += e
    end
finally
    rmprocs(wp)
end
ref = FTShadows.csv_rows(input_file("bb72_qr.csv"))
rows = []
for (i, (P, es)) in enumerate(cands)
    n = (shots ÷ CHUNK) * CHUNK
    lo, hi = wilson(errs[i], n)
    exp = FTShadows.small_graph_q_experiment(P[1], es; rounds = R)
    d = first(x for x in ref if x["P"] == P && x["r"] == R && x["state"] == "+1")
    push!(rows, FTShadows.OrderedDict{String,Any}("P" => P, "edges" => FTShadows.pyjson([collect(e) for e in es]),
        "E" => length(es), "extra_qubits" => 2 * length(es) + 1, "max_check_weight" => exp.graph_info["max_check_weight"],
        "r" => R, "p" => P_PHYS, "shots" => n, "errors" => errs[i], "q" => errs[i] / n, "ci_lo" => lo, "ci_hi" => hi,
        "default_q" => d["q"], "default_ci_lo" => d["ci_lo"], "default_ci_hi" => d["ci_hi"],
        "default_shots" => d["shots"], "default_extra_qubits" => d["cost_extra_qubits"]))
    println(FTShadows.pyjson(rows[end]))
    flush(stdout)
end
FTShadows.write_frame_csv(joinpath(OUT, "bb72_small_graph_q.csv"), rows)
