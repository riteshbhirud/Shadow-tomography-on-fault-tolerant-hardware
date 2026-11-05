# Screen seeded Gauss-check gate orders for BB72 gauging: refuted if stim finds distance < 6, null if the search times out.
# Run: julia --project=. scripts/taskA_bb_gauss_order_screen.jl P [KWARGS_JSON] [N_SEEDS]
# Output: results/taskA_julia/bb72_gauss_order_screen_{P}.json
const REPO = dirname(@__DIR__)
get!(ENV, "JULIA_PYTHONCALL_EXE", joinpath(REPO, ".venv", "bin", "python"))
get!(ENV, "JULIA_CONDAPKG_BACKEND", "Null")
using Distributed, FTShadows, JSON

const OUT = joinpath(REPO, "results", "taskA_julia")
const BATCH, LIMIT = 8, 90

# One worker process per seed, all started at once; after the time limit the ones still searching are killed.
function run_batch(P, seeds, kw, limit)
    wp = addprocs(length(seeds); exeflags = "--project=$(Base.active_project())")
    res = FTShadows.OrderedDict{Int,Any}()
    try
        Distributed.remotecall_eval(Main, wp, :(using FTShadows))
        # compile the builder first so the time limit covers only the search
        Distributed.remotecall_eval(Main, wp, :(FTShadows.Circuits.build(FTShadows.Circuits.bb_gauging('X', 1), 1e-3); nothing))
        # each result is fetched by a local task: asking a worker isready() would block while it sits in stim
        tasks = [@async remotecall_fetch(FTShadows.gauss_order_distance, w, P, sd, kw) for (sd, w) in zip(seeds, wp)]
        t0 = time()
        for (sd, tk) in zip(seeds, tasks)
            timedwait(() -> istaskdone(tk), max(0.0, limit - (time() - t0)); pollint = 0.2)
            res[sd] = istaskdone(tk) ? fetch(tk) : nothing
        end
    finally
        for w in wp
            p = Distributed.worker_from_id(w).config.process
            p === nothing || kill(p)
        end
        try
            rmprocs(wp; waitfor = 5)
        catch err
            @warn "removing killed workers" err
        end
    end
    return res
end

P = ARGS[1]
kw = length(ARGS) > 1 ? JSON.parse(ARGS[2]) : FTShadows.OrderedDict{String,Any}()
seeds = collect(0:(length(ARGS) > 2 ? parse(Int, ARGS[3]) : 24)-1)
t = time()
res = FTShadows.OrderedDict{Int,Any}()
for i in 1:BATCH:length(seeds)
    merge!(res, run_batch(P[1], seeds[i:min(i + BATCH - 1, end)], Dict{String,Any}(kw), LIMIT))
    println(res)
    flush(stdout)
end
mkpath(OUT)
out = FTShadows.OrderedDict("P" => P, "graph_kw" => kw, "results" => FTShadows.OrderedDict(string(k) => v for (k, v) in res),
                            "runtime_s" => time() - t)
FTShadows.write_json(joinpath(OUT, "bb72_gauss_order_screen_$P.json"), out)
println(FTShadows.pyjson(FTShadows.jsonable(out)))
