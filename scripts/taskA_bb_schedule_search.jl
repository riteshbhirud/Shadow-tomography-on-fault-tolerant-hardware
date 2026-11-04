# Search translation-invariant CNOT orders for BB72 syndrome extraction that keep circuit distance 6.
# The memory experiment it uses is
# FTShadows.Circuits.memory_experiment, also behind taskA_bb72_memory_sanity.jl.
# Run: julia --project=. scripts/taskA_bb_schedule_search.jl [WORKERS]    Output: results/taskA_julia/bb72_schedule_search.json
const REPO = dirname(@__DIR__)
get!(ENV, "JULIA_PYTHONCALL_EXE", joinpath(REPO, ".venv", "bin", "python"))
get!(ENV, "JULIA_CONDAPKG_BACKEND", "Null")
using Distributed, FTShadows

const OUT = joinpath(REPO, "results", "taskA_julia")
mkpath(OUT)
t = time()
jobs = [(k, p) for k in "XZ" for p in FTShadows.py_permutations(6)]
wp = addprocs(isempty(ARGS) ? 4 : parse(Int, ARGS[1]); exeflags = "--project=$(Base.active_project())")
res = FTShadows.OrderedDict("X" => FTShadows.OrderedDict{String,Int}(), "Z" => FTShadows.OrderedDict{String,Int}())
try
    @everywhere wp using FTShadows
    # results are kept in job order (the reference stores them in completion order)
    for (kind, perm, d) in pmap(j -> FTShadows.schedule_distance(j...), WorkerPool(wp), jobs; batch_size = 8)
        res[string(kind)][join(perm)] = d
    end
finally
    rmprocs(wp)
end
summary = FTShadows.schedule_summary(res)
FTShadows.write_json(joinpath(OUT, "bb72_schedule_search.json"),
                     FTShadows.OrderedDict("summary" => summary, "all" => res, "runtime_s" => time() - t))
println(FTShadows.pyjson(summary; indent = 1))
