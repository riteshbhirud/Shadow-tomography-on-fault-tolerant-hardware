# Decoder sanity for BB72: Z- and X-basis memory with our syndrome schedule, r = 6, p = 1e-3, Relay-BP.
# Compared with Bravyi et al.'s ~7e-5
# per cycle only to within a small factor (different schedule, noise and decoder).
# Run: julia --project=. scripts/taskA_bb72_memory_sanity.jl [BASES] [--workers 12]
# Output: results/taskA_julia/bb72_memory_sanity.json
const REPO = dirname(@__DIR__)
get!(ENV, "JULIA_PYTHONCALL_EXE", joinpath(REPO, ".venv", "bin", "python"))
get!(ENV, "JULIA_CONDAPKG_BACKEND", "Null")
using Distributed, FTShadows

const OUT = joinpath(REPO, "results", "taskA_julia")
const R, P, SEEDS, SHOTS = 6, 1e-3, 24, 20000

argval(flag, default) = (i = findfirst(==(flag), ARGS); i === nothing ? default : ARGS[i+1])
pos = filter(a -> !startswith(a, "--") && (i = findfirst(==(a), ARGS); i == 1 || ARGS[i-1] != "--workers"), ARGS)
bases = isempty(pos) ? "ZX" : pos[1]
t = time()
mkpath(OUT)
fout = joinpath(OUT, "bb72_memory_sanity.json")
out = isfile(fout) && bases != "ZX" ? FTShadows.read_json_or_empty(fout) : FTShadows.OrderedDict{String,Any}()
wp = addprocs(parse(Int, argval("--workers", "12")); exeflags = "--project=$(Base.active_project())")
try
    @everywhere wp using FTShadows
    for noise in ["sd6", "stim_generated"], bk in bases
        jobs = [(seed_for("bb72_memory", noise, string(bk), s), SHOTS, noise, string(bk)) for s in 0:SEEDS-1]
        res = pmap(FTShadows.bb72_memory_failures, WorkerPool(wp), jobs)
        k, n = sum(first, res), sum(last, res)
        lo, hi = wilson(k, n)
        out["$(noise)_$(bk)memory"] = FTShadows.OrderedDict{String,Any}("shots" => n, "failures" => k, "P_L" => k / n,
            "ci" => [lo, hi], "per_round" => 1 - FTShadows.cpow(1 - k / n, 1 / R), "rounds" => R, "p" => P)
        println(noise, " ", bk, " ", FTShadows.pyjson(out["$(noise)_$(bk)memory"]))
        flush(stdout)
    end
finally
    rmprocs(wp)
end
out["reference"] = "Bravyi et al. 2308.07915: BB72 ~7e-5 per cycle at p=1e-3 (X and Z failures)"
out["runtime_s_$bases"] = time() - t
FTShadows.write_json(fout, out)
