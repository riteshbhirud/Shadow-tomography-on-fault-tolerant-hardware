# La-cross [[65,9,4]] feasibility scan: code (checks, distance, folds), circuit (circuit distance), sim (BP-OSD memory), summary.
# Run: julia --project=. scripts/b3_lacross.jl {code|circuit|sim|summary} [options] [--out DIR]
# Output: results/b3_julia/ unless --out is given.
const REPO = dirname(@__DIR__)
get!(ENV, "JULIA_PYTHONCALL_EXE", joinpath(REPO, ".venv", "bin", "python"))
get!(ENV, "JULIA_CONDAPKG_BACKEND", "Null")
using Distributed, FTShadows

const USAGE = "usage: b3_lacross.jl {code|circuit|sim|summary} [options] [--out DIR]"
isempty(ARGS) && (println(USAGE); exit(2))
cmd = ARGS[1]
argval(flag, default) = (i = findfirst(==(flag), ARGS); i === nothing ? default : ARGS[i+1])
const OUT = argval("--out", joinpath(REPO, "results", "b3_julia"))
mkpath(OUT)
writejson(name, x) = write(joinpath(OUT, name), FTShadows.pyjson(x; indent = 1))

if cmd == "code"
    rec, arrays = FTShadows.lacross_code_record()
    np = FTShadows.pymodule(:numpy)
    np.savez(joinpath(OUT, "lacross_65_9_4_checks.npz"); (Symbol(k) => FTShadows.to_numpy_u8(v) for (k, v) in arrays)...)
    writejson("code.json", rec)
    println(FTShadows.pyjson(filter(kv -> !(kv[1] in ("seed_H", "seed_G")), rec); indent = 1))
elseif cmd == "circuit"
    rec = FTShadows.lacross_circuit_record(; bases = argval("--bases", "Z"), max_dets = argval("--max_dets", "3,4"),
                                           max_degree = parse(Int, argval("--max_degree", "6")),
                                           save = "--save" in ARGS, out = OUT)
    writejson("circuit.json", rec)
elseif cmd == "sim"
    p = parse(Float64, argval("--p", "1e-3"))
    basis, noise = argval("--basis", "Z"), argval("--noise", "uniform")
    batch, batches = parse(Int, argval("--batch", "2000")), parse(Int, argval("--batches", "100"))
    workers_ = min(3, parse(Int, argval("--workers", "3")))
    seed = parse(Int, argval("--seed", "1"))
    dec = Dict("max_iter" => parse(Int, argval("--max_iter", "100")), "ms" => parse(Float64, argval("--ms", "0.625")),
               "osd_order" => parse(Int, argval("--osd_order", "7")))
    target_fails = parse(Int, argval("--target_fails", "200"))
    max_seconds = parse(Float64, argval("--max_seconds", "1200"))
    pf = FTShadows.pyfloat
    key = "$(noise)_$(basis)_p$(pf(p))_it$(dec["max_iter"])_ms$(pf(dec["ms"]))_osd$(dec["osd_order"])"
    t0 = time()
    wp = addprocs(workers_; exeflags = "--project=$(Base.active_project())")
    @everywhere wp using FTShadows
    # batches are handed out as workers free up and counted in the order they finish, until the time or
    # failure budget runs out, then the rest are dropped
    function run_batches(wp)
        fails = shots = 0
        next, stop = Ref(0), Ref(false)
        @sync for w in wp
            @async while !stop[] && next[] < batches
                b = next[]
                next[] += 1
                f, s = remotecall_fetch(FTShadows.lacross_sim_batch, w, key, p, basis, noise, batch, 1000 * seed + b, dec)
                stop[] && break
                fails += f
                shots += s
                (time() - t0 > max_seconds || fails >= target_fails) && (stop[] = true)
            end
        end
        return fails, shots
    end
    fails, shots = try
        run_batches(wp)
    finally
        rmprocs(wp)
    end
    rec = FTShadows.lacross_sim_record(key, p, basis, noise, fails, shots, dec, seed,
                                       FTShadows.py_round(time() - t0, 1), workers_)
    println(FTShadows.pyjson(rec))
    open(io -> println(io, FTShadows.pyjson(rec)), joinpath(OUT, "sim_results.jsonl"), "a")
elseif cmd == "summary"
    writejson("lacross_summary.json", FTShadows.lacross_summary(OUT))
    println("wrote ", joinpath(OUT, "lacross_summary.json"))
else
    println(USAGE)
    exit(2)
end
