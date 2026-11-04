# BB72 gauging q(r) for a single logical X, Z, Y.
# q(r) = P(decoded outcome != eigenvalue) on +1 eigenstates, plain Relay-BP; state independence on -1
# eigenstates. --small-graphs: X and Z with the smaller 6-edge graphs (output bb72_qr_E6.csv).
# Chunks are cached in results/taskA_julia/bb72_qr_chunks*/ (same .npz files as the reference), so a run resumes.
# Run: julia --project=. scripts/taskA_bb72_calibration.jl [--quick] [--small-graphs] [--workers 12]
# Output: results/taskA_julia/bb72_qr.csv
const REPO = dirname(@__DIR__)
get!(ENV, "JULIA_PYTHONCALL_EXE", joinpath(REPO, ".venv", "bin", "python"))
get!(ENV, "JULIA_CONDAPKG_BACKEND", "Null")
using Distributed, FTShadows

const OUT = joinpath(REPO, "results", "taskA_julia")
const CHUNK = 5000

argval(flag, default) = (i = findfirst(==(flag), ARGS); i === nothing ? default : ARGS[i+1])
quick, small = "--quick" in ARGS, "--small-graphs" in ARGS
target, mx = quick ? (30, 3e4) : (200, 6e4)
t = time()
ps = small ? "XZ" : "XZY"
cdir = joinpath(OUT, small ? "bb72_qr_chunks_E6" : "bb72_qr_chunks")
mkpath(cdir)
cfgs = FTShadows.bb72_qr_configs(ps)  # the schedule the seed-uniqueness test checks
wp = addprocs(parse(Int, argval("--workers", "12")); exeflags = "--project=$(Base.active_project())")
rows = try
    @everywhere wp using FTShadows
    FTShadows.bb72_qr_adaptive(cfgs, target, mx; small, chunk_dir = cdir, chunk = CHUNK,
                               pmap = (f, xs) -> pmap(f, WorkerPool(wp), xs))
finally
    rmprocs(wp)
end
name = small ? "bb72_qr_E6" : "bb72_qr"
FTShadows.write_frame_csv(joinpath(OUT, "$name.csv"), rows)
for r in rows
    println(join((r[c] for c in ("P", "r", "state", "shots", "errors", "q", "ci_lo", "ci_hi", "cost_extra_qubits")), "  "))
end
write(joinpath(OUT, replace(name, "qr" => "calibration") * "_meta.json"),
      FTShadows.pyjson(FTShadows.OrderedDict("runtime_s" => time() - t, "quick" => quick)))
