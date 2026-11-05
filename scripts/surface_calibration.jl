# Surface-code lattice-surgery flip rate q(r) per basis and state independence (bravyi_cal noise; sd6, stim_generated variants).
# Run: julia --project=. scripts/surface_calibration.jl [--quick] [--workers 12]
# Output: results/taskA_julia/surface_qr.csv, surface_state_independence.csv, surface_calibration_meta.json
const REPO = dirname(@__DIR__)
get!(ENV, "JULIA_PYTHONCALL_EXE", joinpath(REPO, ".venv", "bin", "python"))
get!(ENV, "JULIA_CONDAPKG_BACKEND", "Null")
using Distributed, FTShadows

argval(flag, default) = (i = findfirst(==(flag), ARGS); i === nothing ? default : ARGS[i+1])
quick = "--quick" in ARGS
out = joinpath(REPO, "results", "taskA_julia")
mkpath(out)
target, max5, max7 = quick ? (50, 2e5, 1e5) : (300, 4e6, 2e6)
t0 = time()
wp = addprocs(parse(Int, argval("--workers", "12")); exeflags = "--project=$(Base.active_project())")
@everywhere wp using FTShadows
try
    pm = (f, xs) -> pmap(f, WorkerPool(wp), xs)
    cfgs = FTShadows.surface_qr_configs()
    rows = []
    for (d, mx) in ((5, max5), (7, max7))
        append!(rows, FTShadows.surface_qr_adaptive(filter(c -> c["d"] == d, cfgs), target, mx; pmap = pm))
        println("d=$d done $(round(Int, time() - t0))s")
    end
    FTShadows.pandas_csv(joinpath(out, "surface_qr.csv"), rows)
    si = FTShadows.surface_state_independence_configs()
    rows2 = FTShadows.surface_qr_adaptive(si, quick ? 50 : 10 * target, 2e6; pmap = pm)
    FTShadows.pandas_csv(joinpath(out, "surface_state_independence.csv"), rows2)
    foreach(println, rows)
    foreach(println, rows2)
    meta = FTShadows.OrderedDict("runtime_s" => time() - t0, "p" => FTShadows.SURFACE_QR_P, "quick" => quick)
    write(joinpath(out, "surface_calibration_meta.json"), FTShadows.pyjson(meta))
finally
    rmprocs(wp)
end
