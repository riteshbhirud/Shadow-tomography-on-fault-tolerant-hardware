# Reproduce Bravyi et al.'s BB72 memory logical error rate (arXiv:2308.07915, Tables 1 and 6) with BP-OSD and Relay-BP.
# Run: julia --project=. scripts/calib_bravyi_reproduction.jl [--quick] [--workers 10] [--scaling 1.0] [--aggregate]
# Output: results/calibration_julia/bravyi_reproduction[_quick].json, chunks in bravyi_chunks_ms<scaling>/
const REPO = dirname(@__DIR__)
get!(ENV, "JULIA_PYTHONCALL_EXE", joinpath(REPO, ".venv", "bin", "python"))
get!(ENV, "JULIA_CONDAPKG_BACKEND", "Null")
using FTShadows
using FTShadows: OrderedDict, pyjson, pyfloat

argval(flag, default) = (i = findfirst(==(flag), ARGS); i === nothing ? default : ARGS[i+1])
quick = "--quick" in ARGS
scaling = parse(Float64, argval("--scaling", "1.0"))
worker = argval("--worker", nothing)
target, mx = quick ? (10, 2e4) : (100, 1e6)
cfgs = FTShadows.bravyi_configs()
outdir = joinpath(REPO, "results", "calibration_julia")
chunk_dir = joinpath(outdir, "bravyi_chunks_ms$(pyfloat(scaling))")
mkpath(chunk_dir)
if worker !== nothing  # crash-proof worker: claim chunks via exclusive file creation
    FTShadows.bravyi_run_worker(cfgs, target, mx, chunk_dir, scaling)
    exit()
end
if !("--aggregate" in ARGS)
    for f in readdir(chunk_dir)  # resume: drop claims whose chunk never finished
        endswith(f, ".claim") && !isfile(joinpath(chunk_dir, f[1:end-6] * ".json")) && rm(joinpath(chunk_dir, f))
    end
    procs = [run(`$(Base.julia_cmd()) --project=$(Base.active_project()) $(PROGRAM_FILE) --worker $i --scaling $(pyfloat(scaling))
                  $(quick ? ["--quick"] : [])`; wait = false) for i in 0:parse(Int, argval("--workers", "10"))-1]
    foreach(wait, procs)
end
rows = FTShadows.bravyi_rows(chunk_dir, cfgs)
out = OrderedDict{String,Any}("rows" => rows, "quick" => quick, "ms_scaling" => scaling,
                              "chunk_dir" => relpath(chunk_dir, REPO), "source" => "github.com/sbravyi/BivariateBicycleCodes@fa77e33")
write(joinpath(outdir, quick ? "bravyi_reproduction_quick.json" : "bravyi_reproduction.json"), pyjson(out; indent = 1))
println(pyjson(rows; indent = 1))
