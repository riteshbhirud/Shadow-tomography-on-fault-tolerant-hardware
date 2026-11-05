# Generate the v1.0/v1.1 datasets: <name>.meta.json, .records.jsonl (one shot per line) and .exact.csv.
# Run: julia --project=. scripts/make_datasets.jl [--only NAME] [--draws N] [--bb-draws N] [--workers 3] [--out DIR]
# Output: handover/data_julia/ by default, never handover/data/.
const REPO = dirname(@__DIR__)
get!(ENV, "JULIA_PYTHONCALL_EXE", joinpath(REPO, ".venv", "bin", "python"))
get!(ENV, "JULIA_CONDAPKG_BACKEND", "Null")
using Distributed, FTShadows

argval(flag, default) = (i = findfirst(==(flag), ARGS); i === nothing ? default : ARGS[i+1])
only_ = argval("--only", "")
draws = parse(Int, argval("--draws", "20000"))
bb_draws = parse(Int, argval("--bb-draws", "10000"))
out = argval("--out", joinpath(REPO, "handover", "data_julia"))
abspath(out) == abspath(joinpath(REPO, "handover", "data")) && error("refusing to write into handover/data")

# PythonCall is not thread-safe: the settings of a dataset run on worker processes
wp = addprocs(parse(Int, argval("--workers", "3")); exeflags = "--project=$(Base.active_project())")
@everywhere wp using FTShadows
version = try
    strip(read(`git -C $REPO rev-parse --short HEAD`, String))
catch
    "unknown"
end
try
    for spec in FTShadows.v1_datasets(draws, bb_draws)
        isempty(only_) || occursin(only_, spec["name"]) || continue
        t0 = time()
        n = FTShadows.generate_v1(spec, out; code_version = version, pmap = (f, xs) -> pmap(f, WorkerPool(wp), xs))
        println(spec["name"], ": ", spec["draws"], " draws, ", n, " learnable targets, ", round(time() - t0; digits = 1), " s")
        flush(stdout)
    end
finally
    rmprocs(wp)
end
