# Generate the Nov-20 Bell / GHZ datasets.
# Output goes to handover/data_julia/nov20/ by default, never to handover/data/.
# Run: julia --project=. scripts/make_local_datasets.jl [--only SUBSTR[,SUBSTR..]] [--workers 6] [--out DIR]
const REPO = dirname(@__DIR__)
get!(ENV, "JULIA_PYTHONCALL_EXE", joinpath(REPO, ".venv", "bin", "python"))
get!(ENV, "JULIA_CONDAPKG_BACKEND", "Null")
using Distributed, FTShadows

argval(flag, default) = (i = findfirst(==(flag), ARGS); i === nothing ? default : ARGS[i+1])
only_ = argval("--only", "")
out = argval("--out", joinpath(REPO, "handover", "data_julia", "nov20"))
abspath(out) == abspath(joinpath(REPO, "handover", "data", "nov20")) && error("refusing to write into handover/data")
nworkers_ = parse(Int, argval("--workers", "6"))

# PythonCall is not thread-safe: chunks run on worker processes
wp = addprocs(nworkers_; exeflags = "--project=$(Base.active_project())")
@everywhere wp using FTShadows
calibration = FTShadows.csv_rows(FTShadows.upstream(REPO, "handover", "calibration_local.csv"))
version = try
    strip(read(`git -C $REPO rev-parse --short HEAD`, String))
catch
    "unknown"
end
try
    for ds in FTShadows.nov20_datasets()
        isempty(only_) || any(s -> occursin(s, ds["name"]), split(only_, ',')) || continue
        t0 = time()
        known = FTShadows.generate_nov20(ds, out; calibration, code_version = version,
                                         pmap = (f, xs) -> pmap(f, WorkerPool(wp), xs))
        println(ds["name"], ": F ideal ", known["fidelity_ideal"], ", expected ", known["fidelity_expected_noisy"],
                " ± ", known["fidelity_expected_noisy_se"], ", ", round(time() - t0; digits = 1), " s")
        flush(stdout)
    end
finally
    rmprocs(wp)
end
