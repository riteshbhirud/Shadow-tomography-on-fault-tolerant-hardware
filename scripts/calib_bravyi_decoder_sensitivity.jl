# Bravyi et al.'s BB72 memory circuit (Z basis, N_c = 6, p = 3e-3, 3000 shots) decoded with each BP variant on the same samples.
# Run: julia --project=. scripts/calib_bravyi_decoder_sensitivity.jl [--workers 6]
# Output: results/calibration_julia/bravyi_decoder_sensitivity.json
const REPO = dirname(@__DIR__)
get!(ENV, "JULIA_PYTHONCALL_EXE", joinpath(REPO, ".venv", "bin", "python"))
get!(ENV, "JULIA_CONDAPKG_BACKEND", "Null")
using Distributed, FTShadows

argval(flag, default) = (i = findfirst(==(flag), ARGS); i === nothing ? default : ARGS[i+1])
out = joinpath(REPO, "results", "calibration_julia", "bravyi_decoder_sensitivity.json")
mkpath(dirname(out))
# one process per variant unless --workers says otherwise
wp = addprocs(parse(Int, argval("--workers", string(length(FTShadows.BRAVYI_SENS_VARIANTS))));
              exeflags = "--project=$(Base.active_project())")
@everywhere wp using FTShadows
try
    rows = pmap(FTShadows.bravyi_sensitivity, WorkerPool(wp), FTShadows.BRAVYI_SENS_VARIANTS)
    for r in rows
        println(r["bp"], " ", something(r["scaling"], "None"), " ", r["fails"], " / ", r["shots"])
    end
    write(out, FTShadows.pyjson(rows))
finally
    rmprocs(wp)
end
