# Flip rates of the compiled BB72 local settings per string and subset product (joint +1 eigenstates, p = 1e-3, Relay-BP).
# Run: julia --project=. scripts/taskA_bb72_local_calibration.jl [--ghz] [--shots 10000] [--shots-r6 5000] [--rounds 3 6] [--ghz-rounds 3] [--workers 8]
# Output: results/taskA_julia/bb72_local_q.csv; resumes from bb72_local_q_chunks/
const REPO = dirname(@__DIR__)
get!(ENV, "JULIA_PYTHONCALL_EXE", joinpath(REPO, ".venv", "bin", "python"))
get!(ENV, "JULIA_CONDAPKG_BACKEND", "Null")
using Distributed, FTShadows

const OUT = joinpath(REPO, "results", "taskA_julia")
const CHUNK = 1000

argval(flag, default) = (i = findfirst(==(flag), ARGS); i === nothing ? default : ARGS[i+1])
function argvals(flag, default)  # every value up to the next flag
    i = findfirst(==(flag), ARGS)
    i === nothing && return default
    j = something(findnext(startswith("--"), ARGS, i + 1), length(ARGS) + 1)
    return parse.(Int, ARGS[i+1:j-1])
end
shots, shots_r6 = parse(Int, argval("--shots", "10000")), parse(Int, argval("--shots-r6", "5000"))
rounds, ghz_rounds = argvals("--rounds", [3, 6]), argvals("--ghz-rounds", [3])
chunks = joinpath(OUT, "bb72_local_q_chunks")
mkpath(chunks)
sets = [(b, r) for b in FTShadows.product_strings("XYZ", 2) for r in rounds]
"--ghz" in ARGS && append!(sets, [(b, r) for b in FTShadows.product_strings("XYZ", 3) for r in ghz_rounds])
nch(r) = (r >= 6 ? shots_r6 : shots) ÷ CHUNK
# biggest circuits first so the pool tail is short
jobs = sort([(b, r, c) for (b, r) in sets for c in 0:nch(r)-1]; by = j -> (-j[2], -count(==('Y'), j[1]), j))
wp = addprocs(parse(Int, argval("--workers", "8")); exeflags = "--project=$(Base.active_project())")
try
    @everywhere wp using FTShadows
    for (i, _) in enumerate(pmap(j -> FTShadows.bb72_local_chunk((j..., chunks, CHUNK)), WorkerPool(wp), jobs))
        (i - 1) % 100 == 0 && println("$(i - 1)/$(length(jobs))")
    end
finally
    rmprocs(wp)
end
rows = []
for (bases, r) in sets
    o = reduce(vcat, [FTShadows.npy_load_u8(joinpath(chunks, "$(bases)_r$(r)_c$(c).npy")) for c in 0:nch(r)-1])
    append!(rows, FTShadows.bb72_local_rows(bases, r, o))
end
FTShadows.write_frame_csv(joinpath(OUT, "bb72_local_q.csv"), rows)
for x in rows
    x["kind"] == "string" && println(join((x[c] for c in ("setting_id", "r", "target", "shots", "flips", "q")), "  "))
end
