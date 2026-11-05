# Deformed-code and r = 6 circuit distances of the BB72 gauging combinations in the Bell / GHZ settings, with and without spectators.
# Run: julia --project=. scripts/taskA_bb72_local_distance.jl [--only 2|3] [--exclude N] [--skip-spectators] [--out PATH]
# Output: results/taskA_julia/bb72_local_distance.json (finished combinations are kept, so a rerun resumes)
const REPO = dirname(@__DIR__)
get!(ENV, "JULIA_PYTHONCALL_EXE", joinpath(REPO, ".venv", "bin", "python"))
get!(ENV, "JULIA_CONDAPKG_BACKEND", "Null")
using FTShadows

const OUT = joinpath(REPO, "results", "taskA_julia")
const R = 6

argval(flag, default) = (i = findfirst(==(flag), ARGS); i === nothing ? default : ARGS[i+1])
out = argval("--out", joinpath(OUT, "bb72_local_distance.json"))
mkpath(dirname(out))
only_ = argval("--only", nothing)
ms = only_ === nothing ? [2, 3] : [parse(Int, only_)]
skip_spect = "--skip-spectators" in ARGS
exclude = argval("--exclude", nothing)
res = FTShadows.read_json_or_empty(out)
skip = exclude === nothing ? Set{String}() : Set(keys(FTShadows.local_distance_combos([parse(Int, exclude)])))
for (key, bases) in sort(collect(FTShadows.local_distance_combos(ms)); by = kv -> (count(==(','), kv[1]), kv[1]))
    key in skip && continue
    row = get(res, key, FTShadows.OrderedDict{String,Any}())
    (haskey(row, "circuit_distance_all") || (skip_spect && haskey(row, "circuit_distance_L_only"))) && continue
    t = time()
    haskey(row, "deformed_distance_stim") || merge!(row, FTShadows.local_deformed_fields(bases, R))
    if !haskey(row, "circuit_distance_L_only")
        row["circuit_distance_L_only"] = FTShadows.local_distance_L(bases, R)
        row["seconds_L_only"] = round(Int, time() - t)
    end
    res[key] = row
    FTShadows.write_json(out, res)
    println(key, " ", FTShadows.pyjson(FTShadows.jsonable(row)))
    flush(stdout)
    skip_spect && continue
    t = time()
    row["n_targets_all"], row["circuit_distance_all"] = FTShadows.local_distance_all(bases, R)
    row["seconds_all"] = round(Int, time() - t)
    res[key] = row
    FTShadows.write_json(out, res)
    println(key, " ", FTShadows.pyjson(FTShadows.jsonable(row)))
    flush(stdout)
end
