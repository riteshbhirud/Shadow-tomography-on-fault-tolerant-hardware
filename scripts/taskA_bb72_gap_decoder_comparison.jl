# Gap decoders for the plain-anchored Relay gap on the same 1000 shots (Xbar, r = 3, p = 3e-3, sd6): s/shot, q, AUC, kept-shot q.
# Run: julia --project=. scripts/taskA_bb72_gap_decoder_comparison.jl [VARIANT ...] [--workers N]   (default all variants)
# Output: results/taskA_julia/bb72_gap_decoder_comparison.json
const REPO = dirname(@__DIR__)
get!(ENV, "JULIA_PYTHONCALL_EXE", joinpath(REPO, ".venv", "bin", "python"))
get!(ENV, "JULIA_CONDAPKG_BACKEND", "Null")
using Distributed, FTShadows
using FTShadows: OrderedDict

const OUT = joinpath(REPO, "results", "taskA_julia")
const SHOTS = 1000

i = findfirst(==("--workers"), ARGS)
only_ = i === nothing ? ARGS : vcat(ARGS[1:i-1], ARGS[i+2:end])  # optional subset of variant names
fout = joinpath(OUT, "bb72_gap_decoder_comparison.json")
names = isempty(only_) ? collect(keys(FTShadows.GAP_VARIANTS)) : only_
wp = addprocs(i === nothing ? length(names) : parse(Int, ARGS[i+1]); exeflags = "--project=$(Base.active_project())")
out = try
    @everywhere wp using FTShadows
    OrderedDict{String,Any}(pmap(n -> FTShadows.gap_variant(n; shots = SHOTS), WorkerPool(wp), names))
finally
    rmprocs(wp)
end
out_save = nothing
# a subset run is merged into the existing JSON (its later 'anchored_relay5' entry uses bravyi_cal noise, so it is not
# comparable with the rest)
if !isempty(only_) && isfile(fout)
    old = FTShadows.read_json_or_empty(fout)
    for (k, v) in out
        k != "anchored_relay100" && (old[k] = v)
    end
    out_save = old
end
base = haskey(out, "anchored_relay100") ? out["anchored_relay100"]["outcomes"] : nothing
for (k, v) in out
    o = pop!(v, "outcomes")
    base === nothing || (v["same_predictions_as_anchored100"] = o == base)
    println(k, " ", FTShadows.pyjson(v))
end
mkpath(OUT)
FTShadows.write_json(fout, out_save === nothing ? out : out_save)
