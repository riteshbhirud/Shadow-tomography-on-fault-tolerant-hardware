# BB72 gap postselection with the plain-anchored Relay gap: Xbar_0, Zbar_0, Ybar_0 gauging at r = 2, 3, p = 1e-3.
# Run: julia --project=. scripts/taskA_bb72_postselection.jl [--workers 6] [--shots 30000]
# Output: results/taskA_julia/bb72_post_{P}_r{r}.npz and bb72_postselection.csv; resumes from bb72_post_chunks/
const REPO = dirname(@__DIR__)
get!(ENV, "JULIA_PYTHONCALL_EXE", joinpath(REPO, ".venv", "bin", "python"))
get!(ENV, "JULIA_CONDAPKG_BACKEND", "Null")
using Distributed, FTShadows

const OUT = joinpath(REPO, "results", "taskA_julia")
const CHUNK = 500

argval(flag, default) = (i = findfirst(==(flag), ARGS); i === nothing ? default : ARGS[i+1])
shots = parse(Int, argval("--shots", "30000"))
t0 = time()
cdir = joinpath(OUT, "bb72_post_chunks")
mkpath(cdir)
cfgs = [(string(Pl), r) for r in (3, 2) for Pl in "XZY"]
jobs = [(Pl, r, c, cdir, CHUNK) for (Pl, r) in cfgs for c in 0:shots÷CHUNK-1]
acc = FTShadows.OrderedDict(c => (Vector{UInt8}[], Vector{Float64}[]) for c in cfgs)
wp = addprocs(parse(Int, argval("--workers", "6")); exeflags = "--project=$(Base.active_project())")
try
    @everywhere wp using FTShadows
    for (k, (key, e, g)) in enumerate(pmap(FTShadows.bb72_post_chunk, WorkerPool(wp), jobs))
        push!(acc[key][1], e)
        push!(acc[key][2], g)
        (k - 1) % 20 == 0 && println(lpad(round(Int, time() - t0), 7), "s $k/$(length(jobs)) chunks")
    end
finally
    rmprocs(wp)
end
rows = []
for ((Pl, r), (es, gs)) in acc
    e, g = reduce(vcat, es), reduce(vcat, gs)
    FTShadows.npz_save_compressed(joinpath(OUT, "bb72_post_$(Pl)_r$(r).npz"); err = e, gap = g)
    push!(rows, FTShadows.postselection_row(Pl, r, e, g))
    println(FTShadows.pyjson(rows[end]))
end
FTShadows.write_frame_csv(joinpath(OUT, "bb72_postselection.csv"), rows)
write(joinpath(OUT, "bb72_postselection_meta.json"),
      FTShadows.pyjson(FTShadows.OrderedDict("runtime_s" => time() - t0, "shots" => shots)))
