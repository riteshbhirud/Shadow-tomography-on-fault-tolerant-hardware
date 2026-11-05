# Smallest gauging graph (2|E| + 1 extra qubits) with circuit distance 6 for the weight-6 Xbar_0 / Zbar_0 of BB72.
# Run: julia --project=. scripts/taskA_bb72_small_graph_search.jl [--spectator-screen] [--workers 3]
# Output: results/taskA_julia/bb72_small_graph_search*.json
const REPO = dirname(@__DIR__)
get!(ENV, "JULIA_PYTHONCALL_EXE", joinpath(REPO, ".venv", "bin", "python"))
get!(ENV, "JULIA_CONDAPKG_BACKEND", "Null")
using Distributed, FTShadows
using PythonCall: pyconvert
using FTShadows: OrderedDict, small_graph_candidates, max_degree_of, Circuits

const OUT = joinpath(REPO, "results", "taskA_julia")
const R = 6
const N_HI, N_LO, N_SPECT = 200, 20, 40

i = findfirst(==("--workers"), ARGS)
workers = i === nothing ? 3 : parse(Int, ARGS[i+1])
mkpath(OUT)
wp = addprocs(workers; exeflags = "--project=$(Base.active_project())")
@everywhere wp using FTShadows
screen(jobs) = pmap(j -> FTShadows.small_graph_job(j...; rounds = R), WorkerPool(wp), jobs)
# numpy's rng.choice(n, k, replace=False), sorted and 1-based, so the samples follow the fixed seeds
pick(rng, n, k) = sort(pyconvert(Vector{Int}, rng.choice(n, k; replace = false))) .+ 1
histogram(rows) = OrderedDict(string(d) => count(r -> r["distance"] == d, rows) for d in sort(unique(r["distance"] for r in rows)))

# Second pass (the L-only screen is too weak: every 5-edge star passes it but has spectator circuit
# distance 4 / 5). For |E| = 6, 7, 8: a seeded sample of 40 h >= 1 graphs, preferring the lowest max
# degree, screened with all readable spectator logicals; stop at the first |E| with a d = 6 graph.
function spectator_screen()
    res = OrderedDict{String,Any}()
    for P in "XZ"
        cands = small_graph_candidates(P)
        rng = Circuits.np_rng(20261031)
        res[string(P)] = []
        for E in (6, 7, 8)
            hi = [es for (e, es, h) in cands if e == E && h >= 1]
            degs = max_degree_of.(hi)
            dmin = minimum(degs)
            pool_ = [es for (es, d) in zip(hi, degs) if d <= dmin + 1]
            sel = pool_[pick(rng, length(pool_), min(N_SPECT, length(pool_)))]
            rows = screen([(P, es, true) for es in sel])
            n6 = count(r -> r["distance"] >= 6, rows)
            dist = histogram(rows)
            println(P, " E ", E, " tested ", length(rows), " of ", length(hi), " h>=1 (max degree <= ", dmin + 1, ") distance histogram ", dist)
            flush(stdout)
            push!(res[string(P)], OrderedDict{String,Any}("E" => E, "n_h_ge_1" => length(hi), "tested" => length(rows),
                "max_degree_le" => dmin + 1, "distance_histogram" => dist, "rows" => rows))
            FTShadows.write_json(joinpath(OUT, "bb72_small_graph_search_spectators.json"), res)
            n6 > 0 && break
        end
    end
end

function main()
    res = OrderedDict{String,Any}()
    for P in "XZ"
        cands = small_graph_candidates(P)
        counts = OrderedDict(string(E) => count(c -> c[1] == E, cands) for E in 5:8)
        hcounts = OrderedDict(string(E) => count(c -> c[1] == E && c[3] >= 1, cands) for E in 5:8)
        println(P, " connected graphs per |E| ", counts, " with h>=1 ", hcounts)
        flush(stdout)
        # Level by level from |E| = 5: all h >= 1 graphs (seeded random sample of 200 if more), plus 20
        # h < 1 graphs as a check of the h >= 1 rule; stop at the first |E| with d = 6.
        rng = Circuits.np_rng(20261030)
        rows = []
        for E in 5:8
            hi = [c for c in cands if c[1] == E && c[3] >= 1]
            lo = [c for c in cands if c[1] == E && c[3] < 1]
            length(hi) > N_HI && (hi = hi[pick(rng, length(hi), N_HI)])
            length(lo) > N_LO && (lo = lo[pick(rng, length(lo), N_LO)])
            level = screen([(P, es, false) for (_, es, _) in vcat(hi, lo)])
            append!(rows, level)
            n6 = count(r -> r["distance"] >= 6, level)
            println(P, " E ", E, " screened ", length(level), " (h>=1: ", length(hi), ") d=6: ", n6, " d=6 among h<1: ",
                    count(r -> r["distance"] >= 6 && r["cheeger"] < 1, level))
            flush(stdout)
            n6 > 0 && break
        end
        good = sort([r for r in rows if r["distance"] >= 6]; by = r -> (r["E"], r["max_check_weight"]), alg = MergeSort)
        res[string(P)] = OrderedDict{String,Any}("counts" => counts, "h_ge_1" => hcounts, "screened" => length(rows),
                                                 "screen_rows" => rows)
        if !isempty(good)
            Emin = good[1]["E"]
            top = [r for r in good if r["E"] == Emin][1:min(3, end)]
            res[string(P)]["smallest_E_with_d6_L_only"] = Emin
            res[string(P)]["spectator_checks"] = screen([(P, [Tuple(e) for e in r["edges"]], true) for r in top])
        end
        println(P, " best ", good[1:min(3, end)])
        FTShadows.write_json(joinpath(OUT, "bb72_small_graph_search.json"), res)
    end
end

try
    "--spectator-screen" in ARGS ? spectator_screen() : main()
finally
    rmprocs(wp)
end
