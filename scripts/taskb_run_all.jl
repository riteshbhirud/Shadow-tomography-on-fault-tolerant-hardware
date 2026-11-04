# Reproduce every Task B result (fixed seeds).
# Steps: build the groups of BB72 and gross, the QuantumClifford sign check, the shadow maps
# (with their plots), weight1_detail, designed_cover, lacross, basis_independence, then the Task B tests.
# Run: julia --project=scripts/oscar scripts/taskb_run_all.jl [--no-tests] [--out DIR]   (default results/taskB_julia/)
# The group orders need Oscar, a weak dependency of FTShadows, so the script runs in the small scripts/oscar
# environment (FTShadows from this repository plus Oscar). The plots need the plotting environment:
# julia --project=scripts/plots -e 'using Pkg; Pkg.instantiate()'.
const REPO = dirname(@__DIR__)
get!(ENV, "JULIA_PYTHONCALL_EXE", joinpath(REPO, ".venv", "bin", "python"))
get!(ENV, "JULIA_CONDAPKG_BACKEND", "Null")
using Distributed, Oscar, FTShadows
import FTShadows.TaskB as TB

argval(flag, default) = (i = findfirst(==(flag), ARGS); i === nothing ? default : ARGS[i+1])
const OUT = argval("--out", joinpath(REPO, "results", "taskB_julia"))
const CODES = ("BB72", "gross")
mkpath(OUT)

# run one step with its printed output going to OUT/log as well, so each step keeps its log
function step(f, title, log)
    println("+ ", title)
    t0 = time()
    path = joinpath(OUT, log)
    open(path, "w") do io
        redirect_stdout(io) do
            f()
        end
    end
    print(read(path, String))
    println("  ($(round(time() - t0; digits = 1)) s)")
end

# build_groups: groups_<code>.json, generators_<code>.npz, julia_export_<code>.json, pauli_corrections_<code>.json
built = Dict{String,Any}()
step("build_groups BB72 gross", "build_groups.log") do
    for name in CODES
        t0 = time()
        rec, kept, ens, folds, code = TB.build_groups(name; timings = true)
        for en in TB.ENSEMBLES[2:end]
            e = rec["ensembles"][en]
            println("[$name] $en: order $(e["order"]) (seed2 $(e["order_second_seed"]), bfs $(something(e["order_bfs_exact"], "None")))")
        end
        sample = TB.crosscheck_sample(ens, folds, kept)
        TB.write_pyjson(joinpath(OUT, "pauli_corrections_$name.json"), TB.pauli_corrections_record(TB.all_elements(ens)))
        TB.write_pyjson(joinpath(OUT, "groups_$name.json"), rec; indent = 1)
        TB.write_generators_npz(joinpath(OUT, "generators_$name.npz"), kept, TB.logical_k(code))
        TB.write_pyjson(joinpath(OUT, "julia_export_$name.json"), TB.export_record(name, code, sample))
        built[name] = (; rec, kept, ens, folds, code, sample)
        println("[$name] done in $(round(Int, time() - t0))s")
    end
end

# the state-level sign check with QuantumClifford (the old julia/taskb_crosscheck.jl)
step("QuantumClifford sign check", "julia_crosscheck.log") do
    ok = true
    for name in CODES
        b = built[name]
        text, allok = TB.crosscheck_report(name, b.code, b.sample)
        print(text)
        ok &= allok
    end
    ok || error("Julia cross-check FAILED")
end

# shadow maps: shadow_<code>.json, shadow_<code>_by_weight.csv, lambda_mix_<code>_<ensemble>.npy
lmix = Dict{String,Any}()
step("shadow_map BB72 gross", "shadow_map.log") do
    for name in CODES
        b = built[name]
        orders = Dict(en => b.rec["ensembles"][en]["order"] for en in TB.ENSEMBLES)
        descr = Dict(en => b.rec["ensembles"][en]["description"] for en in TB.ENSEMBLES)
        rec, rows, lm = TB.shadow_map(name, 12, b.rec["basis"], b.kept, orders, descr)
        for en in TB.ENSEMBLES
            TB.write_npy(joinpath(OUT, "lambda_mix_$(name)_$en.npy"), lm[en])
            e = rec["ensembles"][en]
            println("[$name] $en: $(e["n_nontrivial_orbits"]) orbits, mix learnable frac ",
                    round(e["readouts"]["mix"]["all_weights"]["frac_learnable"]; sigdigits = 4))
        end
        TB.write_pyjson(joinpath(OUT, "shadow_$name.json"), rec; indent = 1)
        TB.write_dictcsv(joinpath(OUT, "shadow_$(name)_by_weight.csv"), rows)
        lmix[name] = lm
    end
    # the plots run in the plotting environment (matplotlib through the bridge)
    plots = joinpath(REPO, "scripts", "plots")
    run(`$(Base.julia_cmd()) --project=$plots $(joinpath(plots, "taskb_shadow_plots.jl")) $OUT $CODES`)
end

step("weight1_detail", "weight1_detail.log") do
    out = FTShadows.OrderedDict{String,Any}()
    for name in CODES, en in ("albert_two_fold", "aut_plus_all_folds")
        out["$name/$en"] = TB.weight1_detail(lmix[name][en])
        println(name, " ", en, " unlearnable w=1: ", out["$name/$en"]["unlearnable_weight1"],
                " learnable w=1 norm2 values: ", out["$name/$en"]["weight1_norm2_values"])
    end
    TB.write_pyjson(joinpath(OUT, "weight1_detail.json"), out; indent = 1)
end

step("designed_cover BB72 gross", "designed_cover.log") do
    for name in CODES
        rec = TB.designed_cover(name, built[name].kept, lmix[name], TB.cover_draws(name))
        for en in TB.ENSEMBLES
            e = rec["ensembles"][en]
            println("[$name] $en: targets $(e["n_learnable_targets"])/$(e["n_targets_total"]), native cover ",
                    "$(e["greedy_cover_size"]) (covered $(e["n_covered"])), local-Pauli cover $(e["local_pauli_cover_same_targets"])")
        end
        TB.write_pyjson(joinpath(OUT, "designed_cover_$name.json"), rec; indent = 1)
    end
end

step("lacross", "lacross.log") do
    rec = TB.lacross(3, 7; timings = true)
    for (en, e) in rec["ensembles"]
        println("[$(rec["code"])] $en: order $(e["order"]) full Sp(2k)=$(e["is_full_Sp"]) orbits $(e["n_nontrivial_orbits"]) ",
                "learnable $(round(e["frac_learnable_mix"]; digits = 4))")
    end
    TB.write_pyjson(joinpath(OUT, "lacross.json"), rec; indent = 1)
end

# basis independence on 3 workers, one (code, basis) per task: basis_independence/<code>_basis_<label>.json, then
# <code>_summary.json and the three CSV tables
step("basis_independence --workers 3", "basis_independence.log") do
    dir = joinpath(OUT, "basis_independence")
    mkpath(dir)
    wp = addprocs(3; exeflags = "--project=$(Base.active_project())")
    try
        Distributed.remotecall_eval(Main, wp, :(using FTShadows))
        tasks = [(name, lb) for name in CODES for lb in TB.basis_labels(name)]
        args = [(name, built[name].kept, lb) for (name, lb) in tasks]
        recs = pmap(a -> FTShadows.TaskB.basis_record(a...), WorkerPool(wp), args)
        for ((name, lb), r) in zip(tasks, recs)
            haskey(r, "failed") && println("[$name] basis $lb: ", r["failed"])
            TB.write_pyjson(joinpath(dir, "$(name)_basis_$lb.json"), r)
        end
        for name in CODES
            labels = TB.basis_labels(name)
            data = Dict(lb => r for ((nm, lb), r) in zip(tasks, recs) if nm == name)
            summ, rows, spread, spec = TB.aggregate_bases(name, data, labels)
            TB.write_pyjson(joinpath(dir, "$(name)_summary.json"), summ; indent = 1)
            TB.write_dictcsv(joinpath(dir, "$(name)_per_weight.csv"), rows)
            TB.write_dictcsv(joinpath(dir, "$(name)_per_weight_spread.csv"), spread)
            TB.write_dictcsv(joinpath(dir, "$(name)_spectrum.csv"), spec)
            println("[$name] invariants identical in all $(length(summ["bases"])) bases: ",
                    summ["all_invariants_identical_all_ensembles"], "; CSS claim: ",
                    "{", join(("'$e': " * (v["holds_in_all_bases"] ? "True" : "False") for (e, v) in summ["css_claim_check"]), ", "), "}")
        end
    finally
        rmprocs(wp)
    end
end

if !("--no-tests" in ARGS)
    step("Task B tests", "test_output.txt") do
        cmd = addenv(`$(Base.julia_cmd()) --project=$REPO -e "using Pkg; Pkg.test()"`, "FTSHADOWS_TEST_FILES" => "taskb")
        run(pipeline(cmd; stdout = stdout, stderr = stdout))
    end
end
