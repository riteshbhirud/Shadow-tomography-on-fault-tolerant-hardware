# Full parity with the saved Task B reference results in results/taskB (local only, not tracked).
# Every item is skipped with a message when that directory is missing. They need Oscar, hence :slow.

@testitem "TaskB parity: generators, group orders, groups records" tags = [:slow] setup = [TaskBRef] begin
    using Oscar, JSON, PythonCall
    import FTShadows.TaskB as TB
    if !TaskBRef.has_results()
        TaskBRef.skip_message("groups_*.json and generators_*.npz")
    else
        np = pyimport("numpy")
        for name in ("BB72", "gross")
            b = TaskBRef.build(name)
            ref = JSON.parsefile(joinpath(TaskBRef.RESULTS, "groups_$name.json"))
            foreach(e -> delete!(e, "seconds"), values(ref["ensembles"]))
            @test TaskBRef.same(b.record, ref)
            npz = np.load(joinpath(TaskBRef.RESULTS, "generators_$name.npz"))
            for en in TB.ENSEMBLES
                a = pyconvert(Array{UInt8,3}, npz[en])
                @test b.kept[en] == [TB.GF2Mat(a[i, :, :] .== 1) for i in 1:size(a, 1)]
            end
            # Oscar's own order agrees with the recog orders (the gross group with all folds is left out: GAP's
            # default method does not finish there)
            for en in TB.ENSEMBLES[2:end]
                name == "gross" && en == "aut_plus_all_folds" && continue
                @test TB.group_order(b.kept[en], 24; method=:oscar) == ref["ensembles"][en]["order"]
            end
        end
    end
end

@testitem "TaskB parity: Pauli corrections and signed checks" tags = [:slow] setup = [TaskBRef] begin
    using Oscar, JSON, Random
    import FTShadows.TaskB as TB
    if !TaskBRef.has_results()
        TaskBRef.skip_message("pauli_corrections_*.json and julia_export_*.json")
    else
        for name in ("BB72", "gross")
            b = TaskBRef.build(name)
            alls = TB.all_elements(b.ens)
            ref = JSON.parsefile(joinpath(TaskBRef.RESULTS, "pauli_corrections_$name.json"))
            @test TaskBRef.same([Dict("label" => e.label, "pauli_correction" => TaskBRef.gate_lists(e.pauli_correction))
                                 for e in alls], ref)
            # every corrected circuit maps every stabilizer to +stabilizer; the bare one flips a sign iff a
            # correction was needed
            for e in alls
                @test !any(TB.stabilizer_sign_flips(b.code.stabs, e.gates, b.code.n)[2])
                @test any(TB.stabilizer_sign_flips(b.code.stabs, e.gates_uncorrected, b.code.n)[2]) ==
                      !isempty(e.pauli_correction)
            end
            exp = JSON.parsefile(joinpath(TaskBRef.RESULTS, "julia_export_$name.json"))
            sample = TB.crosscheck_sample(b.ens, b.folds, b.kept)
            @test [e.label for e in sample] == [d["label"] for d in exp["elements"]]
            rng = Xoshiro(20261030)
            for (e, d) in zip(sample, exp["elements"])
                taus = TB.logical_signs(b.code, e)
                @test taus == d["logical_signs"]
                @test TaskBRef.gate_lists(e.gates) == d["gates"]
                @test e.logical == TB.GF2Mat(reduce(vcat, permutedims.([Bool.(r) for r in d["logical"]])))
                @test all(TB.signed_crosscheck(b.code, e, taus, rng))
            end
        end
    end
end

@testitem "TaskB parity: λ tables and shadow records" tags = [:slow] setup = [TaskBRef] begin
    using Oscar, JSON
    import FTShadows.TaskB as TB
    if !TaskBRef.has_results()
        TaskBRef.skip_message("lambda_mix_*.npy, shadow_*, weight1_detail.json and designed_cover_*.json")
    else
        w1 = Dict{String,Any}()
        for name in ("BB72", "gross")
            b = TaskBRef.build(name)
            orders = Dict(en => b.record["ensembles"][en]["order"] for en in TB.ENSEMBLES)
            rec, rows, lmix = TB.shadow_map(name, 12, b.record["basis"], b.kept, orders, TB.ENSEMBLE_DESCRIPTIONS)
            @test TaskBRef.same(rec, JSON.parsefile(joinpath(TaskBRef.RESULTS, "shadow_$name.json")))
            @test TaskBRef.csv_matches(rows, joinpath(TaskBRef.RESULTS, "shadow_$(name)_by_weight.csv"))
            for en in TB.ENSEMBLES
                ref = TaskBRef.read_npy(joinpath(TaskBRef.RESULTS, "lambda_mix_$(name)_$en.npy"))
                @test length(ref) == 4^12
                @test reinterpret(UInt32, lmix[en]) == reinterpret(UInt32, ref)
            end
            for en in ("albert_two_fold", "aut_plus_all_folds")
                w1["$name/$en"] = TB.weight1_detail(lmix[en])
            end
            cover = TB.designed_cover(name, b.kept, lmix, TB.cover_draws(name))
            @test TaskBRef.same(cover, JSON.parsefile(joinpath(TaskBRef.RESULTS, "designed_cover_$name.json")))
        end
        @test TaskBRef.same(w1, JSON.parsefile(joinpath(TaskBRef.RESULTS, "weight1_detail.json")))
    end
end

@testitem "TaskB parity: La-cross" tags = [:slow] setup = [TaskBRef] begin
    using Oscar, JSON
    import FTShadows.TaskB as TB
    if !TaskBRef.has_results()
        TaskBRef.skip_message("lacross.json")
    else
        ref = JSON.parsefile(joinpath(TaskBRef.RESULTS, "lacross.json"))
        delete!(ref, "seconds")
        rec = TB.lacross(3, 7)
        @test TaskBRef.same(rec, ref)
        @test BigInt(order(symplectic_group(18, GF(2)))) == ref["Sp_order"]
    end
end

@testitem "TaskB parity: basis independence" tags = [:slow] setup = [TaskBRef] begin
    using Oscar, JSON
    import FTShadows.TaskB as TB
    dir = joinpath(TaskBRef.RESULTS, "basis_independence")
    if !isdir(dir)
        TaskBRef.skip_message("basis_independence/")
    else
        strip_times(d) = (d = deepcopy(d); delete!(d, "seconds_total");
                          haskey(d, "ensembles") && foreach(e -> delete!(e, "seconds"), values(d["ensembles"])); d)
        for name in ("BB72", "gross")
            b = TaskBRef.build(name)
            kept = TB.kept_elements(b.ens, b.kept)
            labels = TB.basis_labels(name)
            data = TB.basis_independence(name, b.code, kept; labels)
            for lb in labels
                @test TaskBRef.same(data[lb], strip_times(JSON.parsefile(joinpath(dir, "$(name)_basis_$lb.json"))))
            end
            summ, rows, spread, spec = TB.aggregate_bases(name, data, labels)
            @test TaskBRef.same(summ, JSON.parsefile(joinpath(dir, "$(name)_summary.json")))
            @test TaskBRef.csv_matches(rows, joinpath(dir, "$(name)_per_weight.csv"))
            @test TaskBRef.csv_matches(spread, joinpath(dir, "$(name)_per_weight_spread.csv"))
            @test TaskBRef.csv_matches(spec, joinpath(dir, "$(name)_spectrum.csv"))
        end
    end
end
