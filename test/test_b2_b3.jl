# The B2 cost comparison and the B3 La-cross scan reproduce the stored results in results/b2 and results/b3.
# Those folders are local (not tracked), so the parity items skip when they are absent.

@testitem "B2 prep-vs-measure table equals the reference" tags = [:slow] begin
    using FTShadows
    REPO = pkgdir(FTShadows)
    ref = joinpath(REPO, "results", "b2", "prep_vs_measure.csv")
    costs = joinpath(REPO, "results", "taskA", "local_setting_costs.csv")
    if !(isfile(ref) && isfile(costs))
        @info "results/b2 or results/taskA not found: skipping the B2 parity check"
    else
        rows = FTShadows.prep_vs_measure_rows(costs, joinpath(REPO, "scripts", "data", "b2_cost_sources.json"))
        out = tempname()
        FTShadows.write_prep_vs_measure(out, rows)
        @test read(out) == read(ref)
    end
end

@testitem "B2 helpers" begin
    using FTShadows
    @test FTShadows.pyrepr_dict(['X' => 9, 'Y' => 76]) == "{'X': 9, 'Y': 76}"
    @test FTShadows.cost_volume(Dict("data_qubits" => 9, "ancilla_qubits" => 8, "merged_rounds" => 3)) == 68
end

@testitem "B3 La-cross code and circuits equal results/b3" tags = [:slow] begin
    using FTShadows, JSON
    dir = joinpath(pkgdir(FTShadows), "results", "b3")
    if !isfile(joinpath(dir, "code.json"))
        @info "results/b3 not found: skipping the B3 parity check"
    else
        # timings are the only fields that differ between runs
        strip!(d) = (delete!(d, "seconds"); delete!(d["distance_stim_code_capacity"], "seconds"); d)
        rec, arrays = FTShadows.lacross_code_record()
        @test FTShadows.pyjson(strip!(rec); indent = 1) ==
              FTShadows.pyjson(strip!(JSON.parsefile(joinpath(dir, "code.json"))); indent = 1)
        circ = FTShadows.lacross_memory(FTShadows.LaCross(), 1e-3, 4; basis = "Z", noise = "uniform")
        @test FTShadows.Circuits.stim_text(circ) * "\n" == read(joinpath(dir, "lacross65_memory_uniform_Z_r4_p1e-3.stim"), String)
    end
end

@testitem "B3 La-cross pieces" tags = [:bridge] begin
    using FTShadows, PythonCall, Random
    c = FTShadows.LaCross()
    @test (c.N, size(c.HX, 1), c.N - 2 * FTShadows.Circuits.gf2rank(c.HX)) == (65, 28, 9)
    @test FTShadows.gamma_order(c) == 168
    @test FTShadows.exhaustive_min_logical(c.HZ, c.HX, 3)[1] === nothing
    @test FTShadows.exhaustive_min_logical(c.HZ, c.HX, 4)[1] == 4
    @test FTShadows.cnot_count_gauss(Bool[0 1; 1 0]) == (3, 0)
    f = FTShadows.fold_checks(c)
    @test all(f[k]["preserves_stabilizer_group"] && f[k]["logical_action_matches_paper_up_to_stabilizers_and_sign"]
              for k in ("hadamard_fold", "phase_fold"))
    # the ISD search equals the original implementation's (same numpy permutations); needs its local
    # copy (see runtests.jl)
    if isdir(joinpath(pkgdir(FTShadows), "reference", "python", "src"))
        gf2 = pyimport("ftoracle.gf2")
        np = pyimport("numpy")
        rng = Xoshiro(7)
        for trial in 1:4
            A = rand(rng, Bool, 6, 14)
            B = rand(rng, Bool, 2, 14)
            w, _ = FTShadows.min_weight_codeword(A, B; iters = 20, seed = trial)
            pw = gf2.min_weight_codeword(np.array(UInt8.(A)), np.array(UInt8.(B)); iters = 20, seed = trial)[0]
            @test w == pyconvert(Int, pw)
        end
    end
    # a noiseless memory circuit is deterministic and silent
    py = FTShadows.Circuits.pycircuit(FTShadows.lacross_memory(c, 1e-3, 2; basis = "X", noise = "agnostic"))
    @test FTShadows.lacross_noiseless_check(py; shots = 50) == (0, 0)
end
