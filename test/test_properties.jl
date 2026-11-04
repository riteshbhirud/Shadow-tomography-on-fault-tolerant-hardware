# Statistical and structural properties of the oracle: noiseless exactness on random states, the bias
# formula under noise, seed independence, confidence kinds and the decoder gap. States, seeds, shots and
# thresholds are fixed; stim samples the same circuits with the same seeds, so the samples are reproducible.

@testmodule Props begin
    using FTShadows, PythonCall
    using FTShadows.Circuits: exact_expectation

    const ALPHA = 1e-3
    const STATES = joinpath(pkgdir(FTShadows), "test", "reference", "circuits", "states.json")

    st(state) = [(String(l), Int(s)) for (l, s) in state]

    """Phase-free product of two Pauli labels (phases are handled by exact_expectation)."""
    function mul(a, b)
        t(p, q) = p == 'I' ? q : q == 'I' ? p : p == q ? 'I' : only(setdiff("XYZ", (p, q)))
        return join(t(p, q) for (p, q) in zip(a, b))
    end

    anticommute(a, b) = isodd(count(x != 'I' && y != 'I' && x != y for (x, y) in zip(a, b)))

    """Nonempty subsets of 0:m-1 in itertools order: by size, then lexicographic."""
    function subsets(m)
        out = Vector{Vector{Int}}()
        function rec(cur, start, r)
            if length(cur) == r
                push!(out, copy(cur))
                return
            end
            for i in start:m-1
                push!(cur, i)
                rec(cur, i + 1, r)
                pop!(cur)
            end
        end
        for r in 1:m
            rec(Int[], 0, r)
        end
        return out
    end

    """Two-sided exact binomial p-value (scipy.stats.binomtest, through the bridge)."""
    binom_p(k, n, p=0.5) = pyconvert(Float64, pyimport("scipy.stats").binomtest(k, n, p).pvalue)

    """Every nonempty subset product of the setting's strings against the exact expectation (first 64)."""
    function check_setting!(state, labels, outcomes, pvals, signs_seen)
        for sub in first(subsets(length(labels)), 64)
            P = foldl(mul, (labels[j+1] for j in sub[2:end]); init=labels[sub[1]+1])
            all(==('I'), P) && continue
            ex = exact_expectation(state, P)
            bits = reduce(.⊻, (outcomes[:, j+1] for j in sub))
            if ex != 0
                all(bits .== (ex == 1 ? 0 : 1)) || return false
                push!(get!(signs_seen, P, Set{Int}()), ex)
            else
                push!(pvals, binom_p(Int(sum(bits)), length(bits)))
            end
        end
        return true
    end

    """Holm check: every sorted p-value above alpha / (m - i)."""
    function holm_ok(pvals; alpha=ALPHA)
        p = sort(pvals)
        m = length(p)
        return all(p[i] > alpha / (m - i + 1) for i in 1:m)
    end

    """Two-sided normal p-value 2 sf(|z|)."""
    norm_p(z) = FTShadows.erfc(abs(z) / sqrt(2))

    """[(u, v, om, weight)] of the graphlike edges of a stim DEM, v = -1 for boundary edges."""
    function weighted_edges(dem)
        E = FTShadows.graphlike_edges(FTShadows.read_dem(dem))
        return [(u, v, om, FTShadows.llr(q)) for ((u, v, om), q) in zip(E.keys, E.p)]
    end
end

@testitem "day 2: surface random states, noiseless, match the exact values" tags = [:properties, :bridge] setup = [Props] begin
    using JSON
    using FTShadows.Circuits: setting_experiment, exact_expectation
    using .Props: st, check_setting!, holm_ok, anticommute
    # the random states and the numpy-drawn local settings are the exported ones
    S = JSON.parsefile(Props.STATES)
    for np in (2, 3)
        pvals, signs_seen, single = Float64[], Dict{String,Set{Int}}(), Dict{String,Set{Int}}()
        for t in 0:7
            state = st(S[iseven(t) ? "stab:$np:$(1000np + t)" : "prod:$np:$(2000np + t)"])
            for setting in vcat(["surgery:ZZ", "surgery:XX"], S["day2_settings:$np"][t+1])
                exp, labels = setting_experiment(setting, 3, np, state, 2)
                res = run_oracle(exp, 0.0, 400; seed=t, setting_labels=labels)
                @test check_setting!(state, labels, res.outcomes, pvals, signs_seen)
                # sign coverage: every deterministic string is rerun on the state conjugated by a single-qubit
                # Pauli that anticommutes with it, which flips its sign
                for lab in labels
                    ex = exact_expectation(state, lab)
                    ex == 0 && continue
                    i = findfirst(!=('I'), lab)
                    Q = join(j == i ? (lab[i] == 'X' ? 'Z' : 'X') : 'I' for j in 1:np)
                    flipped = [(g, sg ⊻ Int(anticommute(g, Q))) for (g, sg) in state]
                    @test exact_expectation(flipped, lab) == -ex
                    exp2, labels2 = setting_experiment(setting, 3, np, flipped, 2)
                    res2 = run_oracle(exp2, 0.0, 400; seed=t, setting_labels=labels2)
                    @test check_setting!(flipped, labels2, res2.outcomes, pvals, signs_seen)
                    union!(get!(single, lab, Set{Int}()), (ex, -ex))
                end
            end
        end
        @test holm_ok(pvals)
        @test !isempty(single) && all(==(Set([1, -1])), values(single))
    end
end

@testitem "day 2: BB72 random states, noiseless, match the exact values" tags = [:properties, :bridge, :slow] setup = [Props] begin
    using JSON
    using FTShadows.Circuits: bb72_readout_experiment, bb_gauging
    using .Props: st, check_setting!, holm_ok
    S = JSON.parsefile(Props.STATES)
    pvals, signs_seen = Float64[], Dict{String,Set{Int}}()
    for t in 0:3
        state = st(S[iseven(t) ? "stab:12:$(77 + t)" : "prod:12:$(88 + t)"])
        for kind in "XZ"
            exp, labels = bb72_readout_experiment(kind, state)
            res = run_oracle(exp, 0.0, 400; seed=t, setting_labels=labels)
            @test check_setting!(state, labels, res.outcomes, pvals, signs_seen)
        end
        for P in "XZY"
            res = run_oracle(bb_gauging(P, 1; state_labels=state), 0.0, 400; decoder="relay", seed=t)
            @test check_setting!(state, [P * "I"^11], res.outcomes, pvals, signs_seen)
        end
    end
    @test holm_ok(pvals)
    @test any(==(Set([1, -1])), values(signs_seen))
end

@testitem "day 2: bias formula under noise" tags = [:properties] setup = [Props] begin
    using Statistics
    using FTShadows.Circuits: setting_experiment, exact_expectation
    # E[(-1)^b] = (1 - 2q) <P> for <P> = +1, -1, 0
    # (ZZ surgery, d = 3, r = 3, p = 2e-3), q calibrated on an independent sample; z-tests, Holm at 1e-3
    p, shots = 2e-3, 40000
    cal = run_oracle(setting_experiment("surgery:ZZ", 3, 2, [("ZZ", 0), ("XX", 0)], 3)[1], p, shots; seed=101, with_gap=false)
    q = mean(cal.outcomes)
    @test 0 < q < 0.05
    pv = Float64[]
    for (state, ev, seed) in (([("ZZ", 0), ("XX", 0)], 1, 202), ([("ZZ", 1), ("XX", 0)], -1, 303), ([("XI", 0), ("IX", 0)], 0, 404))
        res = run_oracle(setting_experiment("surgery:ZZ", 3, 2, state, 3)[1], p, shots; seed, with_gap=false)
        x = 1 .- 2 .* Float64.(res.outcomes[:, 1])
        @test exact_expectation(state, "ZZ") == ev
        mu = (1 - 2q) * ev
        var_mu = (4q * (1 - q) / shots) * ev^2  # uncertainty of the calibrated q
        se = sqrt(var(x) / shots + var_mu)
        push!(pv, Props.norm_p((mean(x) - mu) / se))
    end
    @test Props.holm_ok(pv)
end

@testitem "day 2: same seed is bit-identical, other seeds independent" tags = [:properties] begin
    using Statistics
    using FTShadows.Circuits: setting_experiment
    exp, _ = setting_experiment("surgery:ZZ", 3, 2, [("XI", 0), ("IX", 0)], 3)
    a = run_oracle(exp, 1e-3, 4000; seed=9)
    b = run_oracle(exp, 1e-3, 4000; seed=9)
    @test a.outcomes == b.outcomes && a.raw == b.raw
    @test isequal(a.confidence, b.confidence)
    c = run_oracle(exp, 1e-3, 4000; seed=10)
    r = cor(Float64.(a.raw[:, 1]), Float64.(c.raw[:, 1]))  # random outcomes (<ZZ> = 0): independent
    @test abs(r) < 4 / sqrt(4000)
end

@testitem "day 2: worker processes give independent streams" tags = [:properties, :slow] begin
    using Distributed, Statistics
    # Distributed workers, as in the scripts (they run chunks on worker processes; PythonCall is not thread-safe)
    seeds = [seed_for("worker_test", i) for i in 0:3]
    @test length(unique(seeds)) == 4
    wp = addprocs(2; exeflags="--project=$(Base.active_project())")
    try
        @everywhere begin
            using FTShadows
            function worker_raw(seed)
                exp, _ = FTShadows.Circuits.setting_experiment("surgery:XX", 3, 2, [("ZI", 0), ("IZ", 0)], 2)
                return run_oracle(exp, 1e-3, 3000; seed, with_gap=false).raw[:, 1]
            end
        end
        f = invokelatest(getglobal, Main, :worker_raw)
        outs = pmap(f, WorkerPool(wp), seeds)
        for i in 1:4, j in i+1:4
            @test abs(cor(Float64.(outs[i]), Float64.(outs[j]))) < 4 / sqrt(3000)
        end
        # the same seed on two different workers gives identical shots
        @test remotecall_fetch(f, wp[1], seeds[1]) == remotecall_fetch(f, wp[2], seeds[1])
    finally
        rmprocs(wp)
    end
end

@testitem "day 2: confidence kinds consistent with values" tags = [:properties] begin
    using FTShadows.Circuits: multi_patch_surgery, bb_gauging
    bell = [("ZZ", 0), ("XX", 0)]
    r1 = run_oracle(multi_patch_surgery(3, 2, "ZZ", bell, 1), 1e-3, 300; seed=1)
    @test Set(r1.confidence_kind) == Set(["unavailable"]) && all(isnan, r1.confidence)
    r3 = run_oracle(multi_patch_surgery(3, 2, "ZZ", bell, 3), 1e-3, 300; seed=1)
    @test Set(r3.confidence_kind) == Set(["exact"]) && all(r3.confidence .>= 0)
    r0 = run_oracle(multi_patch_surgery(3, 2, "ZZ", bell, 3), 0.0, 50; seed=1)
    @test Set(r0.confidence_kind) == Set(["noiseless"]) && all(isinf, r0.confidence)
    rb = run_oracle(bb_gauging('X', 1; state_kind="+1"), 2e-3, 80; decoder="relay", seed=1)
    k, c = rb.confidence_kind, rb.confidence
    @test issubset(Set(k), Set(["approx", "not_found", "decoder_failed"]))
    @test all(isinf, c[k .== "not_found"]) && all(isfinite, c[k .== "approx"])
    @test all(c[k .== "decoder_failed"] .== 0)
    rec = records(r3)[1]
    @test rec["confidence_kind"] == "exact" && haskey(rec, "meta")
end

@testitem "day 2: anchored-gap predictions equal plain Relay-BP" tags = [:properties, :slow] begin
    using FTShadows.Circuits: bb_gauging
    for P in "XY"
        exp = bb_gauging(P, 2; state_kind="+1")
        a = run_oracle(exp, 2e-3, 300; decoder="relay", seed=11, with_gap=true)
        b = run_oracle(exp, 2e-3, 300; decoder="relay", seed=11, with_gap=false)
        @test sum(a.outcomes[:, 1]) > 0  # the comparison is not vacuous
        @test a.outcomes == b.outcomes
    end
end

@testitem "day 3: the surface Y measurement has an exact gap" tags = [:properties] begin
    using FTShadows.Circuits: surface_local_experiment
    exp, lab = surface_local_experiment(3, "YZ", [("YI", 0), ("IZ", 0)])
    r = run_oracle(exp, 1e-3, 2000; seed=5, setting_labels=lab)
    @test Set(r.confidence_kind) == Set(["exact"])
    @test all(r.confidence .>= 0)
end

@testitem "day 3: logical channel and known answers" tags = [:properties] begin
    using FTShadows: NOV20_STATES, flipped_state, nov20_expectations, product_strings
    @test flipped_state(NOV20_STATES["bell"], 0) == [("XX", 1), ("ZZ", 0)]
    sids = product_strings("XYZ", 2)
    ds = Dict{String,Any}("code" => "surface", "d" => 3, "state" => "bell", "m" => 2, "p" => 0.0, "r" => 0,
                          "pL" => 0.05, "pL_qubit" => 0)
    rows, known = nov20_expectations(ds, sids, [])
    @test known["fidelity_ideal"] == 1.0
    @test abs(known["fidelity_with_logical_noise"] - 0.95) < 1e-12
    by = Dict(r["pauli"] => r for r in rows)
    @test by["YY"]["exact"] == -1 && abs(by["YY"]["exact_with_logical_noise"] + 0.9) < 1e-12
    @test by["XX"]["lambda_P"] ≈ 1 / 9 && by["XI"]["lambda_P"] ≈ 1 / 3
end

@testitem "day 3: the handover example runs clean" tags = [:properties, :handover] begin
    using FTShadows: generate_nov20, csv_rows, load_dataset
    using JSON
    # generate a dataset with generate_nov20, then run the example that estimates the fidelity.
    # Skipped when handover/ (not part of the repository) is absent.
    repo = pkgdir(FTShadows)
    ex = joinpath(repo, "handover", "examples", "bell_fidelity.jl")
    cal = joinpath(repo, "handover", "calibration_local.csv")
    if !(isfile(ex) && isfile(cal))
        @info "skipping: handover/examples not found"
    else
        spec = Dict{String,Any}("name" => "example_surface_d3_bell_random", "code" => "surface", "d" => 3,
                                "state" => "bell", "m" => 2, "p" => 1e-3, "r" => 3, "kind" => "random",
                                "shots" => 2000, "pL" => 0.0)
        mktempdir() do out
            generate_nov20(spec, out; calibration=csv_rows(cal))
            meta = JSON.parsefile(joinpath(out, "example_surface_d3_bell_random.meta.json"))
            @test meta["format_version"] == "1.1" && meta["draws"] == 2000
            julia = joinpath(Sys.BINDIR, Base.julia_exename())
            cmd = Cmd(`$julia --project=. ../examples/bell_fidelity.jl $(joinpath(out, "example_surface_d3_bell_random"))`;
                      dir=joinpath(repo, "handover", "julia"))
            txt = read(cmd, String)  # throws if the example exits with an error
            F = parse(Float64, split(split(txt, "F = ")[2])[1])
            se = parse(Float64, split(split(txt, "±")[2])[1])
            @test abs(F - meta["known_answers"]["fidelity_expected_noisy"]) < 4se
        end
    end
end

@testitem "surface: gap decoder consistent with PyMatching (r = 3)" tags = [:properties, :bridge] begin
    using PythonCall, Statistics
    using FTShadows.Circuits: multi_patch_surgery, build, pycircuit
    for (kind, n) in (("ZZ", 2), ("XX", 2), ("ZZ", 3))
        c = pycircuit(build(multi_patch_surgery(3, n, kind, [], 3), 5e-3; eigen_targets=true).circuit)
        dec = MatchingGapDecoder(c.detector_error_model(; decompose_errors=true))
        s = c.compile_detector_sampler(; seed=2).sample(3000; separate_observables=true)
        obs = pyconvert(Matrix{Bool}, s[1])
        pred, gap, _ = decode(dec, pyconvert(Matrix{Bool}, s[0]))
        @test dec.max_weight_mismatch < 1e-6
        @test dec.pred_mismatches == 0
        @test all(gap .>= 0)
        # higher confidence, lower error: errors sit at small gaps
        wrong = vec(any(pred .!= obs; dims=2))
        count(wrong) > 20 && @test median(gap[wrong]) < median(gap[.!wrong])
    end
end

@testitem "surface: memory experiment matches stim's generated circuit" tags = [:properties, :bridge] begin
    using PythonCall, Statistics
    using FTShadows.Circuits: single_patch_readout
    using FTShadows: pymodule
    # same channels as stim's rotated memory circuit ("stim_generated" noise); agreement within a factor 2
    d, r, p, shots = 3, 3, 3e-3, 200_000
    res = run_oracle(single_patch_readout(3, 'Z', [("Z", 0)], r), p, shots; seed=5, noise="stim_generated")
    ours = mean(res.outcomes)
    stim = pymodule(:stim)
    ref = stim.Circuit.generated("surface_code:rotated_memory_z"; distance=d, rounds=r, after_clifford_depolarization=p,
        before_measure_flip_probability=p, after_reset_flip_probability=p, before_round_data_depolarization=p)
    s = ref.compile_detector_sampler(; seed=5).sample(shots; separate_observables=true)
    m = pyimport("pymatching").Matching.from_detector_error_model(ref.detector_error_model(; decompose_errors=true))
    theirs = mean(any(pyconvert(Matrix{UInt8}, m.decode_batch(s[0])) .!= pyconvert(Matrix{Bool}, s[1]); dims=2))
    @test ours > 0 && theirs > 0
    @test 0.5 < ours / theirs < 2
end

@testitem "audit: gap respects parallel edges with different observables (r = 2)" tags = [:properties, :bridge] setup = [Props] begin
    using PythonCall
    using FTShadows.Circuits: multi_patch_surgery, build
    using FTShadows: detector_error_model, to_numpy_u8
    # p = 1e-3: for every single-detector syndrome, U_o = lightest correction of <= 2 graphlike mechanisms with observable
    # parity o; then w_best + gap <= max(U_0, U_1)
    c = build(multi_patch_surgery(3, 2, "ZZ", [], 2, "conjugate"), 1e-3; eigen_targets=true).circuit
    dem = detector_error_model(FTShadows.Circuits.stim_text(c); decompose_errors=true)
    edges = Props.weighted_edges(dem)
    U = Dict{Set{Int},Vector{Float64}}()
    for i in eachindex(edges), j in [[nothing]; collect(i+1:length(edges))]
        comb = j === nothing ? (edges[i],) : (edges[i], edges[j])
        s, o, w = Set{Int}(), UInt64(0), 0.0
        for (u, v, om, wt) in comb
            symdiff!(s, v == -1 ? Set([u]) : Set([u, v]))
            o ⊻= om
            w += wt
        end
        length(s) == 1 || continue
        d = get!(U, s, [Inf, Inf])
        d[o+1] = min(d[o+1], w)
    end
    dec = MatchingGapDecoder(dem)
    nd = pyconvert(Int, dem.num_detectors)
    worst = Float64[]
    for (s, (u0, u1)) in U
        (isfinite(u0) && isfinite(u1)) || continue
        det = zeros(UInt8, 1, nd)
        det[1, first(s)+1] = 1
        wb = pyconvert(Vector{Float64}, dec.base.decode_batch(to_numpy_u8(det); return_weights=true)[1])[1]
        push!(worst, wb + decode(dec, det)[2][1] - max(u0, u1))
    end
    @test maximum(worst) <= 1e-6
end

@testitem "audit: gap is calibrated for r = 2, conjugate" tags = [:properties] begin
    using Statistics
    using FTShadows.Circuits: multi_patch_surgery
    # confidence >= 6 means an error rate
    # well below 1%
    res = run_oracle(multi_patch_surgery(3, 2, "ZZ", [("ZI", 0), ("IZ", 0)], 2), 2e-3, 60000; seed=9)
    hi = res.confidence .>= 6
    @test count(hi) > 1000
    @test mean(res.outcomes[hi, 1] .!= 0) < 0.01
end

@testitem "audit: gap accounts for undetectable observable flips (r = 1)" tags = [:properties, :bridge] begin
    using PythonCall
    using FTShadows.Circuits: multi_patch_surgery, build, stim_text
    using FTShadows: detector_error_model, read_dem, llr
    # p = 1e-3
    c = build(multi_patch_surgery(3, 2, "ZZ", [], 1), 1e-3; eigen_targets=true).circuit
    dem = detector_error_model(stim_text(c); decompose_errors=true)
    pure = [e.p for e in read_dem(dem).errors if all(isempty(comp[1]) for comp in e.components)]
    @test !isempty(pure)  # a detector-free mechanism that flips the observable
    bound = minimum(llr.(pure))
    dec = MatchingGapDecoder(dem)
    _, gap, _ = decode(dec, zeros(UInt8, 1, pyconvert(Int, dem.num_detectors)))
    # r = 1 has an odd bulk cycle: no exact gap (NaN), never overconfident
    @test isnan(gap[1]) || gap[1] <= bound + 1e-6
    @test dec.undetectable_w[1] <= bound + 1e-9
end

@testitem "audit: XX surgery at r = 1 with the gap does not crash" tags = [:properties] begin
    using FTShadows.Circuits: multi_patch_surgery
    for n in (2, 3)
        state = n == 2 ? [("ZZ", 0), ("XX", 0)] : [("XXI", 0), ("IXX", 0), ("ZZZ", 0)]
        c = run_oracle(multi_patch_surgery(3, n, "XX", state, 1), 1e-3, 500; seed=1).confidence
        @test all(isnan.(c) .| (c .>= 0))  # NaN by policy at r = 1
    end
end

@testitem "BB72: Relay-BP gap decoder beats raw readout" tags = [:properties, :slow] begin
    using Statistics
    using FTShadows.Circuits: bb_gauging
    res = run_oracle(bb_gauging('X', 3; state_kind="+1"), 2e-3, 600; decoder="relay", seed=3)
    @test mean(res.outcomes) < mean(res.raw)
    @test all(res.confidence .>= 0)
end

@testitem "audit: BP-OSD multi-observable confidence is non-negative" tags = [:properties, :slow] begin
    using FTShadows.Circuits: bb_gauging
    exp = bb_gauging('X', 1; state_kind="+1", extra_targets=true)
    res = run_oracle(exp, 4e-3, 40; decoder="bposd", seed=5, decoder_kwargs=(osd_order=2, max_iter=30))
    @test size(res.outcomes, 2) > 1
    @test all(res.confidence .>= 0)
end
