# Circuits: (a) every stored configuration builds to the stored circuit text and build record;
# (b) noiseless circuits give the exact logical outcomes; (c) circuit distances. Decoding is tested with the oracle.

@testitem "IR: stim text, fusion and parsing" tags = [:circuits, :bridge] setup = [CircRef] begin
    using .CircRef, PythonCall
    using FTShadows.Circuits: Circuit, parse_stim, num_measurements, rec, qubit
    c = Circuit()
    append!(c, "CX", [0, 1])
    append!(c, "CX", [2, 3])
    append!(c, "DETECTOR", [rec(-1)], [0, 1])
    append!(c, "DETECTOR", [rec(-2)], [0, 1])
    append!(c, "X_ERROR", [3], [1e-3 / 3])
    append!(c, "TICK")
    append!(c, "TICK")
    append!(c, "QUBIT_COORDS", [1], [1.5, 2.0])
    @test stim_text(c) == "CX 0 1 2 3\nDETECTOR(0, 1) rec[-1]\nDETECTOR(0, 1) rec[-2]\nX_ERROR(0.000333333) 3\nTICK\nTICK\nQUBIT_COORDS(1.5, 2) 1"
    for f in readdir(joinpath(CIRC_REF, "texts"); join=true)
        text = read(f, String)
        p = parse_stim(text)
        @test stim_text(p) == text == normalized_text(text)
        @test num_measurements(p) == pyconvert(Int, pycircuit(p).num_measurements)
    end
end

@testitem "tracker rules" tags = [:circuits] begin
    using FTShadows.Circuits: Tracker, PauliAlgebra, from_dict, measure!, query, reset!, bits_of, Bits
    alg = PauliAlgebra(3)
    t = Tracker(3)
    zz01, zz12, zz02 = from_dict(alg, ((0, 'Z'), (1, 'Z'))), from_dict(alg, ((1, 'Z'), (2, 'Z'))), from_dict(alg, ((0, 'Z'), (2, 'Z')))
    @test measure!(t, zz01, 0) === nothing  # unknown state: random
    @test measure!(t, zz12, 1) === nothing
    @test measure!(t, zz02, 2) == (bits_of([0, 1]), Bits())  # product of the first two
    @test measure!(t, from_dict(alg, ((0, 'X'),)), 3) === nothing  # anticommutes
    @test query(t, zz01) === nothing
    @test query(t, zz12) !== nothing
    alg = PauliAlgebra(2)
    t = Tracker(2)
    measure!(t, from_dict(alg, ((0, 'Z'), (1, 'Z'))), 0)
    reset!(t, 0, 'Z')
    @test query(t, from_dict(alg, ((0, 'Z'),))) == (Bits(), Bits())  # known constant
    q = query(t, from_dict(alg, ((1, 'Z'),)))
    @test q === nothing || !iszero(q[2])  # Z0Z1 information destroyed by the conditional flip
end

@testitem "surface surgery, readout and memory circuits match the stored text" tags = [:circuits] setup = [CircRef] begin
    using .CircRef
    fast = ("single_patch_readout", "multi_patch_surgery", "setting_experiment", "multi_patch_readout",
            "bravyi_memory", "bb72_readout", "multi_gauging")
    recs = [r for r in circ_configs() if r["spec"]["kind"] in fast]
    for rec in recs
        ok, msg = compare_config(rec; normalize=false)
        ok || @info rec["name"] msg
        @test ok
    end
    @test length(recs) == 548
end

@testitem "all stored circuit configurations build identically" tags = [:circuits, :slow, :bridge] setup = [CircRef] begin
    using .CircRef
    # every stored configuration: our text, normalized by stim, equals the stored text, and the build
    # records (records, detectors, targets, cost) agree
    counts = Dict{String,Int}()
    for rec in circ_configs()
        ok, msg = compare_config(rec)
        if ok && !haskey(rec, "error")
            c = first(build_config(rec["spec"]))
            ok = sha_hex(normalized_text(stim_text(c))) == rec["sha256"]
        end
        ok || @info rec["name"] msg
        @test ok
        k = rec["spec"]["kind"]
        counts[k] = get(counts, k, 0) + 1
    end
    @test sum(values(counts)) == 1772
    @test counts["surface_local"] == 759 && counts["bb_local"] == 303 && counts["bb_gauging"] == 144
end

@testitem "noiseless surface surgery and readout are exact" tags = [:circuits, :bridge] setup = [CircRef] begin
    using .CircRef, PythonCall
    using FTShadows.Circuits: logical_from_label
    # surface eigenstate cases
    cases = [
        (("ZZ", 2), [("ZI", 0), ("IZ", 0)], [0]), (("ZZ", 2), [("ZI", 0), ("IZ", 1)], [1]),
        (("ZZ", 2), [("ZI", 1), ("IZ", 1)], [0]), (("ZZ", 2), [("ZZ", 0), ("XX", 0)], [0]),
        (("ZZ", 2), [("ZZ", 1), ("XX", 0)], [1]), (("XX", 2), [("XI", 0), ("IX", 1)], [1]),
        (("XX", 2), [("ZZ", 0), ("XX", 0)], [0]), (("XX", 2), [("ZZ", 0), ("XX", 1)], [1]),
        (("ZZ", 3), [("ZZI", 0), ("IZZ", 1), ("XXX", 0)], [0, 1]), (("XX", 3), [("XXI", 1), ("IXX", 0), ("ZZZ", 0)], [1, 0]),
    ]
    for term in ("conjugate", "same"), ((kind, np), state, expected) in cases
        out = noiseless_outcomes(build(multi_patch_surgery(3, np, kind, state, 2, term), 0.0), 200; seed=7)
        @test all(out .== permutedims(expected))
    end
    for (basis, state, bitv) in (('Z', [("Z", 0)], 0), ('Z', [("Z", 1)], 1), ('X', [("X", 0)], 0), ('X', [("X", 1)], 1)),
        rounds in (0, 2)
        @test all(noiseless_outcomes(build(single_patch_readout(3, basis, state, rounds), 0.0), 100; seed=3) .== bitv)
    end
    fair(k, n) = abs(k - n / 2) / sqrt(n / 4) < 3.29  # two-sided p > 1e-3
    for (kind, state) in (("ZZ", [("XI", 0), ("IX", 0)]), ("XX", [("ZI", 0), ("IZ", 0)]), ("ZZ", [("XI", 0), ("IZ", 0)]))
        out = noiseless_outcomes(build(multi_patch_surgery(3, 2, kind, state, 2), 0.0), 4000; seed=11)
        @test fair(count(out), 4000)
    end
    # after ZZ on |++>, conjugate termination hides the ZZ outcome, same-basis termination repeats it
    for (term, corr) in (("conjugate", false), ("same", true))
        e = multi_patch_surgery(3, 2, "ZZ", [("XI", 0), ("IX", 0)], 2, term)
        b = build(e, 0.0)
        M = pyconvert(Matrix{Bool}, pycircuit(b.circuit).compile_sampler(; seed=5).sample(2000))
        zz = reduce(.⊻, (M[:, r+1] for r in b.target_records[1]))
        ro = M[:, end-e.n_data+1:end]
        P = term == "same" ? "ZZ" : "XX"
        qs = sort(collect(keys(logical_from_label(P, e.logical_basis, e.n_data))))
        term == "conjugate" && append!(qs, [q for (q, xy) in e.coords if xy == [3.0, 0.0]])
        final = reduce(.⊻, (ro[:, q+1] for q in sort(qs)))
        if corr
            @test final == zz
        else
            @test abs(sum(final .== zz) / 2000 - 0.5) < 0.06
            @test !any(final)
        end
    end
    # day 2: seeded random and product states, every subset product of every setting
    S = circ_states()
    for np in (2, 3), t in 0:7
        state = st(S[iseven(t) ? "stab:$np:$(1000np + t)" : "prod:$np:$(2000np + t)"])
        for setting in vcat(["surgery:ZZ", "surgery:XX"], S["day2_settings:$np"][t+1])
            e, labels = setting_experiment(setting, 3, np, state, 2)
            out = noiseless_outcomes(build(e, 0.0), 400; seed=t)
            bases = startswith(setting, "local") ? setting[7:end] : nothing
            if bases !== nothing
                @test check_products(state, bases, out) == 0
            else
                for (j, lab) in enumerate(labels)
                    ex = exact_expectation(state, lab)
                    ex == 0 ? (@test length(unique(out[:, j])) == 2) : (@test all(out[:, j] .== (ex == 1 ? 0 : 1)))
                end
            end
        end
    end
    # audit: unspecified logical qubits are maximally mixed; anticommuting generators are rejected
    e, _ = setting_experiment("local:ZX", 3, 2, [("IZ", 0)], 0)
    m = sum(noiseless_outcomes(build(e, 0.0), 2000; seed=3)[:, 1]) / 2000
    @test 0.4 < m < 0.6
    @test_throws ArgumentError build(multi_patch_surgery(3, 2, "ZZ", [("ZI", 0), ("XI", 0)], 2), 0.0)
end

@testitem "noiseless surface local settings (Gidney fragments) are exact" tags = [:circuits, :bridge] setup = [CircRef] begin
    using .CircRef
    # surface Y eigenstates, and determinism
    for d in (3, 5), sign in (0, 1)
        e, _ = surface_local_experiment(d, "Y", [("Y", sign)])
        pycircuit(build(e, 0.0; eigen_targets=true).circuit).detector_error_model()
        @test all(noiseless_outcomes(build(e, 0.0), 200; seed=sign)[:, 1] .== sign)
    end
    S = circ_states()
    for mem in (0, 2), t in 0:2, bases in ("YXZ", "YYY", "ZYX")
        state = st(S["stab:3:$(700 + t)"])
        e, _ = surface_local_experiment(3, bases, state; memory_rounds=mem)
        @test check_products(state, bases, noiseless_outcomes(build(e, 0.0), 300; seed=t)) == 0
    end
    # all Bell (9) and GHZ (27) settings, both signs, as in the Nov-20 datasets
    for (state, m) in ((BELL, 2), (flipped(BELL, 0), 2), (GHZ, 3), (flipped(GHZ, 0), 3)), mem in (0, 3)
        for bases in vec(join.(Iterators.product(fill("XYZ", m)...)))
            e, _ = surface_local_experiment(3, bases, state; memory_rounds=mem)
            @test check_products(state, bases, noiseless_outcomes(build(e, 0.0), 200; seed=1)) == 0
        end
    end
end

@testitem "noiseless surface local settings at d = 5" tags = [:circuits, :bridge, :slow] setup = [CircRef] begin
    using .CircRef
    for (state, m) in ((BELL, 2), (flipped(BELL, 0), 2), (GHZ, 3), (flipped(GHZ, 0), 3)), mem in (0, 3)
        for bases in vec(join.(Iterators.product(fill("XYZ", m)...)))
            e, _ = surface_local_experiment(5, bases, state; memory_rounds=mem)
            @test check_products(state, bases, noiseless_outcomes(build(e, 0.0), 200; seed=1)) == 0
        end
    end
    for sign in (0, 1)
        e, _ = surface_local_experiment(7, "Y", [("Y", sign)])
        @test all(noiseless_outcomes(build(e, 0.0), 200; seed=sign)[:, 1] .== sign)
    end
end

@testitem "noiseless BB72 gauging and readout are exact" tags = [:circuits, :bridge] setup = [CircRef] begin
    using .CircRef
    using FTShadows.Circuits: group_checks, pauli_product
    fair(k, n) = abs(k - n / 2) / sqrt(n / 4) < 3.29
    for P in "XZY"
        # noiseless outcomes exact, detectors deterministic
        for (sk, bitv) in (("+1", 0), ("-1", 1))
            @test all(noiseless_outcomes(build(bb_gauging(P, 2; state_kind=sk), 0.0), 100; seed=4) .== bitv)
        end
        out = noiseless_outcomes(build(bb_gauging(P, 2; state_kind="random"), 0.0), 2000; seed=5)
        @test fair(count(out), 2000)
        e = bb_gauging(P, 2; state_kind="+1")
        noiseless_outcomes(build(e, 0.0; eigen_targets=true), 100; seed=1)
        pycircuit(build(e, 1e-3; eigen_targets=true).circuit).detector_error_model()
        # audit AUDIT-2: the raw outcome is the product of the round-1 Gauss records only
        b = build(e, 0.0)
        rid, k = Dict{Int,Int}(), b.n_prep
        for g in group_checks(e.epochs[1].checks)
            foreach(((j, i),) -> rid[i] = k + j - 1, enumerate(g))
            k += length(g)
        end
        @test Set(b.target_records[1]) == Set(rid[i - 1] for (i, c) in enumerate(e.epochs[1].checks) if startswith(c.name, "A_"))
        # BB72 eigenstates of both signs
        S = circ_states()
        for sign in (0, 1)
            state = vcat([(P * "I"^11, sign)], [("I" * l, s) for (l, s) in st(S["stab:11:$(5 + sign)"])])
            @test all(noiseless_outcomes(build(bb_gauging(P, 1; state_labels=state), 0.0), 200; seed=sign)[:, 1] .== sign)
            if P != 'Y'
                er, _ = bb72_readout_experiment(P, state)
                @test all(noiseless_outcomes(build(er, 0.0), 200; seed=sign)[:, 1] .== sign)
            end
        end
    end
    # audit AUDIT-1: +1 eigenstate of the phase-free product of the measured X and Z logicals gives Y outcome 0
    s = bb_setup()
    e = bb_gauging('Y', 1; state_kind="+1")
    e.logical_state = [(pauli_product(72, [s.L['X'], s.L['Z']]), 0)]
    @test all(noiseless_outcomes(build(e, 0.0), 200; seed=3)[:, 1] .== 0)
    # BB72 local settings, noiseless
    S = circ_states()
    for bases in ("XY", "YY", "ZX", "XYZ")
        m = length(bases)
        state = st(S["stab:$m:$(900 + m)"])
        ex, labels, _ = local_experiment(bases, state, 1)
        @test labels == [join(j == i ? b : 'I' for j in 1:m) for (i, b) in enumerate(bases)]
        @test check_products(state, bases, noiseless_outcomes(build(ex, 0.0), 200; seed=1)) == 0
    end
end

@testitem "noiseless BB72: random states and all Bell/GHZ local settings" tags = [:circuits, :bridge, :slow] setup = [CircRef] begin
    using .CircRef
    S = circ_states()
    # BB72 random states, noiseless, match the exact values
    for t in 0:3
        state = st(S[iseven(t) ? "stab:12:$(77 + t)" : "prod:12:$(88 + t)"])
        for kind in "XZ"
            e, labels = bb72_readout_experiment(kind, state)
            out = noiseless_outcomes(build(e, 0.0), 400; seed=t)
            for (j, lab) in enumerate(labels)
                ex = exact_expectation(state, lab)
                ex == 0 ? (@test length(unique(out[:, j])) == 2) : (@test all(out[:, j] .== (ex == 1 ? 0 : 1)))
            end
        end
        for P in "XZY"
            out = noiseless_outcomes(build(bb_gauging(P, 1; state_labels=state), 0.0), 400; seed=t)
            ex = exact_expectation(state, P * "I"^11)
            ex == 0 ? (@test length(unique(out[:, 1])) == 2) : (@test all(out[:, 1] .== (ex == 1 ? 0 : 1)))
        end
    end
    # the Nov-20 BB72 settings (noiseless datasets use the real r = 3 circuits), both signs
    for (state, m) in ((BELL, 2), (flipped(BELL, 0), 2), (GHZ, 3), (flipped(GHZ, 0), 3))
        for bases in vec(join.(Iterators.product(fill("XYZ", m)...)))
            e, _, _ = local_experiment(bases, state, 3)
            @test check_products(state, bases, noiseless_outcomes(build(e, 0.0), 100; seed=1)) == 0
        end
    end
end

@testitem "plug-ins build the same circuits as direct calls" tags = [:circuits] setup = [CircRef] begin
    using .CircRef
    using FTShadows.Circuits: build_setting, REGISTRY, op, gauge_spec, gauging_experiment, multi_gauging, local_basis
    st2 = [("XX", 0), ("ZZ", 0)]
    e1, l1, d1 = build_setting("surface", "local:XY", st2; d=3, memory_rounds=3)
    e2, l2 = surface_local_experiment(3, "XY", st2; memory_rounds=3)
    @test stim_text(build(e1, 1e-3).circuit) == stim_text(build(e2, 1e-3).circuit) && l1 == l2 && d1 == "pymatching"
    e1, l1, d1 = build_setting("bb72", "local:YZ", st2; rounds=2)
    e2, l2, _ = local_experiment("YZ", st2, 2)
    @test stim_text(build(e1, 1e-3).circuit) == stim_text(build(e2, 1e-3).circuit) && l1 == l2 && d1 == "relay"
    e1, l1, _ = build_setting("bb72", "gauge:Y0", [("Y" * "I"^11, 0)]; rounds=2)
    @test l1 == ["Y" * "I"^11] && stim_text(build(e1, 0.0).circuit) ==
                                  stim_text(build(bb_gauging('Y', 2; state_labels=[("Y" * "I"^11, 0)]), 0.0).circuit)
    @test_throws ArgumentError build_setting("bb72", "gauge:X1", []; rounds=2)
    @test collect(keys(REGISTRY)) == ["surface", "bb72"]
    @test op(REGISTRY["bb72"], "gauge:Y0").kind == "gauging"
    @test_throws KeyError op(REGISTRY["surface"], "gauge:X0")
    # single-logical BB72 gauging builds the same circuit as the gauging experiment
    code = bb_setup().code
    for P in "XYZ"
        L, G, pr, os_ = gauge_spec(0, P)
        rb = Dict('X' => 'Z', 'Z' => 'X', 'Y' => 'Z')[P]
        a = gauging_experiment(code.n, ordered_checks(code), L, 2, rb; logical_basis=local_basis(), G, pairings=pr, order_seed=os_)
        b = multi_gauging([gauge_spec(0, P)], 2, rb; logical_basis=local_basis())
        @test stim_text(build(a, 1e-3).circuit) == stim_text(build(b, 1e-3).circuit)
    end
end

@testitem "circuit distances: surface code" tags = [:circuits, :bridge, :slow, :distance] setup = [CircRef] begin
    using .CircRef, PythonCall
    using FTShadows.Circuits: stim_search
    graphlike(c) = pyconvert(Int, pybuiltins.len(pycircuit(c).detector_error_model(; decompose_errors=true).shortest_graphlike_error()))
    # deterministic detectors and circuit distance, two and three patches
    for d in (3, 5), kind in ("ZZ", "XX"), term in ("conjugate", "same")
        e = multi_patch_surgery(d, 2, kind, [], d, term)
        @test !any(noiseless_outcomes(build(e, 0.0; eigen_targets=true), 50; seed=1))
        @test graphlike(build(e, 1e-3; eigen_targets=true).circuit) == d
    end
    for d in (3, 5)
        @test graphlike(build(multi_patch_surgery(d, 3, "ZZ", [], d), 1e-3; eigen_targets=true).circuit) == d
    end
    # the hypergraph search finds the same distance (d = 3)
    @test stim_search(build(multi_patch_surgery(3, 2, "ZZ", [], 3), 1e-3; eigen_targets=true).circuit; max_dets=4, max_degree=4) == 3
    # detectors stay local (d = 5, 7)
    for d in (5, 7), kind in ("ZZ", "XX"), term in ("conjugate", "same")
        c = build(multi_patch_surgery(d, 2, kind, [], 2, term), 1e-3; eigen_targets=true).circuit
        @test maximum(length(i.targets) for i in c.ops if i.name == "DETECTOR") <= 5
        pycircuit(c).detector_error_model(; decompose_errors=true)
    end
    @test graphlike(build(multi_patch_surgery(7, 2, "ZZ", [], 7), 1e-3; eigen_targets=true).circuit) == 7
    # Gidney Y measurement: distance d needs b = 1, 2, 4 boundary rounds at d = 3, 5, 7
    for (d, b, expect, dets) in ((3, 1, 3, 29), (5, 2, 5, 110), (5, 1, 4, nothing), (7, 3, 6, nothing), (7, 4, 7, 315))
        e, _ = surface_local_experiment(d, "Y", [("Y", 0)]; boundary_rounds=b)
        r = build(e, 1e-3; eigen_targets=true)
        dets === nothing || @test r.detector_count == dets
        @test stim_search(r.circuit; max_dets=6, max_degree=8) == expect
    end
end

@testitem "circuit distances: BB72 gauging (r = 6, measured logicals)" tags = [:circuits, :bridge, :slow, :distance] setup = [CircRef] begin
    using .CircRef
    using FTShadows.Circuits: stim_search
    # Each search takes minutes and several GB of memory, so only these three run here; the other
    # gauging circuits (spectators, local combinations) are covered by the stored-circuit identity above.
    for P in "XZY"
        @test stim_search(build(bb_gauging(P, 6; state_kind="+1"), 1e-3; eigen_targets=true).circuit) == 6
        GC.gc()
    end
    for bases in ("XZ", "YX")  # combinations Z:X0 and X:Y0
        e, _, _ = local_experiment(bases, [], 6)
        @test stim_search(build(e, 1e-3; eigen_targets=true).circuit) == 6
        GC.gc()
    end
end

@testitem "Gidney's generator reproduces his archived circuits (bridge)" tags = [:circuits, :bridge] begin
    using PythonCall
    using FTShadows.Circuits: VENDOR
    # Gidney's archived circuits (test/reference/snapshots/gidney_reference) are reproduced byte for byte
    sys = pyimport("sys")
    pyin(VENDOR, sys.path) || sys.path.insert(0, VENDOR)
    make_circuit = pyimport("midout._make_circuit").make_circuit
    gen = pyimport("midout.gen")
    dir = joinpath(pkgdir(FTShadows), "test", "reference", "snapshots", "gidney_reference")
    for d in (3, 5)
        c = make_circuit(; basis="Y", distance=d, noise=gen.NoiseModel.si1000(1e-3), boundary_rounds=d ÷ 2, memory_rounds=d)
        f = only(filter(n -> startswith(n, "r=$d,d=$d,p=0.001,noise=SI1000,b=Y,rb=$(d ÷ 2),q="), readdir(dir)))
        @test strip(pyconvert(String, pybuiltins.str(c))) == strip(read(joinpath(dir, f), String))
    end
end
