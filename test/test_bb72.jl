# The library parts behind the BB72 study scripts against reference values
# (test/reference/bb72/bb72.json, exported from the original implementation).

@testmodule BB72Ref begin
    using JSON, FTShadows
    const REF = JSON.parsefile(joinpath(pkgdir(FTShadows), "test", "reference", "bb72", "bb72.json"))
    support(d) = [[q, string(p)] for (q, p) in sort!(collect(d); by = first)]
end

@testitem "seeded numpy draws: random states, Gauss orders, Cheeger padding" tags = [:bb72, :circuits] setup = [BB72Ref] begin
    using FTShadows.Circuits: random_stabilizer_state, bb_setup, ordered_checks, build_graph, bb_gauging,
                              gauss_order_numpy, python_rng, sorted_edges
    for c in BB72Ref.REF["random_stabilizer_state"]
        @test random_stabilizer_state(c["k"], c["seed"]) == [(String(l), Int(s)) for (l, s) in c["gens"]]
    end
    # the exported Gauss permutations are numpy's
    for (key, perm) in python_rng()["gauss_perms"]
        seed, len = parse.(Int, split(key, ':'))
        @test gauss_order_numpy(seed, 0:len-1) == perm
    end
    s = bb_setup()
    G, _ = build_graph([c.pauli for c in ordered_checks(s.code)], s.L['Y']; target_cheeger = 1.5, max_degree = 4)
    @test [collect(e) for e in sorted_edges(G)] == BB72Ref.REF["cheeger_padding_Y_1.5"]
    for (seed, orders) in BB72Ref.REF["gauss_orders_Y_1.5"]
        e = bb_gauging('Y', 1; order_seed = parse(Int, seed), target_cheeger = 1.5)
        @test [c.order for c in e.epochs[1].checks if startswith(c.name, "A_")] == orders
    end
end

@testitem "ISD distance bound, CSS basis and memory circuits" tags = [:bb72, :circuits] setup = [BB72Ref] begin
    using SHA
    using FTShadows.Circuits: bb_gauging, code_distance_isd, css_basis, BBCode, memory_experiment, ordered_checks, build,
                              stim_text
    for c in BB72Ref.REF["isd"]
        kw = isempty(c["kw"]) ? (;) : (use_defaults = false,)
        e = bb_gauging(c["P"][1], 1; kw...)
        w, v = code_distance_isd(e.n_data, [ch.pauli for ch in e.epochs[1].checks]; iters = 20, seed = c["seed"])
        @test w == c["weight"]
        @test BB72Ref.support(v) == c["logical"]
    end
    code = BBCode("bb72")
    basis = css_basis(code)
    @test [[sort(collect(keys(x))), sort(collect(keys(z)))] for (x, z) in basis] == BB72Ref.REF["css_basis"]
    for (key, sha) in BB72Ref.REF["memory_circuit_sha256"]
        bk, perm = split(key, ':')
        exp = memory_experiment(code, ordered_checks(code; order_x = parse.(Int, collect(perm))), bk[1], 3, basis)
        @test bytes2hex(sha256(stim_text(build(exp, 1e-3; eigen_targets = true).circuit))) == sha
    end
end

@testitem "two-Xbar pair, 6-edge candidates, local combinations and AUC" tags = [:bb72] setup = [BB72Ref] begin
    F = FTShadows
    code, Ls = F.multi_pair()
    @test [sort(collect(keys(L))) for L in Ls] == BB72Ref.REF["multi_pair"]
    @test [sort(collect(keys(d))) for d in F.multi_spectators(code, Ls)] == BB72Ref.REF["multi_spectators"]
    for (P, ref) in BB72Ref.REF["small_graph_candidates"]
        c = F.small_graph_candidates(P[1])
        @test length(c) == ref["n"]
        @test [[collect(e) for e in c[1][2]], c[1][3]] == ref["first"]
        @test [[collect(e) for e in c[end][2]], c[end][3]] == ref["last"]
        @test Float64(sum(Rational{BigInt}(x[3]) for x in c)) == ref["cheeger_sum"]  # Python 3.12 sums floats exactly
    end
    @test collect(F.local_distance_combos([2, 3])) == [String(k) => String(v) for (k, v) in BB72Ref.REF["local_distance_combos"]]
    a = BB72Ref.REF["auc"]
    @test F.auc(Float64.(a["score"]), Bool.(a["err"])) === a["auc"]
    @test F.py_permutations(3) == [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]]
    @test F.py_combinations(collect(1:4), 2) == [[1, 2], [1, 3], [1, 4], [2, 3], [2, 4], [3, 4]]
end

@testitem "pandas-style CSV and numpy chunk files" tags = [:bb72] begin
    F = FTShadows
    using FTShadows: OrderedDict
    dir = mktempdir()
    rows = [OrderedDict{String,Any}("P" => "X", "q" => 0.1, "c" => NaN, "ok" => true, "g" => nothing),
            OrderedDict{String,Any}("P" => "Z", "q" => 1e-5, "c" => 0.5, "ok" => false, "g" => 1.0)]
    F.write_frame_csv(joinpath(dir, "a.csv"), rows)
    @test read(joinpath(dir, "a.csv"), String) == "P,q,c,ok,g\nX,0.1,,True,\nZ,1e-05,0.5,False,1.0\n"
    o = UInt8[0 1; 1 1; 0 0]
    F.npy_save(joinpath(dir, "o.npy"), o)
    @test F.npy_load_u8(joinpath(dir, "o.npy")) == o
    F.npz_save_compressed(joinpath(dir, "z.npz"); err = UInt8[0, 1], gap = [0.5, Inf], cost = OrderedDict("a" => 1, "b" => 2))
    z = F.npz_load(joinpath(dir, "z.npz"))
    @test z["err"] == UInt8[0, 1] && z["gap"] == [0.5, Inf] && collect(z["cost"]) == ["a" => 1, "b" => 2]
end

@testitem "BB72 calibration seeds are unique across configurations" tags = [:bb72] begin
    # Run the script's adaptive schedule (default, not --quick) with a fake pool that records the jobs and
    # returns no errors, so every configuration runs to the shot cap; no seed may belong to two configurations.
    F = FTShadows
    jobs = []
    fake_pmap(f, js) = (append!(jobs, js);
                        [(zeros(UInt8, j[5]), zeros(j[5]), F.OrderedDict("extra_qubits" => 0), F.OrderedDict("n" => 0)) for j in js])
    F.bb72_qr_adaptive(F.bb72_qr_configs("XZY"), 200, 6e4; chunk_dir = mktempdir(), pmap = fake_pmap)
    owners = Dict{Int,Set{Tuple}}()
    for (Pl, r, st, seed, _...) in jobs
        push!(get!(owners, seed, Set{Tuple}()), (Pl, r, st, false))
    end
    clash = filter(kv -> length(kv[2]) > 1, owners)
    @test isempty(clash)
    @test length(jobs) == 24 * 12  # 24 configurations, 12 chunks of 5000 to the 60000-shot cap
end
