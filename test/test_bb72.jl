# Library pieces behind the BB72 study scripts, against the stored values in test/reference/bb72/bb72.json.

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
    # the stored Gauss permutations are numpy draws
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

