# Codes, check lists, logical bases, gauging graphs and distances against the stored values in
# test/reference/circuits/codes.json.

@testitem "Python iteration orders are replayed" tags = [:circuits, :bridge] begin
    using PythonCall
    using FTShadows.Circuits: PyRandom, pyshuffle!, pyset_order
    pyrandom = pyimport("random")
    for N in [1, 2, 3, 7, 8, 63, 64, 65, 100, 160, 333]
        r = PyRandom(0)
        py = pyrandom.Random(0)
        for _ in 1:3
            x = pylist(0:N-1)
            py.shuffle(x)
            @test pyshuffle!(r, collect(0:N-1)) == pyconvert(Vector{Int}, x)
        end
    end
    for seed in (1, 12345, 2^40 + 7)
        x = pylist(0:49)
        pyrandom.Random(seed).shuffle(x)
        @test pyshuffle!(PyRandom(seed), collect(0:49)) == pyconvert(Vector{Int}, x)
    end
    tolist = pyeval("lambda xs: list(set(xs))", Main)
    rng = PyRandom(7)
    for trial in 1:400
        len = 1 + trial % 40
        xs = [pyshuffle!(rng, collect(0:199))[i] for i in 1:len]
        @test pyset_order(xs) == pyconvert(Vector{Int}, tolist(pylist(xs)))
    end
end

@testitem "GF(2) helpers" tags = [:circuits] begin
    using FTShadows.Circuits: rref, gf2rank, nullspace_gf2, gf2_solve, logical_basis, symp_inner, css_stabilizers
    M = Bool[1 1 0 1; 0 1 1 0; 1 0 1 1]
    R, piv = rref(M)
    @test piv == [1, 2] && gf2rank(M) == 2
    N = nullspace_gf2(M)
    @test size(N, 1) == 2 && all(iszero, isodd.(Int.(M) * Int.(permutedims(N))))
    x = gf2_solve(M, Bool[1, 1, 0])
    @test isodd.(Int.(M) * Int.(x)) == [1, 1, 0]
    @test gf2_solve(Bool[1 1; 1 1], Bool[1, 0]) === nothing
    S = css_stabilizers(Bool[1 1 1 1], Bool[1 1 1 1])  # [[4,2,2]]
    pairs = logical_basis(S)
    @test length(pairs) == 2
    for (i, (a, _)) in enumerate(pairs), (j, (_, b)) in enumerate(pairs)
        @test symp_inner(a, b) == (i == j)
    end
end

@testitem "surface patches and merged codes (d = 3, 5, 7)" tags = [:circuits] setup = [CircRef] begin
    using .CircRef
    using FTShadows.Circuits: Layout, add_rect!, rect_checks, patch_logicals, PauliAlgebra, from_dict, anticommutes,
                              dict_to_vec, gf2rank
    ref = circ_codes()
    chk(c) = Dict("pauli" => [[q, string(s)] for (q, s) in c.pauli],
                  "layer_of" => c.layer_of === nothing ? nothing : [[q, L] for (q, L) in c.layer_of],
                  "group" => c.group, "name" => c.name, "order" => c.order)
    pdl(p) = [[q, string(s)] for (q, s) in p]
    for d in (3, 5, 7)
        L = Layout()
        add_rect!(L, (0, d - 1), (0, d - 1))
        checks = rect_checks((0, d - 1), (0, d - 1), L.index)
        @test chk.(checks) == ref["surface_d$d"]
        n = length(L.index)
        alg = PauliAlgebra(n)
        vecs = [from_dict(alg, c.pauli) for c in checks]
        @test !any(anticommutes(alg, a, b) for a in vecs, b in vecs)
        @test gf2rank(permutedims(reduce(hcat, [dict_to_vec(c.pauli, n) for c in checks]))) == n - 1  # k = 1
        X, Z = patch_logicals((0, d - 1), (0, d - 1), L.index)
        xv, zv = from_dict(alg, X), from_dict(alg, Z)
        @test anticommutes(alg, xv, zv)
        @test !any(anticommutes(alg, xv, v) || anticommutes(alg, zv, v) for v in vecs)
        for kind in ("ZZ", "XX"), np in (2, 3)
            e = multi_patch_surgery(d, np, kind, [], 1)
            r = ref["surface_merged_$(kind)_$(np)p_d$d"]
            @test pdl.(e.prep_checks) == r["prep"]
            @test chk.(e.epochs[1].checks) == r["merged"]
            @test [[pdl(x), pdl(z)] for (x, z) in e.logical_basis] == r["basis"]
            @test pdl.(e.targets) == r["targets"]
            @test [[q, xy] for (q, xy) in sort!(collect(e.coords); by=first)] == r["coords"]
        end
    end
end

@testitem "BB72 and gross codes equal QuantumClifford's" tags = [:circuits] setup = [CircRef] begin
    using .CircRef
    using QuantumClifford, QuantumClifford.ECC
    using FTShadows.Circuits: BBCode, gf2rank
    ref = circ_codes()
    bits(M) = [join(Int.(M[i, :])) for i in axes(M, 1)]
    for (name, l) in (("bb72", 6), ("gross", 12))
        c = BBCode(name)
        @test bits(c.HX) == ref[name]["HX"] && bits(c.HZ) == ref[name]["HZ"]
        @test c.k == ref[name]["k"] == 12
        @test !any(isodd.(Int.(c.HX) * Int.(permutedims(c.HZ))))
        qc = BivariateBicycleViaCirculantMat(l, 6, [(:x, 3), (:y, 1), (:y, 2)], [(:y, 3), (:x, 1), (:x, 2)])
        @test Bool.(parity_matrix_x(qc)) == c.HX
        @test Bool.(parity_matrix_z(qc)) == c.HZ
        @test code_n(qc) == c.n && code_k(qc) == c.k
        chks = ordered_checks(c)
        @test length(chks) == 2 * l * 6
        @test [[[q, string(s)] for (q, s) in x.pauli] for x in chks] == [r["pauli"] for r in ref[name]["checks"]]
        @test [x.order for x in chks] == [r["order"] for r in ref[name]["checks"]]
        @test [(x.group, x.name) for x in chks] == [(r["group"], r["name"]) for r in ref[name]["checks"]]
    end
    c = BBCode("bb72")
    @test gf2rank(c.HX) == 30 && gf2rank(c.HZ) == 30
end

@testitem "BB72 oracle and local bases" tags = [:circuits] setup = [CircRef] begin
    using .CircRef
    using FTShadows.Circuits: bb_setup, local_basis, compile_setting, extra_qubits_of, PauliAlgebra, from_dict,
                              anticommutes, check_dicts, pauli_product, symp_inner, dict_to_vec, RowSpace, addrow!
    ref = circ_codes()
    pdl(p) = [[q, string(s)] for (q, s) in p]
    s = bb_setup()
    @test [[pdl(x), pdl(z)] for (x, z) in s.basis] == ref["bb72_oracle_basis"]
    @test [[pdl(x), pdl(z)] for (x, z) in local_basis()] == ref["bb72_local_basis"]
    @test Dict(string(P) => pdl(s.L[P]) for P in "XYZ") == ref["bb72_L"]
    # the bases handed to the estimator side
    for (file, B) in (("bb72_oracle_basis.json", s.basis), ("bb72_local_basis.json", local_basis()))
        h = JSON.parsefile(joinpath(CIRC_REF, "bases", file))
        @test h["xbar_support"] == [sort(collect(keys(x))) for (x, _) in B]
        @test h["zbar_support"] == [sort(collect(keys(z))) for (_, z) in B]
        @test all(all(==('X'), values(x)) && all(==('Z'), values(z)) for (x, z) in B)
    end
    # the chosen BB72 logicals
    @test s.weights['X'] == 6 && s.weights['Z'] == 6
    code = s.code
    alg = PauliAlgebra(code.n)
    stabs = [from_dict(alg, c) for c in check_dicts(code)]
    for P in "XYZ"
        L = from_dict(alg, s.L[P])
        @test !any(anticommutes(alg, L, v) for v in stabs)
        for (i, (Xi, Zi)) in enumerate(s.bases[P]), (j, (Xj, Zj)) in enumerate(s.bases[P])
            @test anticommutes(alg, from_dict(alg, Xi), from_dict(alg, Zj)) == (i == j)
            @test !anticommutes(alg, from_dict(alg, Xi), from_dict(alg, Xj))
            @test !anticommutes(alg, from_dict(alg, Zi), from_dict(alg, Zj))
        end
    end
    # audit AUDIT-1: the Y setting acts on the qubit measured by the X and Z settings
    yxz = pauli_product(code.n, [s.L['X'], s.L['Z']])
    @test s.L['Y'] == yxz && mod(length(intersect(keys(s.L['X']), keys(s.L['Z']))), 4) == 1
    # the BB72 local basis
    B = local_basis()
    @test length(B) == 12 && B[1] == s.Ypair
    for i in 1:3
        X, Z = B[i]
        @test length(X) == length(Z) == 6 && length(intersect(keys(X), keys(Z))) == 1
        for j in 1:i-1
            @test isempty(intersect(union(keys(X), keys(Z)), union(keys(B[j][1]), keys(B[j][2]))))
        end
    end
    for (i, (X, _)) in enumerate(B), (j, (_, Z)) in enumerate(B)
        @test symp_inner(dict_to_vec(X, 72), dict_to_vec(Z, 72)) == (i == j)
    end
    for (X, Z) in B, v in (X, Z)
        @test !any(symp_inner(dict_to_vec(v, 72), code.S[r, :]) for r in axes(code.S, 1))
    end
    # the BB72 compile rule and the full compile table
    @test compile_setting("XX") == (readout='X', gauged=Tuple{Int,Char}[], extra_qubits_estimate=0)
    @test compile_setting("XZ").gauged == [(0, 'X')] && compile_setting("XZ").readout == 'Z'
    @test compile_setting("XY") == (readout='X', gauged=[(1, 'Y')], extra_qubits_estimate=41)
    @test compile_setting("YY").gauged == [(0, 'Y'), (1, 'Y')]
    for (b, r) in ref["bb72_compile"]
        c = compile_setting(b)
        @test string(c.readout) == r["readout"] && c.extra_qubits_estimate == r["extra_qubits_estimate"]
        @test [[i, string(P)] for (i, P) in c.gauged] == r["gauged"]
    end
    @test Dict("$P$i" => extra_qubits_of(i, P) for i in 0:2 for P in "XYZ") == ref["bb72_extra_qubits"]
end

@testitem "gauging graphs and deformed codes" tags = [:circuits] setup = [CircRef] begin
    using .CircRef
    using FTShadows.Circuits: gauge_spec, cheeger, cycle_basis_edges, group_checks, PauliAlgebra, from_dict,
                              anticommutes, SMALL_GRAPHS, xor, zero_pvec
    ref = circ_codes()
    pdl(p) = [[q, string(s)] for (q, s) in p]
    graph(G) = Dict("nodes" => G.nodes, "adj" => [[u, G.adj[u]] for u in G.nodes])
    for i in 0:2, P in "XYZ"
        L, G, prs, os_ = gauge_spec(i, P)
        r = ref["bb72_gauge_specs"]["$P$i"]
        @test pdl(L) == r["L"]
        @test graph(G) == r["graph"]  # adjacency order included
        @test prs == r["pairings"]
        @test os_ == r["order_seed"]
        @test cheeger(G) == r["cheeger"]
        @test [[collect(e) for e in c] for c in cycle_basis_edges(G)] == r["cycle_basis"]
    end
    chk(c) = Dict("pauli" => pdl(c.pauli), "layer_of" => nothing, "group" => c.group, "name" => c.name, "order" => c.order)
    for (key, kw) in (("X", (;)), ("Z", (;)), ("Y", (;)), ("Y_minimal", (use_defaults=false,)),
                      ("X_small", SMALL_GRAPHS['X']), ("Z_small", SMALL_GRAPHS['Z']))
        e = bb_gauging(key[1], 2; state_kind="+1", extra_targets=true, kw...)
        r = ref["bb72_deformed"][key]
        ch = e.epochs[1].checks
        @test e.n_data == r["n"]
        @test chk.(ch) == r["checks"]
        @test group_checks(ch) == r["groups"]
        @test graph(e.graph) == r["graph"]
        @test all(e.graph_info[k] == v for (k, v) in r["graph_info"])
        @test pdl.(e.targets) == r["targets"]
        n = e.n_data
        S = permutedims(reduce(hcat, [dict_to_vec(c.pauli, n) for c in ch]))
        pairs = logical_basis(Matrix{Bool}(S))
        @test length(pairs) == r["k"] == 11
        @test [[pdl(vec_to_dict(a)), pdl(vec_to_dict(b))] for (a, b) in pairs] == r["logical_basis"]
        # the deformed code commutes and measures L
        alg = PauliAlgebra(n)
        v = [from_dict(alg, c.pauli) for c in ch]
        @test !any(anticommutes(alg, a, b) for a in v, b in v)
        prod = foldl(xor, (from_dict(alg, c.pauli) for c in ch if startswith(c.name, "A_")); init=zero_pvec(alg))
        @test prod == from_dict(alg, e.targets[1])
    end
    for (bases, r) in ref["bb72_local_experiments"]
        e, labels, comp = local_experiment(bases, [], 2)
        @test labels == r["labels"]
        @test pdl.(e.targets) == r["targets"]
        @test chk.(isempty(e.epochs) ? FTShadows.Circuits.Check[] : e.epochs[1].checks) == r["checks"]
        @test all(e.graph_info[k] == v for (k, v) in r["graph_info"])
    end
end

@testitem "code distances" tags = [:circuits, :slow, :bridge] setup = [CircRef] begin
    using .CircRef
    using QuantumClifford, QuantumClifford.ECC
    import JuMP, HiGHS
    using FTShadows.Circuits: BBCode, Layout, add_rect!, rect_checks, code_distance_stim, SMALL_GRAPHS

    # Exact distance of any stabilizer code: the lightest Pauli that commutes with every check and
    # anticommutes with one basis logical, minimised over the basis (one MIP per logical).
    function mip_distance(checks, n)
        S = permutedims(reduce(hcat, [dict_to_vec(c, n) for c in checks]))
        best = typemax(Int)
        for (a, b) in logical_basis(Matrix{Bool}(S)), L in (a, b)
            m = JuMP.Model(HiGHS.Optimizer)
            JuMP.set_silent(m)
            JuMP.@variable(m, x[1:n], Bin)
            JuMP.@variable(m, z[1:n], Bin)
            JuMP.@variable(m, w[1:n], Bin)
            JuMP.@constraint(m, [i in 1:n], w[i] >= x[i])
            JuMP.@constraint(m, [i in 1:n], w[i] >= z[i])
            sp(v) = JuMP.@expression(m, sum(x[i] for i in 1:n if v[n+i]) + sum(z[i] for i in 1:n if v[i]))
            for r in axes(S, 1)
                t = JuMP.@variable(m, integer = true, lower_bound = 0)
                JuMP.@constraint(m, sp(S[r, :]) == 2t)
            end
            t = JuMP.@variable(m, integer = true, lower_bound = 0)
            JuMP.@constraint(m, sp(L) == 2t + 1)
            JuMP.@objective(m, Min, sum(w))
            JuMP.optimize!(m)
            best = min(best, round(Int, JuMP.objective_value(m)))
        end
        return best
    end
    css(xs, zs, n) = CSS(Bool[haskey(c, q) for c in xs, q in 0:n-1], Bool[haskey(c, q) for c in zs, q in 0:n-1])
    mipalg = DistanceMIPAlgorithm(solver=HiGHS)

    for d in (3, 5, 7)
        L = Layout()
        add_rect!(L, (0, d - 1), (0, d - 1))
        ch = [c.pauli for c in rect_checks((0, d - 1), (0, d - 1), L.index)]
        n = d^2
        q = css([c for c in ch if 'X' in values(c)], [c for c in ch if 'Z' in values(c)], n)
        @test code_k(q) == 1
        @test distance(q, mipalg) == d
    end
    bb = BBCode("bb72")
    qbb = BivariateBicycleViaCirculantMat(6, 6, [(:x, 3), (:y, 1), (:y, 2)], [(:y, 3), (:x, 1), (:x, 2)])
    @test distance(qbb, mipalg) == 6
    # deformed gauging codes: k = 11, d = 6 (QuantumClifford's MIP for the CSS ones, our MIP for all,
    # and stim's search)
    for (P, kw) in (('X', (;)), ('Z', (;)), ('Y', (;)), ('X', SMALL_GRAPHS['X']), ('Z', SMALL_GRAPHS['Z']))
        e = bb_gauging(P, 2; kw...)
        ch = [c.pauli for c in e.epochs[1].checks]
        n = e.n_data
        # Z gauging has Gauss checks Z_v X_e, so only the X-gauging codes are CSS
        if all(c -> length(unique(values(c))) == 1, ch)
            q = css([c for c in ch if all(==('X'), values(c))], [c for c in ch if all(==('Z'), values(c))], n)
            @test code_k(q) == 11
            @test distance(q, mipalg) == 6
        end
        @test mip_distance(ch, n) == 6
        @test code_distance_stim(n, ch) == (6, 11)
    end
end

@testitem "logical labels: phase-tracked and phase-free maps agree" tags = [:circuits] setup = [CircRef] begin
    using .CircRef
    using FTShadows.Circuits: logical_to_physical, parse_state
    # logical-to-physical maps agree with the surface parser, on every label
    e, _ = setting_experiment("local:ZZZ", 3, 3, [], 0)
    for lab in vec(join.(Iterators.product(fill("IXYZ", 3)...)))
        lab == "III" && continue
        d1, s1 = logical_to_physical(lab, e.logical_basis, e.n_data)
        ((d2, _),) = parse_state([(lab, 0)], e.logical_basis, e.n_data)
        @test d1 == d2 && s1 == 0
    end
    @test exact_expectation([("ZZ", 0), ("XX", 0)], "YY") == -1
    @test exact_expectation([("XX", 0)], "ZZ") == 0
    @test_throws ArgumentError logical_to_physical("XQ", e.logical_basis[1:2], e.n_data)
end
