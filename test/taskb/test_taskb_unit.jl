# Unit tests of Task B: GF(2) helpers, BB codes, ensembles, groups, shadow maps and earlier regressions. Items
# that need Oscar (group orders, generator reduction) or take more than a few seconds are tagged :slow.

@testmodule TaskBToys begin
    using FTShadows
    import FTShadows.TaskB as TB
    G(name, qs...) = TB.Gate(name, collect(qs))
    # the k = 3 toy groups; SWAP gates are allowed in gate lists that only go to QuantumClifford
    const TOY_GENS = Dict(
        "phase_type" => [[G(:H, 0), G(:H, 1), G(:H, 2)], [G(:SWAP, 0, 1), G(:SWAP, 1, 2)], [G(:CZ, 0, 1), G(:S, 2)]],
        "css_only" => [[G(:H, 0), G(:H, 1), G(:H, 2)], [G(:SWAP, 0, 1), G(:SWAP, 1, 2)], [G(:CX, 0, 1)]],
    )
    toy_symplectic(name) = [TB.tableau_symplectic(TB.circuit_clifford(c, 3)) for c in TOY_GENS[name]]
    sp_gens(k) = [TB.gates_to_symplectic(g, k) for g in
                  [[G(:H, 0)], [G(:S, 0)], (k > 1 ? [[G(:CX, 0, 1)], [TB.Gate(:SWAPPERM, [mod(i + 1, k) for i in 0:k-1])]] : [])...]]
end

@testitem "TaskB: GF(2) inverse and nullspace" begin
    using PythonCall
    import FTShadows.TaskB as TB
    # matrices from numpy's default_rng(0), through the bridge
    np = pyimport("numpy")
    rng = np.random.default_rng(0)
    draw(r, c) = Bool.(pyconvert(Matrix{UInt8}, rng.integers(0, 2, (r, c); dtype=np.uint8)))
    for _ in 1:20
        m = draw(8, 8)
        while TB.gf2rank(m) < 8
            m = draw(8, 8)
        end
        @test TB.gf2mul(m, TB.gf2inv(m)) == TB.GF2Mat(TB.I(8))
        a = draw(5, 9)
        ns = TB.nullspace_gf2(a)
        @test size(ns, 1) == 9 - TB.gf2rank(a)
        @test !any(TB.gf2mul(a, permutedims(ns)))
    end
end

@testitem "TaskB: BB parameters" begin
    import FTShadows.TaskB as TB
    for (name, n) in (("BB72", 72), ("gross", 144))
        _, _, _, hx, hz = TB.bb_matrices(TB.BB_PARAMS[name]...)
        @test size(hx, 2) == n
        @test !any(TB.gf2mul(hx, permutedims(hz)))
        @test n - TB.gf2rank(hx) - TB.gf2rank(hz) == 12
        # QuantumClifford's BB matrices use the convention x = S_l ⊗ I_m, S_l[i, i+1] = 1
        @test (hx, hz) == TB.bb_matrices_direct(TB.BB_PARAMS[name]...)[4:5]
    end
end

@testitem "TaskB: BB72 distance is 6 exactly" tags = [:slow] begin
    import FTShadows.TaskB as TB
    code = TB.make_bb_code("BB72", "generic")
    n = code.n
    hx, hz = code.meta[:hx], code.meta[:hz]
    @test TB.min_weight_logical_upto(hz, code.lz[:, n+1:end], 5) == (nothing, nothing)
    @test TB.min_weight_logical_upto(hx, code.lx[:, 1:n], 5) == (nothing, nothing)
    @test TB.min_weight_logical_upto(hz, code.lz[:, n+1:end], 6)[1] == 6
    @test TB.min_weight_logical_upto(hx, code.lx[:, 1:n], 6)[1] == 6
end

@testitem "TaskB: gross published basis" begin
    import FTShadows.TaskB as TB
    code = TB.make_bb_code("gross", "fgh")
    n = code.n
    @test TB.gf2mul(code.lx[:, 1:n], permutedims(code.lz[:, n+1:end])) == TB.I(12)
    @test Set(sum(code.lx; dims=2)) == Set(sum(code.lz; dims=2)) == Set([12])
end

@testitem "TaskB: La-cross [[65,9]]" begin
    c = FTShadows.TaskB.make_lacross_code(3, 7)
    @test (c.n, FTShadows.TaskB.logical_k(c)) == (65, 9)
end

@testitem "TaskB: automorphism counts and small group orders" begin
    import FTShadows.TaskB as TB
    # Sayginel et al. 2409.18175 Table VIII: |Aut| = 864 (BB72), 288 (gross)
    for (name, aut, orders) in (("BB72", 864, (36, 72, 864)), ("gross", 288, (36, 72, 144)))
        code = TB.make_bb_code(name, "fgh")
        fam, auts, duals, folds, info = TB.bb_families(code)
        @test info["n_perm_auts"] + info["n_zx_dualities"] == aut
        sh = [fam["shift_x"].logical, fam["shift_y"].logical]
        @test TB.closure(sh)[1] == orders[1]
        @test TB.closure([sh..., fam["zx_bravyi"].logical])[1] == orders[2]
        @test TB.closure([sh; [e.logical for e in TB.dedupe([auts; duals])]])[1] == orders[3]
    end
end

@testitem "TaskB: gross x^6 is the logical identity" begin
    import FTShadows.TaskB as TB
    code = TB.make_bb_code("gross", "fgh")
    @test TB.logical_action(code, TB.perm_op(TB.shift_perm(12, 6, 6, 0))) == TB.I(24)
end

@testitem "TaskB: elements are logical Cliffords and QuantumClifford agrees" begin
    using QuantumClifford
    import FTShadows.TaskB as TB
    code = TB.make_bb_code("BB72", "fgh")
    fam, auts, duals, folds, _ = TB.bb_families(code)
    n = code.n
    L = TB.logicals(code)
    for e in [fam["shift_x"], fam["shift_y"], fam["zx_bravyi"], auts[1:4]..., duals[1:4]..., folds[1:12:end]...]
        @test TB.is_symplectic(e.logical)
        @test TB.preserves_stabilizer_group(code, e.phys)
        @test TB.gates_to_symplectic(e.gates, n) == e.phys
        st = TB.apply_gates!(TB._rows_tableau(L, n), e.gates)
        imgs = reduce(vcat, [permutedims([xbit(st[i]); zbit(st[i])]) for i in axes(L, 1)])
        @test TB.logical_coords(code, imgs) == e.logical
    end
end

@testitem "TaskB: gross fold circuits agree with QuantumClifford" begin
    using QuantumClifford
    import FTShadows.TaskB as TB
    # the extra gross fold circuits are valid (the circuit part; the orders are in "Albert gross group order
    # 460800"): every 7th fold, conjugated independently
    code = TB.make_bb_code("gross", "fgh")
    _, _, _, folds, _ = TB.bb_families(code)
    @test length(folds) == 168
    n = code.n
    L = TB.logicals(code)
    for e in folds[1:7:end]
        st = TB.apply_gates!(TB._rows_tableau(L, n), e.gates)
        imgs = reduce(vcat, [permutedims([xbit(st[i]); zbit(st[i])]) for i in axes(L, 1)])
        @test TB.logical_coords(code, imgs) == e.logical
    end
end

@testitem "TaskB: invalid operations are rejected" begin
    using PythonCall
    import FTShadows.TaskB as TB
    code = TB.make_bb_code("BB72", "fgh")
    # a random qubit permutation (numpy's default_rng(3)) is almost surely not an automorphism
    perm = pyconvert(Vector{Int}, pyimport("numpy").random.default_rng(3).permutation(code.n))
    @test !TB.preserves_stabilizer_group(code, TB.perm_op(perm))
    s = TB.matching_matrix([i ⊻ 1 for i in 0:code.n-1])
    @test !TB.preserves_stabilizer_group(code, TB.z_diag(s))
end

@testitem "TaskB: orbits match single-orbit BFS" setup = [TaskBToys] begin
    import FTShadows.TaskB as TB
    gens = [TB.tableau_symplectic(TB.circuit_clifford(c, 3)) for c in
            [[TaskBToys.G(:H, 0), TaskBToys.G(:H, 1), TaskBToys.G(:H, 2)], [TaskBToys.G(:CX, 0, 1)], [TaskBToys.G(:SWAP, 1, 2)]]]
    lab = TB.orbits(gens, 6)
    for v in (1, 5, 17, 63)
        @test Set(findall(==(lab[v+1]), lab) .- 1) == Set(TB.orbit_of(gens, v, 6))
    end
end

@testitem "TaskB: orbit formula equals group enumeration" tags = [:slow] setup = [TaskBToys] begin
    import FTShadows.TaskB as TB
    for toy in ("phase_type", "css_only")
        group, _ = TB.enumerate_group(TaskBToys.TOY_GENS[toy], 3)
        lz, lx, _, _ = TB.channel_for_group(TaskBToys.toy_symplectic(toy), 3)
        @test lz ≈ TB.exact_lambda(group, 3, "Z")
        @test lx ≈ TB.exact_lambda(group, 3, "X")
    end
end

@testitem "TaskB: brute-force estimator mean and second moment" tags = [:slow] setup = [TaskBToys] begin
    using Statistics
    import FTShadows.TaskB as TB
    k = 3
    circs = TaskBToys.TOY_GENS["phase_type"]
    group, words = TB.enumerate_group(circs, k)
    lz, lx, _, _ = TB.channel_for_group(TaskBToys.toy_symplectic("phase_type"), k)
    for readout in ("Z", "X", "mix")
        rng = TB.numpy_rng(2026)
        lam = Dict("Z" => lz, "X" => lx, "mix" => 0.5 .* (lz .+ lx))[readout]
        psi = TB.haar_state(k, rng)
        paulis = Int[]
        for val in sort(unique(round.(lam[2:end]; digits=9)))
            val > 0 && push!(paulis, findfirst(i -> i > 1 && isapprox(lam[i], val), eachindex(lam)) - 1)
        end
        paulis = sort(unique([paulis; 0b000001; 0b001001; 0b111111]))
        shots = 30000
        est = TB.simulate_estimator(group, words, circs, psi, paulis, lam, k, shots, rng; readout)
        for (j, v) in enumerate(paulis)
            truth = TB.expval(psi, TB.pauli_from_vec(v, k))
            se = std(est[j, :]) / sqrt(shots)
            @test abs(mean(est[j, :]) - truth) < 5se
            se2 = std(est[j, :] .^ 2) / sqrt(shots)
            @test abs(mean(est[j, :] .^ 2) - 1 / lam[v+1]) < 5se2 + 1e-9
        end
    end
end

@testitem "TaskB: unlearnable Pauli is invisible to the ensemble" setup = [TaskBToys] begin
    using QuantumClifford
    import FTShadows.TaskB as TB
    k = 3
    circs = TaskBToys.TOY_GENS["css_only"]
    group, words = TB.enumerate_group(circs, k)
    lz, lx, _, _ = TB.channel_for_group(TaskBToys.toy_symplectic("css_only"), k)
    v = 1 | (2 << 3)   # X on qubit 0, Z on qubit 1: mixed
    @test lz[v+1] == 0 && lx[v+1] == 0
    p = TB.circuit_unitary([TaskBToys.G(:X, 0), TaskBToys.G(:Z, 1)], k)
    rho0 = Matrix{ComplexF64}(TB.I(8)) / 8
    rho1 = rho0 + 0.1 * p
    hall = TB.circuit_unitary([TaskBToys.G(:H, q) for q in 0:k-1], k)
    for gi in 1:max(1, length(group) ÷ 200):length(group)
        u = TB.circuit_unitary(reduce(vcat, [circs[i] for i in words[gi]]; init=TB.Gate[]), k)
        for basis in (TB.I(8), hall)
            uu = basis * u
            @test real.([(uu*rho0*uu')[i, i] for i in 1:8]) ≈ real.([(uu*rho1*uu')[i, i] for i in 1:8])
        end
    end
end

@testitem "TaskB: sum rule and global Clifford" setup = [TaskBToys] begin
    import FTShadows.TaskB as TB
    for k in (2, 3)
        lz, _, _, _ = TB.channel_for_group(TaskBToys.sp_gens(k), k)
        @test sum(lz[2:end]) ≈ 2^k - 1
        @test all(1 ./ lz[2:end] .≈ 2^k + 1)
    end
    lz2, _, _, _ = TB.channel_for_group(TaskBToys.toy_symplectic("css_only"), 3)
    @test sum(lz2[2:end]) ≈ 2^3 - 1
end

@testitem "TaskB: local Clifford gives 3^w" setup = [TaskBToys] begin
    import FTShadows.TaskB as TB
    for k in (2, 3, 4)
        gens = [TB.gates_to_symplectic([TaskBToys.G(g, q)], k) for q in 0:k-1 for g in (:H, :S)]
        lz, _, _, _ = TB.channel_for_group(gens, k)
        @test lz[2:end] ≈ 3.0 .^ (-Float64.(TB.logical_weights(k)[2:end]))
    end
end

@testitem "TaskB: CSS-type groups cannot learn mixed Paulis" begin
    import FTShadows.TaskB as TB
    code = TB.make_bb_code("BB72", "fgh")
    fam, auts, duals, _, _ = TB.bb_families(code)
    # no generator reduction here: the group, hence every λ, is the same without it
    lz, lx, _, _ = TB.channel_for_group([fam["shift_x"].logical, fam["shift_y"].logical,
                                         [e.logical for e in TB.dedupe([auts; duals])][1:6]...], 12)
    v = UInt32(0):UInt32((1 << 24) - 1)
    mixed = ((v .& 0xFFF) .!= 0) .& ((v .>> 12) .!= 0)
    @test !any(lz[mixed] .!= 0) && !any(lx[mixed] .!= 0)
    pure = (((v .& 0xFFF) .== 0) .⊻ ((v .>> 12) .== 0)) .& (v .!= 0)
    @test all(0.5 .* (lz .+ lx)[pure] .== 0.5)
end

@testitem "TaskB: Schreier-Sims replacement gives |Sp(2k,2)| and BFS orders" tags = [:slow] setup = [TaskBToys] begin
    using Oscar, PythonCall
    import FTShadows.TaskB as TB
    for k in 1:5
        @test TB.group_order(TaskBToys.sp_gens(k), 2k) == TB.symplectic_group_order(k)
        @test TB.group_order(TaskBToys.sp_gens(k), 2k; method=:oscar) == TB.symplectic_group_order(k)
    end
    # random words in the generators, drawn from numpy's default_rng(5) in a fixed order
    rng = pyimport("numpy").random.default_rng(5)
    for k in (2, 3), _ in 1:5
        allg = TaskBToys.sp_gens(k)
        els = map(1:2) do _
            m = TB.GF2Mat(TB.I(2k))
            for _ in 1:pyconvert(Int, rng.integers(1, 6))
                m = TB.gf2mul(m, allg[pyconvert(Int, rng.integers(length(allg)))+1])
            end
            m
        end
        @test TB.group_order(els, 2k) == TB.closure(els)[1]
    end
end

@testitem "TaskB: Albert gross group order 460800" tags = [:slow] begin
    using Oscar
    import FTShadows.TaskB as TB
    code = TB.make_bb_code("gross", "fgh")
    fam, auts, duals, folds, _ = TB.bb_families(code)
    sig = TB.bravyi_zx_perm(12, 6)
    mats, zmats = TB.GF2Mat[], TB.GF2Mat[]
    for a in 0:11, b in 0:5
        s = TB.matching_matrix(TB.shift_perm(12, 6, a, b)[sig.+1])
        zl = TB.logical_action(code, TB.z_diag(s))
        push!(zmats, zl)
        append!(mats, [zl, TB.logical_action(code, TB.x_diag(s))])
    end
    base = [mats..., fam["zx_bravyi"].logical, [e.logical for e in [auts; duals]]...]
    _, order = TB.reduce_generators(base, 24)
    @test order == 460800
    # the depth-one Z-diagonal family spans a 6-dimensional space of symmetric 12 × 12 matrices
    @test TB.gf2rank(reduce(vcat, [permutedims(vec(permutedims(z[1:12, 13:24]))) for z in zmats])) == 6
    # the second involutive-duality family: valid folds that strictly enlarge Albert's group
    @test length(folds) == 168
    _, o_all = TB.reduce_generators([base; [e.logical for e in folds]], 24)
    @test o_all == 11324620800 && o_all % order == 0
end

@testitem "TaskB: La-cross [[65,9]] generates the full Clifford group" tags = [:slow] begin
    using Oscar
    import FTShadows.TaskB as TB
    r = TB.lacross(3, 7)
    @test r["n_seed_automorphisms"] == 168     # |GL(3,2)|
    @test r["automorphism_lift_failures"] == 0
    @test r["ensembles"]["automorphisms"]["order"] == 168 * 168
    full = r["ensembles"]["aut+Hfold+Pfold"]
    @test full["is_full_Sp"] && full["n_nontrivial_orbits"] == 1
    @test full["lambdaZ_values"] == [TB.py_round(1 / 513, 12)]
    @test r["Sp_order"] == BigInt(order(symplectic_group(18, GF(2))))
end

@testitem "TaskB: fold gate lists preserve stabilizer signs" begin
    import FTShadows.TaskB as TB
    code = TB.make_bb_code("BB72", "fgh")
    _, _, _, folds, _ = TB.bb_families(code)
    flipped(gates) = any(TB.stabilizer_sign_flips(code.stabs, gates, code.n)[2])
    @test !any(flipped(e.gates) for e in folds)
    bare = count(flipped(e.gates_uncorrected) for e in folds)
    @test bare == count(e -> !isempty(e.pauli_correction), folds) == 126
    for e in folds
        @test TB.gates_to_symplectic(e.gates, code.n) == e.phys
        @test all(g -> g.name in (:X, :Y, :Z), e.pauli_correction)
    end
end

@testitem "TaskB: basis-invariant statistics and pinned per-weight counts" begin
    import FTShadows.TaskB as TB
    function audit_basis(code, seed)
        m = TB.frombits(TB.rng_inputs()["audit_random_css"][string(seed)])
        meta = filter(kv -> kv[1] != :_stab_row_basis, copy(code.meta))
        c = TB.StabilizerCode("BB72_alt_css_basis_$seed", code.n, code.stabs, TB.gf2mul(m, code.lx),
                              TB.gf2mul(permutedims(TB.gf2inv(m)), code.lz), meta)
        TB.validate(c)
        return c
    end
    code = TB.make_bb_code("BB72", "fgh")
    fam, auts, duals, folds, _ = TB.bb_families(code)
    k = 12
    w = TB.logical_weights(k)
    els = [fam["shift_x"], fam["shift_y"], fam["zx_bravyi"], folds[1:6]...]
    out = map((code, audit_basis(code, 11))) do c
        lz, lx, _, _ = TB.channel_for_group([TB.logical_action(c, e.phys) for e in els], k)
        mix = 0.5 .* (lz .+ lx)
        (spec=(sort(lz), sort(lx), sort(mix)), n_learn=count(>(0), mix[2:end]), split=TB.type_breakdown(mix, k),
         per_w=[r["n_learnable"] for r in TB.weight_table(mix, w, k; max_w=2)])
    end
    a, b = out
    @test a.spec == b.spec
    @test a.n_learn == b.n_learn
    @test a.split == b.split
    @test a.per_w == [24, 132]
    @test b.per_w == [24, 135]
    # the automorphism + ZX-duality group learns exactly the pure-X and pure-Z Paulis in further CSS bases
    els2 = [fam["shift_x"], fam["shift_y"], auts..., duals...]
    for seed in (21, 22)
        c = audit_basis(code, seed)
        lz, lx, _, _ = TB.channel_for_group([TB.logical_action(c, e.phys) for e in els2], k)
        tb = TB.type_breakdown(0.5 .* (lz .+ lx), k)
        @test tb["mixed"]["n_learnable"] == 0
        @test tb["pure_X"]["n_learnable"] == tb["pure_X"]["n"] == 2^k - 1
        @test tb["pure_Z"]["n_learnable"] == tb["pure_Z"]["n"]
    end
end

@testitem "TaskB: non-CSS logical basis is rejected" begin
    import FTShadows.TaskB as TB
    code = TB.make_bb_code("BB72", "fgh")
    lz2 = copy(code.lz)
    lz2[1, :] .⊻= code.lx[1, :]   # Zbar'_1 = Xbar_1 Zbar_1: symplectic, commutes with S, not physically Z-type
    c2 = TB.StabilizerCode("BB72_nonCSS", code.n, code.stabs, code.lx, lz2, code.meta)
    @test_throws r"not CSS" TB.validate(c2)
    @test TB.validate(c2; require_css=false) === nothing
end

@testitem "TaskB: brute-force samples equal the reference bit for bit" tags = [:slow, :bridge, :pyref] setup = [TaskBToys] begin
    using PythonCall
    import FTShadows.TaskB as TB
    # the reference setting: seed 2026, 30000 shots, the phase-type toy group, readouts Z, X and mix
    bf, sh, st, np = pyimport("taskb.bruteforce"), pyimport("taskb.shadow"), pyimport("stim"), pyimport("numpy")
    pycirc(c) = st.Circuit(c * "\nI 0 1 2")
    pyc = [pycirc(c) for c in ["H 0 1 2", "SWAP 0 1\nSWAP 1 2", "CZ 0 1\nS 2"]]
    pygroup = bf.enumerate_group(pyc, 3)
    pylam = sh.channel_for_group(pylist([bf.tableau_symplectic(st.Tableau.from_circuit(c)) for c in pyc]), 3)
    k = 3
    circs = TaskBToys.TOY_GENS["phase_type"]
    group, words = TB.enumerate_group(circs, k)
    lz, lx, _, _ = TB.channel_for_group(TaskBToys.toy_symplectic("phase_type"), k)
    @test lz == pyconvert(Vector{Float64}, pylam[0]) && lx == pyconvert(Vector{Float64}, pylam[1])
    for readout in ("Z", "X", "mix")
        lam = Dict("Z" => lz, "X" => lx, "mix" => 0.5 .* (lz .+ lx))[readout]
        paulis = [1, 9, 63, 3, 24]
        rng, prng = TB.numpy_rng(2026), np.random.default_rng(2026)
        psi = TB.haar_state(k, rng)
        ppsi = bf.haar_state(k, prng)
        @test psi == pyconvert(Vector{ComplexF64}, ppsi)
        shots = 30000
        est = TB.simulate_estimator(group, words, circs, psi, paulis, lam, k, shots, rng; readout)
        pest = bf.simulate_estimator(pygroup, ppsi, pylist(paulis), np.asarray(lam), k, shots, prng; readout)
        @test est == pyconvert(Matrix{Float64}, pest)
        # both generators are left in the same state
        @test pyconvert(Float64, rng.random()) == pyconvert(Float64, prng.random())
    end
end
