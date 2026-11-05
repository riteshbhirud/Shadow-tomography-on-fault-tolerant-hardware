# Fast checks against the small fixtures in test/reference/taskb. They do not need Oscar.

@testitem "TaskB fixtures: GF(2) routines match stored results" setup = [TaskBRef] begin
    import FTShadows.TaskB as TB
    ref = TaskBRef.small_refs()
    for g in ref["gf2"]
        m = TB.frombits(g["m"])
        r, piv = TB.rref(m)
        @test TB.bitrows(r) == g["rref"]
        @test piv .- 1 == g["pivots"]
        ns = TB.nullspace_gf2(m)
        @test isempty(g["nullspace"]) ? size(ns, 1) == 0 : TB.bitrows(ns) == g["nullspace"]
        x = TB.solve_left(m, TB.frombits(g["b"]))
        @test x === nothing ? g["solve_left_b"] === nothing : TB.bitrows(x) == g["solve_left_b"]
        @test TB.bitrows(TB.solve_left(m, TB.frombits(g["consistent_rhs"]))) == g["solve_left_consistent"]
    end
end

@testitem "TaskB fixtures: logical bases match the stored ones" setup = [TaskBRef] begin
    import FTShadows.TaskB as TB
    ref = TaskBRef.small_refs()["bases"]
    for key in ("BB72/fgh", "BB72/generic", "gross/fgh", "gross/generic")
        name, basis = split(key, "/")
        c = TB.make_bb_code(name, basis)
        @test TB.bitrows(c.lx[:, 1:c.n]) == ref[key]["xl"]
        @test TB.bitrows(c.lz[:, c.n+1:end]) == ref[key]["zl"]
    end
end

@testitem "TaskB fixtures: Bravyi-style search seeds" begin
    import FTShadows.TaskB as TB
    # only f, g, h come from the stored RNG draws; the greedy completion (and whether it fails) is recomputed here
    for (s, d) in TB.rng_inputs()["bb72_fgh_search"]
        res = TB.bb72_fgh_basis(parse(Int, s))
        @test (res !== nothing) == d["python_found_basis"]
        if res !== nothing
            @test TB.bitrows(res[1]) == d["xl"]
            @test TB.bitrows(res[2]) == d["zl"]
        end
    end
end

@testitem "TaskB fixtures: toy shadow channels match stored values" setup = [TaskBRef] begin
    import FTShadows.TaskB as TB
    for (name, d) in TaskBRef.small_refs()["toy_lambda"]
        gens = [TB.frombits(g) for g in d["generators"]]
        lz, lx, size_, lab = TB.channel_for_group(gens, 3)
        @test TaskBRef.same(lz, d["lambda_Z"])
        @test TaskBRef.same(lx, d["lambda_X"])
        @test size_ == d["orbit_size"]
        @test lab == d["labels"]
        @test TaskBRef.same(TB.weight_table(0.5 .* (lz .+ lx), TB.logical_weights(3), 3), d["weight_table_mix"])
    end
end

@testitem "TaskB fixtures: BB72 Pauli corrections and λ_mix sample" setup = [TaskBRef] begin
    import FTShadows.TaskB as TB
    ref = TaskBRef.small_refs()
    code = TB.make_bb_code("BB72", "fgh")
    ens, folds, _ = TB.build_ensembles(code)
    @test TaskBRef.same([Dict("label" => e.label, "pauli_correction" => TaskBRef.gate_lists(e.pauli_correction))
                         for e in folds[1:12:end]], ref["bb72_fold_corrections"])
    lam_ref = ref["bb72_lambda_mix_weight_le_2"]
    idx = lam_ref["indices"] .+ 1
    for ename in ("aut", "albert_two_fold")
        lz, lx, _, _ = TB.channel_for_group([e.logical for e in ens[ename]], 12)
        @test TaskBRef.same((0.5 .* (lz .+ lx))[idx], lam_ref[ename])
    end
end

@testitem "TaskB fixtures: the record comparison is strict" setup = [TaskBRef] begin
    py = Dict("a" => Any[1, 2.5])
    @test TaskBRef.same(Dict("a" => Any[1, 2.5]), py)
    @test !TaskBRef.same(Dict("a" => [1.0, 2.5]), py)       # int vs float
    @test !TaskBRef.same(Dict("a" => Any[1, nextfloat(2.5)]), py)
    @test !TaskBRef.same(Dict("a" => Any[1, 2.5], "b" => 2), py)
end
