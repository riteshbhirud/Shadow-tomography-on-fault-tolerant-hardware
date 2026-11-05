# La-cross codes (Di Fini et al. 2610.12685): cheap logical Cliffords built physically, and their shadow map.

"""
    seed_automorphisms(h)

Bit permutations P (as matrices) with rowspace(H P) = rowspace(H), each with the Q solving Q H = H P, in
lexicographic order of the permutation.
"""
function seed_automorphisms(h::AbstractMatrix)
    n = size(h, 2)
    rb = row_basis(h)
    out = Tuple{GF2Mat,GF2Mat}[]
    for perm in _permutations(n)
        p = zeros(Bool, n, n)
        for i in 1:n
            p[i, perm[i]] = true
        end
        hp = gf2mul(h, p)
        gf2rank(vcat(rb, hp)) == size(rb, 1) || continue
        q = solve_left(h, hp)
        @assert q !== nothing
        push!(out, (p, q))
    end
    return out
end

# all permutations of 1:n in lexicographic order
function _permutations(n::Integer)
    out = Vector{Int}[]
    p = collect(1:n)
    while true
        push!(out, copy(p))
        i = n - 1
        while i >= 1 && p[i] > p[i+1]
            i -= 1
        end
        i == 0 && return out
        j = n
        while p[j] < p[i]
            j -= 1
        end
        p[i], p[j] = p[j], p[i]
        reverse!(p, i + 1, n)
    end
end

"""
    lift_candidates(p1, q1, p2, q2)

Candidate X-part maps on the check sector, tried in this order; the first one that preserves the code is used.
"""
function lift_candidates(p1, q1, p2, q2)
    qi1, qi2 = gf2inv(q1), gf2inv(q2)
    t = permutedims
    return [
        "Q1^-T (x) Q2^-T" => gf2(kron(t(qi1), t(qi2))),
        "Q1^-T (x) Q2" => gf2(kron(t(qi1), q2)),
        "Q1 (x) Q2^-T" => gf2(kron(q1, t(qi2))),
        "Q1 (x) Q2" => gf2(kron(q1, q2)),
        "Q1^T (x) Q2^T" => gf2(kron(t(q1), t(q2))),
    ]
end

"""
    x_linear_op(lx)

Physical symplectic of the CNOT circuit x -> x Lx, z -> z Lx^{-T}.
"""
function x_linear_op(lx::AbstractMatrix)
    n = size(lx, 1)
    return GF2Mat([lx zeros(Bool, n, n); zeros(Bool, n, n) permutedims(gf2inv(lx))])
end

"""
    transpose_perm(n1, m1)

Qubit permutation (0-based) that transposes both hypergraph-product sectors.
"""
transpose_perm(n1::Integer, m1::Integer) =
    [[j * n1 + i for i in 0:n1-1 for j in 0:n1-1]; [n1 * n1 + b * m1 + a for a in 0:m1-1 for b in 0:m1-1]]

"""
    lacross(kk=3, nn=7; orders=true, timings=false)

Generators, group orders and shadow statistics of the La-cross code H_{kk,nn} (the `lacross.json` record).
Group orders need `using Oscar`; `orders=false` skips them. `timings=true` adds the run time.
"""
function lacross(kk::Integer=3, nn::Integer=7; orders::Bool=true, timings::Bool=false)
    t0 = time()
    code = make_lacross_code(kk, nn)
    h = lacross_seed(kk, nn)
    m1, n1 = size(h)
    k = logical_k(code)
    auts = seed_automorphisms(h)
    res = OrderedDict{String,Any}("code" => code.name, "n" => code.n, "k" => k, "n_seed_automorphisms" => length(auts),
                           "automorphism_lift_failures" => 0)
    gens, labels = GF2Mat[], String[]
    lift_used = OrderedDict{String,Int}()
    eye_n, eye_q = GF2Mat(I(n1)), GF2Mat(I(m1))
    for which in (0, 1), (p, q) in auts[2:end]     # all non-identity ones
        p1, q1, p2, q2 = which == 0 ? (p, q, eye_n, eye_q) : (eye_n, eye_q, p, q)
        ok = false
        for (name, r) in lift_candidates(p1, q1, p2, q2)
            lx = zeros(Bool, code.n, code.n)
            lx[1:n1*n1, 1:n1*n1] = gf2(kron(p1, p2))
            lx[n1*n1+1:end, n1*n1+1:end] = r
            f = x_linear_op(lx)
            if preserves_stabilizer_group(code, f)
                push!(gens, logical_action(code, f))
                push!(labels, "automorphism lift side=$which ($name)")
                lift_used[name] = get(lift_used, name, 0) + 1
                ok = true
                break
            end
        end
        res["automorphism_lift_failures"] += !ok
    end
    tp = transpose_perm(n1, m1)
    f_h = gf2mul(perm_op(tp), h_all(code.n))
    preserves_stabilizer_group(code, f_h) || error("Hadamard fold invalid")
    push!(gens, logical_action(code, f_h)); push!(labels, "Hadamard fold F_tau")
    f_p = z_diag(matching_matrix(tp))
    preserves_stabilizer_group(code, f_p) || error("phase fold invalid")
    push!(gens, logical_action(code, f_p)); push!(labels, "phase fold P_tau")
    res["lift_variants_used"] = lift_used
    res["generator_labels"] = labels
    d = 2k
    sets = [
        "automorphisms" => [g for (g, lb) in zip(gens, labels) if startswith(lb, "automorphism")],
        "aut+Hfold" => [g for (g, lb) in zip(gens, labels) if !startswith(lb, "phase")],
        "aut+Hfold+Pfold" => gens,
    ]
    res["ensembles"] = OrderedDict{String,Any}()
    wts = logical_weights(k)
    sp = symplectic_group_order(k)
    for (name, gs) in sets
        lab = orbits(gs, d)
        lz, lx, _ = lambdas_from_labels(lab, k)
        mix = 0.5 .* (lz .+ lx)
        e = OrderedDict{String,Any}()
        if orders
            km, o = reduce_generators(gs, d)   # GAP is much faster on the reduced generator set
            e["order"] = o
            e["order_second_seed"] = group_order(km, d)   # an independent recognition run
            e["is_full_Sp"] = o == sp
        end
        merge!(e, OrderedDict{String,Any}(
            "n_nontrivial_orbits" => length(unique(lab)) - 1,
            "frac_learnable_mix" => np_mean(mix[2:end] .> 0),
            "by_weight_mix" => weight_table(mix, wts, k; max_w=4),
            "lambdaZ_values" => sort!(unique(py_round.(lz[2:end], 12)))[1:min(end, 10)]))
        res["ensembles"][name] = e
    end
    res["Sp_order"] = sp
    res["basis"] = "generic CSS basis (Gaussian elimination); irrelevant when G = Sp(2k,2) (lambda uniform)"
    timings && (res["seconds"] = py_round(time() - t0, 1))
    return res
end
