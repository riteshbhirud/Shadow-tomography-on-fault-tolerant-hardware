# Code distance and circuit distance. The circuit distance comes from stim's
# undetectable-logical-error search through the bridge; there is no Julia equivalent.
# The randomized ISD upper bound draws its column permutations from numpy through the bridge,
# so a given seed finds the same vectors; the tests also check code distances with an exact MIP.

"""
    code_capacity_circuit(n, checks, observables; p=0.01)

Ideal MPP of checks and observables, one DEPOLARIZE1 layer, ideal MPP again.
"""
function code_capacity_circuit(n::Integer, checks, observables; p::Real=0.01)
    c = Circuit()
    foreach(ch -> append!(c, "MPP", mpp_targets(ch)), checks)
    foreach(ob -> append!(c, "MPP", mpp_targets(ob)), observables)
    append!(c, "DEPOLARIZE1", 0:n-1, [p])
    foreach(ch -> append!(c, "MPP", mpp_targets(ch)), checks)
    foreach(ob -> append!(c, "MPP", mpp_targets(ob)), observables)
    nc, no = length(checks), length(observables)
    tot = 2 * (nc + no)
    for i in 0:nc-1
        append!(c, "DETECTOR", [rec(i - tot), rec(nc + no + i - tot)])
    end
    for j in 0:no-1
        append!(c, "OBSERVABLE_INCLUDE", [rec(nc + j - tot), rec(2nc + no + j - tot)], [j])
    end
    return c
end

"""
    stim_search(circuit; max_dets=6, max_degree=8)

Weight of the smallest undetectable error that flips an observable, from stim's
`search_for_undetectable_logical_errors` within its exploration bounds (bridge).
"""
function stim_search(circuit::Circuit; max_dets::Integer=6, max_degree::Integer=8)
    errs = pycircuit(circuit).search_for_undetectable_logical_errors(;
        dont_explore_detection_event_sets_with_size_above=max_dets,
        dont_explore_edges_with_degree_above=max_degree,
        dont_explore_edges_increasing_symptom_degree=false,
        canonicalize_circuit_errors=true)
    return pyconvert(Int, pybuiltins.len(errs))
end

"""
    code_distance_stim(n, checks; max_dets=6, max_degree=8)

`(code distance, k)`: stim search on each half (X-like, Z-like) of a symplectic logical basis.
"""
function code_distance_stim(n::Integer, checks; max_dets::Integer=6, max_degree::Integer=8)
    S = Matrix{Bool}(permutedims(reduce(hcat, [dict_to_vec(c, n) for c in checks])))
    pairs = logical_basis(S)
    best = typemax(Int)
    for half in 1:2
        obs = [vec_to_dict(pr[half]) for pr in pairs]
        best = min(best, stim_search(code_capacity_circuit(n, checks, obs); max_dets, max_degree))
    end
    return best, length(pairs)
end

"""
    min_weight_codeword(G, Gtriv; iters=2000, seed=0, weight=count)

Randomized information-set search for a low-weight vector in the row span
of `G` that is not in the row span of `Gtriv`: per iteration, row-reduce `G` with its columns in the
order `default_rng(seed).permutation(n)` and try every row and the sums of rows less than 6 apart.
Returns `(weight, vector)`, an upper bound on the minimum, or `(10^9, nothing)`.
"""
function min_weight_codeword(G::AbstractMatrix, Gtriv::AbstractMatrix; iters::Integer=2000, seed::Integer=0, weight=count)
    rng = pymodule(:numpy).random.default_rng(seed)
    n = size(G, 2)
    triv = RowSpace(Gtriv)
    best = (10^9, nothing)
    for _ in 1:iters
        perm = pyconvert(Vector{Int}, rng.permutation(n)) .+ 1
        R, _ = rref(G[:, perm])
        R = R[:, invperm(perm)]
        cands = [R[i, :] for i in axes(R, 1)]
        for i in axes(R, 1), j in i+1:min(size(R, 1), i + 5)
            push!(cands, R[i, :] .⊻ R[j, :])
        end
        for v in cands
            w = weight(v)
            if 0 < w < best[1] && (length(triv) == 0 || any(residual(triv, v)))
                best = (w, copy(v))
            end
        end
    end
    return best
end

"""
    code_distance_isd(n, checks; iters=300, seed=0)

`(weight, Pauli dict)` of the lightest nontrivial logical found by [`min_weight_codeword`](@ref) over
the normalizer of `checks`: an upper bound on the code distance (distance.code_distance_isd).
"""
function code_distance_isd(n::Integer, checks; iters::Integer=300, seed::Integer=0)
    S = Matrix{Bool}(permutedims(reduce(hcat, [dict_to_vec(c, n) for c in checks])))
    w, v = min_weight_codeword(normalizer(S), S; iters, seed, weight=v -> count(v[1:n] .| v[n+1:end]))
    return w, v === nothing ? nothing : vec_to_dict(v)
end
