# PyMatching prediction with the exact complementary gap.

"""
    MatchingGapDecoder(dem)

PyMatching prediction plus the exact complementary gap, for a decomposed (graphlike) DEM.

For observable j, a detector potential s with s_u ⊻ s_v = label on every bulk edge moves the
observable onto the boundary: l'(u) = l(e) ⊻ s_u. Boundary edges with l' = 1 go to an extra node A,
the others to B, and hypothesis o is one matching with A a defect iff o ⊻ Σ_T s_t = 1. The smaller of
the two weights is checked against the plain matching on every call.

If a bulk cycle flips j (e.g. conjugate termination at r = 1) there is no potential, and the gap is
NaN with `gap_exact[j] = false`. Detector-free mechanisms that flip j cap the gap at their weight.
The split graph keeps all parallel mechanisms (the prediction graph only the most likely one), and a
hypothesis with no valid matching has weight `Inf`.
"""
mutable struct MatchingGapDecoder
    nd::Int
    no::Int
    undetectable_w::Vector{Float64}
    base::Py
    pot::Vector{Union{Nothing,Vector{UInt8}}}
    split::Vector{Union{Nothing,Py}}
    gap_exact::Vector{Bool}
    max_weight_mismatch::Float64
    pred_mismatches::Int
end

MatchingGapDecoder(dem::Py) = MatchingGapDecoder(read_dem(dem))

function MatchingGapDecoder(dem::Dem)
    pm = pymodule(:pymatching)
    nd, no = dem.num_detectors, dem.num_observables
    undet = _MergedMechanisms{UInt64}()
    edges = graphlike_edges(dem; undetectable = undet)
    undetectable_w = fill(Inf, no)
    for (om, pr) in zip(undet.keys, undet.p), j in 0:(no - 1)
        (om >> j) & 1 == 1 && (undetectable_w[j + 1] = min(undetectable_w[j + 1], llr(pr)))
    end
    # most likely mechanism per detector pair, kept at the position the pair first appeared
    best = Dict{Tuple{Int,Int},Int}()
    order = Tuple{Int,Int}[]
    choice = Tuple{UInt64,Float64}[]
    for ((u, v, om), pr) in zip(edges.keys, edges.p)
        i = get(best, (u, v), 0)
        if i == 0
            push!(order, (u, v))
            push!(choice, (om, pr))
            best[(u, v)] = length(order)
        elseif pr > choice[i][2]
            choice[i] = (om, pr)
        end
    end
    base = pm.Matching()
    for ((u, v), (om, pr)) in zip(order, choice)
        fid = pybuiltins.set(pylist([j for j in 0:(no - 1) if (om >> j) & 1 == 1]))
        if v == -1
            base.add_boundary_edge(u; fault_ids = fid, weight = llr(pr), error_probability = pr)
        else
            base.add_edge(u, v; fault_ids = fid, weight = llr(pr), error_probability = pr)
        end
    end
    pot = Union{Nothing,Vector{UInt8}}[]
    split = Union{Nothing,Py}[]
    gap_exact = Bool[]
    A, B = nd, nd + 1
    for j in 0:(no - 1)
        s = _potential(edges, nd, j)
        push!(gap_exact, s !== nothing)
        push!(pot, s)
        if s === nothing
            push!(split, nothing)
            continue
        end
        m = pm.Matching()
        for ((u, v, om), pr) in zip(edges.keys, edges.p)  # every mechanism, not only the most likely
            lab = UInt8((om >> j) & 1)
            m.add_edge(u, v == -1 ? (lab ⊻ s[u + 1] == 1 ? A : B) : v; weight = llr(pr),
                merge_strategy = "smallest-weight")
        end
        push!(split, m)
    end
    return MatchingGapDecoder(nd, no, undetectable_w, base, pot, split, gap_exact, 0.0, 0)
end

# Potential s with s_u ⊻ s_v = label_j on every bulk edge (breadth-first, roots in detector order),
# or `nothing` if a bulk cycle flips observable j.
function _potential(edges::_MergedMechanisms, nd::Int, j::Int)
    adj = [Tuple{Int,UInt8}[] for _ in 1:nd]
    for (u, v, om) in edges.keys
        v == -1 && continue
        lab = UInt8((om >> j) & 1)
        push!(adj[u + 1], (v, lab))
        push!(adj[v + 1], (u, lab))
    end
    s = zeros(UInt8, nd)
    seen = falses(nd)
    queue = Int[]
    for root in 0:(nd - 1)
        seen[root + 1] && continue
        seen[root + 1] = true
        push!(queue, root)
        while !isempty(queue)
            u = popfirst!(queue)
            for (v, lab) in adj[u + 1]
                if !seen[v + 1]
                    seen[v + 1] = true
                    s[v + 1] = s[u + 1] ⊻ lab
                    push!(queue, v)
                elseif s[v + 1] != s[u + 1] ⊻ lab
                    return nothing
                end
            end
        end
    end
    return s
end

"""
    decode(dec, dets; with_gap=true) -> (pred, gap, kind)

Decode a shots × detectors 0/1 matrix. Returns predicted observable flips (shots × observables,
`UInt8`), the per-shot confidence (complementary gap; minimum over observables) and its kind (see
[`CONF_KINDS`](@ref)).
"""
function decode(dec::MatchingGapDecoder, dets::AbstractMatrix; with_gap::Bool = true)
    n = size(dets, 1)
    D = convert(Matrix{UInt8}, dets)
    res = dec.base.decode_batch(to_numpy_u8(D); return_weights = true)
    pred = reshape(_from_numpy_u8(res[0]), n, :)
    w_base = pyconvert(Vector{Float64}, res[1])
    with_gap || return pred, fill(NaN, n), fill("unavailable", n)
    par = [sum(Int, view(D, i, :)) % 2 for i in 1:n]
    gaps = fill(Inf, n)
    for j in 1:dec.no
        m = dec.split[j]
        if m === nothing
            gaps = fill(NaN, n)
            continue
        end
        sT = (Int.(D) * Int.(dec.pot[j])) .% 2
        w = map((0, 1)) do o
            a = UInt8.(sT .⊻ o)
            bb = UInt8.((par .+ a) .% 2)
            _matching_weights(m, hcat(D, a, bb))
        end
        if n > 0
            mm = maximum(abs.(min.(w[1], w[2]) .- w_base))
            dec.max_weight_mismatch = _pymax(dec.max_weight_mismatch, mm)
            mm > 1e-6 * _pymax(1.0, maximum(w_base)) &&
                error("gap decoder inconsistent with plain matching (|dW| = $mm)")
        end
        bad = count(i -> UInt8(w[2][i] < w[1][i]) != pred[i, j] && !_isclose(w[2][i], w[1][i]), 1:n)
        dec.pred_mismatches += bad
        bad > 0 && error("gap decoder prediction differs from plain matching on $bad shots")
        g = abs.(w[2] .- w[1])
        gaps = _npmin.(gaps, _npmin.(g, dec.undetectable_w[j]))
    end
    return pred, gaps, [isnan(g) ? "unavailable" : "exact" for g in gaps]
end

# Matching weights; rows without a valid matching (impossible hypothesis) get Inf.
function _matching_weights(m::Py, ext::Matrix{UInt8})
    try
        return pyconvert(Vector{Float64}, m.decode_batch(to_numpy_u8(ext); return_weights = true)[1])
    catch e
        _is_value_error(e) || rethrow()
        return map(1:size(ext, 1)) do i
            try
                pyconvert(Float64, m.decode(to_numpy_u8(ext[i, :]); return_weight = true)[1])
            catch e2
                _is_value_error(e2) || rethrow()
                Inf
            end
        end
    end
end

_is_value_error(e) = e isa PyException && pyisinstance(e.v, pybuiltins.ValueError)
