# Mixed-state stabilizer tracker. Each known stabilizer generator carries the
# record bitset that fixes its sign in a noiseless run and a bitset of hidden random bits (the
# implicit measurement inside a reset). Rows are kept in echelon form keyed by pivot; the
# elimination order is fixed, because it decides which record combinations become detectors.
#
# QuantumClifford's MixedDestabilizer could track the group, but not the record parities that
# go with each generator, so the tracker is written out here.

mutable struct Row
    vec::PVec
    rec::Bits
    hid::Bits
end

mutable struct Tracker
    alg::PauliAlgebra
    rows::Dict{Int,Row}  # pivot => row
    n_hidden::Int
end

Tracker(n::Integer) = Tracker(PauliAlgebra(n), Dict{Int,Row}(), 0)

function _reduce(t::Tracker, vec::PVec)
    rec = hid = Bits()
    top = nothing
    while !iszero(vec)
        p = pivot(t.alg, vec)
        row = get(t.rows, p, nothing)
        row === nothing && return vec, rec, hid, top
        top === nothing && (top = p)
        vec = xor(vec, row.vec)
        rec = xor(rec, row.rec)
        hid = xor(hid, row.hid)
    end
    return vec, rec, hid, top
end

function _insert!(t::Tracker, vec::PVec, rec::Bits, hid::Bits)
    while !iszero(vec)
        p = pivot(t.alg, vec)
        row = get(t.rows, p, nothing)
        if row === nothing
            t.rows[p] = Row(vec, rec, hid)
            return
        end
        vec = xor(vec, row.vec)
        rec = xor(rec, row.rec)
        hid = xor(hid, row.hid)
    end
end

"""
    query(t, vec)

`(rec, hid)` if `vec` is in the known group, else `nothing`.
"""
function query(t::Tracker, vec::PVec)
    left, rec, hid, _ = _reduce(t, vec)
    return iszero(left) ? (rec, hid) : nothing
end

"""
    measure!(t, vec, rec_id; hidden=false)

Projective measurement of `vec` with outcome record `rec_id` (or a fresh hidden bit). Returns
`nothing` if the outcome is random, else `(rec, hid)` with outcome = parity(rec) ⊻ parity(hid).
"""
function measure!(t::Tracker, vec::PVec, rec_id::Union{Int,Nothing}=nothing; hidden::Bool=false)
    if hidden
        my_rec, my_hid = Bits(), bit(t.n_hidden)
        t.n_hidden += 1
    else
        my_rec, my_hid = bit(rec_id), Bits()
    end
    anti = [p for (p, r) in t.rows if anticommutes(t.alg, r.vec, vec)]
    if !isempty(anti)
        p0 = minimum(anti)  # smallest pivot keeps the echelon form when multiplied in
        g0 = pop!(t.rows, p0)
        for p in anti
            p == p0 && continue
            r = t.rows[p]
            r.vec = xor(r.vec, g0.vec)
            r.rec = xor(r.rec, g0.rec)
            r.hid = xor(r.hid, g0.hid)
        end
        _insert!(t, vec, my_rec, my_hid)
        return nothing
    end
    left, rec, hid, top = _reduce(t, vec)
    if !iszero(left)
        _insert!(t, vec, my_rec, my_hid)
        return nothing
    end
    # replacing the top row keeps its pivot, so the echelon form survives
    t.rows[top] = Row(vec, my_rec, my_hid)
    return rec, hid
end

"""
    reset!(t, q, basis)

Reset data qubit `q` to |0> (basis 'Z') or |+> (basis 'X').
"""
function reset!(t::Tracker, q::Integer, basis::Char)
    a = t.alg
    meas = from_dict(a, ((q, basis),))
    flip = from_dict(a, ((q, basis == 'Z' ? 'X' : 'Z'),))
    measure!(t, meas; hidden=true)
    h = bit(t.n_hidden - 1)
    for r in values(t.rows)
        anticommutes(a, r.vec, flip) && (r.hid = xor(r.hid, h))
    end
    return t
end
