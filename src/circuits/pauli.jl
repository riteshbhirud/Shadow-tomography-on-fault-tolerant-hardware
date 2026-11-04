# Phase-free Paulis. `PauliDict` is the {qubit => 'X'|'Y'|'Z'} dict used throughout; it is
# insertion-ordered because that order reaches the circuits (gate order of checks without an
# explicit order, local detectors). `PVec` is the tracker's bit encoding.

"""
    PauliDict

Phase-free Pauli as an insertion-ordered `OrderedDict{Int,Char}` from 0-based qubit to 'X', 'Y' or 'Z'.
"""
const PauliDict = OrderedDict{Int,Char}

pauli(ps) = PauliDict(Int(q) => Char(s) for (q, s) in ps)
pauli(d::PauliDict) = d

"""
    sortedpauli(p)

The pairs of `p` sorted by qubit (Python's `sorted(p.items())`).
"""
sortedpauli(p) = sort!(collect(p); by=first)

# x bits at 0..n-1 and z bits at W..W+n-1 with W a multiple of 64, so the two halves are
# chunk-aligned. A single integer x | z << n would give the same pivot order (z bits above x bits, each by qubit) is
# order, and pivots are only compared, never used as numbers.
struct PVec
    v::Vector{UInt64}
end
Base.:(==)(a::PVec, b::PVec) = a.v == b.v
Base.hash(a::PVec, h::UInt) = hash(a.v, hash(:PVec, h))
Base.xor(a::PVec, b::PVec) = PVec(a.v .⊻ b.v)
Base.iszero(a::PVec) = all(iszero, a.v)

struct PauliAlgebra
    n::Int
    nc::Int  # chunks per half
end
PauliAlgebra(n::Integer) = PauliAlgebra(n, max(1, cld(n, 64)))

zero_pvec(alg::PauliAlgebra) = PVec(zeros(UInt64, 2 * alg.nc))

function from_dict(alg::PauliAlgebra, p)
    v = zeros(UInt64, 2 * alg.nc)
    for (q, s) in p
        0 <= q < alg.n || throw(ArgumentError("qubit $q out of range"))
        s in ('X', 'Y', 'Z') || throw(ArgumentError("bad Pauli letter $s"))
        w, b = (q >> 6) + 1, UInt64(1) << (q & 63)
        s in ('X', 'Y') && (v[w] |= b)
        s in ('Z', 'Y') && (v[alg.nc+w] |= b)
    end
    return PVec(v)
end

function to_dict(alg::PauliAlgebra, a::PVec)
    out = PauliDict()
    for q in 0:alg.n-1
        w, s = (q >> 6) + 1, q & 63
        xb = (a.v[w] >> s) & 1 == 1
        zb = (a.v[alg.nc+w] >> s) & 1 == 1
        (xb || zb) && (out[q] = xb && zb ? 'Y' : (xb ? 'X' : 'Z'))
    end
    return out
end

function anticommutes(alg::PauliAlgebra, a::PVec, b::PVec)
    nc = alg.nc
    s = 0
    @inbounds for i in 1:nc
        s += count_ones((a.v[i] & b.v[nc+i]) ⊻ (a.v[nc+i] & b.v[i]))
    end
    return isodd(s)
end

"""
    pivot(alg, a)

Highest set bit of `a` in the padded layout, or -1 for the identity.
"""
function pivot(alg::PauliAlgebra, a::PVec)
    @inbounds for i in length(a.v):-1:1
        a.v[i] != 0 && return 64 * (i - 1) + 63 - leading_zeros(a.v[i])
    end
    return -1
end

"""
    support_overlaps(alg, a, b)

Whether the qubit supports of `a` and `b` intersect.
"""
function support_overlaps(alg::PauliAlgebra, a::PVec, b::PVec)
    nc = alg.nc
    @inbounds for i in 1:nc
        ((a.v[i] | a.v[nc+i]) & (b.v[i] | b.v[nc+i])) != 0 && return true
    end
    return false
end

"""
    pauli_product(n, dicts)

Phase-free product of Pauli dicts on `n` qubits (surface.pauli_product); ascending qubit order.
"""
function pauli_product(n::Integer, dicts)
    alg = PauliAlgebra(n)
    v = zero_pvec(alg)
    for d in dicts
        v = xor(v, from_dict(alg, d))
    end
    return to_dict(alg, v)
end

function vec_to_dict(v::AbstractVector)
    n = length(v) ÷ 2
    out = PauliDict()
    for q in 0:n-1
        x, z = v[q+1], v[n+q+1]
        (x || z) && (out[q] = x && z ? 'Y' : (x ? 'X' : 'Z'))
    end
    return out
end

function dict_to_vec(p, n::Integer)
    v = falses(2n)
    for (q, s) in p
        s in ('X', 'Y') && (v[q+1] = true)
        s in ('Z', 'Y') && (v[n+q+1] = true)
    end
    return Vector{Bool}(v)
end
