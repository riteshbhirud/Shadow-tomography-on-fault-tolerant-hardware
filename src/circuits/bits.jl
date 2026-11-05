# Immutable bitsets over measurement-record ids: xor, highest set bit, popcount. Stored as
# little-endian UInt64 chunks without trailing zero chunks, so == and hash compare the vectors.

struct Bits
    c::Vector{UInt64}
end

Bits() = Bits(UInt64[])

function _trim!(c::Vector{UInt64})
    while !isempty(c) && iszero(c[end])
        pop!(c)
    end
    return c
end

"""
    bit(i)

The bitset with only bit `i` (0-based) set.
"""
function bit(i::Integer)
    c = zeros(UInt64, (i >> 6) + 1)
    c[end] = UInt64(1) << (i & 63)
    return Bits(c)
end

function Base.xor(a::Bits, b::Bits)
    la, lb = length(a.c), length(b.c)
    la < lb && return xor(b, a)
    c = copy(a.c)
    @inbounds for i in 1:lb
        c[i] ⊻= b.c[i]
    end
    la == lb && _trim!(c)
    return Bits(c)
end

function Base.:&(a::Bits, b::Bits)
    l = min(length(a.c), length(b.c))
    return Bits(_trim!([a.c[i] & b.c[i] for i in 1:l]))
end

Base.iszero(a::Bits) = isempty(a.c)
Base.:(==)(a::Bits, b::Bits) = a.c == b.c
Base.hash(a::Bits, h::UInt) = hash(a.c, hash(:Bits, h))
Base.count_ones(a::Bits) = sum(count_ones, a.c; init=0)

"""
    topbit(a)

Index of the highest set bit, or -1 for zero.
"""
topbit(a::Bits) = iszero(a) ? -1 : 64 * (length(a.c) - 1) + 63 - leading_zeros(a.c[end])

isset(a::Bits, i::Integer) = (i >> 6) < length(a.c) && ((a.c[(i >> 6) + 1] >> (i & 63)) & 1) == 1

overlaps(a::Bits, b::Bits) = any(i -> (a.c[i] & b.c[i]) != 0, 1:min(length(a.c), length(b.c)))

"""
    setbits(a)

The set bit indices of `a`, ascending (0-based).
"""
function setbits(a::Bits)
    out = Int[]
    for (k, w) in enumerate(a.c)
        while w != 0
            push!(out, 64 * (k - 1) + trailing_zeros(w))
            w &= w - 1
        end
    end
    return out
end

bits_of(ids) = foldl(xor, (bit(i) for i in ids); init=Bits())

"""
    rank_bits(vecs)

GF(2) rank of a list of bitsets.
"""
function rank_bits(vecs)
    rows = Dict{Int,Bits}()
    r = 0
    for v in vecs
        while !iszero(v)
            p = topbit(v)
            if haskey(rows, p)
                v = xor(v, rows[p])
            else
                rows[p] = v
                r += 1
                break
            end
        end
    end
    return r
end
