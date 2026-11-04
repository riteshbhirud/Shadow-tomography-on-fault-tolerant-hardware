# GF(2) linear algebra on Bool matrices. Symplectic vectors have length 2n, laid out (x | z).
# Column indices are 1-based here; the Pauli dicts built from them use 0-based qubit ids.

"""
    rref(M)

Row-reduced echelon form over GF(2), eliminating above and below each pivot, with the first
nonzero row of a column as pivot row (as gf2.rref). Returns `(R, pivot_columns)`.
"""
function rref(M::AbstractMatrix)
    R = Matrix{Bool}(isodd.(M))
    m, n = size(R)
    piv = Int[]
    r = 1
    for c in 1:n
        r > m && break
        p = findfirst(view(R, r:m, c))
        p === nothing && continue
        p += r - 1
        if p != r
            R[[r, p], :] = R[[p, r], :]
        end
        for i in 1:m
            if i != r && R[i, c]
                @views R[i, :] .= R[i, :] .⊻ R[r, :]
            end
        end
        push!(piv, c)
        r += 1
    end
    return R[1:r-1, :], piv
end

gf2rank(M::AbstractMatrix) = isempty(M) ? 0 : length(rref(M)[2])

"""
    nullspace_gf2(M)

Basis of {x : M x = 0} as rows, one per free column in increasing order (gf2.nullspace).
"""
function nullspace_gf2(M::AbstractMatrix)
    m, n = size(M)
    R, piv = rref(M)
    pset = Set(piv)
    free = [c for c in 1:n if !(c in pset)]
    out = falses(length(free), n)
    for (i, f) in enumerate(free)
        out[i, f] = true
        for (r, pc) in enumerate(piv)
            out[i, pc] = R[r, f]
        end
    end
    return Matrix{Bool}(out)
end

function symp_inner(a::AbstractVector, b::AbstractVector)
    n = length(a) ÷ 2
    s = false
    @inbounds for i in 1:n
        s ⊻= (a[i] & b[n+i]) ⊻ (a[n+i] & b[i])
    end
    return s
end

css_stabilizers(HX, HZ) = Matrix{Bool}([HX falses(size(HX)); falses(size(HZ)) HZ])

"""
    normalizer(S)

Basis of the Paulis commuting with every row of `S` (gf2.normalizer).
"""
function normalizer(S::AbstractMatrix)
    n = size(S, 2) ÷ 2
    return nullspace_gf2([S[:, n+1:end] S[:, 1:n]])
end

# Incremental row space, for repeated tests of `rank(vstack(cur, v)) > rank(cur)`.
struct RowSpace
    rows::Dict{Int,BitVector}  # leading index -> row
end
RowSpace() = RowSpace(Dict{Int,BitVector}())

function residual(R::RowSpace, v::AbstractVector)
    v = BitVector(v)
    while true
        f = findfirst(v)
        (f === nothing || !haskey(R.rows, f)) && return v
        v .⊻= R.rows[f]
    end
end

"""
    addrow!(R, v)

Add `v` to the row space; returns whether the rank grew.
"""
function addrow!(R::RowSpace, v::AbstractVector)
    w = residual(R, v)
    f = findfirst(w)
    f === nothing && return false
    R.rows[f] = w
    return true
end

function RowSpace(M::AbstractMatrix)
    R = RowSpace()
    for i in axes(M, 1)
        addrow!(R, view(M, i, :))
    end
    return R
end

Base.length(R::RowSpace) = length(R.rows)

"""
    logical_basis(S; first=[])

Symplectic logical basis `[(X_i, Z_i)]` of the stabilizer code with rows `S` (gf2.logical_basis):
independent representatives from `first` then the normalizer, paired by symplectic Gram-Schmidt.
"""
function logical_basis(S::AbstractMatrix; first=Vector{Bool}[])
    N = normalizer(S)
    pool = vcat([Vector{Bool}(v) for v in first], [N[i, :] for i in axes(N, 1)])
    cur = RowSpace(S)
    rS = length(cur)
    reps = Vector{Bool}[]
    for v in pool
        addrow!(cur, v) && push!(reps, Vector{Bool}(v))
    end
    iseven(length(cur) - rS) || error("odd number of logical representatives")
    pairs = Tuple{Vector{Bool},Vector{Bool}}[]
    while !isempty(reps)
        a = popfirst!(reps)
        j = findfirst(b -> symp_inner(a, b), reps)
        j === nothing && throw(ArgumentError("degenerate logical set (vector commutes with all others)"))
        b = popat!(reps, j)
        reps = map(reps) do c
            symp_inner(c, b) && (c = c .⊻ a)
            symp_inner(c, a) && (c = c .⊻ b)
            c
        end
        push!(pairs, (a, b))
    end
    return pairs
end

"""
    gf2_solve(A, b)

One solution of `A x = b` over GF(2), or `nothing` (builder.gf2_solve).
"""
function gf2_solve(A::AbstractMatrix, b::AbstractVector)
    m, n = size(A)
    M = Matrix{Bool}([isodd.(A) isodd.(b)])
    piv = Int[]
    r = 1
    for col in 1:n
        r > m && break
        p = findfirst(view(M, r:m, col))
        p === nothing && continue
        p += r - 1
        M[[r, p], :] = M[[p, r], :]
        for rr in 1:m
            if rr != r && M[rr, col]
                @views M[rr, :] .= M[rr, :] .⊻ M[r, :]
            end
        end
        push!(piv, col)
        r += 1
    end
    any(view(M, r:m, n + 1)) && return nothing
    x = falses(n)
    for (i, col) in enumerate(piv)
        x[col] = M[i, n+1]
    end
    return Vector{Bool}(x)
end
