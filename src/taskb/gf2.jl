# GF(2) linear algebra on Matrix{Bool}.
#
# Conventions for all of Task B: a Pauli on n qubits is a row vector (x | z) of length 2n,
# the symplectic form is w(u, v) = u_x . v_z + u_z . v_x, and a Clifford acts on Pauli rows from the right,
# v -> v F. Phases are ignored at this level.

const GF2Mat = Matrix{Bool}

"""
    gf2(a)

`a` reduced mod 2 as a `Matrix{Bool}` (vectors become 1-row matrices).
"""
gf2(a::AbstractMatrix) = GF2Mat(isodd.(a))
gf2(a::AbstractVector) = GF2Mat(reshape(isodd.(a), 1, :))

"""
    gf2mul(a, b)

Matrix product over GF(2).
"""
gf2mul(a::AbstractVecOrMat, b::AbstractVecOrMat) = GF2Mat(isodd.(Int.(a) * Int.(b)))

"""
    rref(m) -> (r, pivots)

Reduced row echelon form over GF(2) and the (1-based) pivot columns. Same elimination order as the
reference, so pivots and the reduced matrix agree entry for entry.
"""
function rref(m::AbstractMatrix)
    r = gf2(m)
    rows, cols = size(r)
    pivots = Int[]
    i = 1
    for j in 1:cols
        i > rows && break
        p = findnext(view(r, :, j), i)
        p === nothing && continue
        if p != i
            for c in 1:cols
                r[i, c], r[p, c] = r[p, c], r[i, c]
            end
        end
        for o in 1:rows
            if o != i && r[o, j]
                @inbounds for c in j:cols
                    r[o, c] ⊻= r[i, c]
                end
            end
        end
        push!(pivots, j)
        i += 1
    end
    return r, pivots
end

"""
    gf2rank(m)

Rank over GF(2).
"""
gf2rank(m::AbstractMatrix) = isempty(m) ? 0 : length(rref(m)[2])

"""
    row_basis(m)

Basis of the row space: the nonzero rows of the RREF.
"""
function row_basis(m::AbstractMatrix)
    r, piv = rref(m)
    return r[1:length(piv), :]
end

"""
    nullspace_gf2(m)

Basis, as rows, of `{v : m v^T = 0}`.
"""
function nullspace_gf2(m::AbstractMatrix)
    cols = size(m, 2)
    r, piv = rref(m)
    free = setdiff(1:cols, piv)
    basis = zeros(Bool, length(free), cols)
    for (t, f) in enumerate(free)
        basis[t, f] = true
        for (i, p) in enumerate(piv)
            basis[t, p] = r[i, f]
        end
    end
    return basis
end

"""
    in_rowspace(rows, v)

Whether the row vector `v` lies in the row space of `rows`.
"""
in_rowspace(rows::AbstractMatrix, v::AbstractVector) = gf2rank(vcat(gf2(rows), gf2(v))) == gf2rank(rows)

"""
    solve_left(a, b)

A solution `x` of `x a = b` over GF(2) (`b` a vector or a matrix of row vectors), or `nothing`. Free variables
are set to zero, as in the reference, so the particular solution is the same.
"""
function solve_left(a::AbstractMatrix, b::AbstractVecOrMat)
    single = b isa AbstractVector
    b2 = gf2(b)
    at = permutedims(gf2(a))
    nvar = size(at, 2)
    r, piv = rref(hcat(at, permutedims(b2)))
    any(>(nvar), piv) && return nothing
    x = zeros(Bool, size(b2, 1), nvar)
    for (i, p) in enumerate(piv)
        x[:, p] = r[i, nvar+1:end]
    end
    return single ? vec(x) : x
end

"""
    gf2inv(m)

Inverse over GF(2); throws an `ArgumentError` if `m` is singular.
"""
function gf2inv(m::AbstractMatrix)
    n = size(m, 1)
    r, piv = rref(hcat(gf2(m), GF2Mat(I(n))))
    piv[1:min(n, end)] == 1:n || throw(ArgumentError("matrix not invertible over GF(2)"))
    return r[:, n+1:end]
end

"""
    omega(nq)

Symplectic form `[0 I; I 0]` for (x|z) row vectors on `nq` qubits.
"""
omega(nq::Integer) = GF2Mat([zeros(Bool, nq, nq) I(nq); I(nq) zeros(Bool, nq, nq)])

"""
    symp_inner(u, v)

Matrix of symplectic products between the rows of `u` and the rows of `v`.
"""
function symp_inner(u::AbstractVecOrMat, v::AbstractVecOrMat)
    u2, v2 = u isa AbstractVector ? gf2(u) : u, v isa AbstractVector ? gf2(v) : v
    nq = size(u2, 2) ÷ 2
    return gf2mul(u2, gf2mul(omega(nq), permutedims(v2)))
end

"""
    is_symplectic(f)

Whether `f` preserves the symplectic form, `f Ω fᵀ = Ω`.
"""
function is_symplectic(f::AbstractMatrix)
    om = omega(size(f, 1) ÷ 2)
    return gf2mul(gf2mul(f, om), permutedims(f)) == om
end

"""
    bitrows(m)

Rows of a 0/1 matrix as strings of '0'/'1' (the format of the exported reference files).
"""
bitrows(m::AbstractMatrix) = [join(Int.(m[i, :])) for i in 1:size(m, 1)]

"""
    frombits(rows)

Inverse of [`bitrows`](@ref).
"""
frombits(rows::AbstractVector) =
    isempty(rows) ? GF2Mat(undef, 0, 0) : GF2Mat(reduce(vcat, [permutedims([c == '1' for c in r]) for r in rows]))
