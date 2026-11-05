# Bivariate bicycle codes in the conventions of Bravyi et al., arXiv:2308.07915.
# x = S_l ⊗ I_m, y = I_l ⊗ S_m, monomial x^a y^b <-> index a*m + b; H_X = [A | B], H_Z = [B^T | A^T].
# Qubits 0..n/2-1 are the L block, n/2..n-1 the R block. The same H comes out of QuantumClifford's
# BivariateBicycleViaCirculantMat (checked in the tests).

"""Cyclic shift S_k: row i has its 1 at column i + 1."""
shiftmat(k::Integer) = Bool[mod(j - i, k) == 1 for i in 0:k-1, j in 0:k-1]

matpow2(M, e) = e == 0 ? Matrix{Bool}([i == j for i in axes(M, 1), j in axes(M, 2)]) : isodd.(Int.(M)^e)

function monomial(l, m, a, b)
    X = kron(matpow2(shiftmat(l), mod(a, l)), I2(m))
    Y = kron(I2(l), matpow2(shiftmat(m), mod(b, m)))
    return isodd.(Int.(X) * Int.(Y))
end

poly(l, m, terms) = reduce((A, B) -> A .⊻ B, (monomial(l, m, a, b) for (a, b) in terms); init=falses(l * m, l * m))

# name => (l, m, A terms, B terms); A = x^3 + y + y^2, B = y^3 + x + x^2
const BB_CODES = Dict(
    "bb72" => (6, 6, [(3, 0), (0, 1), (0, 2)], [(0, 3), (1, 0), (2, 0)]),
    "gross" => (12, 6, [(3, 0), (0, 1), (0, 2)], [(0, 3), (1, 0), (2, 0)]),
)

"""
    BBCode(name="bb72")

Check matrices `HX`, `HZ`, `n`, `k` and symplectic stabilizer matrix `S` of a code in `BB_CODES`.
"""
struct BBCode
    name::String
    l::Int
    m::Int
    A::Matrix{Bool}
    B::Matrix{Bool}
    HX::Matrix{Bool}
    HZ::Matrix{Bool}
    n::Int
    k::Int
    S::Matrix{Bool}
end

function BBCode(name::AbstractString="bb72")
    l, m, at, bt = BB_CODES[name]
    A, B = Matrix{Bool}(poly(l, m, at)), Matrix{Bool}(poly(l, m, bt))
    HX = [A B]
    HZ = [permutedims(B) permutedims(A)]
    @assert !any(isodd.(Int.(HX) * Int.(permutedims(HZ))))
    n = 2 * l * m
    k = n - gf2rank(HX) - gf2rank(HZ)
    return BBCode(name, l, m, A, B, HX, HZ, n, k, css_stabilizers(HX, HZ))
end

x_checks(c::BBCode) = [PauliDict(q - 1 => 'X' for q in findall(c.HX[i, :])) for i in axes(c.HX, 1)]
z_checks(c::BBCode) = [PauliDict(q - 1 => 'Z' for q in findall(c.HZ[i, :])) for i in axes(c.HZ, 1)]
check_dicts(c::BBCode) = vcat(x_checks(c), z_checks(c))

"""
    css_logical_basis_with(code, x0, z0)

CSS symplectic logical basis `[(Xbar_i, Zbar_i)]` (symplectic vectors) whose first pair is the
anticommuting pair `(x0, z0)`.
"""
function css_logical_basis_with(code::BBCode, x0::Vector{Bool}, z0::Vector{Bool})
    n = code.n
    symp_inner(x0, z0) || throw(ArgumentError("x0 and z0 must anticommute"))
    NZ, NX = nullspace_gf2(code.HZ), nullspace_gf2(code.HX)
    KX = [vcat(NZ[i, :], falses(n)) for i in axes(NZ, 1)]
    KZ = [vcat(falses(n), NX[i, :]) for i in axes(NX, 1)]
    function indep(first, other, vs)
        out = [first]
        cur = RowSpace(code.S)
        addrow!(cur, first)
        for v in vs
            v = symp_inner(v, other) ? v .⊻ first : v
            addrow!(cur, v) && push!(out, Vector{Bool}(v))
        end
        return out
    end
    Xs, Zs = indep(x0, z0, KX), indep(z0, x0, KZ)
    @assert length(Xs) == length(Zs) == code.k
    return _pair_rest([(x0, z0)], Xs[2:end], Zs[2:end])
end

# Make the remaining X and Z representatives dual: Z' = M^{-T} Z with M_ab = <X_a, Z_b>.
function _pair_rest(firsts, Xr, Zr)
    k1 = length(Xr)
    M = Bool[symp_inner(a, b) for a in Xr, b in Zr]
    R, piv = rref([M I2(k1)])
    piv[1:k1] == 1:k1 || error("logical representatives are not independent")
    Minv = R[:, k1+1:end]
    Zm = permutedims(reduce(hcat, Zr))
    Zp = isodd.(Int.(permutedims(Minv)) * Int.(Zm))
    pairs = vcat(firsts, [(x, Vector{Bool}(Zp[i, :])) for (i, x) in enumerate(Xr)])
    for (i, (a, _)) in enumerate(pairs), (j, (_, b)) in enumerate(pairs)
        @assert symp_inner(a, b) == (i == j)
    end
    return pairs
end

logical_basis_with(code::BBCode, first_vecs=Vector{Bool}[]) = logical_basis(code.S; first=first_vecs)

"""
    term_perms(l, m, terms)

For each monomial term, the permutation π with row α -> column π[α] (0-based).
"""
term_perms(l, m, terms) = [[findfirst(monomial(l, m, a, b)[r, :]) - 1 for r in 1:l*m] for (a, b) in terms]

"""
    ordered_checks(code; order_x=0:5, order_z=0:5)

Checks with translation-invariant gate orders: X checks in group 0, Z checks in group 1.
Slots 0-2 are the A (X check) or B^T (Z check) terms on the L block, slots 3-5 the B or A^T terms
on the R block; `order_x`, `order_z` permute these slots.
"""
function ordered_checks(code::BBCode; order_x=0:5, order_z=0:5)
    l, m, at, bt = BB_CODES[code.name]
    half = l * m
    pA, pB = term_perms(l, m, at), term_perms(l, m, bt)
    invp(p) = invperm(p .+ 1) .- 1  # a transposed permutation matrix has row α's 1 at π^{-1}(α)
    ipA, ipB = invp.(pA), invp.(pB)
    xs, zs = Check[], Check[]
    for a in 0:half-1
        slots_x = vcat([p[a+1] for p in pA], [half + p[a+1] for p in pB])
        slots_z = vcat([p[a+1] for p in ipB], [half + p[a+1] for p in ipA])
        push!(xs, Check(PauliDict(q => 'X' for q in slots_x); order=[slots_x[t+1] for t in order_x], group=0, name="X$a"))
        push!(zs, Check(PauliDict(q => 'Z' for q in slots_z); order=[slots_z[t+1] for t in order_z], group=1, name="Z$a"))
    end
    for (c, i) in zip(xs, axes(code.HX, 1))
        @assert sort(collect(keys(c.pauli))) == findall(code.HX[i, :]) .- 1
    end
    for (c, i) in zip(zs, axes(code.HZ, 1))
        @assert sort(collect(keys(c.pauli))) == findall(code.HZ[i, :]) .- 1
    end
    return vcat(xs, zs)
end
