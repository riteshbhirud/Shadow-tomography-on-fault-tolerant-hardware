# Physical Clifford layers as 2n × 2n symplectic matrices (v -> v F), and gate lists that can be replayed on
# QuantumClifford tableaux to get the signs.
#
# Qubit indices in gate lists and permutations are 0-based, exactly as in the saved
# JSON files. "SWAPPERM" carries a whole permutation: the Pauli on qubit i moves to qubit perm[i].

"""
    Gate(name, qubits)

One gate of a physical gate list. `name` is one of `:SWAPPERM, :H, :S, :CZ, :CX, :X, :Y, :Z`; `qubits` are
0-based (for `:SWAPPERM` the full permutation).
"""
struct Gate
    name::Symbol
    qubits::Vector{Int}
end

Base.:(==)(a::Gate, b::Gate) = a.name == b.name && a.qubits == b.qubits
Base.hash(g::Gate, h::UInt) = hash(g.qubits, hash(g.name, h))

"""
    perm_op(perm)

Symplectic matrix of the qubit permutation i -> perm[i] (0-based).
"""
function perm_op(perm::AbstractVector{<:Integer})
    n = length(perm)
    p = zeros(Bool, n, n)
    for i in 1:n
        p[i, perm[i]+1] = true
    end
    return GF2Mat([p zeros(Bool, n, n); zeros(Bool, n, n) p])
end

"""
    h_all(n)

Symplectic matrix of transversal H.
"""
h_all(n::Integer) = omega(n)

"""
    z_diag(s)

`[I S; 0 I]`: CZ on the off-diagonal pairs of the symmetric `s`, S (√Z) on its diagonal.
"""
function z_diag(s::AbstractMatrix)
    s = gf2(s)
    @assert s == permutedims(s)
    n = size(s, 1)
    return GF2Mat([I(n) s; zeros(Bool, n, n) I(n)])
end

"""
    x_diag(t)

`[I 0; T I]` = Hⁿ U(T) Hⁿ.
"""
function x_diag(t::AbstractMatrix)
    t = gf2(t)
    @assert t == permutedims(t)
    n = size(t, 1)
    return GF2Mat([I(n) zeros(Bool, n, n); t I(n)])
end

"""
    matching_matrix(perm)

Symmetric S with S[i, perm[i]] = 1 for an involution: CZ across 2-cycles, √Z on fixed points.
"""
function matching_matrix(perm::AbstractVector{<:Integer})
    n = length(perm)
    @assert all(perm[perm[i]+1] == i - 1 for i in 1:n) "not an involution"
    s = zeros(Bool, n, n)
    for i in 1:n
        s[i, perm[i]+1] = true
    end
    return s
end

gates_perm(perm) = [Gate(:SWAPPERM, collect(Int, perm))]
gates_h_all(n) = [Gate(:H, [q]) for q in 0:n-1]

"""
    gates_z_diag(s)

Gate list of [`z_diag`](@ref): S on diagonal entries, CZ on off-diagonal pairs (row-major upper triangle).
"""
function gates_z_diag(s::AbstractMatrix)
    n = size(s, 1)
    g = Gate[]
    for i in 1:n
        s[i, i] && push!(g, Gate(:S, [i - 1]))
        for j in i+1:n
            s[i, j] && push!(g, Gate(:CZ, [i - 1, j - 1]))
        end
    end
    return g
end

gates_x_diag(t) = [gates_h_all(size(t, 1)); gates_z_diag(t); gates_h_all(size(t, 1))]

"""
    gates_to_symplectic(gates, n)

The symplectic matrix of a gate list, derived gate by gate (independent of how the list was built).
"""
function gates_to_symplectic(gates, n::Integer)
    f = GF2Mat(I(2n))
    for g in gates
        m = GF2Mat(I(2n))
        q = g.qubits .+ 1
        if g.name == :SWAPPERM
            m = perm_op(g.qubits)
        elseif g.name == :H
            m[q[1], q[1]] = m[n+q[1], n+q[1]] = false
            m[q[1], n+q[1]] = m[n+q[1], q[1]] = true
        elseif g.name == :S
            m[q[1], n+q[1]] = true
        elseif g.name == :CZ
            m[q[1], n+q[2]] = true
            m[q[2], n+q[1]] = true
        elseif g.name == :CX
            m[q[1], q[2]] = true          # X_c -> X_c X_t
            m[n+q[2], n+q[1]] = true      # Z_t -> Z_c Z_t
        elseif g.name in (:X, :Y, :Z)
        else
            error("unknown gate $(g.name)")
        end
        f = gf2mul(f, m)
    end
    return f
end

"""
    group_automorphisms(l, m)

All automorphisms of Z_l × Z_m, as pairs (image of x, image of y), in the reference enumeration order.
"""
function group_automorphisms(l::Integer, m::Integer)
    elems = [(a, b) for a in 0:l-1 for b in 0:m-1]
    out = Tuple{Tuple{Int,Int},Tuple{Int,Int}}[]
    for (xa, xb) in elems
        mod(l * xb, m) == 0 || continue
        for (ya, yb) in elems
            mod(m * ya, l) == 0 || continue
            img = Set((mod(xa * a + ya * b, l), mod(xb * a + yb * b, m)) for (a, b) in elems)
            length(img) == l * m && push!(out, ((xa, xb), (ya, yb)))
        end
    end
    return out
end

"""
    affine_perm(l, m, phi, c, cprime, swap)

Qubit permutation L(α) -> L(φ(α) + c), R(α) -> R(φ(α) + c′); with `swap`, L goes to R and R to L.
"""
function affine_perm(l, m, phi, c, cprime, swap::Bool)
    (xa, xb), (ya, yb) = phi
    N = l * m
    perm = zeros(Int, 2N)
    for a in 0:l-1, b in 0:m-1
        pa = mod(xa * a + ya * b, l)
        pb = mod(xb * a + yb * b, m)
        iL = mod(pa + c[1], l) * m + mod(pb + c[2], m)
        iR = mod(pa + cprime[1], l) * m + mod(pb + cprime[2], m)
        src = a * m + b
        if swap
            perm[src+1], perm[N+src+1] = N + iL, iR
        else
            perm[src+1], perm[N+src+1] = iL, N + iR
        end
    end
    return perm
end

_rowset(h::AbstractMatrix) = Set(BitVector(h[i, :]) for i in axes(h, 1))

"""
    search_bb_automorphisms(hx, hz, l, m) -> (perm_auts, zx_dualities)

Exhaustive search over the affine family with c = 0 (the l·m shifts give the rest). Automorphisms satisfy
H_X P ~ H_X and H_Z P ~ H_Z as row sets; ZX-dualities swap the two.
"""
function search_bb_automorphisms(hx, hz, l::Integer, m::Integer)
    hx_set, hz_set = _rowset(hx), _rowset(hz)
    auts, duals = Vector{Int}[], Vector{Int}[]
    elems = [(a, b) for a in 0:l-1 for b in 0:m-1]
    for phi in group_automorphisms(l, m), swap in (false, true), cp in elems
        perm = affine_perm(l, m, phi, (0, 0), cp, swap)
        inv = sortperm(perm)   # column j of the permuted matrix is old column perm^{-1}(j)
        sx = _rowset(hx[:, inv])
        if sx == hx_set
            _rowset(hz[:, inv]) == hz_set && push!(auts, perm)
        elseif sx == hz_set
            _rowset(hz[:, inv]) == hx_set && push!(duals, perm)
        end
    end
    return auts, duals
end

# --- signs through QuantumClifford -----------------------------------------------------------------------------

"""
    expand_gates(gates)

The gate list with every `:SWAPPERM` written as a SWAP network (one SWAP per transposition of each
cycle).
"""
function expand_gates(gates)
    out = Gate[]
    for g in gates
        if g.name != :SWAPPERM
            push!(out, g)
            continue
        end
        perm = g.qubits
        seen = falses(length(perm))
        for c0 in 0:length(perm)-1
            if seen[c0+1] || perm[c0+1] == c0
                seen[c0+1] = true
                continue
            end
            cyc = [c0]
            seen[c0+1] = true
            j = perm[c0+1]
            while j != c0
                push!(cyc, j)
                seen[j+1] = true
                j = perm[j+1]
            end
            append!(out, [Gate(:SWAP, [c0, cj]) for cj in cyc[2:end]])
        end
    end
    return out
end

function _qc_op(g::Gate)
    q = g.qubits .+ 1
    g.name == :SWAP && return sSWAP(q[1], q[2])
    g.name == :H && return sHadamard(q[1])
    g.name == :S && return sPhase(q[1])
    g.name == :CZ && return sCPHASE(q[1], q[2])
    g.name == :CX && return sCNOT(q[1], q[2])
    g.name == :X && return sX(q[1])
    g.name == :Y && return sY(q[1])
    g.name == :Z && return sZ(q[1])
    error("unknown gate $(g.name)")
end

"""
    apply_gates!(state, gates)

Apply a gate list to a QuantumClifford tableau (`Stabilizer`, `MixedDestabilizer`, ...).
"""
function apply_gates!(st, gates)
    for g in expand_gates(gates)
        apply!(st, _qc_op(g))
    end
    return st
end

"""
    rowpauli(bits, n)

The Hermitian Pauli with (x|z) bits `bits` and sign +, as a QuantumClifford `PauliOperator` (Y where x = z = 1).
"""
rowpauli(bits::AbstractVector, n::Integer) = PauliOperator(0x0, Vector{Bool}(bits[1:n]), Vector{Bool}(bits[n+1:2n]))

_rows_tableau(rows::AbstractMatrix, n) =
    Stabilizer(zeros(UInt8, size(rows, 1)), Matrix{Bool}(rows[:, 1:n]), Matrix{Bool}(rows[:, n+1:2n]))

_sign(p::PauliOperator) = p.phase[] == 0x0 ? 1 : p.phase[] == 0x2 ? -1 : error("non-Hermitian Pauli")

"""
    css_group_sign(bits, n)

Sign of the element of the CSS stabilizer group ⟨+X-type rows, +Z-type rows⟩ with bits `bits`, relative to the
Hermitian Pauli with those bits: X(x) Z(z) = (-i)^{|x∧z|} Y-string, and |x∧z| is even, so (-1)^{|x∧z|/2}.
"""
function css_group_sign(bits::AbstractVector, n::Integer)
    m = count(bits[1:n] .& bits[n+1:2n])
    iseven(m) || error("not an element of a CSS stabilizer group (X and Z parts anticommute)")
    return isodd(m ÷ 2) ? -1 : 1
end

"""
    stabilizer_sign_flips(stabs, gates, n) -> (images, flip)

Conjugate every (sign +) stabilizer row through `gates` with QuantumClifford; `flip[i]` is true iff the image is
minus an element of S.
"""
function stabilizer_sign_flips(stabs::AbstractMatrix, gates, n::Integer)
    st = apply_gates!(_rows_tableau(stabs, n), gates)
    r = size(stabs, 1)
    imgs = zeros(Bool, r, 2n)
    flip = zeros(Bool, r)
    for i in 1:r
        p = st[i]
        b = [xbit(p); zbit(p)]
        imgs[i, :] = b
        flip[i] = _sign(p) != css_group_sign(b, n)
    end
    return imgs, flip
end

"""
    pauli_correction(stabs, gates, n) -> (gates, c)

A Pauli C (as X/Y/Z gates) such that C·U maps every stabilizer to +stabilizer, from ω(c, imgᵢ) = flipᵢ solved
over GF(2) with the same particular solution as the reference.
"""
function pauli_correction(stabs::AbstractMatrix, gates, n::Integer)
    imgs, flip = stabilizer_sign_flips(stabs, gates, n)
    any(flip) || return Gate[], zeros(Bool, 2n)
    sw = hcat(imgs[:, n+1:end], imgs[:, 1:n])      # omega(c, r) = c . swap(r)
    c = solve_left(permutedims(sw), flip)
    c === nothing && error("inconsistent sign pattern: no Pauli correction exists")
    out = Gate[]
    for q in 1:n
        cx, cz = c[q], c[n+q]
        (cx || cz) && push!(out, Gate(cx && cz ? :Y : cx ? :X : :Z, [q - 1]))
    end
    return out, c
end
