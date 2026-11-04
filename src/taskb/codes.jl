# Stabilizer codes for Task B: BB codes in the Bravyi et al. (2308.07915) main-text convention, La-cross,
# logical bases, and the logical (symplectic) action of a physical Clifford.

"""
    StabilizerCode

Stabilizer code with a symplectic logical basis. All rows are (x|z) over GF(2).

- `stabs`: r × 2n stabilizer rows (may be redundant);
- `lx`, `lz`: k × 2n logical X̄ᵢ and Z̄ᵢ;
- `meta`: construction data (`:hx`, `:hz`, `:l`, `:m`, `:basis`, ...).
"""
struct StabilizerCode
    name::String
    n::Int
    stabs::GF2Mat
    lx::GF2Mat
    lz::GF2Mat
    meta::OrderedDict{Symbol,Any}
end

logical_k(c::StabilizerCode) = size(c.lx, 1)

"""
    logicals(code)

The 2k × 2n matrix X̄₁..X̄ₖ, Z̄₁..Z̄ₖ; row j is logical coordinate j.
"""
logicals(c::StabilizerCode) = vcat(c.lx, c.lz)

"""
    is_css_basis(code)

CSS stabilizers (every row pure X or pure Z) and a CSS logical basis (X̄ᵢ physically X-type, Z̄ᵢ Z-type).
"""
function is_css_basis(c::StabilizerCode)
    n = c.n
    s = c.stabs
    rows_pure = !any(any(s[i, 1:n]) && any(s[i, n+1:end]) for i in axes(s, 1))
    return rows_pure && !any(c.lx[:, n+1:end]) && !any(c.lz[:, 1:n])
end

"""
    validate(code; require_css=true)

Check commutation, the rank count r + k = n, and that the logical basis is symplectic and independent mod S.
The shadow pipeline needs a CSS logical basis, so a non-CSS one is rejected unless `require_css=false`.
"""
function validate(c::StabilizerCode; require_css::Bool=true)
    n, k = c.n, logical_k(c)
    s = c.stabs
    require_css && !is_css_basis(c) &&
        error("$(c.name): logical basis (or stabilizer set) is not CSS: every Xbar_i must be physically X-type " *
              "and every Zbar_i physically Z-type")
    any(symp_inner(s, s)) && error("stabilizers do not commute")
    r = gf2rank(s)
    r + k == n || error("rank(S)=$r, k=$k, n=$n")
    L = logicals(c)
    any(symp_inner(L, s)) && error("logicals do not commute with stabilizers")
    symp_inner(L, L) == omega(k) || error("logical basis not symplectic")
    gf2rank(vcat(s, L)) == r + 2k || error("logicals not independent mod S")
    return nothing
end

"""
    logical_coords(code, p; check=true)

Coordinates (a|b) ∈ GF(2)^{2k} of normalizer elements `p` (rows), with P ~ ∏ X̄^a Z̄^b.
"""
function logical_coords(c::StabilizerCode, p::AbstractMatrix; check::Bool=true)
    p = gf2(p)
    coords = hcat(symp_inner(p, c.lz), symp_inner(p, c.lx))  # coefficient on Xbar_i is w(p, Zbar_i)
    if check
        any(symp_inner(p, c.stabs)) && error("not in normalizer")
        resid = p .⊻ gf2mul(coords, logicals(c))
        srow = row_basis(c.stabs)
        rs = gf2rank(srow)
        for i in axes(resid, 1)
            gf2rank(vcat(srow, resid[i:i, :])) == rs || error("residual not a stabilizer")
        end
    end
    return coords
end

"""
    preserves_stabilizer_group(code, f)

Whether the physical symplectic `f` maps the stabilizer group onto itself.
"""
function preserves_stabilizer_group(c::StabilizerCode, f::AbstractMatrix)
    srow = get!(() -> row_basis(c.stabs), c.meta, :_stab_row_basis)
    return gf2rank(vcat(srow, gf2mul(c.stabs, f))) == size(srow, 1)
end

"""
    logical_action(code, f; check=true)

The 2k × 2k symplectic matrix Λ with coords ↦ coords Λ for the physical symplectic `f`.
"""
function logical_action(c::StabilizerCode, f::AbstractMatrix; check::Bool=true)
    check && !preserves_stabilizer_group(c, f) && error("physical op does not preserve the stabilizer group")
    lam = logical_coords(c, gf2mul(logicals(c), f); check)
    check && !is_symplectic(lam) && error("logical action not symplectic")
    return lam
end

"""
    css_stabs(hx, hz)

Stabilizer rows (x|z) of the CSS code with check matrices `hx`, `hz`.
"""
function css_stabs(hx::AbstractMatrix, hz::AbstractMatrix)
    n = size(hx, 2)
    return GF2Mat([gf2(hx) zeros(Bool, size(hx, 1), n); zeros(Bool, size(hz, 1), n) gf2(hz)])
end

"""
    generic_css_logicals(hx, hz) -> (xl, zl)

Symplectic CSS logical basis by Gaussian elimination: X̄ from ker(H_Z) mod row(H_X), Z̄ from ker(H_X) mod
row(H_Z), then Z̄ rotated so that X̄ Z̄ᵀ = I.
"""
function generic_css_logicals(hx::AbstractMatrix, hz::AbstractMatrix)
    function complement(kernel, rows)
        cur = row_basis(rows)
        chosen = Vector{Bool}[]
        kb = row_basis(kernel)
        for i in axes(kb, 1)
            t = vcat(cur, kb[i:i, :])
            if gf2rank(t) > size(cur, 1)
                cur = row_basis(t)
                push!(chosen, kb[i, :])
            end
        end
        return GF2Mat(reduce(vcat, permutedims.(chosen)))
    end
    xl = complement(nullspace_gf2(hz), hx)
    zl = complement(nullspace_gf2(hx), hz)
    zl = gf2mul(permutedims(gf2inv(gf2mul(xl, permutedims(zl)))), zl)
    @assert gf2mul(xl, permutedims(zl)) == I
    return xl, zl
end

"""
    make_css_code(name, hx, hz, xl=nothing, zl=nothing; meta...)

Validated CSS code; the generic logical basis is used when `xl`, `zl` are not given.
"""
function make_css_code(name, hx, hz, xl=nothing, zl=nothing; meta...)
    hx, hz = gf2(hx), gf2(hz)
    n = size(hx, 2)
    if xl === nothing
        xl, zl = generic_css_logicals(hx, hz)
    end
    k = size(xl, 1)
    lx = GF2Mat(hcat(xl, zeros(Bool, k, n)))
    lz = GF2Mat(hcat(zeros(Bool, k, n), zl))
    code = StabilizerCode(name, n, css_stabs(hx, hz), lx, lz, OrderedDict{Symbol,Any}(meta..., :hx => hx, :hz => hz))
    validate(code)
    return code
end

# --- BB codes -------------------------------------------------------------------------------------------------

"""
    BBRing(l, m)

The group algebra F₂[Z_l × Z_m]. Monomial x^a y^b has (0-based) index a*m + b.
"""
struct BBRing
    l::Int
    m::Int
end

ring_size(r::BBRing) = r.l * r.m
mono_index(r::BBRing, a, b) = mod(a, r.l) * r.m + mod(b, r.m)   # 0-based
mono(r::BBRing, i) = divrem(i, r.m)                              # 0-based index -> (a, b)

function ring_vec(r::BBRing, monos)
    v = zeros(Bool, ring_size(r))
    for (a, b) in monos
        v[mono_index(r, a, b)+1] ⊻= true
    end
    return v
end

"""
    ring_shift(r, v, a, b)

Multiplication of `v` by the monomial x^a y^b.
"""
ring_shift(r::BBRing, v::AbstractVector, a, b) = vec(permutedims(circshift(permutedims(reshape(v, r.m, r.l)), (a, b))))

function ring_transpose(r::BBRing, v::AbstractVector)
    out = zeros(Bool, length(v))
    for i in findall(v)
        a, b = mono(r, i - 1)
        out[mono_index(r, -a, -b)+1] = true
    end
    return out
end

"""
    BB_PARAMS

(l, m, A monomials, B monomials) of BB72 and the gross code; A = x³ + y + y², B = y³ + x + x².
"""
const BB_PARAMS = OrderedDict(
    "BB72" => (6, 6, [(3, 0), (0, 1), (0, 2)], [(0, 3), (1, 0), (2, 0)]),
    "gross" => (12, 6, [(3, 0), (0, 1), (0, 2)], [(0, 3), (1, 0), (2, 0)]),
)

# QuantumClifford's term format: a monomial is a pure power of x or of y
_qc_terms(monos) = [a == 0 ? (:y, b) : (@assert(b == 0); (:x, a)) for (a, b) in monos]

"""
    bb_matrices(l, m, amon, bmon) -> (ring, A, B, H_X, H_Z)

Check matrices H_X = [A|B], H_Z = [Bᵀ|Aᵀ] with x = S_l ⊗ I_m. The matrices come from QuantumClifford's
`BivariateBicycleViaCirculantMat`, which uses the same convention (checked against our own construction in
the tests).
"""
function bb_matrices(l, m, amon, bmon)
    c = QuantumClifford.ECC.BivariateBicycleViaCirculantMat(l, m, _qc_terms(amon), _qc_terms(bmon))
    hx = gf2(QuantumClifford.ECC.parity_matrix_x(c))
    hz = gf2(QuantumClifford.ECC.parity_matrix_z(c))
    N = l * m
    return BBRing(l, m), hx[:, 1:N], hx[:, N+1:end], hx, hz
end

"""
    bb_matrices_direct(l, m, amon, bmon)

Same as [`bb_matrices`](@ref) but built directly from shift matrices, S_l[i, i+1] = 1 (used to check the
QuantumClifford convention).
"""
function bb_matrices_direct(l, m, amon, bmon)
    shift(n) = GF2Mat([j == mod(i + 1, n) for i in 0:n-1, j in 0:n-1])
    polymat(monos) = reduce((p, q) -> p .⊻ q,
        [gf2(kron(Int.(shift(l))^mod(a, l), Int.(shift(m))^mod(b, m))) for (a, b) in monos])
    a, b = polymat(amon), polymat(bmon)
    return BBRing(l, m), a, b, hcat(a, b), hcat(permutedims(b), permutedims(a))
end

"""
    gross_fgh_logicals(ring)

The published Bravyi et al. (Supp.) f, g, h logical basis of the gross code; qubits 1–6 unprimed, 7–12 primed.
"""
function gross_fgh_logicals(ring::BBRing)
    f = ring_vec(ring, [(0, 0), (1, 0), (2, 0), (3, 0), (6, 0), (7, 0), (8, 0), (9, 0), (1, 3), (5, 3), (7, 3), (11, 3)])
    g = ring_vec(ring, [(1, 0), (2, 1), (0, 2), (1, 2), (2, 3), (0, 4)])
    h = ring_vec(ring, [(0, 0), (0, 1), (1, 1), (0, 2), (0, 3), (1, 3)])
    ns = [(0, 0), (0, 1), (2, 1), (2, 5), (3, 2), (4, 0)]
    ms = [(0, 1), (0, 5), (1, 1), (0, 0), (4, 0), (5, 2)]
    return _fgh_basis(ring, f, g, h, ns, ms)
end

function _fgh_basis(ring, f, g, h, ns, ms)
    zero = zeros(Bool, ring_size(ring))
    xl = [[ring_shift(ring, f, a, b); zero] for (a, b) in ns]                                   # X(alpha f, 0)
    append!(xl, [[ring_shift(ring, g, a, b); ring_shift(ring, h, a, b)] for (a, b) in ns])     # X(alpha g, alpha h)
    ht, gt, ft = ring_transpose(ring, h), ring_transpose(ring, g), ring_transpose(ring, f)
    zl = [[ring_shift(ring, ht, a, b); ring_shift(ring, gt, a, b)] for (a, b) in ms]           # Z(beta h^T, beta g^T)
    append!(zl, [[zero; ring_shift(ring, ft, a, b)] for (a, b) in ms])                          # Z(0, beta f^T)
    return GF2Mat(reduce(vcat, permutedims.(xl))), GF2Mat(reduce(vcat, permutedims.(zl)))
end

function _greedy_independent(vecs, base, count)
    cur = row_basis(base)
    chosen = Int[]
    for (i, v) in enumerate(vecs)
        t = vcat(cur, gf2(v))
        if gf2rank(t) > size(cur, 1)
            cur = row_basis(t)
            push!(chosen, i)
            length(chosen) == count && return chosen
        end
    end
    return nothing
end

"""
    fgh_basis_from(ring, A, B, f, g, h) -> (xl, zl, info) or nothing

The deterministic part of the Bravyi-style block-basis search for given f and (g, h): pick monomial
shifts greedily so each block has rank k/2, then rotate the Z̄ so that the commutation matrix is I. The random
search for f, g, h itself used numpy's RNG; its result is loaded from `test/reference/taskb/rng_inputs.json`.
"""
function fgh_basis_from(ring::BBRing, a, b, f, g, h)
    N = ring_size(ring)
    hx = hcat(a, b)
    hz = hcat(permutedims(b), permutedims(a))
    k = 2N - gf2rank(hx) - gf2rank(hz)
    monos = [mono(ring, i) for i in 0:N-1]
    zero = zeros(Bool, N)
    xs_un = [[ring_shift(ring, f, al...); zero] for al in monos]
    xs_pr = [[ring_shift(ring, g, al...); ring_shift(ring, h, al...)] for al in monos]
    ns = _greedy_independent(xs_un, hx, k ÷ 2)
    ns === nothing && return nothing
    ns_p = _greedy_independent(xs_pr, vcat(hx, gf2.(xs_un[ns])...), k ÷ 2)
    ns_p === nothing && return nothing
    ht, gt, ft = ring_transpose(ring, h), ring_transpose(ring, g), ring_transpose(ring, f)
    zs_un = [[ring_shift(ring, ht, al...); ring_shift(ring, gt, al...)] for al in monos]
    zs_pr = [[zero; ring_shift(ring, ft, al...)] for al in monos]
    ms = _greedy_independent(zs_un, hz, k ÷ 2)
    ms === nothing && return nothing
    ms_p = _greedy_independent(zs_pr, vcat(hz, gf2.(zs_un[ms])...), k ÷ 2)
    ms_p === nothing && return nothing
    xl = GF2Mat(reduce(vcat, permutedims.([xs_un[ns]; xs_pr[ns_p]])))
    zl = GF2Mat(reduce(vcat, permutedims.([zs_un[ms]; zs_pr[ms_p]])))
    c = gf2mul(xl, permutedims(zl))
    zl = gf2mul(permutedims(gf2inv(c)), zl)   # c is block diagonal, so the blocks survive the rotation
    h2 = k ÷ 2
    info = OrderedDict{String,Any}(
        "f" => [collect(mono(ring, i - 1)) for i in findall(f)],
        "g" => [collect(mono(ring, i - 1)) for i in findall(g)],
        "h" => [collect(mono(ring, i - 1)) for i in findall(h)],
        "n" => [collect(monos[i]) for i in ns],
        "n_primed" => [collect(monos[i]) for i in ns_p],
        "m" => [collect(monos[i]) for i in ms],
        "m_primed" => [collect(monos[i]) for i in ms_p],
        "commutation_block_diagonal" => !any(c[1:h2, h2+1:end]) && !any(c[h2+1:end, 1:h2]),
    )
    return xl, zl, info
end

const TASKB_REFDIR = joinpath(dirname(dirname(@__DIR__)), "test", "reference", "taskb")

"""
    rng_inputs()

The stored RNG-derived inputs (exported from the original implementation).
"""
rng_inputs() = JSON.parsefile(joinpath(TASKB_REFDIR, "rng_inputs.json"))

"""
    bb72_fgh_search(seed)

f, g, h found by the seeded search with `seed` (seed 0 gives the named BB72 basis), or `nothing` where
the search found no f or (g, h).
"""
function bb72_fgh_search(seed::Integer)
    d = rng_inputs()["bb72_fgh_search"][string(seed)]
    v(s) = s === nothing ? nothing : [c == '1' for c in s]
    return v(d["f"]), v(d["g"]), v(d["h"])
end

"""
    bb72_fgh_basis(seed=0)

The Bravyi-style BB72 block basis for search seed `seed`, or `nothing` if the search fails.
"""
function bb72_fgh_basis(seed::Integer=0)
    l, m, amon, bmon = BB_PARAMS["BB72"]
    ring, a, b, _, _ = bb_matrices(l, m, amon, bmon)
    f, g, h = bb72_fgh_search(seed)
    (f === nothing || g === nothing) && return nothing
    return fgh_basis_from(ring, a, b, f, g, h)
end

"""
    make_bb_code(name, basis="fgh")

BB72 or the gross code with the f,g,h block basis (default) or the generic CSS basis.
"""
function make_bb_code(name::AbstractString, basis::AbstractString="fgh")
    l, m, amon, bmon = BB_PARAMS[name]
    ring, a, b, hx, hz = bb_matrices(l, m, amon, bmon)
    meta = OrderedDict{Symbol,Any}(:l => l, :m => m, :A => amon, :B => bmon, :ring => ring, :Amat => a, :Bmat => b)
    if basis == "fgh" && name == "gross"
        xl, zl = gross_fgh_logicals(ring)
        meta[:basis] = "Bravyi et al. f,g,h basis (Supp. eq:logicalpaulis), qubits 1-6 unprimed, 7-12 primed"
    elseif basis == "fgh"
        res = bb72_fgh_basis(0)
        res === nothing && error("fgh basis search failed")
        xl, zl, info = res
        meta[:basis] = "computed Bravyi-style f,g,h block basis (qubits 1-6 unprimed, 7-12 primed)"
        meta[:fgh_info] = info
    else
        xl, zl = generic_css_logicals(hx, hz)
        meta[:basis] = "generic CSS basis (Gaussian elimination)"
    end
    return make_css_code(name, hx, hz, xl, zl; meta...)
end

# --- La-cross -------------------------------------------------------------------------------------------------

"""
    lacross_seed(kk, nn)

Open-boundary La-cross seed matrix H_{k,n} (Di Fini et al. 2610.12685): row a has ones at a, a+1, a+k.
"""
function lacross_seed(kk::Integer, nn::Integer)
    mm = nn - kk
    h = zeros(Bool, mm, nn)
    for a in 0:mm-1, j in (a, a + 1, a + kk)
        h[a+1, j+1] ⊻= true
    end
    return h
end

"""
    make_lacross_code(kk, nn)

La-cross code: hypergraph product of the seed H_{k,n} with itself, generic CSS logical basis.
"""
function make_lacross_code(kk::Integer, nn::Integer)
    h = lacross_seed(kk, nn)
    mm, n1 = size(h)
    i_n, i_m = I(n1), I(mm)
    hx = hcat(kron(h, i_n), kron(i_m, permutedims(h)))
    hz = hcat(kron(i_n, h), kron(permutedims(h), i_m))
    return make_css_code("LaCross[k=$kk,n=$nn]", hx, hz; seed=h, kk=kk, nn=nn)
end

"""
    min_weight_logical_upto(h, other_logicals, wmax=6) -> (w, support) or (nothing, nothing)

Smallest weight w ≤ `wmax` of x with h x = 0 and x · other_logicalsᵀ ≠ 0, by meet-in-the-middle over the
syndromes of column subsets of size ≤ ⌈wmax/2⌉. Exhaustive. The support is 0-based.
"""
function min_weight_logical_upto(h::AbstractMatrix, other::AbstractMatrix, wmax::Integer=6)
    h = gf2(h)
    n = size(h, 2)
    half = (wmax + 1) ÷ 2
    cols = [BitVector(h[:, j]) for j in 1:n]
    lg = gf2(other)
    buckets = OrderedDict{BitVector,Vector{Vector{Int}}}()
    order = BitVector[]
    for sz in 0:half, sub in _combinations(n, sz)
        s = falses(size(h, 1))
        for c in sub
            s .⊻= cols[c+1]
        end
        haskey(buckets, s) || push!(order, s)
        push!(get!(buckets, s, Vector{Int}[]), sub)
    end
    best_w, best_s = nothing, nothing
    for key in order
        subs = buckets[key]
        length(subs) < 2 && continue
        for i in eachindex(subs), j in i+1:lastindex(subs)
            supp = symdiff(subs[i], subs[j])
            w = length(supp)
            (w == 0 || w > wmax || (best_w !== nothing && w >= best_w)) && continue
            x = zeros(Bool, n)
            x[supp.+1] .= true
            if any(gf2mul(permutedims(x), permutedims(lg)))
                best_w, best_s = w, sort(supp)
            end
        end
    end
    return best_w, best_s
end

# 0-based combinations of 0:n-1 of size k in lexicographic order (itertools.combinations)
function _combinations(n::Integer, k::Integer)
    k == 0 && return [Int[]]
    out = Vector{Int}[]
    c = collect(0:k-1)
    while true
        push!(out, copy(c))
        i = k
        while i >= 1 && c[i] == n - k + i - 1
            i -= 1
        end
        i == 0 && return out
        c[i] += 1
        for j in i+1:k
            c[j] = c[j-1] + 1
        end
    end
end
