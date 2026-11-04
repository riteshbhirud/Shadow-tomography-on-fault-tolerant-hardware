# Single-logical measurements on BB codes by gauging.
# `bb_setup` fixes the logical basis shared by every BB72 oracle setting
# (handover/bases/bb72_oracle_basis.json).

"""Qubit permutation (0-based images) of the shift automorphism x^a y^b on both blocks."""
function shift_perm(code::BBCode, a::Integer, b::Integer)
    l, m = code.l, code.m
    half = l * m
    perm = zeros(Int, code.n)
    for blk in 0:1, i in 0:l-1, j in 0:m-1
        perm[blk*half+i*m+j+1] = blk * half + mod(i + a, l) * m + mod(j + b, m)
    end
    return perm
end

"""
    orbit(code, qubits)

Distinct supports of `qubits` under all shift automorphisms, sorted.
"""
function orbit(code::BBCode, qubits)
    out = Set{Vector{Int}}()
    for a in 0:code.l-1, b in 0:code.m-1
        p = shift_perm(code, a, b)
        push!(out, sort([p[q+1] for q in qubits]))
    end
    return sort!(collect(out))
end

const _BB_SETUP = Dict{String,Any}()

"""
    bb_setup(name="bb72")

Code, measured logicals and common logical basis: weight-6 X and Z logicals overlapping
once with minimal union are (Xbar_0, Zbar_0), and Ybar_0 is their phase-free product. Returns a
named tuple `(code, L, bases, basis, weights, Ypair, min_weight_found)`.

The starting low-weight logicals come from a seeded random search (stored in
`data/python_rng.json`); the orbit
search and basis completion are recomputed.
"""
function bb_setup(name::AbstractString="bb72")
    get!(_BB_SETUP, name) do
        name == "bb72" || error("only bb72 has exported starting logicals")
        code = BBCode(name)
        n = code.n
        mw = python_rng()["bb72_min_weight_logicals"]
        sx, sz = Vector{Int}(mw["X"]["support"]), Vector{Int}(mw["Z"]["support"])
        @assert _is_logical(code, PauliDict(q => 'X' for q in sx)) && _is_logical(code, PauliDict(q => 'Z' for q in sz))
        xs, zs = orbit(code, sx), orbit(code, sz)
        best = nothing
        for X in xs, Z in zs
            if isodd(length(X ∩ Z))
                u = length(X ∪ Z)
                (best === nothing || u < best[1]) && (best = (u, X, Z))
            end
        end
        best === nothing && error("no anticommuting pair among the shift orbits")
        _, Xy, Zy = best
        # Overlap = 1 mod 4 makes the phase-free product equal i Xbar_0 Zbar_0 = Ybar_0.
        mod(length(Xy ∩ Zy), 4) == 1 || error("X/Z overlap is not 1 mod 4")
        Xyd = PauliDict(q => 'X' for q in Xy)
        Zyd = PauliDict(q => 'Z' for q in Zy)
        Ldict = Dict('X' => Xyd, 'Z' => Zyd, 'Y' => pauli_product(n, [Xyd, Zyd]))
        pairs = css_logical_basis_with(code, dict_to_vec(Xyd, n), dict_to_vec(Zyd, n))
        basis = LogicalBasis([(vec_to_dict(a), vec_to_dict(b)) for (a, b) in pairs])
        @assert basis[1] == (Xyd, Zyd)
        (code=code, L=Ldict, bases=Dict(P => basis for P in "XYZ"), basis=basis,
         weights=Dict(P => length(Ldict[P]) for P in "XYZ"), Ypair=(Xyd, Zyd),
         min_weight_found=Dict('X' => mw["X"]["weight"], 'Z' => mw["Z"]["weight"]))
    end
end

function _is_logical(code::BBCode, p::PauliDict)
    alg = PauliAlgebra(code.n)
    v = from_dict(alg, p)
    return !any(anticommutes(alg, v, from_dict(alg, s)) for s in check_dicts(code)) &&
           addrow!(RowSpace(code.S), dict_to_vec(p, code.n))
end

const READOUT = Dict('X' => 'Z', 'Z' => 'X', 'Y' => 'Z')  # conjugate-basis termination

# The Y graph needs Cheeger constant >= 1 for circuit distance 6.
const DEFAULTS = Dict('X' => (;), 'Z' => (;), 'Y' => (target_cheeger=1.0, max_degree=4, order_seed=0))
# Verified 6-edge graphs: an option, not the default.
const SMALL_GRAPHS = Dict(
    'X' => (edges=((36, 37), (36, 38), (36, 39), (37, 38), (37, 55), (38, 59)),),
    'Z' => (edges=((36, 42), (36, 48), (36, 66), (42, 48), (42, 57), (45, 48)),),
)

"""
    bb_gauging(P, rounds; state_kind="+1", name="bb72", termination="conjugate", extra_targets=false,
               order_seed=nothing, separate_gauss=false, use_defaults=true, flux_keep=nothing,
               state_labels=nothing, graph_kw...)

Experiment measuring Pbar_0 (P ∈ 'XYZ') by gauging for `rounds` rounds, then destructive readout.
`state_kind` is "+1"/"-1" (eigenstate of Pbar_0) or "random"; other logical qubits are Zbar
eigenstates. `state_labels` gives an arbitrary logical state `[(12-letter label, sign bit)]`.
`extra_targets` also observes every other logical the final readout reads (distance checks).
"""
function bb_gauging(P::Char, rounds::Integer; state_kind="+1", name="bb72", termination="conjugate",
                    extra_targets=false, order_seed=nothing, separate_gauss=false, use_defaults=true,
                    flux_keep=nothing, state_labels=nothing, graph_kw...)
    gkw = Dict{Symbol,Any}(pairs(graph_kw))
    if use_defaults
        kw = Dict{Symbol,Any}(pairs(DEFAULTS[P]))
        merge!(kw, gkw)
        if order_seed === nothing
            order_seed = pop!(kw, :order_seed, nothing)
        else
            pop!(kw, :order_seed, nothing)
        end
        gkw = kw
    end
    s = bb_setup(name)
    code, L, basis = s.code, s.L[P], s.bases[P]
    n = code.n
    X0, Z0 = basis[1]
    if state_labels !== nothing
        gens = physical_state(state_labels, basis, n)
    elseif state_kind in ("+1", "-1")
        gens = StateGens([(L, state_kind == "+1" ? 0 : 1)])
    elseif state_kind == "random"
        gens = StateGens([(P in "XY" ? Z0 : X0, 0)])
    else
        throw(ArgumentError(state_kind))
    end
    if state_labels === nothing
        for (_, Zi) in basis[2:end]
            push!(gens, (Zi, 0))
        end
    end
    rb = termination == "conjugate" ? READOUT[P] : P
    checks = ordered_checks(code)
    G, prs = build_graph([c.pauli for c in checks], L; gkw...)
    extra = PauliDict[]
    if extra_targets
        for (Xi, Zi) in basis[2:end], op in (Xi, Zi)
            all(==(rb), values(op)) && push!(extra, deform_operator(op, L, G, n))
        end
    end
    return gauging_experiment(n, checks, L, rounds, rb; state=gens, logical_basis=basis, targets_extra=extra,
                              G, pairings=prs, order_seed, separate_gauss, flux_keep, name="$(name)_gauge_$(P)_r$(rounds)")
end
