# Cheap logical Clifford generators of BB codes: (physical gate list, physical symplectic, logical symplectic).
#
# shifts : x^a y^b on both blocks
# ibm    : shifts + the Bravyi ZX-duality L(α) <-> R(αᵀ), then Hⁿ
# aut    : every affine permutation automorphism and ZX-duality
# folds  : Z- and X-diagonal matching circuits of every involutive ZX-duality that preserves the code

"""
    Element

A verified logical Clifford. `gates` is the physical gate list including the Pauli correction that makes it map
every stabilizer to +stabilizer; `gates_uncorrected` is the bare circuit. The correction does not change `phys`
or `logical`.
"""
struct Element
    label::String
    gates::Vector{Gate}
    phys::GF2Mat
    logical::GF2Mat
    pauli_correction::Vector{Gate}
    gates_uncorrected::Vector{Gate}
end

"""
    shift_perm(l, m, a, b)

Qubit permutation of the shift x^a y^b on both blocks.
"""
shift_perm(l, m, a, b) = affine_perm(l, m, ((1, 0), (0, 1)), (a, b), (a, b), false)

"""
    bravyi_zx_perm(l, m)

The permutation q(L, α) <-> q(R, αᵀ).
"""
function bravyi_zx_perm(l, m)
    N = l * m
    perm = zeros(Int, 2N)
    for a in 0:l-1, b in 0:m-1
        i = a * m + b
        it = mod(-a, l) * m + mod(-b, m)
        perm[i+1] = N + it
        perm[N+i+1] = it
    end
    return perm
end

"""
    make_element(code, label, gates, phys)

Verified [`Element`](@ref) with the sign-fixing Pauli correction appended.
"""
function make_element(code::StabilizerCode, label, gates, phys)
    lam = logical_action(code, phys)
    corr, _ = pauli_correction(code.stabs, gates, code.n)
    full = [gates; corr]
    _, flip = stabilizer_sign_flips(code.stabs, full, code.n)
    any(flip) && error("$label: Pauli correction failed to fix stabilizer signs")
    return Element(label, full, phys, lam, corr, collect(gates))
end

perm_element(code, label, perm) = make_element(code, label, gates_perm(perm), perm_op(perm))
zx_element(code, label, perm) =
    make_element(code, label, [gates_perm(perm); gates_h_all(code.n)], gf2mul(perm_op(perm), h_all(code.n)))

"""
    logical_signs(code, e)

τᵢ ∈ {±1} for each logical basis row Lᵢ: U Lᵢ U† = τᵢ Qᵢ sᵢ with Qᵢ the Hermitian Pauli with bits
(row i of Λ)·L and sᵢ a +element of S. Computed with QuantumClifford.
"""
function logical_signs(code::StabilizerCode, e::Element)
    n = code.n
    L = logicals(code)
    st = apply_gates!(_rows_tableau(L, n), e.gates)
    q = gf2mul(e.logical, L)
    taus = Int[]
    for i in axes(L, 1)
        r = st[i] * rowpauli(q[i, :], n)
        b = [xbit(r); zbit(r)]
        in_rowspace(code.stabs, b) || error("residual of logical image is not a stabilizer")
        push!(taus, _sign(r) * css_group_sign(b, n))
    end
    return taus
end

"""
    bb_families(code) -> (fam, auts, duals, folds, info)

Named families, automorphisms, ZX-dualities and valid fold circuits of a BB code.
"""
function bb_families(code::StabilizerCode)
    l, m = code.meta[:l], code.meta[:m]
    n = code.n
    N = l * m
    hx, hz = code.meta[:hx], code.meta[:hz]
    fam = OrderedDict{String,Element}()
    fam["shift_x"] = perm_element(code, "shift x", shift_perm(l, m, 1, 0))
    fam["shift_y"] = perm_element(code, "shift y", shift_perm(l, m, 0, 1))
    fam["zx_bravyi"] = zx_element(code, "ZX duality L(a)<->R(a^T), H^n", bravyi_zx_perm(l, m))
    auts, duals = search_bb_automorphisms(hx, hz, l, m)
    aut_elems = [perm_element(code, "aut#$(i-1)", p) for (i, p) in enumerate(auts)]
    dual_elems = [zx_element(code, "zxdual#$(i-1)", p) for (i, p) in enumerate(duals)]
    shifts = [shift_perm(l, m, a, b) for a in 0:l-1 for b in 0:m-1]
    all_dual_perms = [s[p.+1] for p in duals for s in shifts]   # p then s: qubit i -> s[p[i]]
    invol = [p for p in all_dual_perms if p[p.+1] == 0:n-1]
    fold_elems = Element[]
    fold_tried = 0
    for (i, p) in enumerate(invol)
        smat = matching_matrix(p)
        for (kind, phys, gates) in (("Zdiag", z_diag(smat), gates_z_diag(smat)),
                                    ("Xdiag", x_diag(smat), gates_x_diag(smat)))
            fold_tried += 1
            if preserves_stabilizer_group(code, phys)
                label = "$kind matching of involutive ZX-duality #$(i-1) (fixed points: $(count(p .== 0:n-1)))"
                push!(fold_elems, make_element(code, label, gates, phys))
            end
        end
    end
    info = OrderedDict{String,Any}(
        "n_perm_auts" => length(auts) * N,
        "n_zx_dualities" => length(duals) * N,
        "n_perm_auts_c0" => length(auts),
        "n_dual_c0" => length(duals),
        "n_involutive_dualities" => length(invol),
        "n_fold_circuits_tried" => fold_tried,
        "n_fold_circuits_valid" => length(fold_elems),
    )
    return fam, aut_elems, dual_elems, fold_elems, info
end

"""
    dedupe(elems)

Drop elements whose logical action repeats an earlier one.
"""
function dedupe(elems)
    seen = Set{GF2Mat}()
    return [e for e in elems if !(e.logical in seen) && (push!(seen, e.logical); true)]
end

"""
    ENSEMBLES

Names of the native ensembles that Task B builds groups and shadow maps for.
"""
const ENSEMBLES = ["readout_only", "shifts", "ibm", "aut", "albert_two_fold", "aut_plus_all_folds"]

const ENSEMBLE_DESCRIPTIONS = OrderedDict(
    "readout_only" => "no unitary; destructive all-Z / all-X readout only",
    "shifts" => "shift automorphisms x^a y^b (logical CNOT circuits)",
    "ibm" => "shifts + Bravyi ZX-duality L(a)<->R(a^T) then H^n  (Sayginel 'Aut_IBM')",
    "aut" => "all SWAP/H automorphisms: affine permutation automorphisms + ZX-dualities (Sayginel 'Aut')",
    "albert_two_fold" => "aut + Z- and X-diagonal matching circuits (CZ perfect matchings) of all translates of the " *
                         "block-swap duality sigma (Albert 2608.05688 Sec VII A generator set)",
    "aut_plus_all_folds" => "aut + Z-/X-diagonal matching circuits of EVERY involutive ZX-duality whose circuit " *
                            "preserves the stabilizer group (CZ on 2-cycles, sqrt(Z)/sqrt(X) on fixed points)",
)

"""
    build_ensembles(code) -> (ensembles, folds, info)

Candidate generator elements of every ensemble. The order matters: it decides which generators the reduction
keeps.
"""
function build_ensembles(code::StabilizerCode)
    l, m = code.meta[:l], code.meta[:m]
    fam, auts, duals, folds, info = bb_families(code)
    shifts = [fam["shift_x"], fam["shift_y"]]
    sig = bravyi_zx_perm(l, m)
    albert_diag = Element[]
    for a in 0:l-1, b in 0:m-1
        s = matching_matrix(shift_perm(l, m, a, b)[sig.+1])
        push!(albert_diag, make_element(code, "Zdiag matching of sigma-translate x^$(a)y^$(b)", gates_z_diag(s), z_diag(s)))
        push!(albert_diag, make_element(code, "Xdiag matching of sigma-translate x^$(a)y^$(b)", gates_x_diag(s), x_diag(s)))
    end
    ens = OrderedDict{String,Vector{Element}}(
        "readout_only" => Element[],
        "shifts" => shifts,
        "ibm" => [shifts; fam["zx_bravyi"]],
        "aut" => [shifts; auts; duals],
        "albert_two_fold" => [albert_diag; shifts; fam["zx_bravyi"]; auts; duals],
        "aut_plus_all_folds" => [shifts; auts; duals; folds],
    )
    return ens, folds, info
end

"""
    all_elements(ensembles)

The distinct `Element` objects over all ensembles, in first-seen order.
"""
function all_elements(ens::AbstractDict)
    seen = IdDict{Element,Nothing}()
    out = Element[]
    for name in ENSEMBLES, e in ens[name]
        haskey(seen, e) || (seen[e] = nothing; push!(out, e))
    end
    return out
end

"""
    logical_support(lam, k)

0-based logical qubits on which `lam` acts non-trivially.
"""
function logical_support(lam::AbstractMatrix, k::Integer)
    d = lam .⊻ I(2k)
    return [q for q in 0:k-1 if any(d[[q + 1, k + q + 1], :]) || any(d[:, [q + 1, k + q + 1]])]
end
