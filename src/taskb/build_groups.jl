# Build and verify the ensembles of BB72 and the gross code, with their group orders.
# Orders and the generator reduction need `using Oscar`.

const PUBLISHED = OrderedDict(
    # Sayginel et al. 2409.18175 Table VIII; Cross et al. 2407.18393 Sec 4.3; Albert 2608.05688 Eq. (39)
    "BB72" => OrderedDict{String,Any}("aut_physical" => 864, "aut_logical" => 864, "ibm_logical" => 72),
    "gross" => OrderedDict{String,Any}("aut_physical" => 288, "aut_logical" => 144, "ibm_logical" => 72,
                                "shifts_logical" => 36, "albert_two_fold_lower_bound" => 460800),
)

fixed_points(label::AbstractString) = parse(Int, rstrip(split(label, "fixed points: ")[2], ')'))

"""
    qc_crosscheck(code, e)

Recompute the logical and physical action of `e` from its gate list with QuantumClifford tableaux and compare
with `e.logical` and `e.phys`.
"""
function qc_crosscheck(code::StabilizerCode, e::Element)
    L = logicals(code)
    st = apply_gates!(_rows_tableau(L, code.n), e.gates)
    imgs = GF2Mat(reduce(vcat, [permutedims([xbit(st[i]); zbit(st[i])]) for i in axes(L, 1)]))
    return logical_coords(code, imgs) == e.logical && gates_to_symplectic(e.gates, code.n) == e.phys
end

"""
    build_groups(name; timings=false) -> (record, kept, ensembles, folds, code)

The `groups_<code>.json` record, the kept generators of every ensemble and the candidate elements. Needs
`using Oscar`.
"""
function build_groups(name::AbstractString; timings::Bool=false)
    code = make_bb_code(name, "fgh")
    k, n = logical_k(code), code.n
    ens, folds, info = build_ensembles(code)
    res = OrderedDict{String,Any}("code" => name, "n" => n, "k" => k, "basis" => code.meta[:basis],
                           "automorphism_search" => info, "published" => PUBLISHED[name],
                           "ensembles" => OrderedDict{String,Any}())
    hx, hz = code.meta[:hx], code.meta[:hz]
    wx = min_weight_logical_upto(hz, code.lz[:, n+1:end], 6)
    if name == "BB72"
        wz = min_weight_logical_upto(hx, code.lx[:, 1:n], 6)
        res["distance_check"] = OrderedDict{String,Any}("method" => "exhaustive meet-in-the-middle up to weight 6",
            "min_x_logical_weight" => wx[1], "min_z_logical_weight" => wz[1])
    else
        res["distance_check"] = OrderedDict{String,Any}(
            "method" => "exhaustive up to weight 6 (no logical found => d >= 7); upper bound from basis representatives",
            "min_x_logical_weight_upto6" => wx[1],
            "min_basis_rep_weight" => min(minimum(sum(code.lx; dims=2)), minimum(sum(code.lz; dims=2))))
    end
    res["basis_weights"] = OrderedDict("X" => vec(sum(code.lx; dims=2)), "Z" => vec(sum(code.lz; dims=2)))
    haskey(code.meta, :fgh_info) && (res["fgh_info"] = code.meta[:fgh_info])
    kept = OrderedDict{String,Vector{GF2Mat}}()
    sp = symplectic_group_order(k)
    for ename in ENSEMBLES
        t = time()
        els = ens[ename]
        if isempty(els)
            res["ensembles"][ename] = OrderedDict{String,Any}("description" => ENSEMBLE_DESCRIPTIONS[ename], "order" => 1,
                                                       "generators" => Any[])
            kept[ename] = GF2Mat[]
            continue
        end
        km, order = reduce_generators([e.logical for e in els], 2k)
        order2 = group_order(km, 2k)              # an independent recognition run
        bfs = closure(km; cap=1_000_000)
        kept_el = [els[findfirst(e -> e.logical == m, els)] for m in km]
        res["ensembles"][ename] = OrderedDict{String,Any}(
            "description" => ENSEMBLE_DESCRIPTIONS[ename],
            "n_candidate_elements" => length(els),
            "n_distinct_logical_actions" => length(unique(e.logical for e in els)),
            "order" => order, "order_second_seed" => order2,
            "order_bfs_exact" => bfs[3] ? bfs[1] : nothing,
            "frac_of_Sp" => Float64(order // sp),
            "generators" => [OrderedDict{String,Any}("label" => e.label, "logical_support" => logical_support(e.logical, k),
                                              "n_physical_gates" => length(e.gates),
                                              "stim_crosscheck" => qc_crosscheck(code, e)) for e in kept_el],
        )
        timings && (res["ensembles"][ename]["seconds"] = py_round(time() - t, 1))
        kept[ename] = km
    end
    alls = all_elements(ens)
    res["all_elements_verified"] = OrderedDict{String,Any}("count" => length(alls),
        "checks" => "physical op preserves stabilizer rowspace; logical images in normalizer " *
                    "and residual in stabilizer group; logical action symplectic")
    res["fold_stim_crosscheck_all"] = all(e -> qc_crosscheck(code, e), folds)
    res["pauli_corrections"] = OrderedDict{String,Any}(
        "n_elements_needing_correction" => count(e -> !isempty(e.pauli_correction), alls),
        "n_elements" => length(alls),
        "n_folds_needing_correction" => count(e -> !isempty(e.pauli_correction), folds),
        "n_folds" => length(folds),
        "note" => "each Element.gates = bare CZ/S/H/SWAP circuit + Pauli correction (X/Y/Z gates); the corrected " *
                  "circuit maps every stabilizer to +stabilizer (exact stim check in ensembles.make_element). " *
                  "Corrections do not change the symplectic/logical action, hence no lambda.")
    return res, kept, ens, folds, code
end

"""
    crosscheck_sample(ens, folds, kept)

Elements for the state-level sign check: those realising the kept generators of the three largest ensembles,
then the first fold with fixed points.
"""
function crosscheck_sample(ens::AbstractDict, folds, kept::AbstractDict)
    alls = all_elements(ens)
    sample = Element[]
    add!(e) = any(s -> s.label == e.label, sample) || push!(sample, e)
    for ename in ("aut", "albert_two_fold", "aut_plus_all_folds"), km in kept[ename]
        add!(alls[findfirst(e -> e.logical == km, alls)])
    end
    i = findfirst(e -> fixed_points(e.label) > 0, folds)
    i === nothing || add!(folds[i])
    return sample
end

"""
    signed_crosscheck(code, e, taus, rng; trials=3) -> (prediction_ok, negative_control_ok, sign_negative_control_ok)

State-level sign check of one element. Random logical stabilizer states are evolved through `e.gates`. Every
code stabilizer and every predicted image (signed with `taus`) must have expectation +1. Negative controls:
the identity prediction fails when Λ ≠ I, and the bare circuit flips some stabilizer sign when the correction
is nonempty.
"""
function signed_crosscheck(code::StabilizerCode, e::Element, taus, rng; trials::Integer=3)
    n, k = code.n, logical_k(code)
    L = logicals(code)
    Lp = [rowpauli(L[i, :], n) for i in 1:2k]
    Q = gf2mul(e.logical, L)
    Qp = [rowpauli(Q[i, :], n) for i in 1:2k]
    stabs = [rowpauli(code.stabs[i, :], n) for i in axes(code.stabs, 1)]
    ok_pred, ok_neg = true, e.logical == I
    ok_sign_neg = isempty(e.pauli_correction)
    for _ in 1:trials
        rs = random_stabilizer(rng, k)
        coords = zeros(Bool, k, 2k)
        for j in 1:k, q in 1:k
            coords[j, q], coords[j, k+q] = rs[j][q]
        end
        phys = gf2mul(coords, L)
        rows = [stabs; [rowpauli(phys[j, :], n) for j in 1:k]]
        md = MixedDestabilizer(Stabilizer(rows))
        rank(md) == n || error("state not pure")
        st0 = copy(stabilizerview(md))
        st = apply_gates!(copy(st0), e.gates)
        all(s -> expect(s, st) == 1, stabs) || (ok_pred = false)
        for j in 1:k
            prodL = zero(PauliOperator, n)
            pred = zero(PauliOperator, n)
            for i in 1:2k
                coords[j, i] || continue
                prodL = prodL * Lp[i]
                pred = pred * Qp[i]
                taus[i] == -1 && (pred.phase[] = (pred.phase[] + 0x2) & 0x3)
            end
            pred.phase[] = (pred.phase[] - prodL.phase[]) & 0x3   # P_j (sign +) = i^{-p} prodL
            if isodd(pred.phase[]) || expect(pred, st) != 1
                ok_pred = false
            end
            e.logical == I || abs(expect(rowpauli(phys[j, :], n), st)) == 1 || (ok_neg = true)
        end
        if !isempty(e.pauli_correction)
            st2 = apply_gates!(copy(st0), e.gates_uncorrected)
            any(s -> expect(s, st2) == -1, stabs) && (ok_sign_neg = true)
        end
    end
    return ok_pred, ok_neg, ok_sign_neg
end
