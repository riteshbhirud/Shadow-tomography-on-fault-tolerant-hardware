# BB72 local Pauli settings on logical qubits 0..m-1, m <= 3.
# Qubit 0 of the local basis is the (Xbar_0, Zbar_0) pair of bb_setup; qubits 1 and 2 are
# weight-6 pairs disjoint from each other and from qubit 0, so any mix of Xbar_i, Zbar_i, Ybar_i
# can be gauged at once.

# (Xbar_i, Zbar_i) supports for i = 1, 2 (found by a seeded search and fixed here).
const EXTRA_PAIRS = [
    ((2, 17, 20, 26, 29, 32), (1, 2, 3, 4, 21, 23)),
    ((10, 13, 19, 22, 25, 31), (12, 14, 30, 31, 34, 35)),
]
const M_LOCAL = 3

const _LOCAL_BASIS = Ref{Any}(nothing)

"""
    local_basis()

The 12 local-basis pairs (Xbar_i, Zbar_i); pairs 3..11 complete a CSS symplectic basis.
"""
function local_basis()
    _LOCAL_BASIS[] === nothing || return _LOCAL_BASIS[]::LogicalBasis
    s = bb_setup()
    code = s.code
    n = code.n
    X0, Z0 = s.Ypair
    first = vcat([(X0, Z0)], [(PauliDict(q => 'X' for q in x), PauliDict(q => 'Z' for q in z)) for (x, z) in EXTRA_PAIRS])
    fx = [dict_to_vec(x, n) for (x, _) in first]
    fz = [dict_to_vec(z, n) for (_, z) in first]
    m = length(first)
    for i in 1:m, j in 1:m
        @assert symp_inner(fx[i], fz[j]) == (i == j)
    end
    NZ, NX = nullspace_gf2(code.HZ), nullspace_gf2(code.HX)
    KX = [vcat(NZ[i, :], falses(n)) for i in axes(NZ, 1)]
    KZ = [vcat(falses(n), NX[i, :]) for i in axes(NX, 1)]
    function extend(firsts, partners, pool)
        out = copy(firsts)
        cur = RowSpace(code.S)
        foreach(f -> addrow!(cur, f), firsts)
        for v in pool
            for (f, pz) in zip(firsts, partners)
                symp_inner(v, pz) && (v = v .⊻ f)
            end
            addrow!(cur, v) && push!(out, Vector{Bool}(v))
        end
        return out
    end
    Xs, Zs = extend(fx, fz, KX), extend(fz, fx, KZ)
    @assert length(Xs) == length(Zs) == code.k
    pairs = _pair_rest(collect(zip(fx, fz)), Xs[m+1:end], Zs[m+1:end])
    basis = LogicalBasis([(vec_to_dict(a), vec_to_dict(b)) for (a, b) in pairs])
    @assert basis[1] == (X0, Z0)
    return _LOCAL_BASIS[] = basis
end

"""
    logical_op(i, P)

Physical Pauli of Pbar_i (0-based i) in the local basis; Ybar_i has no sign (overlap 1).
"""
function logical_op(i::Integer, P::Char)
    basis = local_basis()
    d, sb = logical_to_physical(join(j == i ? P : 'I' for j in 0:length(basis)-1), basis, bb_setup().code.n)
    @assert sb == 0
    return d
end

const _GAUGE_SPECS = Dict{Tuple{Int,Char},Any}()

"""
    gauge_spec(i, P)

`(L, G, pairings, order_seed)` of the default gauging graph for Pbar_i.
"""
function gauge_spec(i::Integer, P::Char)
    get!(_GAUGE_SPECS, (i, P)) do
        L = logical_op(i, P)
        kw = Dict{Symbol,Any}(pairs(DEFAULTS[P]))
        order_seed = pop!(kw, :order_seed, nothing)
        G, prs = build_graph([c.pauli for c in ordered_checks(bb_setup().code)], L; kw...)
        (L, G, prs, order_seed)
    end
end

"""Vertex, edge and flux-check qubits added by gauging Pbar_i."""
function extra_qubits_of(i::Integer, P::Char)
    _, G, _, _ = gauge_spec(i, P)
    return number_of_nodes(G) + number_of_edges(G) + length(cycle_basis_edges(G))
end

"""
    compile_setting(bases)

Readout basis R ∈ {X, Z} and gauged operators `[(i, P)]` for `bases`: R minimises the
extra qubits, ties go to Z. Returns `(readout, gauged, extra_qubits_estimate)`.
"""
function compile_setting(bases::AbstractString)
    length(bases) <= M_LOCAL || throw(ArgumentError("local settings are supported on logical qubits 0..$(M_LOCAL - 1)"))
    best = nothing
    for R in "ZX"
        g = [(i - 1, b) for (i, b) in enumerate(bases) if b != R]
        cost = sum((extra_qubits_of(i, b) for (i, b) in g); init=0)
        (best === nothing || cost < best[1]) && (best = (cost, R, g))
    end
    return (readout=best[2], gauged=best[3], extra_qubits_estimate=best[1])
end

"""
    multi_gauging(specs, rounds, readout_basis; state=(), logical_basis=nothing, extra_ops=(), name="bb72_local")

Simultaneous gauging of logicals with disjoint supports, `specs = [(L, G, pairings, order_seed)]`.
For one spec the circuit equals `gauging_experiment`'s. Targets are the gauged logicals, then
`extra_ops` deformed through every graph.
"""
function multi_gauging(specs, rounds::Integer, readout_basis::Char; state=(), logical_basis=nothing, extra_ops=(),
                       name="bb72_local")
    code = bb_setup().code
    n_code = code.n
    code_checks = ordered_checks(code)
    supp = [Set(keys(s[1])) for s in specs]
    for a in eachindex(supp), b in a+1:length(supp)
        isempty(supp[a] ∩ supp[b]) || throw(ArgumentError("gauged logicals must have disjoint supports"))
    end
    all_edges = Tuple{Int,Tuple{Int,Int}}[]
    for (gi, (_, G, _, _)) in enumerate(specs)
        append!(all_edges, [(gi, e) for e in sorted_edges(G)])
    end
    eidx = Dict(ge => n_code + i - 1 for (i, ge) in enumerate(all_edges))
    n = n_code + length(all_edges)
    checks = Check[]
    for (gi, (L, G, _, order_seed)) in enumerate(specs)
        edges = sorted_edges(G)
        for v in sort(collect(keys(L)))
            p = PauliDict(v => L[v])
            for e in edges
                v in e && (p[eidx[(gi, e)]] = 'X')
            end
            od = order_seed === nothing ? nothing : _gauss_order(order_seed, v, p)
            push!(checks, Check(p; order=od, name="A$(gi-1)_$v"))
        end
    end
    for (gi, (_, G, _, _)) in enumerate(specs), c in cycle_basis_edges(G)
        push!(checks, Check(PauliDict(eidx[(gi, e)] => 'Z' for e in c); name="B$(gi-1)"))
    end
    for (ci, s) in enumerate(code_checks)
        p = copy(s.pauli)
        deformed = false
        for (gi, (_, _, prs, _)) in enumerate(specs), path in prs[ci]
            deformed = true
            for e in _path_edges(path)
                _toggle_z!(p, eidx[(gi, e)])
            end
        end
        od = s.order === nothing ? nothing : vcat(s.order, sort([q for q in keys(p) if q >= n_code]))
        push!(checks, Check(p; order=od, name=deformed ? "deformed" : "orig"))
    end
    edge_qubits = [eidx[ge] for ge in all_edges]
    readout = OrderedDict{Int,Char}(q => readout_basis for q in 0:n_code-1)
    foreach(q -> readout[q] = 'Z', edge_qubits)
    geidx = [(L, G, Dict(e => eidx[(gi, e)] for e in sorted_edges(G))) for (gi, (L, G, _, _)) in enumerate(specs)]
    extra = [_deform(PauliDict(op), geidx) for op in extra_ops]
    exp = Experiment(; n_data=n, active0=0:n_code-1, prep_checks=[c.pauli for c in code_checks], logical_state=collect(state),
                     epochs=[Epoch(checks, rounds; new_qubits=OrderedDict(q => 'Z' for q in edge_qubits), name="gauged")],
                     readout, targets=vcat([copy(s[1]) for s in specs], extra), name, logical_basis)
    exp.graph_info = OrderedDict{String,Any}(
        "V" => sum(length(s[1]) for s in specs), "E" => length(all_edges),
        "flux" => count(c -> startswith(c.name, "B"), checks), "cheeger" => [cheeger(s[2]) for s in specs],
        "max_check_weight" => maximum(length(c.pauli) for c in checks))
    return exp
end

"""
    local_experiment(bases, state_labels, rounds; spectators=false)

Experiment for the local setting `bases` on logical qubits 0..m-1. `state_labels` are logical
generators on all 12 qubits or on the m used ones (the rest is then maximally mixed). With
`spectators`, Rbar_j of every unused qubit is observed too. Returns `(experiment, labels, compile_setting result)`.
"""
function local_experiment(bases::AbstractString, state_labels, rounds::Integer; spectators::Bool=false)
    basis = local_basis()
    code = bb_setup().code
    n, k, m = code.n, length(basis), length(bases)
    st = [(lab * "I"^(k - length(lab)), s) for (lab, s) in state_labels]
    gens = physical_state(st, basis, n)
    comp = compile_setting(bases)
    R, gauged = comp.readout, comp.gauged
    gi = Set(i for (i, _) in gauged)
    readout_targets = [i for i in 0:m-1 if !(i in gi)]
    labels_g = [join(j == i ? P : 'I' for j in 0:m-1) for (i, P) in gauged]
    labels_r = [join(j == i ? R : 'I' for j in 0:m-1) for i in readout_targets]
    ops_r = [logical_op(i, R) for i in readout_targets]
    spectators && append!(ops_r, [logical_op(j, R) for j in m:k-1])
    if isempty(gauged)
        exp = Experiment(; n_data=n, active0=0:n-1, prep_checks=check_dicts(code), logical_state=gens, epochs=Epoch[],
                         readout=OrderedDict(q => R for q in 0:n-1), targets=ops_r, name="bb72_local_$bases",
                         logical_basis=basis)
        exp.graph_info = OrderedDict{String,Any}("V" => 0, "E" => 0, "flux" => 0)
    else
        specs = [gauge_spec(i, P) for (i, P) in gauged]
        exp = multi_gauging(specs, rounds, R; state=gens, logical_basis=basis, extra_ops=ops_r,
                            name="bb72_local_$(bases)_r$(rounds)")
    end
    # outcome i belongs to qubit i; spectator targets follow
    order_q = vcat([i for (i, _) in gauged], readout_targets)
    labels_all = vcat(labels_g, labels_r)
    perm = sort(1:m; by=j -> order_q[j])
    exp.targets = vcat(exp.targets[perm], exp.targets[m+1:end])
    exp.n_reported = m
    return exp, labels_all[perm], comp
end
