# Epoch circuit builder. An experiment is ideal preparation,
# epochs of checks measured for r noisy rounds, then destructive readout. Detectors and
# observables are derived by the tracker, never declared.

const PAULI_GATE = Dict('X' => "CX", 'Y' => "CY", 'Z' => "CZ")

"""
    Check(pauli; layer_of=nothing, group=nothing, name="", order=nothing)

One stabilizer check with optional schedule hints: explicit gate layers per qubit, a group
index, or a gate order over the check's qubits (hook control).
"""
mutable struct Check
    pauli::PauliDict
    layer_of::Union{Nothing,OrderedDict{Int,Int}}
    group::Union{Nothing,Int}
    name::String
    order::Union{Nothing,Vector{Int}}
end
Check(pauli; layer_of=nothing, group=nothing, name="", order=nothing) =
    Check(PauliDict(pauli), layer_of, group, name, order === nothing ? nothing : collect(Int, order))

"""
    Epoch(checks, rounds; new_qubits=PauliDict(), name="")

Checks measured for `rounds` noisy rounds after resetting `new_qubits` (qubit => 'Z' for |0>, 'X' for |+>).
"""
struct Epoch
    checks::Vector{Check}
    rounds::Int
    new_qubits::OrderedDict{Int,Char}
    name::String
end
Epoch(checks, rounds; new_qubits=OrderedDict{Int,Char}(), name="") =
    Epoch(collect(Check, checks), rounds, OrderedDict{Int,Char}(new_qubits), name)

const LogicalBasis = Vector{Tuple{PauliDict,PauliDict}}
const StateGens = Vector{Tuple{PauliDict,Int}}

"""
    Experiment

A code-deformation experiment (builder.Experiment): `n_data` qubits, `active0` initially active,
ideal `prep_checks` and `logical_state` (Pauli dict, sign bit) generators, `epochs`, the final
`readout` basis per qubit and the logical `targets`. Unspecified logical qubits of
`logical_basis` are prepared maximally mixed.
"""
mutable struct Experiment
    n_data::Int
    active0::Vector{Int}
    prep_checks::Vector{PauliDict}
    logical_state::StateGens
    destabilizers::Union{Nothing,Vector{PauliDict}}
    epochs::Vector{Epoch}
    readout::OrderedDict{Int,Char}
    targets::Vector{PauliDict}
    coords::Union{Nothing,OrderedDict{Int,Vector{Float64}}}
    name::String
    logical_basis::Union{Nothing,LogicalBasis}
    tracker_readout_fallback::Bool
    graph::Any
    graph_info::Union{Nothing,OrderedDict{String,Any}}
    n_reported::Union{Nothing,Int}
end

function Experiment(; n_data, active0, prep_checks, logical_state, destabilizers=nothing, epochs, readout, targets,
                    coords=nothing, name="", logical_basis=nothing)
    Experiment(n_data, collect(Int, active0), PauliDict[PauliDict(p) for p in prep_checks],
               StateGens([(PauliDict(p), Int(s)) for (p, s) in logical_state]),
               destabilizers === nothing ? nothing : PauliDict[PauliDict(d) for d in destabilizers],
               collect(Epoch, epochs), OrderedDict{Int,Char}(readout), PauliDict[PauliDict(t) for t in targets],
               coords, name, logical_basis === nothing ? nothing : LogicalBasis([(PauliDict(x), PauliDict(z)) for (x, z) in logical_basis]),
               false, nothing, nothing, nothing)
end

"""
    BuildResult

Circuit plus record bookkeeping and cost of one build: `n_meas` records, per target the absolute
record ids whose parity is its raw outcome, `detector_count`, `cost`, `epoch_info`, and `n_prep`
noiseless preparation records.
"""
struct BuildResult
    circuit::Circuit
    n_meas::Int
    target_records::Vector{Vector{Int}}
    detector_count::Int
    cost::OrderedDict{String,Int}
    epoch_info::Vector{OrderedDict{String,Any}}
    n_prep::Int
end

"""
    NOISE_MODELS

The circuit noise models by name: multipliers of p per channel (reset flip, measurement flip,
two-qubit depolarizing, idle depolarizing per layer, data depolarizing per round) and whether
ancillas get idle noise.
"""
const NOISE_MODELS = Dict(
    "sd6" => (reset=1, meas=1, gate2=1, idle=1, round_data=0, idle_ancilla=1),
    "stim_generated" => (reset=1, meas=1, gate2=1, idle=0, round_data=1, idle_ancilla=1),
    "bravyi_cal" => (reset=1, meas=1, gate2=1, idle=1, round_data=0, idle_ancilla=0),
)
"""
    DEFAULT_NOISE

The noise model used when none is given: `"bravyi_cal"`, calibrated against Bravyi et al.'s BB72 memory
error rates.
"""
const DEFAULT_NOISE = "bravyi_cal"

"""
    group_checks(checks)

Partition checks into groups whose members act identically on every shared qubit. Explicit
groups are honoured; otherwise greedy colouring over the natural, reversed, by-weight and 20
`random.Random(0)` orders, keeping the first with the fewest groups.
"""
function group_checks(checks::Vector{Check})
    if all(c -> c.group !== nothing, checks)
        groups = OrderedDict{Int,Vector{Int}}()
        for (i, c) in enumerate(checks)
            push!(get!(groups, c.group, Int[]), i - 1)
        end
        out = [groups[g] for g in sort!(collect(keys(groups)))]
    else
        N = length(checks)
        orders = [collect(0:N-1), collect(N-1:-1:0), sort(0:N-1; by=i -> -length(checks[i+1].pauli), alg=MergeSort)]
        rng = PyRandom(0)
        for _ in 1:20
            push!(orders, pyshuffle!(rng, collect(0:N-1)))
        end
        cands = [_greedy_groups(checks, o) for o in orders]
        out = cands[argmin(length.(cands))]
        out = [sort(g) for g in out]
    end
    for g in out, a in eachindex(g), b in a+1:length(g)
        pa, pb = checks[g[a]+1].pauli, checks[g[b]+1].pauli
        for (q, s) in pa
            haskey(pb, q) && pb[q] != s && throw(ArgumentError("group contains qubit-wise non-commuting checks"))
        end
    end
    return out
end

function _greedy_groups(checks, order)
    out = Vector{Int}[]
    for i in order
        c = checks[i+1]
        placed = false
        for g in out
            if all(get(checks[j+1].pauli, q, s) == s for j in g for (q, s) in c.pauli)
                push!(g, i)
                placed = true
                break
            end
        end
        placed || push!(out, [i])
    end
    return out
end

"""
    schedule_group(checks, idxs)

Gate layers `[(check index, data qubit), ...]` of one group (0-based check indices), no qubit
used twice per layer. The k-th gate of every check is placed before any (k+1)-th gate.
"""
function schedule_group(checks::Vector{Check}, idxs::Vector{Int})
    layers = Dict{Int,Vector{Tuple{Int,Int}}}()
    if all(i -> checks[i+1].layer_of !== nothing, idxs)
        for i in idxs, (q, L) in checks[i+1].layer_of
            push!(get!(layers, L, Tuple{Int,Int}[]), (i, q))
        end
        for (L, gates) in layers
            used = [q for (_, q) in gates]
            allunique(used) || throw(ArgumentError("explicit schedule conflict in layer $L"))
        end
        return [layers[L] for L in sort!(collect(keys(layers)))]
    end
    busy_anc = Dict{Int,Set{Int}}()
    busy_dat = Dict{Int,Set{Int}}()
    seqs = Dict(i => (checks[i+1].order !== nothing ? checks[i+1].order : collect(keys(checks[i+1].pauli))) for i in idxs)
    for i in idxs
        sort(seqs[i]) == sort(collect(keys(checks[i+1].pauli))) || throw(ArgumentError("order must list the check's qubits"))
    end
    prev = Dict(i => -1 for i in idxs)
    maxlen = maximum(length(seqs[i]) for i in idxs)
    for k in 1:maxlen, i in idxs
        k > length(seqs[i]) && continue
        q = seqs[i][k]
        L = checks[i+1].order !== nothing ? prev[i] + 1 : 0
        while i in get(busy_anc, L, ()) || q in get(busy_dat, L, ())
            L += 1
        end
        push!(get!(busy_anc, L, Set{Int}()), i)
        push!(get!(busy_dat, L, Set{Int}()), q)
        push!(get!(layers, L, Tuple{Int,Int}[]), (i, q))
        prev[i] = L
    end
    return [layers[L] for L in sort!(collect(keys(layers)))]
end

"""
    build(exp, p; eigen_targets=false, noise=DEFAULT_NOISE)

Circuit of `exp` at noise strength `p`. With `eigen_targets` the +1 eigenstates of the targets are
prepared instead of `exp.logical_state`, so observables are deterministic (for the DEM).
"""
function build(exp::Experiment, p::Real; eigen_targets::Bool=false, noise::AbstractString=DEFAULT_NOISE)
    p = Float64(p)
    n = exp.n_data
    alg = PauliAlgebra(n)
    tr = Tracker(n)
    c = Circuit()
    n_meas = 0
    det_count = 0
    if exp.coords !== nothing && !isempty(exp.coords)
        for (q, xy) in sort!(collect(exp.coords); by=first)
            append!(c, "QUBIT_COORDS", [q], xy)
        end
    end
    nm = NOISE_MODELS[noise]
    add_noise(name, qs, prob) = (prob > 0 && !isempty(qs) && append!(c, name, sort!(collect(qs)), [prob]); nothing)
    p_reset, p_meas, p_2q = nm.reset * p, nm.meas * p, nm.gate2 * p
    p_idle, p_round = nm.idle * p, nm.round_data * p

    emitted = Bits[]
    function emit_detector(rec_bits::Bits, coord)
        push!(emitted, rec_bits)
        append!(c, "DETECTOR", [rec(i - n_meas) for i in setbits(rec_bits)], coord)
        det_count += 1
    end

    # Pauli => record id of its latest measurement (nothing: a fresh reset, known +1), valid while
    # nothing anticommuting has been measured or reset since. Insertion order matters.
    last_meas = OrderedDict{PVec,Union{Int,Nothing}}()
    invalidate(op) = foreach(v -> delete!(last_meas, v), [v for v in keys(last_meas) if anticommutes(alg, v, op)])

    function local_combo(vec::PVec)
        cands = [(v, r) for (v, r) in last_meas if support_overlaps(alg, v, vec)]
        isempty(cands) && return nothing
        rows = Dict{Int,Tuple{PVec,Bits}}()
        for (v, r) in cands
            rb = r === nothing ? Bits() : bit(r)
            while !iszero(v)
                p_ = pivot(alg, v)
                if haskey(rows, p_)
                    v = xor(v, rows[p_][1])
                    rb = xor(rb, rows[p_][2])
                else
                    rows[p_] = (v, rb)
                    break
                end
            end
        end
        t, rb = vec, Bits()
        while !iszero(t)
            p_ = pivot(alg, t)
            haskey(rows, p_) || return nothing
            t = xor(t, rows[p_][1])
            rb = xor(rb, rows[p_][2])
        end
        return rb
    end

    # Prefer the same check's previous record, then a local combination, then the tracker's.
    function t_measure(vec::PVec, rid::Int, snapshot=nothing)
        loc = local_combo(vec)
        loc === nothing && (loc = snapshot)
        res = measure!(tr, vec, rid)
        absent = !haskey(last_meas, vec)
        prev = absent ? nothing : last_meas[vec]
        invalidate(vec)
        last_meas[vec] = rid
        (res === nothing || !iszero(res[2])) && return nothing
        !absent && prev !== nothing && return xor(bit(prev), bit(rid))
        loc !== nothing && return xor(loc, bit(rid))
        return xor(res[1], bit(rid))
    end

    function t_reset(q, b::Char)
        reset!(tr, q, b)
        invalidate(from_dict(alg, ((q, b),)))
        invalidate(from_dict(alg, ((q, b == 'Z' ? 'X' : 'Z'),)))
        last_meas[from_dict(alg, ((q, b),))] = nothing
    end

    active = Set(exp.active0)
    append!(c, "R", sort!(collect(active)))
    for chk in exp.prep_checks
        append!(c, "MPP", mpp_targets(chk))
        t_measure(from_dict(alg, chk), n_meas)
        n_meas += 1
    end
    if eigen_targets
        lstate = StateGens([(t, 0) for t in exp.targets])
        destabs = _destabilizers_for(alg, lstate, exp)
    else
        lstate, destabs = exp.logical_state, exp.destabilizers
        destabs === nothing && (destabs = isempty(lstate) ? PauliDict[] : _destabilizers_for(alg, lstate, exp))
    end
    _check_state_generators(alg, lstate)
    lrec = Int[]
    for (lp, _) in lstate
        append!(c, "MPP", mpp_targets(lp))
        push!(lrec, n_meas)
        n_meas += 1
    end
    for ((lp, sign), dstab, r) in zip(lstate, destabs, lrec)
        for (q, s) in sortedpauli(dstab)
            append!(c, PAULI_GATE[s], [rec(r - n_meas), qubit(q)])
            sign != 0 && append!(c, string(s), [q])
        end
    end
    # Unspecified logical qubits are maximally mixed: MPP of destabilizer d_j then of h_j
    # gives h_j a random sign per shot; the records are unused and the tracker never sees them.
    if exp.logical_basis !== nothing && !isempty(exp.logical_basis)
        gens = [g for (g, _) in lstate]
        extra = complete_generators(alg, gens, exp.logical_basis)
        if !isempty(extra)
            dmix = destabilizers_from_basis(alg, vcat(gens, extra), exp.logical_basis)[length(gens)+1:end]
            for (h, d) in zip(extra, dmix), op in (d, h)
                append!(c, "MPP", mpp_targets(op))
                n_meas += 1
            end
        end
    end
    append!(c, "TICK")
    n_prep = n_meas

    anc_base = n
    target_records = Union{Nothing,Bits}[nothing for _ in exp.targets]
    epoch_info = OrderedDict{String,Any}[]
    total_layers = 0
    max_anc = 0
    for (ep_idx, ep) in enumerate(exp.epochs)
        zq = [q for (q, b) in ep.new_qubits if b == 'Z']
        xq = [q for (q, b) in ep.new_qubits if b == 'X']
        if !isempty(zq)
            append!(c, "R", sort(zq))
            add_noise("X_ERROR", zq, p_reset)
        end
        if !isempty(xq)
            append!(c, "RX", sort(xq))
            add_noise("Z_ERROR", xq, p_reset)
        end
        isempty(zq) && isempty(xq) || add_noise("DEPOLARIZE1", active, p_idle)
        for (q, b) in ep.new_qubits
            t_reset(q, b)
        end
        union!(active, keys(ep.new_qubits))
        if !(isempty(zq) && isempty(xq))
            append!(c, "TICK")
            total_layers += 1
        end
        snap = [local_combo(from_dict(alg, ch.pauli)) for ch in ep.checks]
        groups = group_checks(ep.checks)
        sched = [schedule_group(ep.checks, g) for g in groups]
        anc(i) = anc_base + i
        max_anc = max(max_anc, length(ep.checks))
        for rnd in 0:ep.rounds-1
            add_noise("DEPOLARIZE1", active, p_round)
            for (g, layers) in zip(groups, sched)
                ancs = [anc(i) for i in g]
                append!(c, "RX", ancs)
                add_noise("Z_ERROR", ancs, p_reset)
                add_noise("DEPOLARIZE1", active, p_idle)
                append!(c, "TICK")
                total_layers += 1
                for layer in layers
                    used = Set{Int}()
                    for (i, q) in layer
                        append!(c, PAULI_GATE[ep.checks[i+1].pauli[q]], [anc(i), q])
                        push!(used, anc(i), q)
                    end
                    p_2q > 0 && append!(c, "DEPOLARIZE2", [x for (i, q) in layer for x in (anc(i), q)], [p_2q])
                    idle_q = nm.idle_ancilla != 0 ? union(active, ancs) : copy(active)
                    add_noise("DEPOLARIZE1", setdiff(idle_q, used), p_idle)
                    append!(c, "TICK")
                    total_layers += 1
                end
                add_noise("Z_ERROR", ancs, p_meas)
                append!(c, "MX", ancs)
                add_noise("DEPOLARIZE1", active, p_idle)
                last = Tuple{Int,Union{Nothing,Bits}}[]
                for (k, i) in enumerate(g)
                    rid = n_meas + k - 1
                    push!(last, (i, t_measure(from_dict(alg, ep.checks[i+1].pauli), rid, rnd == 0 ? snap[i+1] : nothing)))
                end
                n_meas += length(g)
                for (i, det) in last
                    det === nothing || emit_detector(det, [ep_idx - 1, rnd, i])
                end
                append!(c, "TICK")
                total_layers += 1
            end
            for (t_idx, t) in enumerate(exp.targets)
                if target_records[t_idx] === nothing
                    qres = query(tr, from_dict(alg, t))
                    qres !== nothing && iszero(qres[2]) && (target_records[t_idx] = qres[1])
                end
            end
        end
        push!(epoch_info, OrderedDict{String,Any}(
            "name" => ep.name, "n_checks" => length(ep.checks), "n_groups" => length(groups),
            "layers_per_round" => sum(2 + length(L) for L in sched), "rounds" => ep.rounds,
            "new_data" => length(ep.new_qubits)))
    end

    qs = sort!(collect(keys(exp.readout)))
    flipname = Dict('X' => "Z_ERROR", 'Z' => "X_ERROR", 'Y' => "X_ERROR")
    for q in qs
        add_noise(flipname[exp.readout[q]], [q], p_meas)
    end
    mname = Dict('X' => "MX", 'Y' => "MY", 'Z' => "M")
    for q in qs
        append!(c, mname[exp.readout[q]], [q])
    end
    # a known check fully read out in its own basis gives a local detector
    rid_of = Dict(q => n_meas + k - 1 for (k, q) in enumerate(qs))
    last_checks = isempty(exp.epochs) ? [Check(pc) for pc in exp.prep_checks] : exp.epochs[end].checks
    readout_dets = Bits[]
    for chk in last_checks
        if all(get(exp.readout, q, nothing) == s for (q, s) in chk.pauli)
            prev = get(last_meas, from_dict(alg, chk.pauli), nothing)
            prev === nothing && continue
            bits = bit(prev)
            for q in keys(chk.pauli)
                bits = xor(bits, bit(rid_of[q]))
            end
            push!(readout_dets, bits)
        end
    end
    n_tracker_det = 0
    tracker_dets = Bits[]
    for (k, q) in enumerate(qs)
        res = measure!(tr, from_dict(alg, ((q, exp.readout[q]),)), n_meas + k - 1)
        if res !== nothing && iszero(res[2])
            n_tracker_det += 1
            push!(tracker_dets, xor(res[1], bit(n_meas + k - 1)))
        end
    end
    n_meas += length(qs)
    # Projected onto the readout records, the local readout detectors must span the tracker's
    # deterministic readout parities.
    ro_mask = bits_of(rid_of[q] for q in qs)
    proj_rank = rank_bits([b & ro_mask for b in readout_dets])
    if proj_rank != n_tracker_det && exp.tracker_readout_fallback
        for bits in tracker_dets
            if rank_bits([b & ro_mask for b in vcat(readout_dets, [bits])]) > proj_rank
                push!(readout_dets, bits)
                proj_rank += 1
            end
        end
    end
    proj_rank == n_tracker_det ||
        error("readout detectors incomplete: projected rank $proj_rank vs tracker $n_tracker_det")
    for (k, bits) in enumerate(readout_dets)
        emit_detector(bits, [length(exp.epochs), 0, k - 1])
    end
    for (t_idx, t) in enumerate(exp.targets)
        if target_records[t_idx] === nothing
            qres = query(tr, from_dict(alg, t))
            qres !== nothing && iszero(qres[2]) && (target_records[t_idx] = qres[1])
        end
        target_records[t_idx] === nothing && error("target $(t_idx - 1) never became determined by records")
    end
    # XOR-ing a detector leaves an observable's noiseless value unchanged
    reduced = [_reduce_by_detectors(b, emitted) for b in target_records]
    tlist = Vector{Int}[]
    for (t_idx, bits) in enumerate(reduced)
        ids = setbits(bits)
        push!(tlist, ids)
        append!(c, "OBSERVABLE_INCLUDE", [rec(i - n_meas) for i in ids], [t_idx - 1])
    end
    cost = OrderedDict{String,Int}(
        "merged_rounds" => sum((ep.rounds for ep in exp.epochs); init=0),
        "data_qubits" => length(active),
        "initial_data_qubits" => length(exp.active0),
        "ancilla_qubits" => max_anc,
        "initial_ancilla_qubits" => length(exp.prep_checks),
        # extra vs the idle code, which already has one ancilla per check
        "extra_qubits" => (length(active) - length(exp.active0)) + max(0, max_anc - length(exp.prep_checks)),
        "circuit_layers" => total_layers,
    )
    return BuildResult(c, n_meas, tlist, det_count, cost, epoch_info, n_prep)
end

function _reduce_by_detectors(bits::Bits, dets::Vector{Bits})
    improved = true
    while improved
        improved = false
        for d in dets
            if overlaps(d, bits) && count_ones(xor(bits, d)) < count_ones(bits)
                bits = xor(bits, d)
                improved = true
            end
        end
    end
    return bits
end

function _check_state_generators(alg::PauliAlgebra, lstate)
    vecs = [from_dict(alg, g) for (g, _) in lstate]
    for i in eachindex(vecs), j in 1:i-1
        anticommutes(alg, vecs[i], vecs[j]) && throw(ArgumentError("state generators $(i-1) and $(j-1) anticommute"))
    end
end

function _logical_coords(alg, v, ops, k)
    return vcat([anticommutes(alg, v, ops[k+i]) for i in 1:k], [anticommutes(alg, v, ops[i]) for i in 1:k])
end

function _combine(alg, d, ops, k)
    v = zero_pvec(alg)
    for i in 1:k
        d[i] && (v = xor(v, ops[i]))
        d[k+i] && (v = xor(v, ops[k+i]))
    end
    return to_dict(alg, v)
end

"""
    complete_generators(alg, gens, basis)

Logical Paulis `h` such that `gens ∪ h` is a full commuting independent set of `length(basis)` generators.
"""
function complete_generators(alg::PauliAlgebra, gens, basis)
    k = length(basis)
    ops = vcat([from_dict(alg, x) for (x, _) in basis], [from_dict(alg, z) for (_, z) in basis])
    rows = [_logical_coords(alg, from_dict(alg, g), ops, k) for g in gens]
    out = PauliDict[]
    while length(rows) < k
        if isempty(rows)
            cand = Matrix{Bool}(I2(2k))
            r0 = 0
        else
            G = Matrix{Bool}(permutedims(reduce(hcat, rows)))
            cand = nullspace_gf2([G[:, k+1:end] G[:, 1:k]])
            r0 = gf2rank(G)
        end
        found = false
        for i in axes(cand, 1)
            v = cand[i, :]
            M = isempty(rows) ? permutedims(v) : [permutedims(reduce(hcat, rows)); permutedims(v)]
            if gf2rank(M) > r0
                push!(rows, Vector{Bool}(v))
                push!(out, _combine(alg, v, ops, k))
                found = true
                break
            end
        end
        found || error("could not complete logical generators")
    end
    return out
end

I2(m) = [i == j for i in 1:m, j in 1:m]

function _destabilizers_for(alg::PauliAlgebra, lstate, exp::Experiment)
    exp.logical_basis === nothing && throw(ArgumentError("eigen_targets requires exp.logical_basis"))
    return destabilizers_from_basis(alg, [g for (g, _) in lstate], exp.logical_basis)
end

"""
    destabilizers_from_basis(alg, gens, basis)

Logical Paulis `d_j` anticommuting with `gens[i]` iff i == j, built from the symplectic basis.
"""
function destabilizers_from_basis(alg::PauliAlgebra, gens, basis)
    k = length(basis)
    ops = vcat([from_dict(alg, x) for (x, _) in basis], [from_dict(alg, z) for (_, z) in basis])
    m = length(gens)
    G = Matrix{Bool}(permutedims(reduce(hcat, [_logical_coords(alg, from_dict(alg, g), ops, k) for g in gens])))
    # solve D Omega G^T = I, Omega swapping the halves
    OmG = [G[:, k+1:end] G[:, 1:k]]
    out = PauliDict[]
    for j in 1:m
        d = gf2_solve(OmG, [i == j for i in 1:m])
        d === nothing && throw(ArgumentError("logical generators not independent"))
        push!(out, _combine(alg, d, ops, k))
    end
    return out
end
