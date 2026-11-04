# Surface-code local Pauli settings, including Y, built from per-patch fragments.
# Each patch is measured by its own fragment (Gidney's in-place Y measurement, arXiv:2302.07395, or
# an X/Z readout) after an ideal multi-patch state preparation. Gidney's generator (`midout`,
# vendored in python/gidney_inplace_y) runs through the bridge; parsing, noise and composition are here.

const VENDOR = joinpath(dirname(dirname(@__DIR__)), "python", "gidney_inplace_y")
# Smallest padding b >= d ÷ 2 with verified circuit distance d; at d = 7, b = 3 gives distance 6.
const Y_DEFAULT_ROUNDS = Dict(3 => 1, 5 => 2, 7 => 4)
const ANNOT = Set(["DETECTOR", "OBSERVABLE_INCLUDE", "SHIFT_COORDS", "QUBIT_COORDS", "TICK"])

"""
    Fragment

Noisy part of a single-patch destructive measurement on local qubit ids. Local records:
0..length(stabs)-1 are the ideal MPP outcomes of `stabs`, then the body's measurements.
`detectors` and `observable` are lists of local records; `basis` is (Xbar, Zbar).
"""
struct Fragment
    n_qubits::Int
    data::Vector{Int}
    coords::OrderedDict{Int,Vector{Float64}}
    stabs::Vector{PauliDict}
    ideal::Vector{Vector{Instr}}
    detectors::Vector{Vector{Int}}
    observable::Vector{Int}
    basis::Tuple{PauliDict,PauliDict}
    kind::Char
    rounds::Int
    name::String
    n_body_meas::Int
end

"""The measured logical (Y: phase-free i Xbar Zbar)."""
function fragment_target(f::Fragment)
    X, Z = f.basis
    f.kind == 'X' && return copy(X)
    f.kind == 'Z' && return copy(Z)
    return _y_of(X, Z)
end

_y_of(X, Z) = PauliDict(q => (haskey(X, q) && haskey(Z, q)) ? 'Y' : (haskey(X, q) ? 'X' : 'Z') for q in sort!(collect(union(keys(X), keys(Z)))))

const _GIDNEY = Dict{Tuple{Char,Int,Int,Int},String}()

"""
    gidney_text(kind, d, boundary_rounds, memory_rounds)

Flattened text of Gidney's noiseless circuit: ideal MPP init of stabilizers and observable, then
the body. Y: `Y_magic_measure`; X / Z: `memory_rounds - 1` standard rounds and a data-measuring round.
"""
function gidney_text(kind::Char, d::Integer, boundary_rounds::Integer, memory_rounds::Integer)
    get!(_GIDNEY, (kind, d, boundary_rounds, memory_rounds)) do
        sys = pyimport("sys")
        pyin(VENDOR, sys.path) || sys.path.insert(0, VENDOR)
        if kind == 'Y'
            make_circuit = pyimport("midout._make_circuit").make_circuit
            c = make_circuit(; basis="Y_magic_measure", distance=d, noise=pybuiltins.None,
                             boundary_rounds, memory_rounds, convert_to_cz=false)
        else
            memory_rounds >= 1 || throw(ArgumentError("X/Z fragments with rounds need memory_rounds >= 1"))
            gen = pyimport("midout.gen")
            patch = pyimport("midout.circuits.steps._patches").make_xtop_qubit_patch(; distance=d)
            # Build the Python set and dict in a fixed insertion order.
            keep = [q for q in patch.data_set if pyconvert(Float64, kind == 'X' ? q.real : q.imag) == 0]
            qs = pybuiltins.set(pylist(keep))
            obs_d = pydict()
            for q in qs
                obs_d[q] = string(kind)
            end
            obs = gen.PauliString(obs_d)
            mem = gen.standard_surface_code_chunk(patch; obs)
            chunks = pylist([mem.magic_init_chunk()])
            memory_rounds > 1 && chunks.append(mem.with_repetitions(memory_rounds - 1))
            chunks.append(gen.standard_surface_code_chunk(patch; measure_data_basis=string(kind), obs))
            c = gen.compile_chunks_into_circuit(chunks)
        end
        pyconvert(String, pybuiltins.str(c.flattened()))
    end
end

function _split_mpp(inst::Instr)
    prods = PauliDict[]
    cur = PauliDict()
    expect_more = false
    for t in inst.targets
        if t.kind == COMBINER
            expect_more = true
            continue
        end
        if !isempty(cur) && !expect_more
            push!(prods, cur)
            cur = PauliDict()
        end
        cur[t.value] = t.kind == PAULI_X ? 'X' : t.kind == PAULI_Y ? 'Y' : 'Z'
        expect_more = false
    end
    isempty(cur) || push!(prods, cur)
    return prods
end

"""
    gidney_patch(kind, d, boundary_rounds, memory_rounds)

Gidney's circuit split into `(coords, stabs, obs_head, moments, detectors, observable, n_body_meas)`,
with local records numbered stabs first, then body (the head observable record is dropped).
"""
function gidney_patch(kind::Char, d::Integer, boundary_rounds::Integer=0, memory_rounds::Integer=0)
    insts = parse_stim(gidney_text(kind, d, boundary_rounds, memory_rounds)).ops
    coords = OrderedDict{Int,Vector{Float64}}()
    k = 1
    while insts[k].name == "QUBIT_COORDS"
        for t in insts[k].targets
            coords[t.value] = copy(insts[k].args)
        end
        k += 1
    end
    head = insts[k]
    head.name == "MPP" || error("Gidney circuit does not start with an MPP")
    prods = _split_mpp(head)
    stabs, obs_head = prods[1:end-1], prods[end]
    k += 1
    insts[k].name == "OBSERVABLE_INCLUDE" || error("expected the head observable")
    k += 1
    n_head, n_stab = length(prods), length(prods) - 1
    meas = n_head
    loc(r) = r < n_stab ? r : r == n_stab ? error("body refers to the head observable record") : r - 1
    moments, cur, dets, obs = Vector{Instr}[], Instr[], Vector{Int}[], Int[]
    for inst in insts[k:end]
        if inst.name == "TICK"
            isempty(cur) || push!(moments, cur)
            cur = Instr[]
            continue
        end
        if inst.name in ("DETECTOR", "OBSERVABLE_INCLUDE")
            recs = sort([loc(meas + t.value) for t in inst.targets])
            inst.name == "DETECTOR" ? push!(dets, recs) : append!(obs, recs)
            continue
        end
        inst.name in ANNOT && continue
        push!(cur, inst)
        gate_kind(inst.name) == :measure && (meas += length(inst.targets))
    end
    isempty(cur) || push!(moments, cur)
    obs = sort([r for r in unique(obs) if isodd(count(==(r), obs))])
    return coords, stabs, obs_head, moments, dets, obs, meas - n_head
end

"""Gidney's logicals: Xbar = X on x = 0, Zbar = Z on y = 0."""
function _patch_basis(coords, d)
    at = Dict((round(Int, x * 2), round(Int, y * 2)) => q for (q, (x, y)) in coords)
    X = PauliDict(at[(0, 2j)] => 'X' for j in 0:d-1)
    Z = PauliDict(at[(2j, 0)] => 'Z' for j in 0:d-1)
    return X, Z
end

_data_qubits(coords) = sort([q for (q, (x, y)) in coords if isinteger(x) && isinteger(y)])
_n_rounds(moments) = count(m -> any(i -> gate_kind(i.name) == :measure, m), moments)

"""
    y_fragment(d; boundary_rounds=nothing, memory_rounds=0)

`memory_rounds` standard rounds, then Gidney's Y measurement (b + 2 rounds).
"""
function y_fragment(d::Integer; boundary_rounds=nothing, memory_rounds::Integer=0)
    b = boundary_rounds === nothing ? get(Y_DEFAULT_ROUNDS, d, d ÷ 2) : boundary_rounds
    coords, stabs, obs_head, moments, dets, obs, n_body = gidney_patch('Y', d, b, memory_rounds)
    X, Z = _patch_basis(coords, d)
    _y_of(X, Z) == obs_head || error("Gidney's observable is not i Xbar Zbar on the parsed logicals")
    @assert _n_rounds(moments) == memory_rounds + b + 2
    return Fragment(maximum(keys(coords)) + 1, _data_qubits(coords), coords, stabs, moments, dets, obs, (X, Z), 'Y',
                    memory_rounds + b + 2, "gidney_Y_d$(d)_b$(b)_m$(memory_rounds)", n_body)
end

"""
    xz_fragment(d, kind; memory_rounds=0)

X / Z measurement: transversal if `memory_rounds == 0`, else that many rounds, the last reading the data.
"""
function xz_fragment(d::Integer, kind::Char; memory_rounds::Integer=0)
    memory_rounds == 0 && return readout_fragment(d, kind)
    coords, stabs, obs_head, moments, dets, obs, n_body = gidney_patch(kind, d, 0, memory_rounds)
    X, Z = _patch_basis(coords, d)
    obs_head == (kind == 'X' ? X : Z) || error("Gidney's observable is not the parsed logical")
    @assert _n_rounds(moments) == memory_rounds
    return Fragment(maximum(keys(coords)) + 1, _data_qubits(coords), coords, stabs, moments, dets, obs, (X, Z), kind,
                    memory_rounds, "gidney_$(kind)_d$(d)_m$(memory_rounds)", n_body)
end

"""
    readout_fragment(d, kind)

Transversal destructive readout of the patch in X or Z (no syndrome rounds).
"""
function readout_fragment(d::Integer, kind::Char)
    coords, stabs, _, _, _, _, _ = gidney_patch('Y', d, d ÷ 2, 0)
    X, Z = _patch_basis(coords, d)
    data = _data_qubits(coords)
    n_stab = length(stabs)
    recof = Dict(q => n_stab + i - 1 for (i, q) in enumerate(data))
    moments = [[Instr(kind == 'X' ? "MX" : "M", qubit.(data), Float64[])]]
    dets = [sort(vcat([j - 1], [recof[q] for q in keys(s)])) for (j, s) in enumerate(stabs) if Set(values(s)) == Set([kind])]
    L = kind == 'X' ? X : Z
    return Fragment(maximum(keys(coords)) + 1, data, coords, stabs, moments, dets, sort([recof[q] for q in keys(L)]),
                    (X, Z), kind, 0, "readout_$(kind)_d$(d)", length(data))
end

"""
    apply_noise(moments, p, data; noise=DEFAULT_NOISE)

Noisy circuit of a fragment body and its number of moments. Per moment: DEPOLARIZE2 after
two-qubit gates, flips after resets and before measurements, DEPOLARIZE1 on untouched qubits and
after single-qubit gates.
"""
function apply_noise(moments, p::Float64, data::Set{Int}; noise::AbstractString=DEFAULT_NOISE)
    nm = NOISE_MODELS[noise]
    nm.round_data != 0 && throw(ArgumentError("noise model $(repr(noise)) (per-round data noise) is not supported for fragments"))
    out = Circuit()
    p_reset, p_meas, p_2q, p_idle = nm.reset * p, nm.meas * p, nm.gate2 * p, nm.idle * p
    allq = Set{Int}()
    for m in moments, inst in m
        union!(allq, (t.value for t in inst.targets if t.kind == QUBIT))
    end
    union!(allq, data)
    for m in moments
        pre, post, body = Circuit(), Circuit(), Circuit()
        touched = Set{Int}()
        for inst in m
            qs = [t.value for t in inst.targets]
            union!(touched, qs)
            append!(body, inst.name, copy(inst.targets), inst.args)
            p == 0 && continue
            k = gate_kind(inst.name)
            if k == :measure
                p_meas > 0 && append!(pre, inst.name == "MX" ? "Z_ERROR" : "X_ERROR", qs, [p_meas])
            elseif k == :reset
                p_reset > 0 && append!(post, inst.name == "RX" ? "Z_ERROR" : "X_ERROR", qs, [p_reset])
            elseif k == :two
                p_2q > 0 && append!(post, "DEPOLARIZE2", qs, [p_2q])
            else
                p_idle > 0 && append!(post, "DEPOLARIZE1", qs, [p_idle])
            end
        end
        if p > 0 && p_idle > 0
            idle = setdiff(nm.idle_ancilla != 0 ? allq : data, touched)
            isempty(idle) || append!(post, "DEPOLARIZE1", sort!(collect(idle)), [p_idle])
        end
        append!(out, pre)
        append!(out, body)
        append!(out, post)
        append!(out, "TICK")
    end
    return out, length(moments)
end

function _shifted(c::Circuit, o::Int)
    o == 0 && return c
    return Circuit([Instr(i.name, [t.kind == QUBIT ? qubit(t.value + o) : t for t in i.targets], i.args) for i in c.ops])
end

"""
    ComposedExperiment(fragments, state; name="composed")

Patches each measured by their own fragment after ideal preparation of `state = [(label, sign bit)]`.
Costs are reported as if the fragments ran in parallel (rounds and layers: max; qubits: sum).
"""
struct ComposedExperiment
    fragments::Vector{Fragment}
    state::Vector{Tuple{String,Int}}
    name::String
    offsets::Vector{Int}
    logical_basis::LogicalBasis
    n_qubits::Int
    targets::Vector{PauliDict}
    logical_state::StateGens
end

function ComposedExperiment(fragments::Vector{Fragment}, state; name="composed")
    off = 0
    offsets, basis = Int[], LogicalBasis()
    sh(dct, o) = PauliDict(q + o => s for (q, s) in dct)
    for f in fragments
        push!(offsets, off)
        push!(basis, (sh(f.basis[1], off), sh(f.basis[2], off)))
        off += f.n_qubits
    end
    targets = [sh(fragment_target(f), o) for (f, o) in zip(fragments, offsets)]
    st = [(String(lab), Int(s)) for (lab, s) in state]
    lstate = isempty(st) ? StateGens() : physical_state(st, basis, off)
    k = length(fragments)
    for (i, f) in enumerate(fragments)
        lab = join(j == i ? f.kind : 'I' for j in 1:k)
        d, sb = logical_to_physical(lab, basis, off)
        (d == targets[i] && sb == 0) || error("fragment $(i - 1) target is not +$lab")
    end
    return ComposedExperiment(fragments, st, name, offsets, basis, off, targets, lstate)
end

function build(ce::ComposedExperiment, p::Real; eigen_targets::Bool=false, noise::AbstractString=DEFAULT_NOISE)
    p = Float64(p)
    n = ce.n_qubits
    alg = PauliAlgebra(n)
    c = Circuit()
    for (i, (f, o)) in enumerate(zip(ce.fragments, ce.offsets))
        shift = (i - 1) * (maximum(xy[1] for xy in values(f.coords)) + 3)
        for (q, xy) in sort!(collect(f.coords); by=first)
            append!(c, "QUBIT_COORDS", [q + o], [xy[1] + shift, xy[2]])
        end
    end
    datas = sort([q + o for (f, o) in zip(ce.fragments, ce.offsets) for q in f.data])
    append!(c, "R", datas)
    n_meas = 0
    stab_rec = Vector{Int}[]
    for (f, o) in zip(ce.fragments, ce.offsets)
        recs = Int[]
        for s in f.stabs
            append!(c, "MPP", mpp_targets(PauliDict(q + o => b for (q, b) in s)))
            push!(recs, n_meas)
            n_meas += 1
        end
        push!(stab_rec, recs)
    end
    lstate = eigen_targets ? StateGens([(t, 0) for t in ce.targets]) : ce.logical_state
    gens = [g for (g, _) in lstate]
    destabs = isempty(gens) ? PauliDict[] : destabilizers_from_basis(alg, gens, ce.logical_basis)
    lrec = Int[]
    for (lp, _) in lstate
        append!(c, "MPP", mpp_targets(lp))
        push!(lrec, n_meas)
        n_meas += 1
    end
    for ((lp, sign), dstab, r) in zip(lstate, destabs, lrec), (q, s) in sortedpauli(dstab)
        append!(c, PAULI_GATE[s], [rec(r - n_meas), qubit(q)])
        sign != 0 && append!(c, string(s), [q])
    end
    extra = complete_generators(alg, gens, ce.logical_basis)
    if !isempty(extra)
        dmix = destabilizers_from_basis(alg, vcat(gens, extra), ce.logical_basis)[length(gens)+1:end]
        for (h, d_) in zip(extra, dmix), op in (d_, h)
            append!(c, "MPP", mpp_targets(op))
            n_meas += 1
        end
    end
    append!(c, "TICK")
    n_prep = n_meas
    dets, targets = Vector{Int}[], Vector{Int}[]
    layers = anc_peak = rounds = 0
    for (i, (f, o)) in enumerate(zip(ce.fragments, ce.offsets))
        body, nmom = apply_noise(f.ideal, p, Set(f.data); noise)
        append!(c, _shifted(body, o))
        ns = length(f.stabs)
        base, sr = n_meas, stab_rec[i]
        g(j) = j < ns ? sr[j+1] : base + j - ns
        append!(dets, [[g(j) for j in dd] for dd in f.detectors])
        push!(targets, sort([g(j) for j in f.observable]))
        n_meas += f.n_body_meas
        layers = max(layers, f.rounds > 0 ? nmom : 0)
        rounds = max(rounds, f.rounds)
        f.rounds > 0 && (anc_peak += f.n_qubits - length(f.data))
    end
    for dd in dets
        append!(c, "DETECTOR", [rec(j - n_meas) for j in dd])
    end
    for (t_idx, recs) in enumerate(targets)
        append!(c, "OBSERVABLE_INCLUDE", [rec(j - n_meas) for j in recs], [t_idx - 1])
    end
    n_data = length(datas)
    cost = OrderedDict{String,Int}(
        "merged_rounds" => rounds, "data_qubits" => n_data, "initial_data_qubits" => n_data,
        "ancilla_qubits" => anc_peak, "initial_ancilla_qubits" => sum(length(f.stabs) for f in ce.fragments),
        "extra_qubits" => sum((max(0, f.n_qubits - length(f.data) - length(f.stabs)) for f in ce.fragments if f.rounds > 0); init=0),
        "circuit_layers" => layers)
    return BuildResult(c, n_meas, targets, length(dets), cost, OrderedDict{String,Any}[], n_prep)
end

"""
    surface_local_experiment(d, bases, state; boundary_rounds=nothing, memory_rounds=0)

Patch i measured in `bases[i]` ∈ {X, Y, Z} after `memory_rounds` noisy rounds; returns `(experiment, labels)`.
"""
function surface_local_experiment(d::Integer, bases::AbstractString, state; boundary_rounds=nothing, memory_rounds::Integer=0)
    frs = [b == 'Y' ? y_fragment(d; boundary_rounds, memory_rounds) : xz_fragment(d, b; memory_rounds) for b in bases]
    k = length(bases)
    labels = [join(j == i ? b : 'I' for j in 1:k) for (i, b) in enumerate(bases)]
    b_used = boundary_rounds === nothing ? get(Y_DEFAULT_ROUNDS, d, d ÷ 2) : boundary_rounds
    return ComposedExperiment(frs, collect(state); name="surface_local_$(bases)_d$(d)_b$(b_used)_m$(memory_rounds)"), labels
end
