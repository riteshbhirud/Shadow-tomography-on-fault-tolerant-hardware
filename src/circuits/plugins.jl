# Code plug-ins: per code, its logical basis and a catalog of native measurement
# operations. `build_setting` turns a setting id into an experiment; the oracle samples and decodes it.

"""
    NativeOp(name, kind, describe, builder, params)

One native operation; `builder(state; setting_id, params...)` returns `(experiment, labels, decoder)`.
"""
struct NativeOp
    name::String
    kind::String  # readout | local | surgery | gauging
    describe::String
    builder::Function
    params::Dict{String,String}
end
NativeOp(name, kind, describe, builder) = NativeOp(name, kind, describe, builder, Dict{String,String}())

"""
    CodePlugin

A code: its logical basis (a function of parameters), schedule description, native operations
and default decoder.
"""
struct CodePlugin
    name::String
    describe::String
    logical_basis::Function
    schedule::String
    ops::OrderedDict{String,NativeOp}
    decoder::String
end

"""
    op(plugin, setting_id)

The operation for a setting id; keys ending in '*' match every id with that prefix.
"""
function op(pl::CodePlugin, name::AbstractString)
    haskey(pl.ops, name) && return pl.ops[name]
    for (key, o) in pl.ops
        endswith(key, "*") && startswith(name, key[1:end-1]) && return o
    end
    throw(KeyError("$(pl.name): no native op $(repr(name)); available: $(sort(collect(keys(pl.ops))))"))
end

_setting_arg(setting_id) = split(setting_id, ':')[2]

_surface_local(state; setting_id, d=3, memory_rounds=0, boundary_rounds=nothing) =
    (surface_local_experiment(d, _setting_arg(setting_id), state; boundary_rounds, memory_rounds)..., "pymatching")

_surface_surgery(state; setting_id, d=3, n_patches=2, rounds=3) =
    (setting_experiment(setting_id, d, n_patches, state, rounds)..., "pymatching")

"""Per-patch (Xbar, Zbar) in the local numbering of the fragment measuring `kind`."""
function _surface_basis(; d=3, n_patches=1, kind='Y', memory_rounds=0)
    f = kind == 'Y' ? y_fragment(d; memory_rounds) : xz_fragment(d, kind; memory_rounds)
    return fill(f.basis, n_patches)
end

_bb_readout(state; setting_id) = (bb72_readout_experiment(setting_id[end], state)..., "relay")

function _bb_gauge(state; setting_id, rounds=3)
    arg = _setting_arg(setting_id)
    P, i = arg[1], parse(Int, arg[2:end])
    i == 0 || throw(ArgumentError("oracle-basis gauging is set up for logical qubit 0; use local:<bases> for qubits 1, 2"))
    return bb_gauging(P, rounds; state_labels=state), [P * "I"^11], "relay"
end

function _bb_local(state; setting_id, rounds=3)
    exp, labels, _ = local_experiment(_setting_arg(setting_id), state, rounds)
    return exp, labels, "relay"
end

_bb_basis(; kind="oracle") = kind == "local" ? local_basis() : bb_setup().basis

const SURFACE = CodePlugin(
    "surface", "rotated surface code, distance d, one patch per logical qubit", _surface_basis,
    "standard rounds (Gidney's generator for local settings; builder groups X then Z checks for surgery)",
    OrderedDict(
        "local:*" => NativeOp("local:<B0B1..>", "local",
                              "patch i measured in X/Y/Z after memory_rounds rounds; Y by Gidney's in-place measurement",
                              _surface_local,
                              Dict("d" => "distance", "memory_rounds" => "rounds before the measurement",
                                   "boundary_rounds" => "Y padding rounds b (default: smallest b >= d ÷ 2 with circuit distance d)")),
        "surgery:*" => NativeOp("surgery:ZZ|XX", "surgery", "lattice surgery of neighbouring patches, conjugate termination",
                                _surface_surgery, Dict("d" => "distance", "n_patches" => "patches", "rounds" => "r")),
    ),
    "pymatching")

const BB72 = CodePlugin(
    "bb72", "[[72,12,6]] bivariate bicycle code (Bravyi et al.), one block of 12 logical qubits", _bb_basis,
    "translation-invariant CNOT orders, X checks then Z checks (ordered_checks); gauging: greedy grouping",
    OrderedDict(
        "readout:*" => NativeOp("readout:X|Z", "readout", "destructive readout of all 12 logicals in X or Z (oracle basis)", _bb_readout),
        "gauge:*" => NativeOp("gauge:<P>0", "gauging", "Williamson-Yoder gauging of Pbar_0 (oracle basis), r rounds", _bb_gauge,
                              Dict("rounds" => "r")),
        "local:*" => NativeOp("local:<B0..>", "local",
                              "local Pauli setting on qubits 0..m-1 (local basis): readout + simultaneous gauging", _bb_local,
                              Dict("rounds" => "r")),
    ),
    "relay")

"""Code plug-ins by name ("surface", "bb72")."""
const REGISTRY = OrderedDict(p.name => p for p in (SURFACE, BB72))

"""
    build_setting(code, setting_id, state; params...)

`(experiment, labels, decoder)` of one native setting on the logical state `[(label, sign bit)]`,
e.g. `build_setting("surface", "local:XY", state; d=3, memory_rounds=3)` or
`build_setting("bb72", "gauge:Y0", state; rounds=3)`.
"""
build_setting(code::AbstractString, setting_id::AbstractString, state; params...) =
    op(REGISTRY[code], setting_id).builder(state; setting_id, params...)
