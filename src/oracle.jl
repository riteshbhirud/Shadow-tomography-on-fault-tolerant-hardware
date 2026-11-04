# Measurement oracle: `run_oracle` samples and decodes one measurement setting on one input state
# Sampling is stim's, so the same circuit and seed give the same shots as the reference.

using JSON

"""
    OracleCircuit(name, circuit, n_prep, target_records, detector_count, cost)

What the oracle needs from one build of an experiment:
- `name`: experiment name (recorded in `meta`; the default setting label);
- `circuit`: stim circuit text, with `DETECTOR`s and one `OBSERVABLE_INCLUDE` per target;
- `n_prep`: number of noiseless preparation measurement records at the start;
- `target_records`: per target, the 0-based absolute measurement-record ids whose parity is the raw outcome;
- `detector_count`;
- `cost`: string-keyed dict, written to the records in its iteration order.
"""
struct OracleCircuit
    name::String
    circuit::String
    n_prep::Int
    target_records::Vector{Vector{Int}}
    detector_count::Int
    cost::AbstractDict
end

"""
    oracle_circuit(exp, p; eigen_targets=false, noise="bravyi_cal") -> OracleCircuit

Build `exp` at noise strength `p`. This is the interface an experiment type implements to be run by
[`run_oracle`](@ref). With `eigen_targets=true` the targets are prepared in eigenstates, so the DEM has
deterministic observables; that build must have the same detectors and target records (relative to its
own preparation records) as the sampled one.
"""
function oracle_circuit end

"""
    PrebuiltExperiment(name, builds)

An experiment given by finished builds: `builds[(p, eigen_targets, noise)]` is an [`OracleCircuit`](@ref).
"""
struct PrebuiltExperiment
    name::String
    builds::Dict{Tuple{Float64,Bool,String},OracleCircuit}
end

function oracle_circuit(exp::PrebuiltExperiment, p; eigen_targets::Bool = false, noise::AbstractString = "bravyi_cal")
    key = (Float64(p), eigen_targets, String(noise))
    haskey(exp.builds, key) || error("$(exp.name): no build for p = $p, eigen_targets = $eigen_targets, noise = $noise")
    return exp.builds[key]
end

"""
    load_prebuilt(path) -> PrebuiltExperiment

Read an exported reference experiment (keys `sample` and, for
p > 0, `eigen`).
"""
function load_prebuilt(path::AbstractString)
    d = JSON.parsefile(path)
    builds = Dict{Tuple{Float64,Bool,String},OracleCircuit}()
    for (key, eigen) in (("sample", false), ("eigen", true))
        haskey(d, key) || continue
        b = d[key]
        builds[(Float64(d["p"]), eigen, String(d["noise"]))] = OracleCircuit(d["experiment_name"], b["circuit"],
            b["n_prep"], [Vector{Int}(t) for t in b["target_records"]], b["detector_count"], b["cost"])
    end
    return PrebuiltExperiment(d["experiment_name"], builds)
end

"""
    OracleResult

Outcome of one [`run_oracle`](@ref) call: `setting` (measured Pauli labels), `outcomes` and `raw`
(shots × targets bits; 0 means eigenvalue +1; `raw` is undecoded), per-shot `confidence` and
`confidence_kind` (see [`CONF_KINDS`](@ref)), the setting's `cost`, and `meta`.
"""
struct OracleResult
    setting::Vector{String}
    outcomes::Matrix{UInt8}
    raw::Matrix{UInt8}
    confidence::Vector{Float64}
    confidence_kind::Vector{String}
    cost::AbstractDict
    meta::JSON.Object{String,Any}
end

"""
    records(res::OracleResult) -> Vector

One ordered dict per shot (keys `setting`, `outcomes`, `confidence`, `confidence_kind`, `cost`, `meta`);
write them with [`write_records_raw`](@ref) or [`pyjson`](@ref). This is not the dataset record format
(no `shot` or `setting_id`, unrounded confidence), and an infinite confidence is written as `Infinity`,
which JSON readers, including [`load_records`](@ref), reject.
"""
function records(res::OracleResult)
    return map(1:size(res.outcomes, 1)) do s
        r = JSON.Object{String,Any}()
        r["setting"] = res.setting
        r["outcomes"] = Int.(res.outcomes[s, :])
        r["confidence"] = res.confidence[s]
        r["confidence_kind"] = res.confidence_kind[s]
        r["cost"] = res.cost
        r["meta"] = res.meta
        r
    end
end

# One decoder is kept between runs; Relay-BP decoders are rebuilt every run.
const _DEC_CACHE = Ref{Any}(nothing)

_relative(b::OracleCircuit) = [[i >= b.n_prep ? i - b.n_prep : -1 - i for i in t] for t in b.target_records]

function _decoder(be::OracleCircuit, decoder::AbstractString, decoder_kwargs)
    key = (be.circuit, decoder, repr(decoder_kwargs))
    cached = _DEC_CACHE[]
    if cached === nothing || cached[1] != key || decoder in STATEFUL_DECODERS
        dem = detector_error_model(be.circuit; decompose_errors = decoder == "pymatching")
        kw = decoder_kwargs === nothing ? (;) : _kwargs(decoder_kwargs)
        _DEC_CACHE[] = cached = (key, DECODERS[decoder](dem; kw...))
    end
    return cached[2]
end

"""
    run_oracle(exp, p, shots; decoder="pymatching", seed=0, setting_labels=nothing, decoder_kwargs=nothing,
               batch=20000, noise="bravyi_cal", with_gap=true) -> OracleResult

Sample `shots` shots of `exp` at noise strength `p` with stim (seed `seed`) and decode them.
The decoder is built from the eigen-target build (see [`oracle_circuit`](@ref));
the sampled circuit uses the requested logical state. `decoder` is a key of [`DECODERS`](@ref);
`decoder_kwargs` (a dict or named tuple) go to its constructor. With `with_gap=false` only predictions
are computed (confidence `NaN`). At `p = 0` no decoder is used and every detector must stay silent.
"""
function run_oracle(exp, p::Real, shots::Integer; decoder::AbstractString = "pymatching", seed::Integer = 0,
                    setting_labels = nothing, decoder_kwargs = nothing, batch::Integer = 20000,
                    noise::AbstractString = "bravyi_cal", with_gap::Bool = true)
    b = oracle_circuit(exp, p; noise)
    dec = nothing
    if p > 0
        be = oracle_circuit(exp, p; eigen_targets = true, noise)
        (be.detector_count != b.detector_count || _relative(be) != _relative(b)) &&
            error("DEM circuit and sampling circuit disagree on detectors/targets")
        dec = _decoder(be, decoder, decoder_kwargs)
    end
    stim = pymodule(:stim)
    circuit = stim.Circuit(b.circuit)
    sampler = circuit.compile_sampler(; seed)
    conv = circuit.compile_m2d_converter()
    outs, raws, confs, kinds = Matrix{UInt8}[], Matrix{UInt8}[], Vector{Float64}[], Vector{String}[]
    left = shots
    while left > 0
        k = min(batch, left)
        left -= k
        m = sampler.sample(k)
        dets = pyconvert(Matrix{Bool}, conv.convert(; measurements = m, separate_observables = true)[0])
        M = pyconvert(Matrix{Bool}, m)
        raw = UInt8[reduce(xor, (M[s, i + 1] for i in recs); init = false) for s in 1:k, recs in b.target_records]
        if p > 0
            pred, conf, kind = decode(dec, dets; with_gap)
        else
            any(dets) && error("noiseless circuit fired a detector")
            pred, conf, kind = zeros(UInt8, size(raw)), fill(Inf, k), fill("noiseless", k)
        end
        push!(outs, raw .⊻ pred)
        push!(raws, raw)
        push!(confs, conf)
        push!(kinds, kind)
    end
    nt = length(b.target_records)
    meta = JSON.Object{String,Any}()
    for (k, v) in (("p", p), ("noise", String(noise)), ("shots", shots), ("decoder", String(decoder)),
                   ("seed", seed), ("detectors", b.detector_count), ("name", b.name))
        meta[k] = v
    end
    setting = setting_labels === nothing || isempty(setting_labels) ? [b.name] : String.(collect(setting_labels))
    return OracleResult(setting, reduce(vcat, outs; init = zeros(UInt8, 0, nt)),
        reduce(vcat, raws; init = zeros(UInt8, 0, nt)), reduce(vcat, confs; init = Float64[]),
        reduce(vcat, kinds; init = String[]), b.cost, meta)
end
