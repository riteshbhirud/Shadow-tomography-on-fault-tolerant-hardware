# Detector error models: read from stim (bridged) into Julia, then turned into check/observable
# matrices (hypergraph decoders) or graphlike edges (matching).

"""
    DemError

One `error(p)` instruction of a flattened DEM: probability `p` and its components (split at `^`),
each a list of detector ids (0-based, in target order) and an observable bit mask.
"""
struct DemError
    p::Float64
    components::Vector{Tuple{Vector{Int},UInt64}}
end

"""
    Dem

A flattened detector error model: `num_detectors`, `num_observables` and the `errors` in file order.
"""
struct Dem
    num_detectors::Int
    num_observables::Int
    errors::Vector{DemError}
end

"""
    detector_error_model(circuit; decompose_errors=false) -> Py

stim's DEM of `circuit` (stim text or a `stim.Circuit`).
"""
function detector_error_model(circuit; decompose_errors::Bool = false)
    c = circuit isa AbstractString ? pymodule(:stim).Circuit(circuit) : circuit
    return c.detector_error_model(; decompose_errors)
end

"""
    read_dem(dem) -> Dem

Read a stim `DetectorErrorModel` (a `Py`, or DEM text). Targets are parsed from the flattened text;
probabilities are taken from stim's API because the text form rounds them.
"""
function read_dem(dem::Py)
    flat = dem.flattened()
    probs = Float64[]
    for inst in flat
        pyconvert(String, inst.type) == "error" && push!(probs, pyconvert(Float64, inst.args_copy()[0]))
    end
    errors = DemError[]
    for line in eachline(IOBuffer(pyconvert(String, pystr(flat))))
        line = strip(line)
        startswith(line, "error") || continue
        comps = Tuple{Vector{Int},UInt64}[]
        dets, om = Int[], UInt64(0)
        for tok in split(line[(findfirst(')', line) + 1):end])
            if tok == "^"
                push!(comps, (dets, om))
                dets, om = Int[], UInt64(0)
            elseif tok[1] == 'D'
                push!(dets, parse(Int, tok[2:end]))
            elseif tok[1] == 'L'
                k = parse(Int, tok[2:end])
                k < 64 || error("more than 64 observables are not supported")
                om ⊻= UInt64(1) << k
            end
        end
        push!(comps, (dets, om))
        push!(errors, DemError(probs[length(errors) + 1], comps))
    end
    length(errors) == length(probs) || error("DEM text and instructions disagree")
    return Dem(pyconvert(Int, dem.num_detectors), pyconvert(Int, dem.num_observables), errors)
end
read_dem(text::AbstractString) = read_dem(pymodule(:stim).DetectorErrorModel(text))

# Insertion-ordered map key => probability, merging repeats as independent mechanisms
# (Python dicts keep insertion order; the decoders are built in that order).
struct _MergedMechanisms{K}
    index::Dict{K,Int}
    keys::Vector{K}
    p::Vector{Float64}
end
_MergedMechanisms{K}() where {K} = _MergedMechanisms{K}(Dict{K,Int}(), K[], Float64[])

function _merge!(m::_MergedMechanisms, key, pr::Float64)
    i = get(m.index, key, 0)
    if i == 0
        push!(m.keys, key)
        push!(m.p, pr)
        m.index[key] = length(m.keys)
    else
        q = m.p[i]
        m.p[i] = q * (1 - pr) + pr * (1 - q)
    end
    return m
end

"""
    dem_to_matrices(dem::Dem) -> (H, L, p)

Check matrix `H` (detectors × mechanisms), observable matrix `L` (observables × mechanisms) and priors
`p`. Mechanisms with identical detector and observable sets are merged, in first-occurrence order.
"""
function dem_to_matrices(dem::Dem)
    cols = _MergedMechanisms{Tuple{Vector{Int},Vector{Int}}}()
    for err in dem.errors
        dets, obs = Set{Int}(), Set{Int}()
        for (ds, om) in err.components
            for d in ds
                d in dets ? delete!(dets, d) : push!(dets, d)
            end
            for k in 0:63
                (om >> k) & 1 == 1 && (k in obs ? delete!(obs, k) : push!(obs, k))
            end
        end
        _merge!(cols, (sort!(collect(dets)), sort!(collect(obs))), err.p)
    end
    n = length(cols.keys)
    Hi, Hj, Li, Lj = Int[], Int[], Int[], Int[]
    for (j, (d, o)) in enumerate(cols.keys)
        append!(Hi, d .+ 1)
        append!(Hj, fill(j, length(d)))
        append!(Li, o .+ 1)
        append!(Lj, fill(j, length(o)))
    end
    H = sparse(Hi, Hj, ones(UInt8, length(Hi)), dem.num_detectors, n)
    L = sparse(Li, Lj, ones(UInt8, length(Li)), dem.num_observables, n)
    return H, L, cols.p
end

"""
    graphlike_edges(dem::Dem; undetectable=nothing)

Graphlike components of a decomposed DEM as an insertion-ordered list of `(u, v, obs_mask) => p`
(`v = -1` for a boundary edge), merging identical components as independent mechanisms. Components
without detectors that flip observables are merged into `undetectable` (obs_mask => p) if given.
"""
function graphlike_edges(dem::Dem; undetectable::Union{Nothing,_MergedMechanisms{UInt64}} = nothing)
    out = _MergedMechanisms{Tuple{Int,Int,UInt64}}()
    for err in dem.errors, (dets, om) in err.components
        if isempty(dets)
            om != 0 && undetectable !== nothing && _merge!(undetectable, om, err.p)
            continue
        end
        length(dets) > 2 && error("non-graphlike component; decompose_errors=true required")
        u, v = dets[1], length(dets) == 2 ? dets[2] : -1
        v != -1 && v < u && ((u, v) = (v, u))
        _merge!(out, (u, v, om), err.p)
    end
    return out
end
