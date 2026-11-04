# BP-type decoders on hypergraph DEMs (BB codes): ldpc BP-OSD and Relay-BP, bridged; the
# hypothesis split, the plain-anchored gap and the fallbacks are done here.

# What every hypergraph decoder keeps about its DEM.
struct _Hyper
    H::SparseMatrixCSC{UInt8,Int}
    L::SparseMatrixCSC{UInt8,Int}
    p::Vector{Float64}
    llr::Vector{Float64}
    Lrows::Vector{Vector{Int}}  # mechanisms flipping each observable
end

function _Hyper(dem::Dem)
    H, L, p = dem_to_matrices(dem)
    Lt = sparse(transpose(L))
    Lrows = [rowvals(Lt)[nzrange(Lt, j)] for j in 1:size(L, 1)]
    return _Hyper(H, L, p, llr(Vector{Float64}(p)), Lrows)
end

_no(h::_Hyper) = size(h.L, 1)
_obs_parity(h::_Hyper, j::Int, e::AbstractVector) = UInt8(count(i -> e[i] != 0, h.Lrows[j]) % 2)
_with_obs(h::_Hyper, j::Int) = vcat(h.H, h.L[j:j, :])
_np_priors(h::_Hyper) = _np().array(h.p; dtype = _np().float64)

# keyword arguments from a dict with string or symbol keys (decoder_kwargs read from JSON) or a named tuple
_kwargs(d) = NamedTuple(Symbol(k) => v for (k, v) in pairs(d))

# BP-type gaps are approximate; a decoder without its own kinds reports NaN as unavailable.
_default_kind(gaps) = [isnan(g) ? "unavailable" : "approx" for g in gaps]

#### BP-OSD

"""
    BPOSDGapDecoder(dem; max_iter=100, osd_order=7, osd_method="osd_cs")

ldpc BP-OSD (min-sum, scaling 0.625) with the hypothesis-split gap; a secondary decoder used as a
cross-check. With one observable both hypotheses are decoded on `[H; L]`; with several, each
observable's hypotheses are decoded on `[H; L_j]` and the gap is the minimum. If both hypotheses of an
observable fail, the prediction comes from a plain decode and the gap is 0.
"""
mutable struct BPOSDGapDecoder
    h::_Hyper
    kw::NamedTuple
    cache::Dict{Any,Py}
end

BPOSDGapDecoder(dem::Py; kw...) = BPOSDGapDecoder(read_dem(dem); kw...)
BPOSDGapDecoder(dem::Dem; max_iter = 100, osd_order = 7, osd_method = "osd_cs") =
    BPOSDGapDecoder(_Hyper(dem), (; max_iter, osd_order, osd_method), Dict{Any,Py}())

function _bposd_decode(dec::BPOSDGapDecoder, key, M::SparseMatrixCSC, syn::Vector{UInt8})
    d = get!(dec.cache, key) do
        pymodule(:ldpc).BpOsdDecoder(_scipy_sparse(M, :csc); error_channel = pylist(dec.h.p),
            max_iter = dec.kw.max_iter, bp_method = "minimum_sum", ms_scaling_factor = 0.625,
            osd_method = dec.kw.osd_method, osd_order = dec.kw.osd_order)
    end
    return _from_numpy_u8(d.decode(to_numpy_u8(syn)))
end

# Correction for (det_row, observable values) on M, or nothing if BP-OSD's output misses the syndrome.
function _bposd_constrained(dec, key, M, syn)
    e = _bposd_decode(dec, key, M, syn)
    return _syndrome_ok(M, e, syn) ? e : nothing
end

function decode(dec::BPOSDGapDecoder, dets::AbstractMatrix; with_gap::Bool = true)
    h = dec.h
    D = convert(Matrix{UInt8}, dets)
    n, no = size(D, 1), _no(h)
    pred = zeros(UInt8, n, no)
    gaps = fill(Inf, n)
    _w(e) = e === nothing ? Inf : _weight(h.llr, e)
    if no == 1
        Hext = vcat(h.H, h.L)
        for s in 1:n
            w = [_w(_bposd_constrained(dec, "ext", Hext, vcat(D[s, :], UInt8(o)))) for o in 0:1]
            pred[s, 1] = w[2] < w[1]
            gaps[s] = abs(w[2] - w[1])
        end
        return pred, gaps, _default_kind(gaps)
    end
    for s in 1:n
        g = Inf
        e0 = nothing
        for j in 1:no
            Hj = _with_obs(h, j)
            w = [_w(_bposd_constrained(dec, ("flip", j - 1), Hj, vcat(D[s, :], UInt8(o)))) for o in 0:1]
            if isinf(w[1]) && isinf(w[2])
                # Python's int(L[j] @ e0 % 2) raises on numpy >= 2.5 here; this is the intended value
                e0 === nothing && (e0 = _bposd_decode(dec, "H", h.H, D[s, :]))
                pred[s, j] = _obs_parity(h, j, e0)
                g = 0.0
            else
                pred[s, j] = w[2] < w[1]
                g = min(g, abs(w[2] - w[1]))
            end
        end
        gaps[s] = g
    end
    return pred, gaps, _default_kind(gaps)
end

#### Relay-BP

const RELAY_DEFAULTS = (gamma0 = 0.125, pre_iter = 80, num_sets = 100, set_max_iter = 60,
    gamma_dist_interval = (-0.24, 0.66), stop_nconv = 5, seed = 0)

function _relay(h::_Hyper, M::SparseMatrixCSC, kw)
    k = merge(RELAY_DEFAULTS, kw)
    return pymodule(:relay_bp).RelayDecoderF64(_scipy_sparse(M, :csr); error_priors = _np_priors(h),
        gamma0 = k.gamma0, pre_iter = k.pre_iter, num_sets = k.num_sets, set_max_iter = k.set_max_iter,
        gamma_dist_interval = pytuple(k.gamma_dist_interval), stop_nconv = k.stop_nconv, seed = k.seed)
end

# (success, correction) per shot of one detailed batch decode
function _relay_detailed(dec::Py, syn::Matrix{UInt8})
    res = dec.decode_detailed_batch(to_numpy_u8(syn))
    return [(pyconvert(Bool, r.success), _from_numpy_u8(r.decoding)) for r in res]
end

_relay_weights(h::_Hyper, dec::Py, syn) = [ok ? _weight(h.llr, e) : Inf for (ok, e) in _relay_detailed(dec, syn)]

function _relay_plain_pred(h::_Hyper, plain::Py, D::Matrix{UInt8}, j::Int)
    E = reshape(_from_numpy_u8(plain.decode_batch(to_numpy_u8(D))), size(D, 1), :)
    return [_obs_parity(h, j, view(E, s, :)) for s in 1:size(D, 1)]
end

"""
    RelayGapDecoder(dem; gamma0=0.125, pre_iter=80, num_sets=100, set_max_iter=60,
                    gamma_dist_interval=(-0.24, 0.66), stop_nconv=5, seed=0)

Relay-BP with the hypothesis-split gap (`relay_split`): for each observable both values are decoded on
`[H; L_j]`; prediction = lighter hypothesis, gap = |w₁ - w₀| (minimum over observables). A hypothesis
that does not converge weighs `Inf`; if both fail, the plain Relay-BP prediction is used and the gap is
0 (counted in `failures`).
"""
mutable struct RelayGapDecoder
    h::_Hyper
    decs::Vector{Py}
    plain::Py
    failures::Int
end

RelayGapDecoder(dem::Py; kw...) = RelayGapDecoder(read_dem(dem); kw...)
function RelayGapDecoder(dem::Dem; kw...)
    h = _Hyper(dem)
    k = NamedTuple(kw)
    return RelayGapDecoder(h, [_relay(h, _with_obs(h, j), k) for j in 1:_no(h)], _relay(h, h.H, k), 0)
end

function decode(dec::RelayGapDecoder, dets::AbstractMatrix; with_gap::Bool = true)
    h = dec.h
    D = convert(Matrix{UInt8}, dets)
    n, no = size(D, 1), _no(h)
    if !with_gap
        E = reshape(_from_numpy_u8(dec.plain.decode_batch(to_numpy_u8(D))), n, :)
        pred = UInt8[_obs_parity(h, j, view(E, s, :)) for s in 1:n, j in 1:no]
        return pred, fill(NaN, n), fill("unavailable", n)
    end
    pred = zeros(UInt8, n, no)
    gaps = fill(Inf, n)
    for j in 1:no
        w = [_relay_weights(h, dec.decs[j], hcat(D, fill(UInt8(o), n))) for o in 0:1]
        pred[:, j] .= w[2] .< w[1]
        both_fail = isinf.(w[1]) .& isinf.(w[2])
        dec.failures += count(both_fail)
        if any(both_fail)  # no preference: use plain Relay-BP rather than defaulting to 0
            pred[both_fail, j] .= _relay_plain_pred(h, dec.plain, D[both_fail, :], j)
        end
        g = abs.(w[2] .- w[1])
        g[both_fail] .= 0.0
        gaps = _npmin.(gaps, g)
    end
    return pred, gaps, _default_kind(gaps)
end

"""
    RelayAnchoredGapDecoder(dem; comp="relay", comp_kw=nothing, kw...)

Plain-anchored complementary gap (`relay`). Prediction and its weight come from one plain
Relay-BP decode, so predictions equal `with_gap=false`. For each observable one more decode on
`[H; L_j]` with the complementary value gives w_comp (`Inf` if it does not converge); the gap is
max(w_comp - w_pred, 0), minimum over observables. A plain decode that does not converge gives gap 0
(`decoder_failed`). The complementary decoder is Relay-BP with `kw` updated by `comp_kw` (default 20
relay sets) or, with `comp="bposd"`, ldpc BP-OSD (min-sum scaling 1.0, OSD-CS order 7, updated by
`comp_kw`). Negative gaps clipped to 0 are counted in `negative`, failed plain decodes in
`plain_failures`.
"""
mutable struct RelayAnchoredGapDecoder
    h::_Hyper
    plain::Py
    decs::Vector{Py}
    bposd::Vector{Py}
    comp::String
    negative::Int
    plain_failures::Int
    failures::Int
end

RelayAnchoredGapDecoder(dem::Py; kw...) = RelayAnchoredGapDecoder(read_dem(dem); kw...)
function RelayAnchoredGapDecoder(dem::Dem; comp::AbstractString = "relay", comp_kw = nothing, kw...)
    h = _Hyper(dem)
    k = NamedTuple(kw)
    plain = _relay(h, h.H, k)
    no = _no(h)
    decs, bposd = Py[], Py[]
    if comp != "bposd"
        # an empty comp_kw keeps the plain parameters (Python: `if comp_kw:`)
        ckw = comp != "relay" ? (;) : comp_kw === nothing ? (num_sets = 20,) : _kwargs(comp_kw)
        decs = [_relay(h, _with_obs(h, j), merge(k, ckw)) for j in 1:no]
    else
        b = merge((max_iter = 100, bp_method = "minimum_sum", ms_scaling_factor = 1.0, osd_method = "osd_cs",
            osd_order = 7), comp_kw === nothing ? (;) : _kwargs(comp_kw))
        bposd = [pymodule(:ldpc).BpOsdDecoder(_scipy_sparse(_with_obs(h, j), :csr); error_channel = pylist(h.p), b...)
                 for j in 1:no]
    end
    return RelayAnchoredGapDecoder(h, plain, decs, bposd, String(comp), 0, 0, 0)
end

function _comp_weights(dec::RelayAnchoredGapDecoder, j::Int, syn::Matrix{UInt8})
    dec.comp != "bposd" && return _relay_weights(dec.h, dec.decs[j], syn)
    Hj = _with_obs(dec.h, j)
    return map(1:size(syn, 1)) do i
        e = _from_numpy_u8(dec.bposd[j].decode(to_numpy_u8(syn[i, :])))
        _syndrome_ok(Hj, e, syn[i, :]) ? _weight(dec.h.llr, e) : Inf
    end
end

function decode(dec::RelayAnchoredGapDecoder, dets::AbstractMatrix; with_gap::Bool = true)
    h = dec.h
    D = convert(Matrix{UInt8}, dets)
    n, no = size(D, 1), _no(h)
    res = _relay_detailed(dec.plain, D)
    ok = [r[1] for r in res]
    pred = UInt8[_obs_parity(h, j, res[s][2]) for s in 1:n, j in 1:no]
    with_gap || return pred, fill(NaN, n), fill("unavailable", n)
    # (E.astype(bool) * llr).sum(axis=1): the zeros take part in numpy's summation order
    w_pred = [npsum(Float64.(r[2] .!= 0) .* h.llr) for r in res]
    gaps = fill(Inf, n)
    for j in 1:no
        w_comp = _comp_weights(dec, j, hcat(D, UInt8(1) .- pred[:, j]))
        g = w_comp .- w_pred
        dec.negative += count(<(0), g)
        gaps = _npmin.(gaps, _npmax.(g, 0.0))
    end
    dec.plain_failures += count(!, ok)
    gaps[.!ok] .= 0.0
    kind = [isinf(g) ? "not_found" : "approx" for g in gaps]
    kind[.!ok] .= "decoder_failed"
    return pred, gaps, kind
end

#### registry

"""
    CONF_KINDS

Meaning of a per-shot confidence: `exact` (exact complementary gap; `Inf` = no correction with
the flipped outcome exists), `approx` (BP-type gap), `not_found` (BP-type decoder found no
complementary correction; value `Inf`), `decoder_failed` (plain decode did not converge; value 0),
`unavailable` (no gap computed or none exists; `NaN`), `noiseless` (p = 0 run; `Inf`).
"""
const CONF_KINDS = ("exact", "approx", "not_found", "decoder_failed", "unavailable", "noiseless")

"""
    DECODERS

Decoder name => constructor taking a stim DEM (`Py`) or a [`Dem`](@ref) and keyword arguments:
`"pymatching"` (needs a decomposed DEM), `"bposd"`, `"relay"` (plain-anchored gap), `"relay_split"`.
"""
const DECODERS = Dict{String,Any}("pymatching" => MatchingGapDecoder, "bposd" => BPOSDGapDecoder,
    "relay" => RelayAnchoredGapDecoder, "relay_split" => RelayGapDecoder)

# Relay-BP carries an RNG that advances with every decode and cannot be reseeded.
const STATEFUL_DECODERS = ("relay", "relay_split")
