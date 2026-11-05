# Bravyi et al.'s BB72 memory logical error rate (arXiv:2308.07915) on their circuit, with their BP-OSD
# settings and with Relay-BP, and the sensitivity of the failure rate to the BP variant. ldpc and Relay-BP
# are called through PythonCall.

const BRAVYI_NC = 6
const BRAVYI_CHUNK = 2000
const BRAVYI_FIT = (11.09, 365.6, -16088)

# pow from the C library and exp from numpy, so the fit values match the stored ones on every platform
# (Julia's `^` and numpy's SIMD exp can differ in the last bit).
_cpow(x::Float64, y::Float64) = @ccall pow(x::Float64, y::Float64)::Float64

"""The paper's fit pL(p) = p^3 exp(c0 + c1 p + c2 p^2)."""
bravyi_fit_pl(p::Float64) = _cpow(p, 3.0) * pyconvert(Float64, pymodule(:numpy).exp(BRAVYI_FIT[1] + BRAVYI_FIT[2] * p + BRAVYI_FIT[3] * p * p))

"""32-bit seed of chunk `chunk` of a configuration tag (SHA-256 of "tag|chunk", first 4 bytes little-endian)."""
bravyi_seed(tag::AbstractString, chunk::Integer) = Int(ltoh(reinterpret(UInt32, sha256("$tag|$chunk")[1:4])[1]))

"""`(p, basis, decoder)` as written in tags and chunk file names (floats as Python's str)."""
bravyi_key(k) = "$(pyfloat(Float64(k[1])))_$(k[2])_$(k[3])"
bravyi_tag(k, scaling) = "bravyi|$(pyfloat(Float64(k[1])))|$(k[2])|$(k[3])|$(pyfloat(Float64(scaling)))"

"""Check matrix, observable matrix and priors of Bravyi's memory circuit (N_c noisy cycles) in `basis`."""
function bravyi_memory(basis::AbstractString, p::Float64)
    c = Circuits.pycircuit(Circuits.memory_circuit("bb72", only(basis), BRAVYI_NC, p))
    return c, dem_to_matrices(read_dem(c.detector_error_model()))...
end

const _BRAVYI_CACHE = Dict{Any,Any}()

"""
    bravyi_work(p, basis, decoder, chunk, scaling) -> fails

One chunk of `BRAVYI_CHUNK` shots: shots where any of the 12 logicals is predicted wrong. BP-OSD is
min-sum with fixed scaling `scaling`, OSD-CS order 7, max_iter 10^4, decoding shot by shot; Relay-BP
decodes the chunk as one batch (its results depend on the batch, relay-bp issue #47).
"""
function bravyi_work(p::Float64, basis::AbstractString, decoder::AbstractString, chunk::Integer, scaling::Float64)
    c, H, L, pr, dec = get!(_BRAVYI_CACHE, (p, basis, decoder, scaling)) do
        c, H, L, pr = bravyi_memory(basis, p)
        dec = if decoder == "bposd"
            pymodule(:ldpc).BpOsdDecoder(_scipy_sparse(H, :csr); error_channel = pylist(pr), max_iter = 10000,
                bp_method = "minimum_sum", ms_scaling_factor = scaling, osd_method = "osd_cs", osd_order = 7)
        else
            pymodule(:relay_bp).RelayDecoderF64(_scipy_sparse(H, :csr); error_priors = _np().array(pr; dtype = _np().float64),
                gamma0 = 0.125, pre_iter = 80, num_sets = 100, set_max_iter = 60,
                gamma_dist_interval = pytuple((-0.24, 0.66)), stop_nconv = 5)
        end
        (c, H, L, pr, dec)
    end
    seed = bravyi_seed(bravyi_tag((p, basis, decoder), scaling), chunk)
    res = c.compile_detector_sampler(; seed).sample(BRAVYI_CHUNK; separate_observables = true)
    det = pyconvert(Matrix{UInt8}, _np().asarray(res[0]; dtype = _np().uint8))
    obs = pyconvert(Matrix{Bool}, res[1])
    E = if decoder == "bposd"
        reduce(vcat, [permutedims(_from_numpy_u8(dec.decode(to_numpy_u8(det[s, :])))) for s in 1:size(det, 1)])
    else
        reshape(_from_numpy_u8(dec.decode_batch(to_numpy_u8(det))), size(det, 1), :)
    end
    pred = (Int.(L) * Int.(transpose(E))) .% 2 .== 1
    return count(s -> any(pred[:, s] .!= obs[s, :]), 1:size(det, 1))
end

"""The 12 configurations (p, basis, decoder), BP-OSD first."""
bravyi_configs() = [(p, string(b), dec) for dec in ("bposd", "relay") for p in (1e-3, 2e-3, 3e-3) for b in "ZX"]

"""(fails, shots) summed over the finished chunk files of configuration `k`."""
function bravyi_counts(chunk_dir::AbstractString, k)
    e = n = 0
    pre = bravyi_key(k) * "_c"
    for f in readdir(chunk_dir)
        startswith(f, pre) && endswith(f, ".json") || continue
        d = JSON.parsefile(joinpath(chunk_dir, f))
        e += d["fails"]
        n += d["shots"]
    end
    return e, n
end

_claimed(chunk_dir, k) = count(f -> startswith(f, bravyi_key(k) * "_c") && endswith(f, ".claim"), readdir(chunk_dir))

# Create `path` only if it does not exist yet; false if it does.
function _claim(path::AbstractString)
    fl = Base.Filesystem.JL_O_WRONLY | Base.Filesystem.JL_O_CREAT | Base.Filesystem.JL_O_EXCL
    try
        close(Base.Filesystem.open(path, fl, 0o644))
        return true
    catch e
        e isa Base.IOError && e.code == Base.UV_EEXIST && return false
        rethrow()
    end
end

"""
    bravyi_run_worker(cfgs, target, mx, chunk_dir, scaling)

Loop: pick the first configuration still below `target` errors and the shot cap `mx`, claim its lowest
unclaimed chunk index by creating `<cfg>_c<i>.claim` exclusively, run it, write `<cfg>_c<i>.json`.
Seeds depend only on (configuration, chunk index), so results do not depend on scheduling.
"""
function bravyi_run_worker(cfgs, target, mx, chunk_dir::AbstractString, scaling::Float64)
    while true
        did = false
        for k in cfgs
            e, n = bravyi_counts(chunk_dir, k)
            (e >= target || n >= mx || _claimed(chunk_dir, k) * BRAVYI_CHUNK >= mx) && continue
            c = 0
            while !_claim(joinpath(chunk_dir, "$(bravyi_key(k))_c$c.claim"))
                c += 1
            end
            fails = bravyi_work(k..., c, scaling)
            tmp = joinpath(chunk_dir, "$(bravyi_key(k))_c$c.tmp")
            d = OrderedDict("fails" => fails, "shots" => BRAVYI_CHUNK, "seed" => bravyi_seed(bravyi_tag(k, scaling), c))
            write(tmp, pyjson(d))
            mv(tmp, joinpath(chunk_dir, "$(bravyi_key(k))_c$c.json"); force = true)
            did = true
            break
        end
        did || return
    end
end

"""
    bravyi_rows(chunk_dir, cfgs) -> rows

Per p and decoder: failures of the Z and X memories, PL = 1 - (1 - PZ)(1 - PX), pL per cycle
1 - (1 - PL)^(1/N_c), its interval from the two Wilson intervals' endpoints, and the paper's fit.
"""
function bravyi_rows(chunk_dir::AbstractString, cfgs)
    st = Dict(k => bravyi_counts(chunk_dir, k) for k in cfgs)
    rows = OrderedDict{String,Any}[]
    inv = 1 / BRAVYI_NC
    for p in (1e-3, 2e-3, 3e-3), dec in ("bposd", "relay")
        haskey(st, (p, "Z", dec)) || continue
        (ez, nz), (ex, nx) = st[(p, "Z", dec)], st[(p, "X", dec)]
        PZ, PX = ez / nz, ex / nx
        PL = 1 - (1 - PZ) * (1 - PX)
        pl = 1 - _cpow(1 - PL, inv)
        # CI: combine the two Wilson intervals conservatively (endpoints)
        zlo, zhi = wilson(ez, nz)
        xlo, xhi = wilson(ex, nx)
        lo = 1 - _cpow(1 - (1 - (1 - zlo) * (1 - xlo)), inv)
        hi = 1 - _cpow(1 - (1 - (1 - zhi) * (1 - xhi)), inv)
        push!(rows, OrderedDict{String,Any}("p" => p, "decoder" => dec, "cycles" => BRAVYI_NC, "shots_Z" => nz,
            "fails_Z" => ez, "shots_X" => nx, "fails_X" => ex, "PL" => PL, "pL_per_cycle" => pl, "pL_ci" => [lo, hi],
            "paper_fit" => bravyi_fit_pl(p), "paper_table1" => p == 1e-3 ? 7e-5 : nothing,
            "ratio_to_fit" => pl / bravyi_fit_pl(p)))
    end
    return rows
end

#### Decoder sensitivity

const BRAVYI_SENS_P, BRAVYI_SENS_SHOTS = 3e-3, 3000
const BRAVYI_SENS_VARIANTS = [("ms", 1.0), ("ms", 0.9), ("ms", 0.75), ("ms", 0.625), ("ms", 0.5), ("ps", nothing)]

"""
    bravyi_sensitivity(v) -> OrderedDict

Failures of one BP-OSD variant `v = (method, scaling)` ("ms" min-sum with scaling, "ps" product-sum) on
the same 3000 Z-memory shots at p = 3e-3, with the per-shot failure vector packed by `packbits`.
"""
function bravyi_sensitivity(v)
    method, sf = v
    c, H, L, pr = bravyi_memory("Z", BRAVYI_SENS_P)
    ldpc = pymodule(:ldpc)
    kw = (error_channel = pylist(pr), max_iter = 10000, osd_method = "osd_cs", osd_order = 7)
    dec = method == "ms" ? ldpc.BpOsdDecoder(_scipy_sparse(H, :csr); bp_method = "minimum_sum", ms_scaling_factor = sf, kw...) :
          ldpc.BpOsdDecoder(_scipy_sparse(H, :csr); bp_method = "product_sum", kw...)
    res = c.compile_detector_sampler(; seed = seed_for("bravyi_sens", BRAVYI_SENS_P)).sample(BRAVYI_SENS_SHOTS;
                                                                                             separate_observables = true)
    det = pyconvert(Matrix{UInt8}, _np().asarray(res[0]; dtype = _np().uint8))
    obs = pyconvert(Matrix{Bool}, res[1])
    Li = Int.(L)
    fails = Bool[]
    for s in 1:size(det, 1)
        e = _from_numpy_u8(dec.decode(to_numpy_u8(det[s, :])))
        push!(fails, any((Li * Int.(e)) .% 2 .!= obs[s, :]))
    end
    return OrderedDict{String,Any}("bp" => method, "scaling" => sf, "fails" => count(fails),
                                   "shots" => BRAVYI_SENS_SHOTS, "fail_vec" => packbits(fails))
end

"""Bool vector packed 8 per byte, most significant bit first, zero-padded (as `np.packbits`)."""
packbits(b::AbstractVector{Bool}) = [sum(b[i+j] << (7 - j) for j in 0:7 if i + j <= length(b); init = 0)
                                     for i in 1:8:length(b)]
