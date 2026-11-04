# Decoders, seeds and statistics against reference values (test/reference/decode/, exported from the
# original implementation).

@testmodule DecodeRef begin
    using FTShadows, JSON, PythonCall
    using FTShadows: pymodule

    const DIR = joinpath(dirname(@__DIR__), "test", "reference", "decode")

    """Rows of '0'/'1' strings as a shots × width UInt8 matrix."""
    function bits(rows, width = 0)
        isempty(rows) && return zeros(UInt8, 0, width)
        return permutedims(reduce(hcat, [UInt8.(collect(r) .== '1') for r in rows]))
    end

    floats(v) = parse.(Float64, v)

    # bitwise equality, NaN equal to NaN (Python's repr drops NaN payloads)
    same(a::Float64, b::Float64) = (isnan(a) && isnan(b)) || a === b
    same(a::AbstractVector{Float64}, b::AbstractVector{Float64}) = length(a) == length(b) && all(same.(a, b))

    case(name) = JSON.parsefile(joinpath(DIR, "decoders", name * ".json"))
    case_names() = sort([splitext(f)[1] for f in readdir(joinpath(DIR, "decoders"))])

    """Rebuild a case's decoder and syndromes as the reference did, decode, return everything to compare."""
    function run_case(d)
        stim = pymodule(:stim)
        if haskey(d, "dem_text")
            dem = stim.DetectorErrorModel(d["dem_text"])
            nd = pyconvert(Int, dem.num_detectors)
            rows = UInt8[(i >> (nd - k)) & 1 for i in 0:(2^nd - 1), k in 1:nd]  # itertools.product order
            dets = reduce(vcat, fill(rows, d["repeat"]))
            obs = zeros(UInt8, size(dets, 1), pyconvert(Int, dem.num_observables))
        else
            circ = JSON.parsefile(joinpath(DIR, "circuits", d["circuit"] * ".json"))
            c = stim.Circuit(circ["eigen"]["circuit"])
            dem = c.detector_error_model(; decompose_errors = d["decoder"] == "pymatching")
            s = c.compile_detector_sampler(; seed = d["seed"]).sample(d["shots"]; separate_observables = true)
            dets = UInt8.(pyconvert(Matrix{Bool}, s[0]))
            obs = UInt8.(pyconvert(Matrix{Bool}, s[1]))
            if get(d, "single", false)
                nd = size(dets, 2)
                dets = vcat(dets, zeros(UInt8, 1, nd), UInt8.([i == j for i in 1:nd, j in 1:nd]))
            end
        end
        kw = haskey(d, "kw") ? FTShadows._kwargs(d["kw"]) : (;)
        dec = DECODERS[d["decoder"]](dem; kw...)
        pred, gap, kind = decode(dec, dets; with_gap = get(d, "with_gap", true))
        return (; dec, obs, pred, gap, kind)
    end

    """Compare one case; returns a NamedTuple of booleans."""
    function check_case(name)
        d = case(name)
        r = run_case(d)
        no = size(r.pred, 2)
        ok_counters = all(begin
            v = getfield(r.dec, Symbol(k))
            k == "max_weight_mismatch" ? same(Float64(v), parse(Float64, ref)) : v == ref
        end for (k, ref) in d["counters"])
        extra = true
        if haskey(d, "gap_exact")
            extra = r.dec.gap_exact == Bool.(d["gap_exact"]) && same(r.dec.undetectable_w, floats(d["undetectable_w"]))
        end
        return (; samples = r.obs == bits(d["obs"], size(r.obs, 2)), pred = r.pred == bits(d["pred"], no),
            gap = same(r.gap, floats(d["gap"])), kind = r.kind == String.(d["kind"]), counters = ok_counters,
            extra, shots = length(r.gap))
    end
end

@testitem "numpy-compatible numerics" tags = [:decode] begin
    using FTShadows: npsum, llr, pymodule
    using PythonCall, Random
    np = pymodule(:numpy)
    rng = Xoshiro(3)
    for n in [0, 1, 7, 8, 9, 127, 128, 129, 255, 300, 1000, 5000]
        a = randn(rng, n) .* exp.(3 .* rand(rng, n))
        @test npsum(a) === pyconvert(Float64, np.asarray(a).sum())
    end
    p = rand(rng, 20000) .* 0.02
    @test llr.(p) == pyconvert(Vector{Float64}, np.log((1 .- np.asarray(p)) ./ np.asarray(p)))
end

@testitem "DEM matrices and graphlike edges match the reference layout" tags = [:decode, :pyref] setup = [DecodeRef] begin
    using FTShadows: read_dem, dem_to_matrices, graphlike_edges, pymodule
    using PythonCall, JSON
    ft = pymodule(Symbol("ftoracle.decode"))
    for (name, decompose) in [("bb_local_XZ", false), ("s_XX_r2_conj", true), ("s_Y_d3", true)]
        circ = JSON.parsefile(joinpath(DecodeRef.DIR, "circuits", name * ".json"))
        dem = pymodule(:stim).Circuit(circ["eigen"]["circuit"]).detector_error_model(; decompose_errors = decompose)
        H, L, p = dem_to_matrices(read_dem(dem))
        r = ft.dem_to_matrices(dem)
        @test Matrix(H) == pyconvert(Matrix{UInt8}, r[0].toarray())
        @test Matrix(L) == pyconvert(Matrix{UInt8}, r[1].toarray())
        @test p == pyconvert(Vector{Float64}, r[2])
        if decompose
            e = graphlike_edges(read_dem(dem))
            pe = ft.graphlike_edges(dem)
            ref = [(pyconvert(Int, k[0]), pyconvert(Any, k[1]) === nothing ? -1 : pyconvert(Int, k[1]),
                    UInt64(pyconvert(Int, k[2])), pyconvert(Float64, v)) for (k, v) in pe.items()]
            @test [(k..., q) for (k, q) in zip(e.keys, e.p)] == ref
        end
    end
end

@testitem "matching gap decoder: shot-by-shot identity (surface + synthetic)" tags = [:decode] setup = [DecodeRef] begin
    names = filter(n -> endswith(n, "_pm") || endswith(n, "_pm_nogap") || endswith(n, "_pymatching"), DecodeRef.case_names())
    @test length(names) == 14
    for name in names
        r = DecodeRef.check_case(name)
        all(r[k] for k in keys(r) if k != :shots) || @error "decoder case differs" name r
        @test r.samples && r.pred && r.gap && r.kind && r.counters && r.extra
    end
end

@testitem "BP-type decoders: shot-by-shot identity (synthetic DEMs)" tags = [:decode] setup = [DecodeRef] begin
    names = filter(n -> startswith(n, "syn_hyp"), DecodeRef.case_names())
    @test length(names) == 5
    for name in names
        r = DecodeRef.check_case(name)
        all(r[k] for k in keys(r) if k != :shots) || @error "decoder case differs" name r
        @test r.pred && r.gap && r.kind && r.counters
    end
end

@testitem "BP-type decoders: shot-by-shot identity on BB72" tags = [:decode, :slow] setup = [DecodeRef] begin
    names = filter(n -> startswith(n, "bb_"), DecodeRef.case_names())
    @test length(names) == 18
    for name in names
        r = DecodeRef.check_case(name)
        all(r[k] for k in keys(r) if k != :shots) || @error "decoder case differs" name r
        @test r.samples && r.pred && r.gap && r.kind && r.counters
    end
end

@testitem "BP-OSD fallback when both hypotheses fail" tags = [:decode] begin
    # the original code raises here (int() of a 1-element array under numpy 2.5), so this path is checked on its own:
    # odd syndromes are outside the column space, so each observable falls back to the plain decode.
    # The weak mechanisms give the check matrix more columns than rows plus the OSD order (7); with
    # fewer, ldpc's OSD-CS setup reads out of bounds (it crashed on Linux).
    obs = ["", " L0", " L1", " L0 L1"]
    weak = ["error(0.0001) $d$o" for (d, o0) in (("D0 D1", " L0"), ("D1 D2", " L1"), ("D0 D2", "")) for o in obs if o != o0]
    dec = BPOSDGapDecoder(FTShadows.read_dem(join(["error(0.1) D0 D1 L0", "error(0.1) D1 D2 L1", "error(0.05) D0 D2", weak...], "\n")))
    dets = UInt8[(i >> (3 - k)) & 1 for i in 0:7, k in 1:3]
    pred, gap, kind = decode(dec, dets)
    odd = isodd.(vec(sum(dets; dims = 2)))
    @test all(gap[odd] .== 0) && all(gap[.!odd] .> 0)
    @test all(kind .== "approx")
    @test pred[.!odd, :] == UInt8[0 0; 0 1; 0 0; 1 0]  # lightest corrections of 000, 011, 101, 110
end

@testitem "seed_for and wilson match the reference values" tags = [:decode] setup = [DecodeRef] begin
    using JSON
    ref = JSON.parsefile(joinpath(DecodeRef.DIR, "seeds_stats.json"))
    pyvalue(x) = x isa AbstractVector ? Any[pyvalue(y) for y in x] : x
    for s in ref["seeds"]
        @test seed_for(pyvalue(s["parts"])...) == s["seed"]
    end
    for w in ref["wilson"]
        lo, hi = wilson(w["k"], w["n"]; conf = w["conf"])
        @test DecodeRef.same(lo, parse(Float64, w["lo"])) && DecodeRef.same(hi, parse(Float64, w["hi"]))
    end
    @test holm([0.01, 0.04, 0.03]; alpha = 0.05) == [true, false, false]
    @test holm([0.001, 0.02, 0.03]; alpha = 0.1) == [true, true, true]
    z, p = two_proportion_test(10, 1000, 10, 1000)
    @test z == 0 && p == 1
end
