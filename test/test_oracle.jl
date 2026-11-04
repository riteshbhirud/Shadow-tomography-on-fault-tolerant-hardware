# run_oracle reproduces the reference runs (outcomes, raw, confidence, kinds, cost, meta), the hashed
# snapshots in test/reference/snapshots/snapshots.json, and the records byte for byte.

@testmodule OracleRef begin
    using FTShadows, JSON, SHA
    using FTShadows: pyjson

    const DIR = joinpath(dirname(@__DIR__), "test", "reference", "decode")
    const SNAPSHOTS = joinpath(DIR, "..", "snapshots", "snapshots.json")  # hashed snapshots of the original implementation

    bits(rows, width) = isempty(rows) ? zeros(UInt8, 0, width) :
                        permutedims(reduce(hcat, [UInt8.(collect(r) .== '1') for r in rows]))
    same(a::Float64, b::Float64) = (isnan(a) && isnan(b)) || a === b

    case_names() = sort([splitext(f)[1] for f in readdir(joinpath(DIR, "oracle"))])

    """The snapshot digest: sha256 over outcomes, raw and the rounded confidences."""
    function digest(res::OracleResult)
        c = [isnan(x) ? -1.0 : isinf(x) ? 1e9 : x for x in res.confidence]
        c = [round(x * 10000.0) / 10000.0 for x in c]  # np.round(c, 4)
        io = IOBuffer()
        write(io, vec(permutedims(res.outcomes)))
        write(io, vec(permutedims(res.raw)))
        write(io, htol.(c))
        h = bytes2hex(sha256(take!(io)))[1:20]
        return Dict("sha256" => h, "detectors" => res.meta["detectors"], "flips" => Int(sum(res.outcomes)),
            "raw_ones" => Int(sum(res.raw)))
    end

    records_text(res) = join(pyjson(r; separators = (",", ":")) * "\n" for r in records(res))

    function run_case(name)
        d = JSON.parsefile(joinpath(DIR, "oracle", name * ".json"))
        exp = load_prebuilt(joinpath(DIR, "circuits", d["circuit"] * ".json"))
        labels = get(d, "labels", nothing)
        res = run_oracle(exp, d["p"], d["shots"]; decoder = d["decoder"], seed = d["seed"],
            setting_labels = labels === nothing ? nothing : String.(labels), decoder_kwargs = get(d, "kw", nothing),
            batch = get(d, "batch", 20000), with_gap = get(d, "with_gap", true))
        return d, res
    end

    """Field-by-field comparison with the reference result; returns a NamedTuple of booleans."""
    function check_case(name)
        d, res = run_case(name)
        nt = size(res.outcomes, 2)
        conf = parse.(Float64, d["confidence"])
        text = records_text(res)
        out = (; outcomes = res.outcomes == bits(d["outcomes"], nt), raw = res.raw == bits(d["raw"], nt),
            confidence = length(conf) == length(res.confidence) && all(same.(conf, res.confidence)),
            kind = res.confidence_kind == String.(d["kind"]), cost = pyjson(res.cost) == pyjson(d["cost"]),
            meta = pyjson(res.meta) == pyjson(d["meta"]), setting = res.setting == String.(d["setting"]),
            digest = digest(res) == Dict(d["digest"]),
            records = bytes2hex(sha256(text)) == d["records_sha256"] &&
                      split(text, '\n')[1:length(d["records_head"])] == d["records_head"])
        if startswith(name, "snap_")
            snap = JSON.parsefile(SNAPSHOTS)[name[6:end]]
            out = merge(out, (; snapshot = digest(res) == Dict(snap)))
        end
        return out, size(res.outcomes, 1)
    end
end

@testitem "oracle.run parity: surface snapshots and extra cases" tags = [:oracle] setup = [OracleRef] begin
    names = filter(n -> !occursin("bb", n), OracleRef.case_names())
    @test length(names) == 12
    for name in names
        r, _ = OracleRef.check_case(name)
        all(r) || @error "oracle case differs" name r
        @test all(r)
    end
end

@testitem "oracle.run parity: BB72 snapshots and Relay/BP-OSD runs" tags = [:oracle, :slow] setup = [OracleRef] begin
    names = filter(n -> occursin("bb", n), OracleRef.case_names())
    @test length(names) == 8
    for name in names
        r, _ = OracleRef.check_case(name)
        all(r) || @error "oracle case differs" name r
        @test all(r)
    end
end

@testitem "records written by write_records_raw equal the reference bytes" tags = [:oracle] setup = [OracleRef] begin
    using SHA
    using FTShadows: write_records_raw
    d, res = OracleRef.run_case("snap_surface_XX_2p_d3_r1_conj_nan_gap")  # NaN confidences
    mktempdir() do tmp
        f = joinpath(tmp, "r.jsonl")
        write_records_raw(f, records(res))
        @test bytes2hex(sha256(read(f))) == d["records_sha256"]
    end
end

@testitem "oracle conventions" tags = [:oracle] setup = [OracleRef] begin
    exp = load_prebuilt(joinpath(OracleRef.DIR, "circuits", "s_ZZ_r3_conj.json"))
    a = run_oracle(exp, 5e-3, 400; seed = 9)
    b = run_oracle(exp, 5e-3, 400; seed = 9)
    @test a.outcomes == b.outcomes && a.raw == b.raw && isequal(a.confidence, b.confidence)
    @test run_oracle(exp, 5e-3, 400; seed = 10).raw != a.raw  # (stim's stream depends on the batch split too)
    c = run_oracle(exp, 5e-3, 400; seed = 9, with_gap = false)
    @test c.outcomes == a.outcomes && all(isnan, c.confidence) && all(==("unavailable"), c.confidence_kind)
    @test all(a.confidence .>= 0) && all(==("exact"), a.confidence_kind)
    @test_throws ErrorException run_oracle(exp, 1e-3, 10)  # no build at this p
    r1 = load_prebuilt(joinpath(OracleRef.DIR, "circuits", "s_XX_r1_conj.json"))
    @test all(==("unavailable"), run_oracle(r1, 5e-3, 200; seed = 1).confidence_kind)
end

# On fresh seeds, q agrees with the stored result tables within sampling error (pooled
# two-proportion z-tests, Holm over all rows). About 45 CPU-minutes, run on worker processes since
# PythonCall is not thread-safe; set FTSHADOWS_STAT=1 to run it (FTSHADOWS_STAT_WORKERS, default 6).
@testitem "oracle q agrees with the reference tables (statistical)" tags = [:oracle, :slow, :statistical] setup = [OracleRef] begin
    using Distributed, JSON, Printf
    if get(ENV, "FTSHADOWS_STAT", "") != "1"
        @info "statistical gate skipped (set FTSHADOWS_STAT=1)"
    else
        ref = JSON.parsefile(joinpath(OracleRef.DIR, "stat_rows.json"))
        circ(name) = joinpath(OracleRef.DIR, "circuits", name * ".json")
        chunk_of(row) = row["decoder"] == "pymatching" ? 50_000 : 2_500
        # shot counts: as many as the stored row for surface; 15k (r = 2) and 20k (r = 3) for BB72
        julia_shots(row) = row["decoder"] == "pymatching" ? row["shots"] : (row["row"]["r"] == "2" ? 15_000 : 20_000)
        jobs = Tuple[]
        for (i, row) in enumerate(ref["qr"]), c in 0:(cld(julia_shots(row), chunk_of(row)) - 1)
            n = min(chunk_of(row), julia_shots(row) - c * chunk_of(row))
            push!(jobs, (:qr, i, circ(row["circuit"]), row["decoder"], n, seed_for("julia_stat_gate", row["circuit"], c), false))
        end
        for (i, row) in enumerate(ref["postselection"]), c in 0:11
            push!(jobs, (:post, i, circ(row["circuit"]), row["decoder"], 100, seed_for("julia_stat_post", row["circuit"], c), true))
        end
        sort!(jobs; by = j -> (j[1] == :post ? 0 : j[4] == "relay" ? 1 : 2))  # slow jobs first
        wp = addprocs(parse(Int, get(ENV, "FTSHADOWS_STAT_WORKERS", "6")); exeflags = "--project=$(Base.active_project())")
        t0 = time()
        out = try
            @everywhere begin
                using FTShadows
                function stat_chunk(job)
                    kind, i, path, decoder, shots, seed, with_gap = job
                    res = run_oracle(load_prebuilt(path), 1e-3, shots; decoder, seed, with_gap)
                    return (kind, i, res.outcomes[:, 1] .!= 0, res.confidence)
                end
            end
            pmap(invokelatest(getglobal, Main, :stat_chunk), WorkerPool(wp), jobs)
        finally
            rmprocs(wp)
        end
        @info @sprintf("statistical gate: %d chunks in %.0f s wall", length(jobs), time() - t0)
        pvals, lines = Float64[], String[]
        for (i, row) in enumerate(ref["qr"])
            e = reduce(vcat, [o[3] for o in out if o[1] == :qr && o[2] == i])
            k, n = count(e), length(e)
            z, p = two_proportion_test(k, n, row["errors"], row["shots"])
            push!(pvals, p)
            push!(lines, @sprintf("%-14s %-16s q_julia = %d/%d = %.3g   q_python = %d/%d = %.3g   z = %+.2f  p = %.3f",
                row["table"], row["circuit"], k, n, k / n, row["errors"], row["shots"], row["errors"] / row["shots"], z, p))
        end
        for (i, row) in enumerate(ref["postselection"])
            parts = [o for o in out if o[1] == :post && o[2] == i]
            e = reduce(vcat, [o[3] for o in parts])
            g = reduce(vcat, [o[4] for o in parts])
            order = sortperm([isnan(x) ? -0.0 : -x for x in g]; lt = <, alg = Base.Sort.DEFAULT_STABLE)  # np.argsort(kind="stable")
            kept = e[order][1:floor(Int, row["keep"] * length(e))]
            for (k, n, k2, n2, what) in ((count(e), length(e), row["errors"], row["shots"], "all"),
                                         (count(kept), length(kept), row["errors_kept"], row["kept"], "keep $(row["keep"])"))
                z, p = two_proportion_test(k, n, k2, n2)
                push!(pvals, p)
                push!(lines, @sprintf("%-14s %-16s %-9s q_julia = %d/%d = %.3g   q_python = %d/%d = %.3g   z = %+.2f  p = %.3f",
                    "postselection", row["circuit"], what, k, n, k / n, k2, n2, k2 / n2, z, p))
            end
        end
        @info "statistical gate\n" * join(lines, "\n")
        @test length(pvals) == 16
        @test !any(holm(pvals; alpha = 0.01))
    end
end

# Decoder properties (gap, confidences, PyMatching agreement) on the exported circuits.
@testitem "oracle properties: surface" tags = [:oracle] setup = [OracleRef] begin
    using FTShadows: read_dem, graphlike_edges, llr, detector_error_model, pymodule
    using PythonCall, Statistics
    load(name) = load_prebuilt(joinpath(OracleRef.DIR, "circuits", name * ".json"))
    eigen(name, p) = oracle_circuit(load(name), p; eigen_targets = true).circuit

    # gap decoder consistent with plain matching; errors concentrate at small gaps
    for name in ("s_ZZ_r3_conj", "s_XX_r2_conj", "s_ZZ_3p_r2")
        c = pymodule(:stim).Circuit(eigen(name, 5e-3))
        dec = MatchingGapDecoder(c.detector_error_model(; decompose_errors = true))
        s = c.compile_detector_sampler(; seed = 2).sample(3000; separate_observables = true)
        obs = pyconvert(Matrix{Bool}, s[1])
        pred, gap, _ = decode(dec, pyconvert(Matrix{Bool}, s[0]))
        @test dec.max_weight_mismatch < 1e-6 && dec.pred_mismatches == 0 && all(gap .>= 0)
        wrong = vec(any(pred .!= obs; dims = 2))
        count(wrong) > 20 && @test median(gap[wrong]) < median(gap[.!wrong])
    end

    # audit: parallel boundary mechanisms with different labels (ZZ, r = 2): the complementary weight
    # never exceeds an explicitly achievable one (corrections of at most 2 mechanisms)
    dem = detector_error_model(eigen("snap_surface_ZZ_2p_d3_r2_conj", 2e-3); decompose_errors = true)
    E = graphlike_edges(read_dem(dem))
    edges = [(u, v, om, llr(q)) for ((u, v, om), q) in zip(E.keys, E.p)]
    U = Dict{Set{Int},Vector{Float64}}()
    for i in eachindex(edges), j in [[nothing]; collect(i+1:length(edges))]
        comb = j === nothing ? (edges[i],) : (edges[i], edges[j])
        s, o, w = Set{Int}(), UInt64(0), 0.0
        for (u, v, om, wt) in comb
            symdiff!(s, v == -1 ? Set([u]) : Set([u, v]))
            o ⊻= om
            w += wt
        end
        length(s) == 1 || continue
        d = get!(U, s, [Inf, Inf])
        d[o + 1] = min(d[o + 1], w)
    end
    dec = MatchingGapDecoder(dem)
    nd = pyconvert(Int, dem.num_detectors)
    function excess(s, u0, u1)  # w_best + gap - max(U_0, U_1) for a single-detector syndrome
        det = zeros(UInt8, 1, nd)
        det[1, first(s) + 1] = 1
        wb = pyconvert(Vector{Float64}, dec.base.decode_batch(FTShadows.to_numpy_u8(det); return_weights = true)[1])[1]
        return wb + decode(dec, det)[2][1] - max(u0, u1)
    end
    worst = maximum(excess(s, u0, u1) for (s, (u0, u1)) in U if isfinite(u0) && isfinite(u1))
    @test worst <= 1e-6

    # audit: r = 1 has a detector-free observable flip; no exact gap (NaN), never overconfident
    dem1 = detector_error_model(eigen("s_ZZ_r1_conj", 5e-3); decompose_errors = true)
    pure = [e.p for e in read_dem(dem1).errors if all(isempty(c[1]) for c in e.components)]
    @test !isempty(pure)
    d1 = MatchingGapDecoder(dem1)
    _, gap, _ = decode(d1, zeros(UInt8, 1, pyconvert(Int, dem1.num_detectors)))
    @test isnan(gap[1]) && d1.undetectable_w[1] <= minimum(llr.(pure)) + 1e-9

    # XX surgery at r = 1 does not crash; confidence NaN by policy
    c = run_oracle(load("s_XX_r1_conj"), 5e-3, 500; seed = 1).confidence
    @test all(isnan.(c) .| (c .>= 0))

    # confidence kinds consistent with values; the Y measurement has an exact gap
    @test all(==("exact"), run_oracle(load("s_ZZ_r3_conj"), 5e-3, 300; seed = 1).confidence_kind)
    r0 = run_oracle(load("s_ZZ_r3_conj_p0"), 0.0, 50; seed = 1)
    @test all(==("noiseless"), r0.confidence_kind) && all(isinf, r0.confidence)
    ry = run_oracle(load("s_Y_d3"), 5e-3, 2000; seed = 5)
    @test all(==("exact"), ry.confidence_kind) && all(ry.confidence .>= 0)
    @test haskey(records(ry)[1], "meta") && all(haskey(records(ry)[1]["meta"], k) for k in ("p", "noise", "decoder", "seed"))
end

@testitem "oracle properties: BB72 Relay-BP and BP-OSD" tags = [:oracle, :slow] setup = [OracleRef] begin
    using Statistics
    load(name) = load_prebuilt(joinpath(OracleRef.DIR, "circuits", name * ".json"))

    # kinds of the anchored gap
    rb = run_oracle(load("bb_X_r1"), 2e-3, 80; decoder = "relay", seed = 1)
    k, c = rb.confidence_kind, rb.confidence
    @test issubset(Set(k), Set(["approx", "not_found", "decoder_failed"]))
    @test all(isinf, c[k .== "not_found"]) && all(isfinite, c[k .== "approx"]) && all(c[k .== "decoder_failed"] .== 0)

    # anchored predictions equal plain Relay-BP decoding of the same shots
    a = run_oracle(load("stat_bb_X_r2"), 1e-3, 200; decoder = "relay", seed = 11)
    b = run_oracle(load("stat_bb_X_r2"), 1e-3, 200; decoder = "relay", seed = 11, with_gap = false)
    @test a.outcomes == b.outcomes

    # crippled Relay-BP: shots without a usable decode keep the plain prediction, not outcome 0 (state -1)
    res = run_oracle(load("bb_Z_r2m_p1e3"), 1e-3, 200; decoder = "relay", seed = 2,
        decoder_kwargs = Dict("num_sets" => 1, "set_max_iter" => 5, "pre_iter" => 5))
    zero = res.confidence .== 0
    @test count(zero) > 20 && mean(res.outcomes[zero, 1] .== 1) > 0.7

    # the decoder beats raw readout at r = 3, with non-negative gaps
    r3 = run_oracle(load("stat_bb_X_r3"), 1e-3, 150; decoder = "relay", seed = 3)
    @test mean(r3.outcomes) < mean(r3.raw) && all(r3.confidence .>= 0)

    # BP-OSD with several observables never reports a negative gap (audit)
    bp = run_oracle(load("bb_local_XZ"), 2e-3, 40; decoder = "bposd", seed = 5,
        decoder_kwargs = (osd_order = 2, max_iter = 30))
    @test size(bp.outcomes, 2) > 1 && all(bp.confidence .>= 0)
end
