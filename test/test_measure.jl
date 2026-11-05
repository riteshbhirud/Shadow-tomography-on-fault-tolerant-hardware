# measure() and run_oracle reproduce the stored results bit for bit (test/reference/measure and the
# snapshot runs in test/reference/decode/oracle).

@testmodule MeasureRef begin
    using FTShadows, JSON, SHA
    using FTShadows: pyjson
    const DIR = joinpath(pkgdir(FTShadows), "test", "reference", "measure")

    bits(rows, width) = isempty(rows) ? zeros(UInt8, 0, width) :
                        permutedims(reduce(hcat, [UInt8.(collect(r) .== '1') for r in rows]))
    same(a::Float64, b::Float64) = (isnan(a) && isnan(b)) || a === b
    records_text(res) = join(pyjson(r; separators = (",", ":")) * "\n" for r in records(res))

    """Field-by-field comparison of a result with a stored one."""
    function compare(res::OracleResult, d)
        nt = size(res.outcomes, 2)
        conf = parse.(Float64, d["confidence"])
        return (; outcomes = res.outcomes == bits(d["outcomes"], nt), raw = res.raw == bits(d["raw"], nt),
            confidence = length(conf) == length(res.confidence) && all(same.(conf, res.confidence)),
            kind = res.confidence_kind == String.(d["kind"]), cost = pyjson(res.cost) == pyjson(d["cost"]),
            meta = pyjson(res.meta) == pyjson(d["meta"]),
            setting = res.setting == String.(get(d, "setting_labels", get(d, "setting", nothing))),
            records = bytes2hex(sha256(records_text(res))) == d["records_sha256"])
    end

    case_names() = sort([splitext(f)[1] for f in readdir(DIR)])

    function run_case(name)
        d = JSON.parsefile(joinpath(DIR, name * ".json"))
        state = [(String(l), Int(s)) for (l, s) in d["state"]]
        params = Dict(Symbol(k) => v for (k, v) in d["params"])
        res = measure(d["code"], d["setting"], state, d["p"], d["shots"]; seed = d["seed"],
            with_gap = get(d, "with_gap", true), params...)
        return compare(res, d)
    end
end

@testitem "measure reproduces stored results: surface code" tags = [:oracle, :circuits] setup = [MeasureRef] begin
    names = filter(startswith("surface"), MeasureRef.case_names())
    @test length(names) == 8
    for name in names
        r = MeasureRef.run_case(name)
        all(r) || @error "measure differs" name r
        @test all(r)
    end
end

@testitem "measure reproduces stored results: BB72" tags = [:oracle, :circuits, :slow] setup = [MeasureRef] begin
    names = filter(startswith("bb72"), MeasureRef.case_names())
    @test length(names) == 9
    for name in names
        r = MeasureRef.run_case(name)
        all(r) || @error "measure differs" name r
        @test all(r)
    end
end

@testitem "measure reproduces direct builders bit for bit" tags = [:oracle, :circuits, :slow] begin
    using FTShadows.Circuits: surface_local_experiment, local_experiment
    # the plug-ins reproduce the direct builders bit for bit
    st = [("XX", 0), ("ZZ", 0)]
    a = measure("surface", "local:XY", st, 1e-3, 3000; seed = 11, d = 3, memory_rounds = 3)
    exp, lab = surface_local_experiment(3, "XY", st; memory_rounds = 3)
    b = run_oracle(exp, 1e-3, 3000; seed = 11, setting_labels = lab)
    @test a.outcomes == b.outcomes && isequal(a.confidence, b.confidence)
    a = measure("bb72", "local:YZ", st, 1e-3, 50; seed = 3, with_gap = false, rounds = 2)
    exp, lab, _ = local_experiment("YZ", st, 2)
    b = run_oracle(exp, 1e-3, 50; decoder = "relay", seed = 3, setting_labels = lab, with_gap = false)
    @test a.outcomes == b.outcomes
end

@testitem "run_oracle on the snapshot circuits reproduces the stored runs" tags = [:oracle, :circuits] setup = [MeasureRef, OracleRef] begin
    using JSON
    using FTShadows.Circuits: single_patch_readout, multi_patch_surgery, setting_experiment, bb_gauging,
                              bb72_readout_experiment
    st12 = [(String(l), Int(s)) for (l, s) in
            JSON.parsefile(joinpath(pkgdir(FTShadows), "test", "reference", "circuits", "states.json"))["stab:12:4242"]]
    bell = [("ZZ", 0), ("XX", 0)]
    # the snapshot configurations
    configs = [
        "surface_single_X_d3_r2" => () -> single_patch_readout(3, 'X', [("X", 0)], 2),
        "surface_single_Z_d3_r2" => () -> single_patch_readout(3, 'Z', [("Z", 1)], 2),
        "surface_ZZ_2p_d3_r2_conj" => () -> multi_patch_surgery(3, 2, "ZZ", bell, 2),
        "surface_ZZ_2p_d3_r3_conj" => () -> multi_patch_surgery(3, 2, "ZZ", bell, 3),
        "surface_XX_2p_d3_r3_conj" => () -> multi_patch_surgery(3, 2, "XX", bell, 3),
        "surface_ZZ_2p_d3_r2_same" => () -> multi_patch_surgery(3, 2, "ZZ", bell, 2, "same"),
        "surface_ZZ_3p_d3_r2_conj" => () -> multi_patch_surgery(3, 3, "ZZ", [("XXX", 0), ("ZZI", 0), ("IZZ", 0)], 2),
        "surface_XX_2p_d3_r1_conj_nan_gap" => () -> multi_patch_surgery(3, 2, "XX", bell, 1),
        "surface_local_XZ" => () -> setting_experiment("local:XZ", 3, 2, bell, 0)[1],
        "bb72_gauge_X_r1" => () -> bb_gauging('X', 1; state_kind = "+1"),
        "bb72_gauge_Z_r1" => () -> bb_gauging('Z', 1; state_kind = "-1"),
        "bb72_gauge_Y_r1" => () -> bb_gauging('Y', 1; state_labels = st12),
        "bb72_readout_Z" => () -> bb72_readout_experiment('Z', st12)[1],
        "bb72_readout_X" => () -> bb72_readout_experiment('X', st12)[1],
    ]
    snaps = JSON.parsefile(OracleRef.SNAPSHOTS)
    @test sort(first.(configs)) == sort(collect(keys(snaps)))
    for (name, make) in configs
        d = JSON.parsefile(joinpath(pkgdir(FTShadows), "test", "reference", "decode", "oracle", "snap_$name.json"))
        res = run_oracle(make(), d["p"], d["shots"]; decoder = d["decoder"], seed = d["seed"])
        # the hashed snapshot digest, and the full result
        r = merge(MeasureRef.compare(res, d), (; snapshot = OracleRef.digest(res) == Dict(snaps[name])))
        all(r) || @error "snapshot run differs" name r
        @test all(r)
    end
end
