# The setting-cost table, the Nov-20 expected values and datasets, and the surface
# calibration rows reproduce the stored files. The data files are local (results/, handover/ are not
# tracked), so these items skip when they are absent.

@testmodule HandoverRef begin
    using FTShadows, JSON
    const REPO = pkgdir(FTShadows)
    const NOV20 = joinpath(REPO, "handover", "data", "nov20")
    const CAL = joinpath(REPO, "handover", "calibration_local.csv")
    const RESULTS = joinpath(REPO, "results", "taskA")
    # meta fields that record provenance, not content
    const PROVENANCE = ("created", "runtime_s", "generator", "code_version", "calibration_source")

    have(paths...) = all(ispath, paths) || (@info "skipping: $(join(filter(!ispath, collect(paths)), ", ")) not found"; false)

    function meta_equal(a, b)
        a, b = JSON.parsefile(a), JSON.parsefile(b)
        foreach(k -> (delete!(a, k); delete!(b, k)), PROVENANCE)
        return FTShadows.pyjson(a; indent = 1) == FTShadows.pyjson(b; indent = 1)
    end
end

@testitem "numpy and Python numerics are reproduced" tags = [:handover, :bridge] begin
    using PythonCall, Random
    using FTShadows: np_mean, np_round, py_round
    np = pyimport("numpy")
    rng = Xoshiro(3)
    for _ in 1:500
        a = 1 .- 2 .* rand(rng, rand(rng, 1:200)) .* 10.0 .^ (-rand(rng, 1:5))
        @test np_mean(a) === pyconvert(Float64, np.mean(pylist(a)))
        x = (rand(rng) - 0.5) * 2
        @test np_round(x, 12) === pyconvert(Float64, pybuiltins.round(np.float64(x), 12))
        @test py_round(x, 6) === pyconvert(Float64, pybuiltins.round(x, 6))
    end
end

@testitem "setting-cost table matches results/taskA" tags = [:handover, :slow] setup = [HandoverRef] begin
    using FTShadows: local_setting_costs, write_setting_costs
    R = HandoverRef.RESULTS
    if HandoverRef.have(joinpath(R, "local_setting_costs.csv"), joinpath(R, "bb72_local_q.csv"))
        rows = local_setting_costs(joinpath(R, "surface_local_q.csv"), joinpath(R, "bb72_local_q.csv"))
        @test length(rows) == 284
        mktempdir() do tmp
            write_setting_costs(joinpath(tmp, "c.csv"), rows)
            @test read(joinpath(tmp, "c.csv"), String) == read(joinpath(R, "local_setting_costs.csv"), String)
        end
    end
end

@testitem "Nov-20 exact and expected values equal the handover files" tags = [:handover] setup = [HandoverRef] begin
    using JSON
    using FTShadows: nov20_datasets, nov20_expectations, write_nov20, csv_rows, product_strings, pyjson
    if HandoverRef.have(HandoverRef.NOV20, HandoverRef.CAL)
        cal = csv_rows(HandoverRef.CAL)
        dss = nov20_datasets()
        @test length(dss) == 36
        for ds in dss
            rows, known = nov20_expectations(ds, product_strings("XYZ", ds["m"]), cal)
            meta = JSON.parsefile(joinpath(HandoverRef.NOV20, ds["name"] * ".meta.json"))
            mktempdir() do tmp
                write_nov20(tmp, ds, Dict(), [], rows)
                for ext in (".exact.csv", ".expected.csv")
                    ok = read(joinpath(tmp, ds["name"] * ext)) == read(joinpath(HandoverRef.NOV20, ds["name"] * ext))
                    ok || @error "differs" ds["name"] ext
                    @test ok
                end
            end
            @test pyjson(known; indent = 1) == pyjson(meta["known_answers"]; indent = 1)
        end
    end
end

@testitem "Nov-20 datasets regenerate byte for byte" tags = [:handover, :slow] setup = [HandoverRef] begin
    using FTShadows: nov20_datasets, generate_nov20, csv_rows
    if HandoverRef.have(HandoverRef.NOV20, HandoverRef.CAL)
        cal = csv_rows(HandoverRef.CAL)
        mktempdir() do tmp
            for ds in nov20_datasets()
                ds["name"] in ("surface_d3_bell_fixed_p1e-3_r3", "surface_d3_ghz_random_p0",
                               "surface_d3_bell_random_p1e-3_r3_logicalZ0_pL0.05", "bb72_bell_fixed_p0") || continue
                generate_nov20(ds, tmp; calibration = cal)
                for ext in (".records.jsonl", ".exact.csv", ".expected.csv")
                    @test read(joinpath(tmp, ds["name"] * ext)) == read(joinpath(HandoverRef.NOV20, ds["name"] * ext))
                end
            end
        end
    end
end

@testitem "surface calibration rows reproduce calibration_local.csv (d = 3)" tags = [:handover, :slow] setup = [HandoverRef] begin
    using FTShadows: surface_calibration, calibration_local_rows, write_csv_table, CSVTable
    if HandoverRef.have(HandoverRef.CAL)
        st = surface_calibration([(3, string(b), r) for b in "XZY" for r in (0, 3, 6)])
        rows = calibration_local_rows(st)
        cols = collect(keys(rows[1]))
        lines = split(read(HandoverRef.CAL, String), '\n')
        mktempdir() do tmp
            write_csv_table(joinpath(tmp, "c.csv"), CSVTable(cols, [Any[r[c] for c in cols] for r in rows], "\n"))
            mine = split(read(joinpath(tmp, "c.csv"), String), '\n')
            @test mine[1] == lines[1]
            @test mine[2:10] == filter(l -> startswith(l, "surface,3,"), lines)
        end
    end
end
