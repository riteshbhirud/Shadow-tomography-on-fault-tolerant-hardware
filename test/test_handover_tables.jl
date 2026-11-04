# The handover tables that join earlier results (calibration.csv, calibration_local.csv,
# native_ensembles/, bases/) equal the stored files.
# The data are local (handover/, results/ are not tracked), so the file checks skip when absent.

@testmodule TablesRef begin
    using FTShadows
    const REPO = pkgdir(FTShadows)
    have(paths...) = all(ispath, paths) || (@info "skipping: $(join(filter(!ispath, collect(paths)), ", ")) not found"; false)
    # write a CSVTable / JSON record into a temp file and compare bytes with the stored file
    function same_file(write_fn, ref)
        mktempdir() do tmp
            f = joinpath(tmp, basename(ref))
            write_fn(f)
            ok = read(f) == read(ref)
            ok || @error "differs from the stored file" ref
            return ok
        end
    end
end

@testitem "Fraction.limit_denominator matches Python" tags = [:handover, :bridge] begin
    using PythonCall, Random
    using FTShadows: limit_denominator, pyfraction
    fr = pyimport("fractions")
    rng = Xoshiro(7)
    xs = vcat([Float64(Float32(a / b)) for a in 0:40 for b in 1:40 if a <= b], rand(rng, 300), rand(rng, 50) .* 1e-6)
    for x in xs, maxden in (100000, 7)
        mine = pyfraction(limit_denominator(Rational{BigInt}(x), maxden))
        @test mine == pyconvert(String, pystr(fr.Fraction(x).limit_denominator(maxden)))
    end
end

@testitem "pandas read/write typing matches pandas" tags = [:handover, :bridge] begin
    using PythonCall
    using FTShadows: pandas_read, pandas_table, write_csv_table
    pd = pyimport("pandas")
    mktempdir() do tmp
        src = joinpath(tmp, "in.csv")
        write(src, "a,b,c,s,e\n1,0.5,,+1,x\n2,3,1e-05,-1,\n3,7,2,00,y\n")
        rows = pandas_read(src; str_cols = ("s",))
        @test rows[1]["a"] === 1 && rows[2]["b"] === 3.0 && isnan(rows[1]["c"]) && rows[3]["s"] == "00"
        @test ismissing(rows[2]["e"])
        # drop a key from one row: pandas fills NaN and turns the int column into floats
        out = [Dict{String,Any}(k => v for (k, v) in r if !(i == 2 && k == "a")) for (i, r) in enumerate(rows)]
        write_csv_table(joinpath(tmp, "jl.csv"), pandas_table(out, ["a", "b", "c", "s", "e"]))
        # the same rows through pandas, built as the scripts build them (DataFrame of dicts)
        pyrows = pylist([pydict(Dict(k => ismissing(v) ? NaN : v for (k, v) in r)) for r in out])
        pd.DataFrame(pyrows).__getitem__(pylist(["a", "b", "c", "s", "e"])).to_csv(joinpath(tmp, "py.csv"), index = false)
        @test read(joinpath(tmp, "jl.csv"), String) == read(joinpath(tmp, "py.csv"), String)
    end
end

@testitem "calibration.csv and calibration_local.csv equal the handover files" tags = [:handover] setup = [TablesRef] begin
    using FTShadows: calibration_rows, calibration_local_join, pandas_table, write_csv_table,
                     CALIBRATION_COLUMNS, CALIBRATION_LOCAL_COLUMNS
    R, res = TablesRef.REPO, joinpath(TablesRef.REPO, "results", "taskA")
    cal, loc = joinpath(R, "handover", "calibration.csv"), joinpath(R, "handover", "calibration_local.csv")
    if TablesRef.have(cal, loc, joinpath(R, "results", "calibration", "dataset_settings_q.csv"))
        rows = calibration_rows(R)
        @test length(rows) == 368
        @test TablesRef.same_file(f -> write_csv_table(f, pandas_table(rows, CALIBRATION_COLUMNS)), cal)
        rows = calibration_local_join(joinpath(res, "surface_local_q.csv"), joinpath(res, "bb72_local_q.csv"))
        @test TablesRef.same_file(f -> write_csv_table(f, pandas_table(rows, CALIBRATION_LOCAL_COLUMNS)), loc)
    end
end

@testitem "native-ensemble λ tables equal the handover files" tags = [:handover, :slow] setup = [TablesRef] begin
    using FTShadows: native_ensemble_table, read_npy, write_csv_table, NATIVE_ENSEMBLES, NATIVE_EXCLUDE
    tb, ne = joinpath(TablesRef.REPO, "results", "taskB"), joinpath(TablesRef.REPO, "handover", "native_ensembles")
    if TablesRef.have(tb, ne)
        for code in ("BB72", "gross")
            lam = Dict(e => read_npy(joinpath(tb, "lambda_mix_$(code)_$e.npy"))
                       for e in NATIVE_ENSEMBLES if !((code, e) in NATIVE_EXCLUDE))
            @test TablesRef.same_file(f -> write_csv_table(f, native_ensemble_table(code, lam)),
                                      joinpath(ne, "$(code)_lambda_w3.csv"))
        end
    end
end

@testitem "basis exports equal handover/bases" tags = [:handover] setup = [TablesRef] begin
    using FTShadows: handover_bases, write_meta
    dir = joinpath(TablesRef.REPO, "handover", "bases")
    bases = handover_bases()
    @test first.(bases) == ["bb72_oracle_basis.json", "bb72_local_basis.json", "bb72_taskB_basis.json", "gross_taskB_basis.json"]
    @test all(length(b["xbar_support"]) == 12 for (_, b) in bases)
    if TablesRef.have(dir)
        for (fname, rec) in bases
            @test TablesRef.same_file(f -> write_meta(f, rec), joinpath(dir, fname))
        end
    end
end
