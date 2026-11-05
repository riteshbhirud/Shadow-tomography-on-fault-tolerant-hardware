@testitem "pyfloat matches Python repr" tags = [:records] begin
    using FTShadows: pyfloat
    for (x, s) in [(0.25, "0.25"), (1e-3, "0.001"), (3.73e-5, "3.73e-05"), (1e-4, "0.0001"), (1e16, "1e+16"),
                   (1e15, "1000000000000000.0"), (12.41, "12.41"), (-0.0, "-0.0"), (1 / 3, "0.3333333333333333"),
                   (19.331544864059214, "19.331544864059214"), (2.5e-310, "2.5e-310")]
        @test pyfloat(x) == s
    end
end

@testitem "pyjson matches Python json" tags = [:records] begin
    using FTShadows: pyjson
    using JSON
    d = JSON.parse("""{"a": [1, 2.0, "x\\u00b1"], "b": {"c": null, "d": true}, "e": []}""")
    @test pyjson(d) == """{"a": [1, 2.0, "x\\u00b1"], "b": {"c": null, "d": true}, "e": []}"""
    @test pyjson(d; separators = (",", ":")) == """{"a":[1,2.0,"x\\u00b1"],"b":{"c":null,"d":true},"e":[]}"""
    @test pyjson(d; indent = 1) == "{\n \"a\": [\n  1,\n  2.0,\n  \"x\\u00b1\"\n ],\n \"b\": {\n  \"c\": null,\n  \"d\": true\n },\n \"e\": []\n}"
end

@testitem "tracked sample datasets round-trip byte for byte" tags = [:records] begin
    using FTShadows: read_meta, write_meta, read_records_raw, write_records_raw, read_csv_table, write_csv_table
    dir = joinpath(@__DIR__, "reference", "records")
    files = readdir(dir)
    @test !isempty(files)
    mktempdir() do tmp
        for f in files
            src, dst = joinpath(dir, f), joinpath(tmp, f)
            if endswith(f, ".meta.json")
                write_meta(dst, read_meta(src))
            elseif endswith(f, ".records.jsonl")
                write_records_raw(dst, read_records_raw(src))
            elseif endswith(f, ".csv")
                write_csv_table(dst, read_csv_table(src))
            end
            @test read(dst) == read(src)
        end
    end
end

@testitem "loader API" tags = [:records] begin
    dir = joinpath(@__DIR__, "reference", "records")
    ds = load_dataset(joinpath(dir, "surface_d3_bell_p0"))  # records trimmed to the first 500 shots
    @test length(ds.records) == 500
    @test all(r.confidence_kind == "noiseless" && isinf(r.confidence) for r in ds.records)
    idx = subset_index(["ZZI", "IZZ"])
    @test idx["ZIZ"] == [1, 2]
    @test_throws ErrorException subset_index(["XX", "ZZ"])
    loc = load_dataset(joinpath(dir, "bb72_bell_fixed_p0"))
    @test length(loc.records) == loc.meta["draws"]
    shots = local_view(loc.records)
    @test all(length(s.bases) == 2 for s in shots)
end
