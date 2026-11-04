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

# Every local handover dataset; skipped in a clone without handover/.
@testitem "every handover dataset round-trips byte for byte" tags = [:records, :slow] begin
    using FTShadows: read_meta, write_meta, read_records_raw, write_records_raw, read_csv_table, write_csv_table
    root = joinpath(dirname(@__DIR__), "handover", "data")
    if !isdir(root)
        @info "handover/data not present; skipping the full round-trip"
    else
        n = Ref(0)
        mktempdir() do tmp
            for (dir, _, files) in walkdir(root), f in files
                src, dst = joinpath(dir, f), joinpath(tmp, f)
                if endswith(f, ".meta.json")
                    write_meta(dst, read_meta(src))
                elseif endswith(f, ".records.jsonl")
                    write_records_raw(dst, read_records_raw(src))
                elseif endswith(f, ".csv")
                    write_csv_table(dst, read_csv_table(src))
                else
                    continue
                end
                n[] += 1
                ok = read(dst) == read(src)
                ok || @error "not byte-identical" src
                @test ok
            end
        end
        @info "round-tripped $(n[]) handover files"
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

# The package loader must give exactly what handover/julia/load_records.jl gives, and pass
# that loader's estimator test, on every dataset present locally.
@testitem "package loader equals the handover loader" tags = [:records, :slow] begin
    using Statistics
    root = joinpath(dirname(@__DIR__), "handover")
    if !isfile(joinpath(root, "julia", "load_records.jl"))
        @info "handover/ not present; skipping the loader comparison"
    else
        include(joinpath(root, "julia", "load_records.jl"))
        prefixes = String[]
        for (dir, _, files) in walkdir(joinpath(root, "data")), f in files
            endswith(f, ".meta.json") && push!(prefixes, joinpath(dir, replace(f, ".meta.json" => "")))
        end
        for p in prefixes
            a, b = load_dataset(p), LoadRecords.load_dataset(p)
            @test a.meta == b.meta
            @test length(a.records) == length(b.records)
            @test all(getfield(x, f) == getfield(y, f) || (isnan(x.confidence) && isnan(y.confidence))
                      for (x, y) in zip(a.records, b.records) for f in fieldnames(Record))
            @test all(x.pauli == y.pauli && x.exact == y.exact && x.lambda == y.lambda && x.settings == y.settings
                      for (x, y) in zip(a.exact, b.exact))
            if all(r -> startswith(r.setting_id, "local:") && all(s -> count(!=('I'), s) == 1, r.setting), a.records)
                @test [(s.bases, s.outcomes) for s in local_view(a.records)] ==
                      [(s.bases, s.outcomes) for s in LoadRecords.local_view(b.records)]
            end
        end
        # the handover test's estimator check (its datasets: data/ without the nov20 subfolder,
        # whose logical-noise sets deliberately differ from the exact values), on the package loader
        for p in filter(p -> dirname(p) == joinpath(root, "data"), prefixes)
            ds = load_dataset(p)
            idx = Dict(s => subset_index(String.(e["strings"])) for (s, e) in ds.meta["ensemble"])
            n = length(ds.records)
            for e in ds.exact
                v = zeros(n)
                for (i, r) in enumerate(ds.records)
                    sub = get(idx[r.setting_id], e.pauli, nothing)
                    sub === nothing || (v[i] = prod(1 .- 2 .* r.outcomes[sub]) / e.lambda)
                end
                se = std(v) / sqrt(n)
                z = se > 0 ? (mean(v) - e.exact) / se : (mean(v) == e.exact ? 0.0 : Inf)
                abs(z) < 5 || @error "estimate far from exact" p e.pauli z
                @test abs(z) < 5
            end
        end
        @info "compared $(length(prefixes)) datasets"
    end
end
