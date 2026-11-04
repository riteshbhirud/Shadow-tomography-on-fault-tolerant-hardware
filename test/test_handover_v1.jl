# The v1.0/v1.1 datasets: seeded random states, the learnable targets and the regenerated files equal
# the stored ones. The data files are local, so those items skip when absent.

@testitem "seeded random stabilizer and CSS states match the reference" tags = [:handover, :bridge, :pyref] begin
    using PythonCall
    using FTShadows.Circuits: random_stabilizer_state, random_css_state
    lg = pyimport("ftoracle.logical")
    conv(x) = [(pyconvert(String, t[0]), pyconvert(Int, t[1])) for t in x]
    for k in 1:4, s in 0:5
        @test random_stabilizer_state(k, s) == conv(lg.random_stabilizer_state(k, s))
    end
    for m in (4, 11), s in 0:3
        @test random_css_state(m, s) == conv(lg.random_css_state(m, s))
    end
end

@testitem "v1 dataset specs and exact tables equal the handover files" tags = [:handover] begin
    using JSON
    using FTShadows: v1_datasets, learnable, product_strings, py_round
    dir = joinpath(pkgdir(FTShadows), "handover", "data")
    specs = v1_datasets()
    @test [s["name"] for s in specs][[1, 5, 10]] == ["surface_d3_bell_p0", "bb72_randomcss12_p0", "bb72_randomcss12_p1e-3"]
    # one setting with strings XI, IX gives XI, IX, XX
    lam = learnable(Dict("local:XX" => ["XI", "IX"]), Dict("local:XX" => 0.5))
    @test sort(collect(keys(lam))) == ["IX", "XI", "XX"]
    if isdir(dir)
        for s in specs
            f = joinpath(dir, s["name"] * ".meta.json")
            isfile(f) || continue
            meta = JSON.parsefile(f)
            @test [(g["generator"], g["sign"] == -1 ? 1 : 0) for g in meta["state_for_testing"]] == s["state"]
            @test sort(collect(keys(meta["ensemble"]))) == sort(collect(keys(s["ensemble"])))
        end
    end
end

@testitem "v1 datasets regenerate byte for byte" tags = [:handover, :slow] begin
    using JSON
    using FTShadows: v1_datasets, generate_v1, pyjson
    dir = joinpath(pkgdir(FTShadows), "handover", "data")
    if isdir(dir)
        mktempdir() do tmp
            for s in v1_datasets()
                s["name"] == "bb72_randomcss12_p1e-3" && continue  # about 10 CPU-minutes
                generate_v1(s, tmp)
                for ext in (".records.jsonl", ".exact.csv")
                    @test read(joinpath(tmp, s["name"] * ext)) == read(joinpath(dir, s["name"] * ext))
                end
                a, b = JSON.parsefile(joinpath(tmp, s["name"] * ".meta.json")), JSON.parsefile(joinpath(dir, s["name"] * ".meta.json"))
                foreach(k -> (delete!(a, k); delete!(b, k)), ("created", "runtime_s", "generator", "code_version"))
                @test pyjson(a; indent = 1) == pyjson(b; indent = 1)
            end
        end
    end
end
