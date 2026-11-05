# Library pieces behind the study scripts: numpy and pandas compatibility, the postselection tables and
# the surface calibration seed schedule.

@testmodule StudiesRef begin
    using FTShadows
    const REPO = pkgdir(FTShadows)
    have(paths...) = all(ispath, paths) || (@info "skipping: $(join(filter(!ispath, collect(paths)), ", ")) not found"; false)
end

@testitem "numpy std and packbits are reproduced" tags = [:studies, :bridge] begin
    using PythonCall, Random
    using FTShadows: np_std, packbits
    np = pyimport("numpy")
    rng = Xoshiro(5)
    for n in (2, 7, 8, 9, 127, 128, 129, 1000, 20001)
        a = randn(rng, n) .* 10.0 .^ (2 .* randn(rng, n))
        @test np_std(a; ddof = 1) === pyconvert(Float64, np.std(np.array(a); ddof = 1))
        b = rand(rng, Bool, n)
        @test packbits(b) == pyconvert(Vector{Int}, np.packbits(np.array(b)).tolist())
    end
end

@testitem "pandas_csv writes what DataFrame.to_csv writes" tags = [:studies, :bridge] begin
    using PythonCall
    using FTShadows: pandas_csv, OrderedDict
    pd = pyimport("pandas")
    rows = [OrderedDict{String,Any}("i" => 1, "in" => 3, "f" => 0.1, "b" => true, "s" => "a;b", "n" => nothing,
                                    "inf" => Inf, "bn" => true, "x" => 1e-5),
            OrderedDict{String,Any}("i" => 2, "in" => nothing, "f" => NaN, "b" => false, "s" => "c", "n" => nothing,
                                    "inf" => 2, "bn" => nothing, "extra" => "only here")]
    pyval(x) = x === nothing ? pybuiltins.None : x
    mktempdir() do tmp
        pandas_csv(joinpath(tmp, "j.csv"), rows)
        pd.DataFrame(pylist([pydict(Dict(k => pyval(v) for (k, v) in r)) for r in rows])).to_csv(joinpath(tmp, "p.csv"); index = false)
        # pandas' column order follows the dicts' order; ours is the OrderedDict order, so compare by column
        j = split.(split(chomp(read(joinpath(tmp, "j.csv"), String)), '\n'), ',')
        p = split.(split(chomp(read(joinpath(tmp, "p.csv"), String)), '\n'), ',')
        @test Dict(zip(j[1], zip(j[2], j[3]))) == Dict(zip(p[1], zip(p[2], p[3])))
    end
end

@testitem "postselection tables equal the stored ones" tags = [:studies] setup = [StudiesRef] begin
    using FTShadows: postselection_analyse, pandas_csv
    R = joinpath(StudiesRef.REPO, "results", "postselection")
    D = joinpath(StudiesRef.REPO, "handover", "data", "nov20")
    if StudiesRef.have(R, D)
        for name in ("surface_d3_bell_random_p1e-3_r3", "surface_d3_bell_random_p1e-3_r3_logicalZ0_pL0.05")
            _, rows, frows = postselection_analyse(joinpath(D, name))
            mktempdir() do tmp
                pandas_csv(joinpath(tmp, "f.csv"), frows)
                @test read(joinpath(tmp, "f.csv"), String) == read(joinpath(R, "$(name)_fidelity.csv"), String)
                if !isempty(rows)
                    pandas_csv(joinpath(tmp, "e.csv"), rows)
                    @test read(joinpath(tmp, "e.csv"), String) == read(joinpath(R, "$name.csv"), String)
                end
            end
        end
    end
end

@testitem "surface calibration seeds are unique across configurations" tags = [:studies] begin
    # the full (non-quick) schedule of scripts/surface_calibration.jl, run with a fake pool that records the
    # jobs and returns no errors, so every config runs to its shot cap
    using FTShadows: surface_qr_configs, surface_state_independence_configs, surface_qr_adaptive
    jobs = []
    fake(f, js) = (append!(jobs, js); [(0, j[7], Dict("merged_rounds" => 0)) for j in js])
    target, max5, max7 = 300, 4e6, 2e6
    cfgs = surface_qr_configs()
    for (d, mx) in ((5, max5), (7, max7))
        surface_qr_adaptive(filter(c -> c["d"] == d, cfgs), target, mx; pmap = fake)
    end
    surface_qr_adaptive(surface_state_independence_configs(), 10 * target, 2e6; pmap = fake)
    owners = Dict{Int,Set{Tuple}}()
    for (d, kind, term, st, r, seed, shots, noise) in jobs
        push!(get!(owners, seed, Set{Tuple}()), (d, kind, term, st, r, noise))
    end
    clash = filter(p -> length(p.second) > 1, owners)
    # the worst case: same circuit apart from the ideal input state -> identical samples
    same_circuit = filter(p -> length(Set((c[1], c[2], c[3], c[5], c[6]) for c in p.second)) < length(p.second), clash)
    @test isempty(same_circuit)
    @test isempty(clash)
    @test length(jobs) > 0
end
