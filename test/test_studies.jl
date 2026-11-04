# The studies behind results/ (surface calibration, Y verification, Gidney reproduction, end-to-end
# check, Bravyi calibration, dataset-setting flip rates, postselection) and the library pieces they
# needed: each is compared with the original implementation through the bridge, using its local copy in
# reference/python (see runtests.jl).

@testmodule StudiesRef begin
    using FTShadows, PythonCall
    const REPO = pkgdir(FTShadows)
    const PYSCRIPTS = joinpath(REPO, "reference", "python", "scripts")

    """A reference script loaded as a Python module (its main() is not run)."""
    function script(name)
        util = pyimport("importlib.util")
        spec = util.spec_from_file_location(name, joinpath(PYSCRIPTS, name * ".py"))
        m = util.module_from_spec(spec)
        spec.loader.exec_module(m)
        return m
    end

    have(paths...) = all(ispath, paths) || (@info "skipping: $(join(filter(!ispath, collect(paths)), ", ")) not found"; false)
end

@testitem "random_stabilizer_state equals logical.random_stabilizer_state" tags = [:studies, :bridge, :pyref] begin
    using PythonCall
    using FTShadows.Circuits: random_stabilizer_state
    lg = pyimport("ftoracle.logical")
    for k in 1:4, seed in (0, 1, 7, 1805903980)
        ref = [(pyconvert(String, l), pyconvert(Int, s)) for (l, s) in lg.random_stabilizer_state(k, seed)]
        @test random_stabilizer_state(k, seed) == ref
    end
    @test random_stabilizer_state(3, 5; depth = 4) ==
          [(pyconvert(String, l), pyconvert(Int, s)) for (l, s) in lg.random_stabilizer_state(3, 5; depth = 4)]
end

@testitem "reference estimator equals reference_estimator.py" tags = [:studies, :bridge, :pyref] begin
    using PythonCall, Random
    using FTShadows: decompose, lambdas, estimate, pauli_mul, OrderedDict, product_strings
    re = pyimport("ftoracle.reference_estimator")
    pyint(x) = pyconvert(Int, x)
    for a in ("XYZI", "ZZXY"), b in ("YXZZ", "IXYZ")
        k, l = re.mul(a, b)
        @test pauli_mul(a, b) == (pyint(k), pyconvert(String, l))
    end
    labs = ["ZZI", "IZZ", "XXX"]
    for P in product_strings("IXYZ", 3)
        r = re.decompose(P, pylist(labs))
        j = decompose(P, labs)
        @test pyis(r, pybuiltins.None) == (j === nothing)
        j === nothing || @test j == (pyconvert(Vector{Int}, r[0]), pyint(r[1]))
    end
    # a small ensemble with random records
    ens = OrderedDict("local:XX" => 0.25, "local:ZZ" => 0.25, "surgery:ZZ" => 0.5)
    labels = Dict("local:XX" => ["XI", "IX"], "local:ZZ" => ["ZI", "IZ"], "surgery:ZZ" => ["ZZ"])
    targets = [t for t in product_strings("IXYZ", 2) if t != "II"]
    rng = Xoshiro(1)
    recs = [(s = rand(rng, collect(keys(ens))); Dict{String,Any}("setting_labels" => labels[s],
            "outcomes" => rand(rng, 0:1, length(labels[s])))) for _ in 1:300]
    lam = lambdas(targets, ens, labels)
    plam = re.lambdas(pylist(targets), pydict(ens), pydict(Dict(k => pylist(v) for (k, v) in labels)))
    @test all(lam[P] === pyconvert(Float64, plam[P]) || (lam[P] == 0 && pyint(plam[P]) == 0) for P in targets)
    pyrecs = pylist([pydict(Dict("setting_labels" => pylist(r["setting_labels"]), "outcomes" => pylist(r["outcomes"])))
                     for r in recs])
    pest = re.estimate(pylist(targets), pyrecs, plam)
    est = estimate(targets, recs, lam)
    same(a, b) = a === b || (isnan(a) && isnan(b))
    for P in targets
        m, se, n = pest[P]
        @test same(est[P][1], pyconvert(Float64, m)) && same(est[P][2], pyconvert(Float64, se)) && est[P][3] == pyint(n)
    end
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

@testitem "surface study jobs equal the reference scripts' jobs" tags = [:studies, :bridge, :pyref] setup = [StudiesRef] begin
    using PythonCall
    using FTShadows: surface_qr_job, surface_y_job, gidney_repro_job, dataset_settings_job, seed_for, pyjson
    # lattice surgery chunk (taskA_surface_calibration.work) with sd6 noise and 'same' termination
    m = StudiesRef.script("taskA_surface_calibration")
    for job in ((5, "ZZ", "same", "01", 2, 11, 3000, "sd6"), (5, "XX", "conjugate", "bell", 1, 12, 3000, "bravyi_cal"))
        e, n, cost = surface_qr_job(job)
        pe, pn, pcost = m.work(pytuple(job))
        @test (e, n) == (pyconvert(Int, pe), pyconvert(Int, pn))
        @test pyjson(cost) == pyconvert(String, pyimport("json").dumps(pcost))
    end
    # in-place Y chunk (taskA_surface_y_calibration.work), with a smaller chunk
    m = StudiesRef.script("taskA_surface_y_calibration")
    m.CHUNK = 4000
    for (d, b, st) in ((3, 0, "+1"), (3, 1, "-1"), (5, 2, "zero"))
        c, e, _ = surface_y_job((d, b, st, 3, 4000))
        pc, pe, _ = m.work(pytuple((d, b, st, 3)))
        @test e == pyconvert(Int, pe)
    end
    # Gidney's circuits with plain matching (taskA_surface_y_gidney_reproduction.work)
    m = StudiesRef.script("taskA_surface_y_gidney_reproduction")
    m.CHUNK = 5000
    for (b, d) in (("Y", 3), ("X", 5))
        @test gidney_repro_job((b, d, 2, 5000))[3] == pyconvert(Int, m.work(pytuple((b, d, 2)))[2])
    end
    # one surface dataset setting (calib_dataset_settings.job), all of its rows
    m = StudiesRef.script("calib_dataset_settings")
    rows = dataset_settings_job(("surface", "surgery:XX", 3, 5000))
    prows = m.job(pytuple(("surface", "surgery:XX", 3, 5000)))
    @test length(rows) == pyconvert(Int, pybuiltins.len(prows))
    for (r, pr) in zip(rows, prows)
        @test pyjson(r) == pyconvert(String, pyimport("json").dumps(pr))
    end
end

@testitem "Bravyi seeds, fit and per-cycle rates equal the reference's" tags = [:studies, :bridge, :pyref] setup = [StudiesRef] begin
    using PythonCall
    using FTShadows: bravyi_seed, bravyi_tag, bravyi_fit_pl, bravyi_rows, bravyi_configs, bravyi_key, pyjson, OrderedDict
    m = StudiesRef.script("calib_bravyi_reproduction")
    for k in bravyi_configs(), c in (0, 17)
        tag = bravyi_tag(k, 1.0)
        @test tag == "bravyi|$(pyconvert(String, pystr(k[1])))|$(k[2])|$(k[3])|1.0"
        @test bravyi_seed(tag, c) == pyconvert(Int, m.seed_of(tag, c))
    end
    for p in (1e-3, 2e-3, 3e-3)
        @test bravyi_fit_pl(p) === pyconvert(Float64, m.fit_pl(p))
    end
    # aggregation of a made-up chunk directory, against the reference's main(--aggregate)
    mktempdir() do tmp
        cdir = joinpath(tmp, "results", "calibration", "bravyi_chunks_ms1.0")
        mkpath(cdir)
        for (i, k) in enumerate(bravyi_configs()), c in 0:2
            write(joinpath(cdir, "$(bravyi_key(k))_c$c.json"),
                  pyjson(OrderedDict("fails" => (7i + 3c) % 11, "shots" => 2000, "seed" => 0)))
        end
        m.ROOT = pyimport("pathlib").Path(tmp)
        sys = pyimport("sys")
        sys.argv = pylist(["calib_bravyi_reproduction.py", "--aggregate"])
        m.main()
        ref = pyimport("json").load(pybuiltins.open(joinpath(tmp, "results", "calibration", "bravyi_reproduction.json")))
        @test pyjson(bravyi_rows(cdir, bravyi_configs()); indent = 1) ==
              pyconvert(String, pyimport("json").dumps(ref["rows"]; indent = 1))
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

@testitem "end-to-end check equals end_to_end_surface.run" tags = [:studies, :bridge, :slow, :pyref] setup = [StudiesRef] begin
    using PythonCall
    using FTShadows: e2e_run, pandas_csv
    m = StudiesRef.script("end_to_end_surface")
    rows = e2e_run("bell", "native", 1e-3, 2000, 9)
    prows = m.run("bell", "native", 1e-3, 2000, 9)
    mktempdir() do tmp
        pandas_csv(joinpath(tmp, "j.csv"), rows)
        pyimport("pandas").DataFrame(prows).to_csv(joinpath(tmp, "p.csv"); index = false)
        @test read(joinpath(tmp, "j.csv"), String) == read(joinpath(tmp, "p.csv"), String)
    end
end

@testitem "surface calibration seeds are unique across configurations" tags = [:studies] begin
    # The full (non-quick) schedule of scripts/surface_calibration.jl run
    # with a fake pool that records the jobs and returns no errors, so every config runs to its shot cap.
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
