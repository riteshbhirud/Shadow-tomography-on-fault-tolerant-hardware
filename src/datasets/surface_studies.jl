# Surface-code studies behind results/taskA and results/e2e: lattice-surgery q(r), the in-place Y
# measurement, Gidney's memory numbers and the end-to-end check. The drivers in scripts/ only parse
# arguments, start workers and write files. Seeds, chunking and stopping rules are fixed here, so a rerun
# of a configuration gives the same counts.

#### Lattice-surgery flip rate q(r)

const SURFACE_QR_P = 1e-3
const SURFACE_QR_CHUNK = 50_000

# name => (state spec, expected outcome bit of the 2-patch product)
const SURFACE_QR_STATES = OrderedDict(
    "ZZ" => OrderedDict("00" => ([("ZI", 0), ("IZ", 0)], 0), "01" => ([("ZI", 0), ("IZ", 1)], 1),
                        "11" => ([("ZI", 1), ("IZ", 1)], 0), "bell" => ([("ZZ", 0), ("XX", 0)], 0)),
    "XX" => OrderedDict("++" => ([("XI", 0), ("IX", 0)], 0), "+-" => ([("XI", 0), ("IX", 1)], 1),
                        "--" => ([("XI", 1), ("IX", 1)], 0), "bell" => ([("ZZ", 0), ("XX", 0)], 0)))

_qr_cfg(d, kind, term, state, r, noise) =
    OrderedDict{String,Any}("d" => d, "kind" => kind, "term" => term, "state" => state, "r" => r, "noise" => noise)

"""
    surface_qr_configs() -> Vector{OrderedDict}

The configurations of surface_qr.csv: conjugate termination (d = 5 up to r = 10, d = 7 up to r = 12), same
termination (d = 5), sd6 noise (d = 5), stim_generated noise (ZZ).
"""
function surface_qr_configs()
    cfgs = OrderedDict{String,Any}[]
    for (kind, sname) in (("ZZ", "00"), ("XX", "++"))
        for (d, rmax) in ((5, 10), (7, 12)), r in 1:rmax
            push!(cfgs, _qr_cfg(d, kind, "conjugate", sname, r, "bravyi_cal"))
        end
        for r in 1:6
            push!(cfgs, _qr_cfg(5, kind, "same", sname, r, "bravyi_cal"))
        end
        for r in 1:10  # labelled pessimistic variant (idle noise on ancillas too)
            push!(cfgs, _qr_cfg(5, kind, "conjugate", sname, r, "sd6"))
        end
    end
    for r in 1:8  # stim-generated-like noise, for comparison
        push!(cfgs, _qr_cfg(5, "ZZ", "conjugate", "00", r, "stim_generated"))
    end
    return cfgs
end

"""The state-independence configurations: every state of both kinds at d = 5, r = 3."""
surface_state_independence_configs() =
    [_qr_cfg(5, k, "conjugate", s, 3, "bravyi_cal") for k in keys(SURFACE_QR_STATES) for s in keys(SURFACE_QR_STATES[k])]

"""
    surface_qr_job(job) -> (errors, shots, cost)

One chunk; `job = (d, kind, term, state name, r, seed, shots, noise)`. Errors are shots whose product
outcome differs from the state's expected bit.
"""
function surface_qr_job(job)
    d, kind, term, sname, r, seed, shots, noise = job
    state, expect = SURFACE_QR_STATES[kind][sname]
    exp = Circuits.multi_patch_surgery(d, 2, kind, state, r, term)
    res = run_oracle(exp, SURFACE_QR_P, shots; seed, with_gap = false, noise)
    return count(!=(expect), res.outcomes[:, 1]), shots, res.cost
end

"""
    surface_qr_adaptive(configs, target_errors, max_shots; pmap=map, chunk=SURFACE_QR_CHUNK) -> rows

Run chunks two at a time per configuration until it has `target_errors` errors or `max_shots` shots. One
row per configuration: its fields, then p, shots, errors, q, the Wilson interval and the setting's cost.
"""
function surface_qr_adaptive(configs, target_errors, max_shots; pmap = map, chunk::Integer = SURFACE_QR_CHUNK)
    st = [Dict{String,Any}("err" => 0, "shots" => 0, "seed" => 0, "cost" => nothing) for _ in configs]
    active = Set(eachindex(configs))
    while !isempty(active)
        jobs, owners = Tuple[], Int[]
        for i in sort!(collect(active))
            c, s = configs[i], st[i]
            for _ in 1:2
                seed = seed_for("surface_calib", c["d"], c["kind"], c["term"], c["state"], c["r"], c["noise"], s["seed"])
                push!(jobs, (c["d"], c["kind"], c["term"], c["state"], c["r"], seed, chunk, c["noise"]))
                s["seed"] += 1
                push!(owners, i)
            end
        end
        for (i, (e, n, cost)) in zip(owners, pmap(surface_qr_job, jobs))
            st[i]["err"] += e
            st[i]["shots"] += n
            st[i]["cost"] = cost
        end
        filter!(i -> st[i]["err"] < target_errors && st[i]["shots"] < max_shots, active)
    end
    rows = OrderedDict{String,Any}[]
    for (c, s) in zip(configs, st)
        lo, hi = wilson(s["err"], s["shots"])
        row = merge(OrderedDict{String,Any}(), c)
        merge!(row, OrderedDict{String,Any}("p" => SURFACE_QR_P, "shots" => s["shots"], "errors" => s["err"],
            "q" => s["err"] / s["shots"], "ci_lo" => lo, "ci_hi" => hi))
        for (k, v) in s["cost"]
            row["cost_$k"] = v
        end
        push!(rows, row)
    end
    return rows
end

#### In-place Y measurement as a single-patch setting

const SURFACE_Y_P = 1e-3
const SURFACE_Y_STATES = Dict("+1" => ([("Y", 0)], 0), "-1" => ([("Y", 1)], 1), "zero" => ([("Z", 0)], nothing))

"""
    surface_y_job(job) -> (config, errors, cost)

`job = (d, b, state, chunk, chunk_shots)`. For the Z eigenstate ("zero") the count is of outcome 1.
"""
function surface_y_job(job)
    d, b, state, chunk, shots = job
    st, expect = SURFACE_Y_STATES[state]
    exp, _ = Circuits.surface_local_experiment(d, "Y", st; boundary_rounds = b)
    res = run_oracle(exp, SURFACE_Y_P, shots; seed = seed_for("surface_y_calib", d, b, state, chunk), with_gap = false)
    o = res.outcomes[:, 1]
    return (d, b, state), expect === nothing ? count(==(1), o) : count(!=(expect), o), res.cost
end

"""
    surface_y_calibration(; pmap=map, chunk=100_000, target=300, max_shots=Dict(3 => 2e6, 5 => 1e7, 7 => 3e7))

Rows of surface_y_qr.csv: q(b) of the +1 eigenstate for d = 3, 5, 7 and b = 0 .. d÷2 + 2, and at the
default b the -1 eigenstate and P(outcome = 1) of a Z eigenstate (two chunks, no adaptive stopping).
"""
function surface_y_calibration(; pmap = map, chunk::Integer = 100_000, target::Integer = 300,
                               max_shots = Dict(3 => 2_000_000, 5 => 10_000_000, 7 => 30_000_000))
    ydef = Circuits.Y_DEFAULT_ROUNDS
    cfgs = Tuple{Int,Int,String}[(d, b, "+1") for d in (3, 5, 7) for b in 0:d÷2+2]
    append!(cfgs, [(d, ydef[d], s) for d in (3, 5, 7) for s in ("-1", "zero")])
    st = Dict(c => Dict{String,Any}("err" => 0, "shots" => 0, "chunk" => 0, "cost" => nothing) for c in cfgs)
    active = Set(cfgs)
    while !isempty(active)
        jobs = Tuple[]
        for c in sort!(collect(active))
            for _ in 1:2
                push!(jobs, (c..., st[c]["chunk"], chunk))
                st[c]["chunk"] += 1
            end
        end
        for (c, e, cost) in pmap(surface_y_job, jobs)
            st[c]["err"] += e
            st[c]["shots"] += chunk
            st[c]["cost"] = cost
        end
        # the <Y> = 0 check needs no adaptive stopping: 2 chunks
        filter!(c -> c[3] != "zero" && st[c]["err"] < target && st[c]["shots"] < max_shots[c[1]], active)
    end
    rows = OrderedDict{String,Any}[]
    for (d, b, s) in cfgs
        x = st[(d, b, s)]
        lo, hi = wilson(x["err"], x["shots"])
        row = OrderedDict{String,Any}("d" => d, "boundary_rounds" => b, "default" => b == ydef[d], "state" => s,
            "noise" => "bravyi_cal", "p" => SURFACE_Y_P, "shots" => x["shots"], "errors" => x["err"],
            "q" => x["err"] / x["shots"], "ci_lo" => lo, "ci_hi" => hi,
            "quantity" => s == "zero" ? "P(outcome=1)" : "flip rate q")
        for (k, v) in x["cost"]
            row["cost_$k"] = v
        end
        push!(rows, row)
        println(row)
        flush(stdout)
    end
    return rows
end

#### Gidney's Y / X / Z memory numbers with our decoder

"""
    gidney_memory_circuit(b, d) -> Py

Gidney's own circuit from `midout`: basis `b`, distance d, SI1000 noise at p = 1e-3, d memory rounds,
d÷2 boundary rounds for Y.
"""
function gidney_memory_circuit(b::AbstractString, d::Integer)
    sys = pyimport("sys")
    pyin(Circuits.VENDOR, sys.path) || sys.path.insert(0, Circuits.VENDOR)
    gen = pyimport("midout.gen")
    make_circuit = pyimport("midout._make_circuit").make_circuit
    return make_circuit(; basis = b, distance = d, noise = gen.NoiseModel.si1000(1e-3),
                        boundary_rounds = b == "Y" ? d ÷ 2 : 0, memory_rounds = d)
end

"""
    gidney_repro_job(job) -> (b, d, errors)

`job = (b, d, chunk, chunk_shots)`: sample a chunk of Gidney's circuit and decode it with plain matching.
"""
function gidney_repro_job(job)
    b, d, chunk, shots = job
    c = gidney_memory_circuit(b, d)
    dec = MatchingGapDecoder(c.detector_error_model(; decompose_errors = true))
    res = c.compile_detector_sampler(; seed = seed_for("gidney_repro", b, d, chunk)).sample(shots; separate_observables = true)
    dets = pyconvert(Matrix{Bool}, res[0])
    obs = pyconvert(Matrix{Bool}, res[1])
    pred, _ = decode(dec, dets; with_gap = false)
    return b, d, count(s -> (pred[s, 1] != 0) != obs[s, 1], 1:shots)
end

"""
    gidney_reference_stats(path) -> Vector{Dict}

Rows of Gidney's sinter stats.csv with the json_metadata fields (b, d, p, r, rb, ...) merged in.
"""
function gidney_reference_stats(path::AbstractString)
    lines = filter(!isempty, split(read(path, String), '\n'))
    header = strip.(_csv_split(lines[1]))
    rows = Dict{String,Any}[]
    for l in lines[2:end]
        r = Dict{String,Any}(zip(header, strip.(_csv_split(l))))
        for k in ("shots", "errors", "discards")
            r[k] = parse(Int, r[k])
        end
        merge!(r, JSON.parse(r["json_metadata"]))
        push!(rows, r)
    end
    return rows
end

"""
    gidney_reproduction(stats; pmap=map, workers=8, chunk=100_000, target=300, max_shots=20_000_000)

The rows of surface_y_gidney_reproduction.json. Every round runs `workers` chunks per active
configuration, so the counts depend on `--workers`.
"""
function gidney_reproduction(stats; pmap = map, workers::Integer = 8, chunk::Integer = 100_000,
                             target::Integer = 300, max_shots::Integer = 20_000_000)
    cfgs = [(string(b), d) for b in "YXZ" for d in (3, 5, 7)]
    st = Dict(c => Dict("err" => 0, "shots" => 0, "chunk" => 0) for c in cfgs)
    active = Set(cfgs)
    while !isempty(active)
        jobs = Tuple[]
        for c in sort!(collect(active))
            for _ in 1:workers
                push!(jobs, (c..., st[c]["chunk"], chunk))
                st[c]["chunk"] += 1
            end
        end
        for (b, d, e) in pmap(gidney_repro_job, jobs)
            st[(b, d)]["err"] += e
            st[(b, d)]["shots"] += chunk
        end
        filter!(c -> st[c]["err"] < target && st[c]["shots"] < max_shots, active)
    end
    rows = OrderedDict{String,Any}[]
    for (b, d) in cfgs
        rb = b == "Y" ? d ÷ 2 : 0
        ref = stats[findfirst(s -> s["b"] == b && s["d"] == d && s["p"] == 0.001 && s["r"] == d && s["rb"] == rb, stats)]
        e, n = st[(b, d)]["err"], st[(b, d)]["shots"]
        lo, hi = wilson(e, n)
        gl = ref["errors"] / ref["shots"]
        push!(rows, OrderedDict{String,Any}("basis" => b, "d" => d, "r" => d, "rb" => rb, "shots" => n, "errors" => e,
            "ler" => e / n, "ci_lo" => lo, "ci_hi" => hi, "gidney_ler" => gl, "gidney_shots" => ref["shots"],
            "gidney_errors" => ref["errors"], "gidney_decoder" => ref["decoder"],
            "ratio_ours_over_gidney" => (e / n) / gl))
        println(rows[end])
        flush(stdout)
    end
    return rows
end

#### Verification of the in-place Y measurement

"""
    y_reproduction(circuit_dir) -> OrderedDict(d => Bool)

Whether Gidney's generator reproduces his published Y circuit (stripped text) for d = 3, 5, 7.
"""
function y_reproduction(circuit_dir::AbstractString)
    out = OrderedDict{Int,Bool}()
    for d in (3, 5, 7)
        c = pyconvert(String, pystr(gidney_memory_circuit("Y", d)))
        pre = "r=$d,d=$d,p=0.001,noise=SI1000,b=Y,rb=$(d ÷ 2),q="
        f = first(sort!(filter(n -> startswith(n, pre) && endswith(n, ".stim"), readdir(circuit_dir))))
        out[d] = strip(c) == strip(read(joinpath(circuit_dir, f), String))
    end
    return out
end

"""Detector count and fired detectors in 10^4 noiseless shots of the Y setting (stim's DEM must build)."""
function y_determinism()
    out = OrderedDict{Int,Any}()
    for d in (3, 5, 7)
        exp, _ = Circuits.surface_local_experiment(d, "Y", [("Y", 0)])
        c = Circuits.pycircuit(Circuits.build(exp, 0.0).circuit)
        c.detector_error_model()  # raises on non-deterministic detectors
        s = pyconvert(Matrix{Bool}, c.compile_detector_sampler(; seed = 1).sample(10000))
        out[d] = OrderedDict("detectors" => pyconvert(Int, c.num_detectors), "fired" => count(s))
    end
    return out
end

"""
Noiseless exactness: Y eigenstates at d = 3, 5, 7, and every product of every subset of six d = 3
settings on 8 seeded random 3-qubit stabilizer states.
"""
function y_noiseless()
    eigen = OrderedDict{String,Any}[]
    rnd = OrderedDict("states" => 0, "settings" => 0, "products" => 0, "mismatches" => 0)
    for d in (3, 5, 7), sign in (0, 1)
        exp, lab = Circuits.surface_local_experiment(d, "Y", [("Y", sign)])
        r = run_oracle(exp, 0.0, 500; seed = sign, setting_labels = lab)
        push!(eigen, OrderedDict("d" => d, "sign" => sign == 1 ? -1 : 1, "ok" => all(==(sign), r.outcomes[:, 1])))
    end
    for t in 0:7
        st = Circuits.random_stabilizer_state(3, seed_for("surface_y_verify_state", t))
        rnd["states"] += 1
        for bases in ("YXZ", "ZYX", "YYY", "XYY", "YZZ", "ZZY")
            exp, lab = Circuits.surface_local_experiment(3, bases, st)
            r = run_oracle(exp, 0.0, 400; seed = t, setting_labels = lab)
            x = 1 .- 2 .* Int.(r.outcomes)
            rnd["settings"] += 1
            for m in 1:7
                idx = [i for i in 0:2 if (m >> i) & 1 == 1]
                P = join(i in idx ? bases[i+1] : 'I' for i in 0:2)
                e = Circuits.exact_expectation(st, P)
                v = vec(prod(x[:, idx .+ 1]; dims = 2))
                rnd["products"] += 1
                bad = e in (1, -1) ? any(!=(e), v) : length(Set(v)) != 2
                rnd["mismatches"] += Int(bad)
            end
        end
    end
    return OrderedDict("eigen" => eigen, "random" => rnd)
end

"""
Circuit distance (stim search) and shortest graphlike error per (d, b), b from d÷2 to the default.
`seconds` is wall time.
"""
function y_distances()
    out = OrderedDict{String,Any}[]
    for d in (3, 5, 7), b in d÷2:Circuits.Y_DEFAULT_ROUNDS[d]
        t = time()
        exp, _ = Circuits.surface_local_experiment(d, "Y", [("Y", 0)]; boundary_rounds = b)
        c = Circuits.build(exp, 1e-3; eigen_targets = true).circuit
        push!(out, OrderedDict{String,Any}("d" => d, "b" => b, "default" => b == Circuits.Y_DEFAULT_ROUNDS[d],
            "circuit_distance_search" => Circuits.stim_search(c; max_dets = 6, max_degree = 8),
            "graphlike" => pyconvert(Int, pybuiltins.len(Circuits.pycircuit(c).shortest_graphlike_error())),
            "seconds" => py_round(time() - t, 1)))
        println(out[end])
        flush(stdout)
    end
    return out
end

"""X / Z fragments against the builder's per-patch readout at d = 3: flips of each string in 2·10^5 shots."""
function y_readout_crosscheck()
    out = OrderedDict{String,Any}[]
    for bases in ("XZ", "ZX", "XX", "ZZ")
        st = [(join(j == i ? b : 'I' for j in 1:2), 0) for (i, b) in enumerate(bases)]
        frag, lab = Circuits.surface_local_experiment(3, bases, st)
        old = Circuits.multi_patch_readout(3, bases, st)
        n = 200_000
        a = run_oracle(frag, 1e-3, n; seed = seed_for("y_xcheck_frag", bases), with_gap = false).outcomes
        b = run_oracle(old, 1e-3, n; seed = seed_for("y_xcheck_builder", bases), with_gap = false).outcomes
        for j in 1:2
            ea, eb = count(==(1), a[:, j]), count(==(1), b[:, j])
            push!(out, OrderedDict{String,Any}("setting" => bases, "string" => lab[j], "shots" => n,
                "fragment_errors" => ea, "builder_errors" => eb, "fragment_ci" => wilson(ea, n), "builder_ci" => wilson(eb, n)))
        end
    end
    return out
end

#### End-to-end check of the estimator

const E2E_D, E2E_R = 3, 3
const E2E_STATES = OrderedDict("bell" => (2, [("ZZ", 0), ("XX", 0)]),
                               "ghz" => (3, [("ZZI", 0), ("IZZ", 0), ("XXX", 0)]))

"""<P> on a stabilizer state from `decompose` on its generator list: ±1 if P is in the group, else 0."""
function e2e_exact(P, state_gens)
    dec = decompose(P, [g for (g, _) in state_gens])
    dec === nothing && return 0.0
    sub, sign = dec
    return Float64(sign * prod(Int[(-1)^state_gens[i+1][2] for i in sub]))
end

"""The "local" (uniform X/Z per patch) and "native" (1/2 local, 1/4 ZZ-chain, 1/4 XX-chain) ensembles."""
function e2e_ensembles(n)
    loc = ["local:" * b for b in product_strings("XZ", n)]
    e_local = OrderedDict{String,Float64}(s => 1 / length(loc) for s in loc)
    e_native = OrderedDict{String,Float64}(s => 0.5 / length(loc) for s in loc)
    e_native["surgery:ZZ"] = 0.25
    e_native["surgery:XX"] = 0.25
    return OrderedDict("local" => e_local, "native" => e_native)
end

"""
    e2e_run(state, ensemble, p, n_draws, seed) -> rows

Draw `n_draws` settings with numpy's `default_rng(seed)`, run each setting's shots through the oracle,
and compare the reference estimator with the exact value for every non-identity target.
"""
function e2e_run(state_name, ens_name, p, n_draws, seed)
    n, gens = E2E_STATES[state_name]
    ens = e2e_ensembles(n)[ens_name]
    names = collect(keys(ens))
    np = pyimport("numpy")
    draws = pyconvert(Vector{Int}, np.random.default_rng(seed).choice(length(names); size = n_draws,
                                                                       p = pylist([ens[s] for s in names])))
    labels_of = Dict{String,Vector{String}}()
    recs = Dict{String,Any}[]
    for (si, s) in enumerate(names)
        cnt = count(==(si - 1), draws)
        exp, labs = Circuits.setting_experiment(s, E2E_D, n, gens, E2E_R)
        labels_of[s] = labs
        cnt == 0 && continue
        res = run_oracle(exp, p, cnt; seed = seed * 100 + si - 1)
        for i in 1:cnt
            push!(recs, Dict{String,Any}("setting_labels" => labs, "outcomes" => Int.(res.outcomes[i, :])))
        end
    end
    targets = [t for t in product_strings("IXYZ", n) if any(!=('I'), t)]
    lam = lambdas(targets, ens, labels_of)
    est = estimate(targets, recs, lam)
    rows = OrderedDict{String,Any}[]
    for P in targets
        m, se, N = est[P]
        ex = e2e_exact(P, gens)
        push!(rows, OrderedDict{String,Any}("state" => state_name, "ensemble" => ens_name, "p" => p, "P" => P,
            "exact" => ex, "lam" => lam[P], "shadow_norm_sq" => lam[P] > 0 ? 1 / lam[P] : Inf, "estimate" => m,
            "se" => se, "z" => lam[P] > 0 && se > 0 ? (m - ex) / se : NaN, "n" => N))
    end
    return rows
end

"""
    e2e_summary(rows) -> Vector{OrderedDict}

Per (state, ensemble, p), sorted: learnable targets, max |z| (NaN ignored) and the fractions within
2 and 3 standard errors, over targets with lambda > 0.
"""
function e2e_summary(rows)
    out = OrderedDict{String,Any}[]
    for key in sort!(unique([(r["state"], r["ensemble"], r["p"]) for r in rows]))
        g = [r for r in rows if (r["state"], r["ensemble"], r["p"]) == key]
        az = [abs(r["z"]) for r in g if r["lam"] > 0]
        nn = filter(!isnan, az)
        push!(out, OrderedDict{String,Any}("state" => key[1], "ensemble" => key[2], "p" => key[3],
            "n_learnable" => length(az), "n_targets" => length(g),
            "max_abs_z" => isempty(nn) ? NaN : maximum(nn),
            "frac_within_2se" => count(<(2), az) / length(az), "frac_within_3se" => count(<(3), az) / length(az)))
    end
    return out
end
