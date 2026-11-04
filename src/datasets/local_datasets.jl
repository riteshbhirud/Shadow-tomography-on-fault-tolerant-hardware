# The Nov-20 Bell / GHZ datasets: surface d = 3, 5 and
# BB72 (local basis); fixed or random local Pauli settings; noiseless, or p = 1e-3 at r = 3, 6; optional
# logical Z noise. Sampling and decoding go through `measure`'s path (`run_oracle`), with the reference's
# seeds and chunks, so the records equal the stored ones.

const FORMAT_VERSION = "1.1"
const NOV20_STATES = Dict("bell" => [("XX", 0), ("ZZ", 0)], "ghz" => [("XXX", 0), ("ZZI", 0), ("IZZ", 0)])

"""Strings over `letters` of length m in itertools.product order (last character fastest)."""
product_strings(letters, m) = [join(reverse(t)) for t in vec(collect(Iterators.product(fill(letters, m)...)))]

"""
    nov20_datasets() -> Vector{OrderedDict}

Specs of the 36 Nov-20 datasets (name, code, d, state, m, p, r, kind, shots, pL, pL_qubit).
"""
function nov20_datasets()
    out = OrderedDict{String,Any}[]
    for (code, sizes) in (("surface", (3, 5)), ("bb72", (nothing,)))
        for d in sizes
            for st in ("bell", "ghz")
                m = length(NOV20_STATES[st][1][1])
                pre = code == "surface" ? "surface_d$d" : "bb72"
                nfix, nrand = code == "surface" ? ((m == 2 ? 2000 : 1000), 50000) : ((m == 2 ? 200 : 100), 2000)
                # noiseless: surface measures right after prep (r = 0); BB72 uses the real r = 3 gauging circuits
                for (p, r) in ((0.0, code == "surface" ? 0 : 3), (1e-3, 3), (1e-3, 6))
                    code == "bb72" && st == "ghz" && r == 6 && continue
                    tag = p == 0 ? "p0" : "p1e-3_r$r"
                    for (kind, shots) in (("fixed", nfix), ("random", nrand))
                        push!(out, OrderedDict{String,Any}("name" => "$(pre)_$(st)_$(kind)_$tag", "code" => code, "d" => d,
                            "state" => st, "m" => m, "p" => p, "r" => r, "kind" => kind, "shots" => shots, "pL" => 0.0))
                    end
                end
            end
            if d in (3, nothing)  # logical-noise option: Z on qubit 0 with p_L = 0.05
                pre = code == "surface" ? "surface_d$d" : "bb72"
                push!(out, OrderedDict{String,Any}("name" => "$(pre)_bell_random_p1e-3_r3_logicalZ0_pL0.05", "code" => code,
                    "d" => d, "state" => "bell", "m" => 2, "p" => 1e-3, "r" => 3, "kind" => "random",
                    "shots" => code == "bb72" ? 4000 : 50000, "pL" => 0.05, "pL_qubit" => 0))
            end
        end
    end
    return out
end

"""Generators after a logical Z on qubit q (signs of generators with X or Y there flip)."""
flipped_state(state, q) = [(l, s ⊻ Int(l[q+1] in "XY")) for (l, s) in state]

"""(experiment, labels, decoder) of the local setting `bases` through the code's plug-in."""
function nov20_setting(ds, bases, state)
    params = ds["code"] == "surface" ? (d = ds["d"], memory_rounds = ds["r"]) : (rounds = ds["r"],)
    return Circuits.build_setting(ds["code"], "local:" * bases, state; params...)
end

"""
    nov20_job(job) -> (idx, bases, labels, outcomes, confidence, kinds, cost)

One chunk of shots of one setting; `job = (ds, bases, flag, chunk, idx)`. Workers run these.
"""
function nov20_job(job)
    ds, bases, flag, chunk, idx = job
    state = NOV20_STATES[ds["state"]]
    flag && (state = flipped_state(state, ds["pL_qubit"]))
    exp, labels, dec = nov20_setting(ds, bases, state)
    res = run_oracle(exp, ds["p"], length(idx); decoder = dec, seed = seed_for("nov20", ds["name"], bases, flag, chunk),
                     setting_labels = labels)
    return idx, bases, labels, res.outcomes, res.confidence, res.confidence_kind, res.cost
end

"""
    nov20_draws(ds) -> (sids, draw, flags)

Setting ids, the setting index of every shot, and the logical-noise flags. The random draws are numpy's
`default_rng(seed)` streams, taken through the bridge.
"""
function nov20_draws(ds)
    sids = product_strings("XYZ", ds["m"])
    np = pyimport("numpy")
    if ds["kind"] == "fixed"
        draw = repeat(0:length(sids)-1; inner = ds["shots"])
    else
        rng = np.random.default_rng(seed_for("nov20_draws", ds["name"]))
        draw = pyconvert(Vector{Int}, rng.integers(length(sids); size = ds["shots"]))
    end
    n = length(draw)
    flags = falses(n)
    if ds["pL"] > 0
        u = pyconvert(Vector{Float64}, np.random.default_rng(seed_for("nov20_logical_noise", ds["name"])).random(n))
        flags = u .< ds["pL"]
    end
    return sids, draw, flags
end

_confidence_cell(c::Float64) = isnan(c) ? nothing : isinf(c) ? "inf" : py_round(c, 6)

"""
    simulate_nov20(ds; pmap=map) -> (sids, records, costs, labels_of, flags)

All shots of a dataset. `pmap(f, jobs)` runs the chunks, e.g. `Distributed.pmap` over worker processes
(PythonCall is not thread-safe).
"""
function simulate_nov20(ds; pmap = map)
    sids, draw, flags = nov20_draws(ds)
    chunk = ds["code"] == "surface" ? 10000 : 100
    jobs = Tuple[]
    for (si, b) in enumerate(sids), fl in (false, true)
        idx = findall((draw .== si - 1) .& (flags .== fl)) .- 1
        for c in 0:chunk:length(idx)-1
            push!(jobs, (ds, b, fl, c ÷ chunk, idx[c+1:min(c + chunk, end)]))
        end
    end
    sort!(jobs; by = j -> -length(j[5]) * (1 + count(==('Y'), j[2])), alg = MergeSort)
    records = Vector{Any}(nothing, length(draw))
    costs, labels_of = Dict{String,Any}(), Dict{String,Vector{String}}()
    for (idx, b, labels, outc, conf, kinds, cost) in pmap(nov20_job, jobs)
        costs[b], labels_of[b] = cost, labels
        for (j, shot) in enumerate(idx)
            r = OrderedDict{String,Any}()
            r["shot"] = shot
            r["setting_id"] = "local:$b"
            r["setting"] = labels
            r["outcomes"] = Int.(outc[j, :])
            r["confidence"] = _confidence_cell(conf[j])
            r["confidence_kind"] = kinds[j]
            r["cost"] = cost
            records[shot+1] = r
        end
    end
    @assert all(!isnothing, records)
    return sids, records, costs, labels_of, flags
end

"""
    nov20_q_params(ds, calibration) -> (params, q_of)

Calibrated flip-rate parameters `key => (q, sigma)`, sigma = 95% Wilson CI width / 3.92, from the rows of
`calibration_local.csv`, and `q_of(P, S, params)`. Surface: one q per patch basis, and patches are independent;
BB72: one q per (setting, target).
"""
function nov20_q_params(ds, calibration)
    ds["p"] == 0 && return OrderedDict{Any,Tuple{Float64,Float64}}(), (P, S, p_) -> 0.0
    prm = OrderedDict{Any,Tuple{Float64,Float64}}()
    if ds["code"] == "surface"
        for row in calibration
            (row["code"] == "surface" && row["distance"] == ds["d"] && row["r"] == ds["r"]) || continue
            prm[row["target"]] = (Float64(row["q"]), (Float64(row["ci_hi"]) - Float64(row["ci_lo"])) / (2 * Z95))
        end
        function q_surface(P, S, p_)
            used = [i for i in 1:length(P) if P[i] != 'I']
            any(!haskey(p_, string(S[i])) for i in used) && return NaN
            return (1 - foldl(*, [1 - 2 * p_[string(S[i])] for i in used])) / 2
        end
        return prm, q_surface
    end
    for row in calibration
        (row["code"] == "bb72" && row["r"] == ds["r"]) || continue
        prm[(row["setting_id"], row["target"])] = (Float64(row["q"]), (Float64(row["ci_hi"]) - Float64(row["ci_lo"])) / (2 * Z95))
    end
    return prm, (P, S, p_) -> get(p_, ("local:$S", P), NaN)
end

"""All 2^m - 1 non-identity stabilizer-group elements of `state` as (label, sign ±1), with Pauli phases."""
function stabilizer_group(state)
    m = length(state[1][1])
    out = Tuple{String,Int}[]
    for mask in 1:2^length(state)-1
        acc = Circuits._label_pauli("I"^m)
        for (j, (lab, s)) in enumerate(state)
            (mask >> (j - 1)) & 1 == 1 && (acc = acc * Circuits._label_pauli(lab, s))
        end
        lab = join((x || z) ? (x && z ? 'Y' : x ? 'X' : 'Z') : 'I' for (x, z) in (acc[i] for i in 1:m))
        push!(out, (lab, acc.phase[] == 0x0 ? 1 : -1))
    end
    return out
end

"""
    nov20_expectations(ds, sids, calibration) -> (rows, known)

Ideal, after-logical-channel and noise-expected ⟨P⟩ for every Pauli on the m qubits, and the fidelities.
Each calibrated q is shifted by its sigma in turn for the uncertainty.
"""
function nov20_expectations(ds, sids, calibration)
    m, state = ds["m"], NOV20_STATES[ds["state"]]
    prm, q_of = nov20_q_params(ds, calibration)
    paulis = filter(P -> any(!=('I'), P), product_strings("IXYZ", m))
    sets = Dict(P => [S for S in sids if all(P[i] in ('I', S[i]) for i in 1:m)] for P in paulis)
    ex = Dict(P => Circuits.exact_expectation(state, P) for P in paulis)
    pL = ds["pL"]
    exL = Dict{String,Real}(P => (pL > 0 && P[ds["pL_qubit"]+1] in "XY") ? ex[P] * (1 - 2 * pL) : ex[P] for P in paulis)
    G = stabilizer_group(state)
    function evaluate(p_)
        e = Dict(P => exL[P] * np_mean([1 - 2 * q_of(P, S, p_) for S in sets[P]]) for P in paulis)
        return e, (1 + foldl(+, (s * e[g] for (g, s) in G); init = 0)) / 2^m
    end
    base = Dict{Any,Float64}(k => v[1] for (k, v) in prm)
    e0, F0 = evaluate(base)
    dE = Dict(P => 0.0 for P in paulis)
    dF = 0.0
    for (k, (q, sig)) in prm
        e1, F1 = evaluate(merge(base, Dict{Any,Float64}(k => q + sig)))
        for P in paulis
            dE[P] += (e1[P] - e0[P])^2
        end
        dF += (F1 - F0)^2
    end
    rows = OrderedDict{String,Any}[]
    for P in paulis
        qs = [q_of(P, S, base) for S in sets[P]]
        push!(rows, OrderedDict{String,Any}("pauli" => P, "weight" => m - count(==('I'), P), "exact" => ex[P],
            "exact_with_logical_noise" => exL[P] isa Integer ? exL[P] : py_round(Float64(exL[P]), 12),
            "expected_noisy" => np_round(e0[P], 12), "expected_noisy_se" => sqrt(dE[P]),
            "lambda_P" => py_round(length(sets[P]) / length(sids), 12),
            "settings" => join(("local:$S" for S in sets[P]), ";"),
            "q_mean" => (1 - np_mean([1 - 2 * q for q in qs])) / 2))
    end
    fid(vals) = (1 + foldl(+, (s * vals[g] for (g, s) in G); init = 0)) / 2^m
    known = OrderedDict{String,Any}(
        "fidelity_ideal" => fid(ex), "fidelity_with_logical_noise" => fid(exL),
        "fidelity_expected_noisy" => isnan(F0) ? nothing : F0,
        "fidelity_expected_noisy_se" => isnan(F0) ? nothing : sqrt(dF),
        "fidelity_formula" => "F = 2^-m sum over the 2^m stabilizer-group elements g of sign(g) <g>; Bell: " *
                              "(1 + <XX> - <YY> + <ZZ>)/4",
        "stabilizer_group" => [OrderedDict{String,Any}("pauli" => g, "sign" => s) for (g, s) in G],
        "note" => "expected_noisy(P) = mean over the settings S measuring P of (1 - 2 q_{P,S}) * " *
                  "exact_with_logical_noise(P); q_mean = (1 - mean_S(1 - 2 q_{P,S}))/2. q from handover/" *
                  "calibration_local.csv, calibrated on +1 eigenstates at the same p and r. se: each calibrated q " *
                  "(sigma = 95% Wilson CI width / 3.92) shifted in turn, contributions added in quadrature, so " *
                  "values sharing a q are fully correlated. null = no calibration for this (code, d, r).")
    return rows, known
end

_od(ps...) = OrderedDict{String,Any}(ps...)

"""
    nov20_meta(ds, sids, records, costs, labels_of, flags, known, runtime; code_version, generator)

The dataset's meta.json content (format v1.1 with the Nov-20 keys).
"""
function nov20_meta(ds, sids, records, costs, labels_of, flags, known, runtime; code_version = "unknown",
                    generator = "scripts/make_local_datasets.jl", created = Libc.strftime("%Y-%m-%d", time()))
    m = ds["m"]
    surface = ds["code"] == "surface"
    code = surface ?
           _od("name" => "rotated surface code", "distance" => ds["d"], "patches" => m,
               "logical_basis" => "per patch: Gidney's patch (src/ftoracle/_vendor/gidney_inplace_y); " *
                                  "Xbar = X on the x = 0 column, Zbar = Z on the y = 0 row, " *
                                  "Ybar = i Xbar Zbar measured in place (arXiv:2302.07395)") :
           _od("name" => "BB72 [[72,12,6]] bivariate bicycle code", "distance" => 6,
               "logical_basis" => "LOCAL basis, handover/bases/bb72_local_basis.json (D-37); dataset qubit i = code " *
                                  "logical qubit i")
    ens = _od(("local:$S" => _od("probability" => 1 / length(sids), "strings" => labels_of[S], "cost" => get(costs, S, nothing))
               for S in sids)...)
    nm = Circuits.NOISE_MODELS[Circuits.DEFAULT_NOISE]
    noisy = ds["p"] > 0
    ydef = surface ? get(Circuits.Y_DEFAULT_ROUNDS, ds["d"], nothing) : nothing
    return _od(
        "format_version" => FORMAT_VERSION, "dataset" => ds["name"], "created" => created, "code" => code,
        "n_logical" => m, "logical_qubit_order" => "character i of every Pauli string = logical qubit i (0-based)",
        "conventions" => _od(
            "outcome_bit" => "0 -> eigenvalue +1, 1 -> eigenvalue -1 of the listed Hermitian string",
            "Y" => "Y_i = i X_i Z_i (the Hermitian logical Y of qubit i)",
            "setting" => "local:B0B1.. measures B_i on qubit i; all strings measured jointly on one shot",
            "confidence" => "complementary gap, minimum over the setting's strings (INTERFACE.md section 2)",
            "state_for_testing" => "'sign' is +1 or -1 (not a sign bit)"),
        "noise" => _od("model" => noisy ? Circuits.DEFAULT_NOISE : "none", "p" => ds["p"],
            "multipliers" => noisy ? _od((string(k) => v for (k, v) in pairs(nm))...) : nothing,
            "description" => noisy ? "circuit-level bravyi_cal channel set (D-29), ideal state preparation" :
                             "noiseless (p = 0); ideal state preparation"),
        "rounds" => _od("value" => ds["r"], "meaning" => surface ?
            "memory rounds on every patch before its destructive measurement; Y adds 1 + b + 1 rounds with " *
            "b = $(ydef === nothing ? "None" : ydef) (D-36)" :
            "gauging rounds r of every setting that needs gauging (readout-only settings: 0)"),
        "rounds_surgery" => ds["r"],
        "decoders" => _od(("local:$S" => surface ? "pymatching" : "relay" for S in sids)...),
        "ensemble" => ens, "draws" => length(records),
        "sampling" => ds["kind"] == "fixed" ?
            "fixed_per_setting: every setting exactly shots_per_setting times, in setting blocks (sorted setting " *
            "order); probability = 1/3^m is the exact fraction" :
            "iid_uniform: each shot draws a uniformly random basis per qubit (random Pauli snapshots)",
        "shots_per_setting" => ds["kind"] == "fixed" ? ds["shots"] : nothing,
        "draw_seed" => ds["kind"] == "random" ? seed_for("nov20_draws", ds["name"]) : nothing,
        "logical_qubits_used" => collect(0:m-1),
        "spectators" => surface ? nothing : "BB72 logical qubits 3..11 (and 2 for Bell) are maximally mixed",
        "logical_noise" => ds["pL"] > 0 ?
            _od("pauli" => "Z", "qubit" => ds["pL_qubit"], "p_L" => ds["pL"],
                "applied" => "independently per shot to the prepared state, before the measurement",
                "shots_affected" => count(flags)) : nothing,
        "known_answers" => known,
        "calibration_source" => noisy ? "handover/calibration_local.csv" : nothing,
        "code_version" => code_version, "generator" => generator,
        "state_for_testing" => [_od("generator" => l, "sign" => s == 1 ? -1 : 1) for (l, s) in NOV20_STATES[ds["state"]]],
        "runtime_s" => py_round(Float64(runtime), 1))
end

const EXPECTED_COLUMNS = ["pauli", "weight", "exact", "exact_with_logical_noise", "expected_noisy", "expected_noisy_se",
    "lambda_P", "settings", "q_mean"]

# pandas writes a column of ints and floats as floats, and NaN as an empty cell
function _pandas_column(vals)
    anyfloat = any(v -> v isa AbstractFloat, vals)
    return [v isa AbstractFloat ? (isnan(v) ? missing : Float64(v)) : anyfloat && v isa Integer ? Float64(v) : v for v in vals]
end

"""
    write_nov20(dir, ds, meta, records, rows)

Write `<name>.meta.json`, `.records.jsonl`, `.exact.csv` (csv module, CRLF) and `.expected.csv` (pandas).
"""
function write_nov20(dir::AbstractString, ds, meta, records, rows)
    mkpath(dir)
    base = joinpath(dir, ds["name"])
    write_meta(base * ".meta.json", meta)
    write_records_raw(base * ".records.jsonl", records)
    write_csv_table(base * ".exact.csv",
        CSVTable(["pauli", "exact", "lambda_P", "settings"], [Any[r["pauli"], r["exact"], r["lambda_P"], r["settings"]] for r in rows], "\r\n"))
    cols = [_pandas_column([r[c] for r in rows]) for c in EXPECTED_COLUMNS]
    write_csv_table(base * ".expected.csv", CSVTable(EXPECTED_COLUMNS, [Any[col[i] for col in cols] for i in eachindex(rows)], "\n"))
end

"""
    generate_nov20(ds, dir; calibration, pmap=map, kw...) -> known answers

Simulate one dataset spec and write its four files into `dir` (never the reference's handover/data).
`calibration` is the rows of calibration_local.csv (`csv_rows`).
"""
function generate_nov20(ds, dir::AbstractString; calibration, pmap = map, kw...)
    t0 = time()
    sids, records, costs, labels_of, flags = simulate_nov20(ds; pmap)
    rows, known = nov20_expectations(ds, sids, calibration)
    meta = nov20_meta(ds, sids, records, costs, labels_of, flags, known, time() - t0; kw...)
    write_nov20(dir, ds, meta, records, rows)
    return known
end
