# The v1.0/v1.1 datasets: surface d = 3 Bell, singlet, GHZ and a random entangled 3-qubit state under a
# mixed local/surgery ensemble, and a BB72 random CSS state under readout plus Ybar_0 gauging, each noiseless
# and at p = 1e-3. Exact values come from the logical state only. Setting draws use numpy's default_rng so
# a regenerated dataset has the same shots as the published one.

const R_SURGERY = 3

"""1/2 random local X/Z readout of every patch (uniform over the 2^k strings), 1/4 ZZ and 1/4 XX surgery."""
function surface_ensemble(k::Integer)
    ens = OrderedDict{String,Float64}("local:$b" => 0.5 / 2^k for b in product_strings("XZ", k))
    ens["surgery:ZZ"] = 0.25
    ens["surgery:XX"] = 0.25
    return ens
end

"""0.4 destructive all-Z readout, 0.4 all-X readout, 0.2 gauging of Ybar_0."""
bb72_ensemble() = OrderedDict{String,Float64}("readout:Z" => 0.4, "readout:X" => 0.4, "gauge:Y0" => 0.2)

"""(experiment, labels, decoder) of setting `sid` of a v1 dataset spec on `state`."""
function v1_setting(spec, sid::AbstractString, state)
    if spec["code"] == "surface"
        exp, labels = Circuits.setting_experiment(sid, 3, spec["k"], state, R_SURGERY)
        return exp, labels, "pymatching"
    elseif startswith(sid, "readout:")
        exp, labels = Circuits.bb72_readout_experiment(sid[end], state)
        return exp, labels, "relay"
    end
    return Circuits.bb_gauging('Y', R_SURGERY; state_labels = state), ["Y" * "I"^11], "relay"
end

const _MUL = Dict{Tuple{Char,Char},Char}()
for x in "IXYZ"
    _MUL[('I', x)] = x
    _MUL[(x, 'I')] = x
end
for x in "XYZ"
    _MUL[(x, x)] = 'I'
end
for (x, y, z) in (('X', 'Y', 'Z'), ('Y', 'Z', 'X'), ('Z', 'X', 'Y'))
    _MUL[(x, y)] = z
    _MUL[(y, x)] = z
end

"""Letter-wise product of two Pauli labels, phase dropped."""
mul_labels(a::AbstractString, b::AbstractString) = join(_MUL[(p, q)] for (p, q) in zip(a, b))

"""
    learnable(ens_labels, ens) -> Dict(P => (λ_P, setting ids))

Every P that is a product of a nonempty subset of one setting's strings. Within a setting each qubit
carries one letter (or I), so the products have no phase. `ens_labels` iterates in setting order.
"""
function learnable(ens_labels, ens)
    out = Dict{String,Tuple{Float64,Vector{String}}}()
    for (sid, labels) in ens_labels
        seen = Set{String}()
        for r in 1:length(labels), sub in _combinations_sorted(length(labels), r)
            P = labels[sub[1]]
            for j in sub[2:end]
                P = mul_labels(P, labels[j])
            end
            (all(==('I'), P) || P in seen) && continue
            push!(seen, P)
            lam, sids = get(out, P, (0.0, String[]))
            out[P] = (lam + ens[sid], vcat(sids, sid))
        end
    end
    return out
end

"""
    v1_datasets(draws=20000, bb_draws=10000) -> Vector{OrderedDict}

Specs of the ten v1 datasets: per p ∈ (0, 1e-3), surface Bell, singlet, GHZ, random3, then BB72
randomcss12. The random3 state is the first seeded random 3-qubit state with no
single-qubit stabilizer element; the BB72 state is a -Ybar_0 eigenstate times a seeded random CSS state.
"""
function v1_datasets(draws::Integer = 20000, bb_draws::Integer = 10000)
    bell = [("XX", 0), ("ZZ", 0)]
    singlet = [("XX", 1), ("ZZ", 1)]  # |Psi->: negative signs on the surface code
    ghz = [("XXX", 0), ("ZZI", 0), ("IZZ", 0)]
    rnd3 = nothing
    for i in Iterators.countfrom(0)
        rnd3 = Circuits.random_stabilizer_state(3, seed_for("handover_state", "random3", i))
        all(Circuits.exact_expectation(rnd3, join(j == q ? P : 'I' for j in 1:3)) == 0 for q in 1:3 for P in "XYZ") && break
    end
    rnd12 = vcat([("Y" * "I"^11, 1)],
                 [("I" * lab, s) for (lab, s) in Circuits.random_css_state(11, seed_for("handover_state", "random12"))])
    out = OrderedDict{String,Any}[]
    for (p, tag) in ((0.0, "p0"), (1e-3, "p1e-3"))
        for (nm, k, st) in (("bell", 2, bell), ("singlet", 2, singlet), ("ghz", 3, ghz), ("random3", 3, rnd3))
            push!(out, OrderedDict{String,Any}("name" => "surface_d3_$(nm)_$tag", "code" => "surface", "k" => k,
                "state" => st, "ensemble" => surface_ensemble(k), "p" => p, "draws" => draws,
                "code_desc" => OrderedDict{String,Any}("name" => "rotated surface code", "distance" => 3, "patches" => k,
                    "logical_basis" => "per patch: Xbar = X on the top row, Zbar = Z on the left column")))
        end
        push!(out, OrderedDict{String,Any}("name" => "bb72_randomcss12_$tag", "code" => "bb72", "k" => 12,
            "state" => rnd12, "ensemble" => bb72_ensemble(), "p" => p, "draws" => bb_draws,
            "code_desc" => OrderedDict{String,Any}("name" => "BB72 [[72,12,6]] bivariate bicycle code", "distance" => 6,
                "logical_basis" => "common CSS basis of ftoracle/bb_surgery.bb_setup: qubit 0 = " *
                                   "weight-6 (Xbar_0, Zbar_0) with overlap 1 (D-24)")))
    end
    return out
end

"""
    v1_job(job) -> (sid, labels, decoder, result or nothing)

Build one setting and, if it was drawn, run all its shots (`job = (spec, sid, nshots)`). Workers run these.
"""
function v1_job(job)
    spec, sid, nshots = job
    exp, labels, dec = v1_setting(spec, sid, spec["state"])
    nshots == 0 && return sid, labels, dec, nothing
    res = run_oracle(exp, spec["p"], nshots; decoder = dec, seed = seed_for("handover", spec["name"], spec["p"], sid),
                     setting_labels = labels)
    return sid, labels, dec, (outcomes = res.outcomes, confidence = res.confidence, kinds = res.confidence_kind, cost = res.cost)
end

"""
    generate_v1(spec, dir; pmap=map, code_version="unknown", generator=..., created=today) -> number of targets

Simulate one v1 dataset and write `<name>.meta.json`, `.records.jsonl` and `.exact.csv` into `dir`.
"""
function generate_v1(spec, dir::AbstractString; pmap = map, code_version = "unknown",
                     generator = "scripts/make_datasets.jl", created = Libc.strftime("%Y-%m-%d", time()))
    t0 = time()
    name, p, ens, draws = spec["name"], spec["p"], spec["ensemble"], spec["draws"]
    np = pyimport("numpy")
    draw_seed = seed_for("handover_draws", name, p)
    sids = sort(collect(keys(ens)))
    probs = np.array(pylist([ens[s] for s in sids]))
    draw = pyconvert(Vector{Int}, np.random.default_rng(draw_seed).choice(length(sids); size = draws, p = probs / probs.sum()))
    idx_of = Dict(sid => findall(==(si - 1), draw) .- 1 for (si, sid) in enumerate(sids))
    records = Vector{Any}(nothing, draws)
    ens_labels, costs, decoders = OrderedDict{String,Vector{String}}(), Dict{String,Any}(), OrderedDict{String,String}()
    results = pmap(v1_job, [(spec, sid, length(idx_of[sid])) for sid in sids])
    for (sid, labels, dec, res) in results
        ens_labels[sid], decoders[sid] = labels, dec
        res === nothing && continue
        costs[sid] = res.cost
        for (j, shot) in enumerate(idx_of[sid])
            records[shot+1] = OrderedDict{String,Any}("shot" => shot, "setting_id" => sid, "setting" => labels,
                "outcomes" => Int.(res.outcomes[j, :]), "confidence" => _confidence_cell(Float64(res.confidence[j])),
                "confidence_kind" => res.kinds[j], "cost" => res.cost)
        end
    end
    lam = learnable(ens_labels, ens)
    rows = [Any[P, Circuits.exact_expectation(spec["state"], P), py_round(lam[P][1], 12), join(lam[P][2], ";")]
            for P in sort(collect(keys(lam)))]
    nm = Circuits.NOISE_MODELS[Circuits.DEFAULT_NOISE]
    meta = _od(
        "format_version" => FORMAT_VERSION, "dataset" => name, "created" => created, "code" => spec["code_desc"],
        "n_logical" => spec["k"], "logical_qubit_order" => "character i of every Pauli string = logical qubit i (0-based)",
        "conventions" => _od(
            "outcome_bit" => "0 -> eigenvalue +1, 1 -> eigenvalue -1 of the listed Hermitian string",
            "Y" => "Y_i = i X_i Z_i (the Hermitian logical Y of qubit i); strings are tensor products, no extra phase",
            "setting" => "all strings of one setting commute and are measured jointly on the same shot",
            "confidence" => "complementary gap ln[P(best)/P(best with this outcome flipped)], >= 0, minimum " *
                            "over the setting's strings; meaning of each value given by confidence_kind",
            "confidence_kind" => "exact | approx | not_found | decoder_failed | unavailable | noiseless " *
                                 "(INTERFACE.md section 2)",
            "state_for_testing" => "'sign' is +1 or -1 (not a sign bit)"),
        "noise" => _od("model" => p > 0 ? Circuits.DEFAULT_NOISE : "none", "p" => p,
            "multipliers" => _od((string(k) => v for (k, v) in pairs(nm))...),
            "description" => "circuit-level (Bravyi et al. 2308.07915 channel set, D-29): reset and " *
                             "measurement flips p, DEPOLARIZE2(p) after every 2-qubit gate, " *
                             "DEPOLARIZE1(p) on idle data qubits in every layer; ideal (noiseless) " *
                             "state preparation; destructive readout flips p"),
        "rounds_surgery" => R_SURGERY, "decoders" => decoders,
        "ensemble" => _od((sid => _od("probability" => ens[sid], "strings" => ens_labels[sid],
                                      "cost" => get(costs, sid, nothing)) for sid in sids)...),
        "draws" => draws, "draw_seed" => draw_seed, "code_version" => code_version, "generator" => generator,
        "state_for_testing" => [_od("generator" => l, "sign" => s == 1 ? -1 : 1) for (l, s) in spec["state"]],
        "runtime_s" => py_round(time() - t0, 1))
    mkpath(dir)
    base = joinpath(dir, name)
    write_meta(base * ".meta.json", meta)
    write_records_raw(base * ".records.jsonl", records)
    write_csv_table(base * ".exact.csv", CSVTable(["pauli", "exact", "lambda_P", "settings"], rows, "\r\n"))
    return length(rows)
end
