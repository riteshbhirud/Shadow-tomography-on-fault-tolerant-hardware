# Basis (in)dependence of the native-ensemble shadow statistics.
#
# A CSS change of logical basis is a bijection of Pauli coordinates that keeps the Z-type and X-type sets, so the
# λ spectrum, orbit sizes, learnable fraction and the pure-X / pure-Z / mixed split are invariant, while
# per-weight fractions are not. Every logical action is recomputed from the physical generators in each basis.

const BASIS_MAXW = 4
const CSS_TYPE_ENSEMBLES = ["readout_only", "shifts", "ibm", "aut"]

"""
    basis_labels(name; n_random=5, bravyi_seeds=1:8)

Labels of the bases analysed for a code.
"""
function basis_labels(name::AbstractString; n_random::Integer=5, bravyi_seeds=1:8)
    labels = vcat(["named", "generic"], ["random_css_$s" for s in 0:n_random-1])
    name == "BB72" && append!(labels, ["bravyi_seed_$s" for s in bravyi_seeds])
    return labels
end

"""
    random_css_basis(code0, seed)

X̄′ = A X̄ + R_x H_X, Z̄′ = A^{-T} Z̄ + R_z H_Z, with the random A, R_x, R_z for `seed` read from `rng_inputs.json`.
"""
function random_css_basis(code0::StabilizerCode, seed::Integer)
    d = rng_inputs()["random_css"][code0.name][string(seed)]
    a, rx, rz = frombits(d["A"]), frombits(d["Rx"]), frombits(d["Rz"])
    k, n = logical_k(code0), code0.n
    hx, hz = code0.meta[:hx], code0.meta[:hz]
    xl = gf2mul(a, code0.lx[:, 1:n]) .⊻ gf2mul(rx, hx)
    zl = gf2mul(permutedims(gf2inv(a)), code0.lz[:, n+1:end]) .⊻ gf2mul(rz, hz)
    z = zeros(Bool, k, n)
    meta = copy(code0.meta)
    meta[:basis] = "random CSS basis seed $seed"
    delete!(meta, :_stab_row_basis)
    c = StabilizerCode("$(code0.name)[random_css_$seed]", n, code0.stabs, hcat(xl, z), hcat(z, zl), meta)
    validate(c)
    return c
end

"""
    make_basis(name, label, code0)

The code in the basis named by `label` (`named`, `generic`, `random_css_<s>`, `bravyi_seed_<s>`), or
`nothing` when the Bravyi-style search for that seed finds no complete basis.
"""
function make_basis(name::AbstractString, label::AbstractString, code0::StabilizerCode)
    label == "named" && return code0
    label == "generic" && return make_bb_code(name, "generic")
    s = parse(Int, split(label, "_")[end])
    startswith(label, "random_css_") && return random_css_basis(code0, s)
    if startswith(label, "bravyi_seed_")
        l, m, amon, bmon = BB_PARAMS[name]
        _, _, _, hx, hz = bb_matrices(l, m, amon, bmon)
        res = bb72_fgh_basis(s)
        res === nothing && return nothing
        xl, zl, info = res
        meta = filter(kv -> !startswith(String(kv[1]), "_"), copy(code0.meta))
        delete!(meta, :hx); delete!(meta, :hz)
        meta[:basis] = "Bravyi-style search seed $s"
        meta[:fgh_info] = info
        return make_css_code("$name[bravyi_seed_$s]", hx, hz, xl, zl; meta...)
    end
    throw(ArgumentError(label))
end

"""
    exact_spectrum(lab, k)

Exact λ spectra (identity excluded) per readout as reduced fractions with multiplicities, plus the orbit-size
histogram.
"""
function exact_spectrum(lab::AbstractVector{UInt32}, k::Integer)
    xmask, zmask = pauli_masks(k)
    reps = sort(unique(lab))
    pos = zeros(Int, length(lab))
    for (i, r) in enumerate(reps)
        pos[r+1] = i
    end
    size_ = zeros(Int, length(reps)); nz = zeros(Int, length(reps)); nx = zeros(Int, length(reps))
    for (i, l) in enumerate(lab)
        v = UInt32(i - 1)
        j = pos[l+1]
        size_[j] += 1
        v != 0 && v & xmask == 0 && (nz[j] += 1)
        v != 0 && v & zmask == 0 && (nx[j] += 1)
    end
    keep = reps .!= lab[1]
    @assert size_[.!keep] == [1]
    size_, nz, nx = size_[keep], nz[keep], nx[keep]
    out = OrderedDict{String,Any}()
    for (ro, num, den) in (("Z", nz, size_), ("X", nx, size_), ("mix", nz .+ nx, 2 .* size_))
        mult = OrderedDict{Tuple{Int,Int},Int}()
        for (a, b, s) in zip(num, den, size_)
            g = gcd(a, b)
            key = a == 0 ? (0, 1) : (a ÷ g, b ÷ g)
            mult[key] = get(mult, key, 0) + s
        end
        out[ro] = [[a, b, mult[(a, b)]] for (a, b) in sort(collect(keys(mult)))]
    end
    hist = OrderedDict{Int,Int}()
    for s in size_
        hist[s] = get(hist, s, 0) + 1
    end
    out["orbit_size_hist"] = [[s, hist[s]] for s in sort(collect(keys(hist)))]
    return out
end

"""
    analyse_basis(code, kept, k, wts)

Spectrum, learnable counts, type split and per-weight table of every ensemble in this basis. `kept` maps
ensemble name to the kept physical `Element`s.
"""
function analyse_basis(code::StabilizerCode, kept::AbstractDict, k::Integer, wts::AbstractVector; timings::Bool=false)
    xmask, zmask = pauli_masks(k)
    nv = 1 << 2k
    has_x = [UInt32(i) & xmask != 0 for i in 0:nv-1]
    has_z = [UInt32(i) & zmask != 0 for i in 0:nv-1]
    classes = ["pure_X" => has_x .& .!has_z, "pure_Z" => has_z .& .!has_x, "mixed" => has_x .& has_z]
    pure = classes[1][2] .| classes[2][2]
    res = OrderedDict{String,Any}("basis" => code.meta[:basis],
        "basis_weights" => OrderedDict("X" => vec(sum(code.lx; dims=2)), "Z" => vec(sum(code.lz; dims=2))),
        "ensembles" => OrderedDict{String,Any}())
    for ename in ENSEMBLES
        t = time()
        gens = [logical_action(code, e.phys) for e in kept[ename]]
        lab = isempty(gens) ? collect(UInt32(0):UInt32(nv - 1)) : orbits(gens, 2k)
        spec = exact_spectrum(lab, k)
        hit = falses(nv)   # λ_mix > 0 iff the orbit contains a pure-X or pure-Z vector
        for i in 1:nv
            pure[i] && (hit[lab[i]+1] = true)
        end
        learn = [hit[lab[i]+1] for i in 1:nv]
        learn[1] = false
        ts = OrderedDict{String,Any}(c => OrderedDict{String,Any}("n" => count(m), "n_learnable" => count(learn .& m)) for (c, m) in classes)
        r = OrderedDict{String,Any}(
            "n_nonidentity" => nv - 1,
            "n_learnable_mix" => count(learn),
            "frac_learnable_mix" => count(learn) / (nv - 1),
            "n_orbits_nonidentity" => length(unique(lab)) - 1,
            "spectrum" => spec,
            "type_split" => ts,
            "by_weight" => [OrderedDict{String,Any}("weight" => w, "n_paulis" => count(==(w), wts),
                    "n_learnable" => count(learn .& (wts .== w)), "frac_learnable" => np_mean(learn[wts.==w]))
                for w in 1:BASIS_MAXW],
        )
        r["learnable_set_is_exactly_pureX_or_pureZ"] = ts["mixed"]["n_learnable"] == 0 &&
            ts["pure_X"]["n_learnable"] == ts["pure_X"]["n"] && ts["pure_Z"]["n_learnable"] == ts["pure_Z"]["n"]
        timings && (r["seconds"] = py_round(time() - t, 1))
        res["ensembles"][ename] = r
    end
    return res
end

"""
    kept_elements(code0, ens, generators)

For every ensemble, the first physical `Element` whose logical action equals each kept generator.
"""
function kept_elements(ens::AbstractDict, generators::AbstractDict)
    out = OrderedDict{String,Vector{Element}}()
    for ename in ENSEMBLES
        out[ename] = map(generators[ename]) do km
            i = findfirst(e -> e.logical == km, ens[ename])
            i === nothing && error("$ename: kept generator not found among the physical elements")
            ens[ename][i]
        end
    end
    return out
end

"""
    basis_independence(name, code0, kept; labels=basis_labels(name), timings=false)

Per-basis records (label => record, or a `failed` record) for every basis of one code. `timings=true` adds
`seconds` per ensemble and `seconds_total` per basis.
"""
function basis_independence(name::AbstractString, code0::StabilizerCode, kept::AbstractDict;
                            labels=basis_labels(name), timings::Bool=false)
    k = logical_k(code0)
    wts = logical_weights(k)
    out = OrderedDict{String,Any}()
    for lb in labels
        code = make_basis(name, lb, code0)
        if code === nothing
            s = split(lb, "_")[end]
            out[lb] = OrderedDict{String,Any}("label" => lb, "failed" => "Bravyi-style search found no complete basis for seed $s")
            continue
        end
        t = time()
        r = analyse_basis(code, kept, k, wts; timings)
        r["label"] = lb
        timings && (r["seconds_total"] = py_round(time() - t, 1))
        out[lb] = r
    end
    return out
end

"""
    basis_record(name, kept_generators, label)

The `<code>_basis_<label>.json` record of one basis, rebuilt from the kept logical generators. Needs no group
orders, so it runs on a worker without Oscar.
"""
function basis_record(name::AbstractString, kept_generators::AbstractDict, label::AbstractString)
    code0 = make_bb_code(name, "fgh")
    ens, _, _ = build_ensembles(code0)
    kept = kept_elements(ens, kept_generators)
    return basis_independence(name, code0, kept; labels=[label], timings=true)[label]
end

"""
    aggregate_bases(name, data, labels) -> (summary, per_weight_rows, spread_rows, spectrum_rows)

The `<code>_summary.json` record and the rows of the three CSV tables, from the per-basis records.
"""
function aggregate_bases(name::AbstractString, data::AbstractDict, labels)
    failed = OrderedDict{String,Any}(lb => data[lb]["failed"] for lb in labels if haskey(data[lb], "failed"))
    labels = [lb for lb in labels if !haskey(failed, lb)]
    rand = [lb for lb in labels if startswith(lb, "random_css_")]
    brav = [lb for lb in labels if startswith(lb, "bravyi_seed_")]
    named = data["named"]
    summ = OrderedDict{String,Any}("code" => name, "named_basis" => named["basis"], "bases" => labels,
        "failed_basis_searches" => failed,
        "basis_descriptions" => OrderedDict{String,Any}(lb => data[lb]["basis"] for lb in labels),
        "basis_rep_weights" => OrderedDict{String,Any}(lb => data[lb]["basis_weights"] for lb in labels),
        "ensembles" => OrderedDict{String,Any}())
    rows, spread_rows, spec_rows = OrderedDict{String,Any}[], OrderedDict{String,Any}[], OrderedDict{String,Any}[]
    for ename in ENSEMBLES
        e0 = named["ensembles"][ename]
        inv_keys = ("n_learnable_mix", "spectrum", "type_split", "n_orbits_nonidentity")
        identical = OrderedDict{String,Any}(key => all(data[lb]["ensembles"][ename][key] == e0[key] for lb in labels)
                                     for key in inv_keys)
        hist = e0["spectrum"]["orbit_size_hist"]
        es = OrderedDict{String,Any}(
            "basis_independent" => OrderedDict{String,Any}(
                "n_nonidentity" => e0["n_nonidentity"], "n_learnable_mix" => e0["n_learnable_mix"],
                "frac_learnable_mix" => e0["frac_learnable_mix"], "n_orbits_nonidentity" => e0["n_orbits_nonidentity"],
                "type_split" => e0["type_split"],
                "spectrum_mix" => [OrderedDict{String,Any}("lam" => "$a/$b", "lam_float" => a / b,
                    "norm2" => a != 0 ? b / a : nothing, "multiplicity" => c) for (a, b, c) in e0["spectrum"]["mix"]],
                "orbit_size_hist" => length(hist) < 200 ? hist : "$(length(hist)) distinct sizes (see per-basis JSON)"),
            "identical_in_all_bases" => identical,
            "all_invariants_identical" => all(values(identical)),
            "css_claim_holds_in_basis" => OrderedDict{String,Any}(
                lb => data[lb]["ensembles"][ename]["learnable_set_is_exactly_pureX_or_pureZ"] for lb in labels),
            "per_weight" => OrderedDict{String,Any}[],
        )
        for wi in 1:BASIS_MAXW
            val = OrderedDict(lb => data[lb]["ensembles"][ename]["by_weight"][wi] for lb in labels)
            frac(lbs, f) = f(val[lb]["frac_learnable"] for lb in lbs)
            row = OrderedDict{String,Any}("ensemble" => ename, "weight" => wi, "n_paulis" => val["named"]["n_paulis"],
                "named_n_learnable" => val["named"]["n_learnable"], "named_frac" => val["named"]["frac_learnable"],
                "generic_frac" => haskey(val, "generic") ? val["generic"]["frac_learnable"] : "",
                "random_css_min" => frac(rand, minimum), "random_css_max" => frac(rand, maximum),
                "n_random_css" => length(rand))
            if !isempty(brav)
                row["bravyi_style_min"] = frac(brav, minimum)
                row["bravyi_style_max"] = frac(brav, maximum)
                row["n_bravyi_style"] = length(brav)
            end
            push!(es["per_weight"], row)
            push!(spread_rows, row)
            for lb in labels
                push!(rows, OrderedDict{String,Any}("code" => name, "ensemble" => ename, "basis" => lb, "weight" => wi,
                    "n_paulis" => val[lb]["n_paulis"], "n_learnable" => val[lb]["n_learnable"],
                    "frac_learnable" => val[lb]["frac_learnable"]))
            end
        end
        for ro in ("Z", "X", "mix"), (a, b, c) in e0["spectrum"][ro]
            push!(spec_rows, OrderedDict{String,Any}("code" => name, "ensemble" => ename, "readout" => ro,
                "lam" => "$a/$b", "lam_float" => a / b, "norm2" => a != 0 ? b / a : "inf",
                "multiplicity" => c, "identical_in_all_bases" => identical["spectrum"]))
        end
        summ["ensembles"][ename] = es
    end
    summ["css_claim_check"] = OrderedDict{String,Any}(ename => OrderedDict{String,Any}(
        "holds_in_all_bases" => all(values(summ["ensembles"][ename]["css_claim_holds_in_basis"])),
        "n_bases" => length(labels)) for ename in ENSEMBLES)
    summ["all_invariants_identical_all_ensembles"] = all(summ["ensembles"][e]["all_invariants_identical"] for e in ENSEMBLES)
    return summ, rows, spread_rows, spec_rows
end
