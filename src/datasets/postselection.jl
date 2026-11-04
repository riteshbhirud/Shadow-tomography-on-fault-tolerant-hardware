# Postselection on the decoder gap for the noisy Nov-20 datasets, a check of our data
# For gap cuts at quantiles 0-0.5: per setting the kept fraction and
# the flip rate of stabilizer-group elements, and the fidelity estimate with lambda recomputed from the
# kept setting frequencies. Quantiles are numpy's (bridged); means and standard errors follow numpy's sums.

const POSTSELECTION_QUANTILES = (0.0, 0.05, 0.1, 0.2, 0.3, 0.5)

"""
    postselection_load(prefix) -> (meta, conf, bases, vals)

Per shot: the gap (`Inf` for "inf", `NaN` for null), the measured bases (the non-identity letter of every
setting string) and the ±1 outcome values.
"""
function postselection_load(prefix::AbstractString)
    meta = read_meta(prefix * ".meta.json")
    recs = read_records_raw(prefix * ".records.jsonl")
    conf = [r["confidence"] === nothing ? NaN : r["confidence"] == "inf" ? Inf : Float64(r["confidence"]) for r in recs]
    bases = [join(first(c for c in s if c != 'I') for s in r["setting"]) for r in recs]
    vals = permutedims(reduce(hcat, [[1 - 2 * b for b in r["outcomes"]] for r in recs]))
    return meta, conf, bases, vals
end

_matches(g, S) = all(g[i] == 'I' || S[i] == g[i] for i in eachindex(g))

"""
    postselection_fidelity(bases, vals, group, m) -> (F, se)

Reference estimator: per shot 2^-m (1 + sum_g sign_g [match] prod / lambda_g), with lambda the matching
fraction of these shots. (NaN, NaN) if some element is never measured.
"""
function postselection_fidelity(bases, vals, group, m)
    n = length(bases)
    v = ones(n)
    for (g, s) in group
        match = [_matches(g, b) for b in bases]
        lam = count(match) / n
        lam == 0 && return NaN, NaN
        idx = [i for i in 1:m if g[i] != 'I']
        pr = [prod(vals[k, idx]) for k in 1:n]
        v = v .+ (s .* match .* pr) ./ lam
    end
    v = v ./ 2^m
    return np_mean(v), np_std(v; ddof = 1) / sqrt(n)
end

"""
    postselection_analyse(prefix) -> (meta, element rows, fidelity rows)

The two tables of one dataset; element rows only for datasets without logical noise.
"""
function postselection_analyse(prefix::AbstractString)
    meta, conf, bases, vals = postselection_load(prefix)
    m = meta["n_logical"]
    ka = meta["known_answers"]
    group = [(g["pauli"], g["sign"]) for g in ka["stabilizer_group"]]
    finite = filter(isfinite, conf)
    np = pyimport("numpy")
    cuts = isempty(finite) ? [0.0] :
           sort!(unique([pyconvert(Float64, np.quantile(np.array(finite), q)) for q in POSTSELECTION_QUANTILES]))
    rows, frows = OrderedDict{String,Any}[], OrderedDict{String,Any}[]
    for cut in cuts
        keep = conf .>= cut  # NaN -> false, Inf -> true
        F, se = postselection_fidelity(bases[keep], vals[keep, :], group, m)
        push!(frows, OrderedDict{String,Any}("cut" => cut, "fraction_kept" => count(keep) / length(keep),
            "shots_kept" => count(keep), "F" => F, "se" => se, "F_expected_all_shots" => ka["fidelity_expected_noisy"],
            "F_with_logical_noise" => ka["fidelity_with_logical_noise"]))
        get(meta, "logical_noise", nothing) in (nothing, false) || isempty(meta["logical_noise"]) || continue
        for S in sort!(unique(bases))
            sel = bases .== S
            k = sel .& keep
            for (g, s) in group
                _matches(g, S) || continue
                idx = [i for i in 1:m if g[i] != 'I']
                flips = [prod(vals[j, idx]) != s for j in findall(k)]
                push!(rows, OrderedDict{String,Any}("cut" => cut, "setting" => "local:$S",
                    "element" => (s > 0 ? "+" : "-") * g, "fraction_kept" => count(k) / count(sel),
                    "shots_kept" => count(k), "flips" => count(flips),
                    "q_kept" => isempty(flips) ? NaN : count(flips) / length(flips)))
            end
        end
    end
    return meta, rows, frows
end
