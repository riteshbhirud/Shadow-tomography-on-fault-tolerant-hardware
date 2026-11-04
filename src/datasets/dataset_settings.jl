# Flip rates of exactly the settings in the v1.0/v1.1 handover datasets, per string and per subset
# product. Products need their own q because flips within one
# setting are correlated. Inputs are +1 eigenstates of every measured string.

const DATASET_SETTINGS_P, DATASET_SETTINGS_R = 1e-3, 3

# Letterwise product of two Pauli labels, phases dropped.
_mul_nophase(a, b) = join(x == 'I' ? y : y == 'I' ? x : x == y ? 'I' : pauli_mul(string(x), string(y))[2][1]
                          for (x, y) in zip(a, b))

"""
    dataset_setting_experiment(code, sid, k) -> (experiment, labels, decoder)

The setting `sid` on `k` logical qubits, on a +1 eigenstate of every string it measures: surface code
d = 3 with r = 3 surgery rounds (local readout or ZZ/XX chain), or BB72 readout / gauging of Ybar_0.
"""
function dataset_setting_experiment(code::AbstractString, sid::AbstractString, k::Integer)
    R = DATASET_SETTINGS_R
    if code == "surface"
        kind, arg = split(sid, ':')
        if kind == "local"
            state = [(join(j == i ? arg[i+1] : 'I' for j in 0:k-1), 0) for i in 0:k-1]
        else
            Pc = arg[1]
            other = Pc == 'Z' ? 'X' : 'Z'
            state = [(join(j in (i, i + 1) ? Pc : 'I' for j in 0:k-1), 0) for i in 0:k-2]
            push!(state, (string(other)^k, 0))
        end
        exp, labels = Circuits.setting_experiment(sid, 3, k, state, R)
        return exp, labels, "pymatching"
    end
    if startswith(sid, "readout:")
        b = sid[end]
        state = [(join(j == i ? b : 'I' for j in 0:11), 0) for i in 0:11]
        exp, labels = Circuits.bb72_readout_experiment(b, state)
        return exp, labels, "relay"
    end
    return Circuits.bb_gauging('Y', R; state_kind = "+1"), ["Y" * "I"^11], "relay"
end

"""
    dataset_settings_job(job) -> rows

`job = (code, sid, k, shots)`. One row per string and per subset product (BB72 readout: strings and
pairs only), then the rate at which any string flips.
"""
function dataset_settings_job(job)
    code, sid, k, shots = job
    exp, labels, dec = dataset_setting_experiment(code, sid, k)
    res = run_oracle(exp, DATASET_SETTINGS_P, shots; decoder = dec, seed = seed_for("dataset_settings", code, sid, k),
                     with_gap = false, setting_labels = labels)
    flips = Int.(res.outcomes)  # eigenvalue +1 for every string -> outcome bit = flip
    n = size(flips, 1)
    m = length(labels)
    # itertools.combinations order, sizes 1..m
    subsets = [c for r in 1:m for c in _combinations(0:m-1, r)]
    m > 4 && filter!(s -> length(s) <= 2, subsets)  # BB72 readout: single strings, pairs and 'any'
    r = startswith(sid, "surgery") || startswith(sid, "gauge") ? DATASET_SETTINGS_R : 0
    rows = OrderedDict{String,Any}[]
    for sub in subsets
        lab = labels[sub[1]+1]
        for j in sub[2:end]
            lab = _mul_nophase(lab, labels[j+1])
        end
        f = [reduce(xor, flips[s, j+1] for j in sub) for s in 1:n]
        nf = count(==(1), f)
        lo, hi = wilson(nf, n)
        qs = [1 - 2 * (count(==(1), view(flips, :, j + 1)) / n) for j in sub]
        push!(rows, OrderedDict{String,Any}("code" => code, "setting_id" => sid, "n_logical" => k, "r" => r,
            "noise" => Circuits.DEFAULT_NOISE, "p" => DATASET_SETTINGS_P, "target" => lab,
            "kind" => length(sub) == 1 ? "string" : "product", "shots" => n, "flips" => nf, "q" => nf / n,
            "ci_lo" => lo, "ci_hi" => hi, "q_if_independent" => (1 - foldl(*, qs)) / 2))
    end
    na = count(s -> any(==(1), view(flips, s, :)), 1:n)
    lo, hi = wilson(na, n)
    push!(rows, OrderedDict{String,Any}("code" => code, "setting_id" => sid, "n_logical" => k, "r" => rows[1]["r"],
        "noise" => Circuits.DEFAULT_NOISE, "p" => DATASET_SETTINGS_P, "target" => "any", "kind" => "any_string_flipped",
        "shots" => n, "flips" => na, "q" => na / n, "ci_lo" => lo, "ci_hi" => hi, "q_if_independent" => NaN))
    return rows
end

# itertools.combinations(xs, r) in its (lexicographic) order
function _combinations(xs, r)
    xs = collect(xs)
    r == 0 && return [Int[]]
    out = Vector{Int}[]
    for (i, x) in enumerate(xs), rest in _combinations(xs[i+1:end], r - 1)
        push!(out, [x; rest])
    end
    return out
end

"""The jobs of dataset_settings_q.csv: surface k = 2, 3 (local X/Z, ZZ and XX chains), BB72 readouts and Y0."""
function dataset_settings_jobs()
    jobs = Tuple[]
    for k in (2, 3)
        for b in product_strings("XZ", k)
            push!(jobs, ("surface", "local:" * b, k, 100_000))
        end
        push!(jobs, ("surface", "surgery:ZZ", k, 400_000), ("surface", "surgery:XX", k, 400_000))
    end
    push!(jobs, ("bb72", "readout:Z", 12, 20_000), ("bb72", "readout:X", 12, 20_000), ("bb72", "gauge:Y0", 12, 20_000))
    return jobs
end
