# The inverse-channel estimator and its use on oracle records. The estimator is only a test
# helper, so it lives here and not in src/.

@testmodule RefEstimator begin
    using Statistics

    # (a, b) => (k, c) with a * b = i^k c
    const PAULI_MUL = Dict(
        ('I', 'I') => (0, 'I'), ('I', 'X') => (0, 'X'), ('I', 'Y') => (0, 'Y'), ('I', 'Z') => (0, 'Z'),
        ('X', 'I') => (0, 'X'), ('Y', 'I') => (0, 'Y'), ('Z', 'I') => (0, 'Z'),
        ('X', 'X') => (0, 'I'), ('Y', 'Y') => (0, 'I'), ('Z', 'Z') => (0, 'I'),
        ('X', 'Y') => (1, 'Z'), ('Y', 'X') => (3, 'Z'), ('Y', 'Z') => (1, 'X'),
        ('Z', 'Y') => (3, 'X'), ('Z', 'X') => (1, 'Y'), ('X', 'Z') => (3, 'Y'))

    """Product of two Pauli labels as (k, label) with a * b = i^k label."""
    function mul(a::AbstractString, b::AbstractString)
        ph, out = 0, Char[]
        for (x, y) in zip(a, b)
            p, r = PAULI_MUL[(x, y)]
            ph += p
            push!(out, r)
        end
        return mod(ph, 4), String(out)
    end

    """(subset, sign) with P = sign * prod of S over the subset (0-based indices), or nothing; brute force."""
    function decompose(P::AbstractString, S)
        m = length(S)
        for mask in 1:(1<<m)-1
            ph, cur = 0, "I"^length(P)
            for i in 0:m-1
                if (mask >> i) & 1 == 1
                    p, cur = mul(cur, S[i+1])
                    ph += p
                end
            end
            if cur == P
                ph = mod(ph, 4)
                @assert ph in (0, 2) "product of commuting Hermitian Paulis must be Hermitian"
                return [i for i in 0:m-1 if (mask >> i) & 1 == 1], ph == 0 ? 1 : -1
            end
        end
        return nothing
    end

    """λ_P for each target, from {setting => probability} and {setting => measured Pauli labels}."""
    lambdas(targets, ensemble, setting_labels) =
        Dict(P => sum((pr for (s, pr) in ensemble if decompose(P, setting_labels[s]) !== nothing); init=0.0)
             for P in targets)

    """{P => (mean, standard error, n)} from records with keys "setting" (or "setting_labels") and "outcomes"."""
    function estimate(targets, records, lam)
        out = Dict{String,Tuple{Float64,Float64,Int}}()
        for P in targets
            if lam[P] == 0
                out[P] = (NaN, NaN, length(records))
                continue
            end
            vals = zeros(length(records))
            for (i, rec) in enumerate(records)
                dec = decompose(P, haskey(rec, "setting") ? rec["setting"] : rec["setting_labels"])
                if dec !== nothing
                    sub, sign = dec
                    vals[i] = sign * prod((1 - 2 * rec["outcomes"][j+1] for j in sub); init=1) / lam[P]
                end
            end
            out[P] = (mean(vals), std(vals) / sqrt(length(vals)), length(vals))
        end
        return out
    end
end

@testitem "inverse-channel estimator: Pauli products and decompositions" tags = [:estimator] setup = [RefEstimator] begin
    using .RefEstimator: mul, decompose
    # Pauli multiplication and decomposition
    @test mul("X", "Z") == (3, "Y")  # XZ = -iY
    @test decompose("ZIZ", ["ZZI", "IZZ"]) == ([0, 1], 1)
    @test decompose("YY", ["XX", "ZZ"]) == ([0, 1], -1)  # XX*ZZ = -YY
    @test decompose("XI", ["ZZ"]) === nothing
end

@testitem "inverse-channel estimator is unbiased on exact samples" tags = [:estimator, :bridge] setup = [RefEstimator] begin
    using PythonCall
    using .RefEstimator: lambdas, estimate
    # the estimator is unbiased on exact samples: Bell state, local X/Z ensemble, outcomes drawn classically
    # from a fixed numpy stream
    np = pyimport("numpy")
    rng = np.random.default_rng(0)
    names = ["XX", "XZ", "ZX", "ZZ"]
    ens = Dict(s => 0.25 for s in names)
    labels = Dict(s => [s[1] * "I", "I" * s[2]] for s in names)
    pynames = pylist(names)
    recs = Dict{String,Any}[]
    for _ in 1:40000
        s = pyconvert(String, pystr(rng.choice(pynames)))
        if s in ("XX", "ZZ")  # perfectly correlated outcomes
            b = pyconvert(Int, rng.integers(2))
            out = [b, b]
        else
            out = pyconvert(Vector{Int}, rng.integers(2; size=2))
        end
        push!(recs, Dict("setting_labels" => labels[s], "outcomes" => out))
    end
    T = ["XX", "ZZ", "XZ", "XI", "YY"]
    lam = lambdas(T, ens, labels)
    @test lam["XX"] == 0.25 && lam["XI"] == 0.5 && lam["YY"] == 0
    est = estimate(T, recs, lam)
    for (P, ex) in (("XX", 1), ("ZZ", 1), ("XZ", 0), ("XI", 0))
        m, se, _ = est[P]
        @test abs(m - ex) < 4se
    end
    # variance check: 1/lambda - <P>^2
    m, se, n = est["XX"]
    @test abs(se^2 * n - (1 / 0.25 - 1)) < 0.3
    @test isnan(est["YY"][1])
end

@testitem "oracle records feed the inverse-channel estimator" tags = [:estimator, :oracle] setup = [RefEstimator] begin
    using FTShadows.Circuits: setting_experiment
    # oracle records feed the estimator: the documented record
    # (key "setting") goes straight into the estimator
    exp, labs = setting_experiment("local:ZZ", 3, 2, [("ZZ", 0), ("XX", 0)], 0)
    res = run_oracle(exp, 0.0, 200; seed=1, setting_labels=labs)
    est = RefEstimator.estimate(["ZZ"], records(res), Dict("ZZ" => 1.0))
    @test est["ZZ"][1] ≈ 1.0
end
