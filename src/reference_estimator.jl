# Reference inverse-channel estimator used only to test the oracle end to end
# (not the project's estimator).
#
# For a target P, lambda_P = Pr[P = ±(product of a subset of the measured strings)], and the one-shot
# estimate 1[P measurable] * sign * prod (-1)^b / lambda_P is unbiased when lambda_P > 0.

# (a, b) -> (k, c) with a * b = i^k c
const _PAULI_MUL = Dict{Tuple{Char,Char},Tuple{Int,Char}}(
    ('I', 'I') => (0, 'I'), ('I', 'X') => (0, 'X'), ('I', 'Y') => (0, 'Y'), ('I', 'Z') => (0, 'Z'),
    ('X', 'I') => (0, 'X'), ('Y', 'I') => (0, 'Y'), ('Z', 'I') => (0, 'Z'),
    ('X', 'X') => (0, 'I'), ('Y', 'Y') => (0, 'I'), ('Z', 'Z') => (0, 'I'),
    ('X', 'Y') => (1, 'Z'), ('Y', 'X') => (3, 'Z'), ('Y', 'Z') => (1, 'X'),
    ('Z', 'Y') => (3, 'X'), ('Z', 'X') => (1, 'Y'), ('X', 'Z') => (3, 'Y'))

"""
    pauli_mul(a, b) -> (k, label)

Product of two Pauli labels with a * b = i^k label.
"""
function pauli_mul(a::AbstractString, b::AbstractString)
    ph = 0
    out = Char[]
    for (x, y) in zip(a, b)
        p, r = _PAULI_MUL[(x, y)]
        ph += p
        push!(out, r)
    end
    return ph % 4, String(out)
end

"""
    decompose(P, S) -> (subset, sign) or nothing

0-based `subset` of the labels `S` and `sign` with P = sign * prod of S over the subset, found by brute
force over the subsets in increasing bit-mask order (as the reference, so the same subset is returned).
"""
function decompose(P::AbstractString, S::AbstractVector)
    m = length(S)
    for mask in 1:(1<<m)-1
        ph, cur = 0, "I"^length(P)
        for i in 0:m-1
            if (mask >> i) & 1 == 1
                p, cur = pauli_mul(cur, S[i+1])
                ph += p
            end
        end
        if cur == P
            ph %= 4
            ph in (0, 2) || error("product of commuting Hermitian Paulis must be Hermitian")
            return [i for i in 0:m-1 if (mask >> i) & 1 == 1], ph == 0 ? 1 : -1
        end
    end
    return nothing
end

"""
    lambdas(targets, ensemble, setting_labels) -> OrderedDict

lambda_P for each target, from the ensemble `{setting => probability}` (an ordered dict; the sum runs in
its order) and `{setting => measured Pauli labels}`.
"""
function lambdas(targets, ensemble::AbstractDict, setting_labels::AbstractDict)
    lam = OrderedDict{String,Any}()
    for P in targets
        # Python's sum starts from the integer 0, so a target nothing measures gets 0, not 0.0
        lam[P] = foldl(+, (pr for (s, pr) in ensemble if decompose(P, setting_labels[s]) !== nothing); init = 0)
    end
    return lam
end

"""
    estimate(targets, records, lam) -> OrderedDict(P => (mean, standard error, n))

From records with keys "setting" (or "setting_labels") and "outcomes". A target with lambda_P = 0 gives
(NaN, NaN, n). Mean and standard error follow numpy's summation order.
"""
function estimate(targets, records, lam::AbstractDict)
    out = OrderedDict{String,Tuple{Float64,Float64,Int}}()
    n = length(records)
    for P in targets
        if lam[P] == 0
            out[P] = (NaN, NaN, n)
            continue
        end
        vals = zeros(n)
        for (i, rec) in enumerate(records)
            dec = decompose(P, haskey(rec, "setting") ? rec["setting"] : rec["setting_labels"])
            if dec !== nothing
                sub, sign = dec
                vals[i] = sign * prod(Int[1 - 2 * rec["outcomes"][j+1] for j in sub]) / lam[P]
            end
        end
        out[P] = (np_mean(vals), np_std(vals; ddof = 1) / sqrt(n), n)
    end
    return out
end
