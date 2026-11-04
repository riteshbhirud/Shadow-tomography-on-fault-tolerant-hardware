# Greedy cover of target logical Paulis by native settings.
#
# A native setting (g ∈ G, readout Z or X) measures every P with g(P) Z-type (X-type). The settings were drawn
# from a seeded numpy generator; those draws are replayed from test/reference/taskb/designed_cover_draws.json.

"""
    CoverDraws

The numpy draws of one `designed_cover.run` call, consumed in order: local-Pauli settings (12 trits each) and
product-replacement steps (pair i ≠ j of pool slots and a coin).
"""
mutable struct CoverDraws
    trits::Vector{UInt8}
    walk::Vector{UInt8}
    pool::Int
    tpos::Int
    wpos::Int
end

"""
    cover_draws(name)

The recorded draws of the designed-cover run for `name`.
"""
function cover_draws(name::AbstractString)
    d = JSON.parsefile(joinpath(TASKB_REFDIR, "designed_cover_draws.json"))[name]
    trits = UInt8[]
    for byte in base64decode(d["local_trits_b64"])
        for t in 0:4
            push!(trits, (byte ÷ UInt8(3)^t) % 3)
        end
    end
    return CoverDraws(trits, base64decode(d["walk_b64"]), d["pool"], 0, 0)
end

function _next_setting!(r::CoverDraws, k)
    s = r.trits[r.tpos+1:r.tpos+k]
    r.tpos += k
    return s
end

function _next_step!(r::CoverDraws)
    r.wpos += 1
    b = Int(r.walk[r.wpos])
    ij, coin = divrem(b, 2)
    i, j = divrem(ij, r.pool)
    return i + 1, j + 1, coin == 1
end

"""
    random_elements(gens, count, draws; burn=200)

`count` group elements from the recorded product-replacement walk, after `burn` steps.
"""
function random_elements(gens, count::Integer, r::CoverDraws; burn::Integer=200)
    pool = [copy(g) for g in gens]
    while length(pool) < 10
        push!(pool, copy(gens[length(pool)%length(gens)+1]))
    end
    length(pool) == r.pool || error("pool size differs from the recorded run")
    acc = GF2Mat(I(size(gens[1], 1)))
    out = GF2Mat[]
    for t in 0:burn+count-1
        i, j, coin = _next_step!(r)
        pool[i] = coin ? gf2mul(pool[i], pool[j]) : gf2mul(pool[j], pool[i])
        acc = gf2mul(acc, pool[i])
        t >= burn && push!(out, copy(acc))
    end
    return out
end

function _image_of(targets::Vector{UInt32}, mat::AbstractMatrix)
    rows = UInt32.(pack_rows(mat))
    out = zeros(UInt32, length(targets))
    for j in 1:size(mat, 1)
        out .⊻= ((targets .>> (j - 1)) .& UInt32(1)) .* rows[j]
    end
    return out
end

"""
    greedy_cover(sets, n_targets) -> (chosen, n_covered)

Greedy set cover (ties go to the first set, as `np.argmax`). Sets hold 1-based target indices.
"""
function greedy_cover(sets, n_targets::Integer)
    covered = falses(n_targets)
    chosen = Int[]
    while !all(covered)
        gains = [count(i -> !covered[i], c) for c in sets]
        best = argmax(gains)
        gains[best] == 0 && break
        push!(chosen, best)
        covered[sets[best]] .= true
    end
    return chosen, count(covered)
end

"""
    local_pauli_cover(targets, k, draws; pool=4000)

Size of the greedy cover of `targets` by recorded random local settings s ∈ {X,Y,Z}^k (P is covered iff
P_q ∈ {I, s_q} on every qubit).
"""
function local_pauli_cover(targets::Vector{UInt32}, k::Integer, r::CoverDraws; pool::Integer=4000)
    xmask, _ = pauli_masks(k)
    tx = targets .& xmask
    tz = targets .>> k
    support = tx .| tz
    sets = Vector{Int}[]
    for _ in 1:pool
        s = _next_setting!(r, k)    # 0 = X, 1 = Y, 2 = Z
        sx = reduce(|, (UInt32(1) << (q - 1) for q in 1:k if s[q] in (0, 1)); init=UInt32(0))
        sz = reduce(|, (UInt32(1) << (q - 1) for q in 1:k if s[q] in (1, 2)); init=UInt32(0))
        push!(sets, findall(((tx .⊻ (support .& sx)) .== 0) .& ((tz .⊻ (support .& sz)) .== 0)))
    end
    chosen, _ = greedy_cover(sets, length(targets))
    return length(chosen)
end

"""
    designed_cover(name, generators, lambda_mix, draws; k=12, wmax=2, pool=6000, seed=7)

The `designed_cover_<code>.json` record: greedy native and local-Pauli covers of the learnable Paulis of weight
≤ `wmax` for every ensemble. `generators` are the kept generators (in the saved order), `lambda_mix` the
`Float32` λ_mix tables; `seed` is only recorded (the draws come from `draws`).
"""
function designed_cover(name, generators::AbstractDict, lambda_mix::AbstractDict, r::CoverDraws;
                        k::Integer=12, wmax::Integer=2, pool::Integer=6000, seed::Integer=7)
    wts = logical_weights(k)
    xmask, zmask = pauli_masks(k)
    res = OrderedDict{String,Any}("code" => name, "wmax" => wmax, "pool" => pool, "seed" => seed,
                           "ensembles" => OrderedDict{String,Any}())
    all_targets = UInt32[i - 1 for i in eachindex(wts) if 1 <= wts[i] <= wmax]
    res["local_pauli_greedy_cover_all_weight_le_wmax"] = local_pauli_cover(all_targets, k, r)
    for ename in ENSEMBLES
        gens = generators[ename]
        lam = lambda_mix[ename]
        targets = [t for t in all_targets if lam[t+1] > 0]
        els = GF2Mat[GF2Mat(I(2k))]
        isempty(gens) || append!(els, random_elements(gens, pool, r))
        sets = Vector{Int}[]
        for m in els
            img = _image_of(targets, m)
            push!(sets, findall(==(0), img .& xmask))   # Z readout
            push!(sets, findall(==(0), img .& zmask))   # X readout
        end
        chosen, ncov = greedy_cover(sets, length(targets))
        res["ensembles"][ename] = OrderedDict{String,Any}(
            "n_learnable_targets" => length(targets), "n_targets_total" => length(all_targets),
            "greedy_cover_size" => length(chosen), "n_covered" => ncov,
            "local_pauli_cover_same_targets" => local_pauli_cover(targets, k, r))
    end
    return res
end
