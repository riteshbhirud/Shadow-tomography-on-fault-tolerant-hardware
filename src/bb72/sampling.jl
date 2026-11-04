# Sampled BB72 studies: q(r) calibration, local-setting flip rates, gap postselection, the 6-edge
# graphs' q, memory sanity, the gap-decoder comparison and noiseless checks (scripts taskA_bb72_*).
# Every chunk uses the reference's seed_for key, so stim draws the reference's shots.

const BB72_P = 1e-3

#### q(r) of Xbar_0, Zbar_0, Ybar_0 gauging (taskA_bb72_calibration)

const BB72_QR_CHUNK = 5000

"""
    bb72_qr_chunk(job) -> (err, gap, cost, graph_info)

One chunk `job = (P, r, state, seed, shots, small, chunk_dir)`, cached in `chunk_dir` as
`{P}_r{r}_{state}_{seed}.npz` (the reference's file), so a restarted run resumes seed-identically.
"""
function bb72_qr_chunk(job)
    Pl, r, state, seed, shots, small, cdir = job
    f = joinpath(cdir, "$(Pl)_r$(r)_$(state)_$(seed).npz")
    if isfile(f)
        z = npz_load(f)
        return z["err"], z["gap"], z["cost"], z["g"]
    end
    exp = bb_gauging(Pl[1], r; state_kind = state, (small ? pairs(Circuits.SMALL_GRAPHS[Pl[1]]) : ())...)
    res = run_oracle(exp, BB72_P, shots; decoder = "relay", seed, with_gap = false)
    err = UInt8.(res.outcomes[:, 1] .!= (state == "+1" ? 0 : 1))
    atomic_save(f, joinpath(cdir, "$(Pl)_r$(r)_$(state)_$(seed).tmp.npz")) do tmp
        npz_save_compressed(tmp; err, gap = res.confidence, cost = res.cost, g = exp.graph_info)
    end
    return err, res.confidence, res.cost, exp.graph_info
end

"""
    bb72_qr_adaptive(cfgs, target, max_shots; small=false, chunk_dir, pmap=map, chunk=BB72_QR_CHUNK)

Three chunks per active configuration `(P, r, state)` at a time until it has `target` errors or
`max_shots` shots (taskA_bb72_calibration.adaptive). Returns the rows of bb72_qr.csv.
"""
function bb72_qr_adaptive(cfgs, target, max_shots; small::Bool = false, chunk_dir, pmap = map, chunk = BB72_QR_CHUNK)
    st = [Dict{String,Any}("err" => Vector{UInt8}[], "seed" => 0, "cost" => nothing, "g" => nothing) for _ in cfgs]
    active = Set(eachindex(cfgs))
    while !isempty(active)
        jobs, own = Tuple[], Int[]
        for i in sort!(collect(active)), _ in 1:3
            push!(jobs, (cfgs[i]..., seed_for("bb72_q", cfgs[i]..., st[i]["seed"]), chunk, small, chunk_dir))
            st[i]["seed"] += 1
            push!(own, i)
        end
        for (i, (e, _, cost, gi)) in zip(own, pmap(bb72_qr_chunk, jobs))
            push!(st[i]["err"], e)
            st[i]["cost"], st[i]["g"] = cost, gi
        end
        errs(i) = sum(x -> count(!iszero, x), st[i]["err"])
        shots(i) = sum(length, st[i]["err"])
        filter!(i -> errs(i) < target && shots(i) < max_shots, active)
        println("  active $(length(active)): ", join(("$(cfgs[i]):$(errs(i))/$(shots(i))" for i in sort!(collect(active))), ", "))
        flush(stdout)
    end
    rows = OrderedDict{String,Any}[]
    for (i, (Pl, r, state)) in enumerate(cfgs)
        e = reduce(vcat, st[i]["err"])
        k, n = count(!iszero, e), length(e)
        lo, hi = wilson(k, n)
        row = OrderedDict{String,Any}("P" => Pl, "r" => r, "state" => state, "p" => BB72_P, "shots" => n, "errors" => k,
            "q" => k / n, "ci_lo" => lo, "ci_hi" => hi)
        foreach(((a, b),) -> row["cost_$a"] = b, st[i]["cost"])
        foreach(((a, b),) -> row["graph_$a"] = b, st[i]["g"])
        push!(rows, row)
    end
    return rows
end

"""
    bb72_qr_configs(ps; rounds=1:6, minus_rounds=(2, 3))

Configurations `(P, r, state)` of taskA_bb72_calibration: every r in `rounds` on +1 eigenstates, then
the r in `minus_rounds` on -1 eigenstates (state independence).
"""
bb72_qr_configs(ps; rounds = 1:6, minus_rounds = (2, 3)) =
    vcat([(string(P), r, "+1") for P in ps for r in rounds], [(string(P), r, "-1") for P in ps for r in minus_rounds])

#### flip rates of the compiled local settings (taskA_bb72_local_calibration)

const BB72_LOCAL_CHUNK = 1000

"""Chunk `c` of setting `bases` at r rounds on the joint +1 eigenstate, saved as `{bases}_r{r}_c{c}.npy`."""
function bb72_local_chunk(job)
    bases, r, c, cdir, chunk = job
    f = joinpath(cdir, "$(bases)_r$(r)_c$(c).npy")
    if !isfile(f)
        st = [(join(j == i ? b : 'I' for j in 1:length(bases)), 0) for (i, b) in enumerate(bases)]
        exp, labels, _ = local_experiment(bases, st, r)
        res = run_oracle(exp, BB72_P, chunk; decoder = "relay", seed = seed_for("bb72_local_q", bases, r, c), with_gap = false,
                         setting_labels = labels)
        atomic_save(tmp -> npy_save(tmp, res.outcomes), f, joinpath(cdir, "$(bases)_r$(r)_c$(c).tmp.npy"))
    end
    return bases, r
end

"""
    bb72_local_rows(bases, r, o) -> Vector

Rows of bb72_local_q.csv for one setting from its outcome bits `o` (shots × m): every string and
every subset product, with Wilson CI, the independent-flip prediction and the pair correlation.
"""
function bb72_local_rows(bases::AbstractString, r::Integer, o::AbstractMatrix{UInt8})
    comp = compile_setting(bases)
    n, m = size(o)
    rows = OrderedDict{String,Any}[]
    for k in 1:m, sub in py_combinations(collect(0:m-1), k)
        lab = join(q in sub ? bases[q+1] : 'I' for q in 0:m-1)
        f = reduce(.⊻, (o[:, j+1] for j in sub))
        nf = count(!iszero, f)
        lo, hi = wilson(nf, n)
        q1 = [bitmean(o[:, j+1]) for j in sub]
        a, b = k == 2 ? (o[:, sub[1]+1], o[:, sub[2]+1]) : (nothing, nothing)
        corr = k == 2 && any(!iszero, a) && any(!iszero, b) ? np_corrcoef(a, b) : NaN
        push!(rows, OrderedDict{String,Any}("code" => "bb72", "setting_id" => "local:$bases", "readout" => string(comp.readout),
            "gauged" => join(("$P_$i" for (i, P_) in comp.gauged), ";"), "r" => isempty(comp.gauged) ? 0 : r,
            "rounds_parameter" => r, "noise" => "bravyi_cal", "p" => BB72_P, "target" => lab,
            "kind" => k == 1 ? "string" : "product", "shots" => n, "flips" => nf, "q" => nf / n, "ci_lo" => lo, "ci_hi" => hi,
            "q_if_independent" => (1 - foldl(*, [1 - 2x for x in q1])) / 2, "flip_correlation" => corr))
    end
    return rows
end

#### gap postselection (taskA_bb72_postselection)

const BB72_POST_CHUNK = 500

"""Chunk of Pbar_0 gauging at r rounds with the plain-anchored Relay gap, saved as `{P}_r{r}_c{c}.npz`."""
function bb72_post_chunk(job)
    Pl, r, c, cdir, chunk = job
    f = joinpath(cdir, "$(Pl)_r$(r)_c$(c).npz")
    if !isfile(f)
        res = run_oracle(bb_gauging(Pl[1], r), BB72_P, chunk; decoder = "relay", seed = seed_for("bb72_post", Pl, r, c))
        atomic_save(f, joinpath(cdir, "$(Pl)_r$(r)_c$(c).tmp.npz")) do tmp
            npz_save_compressed(tmp; err = UInt8.(res.outcomes[:, 1]), gap = res.confidence)
        end
    end
    z = npz_load(f)
    return (Pl, r), z["err"], z["gap"]
end

"""Order of shots by decreasing gap (NaN as 0), stable: the postselection order."""
gap_order(g) = argsort_stable(.-ifelse.(isnan.(g), 0.0, g))

"""Row of bb72_postselection.csv: q overall and among the shots kept at several fractions."""
function postselection_row(Pl, r, e::AbstractVector, g::AbstractVector)
    o = gap_order(g)
    n = length(e)
    lo, hi = wilson(count(!iszero, e), n)
    row = OrderedDict{String,Any}("P" => Pl, "r" => r, "shots" => n, "errors" => count(!iszero, e), "q" => bitmean(e),
        "ci_lo" => lo, "ci_hi" => hi, "inf_gap_frac" => count(isinf, g) / n)
    for f in (0.99, 0.95, 0.9, 0.75, 0.5)
        kept = e[o][1:trunc(Int, f * n)]
        row["q_keep$(trunc(Int, f * 100))"] = isempty(kept) ? NaN : bitmean(kept)
        row["errors_keep$(trunc(Int, f * 100))"] = count(!iszero, kept)
    end
    return row
end

#### q of the 6-edge graphs (taskA_bb72_small_graph_q)

const SMALL_Q_CHUNK = 2500

"""Pbar_0 gauging on a given graph for r = 3 rounds on the +1 eigenstate (taskA_bb72_small_graph_q.experiment)."""
function small_graph_q_experiment(P::Char, edges; rounds::Integer = 3)
    s = bb_setup()
    code, L, basis = s.code, s.L[P], s.bases[P]
    gens = vcat([(L, 0)], [(Zi, 0) for (_, Zi) in basis[2:end]])
    G = edge_graph(L, edges)
    chk = ordered_checks(code)
    pair = routed_pairings(G, [c.pauli for c in chk], L)
    return gauging_experiment(code.n, chk, L, rounds, READOUT[P]; state = gens, logical_basis = basis, G, pairings = pair,
                              name = "small_$(P)_r$(rounds)")
end

"""Graphs with spectator circuit distance 6 in the spectator screen, as (P, edge tuple)."""
function small_graph_q_candidates(screen)
    out = Tuple{String,Tuple}[]
    for (P, blocks) in screen, blk in blocks, row in blk["rows"]
        row["distance"] == 6 && push!(out, (String(P), Tuple(Tuple(Int.(e)) for e in row["edges"])))
    end
    return out
end

"""Chunk `c` of candidate `i` (0-based): errors of Pbar_0, saved as `{i}_c{c}.npy`."""
function small_graph_q_chunk(job)
    i, P, edges, c, cdir, chunk = job
    f = joinpath(cdir, "$(i)_c$(c).npy")
    if !isfile(f)
        res = run_oracle(small_graph_q_experiment(P[1], edges), BB72_P, chunk; decoder = "relay",
                         seed = seed_for("bb72_small_q", P, edges, c), with_gap = false)
        atomic_save(tmp -> npy_save(tmp, UInt8.(res.outcomes[:, 1])), f, joinpath(cdir, "$(i)_c$(c).tmp.npy"))
    end
    return i, count(!iszero, pyconvert(Vector{UInt8}, pymodule(:numpy).load(f)))
end

#### memory sanity (taskA_bb72_memory_sanity)

"""Failures in `shots` shots of BB72 memory (r = 6, p = 1e-3, Relay-BP) under `noise` in basis `bk`."""
function bb72_memory_failures(job)
    seed, shots, noise, bk = job
    code, basis = bb72_css()
    exp = memory_experiment(code, ordered_checks(code), bk[1], 6, basis)
    res = run_oracle(exp, BB72_P, shots; decoder = "relay", seed, with_gap = false, noise)
    return count(any(!iszero, res.outcomes[s, :]) for s in axes(res.outcomes, 1)), shots
end

#### gap decoders (taskA_bb72_gap_decoder_comparison)

const GAP_P = 3e-3
const GAP_VARIANTS = OrderedDict{String,Tuple{String,Any}}(
    "anchored_relay100" => ("relay", (;)),
    "anchored_relay20" => ("relay", (comp_kw = (num_sets = 20,),)),
    "anchored_bposd" => ("relay", (comp = "bposd",)),
    "split_relay20" => ("relay_split", (num_sets = 20,)),
    "anchored_relay5" => ("relay", (comp_kw = (num_sets = 5,),)),
)

"""
    gap_variant(name; shots=1000) -> (name, row)

The same Xbar_0 gauging shots (r = 3, p = 3e-3, sd6 noise) decoded with one gap decoder: time per
shot, q, AUC of the gap against errors, q of the 90% and 75% largest-gap shots, fraction of infinite
gaps, and the outcomes (for the same-prediction check).
"""
function gap_variant(name::AbstractString; shots::Integer = 1000)
    dec, kw = GAP_VARIANTS[name]
    t = time()
    res = run_oracle(bb_gauging('X', 3), GAP_P, shots; decoder = dec, seed = seed_for("gapcmp", GAP_P), decoder_kwargs = kw,
                     noise = "sd6")
    dt = (time() - t) / shots
    err = res.outcomes[:, 1] .!= 0
    g = ifelse.(isnan.(res.confidence), 0.0, res.confidence)
    o = argsort_stable(.-g)
    keep(f) = (k = err[o][1:trunc(Int, f * shots)]; isempty(k) ? NaN : bitmean(k))
    return name, OrderedDict{String,Any}("s_per_shot" => dt, "q" => bitmean(err), "errors" => count(err),
        "auc" => auc(.-g, err), "q_keep90" => keep(0.9), "q_keep75" => keep(0.75), "inf_frac" => count(isinf, g) / shots,
        "outcomes" => Int.(res.outcomes[:, 1]))
end

#### noiseless exactness of the local settings (taskA_bb72_local_noiseless)

"""
    local_noiseless_check(bases, st, seed) -> (products checked, mismatches, compiled setting)

Noiseless r = 1 run of the setting on state `st`; every subset product of its outcomes must equal
the exact expectation (±1), or take both values when that is 0.
"""
function local_noiseless_check(bases::AbstractString, st, seed::Integer)
    exp, labels, comp = local_experiment(bases, st, 1)
    r = run_oracle(exp, 0.0, 200; decoder = "relay", seed, setting_labels = labels)
    x = 1 .- 2 .* Int.(r.outcomes)
    m, bad, n = length(bases), 0, 0
    for k in 1:length(labels), sub in py_combinations(collect(1:length(labels)), k)
        P = join(_first_non_identity([labels[j][q] for j in sub]) for q in 1:m)
        e = exact_expectation(st, P)
        v = [prod(x[s, sub]) for s in axes(x, 1)]
        n += 1
        bad += e in (1, -1) ? Int(any(v .!= e)) : Int(length(unique(v)) != 2)
    end
    return n, bad, comp
end

_first_non_identity(cs) = (i = findfirst(!=('I'), cs); i === nothing ? 'I' : cs[i])

"""Rows of bb72_local_noiseless.json for m = 2 and 3: 3 random states and the joint eigenstate per setting."""
function local_noiseless_row(bases::AbstractString)
    m = length(bases)
    states = [random_stabilizer_state(m, seed_for("bb72_local_state", m, t)) for t in 0:2]
    eig = [(join(j == i ? b : 'I' for j in 1:m), 0) for (i, b) in enumerate(bases)]
    tot = bad = 0
    comp = nothing
    for (t, st) in enumerate(vcat(states, [eig]))
        n, b, comp = local_noiseless_check(bases, st, t - 1)
        tot += n
        bad += b
    end
    return OrderedDict{String,Any}("setting" => bases, "readout" => string(comp.readout), "gauged" => jsonable(comp.gauged),
        "products_checked" => tot, "mismatches" => bad)
end
