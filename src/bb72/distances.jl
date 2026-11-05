# Deformed-code and circuit distances of BB72 gauging circuits (used by the scripts/taskA_bb*_distance,
# schedule, flux, multi and small_graph scripts). Circuit distances come from stim's search with at most
# 6 detection events and edge degree 8.

using .Circuits: PauliDict, bb_setup, bb_gauging, build, stim_search, ordered_checks, cheeger,
                 cycle_basis_edges, code_distance_isd, code_distance_stim, PyGraph, add_node!, add_edge!,
                 is_connected, number_of_nodes, degree, routed_pairings, gauging_experiment, deform_operator,
                 multi_gauging_experiment, local_experiment, compile_setting, memory_experiment, css_basis,
                 dict_to_vec, vec_to_dict, RowSpace, addrow!, nullspace_gf2, logical_basis_with, orbit,
                 random_stabilizer_state, exact_expectation, symp_inner, READOUT

"""Circuit distance of `exp` at r rounds: stim search on the eigen-target build at p = 1e-3."""
circuit_distance(exp) = stim_search(build(exp, 1e-3; eigen_targets = true).circuit; max_dets = 6, max_degree = 8)

"""Plain checks of the (only) epoch of a gauging experiment."""
epoch_checks(exp) = [c.pauli for c in exp.epochs[1].checks]

"""(ISD upper bound, stim distance, k) of the deformed code of a gauging experiment."""
function deformed_distances(exp)
    ch = epoch_checks(exp)
    d_isd, _ = code_distance_isd(exp.n_data, ch; iters = 300)
    d_stim, k = code_distance_stim(exp.n_data, ch)
    return d_isd, d_stim, k
end

"""
    bb_distance_entry(name; rounds=6) -> (name, row)

One entry of bb72_distances.json: "X", "Z", "Y" (default graphs) or "Y_minimal" (the h = 0.6 graph
without the Cheeger padding), with all spectator logicals and with the measured logical only.
"""
function bb_distance_entry(name::AbstractString; rounds::Integer = 6)
    t = time()
    P = name[1]
    kw = name == "Y_minimal" ? (use_defaults = false,) : (;)
    exp = bb_gauging(P, rounds; extra_targets = true, kw...)
    d_isd, d_stim, k = deformed_distances(exp)
    d_circ = circuit_distance(exp)
    d_L = circuit_distance(bb_gauging(P, rounds; extra_targets = false, kw...))
    return name, OrderedDict{String,Any}("graph" => exp.graph_info, "n_targets" => length(exp.targets),
        "deformed_k" => k, "deformed_distance_isd_upper" => d_isd, "deformed_distance_stim" => d_stim,
        "circuit_distance_all_targets" => d_circ, "circuit_distance_L_only" => d_L, "rounds" => rounds,
        "seconds" => time() - t)
end

"""Circuit distance of Pbar_0 gauging at r = 6 with Gauss gate orders from `seed` (< 0: the default)."""
function gauss_order_distance(P::Char, seed::Integer, graph_kw)
    kw = (Symbol(k) => v for (k, v) in graph_kw)
    return circuit_distance(bb_gauging(P, 6; order_seed = seed < 0 ? nothing : seed, kw...))
end

#### syndrome schedules and memory

const _BB72_CSS = Ref{Any}(nothing)
"""BB72 and its CSS logical basis (cached)."""
function bb72_css()
    if _BB72_CSS[] === nothing
        code = Circuits.BBCode("bb72")
        _BB72_CSS[] = (code, css_basis(code))
    end
    return _BB72_CSS[]
end

"""
    schedule_distance(kind, perm) -> (kind, perm, distance)

Circuit distance of a BB72 memory at r = 3 when the X (`kind = 'X'`) or Z check gate slots are
permuted by `perm`; X-check hooks are X errors, so they are tested in Z memory and vice versa.
"""
function schedule_distance(kind::Char, perm)
    code, basis = bb72_css()
    ox = kind == 'X' ? perm : 0:5
    oz = kind == 'Z' ? perm : 0:5
    exp = memory_experiment(code, ordered_checks(code; order_x = ox, order_z = oz), kind == 'X' ? 'Z' : 'X', 3, basis)
    return kind, perm, circuit_distance(exp)
end

"""Summary of a schedule search, per kind: count, best distance, how many reach it, first best, histogram."""
function schedule_summary(res)
    out = OrderedDict{String,Any}()
    for (k, v) in res
        best = maximum(values(v))
        out[k] = OrderedDict{String,Any}("n" => length(v), "best" => best, "n_best" => count(==(best), values(v)),
            "first_best" => minimum(p for (p, x) in v if x == best),
            "histogram" => OrderedDict(string(d) => count(==(d), values(v)) for d in sort(unique(values(v)))))
    end
    return out
end

#### local settings

"""
    local_distance_combos(ms) -> OrderedDict(key => example setting)

Every distinct gauging combination of the local settings on m ∈ ms qubits, keyed
"R:P0,Q1" (readout, then gauged logicals), with the first setting that uses it.
"""
function local_distance_combos(ms)
    out = OrderedDict{String,String}()
    for m in ms, bases in product_strings("XYZ", m)
        c = compile_setting(bases)
        isempty(c.gauged) && continue
        key = string(c.readout) * ":" * join(("$P$i" for (i, P) in c.gauged), ",")
        haskey(out, key) || (out[key] = bases)
    end
    return out
end

"""Fields of a local-distance row from the deformed code."""
function local_deformed_fields(bases, rounds)
    exp, _, comp = local_experiment(bases, [], rounds)
    d_isd, d_stim, k = deformed_distances(exp)
    return OrderedDict{String,Any}("example_setting" => bases, "readout" => string(comp.readout),
        "gauged" => jsonable(comp.gauged), "graph" => exp.graph_info, "deformed_k" => k,
        "deformed_distance_stim" => d_stim, "deformed_distance_isd_upper" => d_isd, "rounds" => rounds)
end

local_distance_L(bases, rounds) = circuit_distance(local_experiment(bases, [], rounds)[1])

function local_distance_all(bases, rounds)
    exp = local_experiment(bases, [], rounds; spectators = true)[1]
    return length(exp.targets), circuit_distance(exp)
end

#### incomplete flux (Ybar_0)

"""Indices (0-based, sorted) of the two shortest cycles, ties to the lower index."""
two_shortest(lens) = sort(argsort_stable(lens)[1:2] .- 1)

"""
    incomplete_flux_rows(rounds) -> Vector

For the adopted E20 and the minimal E18 Ybar graphs, keep the two shortest flux checks, three random pairs
(numpy `default_rng(20261030)`, so the pairs match the stored rows) or all of them, and search the circuit
distance with the measured logical only.
"""
function incomplete_flux_rows(rounds::Integer)
    out = OrderedDict{String,Any}[]
    for (label, kw) in (("adopted_E20_h1", (;)), ("minimal_E18_h0.6", (use_defaults = false,)))
        full = bb_gauging('Y', rounds; kw...)
        cyc = cycle_basis_edges(full.graph)
        lens = length.(cyc)
        rng = Circuits.np_rng(20261030)
        choices = Any[("shortest2", two_shortest(lens))]
        for t in 0:2
            push!(choices, ("random2_$t", sort(pyconvert(Vector{Int}, rng.choice(length(cyc), 2; replace = false)))))
        end
        push!(choices, ("all", nothing))
        for (cname, keep) in choices
            exp = bb_gauging('Y', rounds; flux_keep = keep, kw...)
            exp.tracker_readout_fallback = true  # unmeasured cycles still give readout detectors
            d = circuit_distance(exp)
            gi = exp.graph_info
            row = OrderedDict{String,Any}("graph" => label, "flux_choice" => cname, "flux_kept" => keep,
                "cycle_lengths" => lens, "V" => gi["V"], "E" => gi["E"], "flux" => gi["flux"],
                "extra_qubits" => gi["V"] + gi["E"] + gi["flux"], "cheeger" => Float64(cheeger(full.graph)),
                "circuit_distance_L_only" => d)
            println(pyjson(row))
            flush(stdout)
            push!(out, row)
        end
    end
    return out
end

"""
    incomplete_flux_spectators(rounds; control=false) -> row

The 33-qubit reconstruction (adopted graph, two shortest flux cycles; all cycles with `control`)
with every readable spectator logical observed.
"""
function incomplete_flux_spectators(rounds::Integer; control::Bool = false)
    full = bb_gauging('Y', rounds)
    keep = control ? nothing : two_shortest(length.(cycle_basis_edges(full.graph)))
    exp = bb_gauging('Y', rounds; flux_keep = keep, extra_targets = true)
    exp.tracker_readout_fallback = true
    d = circuit_distance(exp)
    return OrderedDict{String,Any}("graph" => "adopted_E20_h1", "flux_kept" => keep, "n_targets" => length(exp.targets),
        "rounds" => rounds, "circuit_distance_all_targets" => d)
end

#### two Xbar at once

"""
    multi_pair() -> (code, [L1, L2])

Xbar_0 and the first weight-6 Xbar in its shift orbit that is disjoint from it and logically
independent (v1, v2 and v1 + v2 outside the stabilizer span).
"""
function multi_pair()
    s = bb_setup()
    code = s.code
    X = sort(collect(keys(s.L['X'])))
    for b in orbit(code, X)
        isempty(b ∩ X) || continue
        v1 = dict_to_vec(PauliDict(q => 'X' for q in X), code.n)
        v2 = dict_to_vec(PauliDict(q => 'X' for q in b), code.n)
        if all(v -> addrow!(RowSpace(code.S), v), (v1, v2, v1 .⊻ v2))
            return code, [PauliDict(q => 'X' for q in X), PauliDict(q => 'X' for q in b)]
        end
    end
    error("no disjoint independent pair")
end

"""Symplectic logical basis (Pauli dicts) starting with the two measured logicals."""
multi_basis(code, Ls) = [(vec_to_dict(a), vec_to_dict(b)) for (a, b) in logical_basis_with(code, [dict_to_vec(L, code.n) for L in Ls])]

"""Independent Z-type logicals commuting with both measured logicals."""
function multi_spectators(code, Ls)
    n = code.n
    lv = [dict_to_vec(L, n) for L in Ls]
    out = PauliDict[]
    cur = RowSpace(code.S)
    NX = nullspace_gf2(code.HX)
    for i in axes(NX, 1)
        z = vcat(falses(n), NX[i, :])
        any(L -> symp_inner(z, L), lv) && continue
        addrow!(cur, z) && push!(out, vec_to_dict(z))
    end
    return out
end

"""Simultaneous gauging of the two Xbar for `rounds` rounds, Z readout."""
multi_experiment(code, Ls, rounds; kw...) =
    multi_gauging_experiment(code.n, ordered_checks(code), Ls, rounds, 'Z'; logical_basis = multi_basis(code, Ls), kw...)

#### 6-edge graphs for Xbar_0 / Zbar_0

"""Gauging graph on the sorted support of L with exactly these edges (nodes added first)."""
function edge_graph(L::PauliDict, edges)
    G = PyGraph()
    foreach(v -> add_node!(G, v), sort(collect(keys(L))))
    foreach(e -> add_edge!(G, e[1], e[2]), edges)
    return G
end

"""Graph of an edge list alone."""
function graph_of_edges(edges)
    G = PyGraph()
    foreach(e -> add_edge!(G, e[1], e[2]), edges)
    return G
end

"""
    small_graph_experiment(P, edges; rounds=6, spectators=false)

Gauging of Pbar_0 on the given graph, checks routed along shortest paths.
With `spectators`, every logical the readout reads is observed too, deformed through the graph.
"""
function small_graph_experiment(P::Char, edges; rounds::Integer = 6, spectators::Bool = false)
    s = bb_setup()
    code, basis, L = s.code, s.basis, s.L[P]
    G = edge_graph(L, edges)
    chk = ordered_checks(code)
    pair = routed_pairings(G, [c.pauli for c in chk], L)
    rb = READOUT[P]
    extra = PauliDict[]
    if spectators
        for (Xi, Zi) in basis[2:end], op in (Xi, Zi)
            all(==(rb), values(op)) && push!(extra, deform_operator(op, L, G, code.n))
        end
    end
    return gauging_experiment(code.n, chk, L, rounds, rb; logical_basis = basis, targets_extra = extra, G, pairings = pair,
                              name = "small_$P")
end

"""One screened graph: its circuit distance and graph data."""
function small_graph_job(P::Char, edges, spectators::Bool; rounds::Integer = 6)
    t = time()
    exp = small_graph_experiment(P, edges; rounds, spectators)
    d = circuit_distance(exp)
    return OrderedDict{String,Any}("P" => string(P), "edges" => [collect(e) for e in edges], "E" => length(edges),
        "spectators" => spectators, "distance" => d, "cheeger" => Float64(cheeger(exp.graph)),
        "max_check_weight" => exp.graph_info["max_check_weight"], "max_degree" => exp.graph_info["max_degree"],
        "extra_qubits" => 2 * length(edges) + 1, "seconds" => time() - t)
end

"""
    small_graph_candidates(P) -> [(E, edges, cheeger)]

Every connected graph with 5..8 edges on the 6 support qubits of Pbar_0, in lexicographic order.
"""
function small_graph_candidates(P::Char)
    V = sort(collect(keys(bb_setup().L[P])))
    allpairs = [Tuple(c) for c in py_combinations(V, 2)]
    out = Tuple{Int,Vector{Tuple{Int,Int}},Float64}[]
    for E in 5:8, es in py_combinations(allpairs, E)
        G = graph_of_edges(es)
        number_of_nodes(G) == length(V) && is_connected(G) && push!(out, (E, es, cheeger(G)))
    end
    return out
end

"""Largest vertex degree of the graph with these edges."""
max_degree_of(es) = (G = graph_of_edges(es); maximum(degree(G, v) for v in G.nodes))

"""
    noiseless_gauging_check(P, edges) -> (states, mismatches)

Noiseless Pbar_0 gauging at r = 1 on this graph against the exact expectation, on 2 eigenstates
and 4 random 12-qubit states.
"""
function noiseless_gauging_check(P::Char, edges)
    bad, n = 0, 0
    lab = string(P) * "I"^11
    states = Vector{Tuple{String,Int}}[]
    for sign in (0, 1)
        rest = random_stabilizer_state(11, 50 + sign)
        push!(states, vcat([(lab, sign)], [("I" * l, s) for (l, s) in rest]))
    end
    append!(states, [random_stabilizer_state(12, 60 + t) for t in 0:3])
    for (t, st) in enumerate(states)
        ex = exact_expectation(st, lab)
        res = run_oracle(bb_gauging(P, 1; state_labels = st, edges), 0.0, 200; decoder = "relay", seed = t - 1)
        x = 1 .- 2 .* Int.(res.outcomes[:, 1])
        n += 1
        bad += ex in (1, -1) ? Int(any(x .!= ex)) : Int(length(unique(x)) != 2)
    end
    return OrderedDict{String,Any}("states" => n, "mismatches" => bad)
end
