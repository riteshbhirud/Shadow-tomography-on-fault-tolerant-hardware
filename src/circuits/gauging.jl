# Gauging logical measurement (Williamson and Yoder, arXiv:2410.02213) for any stabilizer code
# Measuring L on support V adds one |0> edge qubit per edge of a graph
# on V, Gauss checks A_v = L_v prod_{e at v} X_e, flux checks on a cycle basis, and deforms every
# original check by Z on edge paths pairing its L-anticommuting vertices.

anti(a::Char, b::Char) = a != b  # both are X, Y or Z

"""
    cheeger(G)

Cheeger constant min over |S| <= |V|/2 of |∂S| / |S|, by enumeration.
"""
function cheeger(G::PyGraph)
    nodes = G.nodes
    n = length(nodes)
    n <= 62 || throw(ArgumentError("graph too large for enumeration"))
    pos = Dict(v => i - 1 for (i, v) in enumerate(nodes))
    es = [(pos[u], pos[v]) for (u, v) in pyedges(G)]
    best = Inf
    for mask in 1:(UInt64(1) << n)-1
        size = count_ones(mask)
        size <= n ÷ 2 || continue
        cut = count(((a, b),) -> isodd((mask >> a) ⊻ (mask >> b)), es)
        best = min(best, cut / size)
    end
    return best
end

_path_edges(path) = [minmax(path[i], path[i+1]) for i in 1:length(path)-1]

function _perfect_matchings(S::Vector{Int})
    isempty(S) && return [Tuple{Int,Int}[]]
    a = S[1]
    out = Vector{Tuple{Int,Int}}[]
    for i in 2:length(S)
        rest = vcat(S[2:i-1], S[i+1:end])
        for m in _perfect_matchings(rest)
            push!(out, vcat([(a, S[i])], m))
        end
    end
    return out
end

function _best_matching(G::PyGraph, S::Vector{Int}, degree_aware::Bool)
    best = nothing
    for m in _perfect_matchings(S)
        new = [e for e in m if !has_edge(G, e...)]
        deg = Dict(v => degree(G, v) for v in G.nodes)
        for (a, b) in new
            deg[a] = get(deg, a, 0) + 1
            deg[b] = get(deg, b, 0) + 1
        end
        key = degree_aware ? (length(new), maximum(get(deg, v, 0) for v in S), sum(get(deg, v, 0) for v in S)) :
              (length(new), 0, 0)
        (best === nothing || key < best[1]) && (best = (key, m))
    end
    return best[2]
end

# Cap vertex degrees by rerouting an excess edge (u, x) as the path u-w-x. The candidate
# w is taken from a Python set, so its iteration order is CPython's.
function _reroute_degree!(G::PyGraph, pairings, max_degree::Int)
    pairings = [[collect(Int, p) for p in pr] for pr in pairings]
    changed = true
    while changed
        changed = false
        for u in sort(G.nodes; by=v -> -degree(G, v), alg=MergeSort)
            degree(G, u) <= max_degree && continue
            nbrs = sort(neighbors(G, u); by=v -> degree(G, v), alg=MergeSort)
            keep = pyset_order(nbrs[1:max_degree])
            for x in nbrs[max_degree+1:end]
                ws = [w for w in keep if w != x &&
                      (has_edge(G, w, x) || (degree(G, w) < max_degree && degree(G, x) <= max_degree))]
                isempty(ws) && continue
                w = ws[argmin([(!has_edge(G, v, x), degree(G, v)) for v in ws])]
                remove_edge!(G, u, x)
                add_edge!(G, w, x)
                for pr in pairings, i in eachindex(pr)
                    path = pr[i]
                    new = Int[]
                    for k in 1:length(path)-1
                        a, b = path[k], path[k+1]
                        push!(new, a)
                        minmax(a, b) == minmax(u, x) && push!(new, w)
                    end
                    push!(new, path[end])
                    pr[i] = new
                end
                changed = true
            end
            changed && break
        end
    end
    return G, pairings
end

"""
    routed_pairings(G, checks, L)

Per check, the matching of its L-anticommuting vertices with minimal total path length in `G`,
as vertex paths.
"""
function routed_pairings(G::PyGraph, checks, L::PauliDict)
    pairings = Vector{Vector{Int}}[]
    for s in checks
        S = sort([v for v in keys(L) if haskey(s, v) && anti(s[v], L[v])])
        if isempty(S)
            push!(pairings, Vector{Int}[])
            continue
        end
        ms = _perfect_matchings(S)
        best = ms[argmin([sum(shortest_path_length(G, a, b) for (a, b) in M) for M in ms])]
        push!(pairings, [shortest_path(G, a, b) for (a, b) in best])
    end
    return pairings
end

# The Cheeger padding breaks ties with numpy's default_rng(seed).random(), which Julia cannot
# replay. The edges used by the default graphs are stored (src/circuits/data/python_rng.json) and
# every imported edge is checked to maximise (Cheeger constant, -degree sum), the deterministic part
# of the choice. Other targets and seeds draw from numpy through the bridge.
const PYTHON_RNG = Ref{Any}(nothing)
function python_rng()
    PYTHON_RNG[] === nothing && (PYTHON_RNG[] = JSON.parsefile(joinpath(@__DIR__, "data", "python_rng.json")))
    return PYTHON_RNG[]
end

function _cheeger_padding!(G::PyGraph, L::PauliDict, target, max_extra, seed, max_degree)
    support = [[q, string(s)] for (q, s) in sort!(collect(L); by=first)]
    entry = nothing
    for e in python_rng()["cheeger_padding"]
        if e["support"] == support && e["target_cheeger"] == target && e["seed"] == seed && e["max_degree"] == max_degree
            entry = e
        end
    end
    V = sort(collect(keys(L)))
    entry === nothing && return _cheeger_padding_numpy!(G, V, target, max_extra, seed)
    h = cheeger(G)
    added = 0
    for (u, v) in entry["edges"]
        (h < target && added < max_extra) || error("exported Cheeger padding is longer than needed")
        keyof(a, b) = (add_edge!(G, a, b); hh = cheeger(G); remove_edge!(G, a, b); (hh, -(degree(G, a) + degree(G, b))))
        best = maximum(keyof(a, b) for a in V for b in V if a < b && !has_edge(G, a, b))
        (!has_edge(G, u, v) && keyof(u, v) == best) || error("exported Cheeger padding edge ($u, $v) is not a best choice")
        add_edge!(G, u, v)
        added += 1
        h = cheeger(G)
    end
    (h >= target || added == max_extra) || error("exported Cheeger padding does not reach the target")
    return G
end

"""numpy's `default_rng(seed)`, through the bridge."""
np_rng(seed::Integer) = pymodule(:numpy).random.default_rng(seed)

# The padding loop of gauging.build_graph, drawing its tie-breaks from numpy: one draw per candidate
# edge, in itertools.combinations order over the sorted support.
function _cheeger_padding_numpy!(G::PyGraph, V, target, max_extra, seed)
    rng = np_rng(seed)
    added = 0
    h = cheeger(G)
    while h < target && added < max_extra
        best = nothing
        for i in eachindex(V), j in i+1:length(V)
            u, v = V[i], V[j]
            has_edge(G, u, v) && continue
            add_edge!(G, u, v)
            hh = cheeger(G)
            remove_edge!(G, u, v)
            key = (hh, -(degree(G, u) + degree(G, v)), pyconvert(Float64, rng.random()))
            (best === nothing || key > best[1]) && (best = (key, (u, v)))
        end
        add_edge!(G, best[2]...)
        added += 1
        h = cheeger(G)
    end
    return G
end

"""
    build_graph(checks, L; extra_edges=(), target_cheeger=nothing, max_extra=50, seed=0,
                degree_aware=true, max_degree=nothing, edges=nothing)

Gauging graph on supp(L) and per check the vertex pairs or paths deforming it. With `edges`, the
graph is exactly that connected graph and checks are routed along shortest paths. Otherwise it
holds the matching edges, made connected, then grown until its Cheeger constant reaches
`target_cheeger`.
"""
function build_graph(checks, L::PauliDict; extra_edges=(), target_cheeger=nothing, max_extra=50, seed=0,
                     degree_aware=true, max_degree=nothing, edges=nothing)
    V = sort(collect(keys(L)))
    if edges !== nothing
        G = PyGraph()
        foreach(v -> add_node!(G, v), V)
        for (a, b) in edges
            add_edge!(G, a, b)
        end
        (Set(G.nodes) == Set(V) && is_connected(G)) ||
            throw(ArgumentError("edges must form a connected graph on the support of L"))
        return G, routed_pairings(G, checks, L)
    end
    G = PyGraph()
    foreach(v -> add_node!(G, v), V)
    pairings = Vector{Vector{Int}}[]
    for s in checks
        S = [v for v in V if haskey(s, v) && anti(s[v], L[v])]
        iseven(length(S)) || throw(ArgumentError("check anticommutes with L"))
        pairs = isempty(S) ? Tuple{Int,Int}[] : _best_matching(G, S, degree_aware)
        for (a, b) in pairs
            add_edge!(G, a, b)
        end
        push!(pairings, [[a, b] for (a, b) in pairs])
    end
    max_degree === nothing || ((G, pairings) = _reroute_degree!(G, pairings, max_degree))
    for (a, b) in extra_edges
        add_edge!(G, a, b)
    end
    comps = connected_components(G)
    while length(comps) > 1
        add_edge!(G, minimum(comps[1]), minimum(comps[2]))
        comps = connected_components(G)
    end
    target_cheeger === nothing || _cheeger_padding!(G, L, target_cheeger, max_extra, seed, max_degree)
    return G, pairings
end

"""
    cycle_basis_edges(G)

Short GF(2) cycle basis as sorted edge lists, shortest first, from shortest-path and BFS-tree
cycles (gauging.cycle_basis_edges).
"""
function cycle_basis_edges(G::PyGraph)
    edges = sorted_edges(G)
    ei = Dict(e => i - 1 for (i, e) in enumerate(edges))
    target = length(edges) - number_of_nodes(G) + length(connected_components(G))
    cands = Vector{Tuple{Int,Int}}[]
    for (u, v) in edges
        H = pycopy(G)
        remove_edge!(H, u, v)
        if has_path(H, u, v)
            push!(cands, vcat([(u, v)], _path_edges(shortest_path(H, u, v))))
        end
    end
    T = bfs_tree(G, minimum(G.nodes))
    for (u, v) in edges
        has_edge(T, u, v) || push!(cands, vcat([(u, v)], _path_edges(shortest_path(T, u, v))))
    end
    sort!(cands; by=length, alg=MergeSort)
    basis = Vector{Tuple{Int,Int}}[]
    rows = Dict{Int,Bits}()
    for c in cands
        w = Bits()
        for e in c
            w = xor(w, bit(ei[e]))
        end
        while !iszero(w)
            p = topbit(w)
            if haskey(rows, p)
                w = xor(w, rows[p])
            else
                rows[p] = w
                push!(basis, sort(unique(c)))
                break
            end
        end
        length(basis) == target && break
    end
    @assert length(basis) == target
    return basis
end

"""
    deform_operator(op, L, G, n_code)

`op` (commuting with L) times Z on shortest edge paths pairing its L-anticommuting vertices, so
that it commutes with every Gauss check and survives the gauged epoch.
"""
function deform_operator(op::PauliDict, L::PauliDict, G::PyGraph, n_code::Integer)
    eidx = Dict(e => n_code + i - 1 for (i, e) in enumerate(sorted_edges(G)))
    return _deform(op, [(L, G, eidx)])
end

function _deform(op::PauliDict, graphs)
    out = copy(op)
    for (L, G, eidx) in graphs
        rem = sort([v for v in keys(L) if haskey(op, v) && anti(op[v], L[v])])
        iseven(length(rem)) || error("operator anticommutes with L")
        while !isempty(rem)
            a = popfirst!(rem)
            b = rem[argmin([shortest_path_length(G, a, w) for w in rem])]
            deleteat!(rem, findfirst(==(b), rem))
            for e in _path_edges(shortest_path(G, a, b))
                _toggle_z!(out, eidx[e])
            end
        end
    end
    return out
end

_toggle_z!(p::PauliDict, q) = haskey(p, q) ? delete!(p, q) : (p[q] = 'Z')

function _gauss_order(order_seed, v, p::PauliDict)
    sp = sort(collect(keys(p)))
    key = "$(order_seed * 1000 + v):$(length(sp))"
    perm = get(python_rng()["gauss_perms"], key, nothing)
    perm === nothing || return sp[perm .+ 1]
    return gauss_order_numpy(order_seed * 1000 + v, sp)
end

"""
    gauss_order_numpy(seed, qubits)

`default_rng(seed).permutation(qubits)`: the Gauss-check gate order of gauging_experiment for seeds
that are not in the exported table.
"""
gauss_order_numpy(seed::Integer, qubits) = pyconvert(Vector{Int}, np_rng(seed).permutation(pylist(collect(Int, qubits))))

"""
    gauging_experiment(n_code, code_checks, L, rounds, readout_basis; state=(), logical_basis=nothing,
                       targets_extra=(), G=nothing, pairings=nothing, name="gauging", order_seed=nothing,
                       separate_gauss=false, flux_keep=nothing, graph_kw...)

Experiment measuring the logical `L` by gauging for `rounds` rounds, then reading the code
qubits out in `readout_basis` and the edge qubits in Z.
"""
function gauging_experiment(n_code::Integer, code_checks, L::PauliDict, rounds::Integer, readout_basis::Char;
                            state=(), logical_basis=nothing, targets_extra=(), G=nothing, pairings=nothing,
                            name="gauging", order_seed=nothing, separate_gauss=false, flux_keep=nothing, graph_kw...)
    code_orders = [c isa Check ? c.order : nothing for c in code_checks]
    plain = [c isa Check ? c.pauli : PauliDict(c) for c in code_checks]
    G === nothing && ((G, pairings) = build_graph(plain, L; graph_kw...))
    edges = sorted_edges(G)
    eidx = Dict(e => n_code + i - 1 for (i, e) in enumerate(edges))
    n = n_code + length(edges)
    V = sort(collect(keys(L)))
    checks = Check[]
    for v in V
        p = PauliDict(v => L[v])
        for e in edges
            v in e && (p[eidx[e]] = 'X')
        end
        od = order_seed === nothing ? nothing : _gauss_order(order_seed, v, p)
        push!(checks, Check(p; order=od, name="A_$v"))
    end
    cyc = cycle_basis_edges(G)
    flux_keep === nothing || (cyc = [cyc[i+1] for i in flux_keep])
    for c in cyc
        push!(checks, Check(PauliDict(eidx[e] => 'Z' for e in c); name="B"))
    end
    for (s, pairs, order) in zip(plain, pairings, code_orders)
        p = copy(s)
        for path in pairs, e in _path_edges(path)
            _toggle_z!(p, eidx[e])
        end
        od = order === nothing ? nothing : vcat(order, sort([q for q in keys(p) if q >= n_code]))
        push!(checks, Check(p; order=od, name=isempty(pairs) ? "orig" : "deformed"))
    end
    if separate_gauss
        gi = [i for (i, c) in enumerate(checks) if startswith(c.name, "A_")]
        rest = [i for (i, c) in enumerate(checks) if !startswith(c.name, "A_")]
        foreach(i -> checks[i].group = 0, gi)
        for (g_no, g) in enumerate(_greedy_groups(checks[rest], collect(0:length(rest)-1))), j in g
            checks[rest[j+1]].group = g_no
        end
    end
    readout = OrderedDict{Int,Char}(q => readout_basis for q in 0:n_code-1)
    for e in edges
        readout[eidx[e]] = 'Z'
    end
    exp = Experiment(; n_data=n, active0=0:n_code-1, prep_checks=plain, logical_state=collect(state),
                     epochs=[Epoch(checks, rounds; new_qubits=OrderedDict(eidx[e] => 'Z' for e in edges), name="gauged")],
                     readout, targets=vcat([copy(L)], collect(targets_extra)), name, logical_basis)
    exp.graph = G
    exp.graph_info = OrderedDict{String,Any}(
        "V" => length(V), "E" => length(edges), "flux" => length(cyc),
        "cheeger" => length(V) <= 22 ? cheeger(G) : nothing,
        "max_check_weight" => maximum(length(c.pauli) for c in checks),
        "max_degree" => maximum(degree(G, v) for v in G.nodes))
    return exp
end

"""
    multi_gauging_experiment(n_code, code_checks, Ls, rounds, readout_basis; state=(), logical_basis=nothing,
                             name="multi_gauging", extra_ops=())

Simultaneous gauging of logicals with pairwise-disjoint supports, one default graph each
(gauging.multi_gauging_experiment). Targets are the logicals, then `extra_ops` deformed through
every graph.
"""
function multi_gauging_experiment(n_code::Integer, code_checks, Ls, rounds::Integer, readout_basis::Char;
                                  state=(), logical_basis=nothing, name="multi_gauging", extra_ops=())
    supp = [Set(keys(L)) for L in Ls]
    for i in eachindex(Ls), j in i+1:length(Ls)
        isempty(supp[i] ∩ supp[j]) || throw(ArgumentError("multi_gauging requires disjoint supports"))
    end
    code_orders = [c isa Check ? c.order : nothing for c in code_checks]
    plain = [c isa Check ? c.pauli : PauliDict(c) for c in code_checks]
    graphs = [build_graph(plain, L) for L in Ls]
    all_edges = Tuple{Int,Tuple{Int,Int}}[]
    for (gi, (G, _)) in enumerate(graphs)
        append!(all_edges, [(gi, minmax(e...)) for e in sort(pyedges(G))])
    end
    eidx = Dict(ge => n_code + i - 1 for (i, ge) in enumerate(all_edges))
    n = n_code + length(all_edges)
    checks = Check[]
    for (gi, (L, (G, _))) in enumerate(zip(Ls, graphs))
        for v in sort(collect(keys(L)))
            p = PauliDict(v => L[v])
            for e in pyedges(G)
                e = minmax(e...)
                v in e && (p[eidx[(gi, e)]] = 'X')
            end
            push!(checks, Check(p; name="A$(gi-1)_$v"))
        end
        for c in cycle_basis_edges(G)
            push!(checks, Check(PauliDict(eidx[(gi, e)] => 'Z' for e in c); name="B$(gi-1)"))
        end
    end
    for (ci, (s, order)) in enumerate(zip(plain, code_orders))
        p = copy(s)
        for (gi, (_, pairings)) in enumerate(graphs), path in pairings[ci], e in _path_edges(path)
            _toggle_z!(p, eidx[(gi, e)])
        end
        od = order === nothing ? nothing : vcat(order, sort([q for q in keys(p) if q >= n_code]))
        push!(checks, Check(p; order=od, name="deformed"))
    end
    readout = OrderedDict{Int,Char}(q => readout_basis for q in 0:n_code-1)
    for (_, q) in sort!(collect(eidx); by=last)
        readout[q] = 'Z'
    end
    geidx = [(L, G, Dict(e => eidx[(gi, e)] for e in sorted_edges(G))) for (gi, (L, (G, _))) in enumerate(zip(Ls, graphs))]
    extra = [_deform(PauliDict(op), geidx) for op in extra_ops]
    exp = Experiment(; n_data=n, active0=0:n_code-1, prep_checks=plain, logical_state=collect(state),
                     epochs=[Epoch(checks, rounds; new_qubits=OrderedDict(q => 'Z' for (_, q) in sort!(collect(eidx); by=last)),
                                   name="gauged")],
                     readout, targets=vcat([copy(L) for L in Ls], extra), name, logical_basis)
    exp.graph_info = OrderedDict{String,Any}("E" => length(all_edges), "V" => sum(length, Ls),
                                             "flux" => count(c -> startswith(c.name, "B"), checks))
    return exp
end
