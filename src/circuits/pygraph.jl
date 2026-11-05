# Gauging graphs. The deformed checks and flux cycles depend on adjacency order through shortest
# paths (bidirectional BFS), BFS trees and edge iteration. The stored circuits were built with
# networkx, so `PyGraph` keeps networkx 3.7's node and adjacency order for the operations used.
# Order-free questions go to Graphs.jl.

"""
    PyGraph()

Undirected graph with networkx's node and adjacency order.
"""
struct PyGraph
    nodes::Vector{Int}
    adj::Dict{Int,Vector{Int}}
end
PyGraph() = PyGraph(Int[], Dict{Int,Vector{Int}}())

function add_node!(G::PyGraph, v::Integer)
    haskey(G.adj, v) && return G
    push!(G.nodes, v)
    G.adj[v] = Int[]
    return G
end

function add_edge!(G::PyGraph, u::Integer, v::Integer)
    add_node!(G, u)
    add_node!(G, v)
    v in G.adj[u] || push!(G.adj[u], v)
    u in G.adj[v] || push!(G.adj[v], u)
    return G
end

function remove_edge!(G::PyGraph, u::Integer, v::Integer)
    has_edge(G, u, v) || throw(ArgumentError("the edge $u-$v is not in the graph"))
    deleteat!(G.adj[u], findfirst(==(v), G.adj[u]))
    u != v && deleteat!(G.adj[v], findfirst(==(u), G.adj[v]))
    return G
end

has_edge(G::PyGraph, u, v) = haskey(G.adj, u) && v in G.adj[u]
degree(G::PyGraph, v) = length(G.adj[v])
neighbors(G::PyGraph, v) = G.adj[v]
number_of_nodes(G::PyGraph) = length(G.nodes)
number_of_edges(G::PyGraph) = sum(length, values(G.adj); init=0) ÷ 2

"""
    pyedges(G)

Edges in networkx's iteration order: by node, then adjacency, each edge once.
"""
function pyedges(G::PyGraph)
    seen = Set{Int}()
    out = Tuple{Int,Int}[]
    for n in G.nodes
        for nbr in G.adj[n]
            nbr in seen || push!(out, (n, nbr))
        end
        push!(seen, n)
    end
    return out
end

"""Sorted edge list with each edge as (min, max)."""
sorted_edges(G::PyGraph) = sort!([minmax(e...) for e in pyedges(G)])

"""
    pycopy(G)

`networkx.Graph.copy`: nodes in order, then every adjacency entry re-added as an edge, which
can reorder adjacency lists.
"""
function pycopy(G::PyGraph)
    H = PyGraph()
    foreach(v -> add_node!(H, v), G.nodes)
    for u in G.nodes, v in G.adj[u]
        add_edge!(H, u, v)
    end
    return H
end

"""
    shortest_path(G, s, t)

`networkx.bidirectional_shortest_path`.
"""
function shortest_path(G::PyGraph, s::Integer, t::Integer)
    s == t && return [s]
    pred = Dict{Int,Union{Int,Nothing}}(s => nothing)
    succ = Dict{Int,Union{Int,Nothing}}(t => nothing)
    fwd, rev = [s], [t]
    meet = nothing
    while !isempty(fwd) && !isempty(rev) && meet === nothing
        if length(fwd) <= length(rev)
            level, fwd = fwd, Int[]
            for v in level
                for w in G.adj[v]
                    if !haskey(pred, w)
                        push!(fwd, w)
                        pred[w] = v
                    end
                    haskey(succ, w) && (meet = w; break)
                end
                meet === nothing || break
            end
        else
            level, rev = rev, Int[]
            for v in level
                for w in G.adj[v]
                    if !haskey(succ, w)
                        succ[w] = v
                        push!(rev, w)
                    end
                    haskey(pred, w) && (meet = w; break)
                end
                meet === nothing || break
            end
        end
    end
    meet === nothing && throw(ArgumentError("no path between $s and $t"))
    path = Int[]
    w = meet
    while w !== nothing
        push!(path, w)
        w = pred[w]
    end
    reverse!(path)
    w = succ[path[end]]
    while w !== nothing
        push!(path, w)
        w = succ[w]
    end
    return path
end

shortest_path_length(G::PyGraph, s, t) = length(shortest_path(G, s, t)) - 1

function has_path(G::PyGraph, s, t)
    try
        shortest_path(G, s, t)
        return true
    catch e
        e isa ArgumentError || rethrow()
        return false
    end
end

"""
    bfs_tree(G, source)

BFS tree of `networkx.bfs_tree` (children in adjacency order), as an undirected `PyGraph`.
"""
function bfs_tree(G::PyGraph, source::Integer)
    T = PyGraph()
    add_node!(T, source)
    seen = Set([source])
    level = [source]
    while !isempty(level)
        nxt = Int[]
        for parent in level, child in G.adj[parent]
            if !(child in seen)
                push!(seen, child)
                push!(nxt, child)
                add_edge!(T, parent, child)
            end
        end
        level = nxt
    end
    return T
end

"""
    simplegraph(G)

`(Graphs.SimpleGraph, vertex labels)` of `G`, vertex i carrying label `G.nodes[i]`.
"""
function simplegraph(G::PyGraph)
    idx = Dict(v => i for (i, v) in enumerate(G.nodes))
    S = Graphs.SimpleGraph(length(G.nodes))
    for (u, v) in pyedges(G)
        Graphs.add_edge!(S, idx[u], idx[v])
    end
    return S, G.nodes
end

"""
    connected_components(G)

Components (vertex label sets) in networkx order: by first node in node order.
"""
function connected_components(G::PyGraph)
    S, labels = simplegraph(G)
    comps = [Set(labels[c]) for c in Graphs.connected_components(S)]
    first_pos = Dict(v => i for (i, v) in enumerate(G.nodes))
    return sort!(comps; by=c -> minimum(first_pos[v] for v in c))
end

is_connected(G::PyGraph) = !isempty(G.nodes) && Graphs.is_connected(first(simplegraph(G)))
