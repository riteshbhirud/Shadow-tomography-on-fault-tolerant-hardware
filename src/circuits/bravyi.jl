# The BB-code memory experiment of Bravyi et al. (arXiv:2308.07915), used to calibrate the noise
# model. Follows github.com/sbravyi/BivariateBicycleCodes at fa77e33:
# depth-8 cycle, their CNOT orders and noise placement, N_c noisy cycles then 2 noiseless ones.

const SX = [:idle, 1, 4, 3, 5, 0, 2]
const SZ = [3, 5, 0, 1, 2, 4, :idle]
const BRAVYI_PARAMS = Dict("bb72" => (6, 6, (3, 1, 2), (3, 1, 2)), "gross" => (12, 6, (3, 1, 2), (3, 1, 2)))

"""
    bb_matrices(name)

`(l, m, (A1, A2, A3), (B1, B2, B3))` in their construction.
"""
function bb_matrices(name::AbstractString="bb72")
    l, m, (a1, a2, a3), (b1, b2, b3) = BRAVYI_PARAMS[name]
    x(i) = monomial(l, m, i, 0)
    y(i) = monomial(l, m, 0, i)
    return l, m, (x(a1), y(a2), y(a3)), (y(b1), x(b2), x(b3))
end

"""
    memory_circuit(name, basis, cycles, p; final_noiseless_cycles=2)

Memory experiment in `basis`: 'Z' prepares |0>, uses Z-check detectors and observes all Zbar; 'X' likewise.
"""
function memory_circuit(name::AbstractString, basis::Char, cycles::Integer, p::Real; final_noiseless_cycles::Integer=2)
    l, m, As, Bs = bb_matrices(name)
    n2 = l * m
    n = 2 * n2
    # qubits: X checks 0..n2-1, left data, right data, Z checks 3n2..4n2-1
    xc(i) = i
    dl(i) = n2 + i
    dr(i) = 2 * n2 + i
    zc(i) = 3 * n2 + i
    data = collect(n2:3n2-1)
    nbs = Dict{Tuple{Char,Int,Int},Int}()
    col(M, i) = findfirst(M[i+1, :]) - 1
    row(M, i) = findfirst(M[:, i+1]) - 1
    for i in 0:n2-1
        for (d, M) in enumerate(As)
            nbs[('X', i, d - 1)] = dl(col(M, i))
            nbs[('Z', i, 3 + d - 1)] = dr(row(M, i))
        end
        for (d, M) in enumerate(Bs)
            nbs[('X', i, 3 + d - 1)] = dr(col(M, i))
            nbs[('Z', i, d - 1)] = dl(row(M, i))
        end
    end
    c = Circuit()
    function cycle(noisy)
        q = noisy ? Float64(p) : 0.0
        function cnots(pairs)
            flat = [t for pr in pairs for t in pr]
            append!(c, "CX", flat)
            q > 0 && append!(c, "DEPOLARIZE2", flat, [q])
        end
        idle(qs) = (q > 0 && !isempty(qs) && append!(c, "DEPOLARIZE1", sort!(collect(qs)), [q]); nothing)
        xs = [xc(i) for i in 0:n2-1]
        zs = [zc(i) for i in 0:n2-1]
        append!(c, "RX", xs)
        q > 0 && append!(c, "Z_ERROR", xs, [q])
        pairs = [(nbs[('Z', i, SZ[1])], zc(i)) for i in 0:n2-1]
        cnots(pairs)
        idle(setdiff(data, first.(pairs)))
        append!(c, "TICK")
        for t in 2:6
            cnots(vcat([(xc(i), nbs[('X', i, SX[t])]) for i in 0:n2-1], [(nbs[('Z', i, SZ[t])], zc(i)) for i in 0:n2-1]))
            append!(c, "TICK")
        end
        q > 0 && append!(c, "X_ERROR", zs, [q])
        append!(c, "M", zs)
        pairs = [(xc(i), nbs[('X', i, SX[7])]) for i in 0:n2-1]
        cnots(pairs)
        idle(setdiff(data, last.(pairs)))
        append!(c, "TICK")
        idle(data)
        q > 0 && append!(c, "Z_ERROR", xs, [q])
        append!(c, "MX", xs)
        append!(c, "R", zs)
        q > 0 && append!(c, "X_ERROR", zs, [q])
        append!(c, "TICK")
    end
    append!(c, "R", [zc(i) for i in 0:n2-1])  # Z-check ancillas start clean (frame = 0)
    append!(c, basis == 'Z' ? "R" : "RX", data)
    append!(c, "TICK")
    per_cycle = 2 * n2  # each cycle measures n2 Z checks, then n2 X checks
    for t in 0:cycles+final_noiseless_cycles-1
        cycle(t < cycles)
        for i in 0:n2-1
            off = basis == 'Z' ? (n2 - i) + n2 : (n2 - i)  # rec offset of check i this cycle
            if t == 0
                append!(c, "DETECTOR", [rec(-off)], [i, t])
            else
                append!(c, "DETECTOR", [rec(-off), rec(-off - per_cycle)], [i, t])
            end
        end
    end
    append!(c, basis == 'Z' ? "M" : "MX", data)
    for (j, L) in enumerate(_bravyi_logicals(name, basis))
        append!(c, "OBSERVABLE_INCLUDE", [rec(-(n - (q - 1))) for q in findall(L)], [j - 1])
    end
    return c
end

function _bravyi_logicals(name, basis)
    _, _, As, Bs = bb_matrices(name)
    A = reduce(.⊻, As)
    B = reduce(.⊻, Bs)
    hx = Matrix{Bool}([A B])
    hz = Matrix{Bool}([permutedims(B) permutedims(A)])
    K, S = basis == 'Z' ? (nullspace_gf2(hx), hz) : (nullspace_gf2(hz), hx)
    cur = RowSpace(S)
    return [K[i, :] for i in axes(K, 1) if addrow!(cur, K[i, :])]
end
