# La-cross [[65,9,4]] feasibility scan (Di Fini et al. arXiv:2610.12685, Pecorari et al.
# arXiv:2404.13010): the code and its folds, a bare-ancilla memory circuit, and a BP-OSD memory Monte Carlo.
# The memory circuit is built in our circuit IR, so its stim text equals the reference one; stim (DEM, sampling,
# distance search) and ldpc's BP-OSD are bridged as everywhere else.

using LinearAlgebra: I
using OrderedCollections: OrderedDict
using QuantumClifford: Stabilizer, apply!, stab_to_gf2, sHadamard, sSWAP, sPhase, sInvPhase, sCPHASE
using SparseArrays: SparseMatrixCSC
import JSON

const LACROSS_K, LACROSS_N = 3, 7  # [[65,9,4]]: k=3, n=7, h_3(x) = 1 + x + x^3 (2610.12685 Sec. 2.1)

"""G = (I_k | A) of the seed code, columns completed by the recurrence c_{a+k} = c_a + c_{a+1}."""
function lacross_seed_G(k::Integer=LACROSS_K, n::Integer=LACROSS_N)
    G = zeros(Bool, k, n)
    G[:, 1:k] .= I(k)
    for a in 1:n-k
        G[:, a+k] .= G[:, a] .⊻ G[:, a+1]
    end
    return G
end

"""
    LaCross(k=3, n=7)

HGP(H, H): H_X = [H (x) I_n | I_m (x) H^T], H_Z = [I_n (x) H | H^T (x) I_m]. Qubit indices are 0-based:
principal (i, j) -> i*n + j, check sector [a, b] -> n^2 + a*m + b. `S` is the symplectic
stabilizer matrix [H_X | 0; 0 | H_Z].
"""
struct LaCross
    k_seed::Int
    n_seed::Int
    m_seed::Int
    H::Matrix{Bool}
    G::Matrix{Bool}
    HX::Matrix{Bool}
    HZ::Matrix{Bool}
    N::Int
    S::Matrix{Bool}
end

function LaCross(k::Integer=LACROSS_K, n::Integer=LACROSS_N)
    H, m = TaskB.lacross_seed(k, n), n - k
    In, Im = I(n), I(m)
    HX = Matrix{Bool}(isodd.(hcat(kron(H, In), kron(Im, permutedims(H)))))
    HZ = Matrix{Bool}(isodd.(hcat(kron(In, H), kron(permutedims(H), Im))))
    N = n * n + m * m
    S = [HX zeros(Bool, size(HX)); zeros(Bool, size(HZ)) HZ]
    return LaCross(k, n, m, H, lacross_seed_G(k, n), HX, HZ, N, S)
end

lc_P(c::LaCross, i, j) = i * c.n_seed + j
lc_C(c::LaCross, a, b) = c.n_seed^2 + a * c.m_seed + b
lc_c_u(c::LaCross, u) = isodd.(Int.(permutedims(c.G)) * Int.(u))

"""
    canonical_logicals(code) -> OrderedDict((i, j) => (x, z))

Eq. eq:canonical-reps: Xbar_ij ~ X(eps_i c_{e_j}^T, 0) (row i), Zbar_ij ~ Z(c_{e_i} eps_j^T, 0) (column j), as
length-N bit vectors, keys 0-based in (i, j) order.
"""
function canonical_logicals(c::LaCross)
    k = c.k_seed
    e(j) = [t == j for t in 0:k-1]
    out = OrderedDict{Tuple{Int,Int},Tuple{Vector{Bool},Vector{Bool}}}()
    for i in 0:k-1, j in 0:k-1
        x, z = zeros(Bool, c.N), zeros(Bool, c.N)
        for col in findall(lc_c_u(c, e(j))) .- 1
            x[lc_P(c, i, col)+1] = true
        end
        for r in findall(lc_c_u(c, e(i))) .- 1
            z[lc_P(c, r, j)+1] = true
        end
        out[(i, j)] = (x, z)
    end
    return out
end

"""Fold involution tau (0-based): (i,j) <-> (j,i), [a,b] <-> [b,a]."""
function transpose_perm(c::LaCross)
    perm = zeros(Int, c.N)
    for i in 0:c.n_seed-1, j in 0:c.n_seed-1
        perm[lc_P(c, i, j)+1] = lc_P(c, j, i)
    end
    for a in 0:c.m_seed-1, b in 0:c.m_seed-1
        perm[lc_C(c, a, b)+1] = lc_C(c, b, a)
    end
    return perm
end

# all w-subsets of 0:n-1 in itertools.combinations order, calling f on each until it returns true
function _each_combination(f, n::Integer, w::Integer)
    comb = collect(0:w-1)
    while true
        f(comb) && return true
        i = w
        while i >= 1 && comb[i] == n - w + i - 1
            i -= 1
        end
        i == 0 && return false
        comb[i] += 1
        for j in i+1:w
            comb[j] = comb[j-1] + 1
        end
    end
end

"""
    exhaustive_min_logical(Hcheck, Hstab, wmax=3) -> (w, support, count_checked)

Smallest-weight v with Hcheck v = 0 and v not in rowspace(Hstab), searched exhaustively up to weight wmax
(`(nothing, nothing, count)` if none, so the exact lower bound is wmax + 1). `support` is 0-based.
"""
function exhaustive_min_logical(Hcheck::AbstractMatrix, Hstab::AbstractMatrix, wmax::Integer=3)
    n = size(Hcheck, 2)
    cols = [BitVector(Hcheck[:, q]) for q in 1:n]
    rS = Circuits.gf2rank(Hstab)
    count_checked = 0
    found = nothing
    for w in 1:wmax
        _each_combination(n, w) do comb
            count_checked += 1
            s = falses(size(Hcheck, 1))
            for q in comb
                s .⊻= cols[q+1]
            end
            any(s) && return false
            v = zeros(Bool, n)
            v[comb.+1] .= true
            Circuits.gf2rank(vcat(Hstab, permutedims(v))) > rS || return false
            found = (w, Tuple(comb))
            return true
        end
        found === nothing || return found[1], found[2], count_checked
    end
    return nothing, nothing, count_checked
end

# the invertible k x k matrices over GF(2), in the reference bit order
function _gl_k2(k::Integer)
    out = Matrix{Bool}[]
    for bits in 0:(1<<(k*k))-1
        A = permutedims(reshape([isodd(bits >> t) for t in 0:k*k-1], k, k))
        Circuits.gf2rank(A) == k && push!(out, A)
    end
    return out
end

"""|Gamma_{k,n}| = #{A in GL(k,2): mu_n(A v) = mu_n(v) for all v} (Eqs. eq:gamma, eq:multiplicity)."""
function gamma_order(c::LaCross)
    k = c.k_seed
    mu = Dict{Vector{Bool},Int}()
    for j in 1:c.n_seed
        col = c.G[:, j]
        mu[col] = get(mu, col, 0) + 1
    end
    vs = [Bool[isodd(b >> (k - 1 - t)) for t in 0:k-1] for b in 0:(1<<k)-1]
    return count(A -> all(v -> get(mu, v, 0) == get(mu, Vector{Bool}(isodd.(Int.(A) * Int.(v))), 0), vs), _gl_k2(k))
end

"""
    cnot_count_gauss(Q) -> (#CNOT, #row swaps)

CNOT count of a naive Gaussian-elimination synthesis of the reversible linear map Q (row operations only; an
UPPER BOUND, not optimal). A missing pivot is fixed with a CNOT, not a swap.
"""
function cnot_count_gauss(Q::AbstractMatrix)
    A = Matrix{Bool}(isodd.(Q))
    m = size(A, 1)
    cn = sw = 0
    for c in 1:m
        if !A[c, c]
            r = findfirst(r -> A[r, c], c+1:m)
            r === nothing && error("cnot_count_gauss: singular matrix")
            A[c, :] .⊻= A[c+r, :]
            cn += 1
        end
        for r in 1:m
            if r != c && A[r, c]
                A[r, :] .⊻= A[c, :]
                cn += 1
            end
        end
    end
    @assert A == I(m)
    return cn, sw
end

# x mod 2 of the rounded float matrix, as numpy's np.rint(...) % 2 / np.round(...) % 2
_round_mod2(M) = mod.(round.(Int, M), 2)

"""
    lacross_seed_automorphisms(code)

All coordinate permutations P of the seed with ker(H P) = ker H; for each, Q with H P = Q H. Reports how many Q
are permutations and the naive Gaussian-elimination CNOT counts of Q and of Q^{-T} (the check-sector map).
"""
function lacross_seed_automorphisms(c::LaCross)
    H, G = Int.(c.H), Int.(c.G)
    n, m = c.n_seed, c.m_seed
    rows = []
    Hinv = inv(Float64.(H[:, 1:m]))  # H's first m columns are upper unitriangular
    for perm in TaskB._permutations(n)
        Pm = zeros(Int, n, n)
        for i in 1:n
            Pm[i, perm[i]] = 1
        end
        GP = mod.(G * Pm, 2)
        any(isodd, H * permutedims(GP)) && continue
        HP = mod.(H * Pm, 2)
        Q = _round_mod2(HP[:, 1:m] * Hinv)  # Q H[:, :m] = HP[:, :m]
        @assert !any(isodd, Q * H .- HP)
        is_perm = all(==(1), sum(Q; dims=1)) && all(==(1), sum(Q; dims=2))
        push!(rows, (perm=perm .- 1, Q=Q, is_perm=is_perm, cnot_Q=cnot_count_gauss(Q)[1]))
    end
    cn = [r.cnot_Q for r in rows]
    cnT = [cnot_count_gauss(permutedims(_round_mod2(inv(Float64.(r.Q)))))[1] for r in rows]
    stats(v) = OrderedDict{String,Any}("min" => minimum(v), "max" => maximum(v), "mean" => np_mean(Float64.(v)))
    return OrderedDict{String,Any}(
        "n_seed_automorphisms" => length(rows),
        "n_with_Q_permutation" => count(r -> r.is_perm, rows),
        "cnot_Q_naive_gauss" => stats(cn),
        "cnot_QinvT_naive_gauss" => stats(cnT),
        "check_sector_map_cnot_upper_bound_note" =>
            "Q1^{-T} (x) Q2 = (Q1^{-T} (x) I_m)(I_m (x) Q2): naive count <= m*c(Q1^{-T}) + m*c(Q2) " *
            "= $m*c1 + $m*c2 on the $(m)x$m check sector; non-permutation Q only; NOT optimized, NOT simulated",
        "worst_case_check_sector_cnots_upper_bound" => m * maximum(cnT) + m * maximum(cn))
end

# images (x|z bits, signs dropped) of the Pauli rows under a list of QuantumClifford gates
function _pauli_images(rows::AbstractMatrix{Bool}, N::Integer, ops)
    st = Stabilizer(zeros(UInt8, size(rows, 1)), Matrix{Bool}(rows[:, 1:N]), Matrix{Bool}(rows[:, N+1:2N]))
    foreach(op -> apply!(st, op), ops)
    return stab_to_gf2(st)
end

"""
    fold_checks(code)

The Hadamard fold (transversal H, then SWAP of mirror pairs) and the phase fold (S on the principal diagonal,
S_DAG on the check diagonal, CZ on mirror pairs): do they preserve the stabilizer group, and is their logical
action the paper's, up to stabilizers and sign? Done on QuantumClifford tableaux.
"""
function fold_checks(c::LaCross)
    N, k, n, m = c.N, c.k_seed, c.n_seed, c.m_seed
    tau = transpose_perm(c)
    pairs = [(q, tau[q+1]) for q in 0:N-1 if tau[q+1] > q]
    hf = vcat([sHadamard(q + 1) for q in 0:N-1], [sSWAP(a + 1, b + 1) for (a, b) in pairs])
    pf = vcat(Any[sPhase(lc_P(c, i, i) + 1) for i in 0:n-1], Any[sInvPhase(lc_C(c, a, a) + 1) for a in 0:m-1],
              Any[sCPHASE(a + 1, b + 1) for (a, b) in pairs])
    S = c.S
    rS = Circuits.gf2rank(S)
    logs = canonical_logicals(c)
    in_stab(v) = Circuits.gf2rank(vcat(S, permutedims(v))) == rS
    function Xbar(M)
        v = zeros(Bool, 2N)
        for ((i, j), (x, _)) in logs
            M[i+1, j+1] && (v[1:N] .⊻= x)
        end
        return v
    end
    function Zbar(M)
        v = zeros(Bool, 2N)
        for ((i, j), (_, z)) in logs
            M[i+1, j+1] && (v[N+1:end] .⊻= z)
        end
        return v
    end
    res = OrderedDict{String,Any}()
    for (name, ops) in (("hadamard_fold", hf), ("phase_fold", pf))
        stab_ok = all(in_stab, eachrow(_pauli_images(S, N, ops)))
        act_ok = true
        for i in 1:k, j in 1:k
            E = zeros(Bool, k, k)
            E[i, j] = true
            img = _pauli_images(permutedims(hcat(Xbar(E), Zbar(E))), N, ops)
            if name == "hadamard_fold"  # (M, N) -> (N^T, M^T): Xbar(M) -> Zbar(M^T), Zbar(N) -> Xbar(N^T)
                expX, expZ = Zbar(permutedims(E)), Xbar(permutedims(E))
            else  # (M, N) -> (M, N + M^T): Xbar(M) -> Xbar(M) Zbar(M^T), Zbar unchanged
                expX, expZ = Xbar(E) .⊻ Zbar(permutedims(E)), Zbar(E)
            end
            act_ok &= in_stab(img[1, :] .⊻ expX) && in_stab(img[2, :] .⊻ expZ)
        end
        res[name] = OrderedDict{String,Any}(
            "preserves_stabilizer_group" => stab_ok,
            "logical_action_matches_paper_up_to_stabilizers_and_sign" => act_ok,
            "n_two_qubit_gates" => length(pairs),
            "n_fixed_qubits" => count(q -> tau[q+1] == q, 0:N-1))
    end
    res["hadamard_fold"]["note"] = "transversal H + SWAP of mirror pairs (SWAPs assumed relabelling/atom moves)"
    res["phase_fold"]["note"] = "S on principal diag, S_DAG on check diag, CZ on mirror pairs (single layer)"
    return res
end

"""
    min_weight_codeword(G_logical, G_trivial; iters=2000, seed=0) -> (weight, vector)

The information-set search of `Circuits.min_weight_codeword` on a Bool copy of
`G_logical`: an upper bound on the weight of the lightest vector in its span outside span(G_trivial).
"""
min_weight_codeword(Glog::AbstractMatrix, Gtriv::AbstractMatrix; iters::Integer=2000, seed::Integer=0) =
    Circuits.min_weight_codeword(Matrix{Bool}(Glog), Gtriv; iters, seed)

"""
    lacross_code_record(; t0=time()) -> (record, arrays)

The `code.json` record of the `code` step and the arrays of `lacross_65_9_4_checks.npz` (HX, HZ, H, G).
"""
function lacross_code_record(; t0::Float64=time())
    c = LaCross()
    H, G = c.H, c.G
    @assert !any(isodd, Int.(H) * permutedims(Int.(G)))
    @assert !any(isodd, Int.(c.HX) * permutedims(Int.(c.HZ)))
    rx, rz = Circuits.gf2rank(c.HX), Circuits.gf2rank(c.HZ)
    K = c.N - rx - rz
    # codewords u G for u in itertools.product((0, 1), repeat=k)
    cw = [isodd.(permutedims(Int.(G)) * [isodd(b >> (c.k_seed - 1 - t)) for t in 0:c.k_seed-1]) for b in 0:(1<<c.k_seed)-1]
    seed_weights = sort([count(x) for x in cw if any(x)])
    seed_T_rank = Circuits.gf2rank(permutedims(H))
    rowlists(M) = [Int.(M[i, :]) for i in axes(M, 1)]
    out = OrderedDict{String,Any}(
        "construction" => "2610.12685 Eq. eq:seed (H_{k,n} rows e_a+e_{a+1}+e_{a+k}, open boundary), " *
                          "Eq. eq:checks (H_X=[H(x)I_n | I_m(x)H^T], H_Z=[I_n(x)H | H^T(x)I_m]); k=3, n=7, h_3=1+x+x^3",
        "seed_H" => rowlists(H), "seed_G" => rowlists(G),
        "seed_codeword_weights" => seed_weights, "seed_distance" => minimum(seed_weights),
        "seed_rank_H" => Circuits.gf2rank(H), "seed_ker_HT_dim" => size(H, 1) - seed_T_rank,
        "N" => c.N, "n_X_checks" => size(c.HX, 1), "n_Z_checks" => size(c.HZ, 1),
        "rank_HX" => rx, "rank_HZ" => rz, "K" => K,
        "max_check_weight" => max(maximum(sum(c.HX; dims=2)), maximum(sum(c.HZ; dims=2))),
        "max_qubit_degree_X" => maximum(sum(c.HX; dims=1)), "max_qubit_degree_Z" => maximum(sum(c.HZ; dims=1)))
    dX = exhaustive_min_logical(c.HZ, c.HX, 3)
    dZ = exhaustive_min_logical(c.HX, c.HZ, 3)
    logs = canonical_logicals(c)
    wX = sort(unique(count(x) for (x, _) in values(logs)))
    wZ = sort(unique(count(z) for (_, z) in values(logs)))
    ok = true
    for (a, (x, _)) in logs, (b, (_, z)) in logs
        ok &= isodd(count(x .& z)) == (a == b)
    end
    for (x, z) in values(logs)
        ok &= !any(isodd, Int.(c.HZ) * Int.(x)) && !any(isodd, Int.(c.HX) * Int.(z))
    end
    out["canonical_logicals_valid_symplectic_basis"] = ok
    out["canonical_logical_weights"] = OrderedDict{String,Any}("X" => wX, "Z" => wZ)
    out["distance_exhaustive"] = OrderedDict{String,Any}(
        "X_type_none_of_weight_le3" => dX[1] === nothing, "Z_type_none_of_weight_le3" => dZ[1] === nothing,
        "n_supports_checked_each" => dX[3], "weight4_logical_exists" => "yes (canonical representatives, weight 4)",
        "conclusion" => "d = 4 exactly (exhaustive lower bound 4, explicit weight-4 logical)")
    wisdX, _ = min_weight_codeword(Circuits.nullspace_gf2(c.HZ), c.HX; iters=200, seed=0)
    wisdZ, _ = min_weight_codeword(Circuits.nullspace_gf2(c.HX), c.HZ; iters=200, seed=1)
    out["distance_isd_upper_bound"] = OrderedDict{String,Any}("X" => wisdX, "Z" => wisdZ, "iters" => 200, "seeds" => [0, 1])
    checks = vcat([Circuits.PauliDict(q - 1 => 'X' for q in findall(c.HX[i, :])) for i in axes(c.HX, 1)],
                  [Circuits.PauliDict(q - 1 => 'Z' for q in findall(c.HZ[i, :])) for i in axes(c.HZ, 1)])
    t1 = time()
    dstim, nl = Circuits.code_distance_stim(c.N, checks; max_dets=4, max_degree=6)
    out["distance_stim_code_capacity"] = OrderedDict{String,Any}("value" => dstim, "n_logical_pairs" => nl,
        "max_dets" => 4, "max_degree" => 6, "seconds" => py_round(time() - t1, 1))
    out["gamma_3_7_order"] = gamma_order(c)
    out["folds"] = fold_checks(c)
    out["seed_automorphisms"] = lacross_seed_automorphisms(c)
    out["seconds"] = py_round(time() - t0, 1)
    arrays = ["HX" => c.HX, "HZ" => c.HZ, "H" => H, "G" => G]
    return out, arrays
end

"""
    lacross_memory(code, p, rounds; basis="Z", noise="uniform", x_order=0:5, z_order=0:5) -> Circuits.Circuit

Bare-ancilla memory circuit (b3_lacross.build_memory). One ancilla per check; each round resets the ancillas,
runs 6 CNOT layers for all X checks, then 6 for all Z checks, and measures the ancillas. Layer t of a block is
slot t of every check (principal offsets 0, 1, k; check-sector offsets 0, 1, k).

`noise="uniform"`: DEPOLARIZE2(p) after every CNOT, DEPOLARIZE1(p) on qubits idle in a CNOT layer and on data
during measure/reset, X/Z flips p after every reset and before every measurement. `noise="agnostic"`: the
hardware-agnostic model of 2404.13010 as described in its text (DEPOLARIZE2 after CNOTs, DEPOLARIZE1 after
ancilla H, no idle, no SPAM); their CNOT schedule is not reproduced.
"""
function lacross_memory(c::LaCross, p::Real, rounds::Integer; basis::AbstractString="Z", noise::AbstractString="uniform",
                        x_order=0:5, z_order=0:5)
    n, m, k, N = c.n_seed, c.m_seed, c.k_seed, c.N
    offs = (0, 1, k)
    nX, nZ = m * n, n * m
    ancX = [N + r for r in 0:nX-1]
    ancZ = [N + nX + r for r in 0:nZ-1]
    allq = 0:N+nX+nZ-1
    X_slots, Z_slots = [Tuple{Int,Int}[] for _ in 1:6], [Tuple{Int,Int}[] for _ in 1:6]
    for a in 0:m-1, j in 0:n-1
        anc = ancX[a*n+j+1]
        for (t, o) in enumerate(offs)
            push!(X_slots[t], (anc, lc_P(c, a + o, j)))
            b = j - o
            0 <= b < m && push!(X_slots[3+t], (anc, lc_C(c, a, b)))
        end
    end
    for i in 0:n-1, b in 0:m-1
        anc = ancZ[i*m+b+1]
        for (t, o) in enumerate(offs)
            push!(Z_slots[t], (lc_P(c, i, b + o), anc))
            a = i - o
            0 <= a < m && push!(Z_slots[3+t], (lc_C(c, a, b), anc))
        end
    end
    # every slot is a matching, and the slots together give back H_X and H_Z
    HXc, HZc = zeros(Bool, size(c.HX)), zeros(Bool, size(c.HZ))
    for sl in X_slots
        used = [q for pr in sl for q in pr]
        @assert allunique(used)
        foreach(((anc, d),) -> HXc[anc-N+1, d+1] ⊻= true, sl)
    end
    for sl in Z_slots
        used = [q for pr in sl for q in pr]
        @assert allunique(used)
        foreach(((d, anc),) -> HZc[anc-N-nX+1, d+1] ⊻= true, sl)
    end
    @assert HXc == c.HX && HZc == c.HZ

    uni = noise == "uniform"
    C = Circuits.Circuit()
    ap!(name, ts, args=Float64[]) = append!(C, name, ts, args)
    for q in 0:N-1
        ap!("QUBIT_COORDS", [q], [q])
    end
    data = collect(0:N-1)
    if basis == "Z"
        ap!("R", data)
        uni && ap!("X_ERROR", data, [p])
    else
        ap!("RX", data)
        uni && ap!("Z_ERROR", data, [p])
    end
    ap!("TICK", Int[])
    function cnot_layer(pairs)
        flat = [q for pr in pairs for q in pr]
        ap!("CX", flat)
        ap!("DEPOLARIZE2", flat, [p])
        if uni
            busy = Set(flat)
            ap!("DEPOLARIZE1", [q for q in allq if !(q in busy)], [p])
        end
        ap!("TICK", Int[])
    end
    nanc = nX + nZ
    for r in 0:rounds-1
        if uni
            ap!("RX", ancX); ap!("Z_ERROR", ancX, [p])
            ap!("R", ancZ); ap!("X_ERROR", ancZ, [p])
            r > 0 && ap!("DEPOLARIZE1", data, [p])
        else
            ap!("R", vcat(ancX, ancZ)); ap!("H", ancX); ap!("DEPOLARIZE1", ancX, [p])
        end
        ap!("TICK", Int[])
        foreach(t -> cnot_layer(X_slots[t+1]), x_order)
        foreach(t -> cnot_layer(Z_slots[t+1]), z_order)
        if uni
            ap!("Z_ERROR", ancX, [p]); ap!("X_ERROR", ancZ, [p]); ap!("DEPOLARIZE1", data, [p])
            ap!("MX", ancX); ap!("M", ancZ)
        else
            ap!("H", ancX); ap!("DEPOLARIZE1", ancX, [p]); ap!("M", vcat(ancX, ancZ))
        end
        ap!("TICK", Int[])
        # detectors: measurement record order ancX (nX) then ancZ (nZ)
        same = basis == "Z" ? (nX:nX+nZ-1) : (0:nX-1)
        for s in 0:nanc-1
            if r == 0
                s in same && ap!("DETECTOR", [Circuits.rec(s - nanc)], [s, r])
            else
                ap!("DETECTOR", [Circuits.rec(s - nanc), Circuits.rec(s - 2nanc)], [s, r])
            end
        end
    end
    if basis == "Z"
        uni && ap!("X_ERROR", data, [p])
        ap!("M", data)
        Hs, off = c.HZ, nX
    else
        uni && ap!("Z_ERROR", data, [p])
        ap!("MX", data)
        Hs, off = c.HX, 0
    end
    for s in 0:size(Hs, 1)-1
        tg = [Circuits.rec(q - 1 - N) for q in findall(Hs[s+1, :])]
        push!(tg, Circuits.rec(off + s - nanc - N))
        ap!("DETECTOR", tg, [off + s, rounds])
    end
    for (idx, ((i, j), (x, z))) in enumerate(sort(collect(canonical_logicals(c)); by=first))
        v = basis == "Z" ? z : x
        ap!("OBSERVABLE_INCLUDE", [Circuits.rec(q - 1 - N) for q in findall(v)], [idx - 1])
    end
    return C
end

"""Noiseless determinism (stim builds a DEM without gauge detectors) and fires in 2000 noiseless shots."""
function lacross_noiseless_check(circ::Py; shots::Integer=2000)
    clean = circ.without_noise()
    clean.detector_error_model(; allow_gauge_detectors=false)
    res = clean.compile_detector_sampler(; seed=0).sample(shots; separate_observables=true)
    return pyconvert(Int, res[0].sum()), pyconvert(Int, res[1].sum())
end

"""
    lacross_circuit_record(; bases="Z", max_dets="3,4", max_degree=6, save=false, out=nothing)

The `circuit.json` record of the `circuit` step: structure and noiseless checks of the r = 4, p = 1e-3
memory circuits, then stim's bounded circuit-distance search. With `out`, the circuits are saved as in the
reference (uniform Z always, all four with `save`).
"""
function lacross_circuit_record(; bases::AbstractString="Z", max_dets::AbstractString="3,4", max_degree::Integer=6,
                                save::Bool=false, out=nothing)
    c = LaCross()
    rec = OrderedDict{String,Any}()
    for noise in ("uniform", "agnostic"), basis in ("Z", "X")
        circ = lacross_memory(c, 1e-3, 4; basis, noise)
        py = Circuits.pycircuit(circ)
        nd, no = lacross_noiseless_check(py)
        dem = py.detector_error_model(; decompose_errors=false)
        key = "$(noise)_$basis"
        rec[key] = OrderedDict{String,Any}(
            "num_qubits" => pyconvert(Int, py.num_qubits), "num_detectors" => pyconvert(Int, py.num_detectors),
            "num_observables" => pyconvert(Int, py.num_observables), "dem_errors" => pyconvert(Int, dem.num_errors),
            "noiseless_detector_fires" => nd, "noiseless_observable_flips" => no,
            "cnot_count" => sum((length(i.targets) ÷ 2 for i in circ.ops if i.name == "CX"); init=0),
            "depth_ticks" => count(i -> i.name == "TICK", circ.ops))
        if out !== nothing && (save || (noise == "uniform" && basis == "Z"))
            write(joinpath(out, "lacross65_memory_$(key)_r4_p1e-3.stim"), Circuits.stim_text(circ) * "\n")
        end
        println(key, " ", pyjson(rec[key]))
    end
    res = OrderedDict{String,Any}[]
    for basis in split(bases, ","), md in parse.(Int, split(max_dets, ","))
        circ = lacross_memory(c, 1e-3, 4; basis, noise="uniform")
        t0 = time()
        fw = Circuits.stim_search(circ; max_dets=md, max_degree)
        r = OrderedDict{String,Any}("basis" => String(basis), "noise" => "uniform", "rounds" => 4, "max_dets" => md,
                                    "max_degree" => max_degree, "found_weight" => fw, "seconds" => py_round(time() - t0, 1))
        println(pyjson(r))
        push!(res, r)
    end
    rec["circuit_distance_search"] = res
    return rec
end

const _LACROSS_DEC = Dict{String,Any}()

"""
    lacross_sim_batch(key, p, basis, noise, shots, seed, dec) -> (fails, shots)

One batch of the BP-OSD memory Monte Carlo (b3_lacross._worker): sample the r = 4 memory circuit with stim at
`seed`, decode every shot with ldpc's BpOsdDecoder (min-sum, OSD-CS) built from the DEM, and count shots where
any of the 9 logical observables is wrong. The decoder is cached per process and `key`.
"""
function lacross_sim_batch(key, p, basis, noise, shots, seed, dec)
    circ, d, L = get!(_LACROSS_DEC, key) do
        py = Circuits.pycircuit(lacross_memory(LaCross(), p, 4; basis, noise))
        H, L, pr = dem_to_matrices(read_dem(py.detector_error_model(; decompose_errors=false)))
        d = pymodule(:ldpc).BpOsdDecoder(_scipy_sparse(H, :csc); error_channel=pylist(pr),
            max_iter=dec["max_iter"], bp_method="minimum_sum", ms_scaling_factor=dec["ms"],
            osd_method="osd_cs", osd_order=dec["osd_order"])
        (py, d, Matrix{Int}(L))
    end
    res = circ.compile_detector_sampler(; seed).sample(shots; separate_observables=true)
    det, obs = pyconvert(Matrix{Bool}, res[0]), pyconvert(Matrix{Bool}, res[1])
    fails = 0
    for s in 1:shots
        if !any(view(det, s, :))
            fails += any(view(obs, s, :))
            continue
        end
        e = _from_numpy_u8(d.decode(to_numpy_u8(det[s, :])))
        pred = isodd.(L * Int.(e))
        fails += pred != view(obs, s, :)
    end
    return fails, shots
end

# Python's float ** (C pow), so the per-round rates are bitwise the reference ones

"""
    lacross_sim_record(key, p, basis, noise, fails, shots, dec, seed_base, seconds, workers)

The line the `sim` step appends to `sim_results.jsonl`.
"""
function lacross_sim_record(key, p, basis, noise, fails, shots, dec, seed_base, seconds, workers)
    R = 4
    pl = fails / shots
    lo, hi = wilson(fails, shots)
    per(x) = 1 - _cpow(1 - x, 1 / R)
    return OrderedDict{String,Any}(
        "key" => key, "p" => p, "basis" => basis, "noise" => noise, "rounds" => R, "shots" => shots, "fails" => fails,
        "p_L_total" => pl, "p_L_total_wilson95" => [lo, hi], "p_L_per_round" => per(pl),
        "p_L_per_round_wilson95" => [per(lo), per(hi)], "failure" => "any of the 9 logical observables wrong",
        "decoder" => OrderedDict{String,Any}("name" => "ldpc BpOsdDecoder", "bp_method" => "minimum_sum",
                                             "osd_method" => "osd_cs", "max_iter" => dec["max_iter"], "ms" => dec["ms"],
                                             "osd_order" => dec["osd_order"]),
        "seed_base" => seed_base, "seconds" => seconds, "workers" => workers)
end

"""
    lacross_summary(dir) -> OrderedDict

`lacross_summary.json` from `code.json`, `circuit.json` and `sim_results.jsonl` in `dir` (each optional), plus
the published comparison and the public-code search notes.
"""
function lacross_summary(dir::AbstractString)
    out = OrderedDict{String,Any}()
    for name in ("code", "circuit")
        fp = joinpath(dir, "$name.json")
        isfile(fp) && (out[name] = JSON.parsefile(fp))
    end
    fp = joinpath(dir, "sim_results.jsonl")
    out["memory"] = isfile(fp) ? [JSON.parse(l) for l in eachline(fp)] : Any[]
    pub(p, pld, err) = OrderedDict{String,Any}("p" => p, "pld" => pld, "pld_err" => err)
    out["published_comparison"] = OrderedDict{String,Any}(
        "source" => "Pecorari et al. Figshare doi:10.6084/m9.figshare.28074533.v2, QEC_data/k3_hardware_agnostic.dat " *
                    "(copied to results/b3/external_pecorari_figshare/)",
        "code" => "[[65,9,4]]", "rounds" => 4,
        "decoder" => "BP+OSD min-sum, scaling 1.0, osd_cs order 1 (file header)",
        "pld_definition" => "1-(1-pl)^(1/d), pl = decoder errors / shots",
        "rows" => [pub(0.0010, 0.00052, 0.000092), pub(0.0020, 0.002585, 0.000507)],
        "noise_model" => "hardware-agnostic: equal gate errors, no SPAM, no idle (2404.13010 text)",
        "caveat" => "memory basis / which observables count is not stated in the paper text")
    out["public_code_search"] = OrderedDict{String,Any}(
        "di_fini_2610_12685" => "no GitHub/Zenodo/code-availability statement in LaTeX source",
        "pecorari_2404_13010" => "data only (Figshare 28074533: 6 .dat files of logical error rates); no circuits or code")
    return out
end
