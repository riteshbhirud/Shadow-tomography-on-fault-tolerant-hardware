# Logical-level helpers. Phases are tracked with QuantumClifford's PauliOperator,
# whose convention is stim's (phase byte p means i^p; X*Z = -iY).

function _qc_pauli(d, n::Integer)
    xs, zs = falses(n), falses(n)
    for (q, s) in d
        s in ('X', 'Y') && (xs[q+1] = true)
        s in ('Z', 'Y') && (zs[q+1] = true)
    end
    return PauliOperator(0x0, xs, zs)
end

"""
    logical_to_physical(label, basis, n)

`(Pauli dict, sign bit)` of a logical label such as "XIY" ('I' or '_' is identity), with
Ybar_i = i Xbar_i Zbar_i. Throws if the product is not Hermitian.
"""
function logical_to_physical(label::AbstractString, basis, n::Integer)
    length(label) == length(basis) ||
        throw(ArgumentError("label $(repr(label)) has length $(length(label)), expected $(length(basis))"))
    acc = _qc_pauli((), n)
    for (ch, (X, Z)) in zip(label, basis)
        ch in ('I', '_') && continue
        op = ch == 'X' ? _qc_pauli(X, n) : ch == 'Z' ? _qc_pauli(Z, n) :
             ch == 'Y' ? im * _qc_pauli(X, n) * _qc_pauli(Z, n) :
             throw(ArgumentError("bad character $(repr(ch)) in $(repr(label))"))
        acc = acc * op
    end
    ph = acc.phase[] & 0x3
    ph in (0x0, 0x2) || throw(ArgumentError("$(repr(label)) maps to a non-Hermitian operator"))
    d = PauliDict()
    for q in 0:n-1
        x, z = acc[q+1]
        (x || z) && (d[q] = x && z ? 'Y' : (x ? 'X' : 'Z'))
    end
    return d, ph == 0x0 ? 0 : 1
end

"""
    physical_state(gens, basis, n)

Map logical generators `[(label, sign bit)]` to `[(Pauli dict, sign bit)]`.
"""
function physical_state(gens, basis, n::Integer)
    out = StateGens()
    for (lab, s) in gens
        d, sb = logical_to_physical(lab, basis, n)
        push!(out, (d, Int(s) ⊻ sb))
    end
    return out
end

function _label_pauli(lab::AbstractString, sign::Integer=0)
    k = length(lab)
    xs, zs = falses(k), falses(k)
    for (i, ch) in enumerate(lab)
        ch in ('X', 'Y') && (xs[i] = true)
        ch in ('Z', 'Y') && (zs[i] = true)
    end
    return PauliOperator(sign == 0 ? 0x0 : 0x2, xs, zs)
end

"""
    exact_expectation(gens, label)

⟨P⟩ ∈ {+1, -1, 0} on the stabilizer state with generators `[(label, sign bit)]` (QuantumClifford `expect`).
"""
exact_expectation(gens, label::AbstractString) =
    Int(expect(_label_pauli(replace(label, '_' => 'I')), Stabilizer([_label_pauli(l, s) for (l, s) in gens])))

"""
    bb72_readout_experiment(kind, state_labels)

Destructive readout of every BB72 data qubit in `kind` ('X' or 'Z'), measuring that logical on all
12 qubits of the oracle basis. Returns `(experiment, labels)`.
"""
function bb72_readout_experiment(kind::Char, state_labels; name="bb72")
    s = bb_setup()
    code, basis = s.code, s.basis
    n, k = code.n, length(basis)
    targets, labels = PauliDict[], String[]
    for i in 1:k
        lab = join(j == i ? kind : 'I' for j in 1:k)
        d, sb = logical_to_physical(lab, basis, n)
        @assert sb == 0 && Set(values(d)) == Set([kind])
        push!(targets, d)
        push!(labels, lab)
    end
    exp = Experiment(; n_data=n, active0=0:n-1, prep_checks=check_dicts(code),
                     logical_state=physical_state(state_labels, basis, n), epochs=Epoch[],
                     readout=OrderedDict(q => kind for q in 0:n-1), targets, name="$(name)_readout_$kind",
                     logical_basis=basis)
    return exp, labels
end

"""
    canonical_stabilizers(stab)

stim's `canonical_stabilizers`: the generators in reduced row echelon form, columns taken qubit by
qubit (X before Z). The result depends only on the stabilizer group, so it equals stim's output.
"""
function canonical_stabilizers(stab::Stabilizer)
    rows = [stab[i] for i in 1:length(stab)]
    n = nqubits(stab)
    piv0 = 1
    for q in 1:n, b in (1, 2)
        bit(P) = b == 1 ? P[q][1] : P[q][2]
        piv = findfirst(i -> bit(rows[i]), piv0:length(rows))
        piv === nothing && continue
        piv += piv0 - 1
        for s in eachindex(rows)
            s != piv && bit(rows[s]) && (rows[s] = rows[s] * rows[piv])
        end
        rows[piv0], rows[piv] = rows[piv], rows[piv0]
        piv0 += 1
    end
    return rows
end

_pauli_label(P::PauliOperator) = join((x && z) ? 'Y' : x ? 'X' : z ? 'Z' : 'I' for (x, z) in (P[i] for i in 1:nqubits(P)))

"""
    random_stabilizer_state(k, seed; depth=40k²)

Canonical generators `[(label, sign bit)]` of a seeded random k-qubit stabilizer state:
`depth` random H / S / CX / X-or-Z gates on |0…0⟩. The gate draws are numpy's `default_rng(seed)`
stream through the bridge; the gates act on a QuantumClifford tableau.
"""
function random_stabilizer_state(k::Integer, seed::Integer; depth = nothing)
    rng = pyimport("numpy").random.default_rng(seed)
    draw(m) = pyconvert(Int, rng.integers(m))
    st = one(Stabilizer, k)
    for _ in 1:(depth === nothing || depth == 0 ? 40 * k * k : depth)
        g = draw(4)
        q = draw(k) + 1
        if g == 0
            apply!(st, sHadamard(q))
        elseif g == 1
            apply!(st, sPhase(q))
        elseif g == 2 && k > 1
            r = draw(k - 1)
            apply!(st, sCNOT(q, r + (r >= q - 1) + 1))
        else
            apply!(st, draw(2) == 1 ? sX(q) : sZ(q))
        end
    end
    return [(_pauli_label(P), P.phase[] == 0x0 ? 0 : 1) for P in canonical_stabilizers(st)]
end

"""
    random_css_state(m, seed; kz=m÷2)

Seeded random CSS state on m qubits: kz Z-type generators (rows of a random full-rank A),
X-type generators spanning ker(A), random signs, all drawn from numpy's `default_rng(seed)`.
"""
function random_css_state(m::Integer, seed::Integer; kz = nothing)
    rng = pyimport("numpy").random.default_rng(seed)
    kz = kz === nothing ? m ÷ 2 : kz
    local A
    while true
        A = pyconvert(Matrix{Int}, rng.integers(0, 2; size = (kz, m)))
        gf2rank(A) == kz && break
    end
    K = nullspace_gf2(A)
    gens = [(join(b == 1 ? 'Z' : 'I' for b in A[i, :]), pyconvert(Int, rng.integers(2))) for i in axes(A, 1)]
    append!(gens, [(join(b ? 'X' : 'I' for b in K[i, :]), pyconvert(Int, rng.integers(2))) for i in axes(K, 1)])
    @assert length(gens) == m
    return gens
end
