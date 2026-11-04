# Rotated surface-code patches and lattice surgery. Data qubit (row i, col j); the
# plaquette with corner (i, j) is X-type if i + j is even. Top/bottom boundaries keep weight-2 Z
# plaquettes, left/right weight-2 X plaquettes, so Zbar is a column and Xbar a row. Gate orders put
# each hook perpendicular to its logical; X checks are group 0, Z checks group 1.

const ORDER = Dict('X' => ["NW", "SW", "NE", "SE"], 'Z' => ["NW", "NE", "SW", "SE"])
const OFFS = OrderedDict("NW" => (0, 0), "NE" => (0, 1), "SW" => (1, 0), "SE" => (1, 1))

"""
    rect_checks(rows, cols, index)

Checks of a rectangular rotated patch; `rows`, `cols` are inclusive ranges and `index` maps
(i, j) to a qubit id.
"""
function rect_checks(rows, cols, index)
    r0, r1 = rows
    c0, c1 = cols
    out = Check[]
    for i in r0-1:r1, j in c0-1:c1
        t = iseven(i + j) ? 'X' : 'Z'
        cells = OrderedDict{String,Tuple{Int,Int}}()
        for (dn, (di, dj)) in OFFS
            a, b = i + di, j + dj
            (r0 <= a <= r1 && c0 <= b <= c1) && (cells[dn] = (a, b))
        end
        if length(cells) == 4
        elseif length(cells) == 2
            if i == r0 - 1 || i == r1
                t == 'Z' || continue
            elseif j == c0 - 1 || j == c1
                t == 'X' || continue
            else
                continue
            end
        else
            continue
        end
        pauli = PauliDict(index[cells[dn]] => t for dn in keys(cells))
        layer_of = OrderedDict(index[cells[dn]] => findfirst(==(dn), ORDER[t]) - 1 for dn in keys(cells))
        push!(out, Check(pauli; layer_of, group=t == 'X' ? 0 : 1, name="$t@($i,$j)"))
    end
    return out
end

"""
    Layout()

Assigns qubit ids to grid positions in the order rectangles are added.
"""
struct Layout
    index::OrderedDict{Tuple{Int,Int},Int}
end
Layout() = Layout(OrderedDict{Tuple{Int,Int},Int}())

function add_rect!(L::Layout, rows, cols)
    for i in rows[1]:rows[2], j in cols[1]:cols[2]
        haskey(L.index, (i, j)) || (L.index[(i, j)] = length(L.index))
    end
    return [L.index[(i, j)] for i in rows[1]:rows[2] for j in cols[1]:cols[2]]
end

coords(L::Layout) = OrderedDict(q => Float64[j, i] for ((i, j), q) in L.index)

"""
    patch_logicals(rows, cols, index)

`(Xbar, Zbar)` of a patch: X on its first row, Z on its first column.
"""
function patch_logicals(rows, cols, index)
    X = PauliDict(index[(rows[1], j)] => 'X' for j in cols[1]:cols[2])
    Z = PauliDict(index[(i, cols[1])] => 'Z' for i in rows[1]:rows[2])
    return X, Z
end

"""
    logical_from_label(label, basis, n)

Phase-free physical Pauli of a logical label such as "XIZ" or "Y_".
"""
function logical_from_label(label::AbstractString, basis, n::Integer)
    length(label) == length(basis) || throw(ArgumentError("label $label has wrong length for $(length(basis)) logical qubits"))
    parts = PauliDict[]
    for (ch, (X, Z)) in zip(label, basis)
        ch in ('I', '_') && continue
        ch in ('X', 'Y') && push!(parts, X)
        ch in ('Z', 'Y') && push!(parts, Z)
    end
    return pauli_product(n, parts)
end

parse_state(state, basis, n) = StateGens([(logical_from_label(lab, basis, n), Int(s)) for (lab, s) in state])

"""
    single_patch_readout(d, basis, state; rounds=0)

Destructive readout of one patch in `basis` ('X' or 'Z') after `rounds` noisy memory rounds.
"""
function single_patch_readout(d::Integer, basis::Char, state, rounds::Integer=0)
    L = Layout()
    rows, cols = (0, d - 1), (0, d - 1)
    add_rect!(L, rows, cols)
    n = length(L.index)
    checks = rect_checks(rows, cols, L.index)
    lb = LogicalBasis([patch_logicals(rows, cols, L.index)])
    epochs = rounds > 0 ? [Epoch(checks, rounds; name="memory")] : Epoch[]
    return Experiment(; n_data=n, active0=0:n-1, prep_checks=[c.pauli for c in checks],
                      logical_state=parse_state(state, lb, n), epochs, readout=OrderedDict(q => basis for q in 0:n-1),
                      targets=[logical_from_label(string(basis), lb, n)], coords=coords(L),
                      name="single_$(basis)_d$(d)_r$(rounds)", logical_basis=lb)
end

"""
    multi_patch_surgery(d, n_patches, kind, state, rounds; termination="conjugate")

Lattice surgery of a row of patches measuring every neighbouring product of type `kind` ("ZZ" or "XX").
"""
function multi_patch_surgery(d::Integer, n_patches::Integer, kind::AbstractString, state, rounds::Integer,
                             termination::AbstractString="conjugate")
    kind in ("ZZ", "XX") || throw(ArgumentError("kind must be ZZ or XX"))
    L = Layout()
    patches = Tuple{Tuple{Int,Int},Tuple{Int,Int}}[]
    seams = Vector{Int}[]
    for a in 0:n_patches-1
        off = a * (d + 1)
        rows, cols = kind == "ZZ" ? ((0, d - 1), (off, off + d - 1)) : ((off, off + d - 1), (0, d - 1))
        add_rect!(L, rows, cols)
        push!(patches, (rows, cols))
    end
    for a in 0:n_patches-2
        off = a * (d + 1) + d
        push!(seams, kind == "ZZ" ? add_rect!(L, (0, d - 1), (off, off)) : add_rect!(L, (off, off), (0, d - 1)))
    end
    n = length(L.index)
    seam_qubits = reduce(vcat, seams; init=Int[])
    sq = Set(seam_qubits)
    active0 = [q for q in 0:n-1 if !(q in sq)]
    prep = PauliDict[]
    for (rows, cols) in patches
        append!(prep, [c.pauli for c in rect_checks(rows, cols, L.index)])
    end
    span = n_patches * (d + 1) - 2
    merged = kind == "ZZ" ? rect_checks((0, d - 1), (0, span), L.index) : rect_checks((0, span), (0, d - 1), L.index)
    lb = LogicalBasis([patch_logicals(r, c, L.index) for (r, c) in patches])
    P = kind[1]
    targets = PauliDict[]
    for a in 1:n_patches-1
        lab = fill('I', n_patches)
        lab[a] = lab[a+1] = P
        push!(targets, logical_from_label(String(lab), lb, n))
    end
    seam_init = kind == "ZZ" ? 'X' : 'Z'
    conj = P == 'Z' ? 'X' : 'Z'
    rb = termination == "conjugate" ? conj : P
    return Experiment(; n_data=n, active0, prep_checks=prep, logical_state=parse_state(state, lb, n),
                      epochs=[Epoch(merged, rounds; new_qubits=OrderedDict(q => seam_init for q in seam_qubits), name="merged")],
                      readout=OrderedDict(q => rb for q in 0:n-1), targets, coords=coords(L),
                      name="surgery_$(kind)_$(n_patches)p_d$(d)_r$(rounds)_$(termination)", logical_basis=lb)
end

"""
    multi_patch_readout(d, bases, state; rounds=0)

Separate patches, `rounds` memory rounds, then destructive readout of patch i in `bases[i]`.
"""
function multi_patch_readout(d::Integer, bases::AbstractString, state, rounds::Integer=0)
    n_p = length(bases)
    L = Layout()
    patches = Tuple{Tuple{Int,Int},Tuple{Int,Int}}[]
    for a in 0:n_p-1
        off = a * (d + 1)
        rows, cols = (0, d - 1), (off, off + d - 1)
        add_rect!(L, rows, cols)
        push!(patches, (rows, cols))
    end
    n = length(L.index)
    checks = Check[]
    readout = OrderedDict{Int,Char}()
    for ((rows, cols), b) in zip(patches, bases)
        append!(checks, rect_checks(rows, cols, L.index))
        for i in rows[1]:rows[2], j in cols[1]:cols[2]
            readout[L.index[(i, j)]] = b
        end
    end
    lb = LogicalBasis([patch_logicals(r, c, L.index) for (r, c) in patches])
    targets = [logical_from_label(join(j == a ? b : 'I' for j in 1:n_p), lb, n) for (a, b) in enumerate(bases)]
    epochs = rounds > 0 ? [Epoch(checks, rounds; name="memory")] : Epoch[]
    return Experiment(; n_data=n, active0=0:n-1, prep_checks=[c.pauli for c in checks],
                      logical_state=parse_state(state, lb, n), epochs, readout, targets, coords=coords(L),
                      name="readout_$(bases)_d$(d)_r$(rounds)", logical_basis=lb)
end

"""
    setting_experiment(setting, d, n_patches, state, rounds)

Experiment and measured logical labels for "local:<bases>", "surgery:ZZ" or "surgery:XX".
"""
function setting_experiment(setting::AbstractString, d::Integer, n_patches::Integer, state, rounds::Integer)
    kind, arg = split(setting, ':')
    if kind == "local"
        length(arg) == n_patches || throw(ArgumentError("setting $setting does not have $n_patches letters"))
        exp = multi_patch_readout(d, arg, state, 0)
        labels = [join(i == a ? b : 'I' for i in 1:n_patches) for (a, b) in enumerate(arg)]
    elseif kind == "surgery"
        exp = multi_patch_surgery(d, n_patches, arg, state, rounds)
        P = arg[1]
        labels = [join(i in (a, a + 1) ? P : 'I' for i in 1:n_patches) for a in 1:n_patches-1]
    else
        throw(ArgumentError(setting))
    end
    return exp, labels
end
