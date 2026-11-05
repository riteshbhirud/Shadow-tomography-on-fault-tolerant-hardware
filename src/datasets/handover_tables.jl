# Tables that join earlier results (calibration, local calibration, native ensembles), and the basis
# exports. The tables are pandas-style CSV, so pandas' read and write typing is reproduced here.

#### pandas read_csv / to_csv typing

# pandas reads a column as int if every cell is an integer, as float if every non-empty cell is a
# number (empty = NaN), and as strings otherwise (empty = NaN). `str_cols` are always strings.
const _PD_INT = r"^[+-]?\d+$"
const _PD_FLOAT = r"^[+-]?(\d+\.?\d*|\d*\.\d+)([eE][+-]?\d+)?$|^[+-]?(nan|NaN|inf|Inf)$"

"""
    pandas_read(path; str_cols=()) -> Vector{Dict{String,Any}}

Rows of a CSV file with the cell types `pandas.read_csv` gives them: `Int`, `Float64` (NaN for an
empty numeric cell), `String`, or `missing` for an empty cell of a string column.
"""
function pandas_read(path::AbstractString; str_cols = ())
    lines = split(chomp(replace(read(path, String), "\r\n" => "\n")), '\n')
    header = _csv_split(lines[1])
    cells = [_csv_split(l) for l in lines[2:end]]
    cols = Dict{String,Vector{Any}}()
    for (j, h) in enumerate(header)
        raw = [r[j] for r in cells]
        nonempty = filter(!isempty, raw)
        cols[h] = if h in str_cols
            Any[isempty(t) ? missing : String(t) for t in raw]
        elseif !isempty(nonempty) && length(nonempty) == length(raw) && all(t -> occursin(_PD_INT, t), raw)
            Any[parse(Int, t) for t in raw]
        elseif all(t -> occursin(_PD_INT, t) || occursin(_PD_FLOAT, t), nonempty)
            Any[isempty(t) ? NaN : parse(Float64, t) for t in raw]
        else
            Any[isempty(t) ? missing : String(t) for t in raw]
        end
    end
    return [Dict{String,Any}(h => cols[h][i] for h in header) for i in eachindex(cells)]
end

# one cell as to_csv writes it: floats by repr, NaN and missing empty
_pd_text(x) = x
_pd_text(x::AbstractFloat) = isnan(x) ? missing : isinf(x) ? (x > 0 ? "inf" : "-inf") : Float64(x)
_pd_text(x::Bool) = x ? "True" : "False"

"""
    pandas_table(rows, cols; newline="\\n") -> CSVTable

Table in pandas `to_csv` format: a key missing from a row is NaN, and a numeric column with any float or
NaN in it is written as floats (an int 11 becomes `11.0`).
"""
function pandas_table(rows, cols; newline = "\n")
    data = [Any[get(r, c, NaN) for r in rows] for c in cols]
    for v in data
        numeric = all(x -> x isa Real && !(x isa Bool), v)
        if numeric && any(x -> x isa AbstractFloat, v)
            v .= Float64.(v)
        end
    end
    return CSVTable(collect(cols), [Any[_pd_text(data[j][i]) for j in eachindex(cols)] for i in eachindex(rows)], newline)
end

#### calibration.csv

const CALIBRATION_COLUMNS = ["section", "code", "distance", "setting_type", "target", "kind", "termination", "r",
    "noise_model", "p", "input_state", "shots", "errors", "q", "ci_lo", "ci_hi", "q_if_independent", "extra_qubits",
    "circuit_layers", "source"]

"""
    calibration_rows(root; input=joinpath) -> Vector{Dict}

The rows of calibration.csv from the result tables under `root`: the dataset-setting flip rates, then the
surface q(r) curves, the state-independence runs and BB72 q(r).
Each table is read from `input(root, "results", ...)` (scripts pass [`upstream`](@ref)).
"""
function calibration_rows(root::AbstractString; input = joinpath)
    res(f) = input(root, "results", "taskA", f)
    rows = Dict{String,Any}[]
    for x in pandas_read(input(root, "results", "calibration", "dataset_settings_q.csv"))
        sid = string(x["setting_id"])
        push!(rows, Dict{String,Any}("section" => "dataset_settings", "code" => x["code"],
            "distance" => x["code"] == "surface" ? 3 : 6, "setting_type" => x["setting_id"], "target" => x["target"],
            "kind" => x["kind"], "termination" => startswith(sid, "surgery") || startswith(sid, "gauge") ? "conjugate" : "",
            "r" => x["r"], "noise_model" => x["noise"], "p" => x["p"], "input_state" => "+1 eigenstate of every string",
            "shots" => x["shots"], "errors" => x["flips"], "q" => x["q"], "ci_lo" => x["ci_lo"], "ci_hi" => x["ci_hi"],
            "q_if_independent" => x["q_if_independent"], "source" => "results/calibration/dataset_settings_q.csv"))
    end
    for f in ("surface_qr.csv", "surface_state_independence.csv")
        for x in pandas_read(res(f); str_cols = ("state",))
            push!(rows, Dict{String,Any}("section" => "scaling", "code" => "surface", "distance" => x["d"],
                "setting_type" => "surgery:$(x["kind"])", "target" => x["kind"], "kind" => "string",
                "termination" => x["term"], "r" => x["r"], "noise_model" => x["noise"], "p" => x["p"],
                "input_state" => x["state"], "shots" => x["shots"], "errors" => x["errors"], "q" => x["q"],
                "ci_lo" => x["ci_lo"], "ci_hi" => x["ci_hi"], "extra_qubits" => x["cost_extra_qubits"],
                "circuit_layers" => x["cost_circuit_layers"], "source" => "results/taskA/$f"))
        end
    end
    for x in pandas_read(res("bb72_qr.csv"); str_cols = ("state",))
        push!(rows, Dict{String,Any}("section" => "scaling", "code" => "bb72", "distance" => 6,
            "setting_type" => "gauge:$(x["P"])0", "target" => "$(x["P"])0", "kind" => "string",
            "termination" => "conjugate", "r" => x["r"], "noise_model" => "bravyi_cal", "p" => x["p"],
            "input_state" => x["state"], "shots" => x["shots"], "errors" => x["errors"], "q" => x["q"],
            "ci_lo" => x["ci_lo"], "ci_hi" => x["ci_hi"], "extra_qubits" => x["cost_extra_qubits"],
            "circuit_layers" => x["cost_circuit_layers"], "source" => "results/taskA/bb72_qr.csv"))
    end
    return rows
end

#### calibration_local.csv

const CALIBRATION_LOCAL_COLUMNS = ["code", "distance", "setting_id", "r", "target", "kind", "shots", "flips", "q",
    "ci_lo", "ci_hi", "noise", "p"]

"""
    calibration_local_join(surface_q, bb72_q) -> Vector{Dict}

The rows of calibration_local.csv from `surface_local_q.csv` and, if it exists, `bb72_local_q.csv` (per
setting and target, at the rounds parameter).
"""
function calibration_local_join(surface_q::AbstractString, bb72_q::AbstractString)
    rows = Dict{String,Any}[]
    for x in pandas_read(surface_q)
        push!(rows, Dict{String,Any}("code" => "surface", "distance" => x["d"], "setting_id" => "patch:$(x["basis"])",
            "r" => x["memory_rounds"], "target" => x["basis"], "kind" => "single_patch", "shots" => x["shots"],
            "flips" => x["errors"], "q" => x["q"], "ci_lo" => x["ci_lo"], "ci_hi" => x["ci_hi"], "noise" => x["noise"],
            "p" => x["p"]))
    end
    if isfile(bb72_q)
        for x in pandas_read(bb72_q)
            push!(rows, Dict{String,Any}("code" => "bb72", "distance" => 6, "setting_id" => x["setting_id"],
                "r" => x["rounds_parameter"], "target" => x["target"], "kind" => x["kind"], "shots" => x["shots"],
                "flips" => x["flips"], "q" => x["q"], "ci_lo" => x["ci_lo"], "ci_hi" => x["ci_hi"],
                "noise" => x["noise"], "p" => x["p"]))
        end
    end
    return rows
end

#### native_ensembles/<code>_lambda_w3.csv

const NATIVE_ENSEMBLES = ["readout_only", "shifts", "ibm", "aut", "albert_two_fold", "aut_plus_all_folds"]
# the gross second fold family is unverified, so it is left out
const NATIVE_EXCLUDE = Set([("gross", "aut_plus_all_folds")])

"""
    read_npy(path) -> Vector{Float32}

A 1-D little-endian float32 `.npy` file (the λ_mix tables of Task B).
"""
function read_npy(path::AbstractString)
    open(path) do io
        read(io, 6) == UInt8[0x93, 'N', 'U', 'M', 'P', 'Y'] || error("not an npy file: $path")
        major = read(io, UInt8)
        read(io, UInt8)
        hlen = major == 1 ? Int(ltoh(read(io, UInt16))) : Int(ltoh(read(io, UInt32)))
        header = String(read(io, hlen))
        occursin("'descr': '<f4'", header) || error("expected float32: $header")
        occursin("'fortran_order': False", header) || error("expected C order: $header")
        n = parse(Int, match(r"'shape': \((\d+),\)", header)[1])
        return ltoh.(read!(io, Vector{Float32}(undef, n)))
    end
end

"""
    limit_denominator(x::Rational, maxden) -> Rational

Closest fraction to `x` with denominator at most `maxden`, by continued fractions. Ties go to the
convergent p1/q1, as in Python's `Fraction.limit_denominator`.
"""
function limit_denominator(x::Rational{BigInt}, maxden::Integer)
    denominator(x) <= maxden && return x
    p0, q0, p1, q1 = big(0), big(1), big(1), big(0)
    n, d = numerator(x), denominator(x)
    while true
        a = fld(n, d)
        q2 = q0 + a * q1
        q2 > maxden && break
        p0, q0, p1, q1 = p1, q1, p0 + a * p1, q2
        n, d = d, n - a * d
    end
    k = fld(maxden - q0, q1)
    return 2 * d * (q0 + k * q1) <= denominator(x) ? p1 // q1 : (p0 + k * p1) // (q0 + k * q1)
end

"""A fraction as "n/d", or "n" for an integer."""
pyfraction(x::Rational) = denominator(x) == 1 ? string(numerator(x)) : "$(numerator(x))/$(denominator(x))"

"""Index of a k-qubit Pauli label in the λ tables: bit i = X on qubit i, bit k+i = Z on qubit i."""
function pauli_index(label::AbstractString, k::Integer)
    v = 0
    for (i, ch) in enumerate(label)
        ch in ('X', 'Y') && (v |= 1 << (i - 1))
        ch in ('Z', 'Y') && (v |= 1 << (k + i - 1))
        ch in "IXYZ" || throw(ArgumentError("bad Pauli letter $ch"))
    end
    return v
end

"""All k-qubit Pauli labels of weight 1..wmax, ordered by weight, then qubit subset, then letters."""
function low_weight_labels(k::Integer, wmax::Integer)
    labels = String[]
    for w in 1:wmax, qs in _combinations_sorted(k, w), letters in product_strings("XYZ", w)
        s = fill('I', k)
        for (q, l) in zip(qs, letters)
            s[q] = l
        end
        push!(labels, String(s))
    end
    return labels
end

# size-w subsets of 1:n in lexicographic order
function _combinations_sorted(n, w)
    out = Vector{Int}[]
    function rec(start, acc)
        length(acc) == w && return push!(out, copy(acc))
        for i in start:n
            push!(acc, i)
            rec(i + 1, acc)
            pop!(acc)
        end
    end
    rec(1, Int[])
    return out
end

"""
    native_ensemble_table(code, tables; k=12) -> CSVTable

Rows of `<code>_lambda_w3.csv`: every Pauli of weight ≤ 3 with its exact λ_P (the float32 value
turned back into a fraction with denominator ≤ 100000) and 1/λ_P per ensemble. `tables` maps
ensemble name => λ_mix vector.
"""
function native_ensemble_table(code::AbstractString, tables::AbstractDict; k::Integer = 12)
    ens = [e for e in NATIVE_ENSEMBLES if !((code, e) in NATIVE_EXCLUDE)]
    header = vcat(["pauli", "weight"], ["lambda_$e" for e in ens], ["norm2_$e" for e in ens])
    rows = Vector{Any}[]
    for lab in low_weight_labels(k, 3)
        v = pauli_index(lab, k)
        ls = [limit_denominator(Rational{BigInt}(Float64(tables[e][v+1])), 100000) for e in ens]
        push!(rows, vcat(Any[lab, count(!=('I'), lab)], pyfraction.(ls), [x > 0 ? pyfraction(1 / x) : "inf" for x in ls]))
    end
    return CSVTable(header, rows, "\r\n")
end

#### bases/

const BASIS_CONVENTION = "physical qubit q = block*l*m + a*m + b for monomial x^a y^b (block 0 = L, 1 = R); " *
                         "H_X = [A|B], H_Z = [B^T|A^T]; logical qubit i is 0-based"

"""Basis file content: sorted 0-based supports of every Xbar_i and Zbar_i."""
basis_record(code_desc, l, m, xs, zs, note) = OrderedDict{String,Any}("code" => code_desc, "l" => l, "m" => m,
    "convention" => BASIS_CONVENTION, "note" => note, "xbar_support" => [sort!(Int.(collect(x))) for x in xs],
    "zbar_support" => [sort!(Int.(collect(z))) for z in zs])

"""
    handover_bases() -> Vector{Pair{String,OrderedDict}}

The four basis files (`file name => content`): the BB72 oracle and local bases, and the Task B f,g,h
bases of BB72 and gross.
"""
function handover_bases()
    bb = "BB72 [[72,12,6]], A = x^3 + y + y^2, B = y^3 + x + x^2"
    out = Pair{String,OrderedDict{String,Any}}[]
    for (fname, basis, note) in (
        ("bb72_oracle_basis.json", Circuits.bb_setup().basis,
         "basis of every BB72 shot record (INTERFACE.md section 4); qubit 0 = the weight-6 pair measured by gauging"),
        ("bb72_local_basis.json", Circuits.local_basis(),
         "basis of the BB72 local-setting datasets (Bell / GHZ on logical qubits 0, 1, 2; D-37): qubits 0-2 have " *
         "weight-6 Xbar/Zbar with overlap 1 and pairwise disjoint supports; qubit 0 equals the oracle basis' qubit 0"))
        @assert all(Set(values(x)) == Set(['X']) && Set(values(z)) == Set(['Z']) for (x, z) in basis)
        push!(out, fname => basis_record(bb, 6, 6, [keys(x) for (x, _) in basis], [keys(z) for (_, z) in basis], note))
    end
    for (code, fname) in (("BB72", "bb72_taskB_basis.json"), ("gross", "gross_taskB_basis.json"))
        c = TaskB.make_bb_code(code, "fgh")
        n = c.n
        @assert !any(!iszero, c.lx[:, n+1:end]) && !any(!iszero, c.lz[:, 1:n])
        xs = [findall(!iszero, c.lx[i, 1:n]) .- 1 for i in axes(c.lx, 1)]
        zs = [findall(!iszero, c.lz[i, n+1:end]) .- 1 for i in axes(c.lz, 1)]
        push!(out, fname => basis_record("$code: $(c.meta[:basis])", c.meta[:l], c.meta[:m], xs, zs,
            "basis of handover/native_ensembles/" * code * "_* per-Pauli files; NOT the oracle basis"))
    end
    return out
end
