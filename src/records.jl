# Shot-record files (format v1.1): reading, the loader API, and writers that reproduce the stored files' bytes exactly.

using JSON

#### Python-compatible JSON

"""
    pyfloat(x)

`x` formatted like Python's `repr(float)`: shortest round-trip digits, scientific notation
below 1e-4 and from 1e16 on, with a two-digit exponent (`1e-05`, not `1e-5`).
"""
function pyfloat(x::Float64)
    isnan(x) && return "NaN"
    isinf(x) && return x > 0 ? "Infinity" : "-Infinity"
    x == 0 && return signbit(x) ? "-0.0" : "0.0"
    s = string(abs(x))  # Julia also prints the shortest round-trip digits
    m, e = occursin('e', s) ? (split(s, 'e')[1], parse(Int, split(s, 'e')[2])) : (s, 0)
    ip, fp = occursin('.', m) ? split(m, '.') : (m, "")
    digits = ip * fp
    exp10 = length(ip) + e  # value = 0.digits × 10^exp10
    nz = findfirst(!=('0'), digits)
    digits = rstrip(digits[nz:end], '0')
    exp10 -= nz - 1
    isempty(digits) && (digits = "0")
    dexp = exp10 - 1
    out = if -4 <= dexp < 16
        if exp10 <= 0
            "0." * "0"^(-exp10) * digits
        elseif exp10 >= length(digits)
            digits * "0"^(exp10 - length(digits)) * ".0"
        else
            digits[1:exp10] * "." * digits[(exp10 + 1):end]
        end
    else
        mant = length(digits) == 1 ? digits : digits[1:1] * "." * digits[2:end]
        mant * "e" * (dexp < 0 ? "-" : "+") * lpad(string(abs(dexp)), 2, '0')
    end
    return signbit(x) ? "-" * out : out
end

# Python's json escapes everything outside printable ASCII (ensure_ascii=True).
function _pyjson_string(io::IO, s::AbstractString)
    print(io, '"')
    for c in s
        u = UInt32(c)
        if c == '"'
            print(io, "\\\"")
        elseif c == '\\'
            print(io, "\\\\")
        elseif c == '\n'
            print(io, "\\n")
        elseif c == '\r'
            print(io, "\\r")
        elseif c == '\t'
            print(io, "\\t")
        elseif c == '\b'
            print(io, "\\b")
        elseif c == '\f'
            print(io, "\\f")
        elseif u < 0x20 || u > 0x7f
            if u > 0xffff  # surrogate pair, as Python writes it
                u -= 0x10000
                print(io, "\\u", string(0xd800 + (u >> 10); base = 16, pad = 4))
                print(io, "\\u", string(0xdc00 + (u & 0x3ff); base = 16, pad = 4))
            else
                print(io, "\\u", string(u; base = 16, pad = 4))
            end
        else
            print(io, c)
        end
    end
    print(io, '"')
end

"""
    pyjson(x; indent=nothing, separators=nothing)

JSON text of `x` exactly as Python's `json.dumps(x, indent=..., separators=...)` writes it.
Dictionaries keep their insertion order, which is the order `JSON.parse` reads.
"""
function pyjson(x; indent::Union{Nothing,Int} = nothing, separators = nothing)
    item, key = separators === nothing ? (indent === nothing ? ", " : ",", ": ") : separators
    io = IOBuffer()
    _pyjson(io, x, indent, item, key, 0)
    return String(take!(io))
end

_pyjson(io, ::Nothing, _...) = print(io, "null")
_pyjson(io, x::Bool, _...) = print(io, x ? "true" : "false")
_pyjson(io, x::Integer, _...) = print(io, x)
_pyjson(io, x::AbstractFloat, _...) = print(io, pyfloat(Float64(x)))
_pyjson(io, x::AbstractString, _...) = _pyjson_string(io, x)
_pyjson(io, x::Symbol, _...) = _pyjson_string(io, String(x))

function _pyjson_items(f, io, items, open, close, indent, item, level)
    isempty(items) && return print(io, open, close)
    print(io, open)
    for (k, it) in enumerate(items)
        k > 1 && print(io, item)
        indent === nothing || print(io, "\n", " "^(indent * (level + 1)))
        f(it)
    end
    indent === nothing || print(io, "\n", " "^(indent * level))
    print(io, close)
end

_pyjson(io, x::Union{AbstractVector,Tuple}, indent, item, key, level) =
    _pyjson_items(v -> _pyjson(io, v, indent, item, key, level + 1), io, x, '[', ']', indent, item, level)

_pyjson(io, x::AbstractDict, indent, item, key, level) = _pyjson_items(io, collect(pairs(x)), '{', '}', indent, item, level) do (k, v)
    _pyjson_string(io, string(k))
    print(io, key)
    _pyjson(io, v, indent, item, key, level + 1)
end

#### Raw files, kept with their types

"""
    read_meta(path) -> JSON.Object

The dataset's `meta.json`, with key order and number types as in the file.
"""
read_meta(path::AbstractString) = JSON.parsefile(path)

"""
    write_meta(path, meta)

Write `meta` like `json.dump(meta, f, indent=1)`: one-space indent, no trailing newline.
"""
write_meta(path::AbstractString, meta) = write(path, pyjson(meta; indent = 1))

"""
    read_records_raw(path) -> Vector

One parsed JSON object per line of a `records.jsonl` file.
"""
read_records_raw(path::AbstractString) = [JSON.parse(l) for l in eachline(path) if !isempty(strip(l))]

"""
    write_records_raw(path, records)

Write one compact JSON object per line, like `json.dumps(r, separators=(",", ":"))`.
"""
function write_records_raw(path::AbstractString, records)
    open(path, "w") do io
        for r in records
            print(io, pyjson(r; separators = (",", ":")), "\n")
        end
    end
end

"""
    CSVTable(header, rows, newline)

A CSV file as stored. Each cell is `missing` (empty), an `Int`, a
`Float64` or a `String`, decided from its text so that writing it back gives the same bytes.
`newline` is `"\\r\\n"` for files from Python's `csv` module (`exact.csv`) and `"\\n"` for
pandas (`expected.csv`).
"""
struct CSVTable
    header::Vector{String}
    rows::Vector{Vector{Any}}
    newline::String
end

function _csv_cell(tok::AbstractString)
    isempty(tok) && return missing
    occursin(r"^-?\d+$", tok) && return parse(Int, tok)
    if occursin(r"^-?(\d+\.\d*|\d*\.\d+|\d+)(e[+-]?\d+)?$", tok) || tok in ("nan", "inf", "-inf")
        x = parse(Float64, tok)
        pyfloat(x) == tok && return x  # otherwise keep the exact text
    end
    return String(tok)
end

# Python's csv module and pandas quote a field only if it contains a separator, quote or newline.
function _csv_split(line::AbstractString)
    fields = String[]
    buf = IOBuffer()
    inquote = false
    i = firstindex(line)
    while i <= lastindex(line)
        c = line[i]
        if inquote
            if c == '"'
                nxt = nextind(line, i)
                if nxt <= lastindex(line) && line[nxt] == '"'
                    print(buf, '"')
                    i = nxt
                else
                    inquote = false
                end
            else
                print(buf, c)
            end
        elseif c == '"'
            inquote = true
        elseif c == ','
            push!(fields, String(take!(buf)))
        else
            print(buf, c)
        end
        i = nextind(line, i)
    end
    push!(fields, String(take!(buf)))
    return fields
end

_csv_text(::Missing) = ""
_csv_text(x::Integer) = string(x)
_csv_text(x::AbstractFloat) = pyfloat(Float64(x))
function _csv_text(x::AbstractString)
    any(c -> c in (',', '"', '\n', '\r'), x) || return String(x)
    return "\"" * replace(x, "\"" => "\"\"") * "\""
end

"""
    read_csv_table(path) -> CSVTable
"""
function read_csv_table(path::AbstractString)
    text = read(path, String)
    newline = occursin("\r\n", text) ? "\r\n" : "\n"
    lines = split(chomp(replace(text, "\r\n" => "\n")), '\n')
    header = _csv_split(lines[1])
    rows = [Any[_csv_cell(t) for t in _csv_split(l)] for l in lines[2:end]]
    return CSVTable(header, rows, newline)
end

"""
    write_csv_table(path, table)
"""
function write_csv_table(path::AbstractString, t::CSVTable)
    open(path, "w") do io
        print(io, join(_csv_text.(t.header), ","), t.newline)
        for r in t.rows
            print(io, join(_csv_text.(r), ","), t.newline)
        end
    end
end

#### Loader API (unchanged from handover/julia/load_records.jl)

"""
    Record

One shot: `setting_id`, the measured `setting` strings, decoded `outcomes` (bit 0 = eigenvalue +1),
`confidence` (`Inf`/`NaN` stand for `"inf"`/`null` in the file) with its `confidence_kind`, and the
setting's `cost`.
"""
struct Record
    shot::Int
    setting_id::String
    setting::Vector{String}
    outcomes::Vector{Int}
    confidence::Float64
    confidence_kind::String
    cost::AbstractDict
end

"""
    ExactRow

A row of `exact.csv`: a learnable logical Pauli, its exact expectation, `lambda` = λ_P, and the
settings that measure it.
"""
struct ExactRow
    pauli::String
    exact::Float64
    lambda::Float64
    settings::Vector{String}
end

_confidence(x) = x === nothing ? NaN : (x == "inf" ? Inf : Float64(x))

"""
    load_records(path) -> Vector{Record}

Read a `records.jsonl` file. Errors if the shots are not in draw order 0, 1, 2, ….
"""
function load_records(path::AbstractString)
    recs = Record[]
    for d in read_records_raw(path)
        push!(recs, Record(d["shot"], d["setting_id"], String.(d["setting"]), Int.(d["outcomes"]),
            _confidence(d["confidence"]), String(get(d, "confidence_kind", "unspecified")), d["cost"]))
    end
    all(r.shot == i - 1 for (i, r) in enumerate(recs)) || error("$path: shots not in draw order 0..N-1")
    return recs
end

"""
    load_exact(path) -> Vector{ExactRow}
"""
function load_exact(path::AbstractString)
    t = read_csv_table(path)
    col = Dict(h => i for (i, h) in enumerate(t.header))
    return [ExactRow(r[col["pauli"]], Float64(r[col["exact"]]), Float64(r[col["lambda_P"]]),
        String.(split(r[col["settings"]], ";"))) for r in t.rows]
end

"""
    load_dataset(prefix) -> (; meta, records, exact)

Load `<prefix>.meta.json`, `<prefix>.records.jsonl` and `<prefix>.exact.csv`.
"""
load_dataset(prefix::AbstractString) = (meta = read_meta(prefix * ".meta.json"),
    records = load_records(prefix * ".records.jsonl"), exact = load_exact(prefix * ".exact.csv"))

const _PAULI_BITS = Dict('I' => (0, 0), 'X' => (1, 0), 'Z' => (0, 1), 'Y' => (1, 1))
const _PAULI_CHAR = Dict(v => k for (k, v) in _PAULI_BITS)
_letterwise(a, b) = String([_PAULI_CHAR[(_PAULI_BITS[x][1] ⊻ _PAULI_BITS[y][1], _PAULI_BITS[x][2] ⊻ _PAULI_BITS[y][2])]
                            for (x, y) in zip(a, b)])

"""
    subset_index(setting) -> Dict(P => indices of the strings whose product is P)

Valid only when every qubit carries one letter type within the setting (true for every v1.1
setting): then the letter-wise product has sign +1. Errors otherwise, since e.g. XX·ZZ = −YY.
"""
function subset_index(setting::Vector{String})
    for q in 1:length(setting[1])
        letters = Set(s[q] for s in setting if s[q] != 'I')
        length(letters) <= 1 || error("setting $setting mixes letters on qubit $(q-1): products carry phases")
    end
    out = Dict{String,Vector{Int}}()
    m = length(setting)
    for mask in 1:(2^m - 1)
        idx = [j for j in 1:m if (mask >> (j - 1)) & 1 == 1]
        P = foldl(_letterwise, setting[idx])
        all(==('I'), P) || haskey(out, P) || (out[P] = idx)
    end
    return out
end

"""
    LocalShot

One shot of a local Pauli setting: `bases[i]` ∈ ('X','Y','Z') for logical qubit i, `outcomes[i]` = ±1,
the decoder `gap` and `gap_kind`, and the setting's `cost`.
"""
struct LocalShot
    bases::Vector{Char}
    outcomes::Vector{Int}
    gap::Float64
    gap_kind::String
    cost::AbstractDict
end

"""
    local_view(records) -> Vector{LocalShot}

The (bases, ±1 outcomes) view of local-setting records. Errors if a record is not a local setting.
"""
function local_view(recs::Vector{Record})
    out = LocalShot[]
    for r in recs
        n = length(r.setting[1])
        bases = fill('?', n)
        vals = zeros(Int, n)
        for (s, b) in zip(r.setting, r.outcomes)
            q = findall(!=('I'), s)
            length(q) == 1 || error("shot $(r.shot): string $s is not single-qubit; not a local setting")
            bases[q[1]] == '?' || error("shot $(r.shot): qubit $(q[1]-1) measured twice")
            bases[q[1]] = s[q[1]]
            vals[q[1]] = 1 - 2b
        end
        any(==('?'), bases) && error("shot $(r.shot): not every qubit measured")
        push!(out, LocalShot(bases, vals, r.confidence, r.confidence_kind, r.cost))
    end
    return out
end
