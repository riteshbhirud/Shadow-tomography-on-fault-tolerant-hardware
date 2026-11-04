# Circuit IR: stim instructions as plain data, emitted as stim text. `append!` fuses an
# instruction into the previous one when stim's Circuit.append would (same gate, same arguments,
# gate not an annotation), so the emitted text equals Python's `str(circuit)`.

@enum TargetKind::UInt8 QUBIT REC PAULI_X PAULI_Y PAULI_Z COMBINER

struct Target
    kind::TargetKind
    value::Int  # qubit id, or the negative record offset of a rec target
end

qubit(q::Integer) = Target(QUBIT, q)
rec(offset::Integer) = (offset < 0 || throw(ArgumentError("rec offsets are negative")); Target(REC, offset))
pauli_target(q::Integer, s::Char) = Target(s == 'X' ? PAULI_X : s == 'Y' ? PAULI_Y : PAULI_Z, q)

struct Instr
    name::String
    targets::Vector{Target}
    args::Vector{Float64}
end

"""
    Circuit

A stim circuit as a list of instructions; build it with `append!` and render with `stim_text`.
"""
struct Circuit
    ops::Vector{Instr}
end
Circuit() = Circuit(Instr[])

const UNFUSABLE = Set(["DETECTOR", "OBSERVABLE_INCLUDE", "TICK", "QUBIT_COORDS", "SHIFT_COORDS", "MPAD", "REPEAT"])

_targets(ts::AbstractVector{Target}) = ts
_targets(ts) = Target[t isa Target ? t : qubit(t) for t in ts]

"""
    append!(c, name, targets, args=Float64[])

Append a stim instruction; integer targets are qubits.
"""
function Base.append!(c::Circuit, name::AbstractString, targets=Target[], args=Float64[])
    ts = _targets(collect(targets))
    as = Float64[args...]
    if !isempty(c.ops)
        last = c.ops[end]
        if last.name == name && last.args == as && !(name in UNFUSABLE)
            append!(last.targets, ts)
            return c
        end
    end
    push!(c.ops, Instr(String(name), ts, as))
    return c
end

Base.append!(c::Circuit, i::Instr) = append!(c, i.name, i.targets, i.args)

function Base.append!(c::Circuit, other::Circuit)
    for i in other.ops
        append!(c, i.name, copy(i.targets), i.args)
    end
    return c
end

"""
    mpp_targets(p)

MPP targets of one Pauli product, factors in qubit order.
"""
function mpp_targets(p)
    out = Target[]
    for (k, (q, s)) in enumerate(sortedpauli(p))
        k > 1 && push!(out, Target(COMBINER, 0))
        push!(out, pauli_target(q, s))
    end
    return out
end

# C++ ostream default formatting of doubles, which stim uses for arguments.
fmtarg(x::Float64) = @sprintf("%g", x)

function _show_target(io::IO, t::Target)
    if t.kind == QUBIT
        print(io, t.value)
    elseif t.kind == REC
        print(io, "rec[", t.value, "]")
    else
        print(io, t.kind == PAULI_X ? 'X' : t.kind == PAULI_Y ? 'Y' : 'Z', t.value)
    end
end

"""
    stim_text(c)

The circuit in stim's text format, as `str(stim.Circuit)` prints it.
"""
function stim_text(c::Circuit)
    io = IOBuffer()
    for (k, op) in enumerate(c.ops)
        k > 1 && print(io, '\n')
        print(io, op.name)
        isempty(op.args) || print(io, '(', join(fmtarg.(op.args), ", "), ')')
        after_combiner = false
        for t in op.targets
            if t.kind == COMBINER
                print(io, '*')
                after_combiner = true
                continue
            end
            after_combiner || print(io, ' ')
            _show_target(io, t)
            after_combiner = false
        end
    end
    return String(take!(io))
end

Base.show(io::IO, ::MIME"text/plain", c::Circuit) = print(io, "Circuit with ", length(c.ops), " instructions")

function _parse_target(s::AbstractString)
    if startswith(s, "rec[")
        return rec(parse(Int, s[5:end-1]))
    elseif s[1] in ('X', 'Y', 'Z')
        return pauli_target(parse(Int, s[2:end]), s[1])
    else
        return qubit(parse(Int, s))
    end
end

"""
    parse_stim(text)

Parse flat stim text (no REPEAT blocks, no inverted targets) into a `Circuit`, fusing like stim.
"""
function parse_stim(text::AbstractString)
    c = Circuit()
    for line in eachline(IOBuffer(text))
        line = strip(line)
        isempty(line) && continue
        m = match(r"^([A-Z_0-9]+)(?:\(([^)]*)\))?\s*(.*)$", line)
        m === nothing && throw(ArgumentError("cannot parse stim line: $line"))
        name = m.captures[1]
        name == "REPEAT" && throw(ArgumentError("REPEAT blocks are not supported; flatten first"))
        args = m.captures[2] === nothing ? Float64[] : [parse(Float64, a) for a in split(m.captures[2], ',')]
        ts = Target[]
        for tok in split(m.captures[3])
            startswith(tok, '!') && throw(ArgumentError("inverted targets are not supported"))
            parts = split(tok, '*')
            for (k, p) in enumerate(parts)
                k > 1 && push!(ts, Target(COMBINER, 0))
                push!(ts, _parse_target(p))
            end
        end
        append!(c, name, ts, args)
    end
    return c
end

# stim.gate_data properties for the gates that appear in the circuits we build or parse.
const GATE_KIND = Dict(
    "M" => :measure, "MX" => :measure, "MY" => :measure, "MPP" => :measure,
    "R" => :reset, "RX" => :reset, "RY" => :reset,
    "CX" => :two, "CY" => :two, "CZ" => :two, "XCY" => :two, "XCZ" => :two, "YCX" => :two, "YCZ" => :two,
    "XCX" => :two, "YCY" => :two, "SWAP" => :two,
    "H" => :one, "S" => :one, "S_DAG" => :one, "SQRT_X" => :one, "SQRT_X_DAG" => :one, "SQRT_Y" => :one,
    "SQRT_Y_DAG" => :one, "X" => :one, "Y" => :one, "Z" => :one, "I" => :one, "H_YZ" => :one, "H_XY" => :one,
    "C_XYZ" => :one, "C_ZYX" => :one,
)

gate_kind(name::AbstractString) =
    get(() -> throw(ArgumentError("unsupported instruction $name")), GATE_KIND, name)

"""
    num_measurements(c)

Number of measurement records the circuit produces.
"""
function num_measurements(c::Circuit)
    n = 0
    for op in c.ops
        get(GATE_KIND, op.name, :none) == :measure || continue
        n += op.name == "MPP" ? count(t -> t.kind != COMBINER, op.targets) - count(t -> t.kind == COMBINER, op.targets) :
             length(op.targets)
    end
    return n
end

"""
    pycircuit(c)

The circuit as a Python `stim.Circuit` (through the bridge).
"""
pycircuit(c::Circuit) = pymodule(:stim).Circuit(stim_text(c))

"""
    normalized_text(text)

`str(stim.Circuit(text))`: stim's own rendering of a circuit text (through the bridge).
"""
normalized_text(text::AbstractString) = pyconvert(String, pybuiltins.str(pymodule(:stim).Circuit(text)))
