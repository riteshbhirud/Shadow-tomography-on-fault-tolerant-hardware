# Writing the Task B result files in their stored formats (scripts/taskb_run_all.jl): JSON as
# Python's json.dump, CSV as csv.DictWriter, the λ tables as .npy, the kept generators as .npz (numpy, bridged).

"""
    write_pyjson(path, x; indent=nothing)

Write `x` as Python's `json.dump(x, fh, indent=indent)` does (no trailing newline).
"""
write_pyjson(path::AbstractString, x; indent=nothing) = write(path, pyjson(x; indent))

# a cell as Python's csv module writes str(value)
function _dictcsv_cell(v)
    v === nothing && return ""
    v isa Bool && return v ? "True" : "False"
    if v isa AbstractFloat
        isnan(v) && return "nan"
        isinf(v) && return v > 0 ? "inf" : "-inf"
        return pyfloat(Float64(v))
    end
    s = string(v)
    any(c -> c in (',', '"', '\n', '\r'), s) && return "\"" * replace(s, "\"" => "\"\"") * "\""
    return s
end

"""
    write_dictcsv(path, rows)

Write dict rows as `csv.DictWriter(fh, fieldnames=<keys of all rows, in first-seen order>)` with `writeheader()`
and `writerows(rows)`: CRLF line ends, missing keys empty.
"""
function write_dictcsv(path::AbstractString, rows)
    cols = String[]
    for r in rows, k in keys(r)
        k in cols || push!(cols, k)
    end
    open(path, "w") do io
        print(io, join(_dictcsv_cell.(cols), ","), "\r\n")
        for r in rows
            print(io, join((_dictcsv_cell(get(r, c, "")) for c in cols), ","), "\r\n")
        end
    end
end

"""
    write_npy(path, v::Vector{Float32})

A 1-D little-endian float32 array in numpy's .npy format, version 1.0 (header padded to 64 bytes), as `np.save`.
"""
function write_npy(path::AbstractString, v::AbstractVector{Float32})
    header = "{'descr': '<f4', 'fortran_order': False, 'shape': ($(length(v)),), }"
    pad = 64 - (10 + length(header) + 1) % 64
    header *= " "^(pad % 64) * "\n"
    open(path, "w") do io
        write(io, UInt8[0x93], "NUMPY", UInt8[0x01, 0x00], htol(UInt16(length(header))), header)
        write(io, htol.(v))
    end
end

"""
    write_generators_npz(path, kept)

The kept generators of every ensemble as `np.savez_compressed` arrays of shape (count, 2k, 2k), uint8, in
[`ENSEMBLES`](@ref) order (numpy through the bridge; the zip timestamps differ from run to run anyway).
"""
function write_generators_npz(path::AbstractString, kept::AbstractDict, k::Integer)
    np = pymodule(:numpy)
    arrs = map(ENSEMBLES) do en
        a = zeros(UInt8, length(kept[en]), 2k, 2k)
        for (i, m) in enumerate(kept[en])
            a[i, :, :] .= m
        end
        Symbol(en) => np.ascontiguousarray(np.asarray(a))
    end
    np.savez_compressed(path; arrs...)
end

"""A gate list in the JSON layout of the stored files: [name, qubits...] (a SWAPPERM carries its list)."""
gate_lists(gates) = [g.name == :SWAPPERM ? Any["SWAPPERM", g.qubits] : Any[string(g.name); g.qubits] for g in gates]

_int_rows(m) = [Int.(m[i, :]) for i in axes(m, 1)]

"""
    pauli_corrections_record(elements)

The `pauli_corrections_<code>.json` list: the Pauli correction of every candidate element.
"""
pauli_corrections_record(alls) =
    [OrderedDict{String,Any}("label" => e.label, "pauli_correction" => gate_lists(e.pauli_correction)) for e in alls]

"""
    export_record(name, code, sample)

The `julia_export_<code>.json` record: code stabilizers and logical basis, and for every sampled element its
gate list, logical action, logical signs and (if corrected) the bare gate list.
"""
function export_record(name, code::StabilizerCode, sample)
    out = OrderedDict{String,Any}("code" => name, "n" => code.n, "k" => logical_k(code),
                                  "stabs" => _int_rows(code.stabs), "logicals" => _int_rows(logicals(code)),
                                  "elements" => Any[])
    for e in sample
        d = OrderedDict{String,Any}("label" => e.label, "gates" => gate_lists(e.gates),
                                    "logical" => _int_rows(e.logical), "logical_signs" => logical_signs(code, e),
                                    "n_pauli_correction_gates" => length(e.pauli_correction))
        isempty(e.pauli_correction) || (d["gates_uncorrected"] = gate_lists(e.gates_uncorrected))
        push!(out["elements"], d)
    end
    return out
end

"""
    crosscheck_report(name, code, sample; trials=3, rng=MersenneTwister(20261030)) -> (text, all_ok)

The state-level sign check of the sampled elements ([`signed_crosscheck`](@ref)), one line per element in the
format of the old `julia/taskb_crosscheck.jl` log.
"""
function crosscheck_report(name, code::StabilizerCode, sample; trials::Integer=3, rng=nothing)
    rng = rng === nothing ? MersenneTwister(20261030) : rng
    io = IOBuffer()
    allok = true
    for e in sample
        ok, neg, sneg = signed_crosscheck(code, e, logical_signs(code, e), rng; trials)
        println(io, rpad(first(e.label, 70), 72), " signed_prediction_ok=", ok, " negative_control_ok=", neg,
                " sign_negative_control_ok=", sneg, " n_pauli_corr=", length(e.pauli_correction))
        allok &= ok & neg & sneg
    end
    println(io, "ALL_OK=", allok, " code=", name, " elements=", length(sample))
    return String(take!(io)), allok
end
