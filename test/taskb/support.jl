# Shared helpers for the Task B parity tests: exact comparison against the saved reference results.
@testmodule TaskBRef begin
    using FTShadows, JSON
    import FTShadows.TaskB as TB
    const RESULTS = joinpath(pkgdir(FTShadows), "results", "taskB")

    """Whether the saved reference results are present (they are local only, not tracked)."""
    has_results() = isdir(RESULTS)

    function skip_message(what)
        @info "results/taskB not found: skipping the full parity check of $what (it holds the stored reference results)"
        return nothing
    end

    """Read a 1-D little-endian .npy file (the λ_mix tables)."""
    function read_npy(path)
        open(path) do io
            magic = read(io, 6)
            magic == UInt8[0x93, 'N', 'U', 'M', 'P', 'Y'] || error("not an npy file")
            major = read(io, UInt8); read(io, UInt8)
            hlen = major == 1 ? Int(read(io, UInt16)) : Int(read(io, UInt32))
            header = String(read(io, hlen))
            occursin("'descr': '<f4'", header) || error("expected float32: $header")
            occursin("'fortran_order': False", header) || error("expected C order")
            n = parse(Int, match(r"'shape': \((\d+),\)", header)[1])
            return ltoh.(read!(io, Vector{Float32}(undef, n)))
        end
    end

    """
        same(julia, python) -> Bool

    Exact structural equality of a Julia record and a parsed reference JSON value: equal keys, equal lengths,
    integers equal as integers, floats bitwise equal as floats (so `1` and `1.0` differ, as they do in the file).
    """
    same(a, b) = (r = first_difference(a, b, "")) === nothing ? true : (@warn "records differ" r; false)

    function first_difference(a, b, path)
        if b isa AbstractDict
            a isa AbstractDict || return "$path: not a dict"
            ka, kb = Set(string.(keys(a))), Set(string.(keys(b)))
            ka == kb || return "$path: keys differ: julia-only $(setdiff(ka, kb)), python-only $(setdiff(kb, ka))"
            for k in kb
                d = first_difference(a[k], b[k], "$path/$k")
                d === nothing || return d
            end
            return nothing
        elseif b isa AbstractVector
            (a isa AbstractVector || a isa Tuple) || return "$path: not a list"
            length(a) == length(b) || return "$path: length $(length(a)) vs $(length(b))"
            for i in eachindex(b)
                d = first_difference(a[i], b[i], "$path[$i]")
                d === nothing || return d
            end
            return nothing
        elseif b isa Bool
            return a isa Bool && a == b ? nothing : "$path: $a vs $b"
        elseif b isa Integer
            return a isa Integer && !(a isa Bool) && a == b ? nothing : "$path: $a vs $b"
        elseif b isa AbstractFloat
            return a isa AbstractFloat && (Float64(a) === Float64(b) || (isnan(a) && isnan(b))) ? nothing : "$path: $a vs $b"
        elseif b === nothing
            return a === nothing ? nothing : "$path: $a vs null"
        else
            return string(a) == string(b) ? nothing : "$path: $a vs $b"
        end
    end

    """
        csv_matches(rows, columns, path)

    Whether the CSV written by Python's csv.DictWriter equals the Julia rows cell by cell: floats are compared
    bitwise after parsing, integers and strings as text, booleans as Python's True/False.
    """
    function csv_matches(rows, path)
        lines = split(chomp(replace(read(path, String), "\r\n" => "\n")), "\n")
        header = split(lines[1], ",")
        length(lines) - 1 == length(rows) || (@warn "row count" length(lines) - 1 length(rows); return false)
        for (row, line) in zip(rows, lines[2:end])
            cells = split(line, ",")
            Set(keys(row)) == Set(header) || (@warn "columns differ" header keys(row); return false)
            for (col, cell) in zip(header, cells)
                v = row[col]
                ok = if v isa Bool
                    cell == (v ? "True" : "False")
                elseif v isa AbstractFloat
                    cell == "inf" ? isinf(v) : (occursin(r"[.eE]|inf|nan", cell) && parse(Float64, cell) === Float64(v))
                else
                    cell == string(v)
                end
                ok || (@warn "cell differs" col cell v line; return false)
            end
        end
        return true
    end

    small_refs() = JSON.parsefile(joinpath(TB.TASKB_REFDIR, "small_references.json"))

    "A gate list in the JSON layout of the reference files: [name, qubits...] (a SWAPPERM carries its list)."
    gate_lists(gates) = [g.name == :SWAPPERM ? Any["SWAPPERM", g.qubits] : Any[string(g.name); g.qubits] for g in gates]

    const CACHE = Dict{Any,Any}()
    cached(f, key) = get!(f, CACHE, key)

    """
        build(name)

    `TaskB.build_groups(name)` once per test process: the groups record, kept generators, candidate elements,
    folds and code. Needs Oscar to be loaded.
    """
    build(name) = cached(("build", name)) do
        record, kept, ens, folds, code = TB.build_groups(name)
        (; record, kept, ens, folds, code)
    end
end
