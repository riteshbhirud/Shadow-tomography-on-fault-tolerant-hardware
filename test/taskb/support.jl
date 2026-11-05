# Helpers for the Task B fixture tests.
@testmodule TaskBRef begin
    using FTShadows, JSON
    import FTShadows.TaskB as TB

    """
        same(julia, stored) -> Bool

    Exact structural equality of a Julia record and a parsed stored JSON value: equal keys, equal lengths,
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

    small_refs() = JSON.parsefile(joinpath(pkgdir(FTShadows), "test", "reference", "taskb", "small_references.json"))

    "A gate list in the JSON layout of the stored files: [name, qubits...] (a SWAPPERM carries its list)."
    gate_lists(gates) = [g.name == :SWAPPERM ? Any["SWAPPERM", g.qubits] : Any[string(g.name); g.qubits] for g in gates]

end
