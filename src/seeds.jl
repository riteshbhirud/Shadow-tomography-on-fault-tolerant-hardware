# Collision-free seeds from a description of the job.

using SHA

"""
    pyrepr(x) -> String

Python's `repr` of a string, integer, float, `Bool`, `nothing` (None), vector (list) or tuple.
"""
pyrepr(x::Bool) = x ? "True" : "False"
pyrepr(x::Integer) = string(x)
pyrepr(::Nothing) = "None"
pyrepr(x::Symbol) = pyrepr(String(x))
function pyrepr(x::AbstractFloat)
    isnan(x) && return "nan"
    isinf(x) && return x > 0 ? "inf" : "-inf"
    return pyfloat(Float64(x))
end
pyrepr(x::AbstractVector) = "[" * join(map(pyrepr, x), ", ") * "]"
pyrepr(x::Tuple) = length(x) == 1 ? "(" * pyrepr(x[1]) * ",)" : "(" * join(map(pyrepr, x), ", ") * ")"

# str.isprintable() is false for these Unicode categories (and true for the ASCII space).
const _NONPRINTABLE = (Base.Unicode.UTF8PROC_CATEGORY_CC, Base.Unicode.UTF8PROC_CATEGORY_CF,
    Base.Unicode.UTF8PROC_CATEGORY_CS, Base.Unicode.UTF8PROC_CATEGORY_CO, Base.Unicode.UTF8PROC_CATEGORY_CN,
    Base.Unicode.UTF8PROC_CATEGORY_ZL, Base.Unicode.UTF8PROC_CATEGORY_ZP, Base.Unicode.UTF8PROC_CATEGORY_ZS)

function pyrepr(s::AbstractString)
    q = occursin('\'', s) && !occursin('"', s) ? '"' : '\''
    io = IOBuffer()
    print(io, q)
    for c in s
        u = UInt32(c)
        if c == q || c == '\\'
            print(io, '\\', c)
        elseif c == '\t'
            print(io, "\\t")
        elseif c == '\n'
            print(io, "\\n")
        elseif c == '\r'
            print(io, "\\r")
        elseif c != ' ' && Base.Unicode.category_code(c) in _NONPRINTABLE
            if u <= 0xff
                print(io, "\\x", string(u; base = 16, pad = 2))
            elseif u <= 0xffff
                print(io, "\\u", string(u; base = 16, pad = 4))
            else
                print(io, "\\U", string(u; base = 16, pad = 8))
            end
        else
            print(io, c)
        end
    end
    print(io, q)
    return String(take!(io))
end

"""
    seed_for(parts...) -> Int

32-bit seed from the SHA-256 of the parts' Python reprs joined by `|` (Python reprs so the seeds
agree with those of the stored datasets).
`seed_for("surface_calib", d, kind, term, state, r, noise, chunk)` differs for every configuration and
chunk, so two configurations never share a noise realisation.
"""
function seed_for(parts...)
    h = sha256(join(map(pyrepr, parts), "|"))
    return Int(ltoh(reinterpret(UInt32, h[1:4])[1]))
end
