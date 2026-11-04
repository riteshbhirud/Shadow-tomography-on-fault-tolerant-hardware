# `pd.DataFrame(rows).to_csv(path, index=False)` for a list of row dicts, byte for byte.
# pandas picks one dtype per column, and that decides how a cell is printed: an integer column
# that also holds None or NaN becomes float ("3.0"), NaN and None print empty, inf prints "inf".

_is_int(x) = x isa Integer && !(x isa Bool)
_is_none(x) = x === nothing || x === missing

# pandas' column type: :bool, :int, :float or :object (strings, mixtures, all-None)
function _pandas_dtype(vals)
    all(_is_none, vals) && return :object
    all(x -> x isa Bool, vals) && return :bool
    all(_is_int, vals) && return :int
    all(x -> _is_none(x) || _is_int(x) || x isa AbstractFloat, vals) && return :float
    return :object
end

function _pandas_cell(x, dtype::Symbol)
    _is_none(x) && return ""
    x isa Bool && return x ? "True" : "False"
    if x isa Real && dtype === :float || x isa AbstractFloat
        y = Float64(x)
        isnan(y) && return ""
        isinf(y) && return y > 0 ? "inf" : "-inf"
        return pyfloat(y)
    end
    return string(x)
end

"""
    pandas_csv(path, rows)

Write `rows` (dicts from column name to value; `nothing` is Python's None) as
`pd.DataFrame(rows).to_csv(path, index=False)`. Columns appear in first-seen order; a key missing from a
row is NaN there.
"""
function pandas_csv(path::AbstractString, rows)
    cols = String[]
    for r in rows, k in keys(r)
        String(k) in cols || push!(cols, String(k))
    end
    vals = [[get(r, c, missing) for r in rows] for c in cols]
    types = [_pandas_dtype(v) for v in vals]
    cells = [Any[_pandas_cell(vals[j][i], types[j]) for j in eachindex(cols)] for i in eachindex(rows)]
    write_csv_table(path, CSVTable(cols, cells, "\n"))
end
