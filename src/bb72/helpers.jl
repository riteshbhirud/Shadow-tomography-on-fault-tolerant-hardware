# Small pieces the BB72 study scripts share: itertools orders, numpy file I/O and statistics that
# have to match numpy/scipy bit for bit, and a pandas-style CSV writer.

"""Size-k combinations of `v` in itertools.combinations order."""
function py_combinations(v::AbstractVector, k::Integer)
    n = length(v)
    out = Vector{Vector{eltype(v)}}()
    k > n && return out
    idx = collect(1:k)
    while true
        push!(out, v[idx])
        i = k
        while i >= 1 && idx[i] == n - k + i
            i -= 1
        end
        i == 0 && return out
        idx[i] += 1
        for j in i+1:k
            idx[j] = idx[j-1] + 1
        end
    end
end

"""All orderings of 0..n-1 in itertools.permutations order (lexicographic)."""
function py_permutations(n::Integer)
    out = Vector{Vector{Int}}()
    p = collect(0:n-1)
    while true
        push!(out, copy(p))
        i = n - 1
        while i >= 1 && p[i] > p[i+1]
            i -= 1
        end
        i < 1 && return out
        j = n
        while p[j] < p[i]
            j -= 1
        end
        p[i], p[j] = p[j], p[i]
        reverse!(p, i + 1, n)
    end
end

"""`np.argsort(x, kind="stable")`, 1-based. Compares with `<`, so -0.0 and 0.0 tie as in numpy."""
argsort_stable(x) = sortperm(x; lt = <, alg = MergeSort)

"""Python's float ** float (the C library's pow)."""
cpow(x::Float64, y::Float64) = @ccall pow(x::Float64, y::Float64)::Float64

"""`scipy.stats.rankdata(x)`: ranks 1..n, ties get their average rank."""
function rankdata(x::AbstractVector)
    n = length(x)
    p = argsort_stable(x)
    r = zeros(n)
    i = 1
    while i <= n
        j = i
        while j < n && x[p[j+1]] == x[p[i]]
            j += 1
        end
        for t in i:j
            r[p[t]] = (i + j) / 2
        end
        i = j + 1
    end
    return r
end

"""
    auc(score, err)

P(score of a random error shot > score of a random correct shot), ties counting 1/2 (Mann-Whitney
from average ranks); `nothing` if either class is empty.
"""
function auc(score::AbstractVector, err::AbstractVector{Bool})
    r = rankdata(score)
    n1, n0 = count(err), count(!, err)
    (n1 > 0 && n0 > 0) || return nothing
    return (np_sum(r[err]) - n1 * (n1 + 1) / 2) / (n1 * n0)
end

"""`np.corrcoef(a, b)[0, 1]` computed by numpy through the bridge (its BLAS dot products decide the last bits)."""
function np_corrcoef(a::AbstractVector, b::AbstractVector)
    np = pymodule(:numpy)
    return pyconvert(Float64, np.corrcoef(to_numpy(Int.(a)), to_numpy(Int.(b)))[0, 1])
end

"""Mean of 0/1 values as numpy computes it (the sum is exact, then one division)."""
bitmean(x) = count(!iszero, x) / length(x)

#### numpy files, written and read by numpy through the bridge so the formats are the reference's

"""A C-ordered numpy copy of a Julia array."""
to_numpy(A::AbstractArray) = (np = pymodule(:numpy); np.ascontiguousarray(np.asarray(A)))

"""`np.save(path, A)`."""
npy_save(path::AbstractString, A::AbstractArray) = (pymodule(:numpy).save(path, to_numpy(A)); nothing)

"""`np.load(path)` of a 2-d uint8 array."""
npy_load_u8(path::AbstractString) = pyconvert(Matrix{UInt8}, pymodule(:numpy).load(path))

"""
    npz_save_compressed(path; arrays...)

`np.savez_compressed(path, ...)`. Julia arrays become numpy arrays; dicts become 0-d object arrays,
as `np.array(d, dtype=object)`.
"""
function npz_save_compressed(path::AbstractString; arrays...)
    np = pymodule(:numpy)
    conv(v::AbstractArray) = to_numpy(v)
    conv(v::AbstractDict) = np.array(pydict(v); dtype = pybuiltins.object)
    np.savez_compressed(path; (k => conv(v) for (k, v) in arrays)...)
    return nothing
end

# Python values of a stored dict back to Julia, keeping key order.
_jl(x::Py) = pyisinstance(x, pybuiltins.dict) ?
             OrderedDict{String,Any}(pyconvert(String, k) => _jl(v) for (k, v) in x.items()) :
             pyis(x, pybuiltins.None) ? nothing : pyisinstance(x, pybuiltins.bool) ? pyconvert(Bool, x) :
             pyisinstance(x, pybuiltins.int) ? pyconvert(Int, x) : pyisinstance(x, pybuiltins.float) ?
             pyconvert(Float64, x) : pyconvert(Any, x)

"""
    npz_load(path) -> Dict

The arrays of an `.npz` file: 1-d arrays as Julia vectors (`uint8` and `float64`), 0-d object arrays
as the dicts they hold.
"""
function npz_load(path::AbstractString)
    z = pymodule(:numpy).load(path; allow_pickle = true)
    out = Dict{String,Any}()
    for k in z.files
        a = z[k]
        key = pyconvert(String, k)
        if pyconvert(Int, a.ndim) == 0
            out[key] = _jl(a.item())
        elseif pyconvert(String, a.dtype.name) == "uint8"
            out[key] = pyconvert(Vector{UInt8}, a)
        else
            out[key] = pyconvert(Vector{Float64}, a)
        end
    end
    z.close()
    return out
end

"""Write `path` through a temporary file and rename it, so a crash never leaves a half-written chunk."""
function atomic_save(save, path::AbstractString, tmp::AbstractString)
    save(tmp)
    mv(tmp, path; force = true)
    return path
end

#### pandas DataFrame(rows).to_csv(index=False)

_frame_cell(x::Bool) = x ? "True" : "False"
_frame_cell(x::AbstractFloat) = isnan(x) ? missing : Float64(x)
_frame_cell(::Nothing) = missing
_frame_cell(x::AbstractChar) = string(x)
_frame_cell(x) = x

"""
    write_frame_csv(path, rows)

Write a vector of ordered dicts as `pd.DataFrame(rows).to_csv(path, index=False)`: columns in order of
first appearance, NaN and None empty, booleans as True/False, floats as Python's repr.
"""
function write_frame_csv(path::AbstractString, rows)
    cols = String[]
    for r in rows, k in keys(r)
        k in cols || push!(cols, k)
    end
    body = [Any[_frame_cell(get(r, c, nothing)) for c in cols] for r in rows]
    write_csv_table(path, CSVTable(cols, body, "\n"))
end

"""Characters as strings and tuples as vectors, for JSON output of compiled settings."""
jsonable(x::AbstractChar) = string(x)
jsonable(x::Tuple) = Any[jsonable(y) for y in x]
jsonable(x::AbstractVector) = Any[jsonable(y) for y in x]
jsonable(x::AbstractDict) = OrderedDict{String,Any}(string(k) => jsonable(v) for (k, v) in x)
jsonable(x) = x

"""Write `x` like `json.dump(x, f, indent=1)`."""
write_json(path::AbstractString, x; indent = 1) = write(path, pyjson(jsonable(x); indent))

"""Read a JSON file keeping key order (`JSON.Object`), or an empty ordered dict if it does not exist."""
read_json_or_empty(path::AbstractString) = isfile(path) ? JSON.parsefile(path) : OrderedDict{String,Any}()
