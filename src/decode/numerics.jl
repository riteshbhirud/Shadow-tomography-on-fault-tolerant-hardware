# Floating-point steps done exactly as numpy does them, so that weights and gaps agree with the
# reference bit for bit, plus the conversions to the numpy/scipy objects the decoders take.

using SparseArrays

# numpy's float64 log calls the C library; Julia's own `log` differs from it in about 1e-4 of cases.
_clog(x::Float64) = @ccall log(x::Float64)::Float64

"""
    llr(p)

Log-likelihood weight `ln((1 - p) / p)` of a mechanism with probability `p`, computed as numpy does.
"""
llr(p::Float64) = _clog((1 - p) / p)

# numpy's pairwise summation (pairwise_sum_DOUBLE, unrolled by 8, blocks of 128).
function _pairwise(a::AbstractVector{Float64}, lo::Int, n::Int)
    if n < 8
        res = 0.0
        for i in 0:(n - 1)
            res += a[lo + i]
        end
        return res
    elseif n <= 128
        r1, r2, r3, r4, r5, r6, r7, r8 = a[lo], a[lo + 1], a[lo + 2], a[lo + 3], a[lo + 4], a[lo + 5], a[lo + 6], a[lo + 7]
        i = 8
        while i < n - n % 8
            r1 += a[lo + i]
            r2 += a[lo + i + 1]
            r3 += a[lo + i + 2]
            r4 += a[lo + i + 3]
            r5 += a[lo + i + 4]
            r6 += a[lo + i + 5]
            r7 += a[lo + i + 6]
            r8 += a[lo + i + 7]
            i += 8
        end
        res = ((r1 + r2) + (r3 + r4)) + ((r5 + r6) + (r7 + r8))
        while i < n
            res += a[lo + i]
            i += 1
        end
        return res
    else
        n2 = n ÷ 2
        n2 -= n2 % 8
        return _pairwise(a, lo, n2) + _pairwise(a, lo + n2, n - n2)
    end
end

"""
    npsum(a)

Sum of a contiguous `Float64` vector with numpy's summation order (`a.sum()`, or one row of a
C-ordered `.sum(axis=1)`).
"""
npsum(a::AbstractVector{Float64}) = 0.0 + _pairwise(a, firstindex(a), length(a))

# np.minimum / np.maximum: NaN propagates, and on ties the first argument is returned.
_npmin(a::Float64, b::Float64) = (a <= b || isnan(a)) ? a : b
_npmax(a::Float64, b::Float64) = (a >= b || isnan(a)) ? a : b

# np.isclose with the default tolerances
_isclose(a::Float64, b::Float64) = isfinite(a) && isfinite(b) ? abs(a - b) <= 1e-8 + 1e-5 * abs(b) : a == b

# Python's builtin max(a, b): b only if it is strictly larger, so a NaN in second place is ignored.
_pymax(a, b) = b > a ? b : a

# Weight of a correction e (0/1 entries): llr[e.astype(bool)].sum()
_weight(llrs::Vector{Float64}, e::AbstractVector) = npsum(llrs[e .!= 0])

#### numpy / scipy objects

_np() = pymodule(:numpy)

"""
    to_numpy_u8(M)

C-ordered numpy `uint8` copy of a 0/1 matrix or vector (the dtype and layout the decoders expect).
"""
to_numpy_u8(M::AbstractArray) = _np().ascontiguousarray(_np().asarray(convert(Array{UInt8}, M)))

_from_numpy_u8(x::Py, dims...) = pyconvert(Array{UInt8}, _np().asarray(x; dtype = _np().uint8))

# The decoders get canonical scipy matrices (sorted indices, int32 index arrays, uint8 ones), as
# `lil -> tocsc` and `sp.vstack(...).tocsr()/.tocsc()` produce.
function _scipy_sparse(M::SparseMatrixCSC, format::Symbol)
    np = _np()
    S = format === :csc ? M : sparse(transpose(M))  # CSR of M = CSC of its transpose
    indptr = np.array(Int32.(S.colptr .- 1); dtype = np.int32)
    indices = np.array(Int32.(S.rowval .- 1); dtype = np.int32)
    data = np.ones(length(S.rowval); dtype = np.uint8)
    ctor = format === :csc ? pymodule(Symbol("scipy.sparse")).csc_matrix : pymodule(Symbol("scipy.sparse")).csr_matrix
    return ctor(pytuple((data, indices, indptr)); shape = pytuple(size(M)))
end

# Parity of H * e (e has 0/1 entries); `syn` is compared against it.
function _syndrome_ok(H::SparseMatrixCSC, e::AbstractVector, syn::AbstractVector)
    s = zeros(UInt8, size(H, 1))
    rows = rowvals(H)
    for j in 1:size(H, 2)
        e[j] == 0 && continue
        for k in nzrange(H, j)
            s[rows[k]] ⊻= 0x01
        end
    end
    return all(s[i] == syn[i] for i in eachindex(s))
end

