# numpy and Python float conventions (pairwise sums, the two rounding rules) used for the datasets' expected
# values, so recomputed values agree with the stored tables to the last bit.

"""Sum in numpy's pairwise order (`np.add.reduce` on float64)."""
function np_sum(a::AbstractVector{Float64})
    n = length(a)
    if n < 8
        r = -0.0
        for x in a
            r += x
        end
        return r
    elseif n <= 128
        r = a[1:8]
        i = 9
        while i <= n - n % 8
            for j in 0:7
                r[j+1] += a[i+j]
            end
            i += 8
        end
        res = ((r[1] + r[2]) + (r[3] + r[4])) + ((r[5] + r[6]) + (r[7] + r[8]))
        while i <= n
            res += a[i]
            i += 1
        end
        return res
    end
    n2 = n ÷ 2
    n2 -= n2 % 8
    return np_sum(a[1:n2]) + np_sum(a[n2+1:end])
end

"""`np.mean` of a list of floats."""
np_mean(a) = np_sum(Float64.(collect(a))) / length(a)

"""`np.std(a, ddof=ddof)`: two-pass variance with pairwise sums, then the square root."""
function np_std(a::AbstractVector{Float64}; ddof::Integer = 0)
    m = np_sum(a) / length(a)
    x = a .- m
    return sqrt(np_sum(x .* x) / max(length(a) - ddof, 0))
end

"""`round(x, k)` on a numpy float64: multiply by 10^k, round half to even (rint), divide."""
np_round(x::Float64, k::Integer) = round(x * 10.0^k) / 10.0^k

"""Python's `round(x, k)` on a float: the correctly rounded decimal."""
py_round(x::Float64, k::Integer) = isfinite(x) ? parse(Float64, @sprintf("%.*f", k, x)) : x

const Z95 = 1.959963984540054
