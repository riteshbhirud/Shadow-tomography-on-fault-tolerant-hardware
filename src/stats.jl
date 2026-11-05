# Small statistics helpers.

using SpecialFunctions: erfc

# Cephes ndtri, with the fused multiply-adds of scipy's build, so the stored Wilson intervals are
# reproduced to the last bit.
const _NDTRI_P0 = (-5.99633501014107895267E1, 9.80010754185999661536E1, -5.66762857469070293439E1,
    1.39312609387279679503E1, -1.23916583867381258016E0)
const _NDTRI_Q0 = (1.95448858338141759834E0, 4.67627912898881538453E0, 8.63602421390890590575E1,
    -2.25462687854119370527E2, 2.00260212380060660359E2, -8.20372256168333339912E1, 1.59056225126211695515E1,
    -1.18331621121330003142E0)
const _NDTRI_P1 = (4.05544892305962419923E0, 3.15251094599893866154E1, 5.71628192246421288162E1,
    4.40805073893200834700E1, 1.46849561928858024014E1, 2.18663306850790267539E0, -1.40256079171354495875E-1,
    -3.50424626827848203418E-2, -8.57456785154685413611E-4)
const _NDTRI_Q1 = (1.57799883256466749731E1, 4.53907635128879210584E1, 4.13172038254672030440E1,
    1.50425385692907503408E1, 2.50464946208309415979E0, -1.42182922854787788574E-1, -3.80806407691578277194E-2,
    -9.33259480895457427372E-4)
const _NDTRI_P2 = (3.23774891776946035970E0, 6.91522889068984211695E0, 3.93881025292474443415E0,
    1.33303460815807542389E0, 2.01485389549179081538E-1, 1.23716634817820021358E-2, 3.01581553508235416007E-4,
    2.65806974686737550832E-6, 6.23974539184983293730E-9)
const _NDTRI_Q2 = (6.02427039364742014255E0, 3.67983563856160859403E0, 1.37702099489081330271E0,
    2.16236993594496635890E-1, 1.34204006088543189037E-2, 3.28014464682127739104E-4, 2.89247864745380683936E-6,
    6.79019408009981274425E-9)

_polevl(x, c) = foldl((a, ci) -> fma(a, x, ci), Base.tail(c); init = c[1])
_p1evl(x, c) = foldl((a, ci) -> fma(a, x, ci), Base.tail(c); init = x + c[1])

"""
    norm_ppf(y)

Standard normal quantile (Cephes `ndtri`).
"""
function norm_ppf(y0::Float64)
    y0 == 0 && return -Inf
    y0 == 1 && return Inf
    0 < y0 < 1 || return NaN
    lower, y = true, y0
    if y > 1.0 - 0.13533528323661269189  # exp(-2)
        y, lower = 1.0 - y, false
    end
    if y > 0.13533528323661269189
        y -= 0.5
        y2 = y * y
        return fma(y, y2 * _polevl(y2, _NDTRI_P0) / _p1evl(y2, _NDTRI_Q0), y) * 2.50662827463100050242E0
    end
    x = sqrt(-2.0 * _clog(y))
    x0 = x - _clog(x) / x
    z = 1.0 / x
    x1 = x < 8.0 ? z * _polevl(z, _NDTRI_P1) / _p1evl(z, _NDTRI_Q1) : z * _polevl(z, _NDTRI_P2) / _p1evl(z, _NDTRI_Q2)
    x = x0 - x1
    return lower ? -x : x
end

"""
    wilson(k, n; conf=0.95) -> (lo, hi)

Wilson score interval for `k` successes in `n` trials; `(NaN, NaN)` if `n = 0`.
"""
function wilson(k::Integer, n::Integer; conf::Real = 0.95)
    n == 0 && return (NaN, NaN)
    z = norm_ppf(0.5 + conf / 2)
    ph = k / n
    den = 1 + z * z / n
    c = (ph + z * z / (2 * n)) / den
    h = z * sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / den
    return (max(0.0, c - h), min(1.0, c + h))
end

"""
    two_proportion_test(k1, n1, k2, n2) -> (z, p)

Pooled two-proportion z-test of k1/n1 against k2/n2; two-sided p-value (1 if both rates are 0 or 1).
"""
function two_proportion_test(k1::Integer, n1::Integer, k2::Integer, n2::Integer)
    pp = (k1 + k2) / (n1 + n2)
    se = sqrt(pp * (1 - pp) * (1 / n1 + 1 / n2))
    se == 0 && return (0.0, 1.0)
    z = (k1 / n1 - k2 / n2) / se
    return (z, erfc(abs(z) / sqrt(2)))
end

"""
    holm(pvals; alpha=0.05) -> Vector{Bool}

Holm step-down rejections at family-wise level `alpha`, in the order of `pvals`.
"""
function holm(pvals::AbstractVector{<:Real}; alpha::Real = 0.05)
    m = length(pvals)
    reject = falses(m)
    for (i, k) in enumerate(sortperm(pvals))
        pvals[k] <= alpha / (m - i + 1) || break
        reject[k] = true
    end
    return reject
end
