# Two pieces of CPython behaviour that decide emitted circuits:
#  - `group_checks` shuffles with `random.Random(0)` (MT19937, CPython's init_by_array seeding
#    and `_randbelow` rejection sampling);
#  - `_reroute_degree!` iterates a `set` of qubit ids, whose order is CPython's hash-table order.
# Both are deterministic and replayed exactly here (checked against CPython through the bridge in the tests).

mutable struct PyRandom
    mt::Vector{UInt32}
    i::Int
end

function _init_genrand(s::UInt32)
    mt = zeros(UInt32, 624)
    mt[1] = s
    for i in 2:624
        mt[i] = UInt32(1812433253) * (mt[i-1] ⊻ (mt[i-1] >> 30)) + UInt32(i - 1)
    end
    return mt
end

"""
    PyRandom(seed)

Python's `random.Random(seed)` for a non-negative integer seed.
"""
function PyRandom(seed::Integer)
    seed >= 0 || throw(ArgumentError("seed must be non-negative"))
    key = UInt32[]
    s = BigInt(seed)
    while s > 0
        push!(key, UInt32(s & 0xffffffff))
        s >>= 32
    end
    isempty(key) && push!(key, 0)
    mt = _init_genrand(UInt32(19650218))
    N = 624
    i, j = 1, 0
    for _ in 1:max(N, length(key))
        mt[i+1] = (mt[i+1] ⊻ ((mt[i] ⊻ (mt[i] >> 30)) * UInt32(1664525))) + key[j+1] + UInt32(j)
        i += 1
        j += 1
        if i >= N
            mt[1] = mt[N]
            i = 1
        end
        j >= length(key) && (j = 0)
    end
    for _ in 1:N-1
        mt[i+1] = (mt[i+1] ⊻ ((mt[i] ⊻ (mt[i] >> 30)) * UInt32(1566083941))) - UInt32(i)
        i += 1
        if i >= N
            mt[1] = mt[N]
            i = 1
        end
    end
    mt[1] = 0x80000000
    return PyRandom(mt, 624)
end

function _genrand!(r::PyRandom)
    mt = r.mt
    if r.i >= 624
        for k in 1:624
            y = (mt[k] & 0x80000000) | (mt[mod1(k + 1, 624)] & 0x7fffffff)
            mt[k] = mt[mod1(k + 397, 624)] ⊻ (y >> 1) ⊻ (isodd(y) ? 0x9908b0df : 0x00000000)
        end
        r.i = 0
    end
    y = mt[r.i+1]
    r.i += 1
    y ⊻= y >> 11
    y ⊻= (y << 7) & 0x9d2c5680
    y ⊻= (y << 15) & 0xefc60000
    y ⊻= y >> 18
    return y
end

function _randbelow!(r::PyRandom, n::Int)
    k = 64 - leading_zeros(UInt64(n))
    k <= 32 || throw(ArgumentError("n too large"))
    x = Int(_genrand!(r) >> (32 - k))
    while x >= n
        x = Int(_genrand!(r) >> (32 - k))
    end
    return x
end

"""
    pyshuffle!(r, x)

`random.Random.shuffle` (Python 3.12).
"""
function pyshuffle!(r::PyRandom, x::AbstractVector)
    for i in length(x)-1:-1:1
        j = _randbelow!(r, i + 1)
        x[i+1], x[j+1] = x[j+1], x[i+1]
    end
    return x
end

"""
    pyset_order(xs)

Iteration order of `set(xs)` for non-negative ints `xs` in CPython 3.12 (hash(i) == i):
open addressing with 9 linear probes, then perturbed probing; table resized to 4x used at 3/5 fill.
"""
function pyset_order(xs)
    table = fill(-1, 8)
    used = 0
    fill_ = 0
    function insert_clean!(tab, key)
        mask = length(tab) - 1
        perturb = UInt(key)
        i = UInt(key) & mask
        while true
            tab[i+1] == -1 && (tab[i+1] = key; return)
            if i + 9 <= mask
                for j in 1:9
                    tab[i+j+1] == -1 && (tab[i+j+1] = key; return)
                end
            end
            perturb >>= 5
            i = (i * 5 + 1 + perturb) & mask
        end
    end
    for key in xs
        key >= 0 || throw(ArgumentError("only non-negative ints"))
        mask = length(table) - 1
        perturb = UInt(key)
        i = UInt(key) & mask
        found = false
        placed = false
        while !found && !placed
            probes = i + 9 <= mask ? 9 : 0
            e = i
            for _ in 0:probes
                if table[e+1] == -1
                    table[e+1] = key
                    placed = true
                    break
                elseif table[e+1] == key
                    found = true
                    break
                end
                e += 1
            end
            if !found && !placed
                perturb >>= 5
                i = (i * 5 + 1 + perturb) & mask
            end
        end
        found && continue
        used += 1
        fill_ += 1
        if fill_ * 5 >= mask * 3
            minused = used > 50000 ? used * 2 : used * 4
            newsize = 8
            while newsize <= minused
                newsize <<= 1
            end
            old = table
            table = fill(-1, newsize)
            for k in old
                k == -1 || insert_clean!(table, k)
            end
            fill_ = used
        end
    end
    return [k for k in table if k != -1]
end
